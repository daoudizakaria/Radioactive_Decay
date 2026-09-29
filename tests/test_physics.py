"""
Checks of the physics in the code.  Run from the repository folder with

    python3 -m unittest discover tests
"""
import math
import os
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nuclear import (simulate_decay, simulate_decay_chain, simulate_complex_decay_chain,
                     simulate_decay_monte_carlo, bateman_daughter)
from nuclear_energy import binding_energy, q_value, parse_nuclide, fit_semf
from nuclides_data import nuclides


class TestDecay(unittest.TestCase):

    def test_half_life(self):
        """After one half-life, half of the nuclei are left."""
        t, _, N, _ = simulate_decay(1000, 5700.0, 10, 5700.0)
        self.assertAlmostEqual(N[-1], 500.0, places=9)

    def test_euler_converges(self):
        """Euler's error at fixed time is proportional to the step size."""
        errors = []
        for L in (100, 200, 400):
            _, N, N_ex, _ = simulate_decay(1.0, 1.0, L, 5.0)
            errors.append(abs(N[-1] - N_ex[-1]) / N_ex[-1])
        self.assertAlmostEqual(errors[0] / errors[1], 2.0, delta=0.05)
        self.assertAlmostEqual(errors[1] / errors[2], 2.0, delta=0.05)

    def test_monte_carlo_average(self):
        """The average of many random simulations follows N0 exp(-λt)."""
        runs = [simulate_decay_monte_carlo(1000, 1.0, 50, 3.0, seed=k)[1] for k in range(400)]
        _, _, N_exp, sigma, _ = simulate_decay_monte_carlo(1000, 1.0, 50, 3.0, seed=0)
        mean = np.mean(runs, axis=0)
        std = np.std(runs, axis=0)
        # the mean of 400 runs is known to within sigma/20
        self.assertTrue(np.all(np.abs(mean - N_exp) <= 5 * sigma / 20 + 1e-9))
        self.assertTrue(np.allclose(std[1:], sigma[1:], rtol=0.15))

    def test_monte_carlo_is_reproducible(self):
        a = simulate_decay_monte_carlo(500, 1.0, 20, 2.0, seed=42)[1]
        b = simulate_decay_monte_carlo(500, 1.0, 20, 2.0, seed=42)[1]
        self.assertTrue(np.array_equal(a, b))


class TestChains(unittest.TestCase):

    def test_bateman_against_small_step_integration(self):
        """The exact solution agrees with a very fine numerical integration."""
        lP, lD = 0.3, 1.1
        t = np.linspace(0, 10, 200001)
        dt = t[1] - t[0]
        P, D = 1.0, 0.0
        for _ in range(len(t) - 1):
            # midpoint (second-order) integration
            P_mid = P - 0.5 * dt * lP * P
            D_mid = D + 0.5 * dt * (lP * P - lD * D)
            P, D = P - dt * lP * P_mid, D + dt * (lP * P_mid - lD * D_mid)
        self.assertAlmostEqual(D, bateman_daughter(1.0, lP, lD, t[-1]), places=8)

    def test_equal_decay_constants(self):
        """The special case λP = λD is the limit of the general formula."""
        t = np.array([0.5, 1.0, 2.0])
        general = bateman_daughter(1.0, 1.0, 1.0 + 1e-7, t)
        special = bateman_daughter(1.0, 1.0, 1.0, t)
        self.assertTrue(np.allclose(general, special, rtol=1e-6))

    def test_secular_equilibrium_uranium(self):
        """U-238 -> Th-234: the two activities become equal."""
        d = nuclides["238U"]
        t, NP, ND, l1, l2 = simulate_decay_chain(1e6, d["half_life"], d["daughter_half_life"],
                                                 1000, 5 * d["half_life"])
        self.assertTrue(np.all(ND >= 0))
        self.assertTrue(np.allclose(l2 * ND[1:], l1 * NP[1:], rtol=1e-9))

    def test_branching_ratios(self):
        """Each daughter receives its share of the parent decays."""
        t, NP, NA, NB, *_ = simulate_complex_decay_chain(1e6, 10.0, 1e-3, 1e-3, 0.3, 0.7, 100, 50.0)
        self.assertTrue(np.allclose(NA[1:] / NB[1:], 0.3 / 0.7))


class TestEnergy(unittest.TestCase):

    def test_parse_names(self):
        self.assertEqual(parse_nuclide("238U"), (92, 238))
        self.assertEqual(parse_nuclide("U-238"), (92, 238))
        self.assertEqual(parse_nuclide("14N"), (7, 14))
        self.assertEqual(parse_nuclide("n"), (0, 1))

    def test_binding_energies(self):
        """Reference values from AME2020."""
        self.assertAlmostEqual(binding_energy("2H"), 2.2246, places=3)
        self.assertAlmostEqual(binding_energy("4He") / 4, 7.0739, places=3)
        self.assertAlmostEqual(binding_energy("56Fe") / 56, 8.7904, places=3)

    def test_q_values(self):
        """Well-known energies released by decays and reactions."""
        self.assertAlmostEqual(q_value(["238U"], ["234Th", "4He"]), 4.270, places=2)
        self.assertAlmostEqual(q_value(["3H"], ["3He"]), 0.01859, places=4)
        self.assertAlmostEqual(q_value(["2H", "3H"], ["4He", "n"]), 17.589, places=2)
        self.assertAlmostEqual(q_value(["1H"] * 4, ["4He"]), 26.731, places=2)

    def test_semf_coefficients(self):
        """The fitted liquid drop coefficients have their usual sizes."""
        coef, rms = fit_semf()
        self.assertTrue(14 < coef["aV"] < 17)
        self.assertTrue(16 < coef["aS"] < 20)
        self.assertTrue(0.6 < coef["aC"] < 0.8)
        self.assertTrue(20 < coef["aA"] < 25)
        self.assertLess(rms, 5.0)


class TestData(unittest.TestCase):

    def test_csv_matches_dictionary(self):
        """nuclides.csv (used by radioactivity.py) agrees with nuclides_data.py."""
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        df = pd.read_csv(os.path.join(here, "nuclides.csv"))
        self.assertEqual(list(df["Nuclide Symbol"]), list(nuclides))
        for row in df.itertuples():
            self.assertAlmostEqual(row[3] / nuclides[row[1]]["half_life"], 1.0, places=3)


if __name__ == "__main__":
    unittest.main()
