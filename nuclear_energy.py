"""
Nuclear masses, binding energies and reaction energies (Q-values).

All masses are atomic mass excesses from AME2020 (data/ame2020_masses.csv):
    Delta = (M_atom - A u) c^2.
Working with atomic masses means that the electrons are counted automatically,
except in beta-plus decay (see q_value).

Run  python3 nuclear_energy.py  to print binding energies, fit the
semi-empirical mass formula and compute the energy released by decays,
fission and fusion.
"""
import os
import re
from functools import lru_cache

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

U_MEV = 931.49410242          # 1 atomic mass unit in MeV/c^2 (CODATA 2018)
M_ELECTRON_MEV = 0.51099895   # electron mass in MeV/c^2 (CODATA 2018)


def load_masses():
    """Measured atomic mass excesses (keV) of AME2020, one row per nuclide (Z, N, A, element)."""
    return pd.read_csv(os.path.join(HERE, "data", "ame2020_masses.csv"), comment="#")


@lru_cache(maxsize=None)
def _tables():
    """Mass excesses in MeV indexed by (Z, A), and the atomic number of each element symbol."""
    df = load_masses()
    masses = {(int(r.Z), int(r.A)): r.mass_excess_keV / 1000.0 for r in df.itertuples()}
    elements = {r.element: int(r.Z) for r in df.itertuples() if r.Z > 0}
    return masses, elements


def parse_nuclide(name):
    """
    Convert a nuclide name into (Z, A).
    Accepted forms: "238U", "U238", "U-238", "4He", "n" (neutron), "p" (hydrogen atom),
    "d" (deuterium), "t" (tritium).
    """
    _, elements = _tables()
    s = name.strip().replace("-", "")
    special = {"n": (0, 1), "p": (1, 1), "d": (1, 2), "t": (1, 3)}
    if s.lower() in special:
        return special[s.lower()]
    m = re.fullmatch(r"(\d+)([A-Za-z]+)", s) or re.fullmatch(r"([A-Za-z]+)(\d+)", s)
    if m is None:
        raise ValueError(f"Cannot read the nuclide name '{name}'")
    a, b = m.groups()
    A, symbol = (int(a), b) if a.isdigit() else (int(b), a)
    symbol = symbol.capitalize()
    if symbol not in elements:
        raise ValueError(f"Unknown element '{symbol}'")
    return elements[symbol], A


def mass_excess(name):
    """Atomic mass excess in MeV."""
    masses, _ = _tables()
    Z, A = parse_nuclide(name)
    if (Z, A) not in masses:
        raise ValueError(f"No measured mass for {name} in AME2020")
    return masses[(Z, A)]


def atomic_mass(name):
    """Atomic mass in atomic mass units (u)."""
    Z, A = parse_nuclide(name)
    return A + mass_excess(name) / U_MEV


def binding_energy(name):
    """
    Total binding energy B (MeV) of the nucleus:
        B = [Z M(1H) + N m_n - M(atom)] c^2
    (the small binding energy of the electrons is neglected).
    """
    Z, A = parse_nuclide(name)
    N = A - Z
    return Z * mass_excess("1H") + N * mass_excess("n") - mass_excess(name)


def binding_energy_table():
    """DataFrame of all measured nuclei with their binding energy per nucleon (MeV)."""
    df = load_masses()
    df = df[df.A > 1].copy()
    df["B_MeV"] = (df.Z * mass_excess("1H") + df.N * mass_excess("n")
                   - df.mass_excess_keV / 1000.0)
    df["B_per_A_MeV"] = df.B_MeV / df.A
    return df


# =============================================================================
# Semi-empirical mass formula (liquid drop model)
# =============================================================================
def semf_terms(Z, A):
    """
    The five terms of the Bethe-Weizsaecker formula, without their coefficients:
        B = aV A - aS A^(2/3) - aC Z(Z-1)/A^(1/3) - aA (A-2Z)^2/A + delta
    with delta = +aP/A^(1/2) (even-even), 0 (odd A), -aP/A^(1/2) (odd-odd).
    """
    Z = np.asarray(Z, dtype=float)
    A = np.asarray(A, dtype=float)
    N = A - Z
    pairing = np.where(A % 2 == 1, 0.0, np.where(Z % 2 == 0, 1.0, -1.0)) / np.sqrt(A)
    return np.stack([A, -A**(2/3), -Z*(Z - 1)/A**(1/3), -(N - Z)**2/A, pairing], axis=-1)


def fit_semf(A_min=20):
    """
    Least-squares fit of the five coefficients (aV, aS, aC, aA, aP) in MeV to the
    measured binding energies of all nuclei with A >= A_min.
    """
    df = binding_energy_table()
    df = df[df.A >= A_min]
    X = semf_terms(df.Z.values, df.A.values)
    coef, *_ = np.linalg.lstsq(X, df.B_MeV.values, rcond=None)
    rms = np.sqrt(np.mean((X @ coef - df.B_MeV.values)**2))
    return dict(zip(["aV", "aS", "aC", "aA", "aP"], coef)), rms


def semf_binding_energy(Z, A, coef):
    """Binding energy (MeV) predicted by the semi-empirical mass formula."""
    c = np.array([coef[k] for k in ["aV", "aS", "aC", "aA", "aP"]])
    return semf_terms(Z, A) @ c


# =============================================================================
# Reaction energies
# =============================================================================
def q_value(reactants, products, beta_plus=0):
    """
    Energy released (MeV) by a nuclear reaction or decay:
        Q = [sum of initial masses - sum of final masses] c^2
    computed with atomic masses. With atomic masses the electrons balance
    automatically in alpha decay, beta-minus decay, electron capture, fission
    and fusion. In beta-plus decay the daughter atom has one electron too many
    and a positron is emitted: each beta-plus decay costs 2 m_e c^2, so give
    the number of beta-plus decays in beta_plus.
    Q > 0: energy is released. Q < 0: energy must be supplied.
    """
    q = sum(mass_excess(x) for x in reactants) - sum(mass_excess(x) for x in products)
    # mass numbers must balance, so the A u terms cancel
    A_in = sum(parse_nuclide(x)[1] for x in reactants)
    A_out = sum(parse_nuclide(x)[1] for x in products)
    if A_in != A_out:
        raise ValueError(f"Mass number is not conserved: {A_in} -> {A_out}")
    return q - 2 * M_ELECTRON_MEV * beta_plus


EXAMPLES = [
    ("Alpha decay of uranium 238", ["238U"], ["234Th", "4He"], 0),
    ("Alpha decay of polonium 210", ["210Po"], ["206Pb", "4He"], 0),
    ("Beta-minus decay of carbon 14", ["14C"], ["14N"], 0),
    ("Beta-minus decay of tritium", ["3H"], ["3He"], 0),
    ("Beta-plus decay of sodium 22", ["22Na"], ["22Ne"], 1),
    ("Fission 235U + n -> 141Ba + 92Kr + 3n", ["235U", "n"], ["141Ba", "92Kr", "n", "n", "n"], 0),
    ("Fusion D + T -> 4He + n", ["2H", "3H"], ["4He", "n"], 0),
    ("Fusion D + D -> 3He + n", ["2H", "2H"], ["3He", "n"], 0),
    ("Proton-proton chain 4 1H -> 4He (+ 2e+ + 2nu)", ["1H"] * 4, ["4He"], 0),
]


def main():
    print("Binding energies (AME2020)")
    print(f"{'nucleus':>8} {'B (MeV)':>10} {'B/A (MeV)':>10}")
    for x in ["2H", "4He", "12C", "16O", "56Fe", "62Ni", "208Pb", "235U", "238U"]:
        Z, A = parse_nuclide(x)
        B = binding_energy(x)
        print(f"{x:>8} {B:10.3f} {B / A:10.4f}")

    coef, rms = fit_semf()
    print("\nSemi-empirical mass formula fitted to all measured nuclei with A >= 20:")
    for k, v in coef.items():
        print(f"  {k} = {v:7.3f} MeV")
    print(f"  rms deviation = {rms:.2f} MeV")

    print("\nEnergy released (Q-value)")
    for label, r, p, bp in EXAMPLES:
        print(f"  {label:48s} Q = {q_value(r, p, bp):10.4f} MeV")

    df = binding_energy_table()
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(df.A, df.B_per_A_MeV, ".", ms=2, color="#0072B2", label="Measured (AME2020)")
    stable = df.groupby("A").apply(lambda g: g.loc[g.B_MeV.idxmax()], include_groups=False)
    A = stable.index.values
    ax.plot(A, semf_binding_energy(stable.Z.values, A, coef) / A, color="#D55E00",
            lw=1.5, label="Semi-empirical mass formula")
    ax.set_xlabel("Mass number A")
    ax.set_ylabel("Binding energy per nucleon B/A (MeV)")
    ax.set_ylim(0, 9.2)
    ax.grid(alpha=0.3)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    main()
