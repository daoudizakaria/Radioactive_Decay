# Radioactive Decay

[![tests](https://github.com/daoudizakaria/Radioactive_Decay/actions/workflows/tests.yml/badge.svg)](https://github.com/daoudizakaria/Radioactive_Decay/actions/workflows/tests.yml)

Python programs that simulate radioactive decay and compute nuclear energies from evaluated nuclear
data. They accompany introductory notes on nuclear physics written for students in the final years of
secondary school and the first year of university,
[`notes/nuclear_physics_notes.pdf`](notes/nuclear_physics_notes.pdf), which cover the structure of the
nucleus, radioactive decay and decay chains, binding energy, fission and fusion.

## Programs

| Program | Description |
|---|---|
| `radioactivity.py` | Decay of a single nuclide. The decay law dN/dt = −λN is integrated with Euler's method and compared with the exact solution N(t) = N₀ exp(−λt). |
| `nuclear.py` | Four interactive simulations: the decay of a nuclide and its activity A = λN; a parent–daughter chain, for example U-238 → Th-234 or Th-232 → Ra-228, which reaches secular equilibrium; a chain with branching; and Monte Carlo decay, in which each nucleus decays at random. The results can be exported to CSV files. |
| `nuclear_energy.py` | Binding energies from measured atomic masses, a least-squares fit of the semi-empirical mass formula, and the energy released (Q-value) in alpha and beta decay, fission and fusion. |
| `make_figures.py` | Reproduces the figures of the notes in `notes/figures/`. |

## Methods

- **Single decay.** Euler's method is first-order: at fixed time its relative error is proportional to
  the time step. The programs compare it with the exact exponential.
- **Decay chains.** The populations are given by the exact Bateman solution. Parent and daughter
  half-lives often differ by many orders of magnitude (4.5 × 10⁹ years for U-238 against 24 days for
  Th-234), and Euler's method becomes unstable as soon as λ_D Δt > 2.
- **Monte Carlo decay.** During a step Δt each remaining nucleus decays independently with probability
  p = 1 − exp(−λΔt), so the number of decays is drawn from a binomial distribution. The mean follows
  N₀ exp(−λt), with fluctuations σ = [N₀ q (1 − q)]^½, where q = exp(−λt).
- **Nuclear energies.** Binding energies and Q-values are computed from atomic mass excesses, so the
  electron masses cancel except in β⁺ decay, where 2mₑc² is subtracted.

## Data

- `nuclides_data.py` and `nuclides.csv`: half-lives of the nuclides offered by the programs, from
  NUBASE2020.
- `data/ame2020_masses.csv`: measured atomic mass excesses of 2550 nuclides from AME2020
  (M. Wang et al., *Chinese Physics C* **45**, 030003 (2021)).
- `data/nubase2020_ground_states.csv`: ground-state half-lives and main decay modes of 3558 nuclides
  from NUBASE2020 (F. G. Kondev et al., *Chinese Physics C* **45**, 030001 (2021)).

## Usage

```bash
pip install -r requirements.txt
python3 radioactivity.py
python3 nuclear.py
python3 nuclear_energy.py
```

In a Jupyter notebook, `nuclear.py` also provides sliders for the single-decay simulation (this requires
the `ipywidgets` package):

```python
from nuclear import interactive_widget
interactive_widget()
```

## Tests

```bash
python3 -m unittest discover tests
```

The tests check the programs against exact results: the half-life, the first-order convergence of
Euler's method, the mean and the fluctuations of the Monte Carlo simulation, the Bateman solution
against a fine numerical integration, secular equilibrium and branching ratios, reference binding
energies and Q-values from AME2020, the fitted coefficients of the semi-empirical mass formula, and the
consistency of the two nuclide tables.

## Notes

The notes are written in LaTeX (`notes/nuclear_physics_notes.tex`). To rebuild the PDF:

```bash
python3 make_figures.py
cd notes && latexmk -pdf nuclear_physics_notes.tex
```

## Licence

GNU General Public License v3.0; see [LICENSE](LICENSE).
