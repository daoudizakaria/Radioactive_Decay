"""
Produce the figures of the notes (notes/figures/*.pdf) from the simulation code.

    python3 make_figures.py
"""
import math
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nuclear import simulate_decay, simulate_decay_monte_carlo, bateman_daughter
from nuclear_energy import binding_energy_table, fit_semf, semf_binding_energy

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "notes", "figures")

# Okabe-Ito colour-blind-safe palette
BLUE, ORANGE, GREEN, VERMILLION, PURPLE, SKY, YELLOW = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442")
GRAY = "#777777"

plt.rcParams.update({
    "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False,
    "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)
    print("  wrote", name)


def fig_decay_law():
    """Exponential decay with the successive half-lives marked."""
    x = np.linspace(0, 5, 400)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.plot(x, 2.0**(-x), color=BLUE, lw=2)
    for n in range(1, 5):
        ax.plot([n, n], [0, 2.0**(-n)], color=GRAY, lw=0.8, ls="--")
        ax.plot([0, n], [2.0**(-n), 2.0**(-n)], color=GRAY, lw=0.8, ls="--")
        ax.annotate(f"$N_0/{2**n}$", (0.05, 2.0**(-n)), xytext=(0, 3),
                    textcoords="offset points", fontsize=9)
    ax.set_xlabel(r"Time $t$ (in units of the half-life $T_{1/2}$)")
    ax.set_ylabel(r"Remaining fraction $N/N_0$")
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.02)
    save(fig, "decay_law")


def fig_euler():
    """Euler's method against the exact solution, and its error."""
    T = 1.0
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    t_ex = np.linspace(0, 5, 400)
    axs[0].plot(t_ex, 2.0**(-t_ex), color="black", lw=1.5, label="Exact")
    for L, c in [(8, VERMILLION), (25, ORANGE), (100, BLUE)]:
        t, N, _, _ = simulate_decay(1.0, T, L, 5 * T)
        axs[0].plot(t, N, "o-", ms=3, lw=1, color=c, label=f"Euler, {L} steps")
    axs[0].set_xlabel(r"$t/T_{1/2}$")
    axs[0].set_ylabel(r"$N/N_0$")
    axs[0].legend()
    axs[0].set_title("(a) Euler's method")

    steps = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
    err = []
    for L in steps:
        t, N, N_ex, _ = simulate_decay(1.0, T, L, 5 * T)
        err.append(abs(N[-1] - N_ex[-1]) / N_ex[-1])
    axs[1].loglog(5 * T / steps, err, "o-", color=BLUE)
    axs[1].loglog(5 * T / steps, 5 * T / steps * err[-1] / (5 * T / steps[-1]),
                  "--", color=GRAY, label=r"$\propto \Delta t$")
    axs[1].set_xlabel(r"Time step $\Delta t / T_{1/2}$")
    axs[1].set_ylabel(r"Relative error at $t = 5\,T_{1/2}$")
    axs[1].legend()
    axs[1].set_title("(b) Error of Euler's method")
    save(fig, "euler")


def fig_monte_carlo():
    """Random decays for a small and a large number of nuclei."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, N0, seed in [(axs[0], 100, 1), (axs[1], 10000, 2)]:
        for k, c in enumerate([BLUE, ORANGE, GREEN, PURPLE, VERMILLION]):
            t, N, N_exp, sigma, _ = simulate_decay_monte_carlo(N0, 1.0, 250, 5.0, seed=seed * 10 + k)
            ax.step(t, N / N0, where="post", color=c, lw=1)
        ax.fill_between(t, (N_exp - sigma) / N0, (N_exp + sigma) / N0, color=GRAY, alpha=0.25,
                        lw=0, label=r"$N_0 e^{-\lambda t} \pm \sigma$")
        ax.plot(t, N_exp / N0, "k--", lw=1.2, label=r"$N_0 e^{-\lambda t}$")
        ax.set_title(f"$N_0 = {N0}$ nuclei, 5 simulations")
        ax.set_xlabel(r"$t/T_{1/2}$")
        ax.legend()
    axs[0].set_ylabel(r"$N/N_0$")
    save(fig, "monte_carlo")


def fig_chains():
    """Secular (Th-232 -> Ra-228) and transient (Ba-140 -> La-140) equilibrium."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))

    # secular: Th-232 (1.40e10 y) -> Ra-228 (5.75 y)
    lP, lD = math.log(2) / 1.40e10, math.log(2) / 5.75
    t = np.linspace(0, 40, 400)
    A_P = lP * np.exp(-lP * t)
    A_D = lD * bateman_daughter(1.0, lP, lD, t)
    axs[0].plot(t, A_P / A_P[0], color=BLUE, lw=2, label=r"$^{232}$Th (parent)")
    axs[0].plot(t, A_D / A_P[0], color=ORANGE, lw=2, label=r"$^{228}$Ra (daughter)")
    axs[0].axvline(5.75, color=GRAY, lw=0.8, ls="--")
    axs[0].text(6.2, 0.1, r"$T_{1/2}(^{228}$Ra$)$", color=GRAY, fontsize=9)
    axs[0].set_xlabel("Time (years)")
    axs[0].set_ylabel("Activity / initial parent activity")
    axs[0].set_title("(a) Secular equilibrium")
    axs[0].legend(loc="center right")

    # transient: Ba-140 (12.7534 d) -> La-140 (40.289 h)
    TP, TD = 12.7534, 40.289 / 24
    lP, lD = math.log(2) / TP, math.log(2) / TD
    t = np.linspace(0, 60, 600)
    A_P = lP * np.exp(-lP * t)
    A_D = lD * bateman_daughter(1.0, lP, lD, t)
    axs[1].semilogy(t, A_P / A_P[0], color=BLUE, lw=2, label=r"$^{140}$Ba (parent)")
    axs[1].semilogy(t, A_D / A_P[0], color=ORANGE, lw=2, label=r"$^{140}$La (daughter)")
    axs[1].set_ylim(0.03, 2)
    axs[1].set_xlabel("Time (days)")
    axs[1].set_title("(b) Transient equilibrium")
    axs[1].legend()
    save(fig, "chains")


def fig_binding_energy():
    """Binding energy per nucleon: data and liquid drop model."""
    df = binding_energy_table()
    coef, _ = fit_semf()
    # most bound isobar for each A: the valley of stability
    best = df.loc[df.groupby("A").B_MeV.idxmax()]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.plot(df.A, df.B_per_A_MeV, ".", ms=2, color=SKY, label="All measured nuclei (AME2020)")
    ax.plot(best.A, best.B_per_A_MeV, ".", ms=4, color=BLUE, label="Most bound nucleus of each $A$")
    A = best.A.values[best.A.values >= 20]
    Z = best.Z.values[best.A.values >= 20]
    ax.plot(A, semf_binding_energy(Z, A, coef) / A, color=VERMILLION, lw=1.8,
            label="Liquid drop model (fit to $A \\geq 20$)")
    for name, Aa, Za, off in [("$^4$He", 4, 2, (8, -4)), ("$^{12}$C", 12, 6, (6, -14)),
                              ("$^{16}$O", 16, 8, (-6, 8)), ("$^{56}$Fe", 56, 26, (-8, 8)),
                              ("$^{238}$U", 238, 92, (-10, 8))]:
        r = df[(df.A == Aa) & (df.Z == Za)].iloc[0]
        ax.plot(Aa, r.B_per_A_MeV, "o", ms=5, mfc="white", mec="black", zorder=5)
        ax.annotate(name, (Aa, r.B_per_A_MeV), xytext=off, textcoords="offset points",
                    zorder=6)
    ax.annotate("", xy=(56, 4.6), xytext=(22, 4.6),
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=2))
    ax.text(22, 4.0, "fusion releases energy", color=GREEN)
    ax.annotate("", xy=(120, 4.6), xytext=(235, 4.6),
                arrowprops=dict(arrowstyle="->", color=PURPLE, lw=2))
    ax.text(128, 4.0, "fission releases energy", color=PURPLE)
    ax.set_xlabel("Mass number $A$")
    ax.set_ylabel("Binding energy per nucleon $B/A$ (MeV)")
    ax.set_ylim(0, 9.3)
    ax.set_xlim(0, 270)
    ax.legend(loc="lower right")
    save(fig, "binding_energy")


def fig_chart():
    """Chart of nuclides coloured by the main decay mode (NUBASE2020)."""
    df = pd.read_csv(os.path.join(HERE, "data", "nubase2020_ground_states.csv"), comment="#")
    groups = [
        ("Stable", ["stable"], "black", 7),
        (r"$\beta^-$", ["B-", "2B-"], SKY, 5),
        (r"$\beta^+$ / electron capture", ["B+", "EC", "2B+", "2EC"], VERMILLION, 5),
        (r"$\alpha$", ["A"], YELLOW, 5),
        ("Spontaneous fission", ["SF"], GREEN, 5),
        ("Proton or neutron emission", ["p", "2p", "n", "2n"], PURPLE, 5),
    ]
    fig, ax = plt.subplots(figsize=(8, 5.6))
    for label, modes, color, size in groups[1:] + groups[:1]:
        d = df[df.main_decay_mode.isin(modes)]
        ax.scatter(d.N, d.Z, s=size, marker="s", color=color, lw=0, label=label)
    ax.plot([0, 130], [0, 130], color=GRAY, lw=1, ls="--")
    ax.text(100, 104, "$N = Z$", color=GRAY)
    ax.set_xlabel("Number of neutrons $N$")
    ax.set_ylabel("Number of protons $Z$")
    ax.set_xlim(0, 180)
    ax.set_ylim(0, 120)
    ax.set_aspect("equal")
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(g[0]) for g in groups]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper left",
              markerscale=2.5, fontsize=9)
    save(fig, "chart")


def fig_carbon():
    """Carbon-14 dating: remaining fraction of C-14 as a function of age."""
    T = 5700.0  # years (NUBASE2020)
    age = np.linspace(0, 40000, 400)
    frac = 2.0**(-age / T)
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.plot(age / 1000, frac, color=BLUE, lw=2)
    a = T * math.log2(4)
    ax.plot([0, a / 1000, a / 1000], [0.25, 0.25, 0], color=VERMILLION, ls="--", lw=1)
    ax.text(a / 1000 + 0.6, 0.27, f"25 % left: age = {a:.0f} years", color=VERMILLION)
    ax.set_xlabel("Age of the sample (thousands of years)")
    ax.set_ylabel(r"Fraction of $^{14}$C left")
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 1.02)
    save(fig, "carbon14")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_decay_law()
    fig_euler()
    fig_monte_carlo()
    fig_chains()
    fig_binding_energy()
    fig_chart()
    fig_carbon()
