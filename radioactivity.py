"""
Radioactive decay of a single nuclide: the simple program.

    python3 radioactivity.py

The decay law dN/dt = -N / tau, with tau = T_1/2 / ln 2, is solved with Euler's method for a
nuclide chosen from nuclides.csv and compared with the exact solution N(t) = N0 exp(-t / tau).
"""
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_STEPS = 100000


def ask_number(prompt, cast, minimum, maximum=None):
    """Ask again until the answer is a valid number in the allowed range."""
    while True:
        answer = input(prompt).strip()
        try:
            value = cast(answer)
        except ValueError:
            print("Please enter a number.")
            continue
        if value < minimum or (maximum is not None and value > maximum):
            print("This value is out of range.")
        else:
            return value


def euler_decay(N0, half_life, total_time, steps):
    """Euler's method for dN/dt = -N / tau. Returns the times, N (Euler) and N (exact)."""
    tau = half_life / math.log(2)
    dt = total_time / steps
    t = dt * np.arange(steps + 1)
    N = np.zeros(steps + 1)
    N[0] = N0
    for i in range(steps):
        N[i + 1] = N[i] - dt * N[i] / tau
    return t, N, N0 * np.exp(-t / tau)


def main():
    df = pd.read_csv("nuclides.csv")
    print("Radioactive decay of a single nuclide")
    print("=" * 37)
    print("\nList of radioactive nuclides")
    print(df)

    while True:
        choice = ask_number("\nChoose the key of a nuclide (or -1 to quit): ", int, -1, len(df) - 1)
        if choice == -1:
            break
        name = df.loc[choice, "Nuclide Name"]
        half_life = float(df.loc[choice, "Half-life (years)"])
        print("=" * 30)
        print(f"You chose {name}. Its half-life is {half_life:.4g} years.")

        N0 = ask_number("Initial number of nuclei: ", int, 1)
        answer = input(f"The number of time steps is {DEFAULT_STEPS} by default. Change it? [yes/no] ")
        if answer.strip().lower() in ("yes", "y"):
            steps = ask_number("New number of steps: ", int, 1)
        else:
            steps = DEFAULT_STEPS
        total_time = ask_number(f"Duration of the simulation in years (try {5 * half_life:.3g}, "
                                "i.e. 5 half-lives): ", float, 1e-300)

        t, N, N_exact = euler_decay(N0, half_life, total_time, steps)
        print(f"After {total_time:.4g} years: N = {N[-1]:.6g} (Euler) and {N_exact[-1]:.6g} (exact).")

        plt.plot(t, N, linewidth=1.5, label="Euler's method")
        plt.plot(t, N_exact, "--", linewidth=1.5, label="Exact solution")
        plt.axhline(N0 / 2, color="gray", linewidth=0.8)
        plt.axvline(half_life, color="gray", linewidth=0.8)
        plt.xlabel("Time (years)")
        plt.ylabel("Number of nuclei")
        plt.title(f"Decay of {name}")
        plt.grid(visible=True)
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
