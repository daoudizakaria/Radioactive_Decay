"""
Radioactive decay, decay chains and Monte Carlo decay: the complete interactive program.

    python3 nuclear.py

Four simulations are offered: (1) decay of a single nuclide and its activity, solved with
Euler's method and compared with the exact solution; (2) a parent-daughter chain; (3) a chain
with branching; (4) Monte Carlo decay, in which every nucleus decays at random. Decay chains
use the exact (Bateman) solution. The results can be exported to CSV files.
"""
import csv
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from nuclides_data import nuclides

# =============================================================================
# Header
# =============================================================================
def display_header():
    print("Radioactive decay and decay chains")
    print("=" * 34)
    print("The simulations illustrate:")
    print("- Radioactive decay, half-life, and decay constant")
    print("- Activity (A = λN)")
    print("- Comparison of numerical (Euler) and analytical solutions")
    print("- Decay chains (parent -> daughter) and branching ratios")
    print("- The random nature of radioactive decay (Monte Carlo)\n")

# =============================================================================
# Nuclide data
# =============================================================================
def load_nuclides_dict():
    """Nuclides of nuclides_data.py as a DataFrame (symbol, name, half-life in years)."""
    df = pd.DataFrame.from_dict(nuclides, orient='index')
    df = df.reset_index().rename(
        columns={
            'index': 'Nuclide Symbol',
            'name': 'Nuclide Name',
            'half_life': 'Half-life (years)'
        }
    )
    return df

def display_nuclides(df):
    """Display the available nuclides with an index for selection."""
    print("List of Radioactive Nuclides:")
    print(df[['Nuclide Symbol', 'Nuclide Name', 'Half-life (years)']])
    print("\n")

def display_decay_chain_candidates():
    """List the nuclides for which a daughter is suggested for the parent-daughter chain."""
    print("Parent nuclides with a suggested daughter:")
    found = False
    for symbol, data in nuclides.items():
        if data.get("daughter") is not None:
            print(f" - {data['name']} ({symbol}) -> {data['daughter']} "
                  f"(daughter half-life: {data['daughter_half_life']:.3e} years)")
            found = True
    if not found:
        print("   None available.")
    print("\n")

# =============================================================================
# Input helpers
# =============================================================================
def ask_yes_no(prompt):
    """Return True for yes/y (any capitalisation), False for anything else."""
    return input(prompt).strip().lower() in ("yes", "y")

def ask_number(prompt, default=None, cast=float, minimum=None, maximum=None):
    """
    Ask for a number until a valid one is given.
    An empty answer returns the default (if there is one).
    """
    while True:
        answer = input(prompt).strip()
        if answer == "" and default is not None:
            return default
        try:
            value = cast(float(answer)) if cast is int else cast(answer)
        except ValueError:
            print("Invalid input. Please enter a numeric value.")
            continue
        if minimum is not None and value < minimum:
            print(f"The value must be at least {minimum}.")
        elif maximum is not None and value > maximum:
            print(f"The value must be at most {maximum}.")
        else:
            return value

def get_nuclide_choice(df):
    """
    Prompt the user to choose a nuclide from the DataFrame.
    Returns: index, nuclide_symbol, nuclide_name, half_life
    """
    choice = ask_number("Choose a key (row index) of a nuclide (or type -1 to exit): ",
                        cast=int, minimum=-1, maximum=len(df) - 1)
    if choice == -1:
        return None, None, None, None
    nuclide_symbol = df.loc[choice, "Nuclide Symbol"]
    nuclide_name = df.loc[choice, "Nuclide Name"]
    half_life = float(df.loc[choice, "Half-life (years)"])
    print("\n" + "=" * 30)
    print(f"You chose: {nuclide_name} ({nuclide_symbol})")
    print(f"Half-life: {half_life:.3e} years")
    return choice, nuclide_symbol, nuclide_name, half_life

def get_simulation_parameters(default_N0, default_L, default_total_time):
    """
    Ask user for simulation parameters with suggested default values.
    """
    print(f"Suggested parameters: N0 = {default_N0}, Steps = {default_L}, "
          f"Total Time = {default_total_time:.4g} years")
    N0 = ask_number(f"Enter the initial number of nuclides (N0) [default: {default_N0}]: ",
                    default=default_N0, cast=int, minimum=1)
    L = ask_number(f"Enter the number of simulation steps [default: {default_L}]: ",
                   default=default_L, cast=int, minimum=1)
    total_time = ask_number(f"Enter the total simulation time (in years) [default: {default_total_time:.4g}]: ",
                            default=default_total_time, cast=float, minimum=1e-300)
    return N0, L, total_time

# =============================================================================
# Single Decay
# =============================================================================
def simulate_decay(N0, half_life, L, total_time):
    """
    Single nuclide decay, dN/dt = -λN.
    Solved with Euler's method and compared with the exact solution N0 exp(-λt).
    Returns: t, N_numerical, N_analytical, lambda_decay
    """
    lambda_decay = math.log(2) / half_life
    dt = total_time / L
    if lambda_decay * dt > 0.1:
        print(f"Warning: λ·dt = {lambda_decay * dt:.3g} is not small, so Euler's method "
              "will be inaccurate. Use more steps to improve it.")
    t = np.linspace(0, total_time, L + 1)
    N_numerical = np.zeros(L + 1)
    N_numerical[0] = N0
    for i in range(L):
        N_numerical[i+1] = N_numerical[i] - lambda_decay * N_numerical[i] * dt
    N_analytical = N0 * np.exp(-lambda_decay * t)
    return t, N_numerical, N_analytical, lambda_decay

def calculate_activity(lambda_decay, N):
    return lambda_decay * N

# =============================================================================
# Monte Carlo decay: every nucleus decays at random
# =============================================================================
def simulate_decay_monte_carlo(N0, half_life, L, total_time, seed=None):
    """
    Radioactive decay as a random process. During a step dt, each remaining
    nucleus decays with probability p = 1 - exp(-λ dt), independently of the
    others, so the number of decays in the step follows a binomial law.
    Returns: t, N (random), N_expected, sigma, lambda_decay
    where N_expected = N0 exp(-λt) is the average and
    sigma = sqrt(N0 q (1-q)), q = exp(-λt), the size of the random fluctuations.
    """
    rng = np.random.default_rng(seed)
    lambda_decay = math.log(2) / half_life
    dt = total_time / L
    p = -math.expm1(-lambda_decay * dt)          # 1 - exp(-λ dt), accurate for small λ dt
    t = np.linspace(0, total_time, L + 1)
    N = np.zeros(L + 1, dtype=np.int64)
    N[0] = N0
    for i in range(L):
        N[i+1] = N[i] - rng.binomial(N[i], p)
    q = np.exp(-lambda_decay * t)
    N_expected = N0 * q
    sigma = np.sqrt(N0 * q * (1 - q))
    return t, N, N_expected, sigma, lambda_decay

def plot_monte_carlo(t, N, N_expected, sigma, nuclide_name, n_runs=1):
    """Random decay compared with the exponential law and its ±1σ band."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.fill_between(t, N_expected - sigma, N_expected + sigma, color='gray', alpha=0.3,
                    label='Exponential law ± 1σ')
    ax.plot(t, N_expected, 'k--', lw=1.5, label='Exponential law N0 exp(-λt)')
    for k, curve in enumerate(np.atleast_2d(N)):
        ax.step(t, curve, where='post', lw=1.2, label='Random decays' if k == 0 else None)
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Number of Nuclides')
    ax.set_title(f'Random decay of {nuclide_name} ({n_runs} simulation{"s" if n_runs > 1 else ""})')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_results(t, N_num, N_ana, lambda_decay, nuclide_name):
    """
    Single decay plot: population and activity.
    """
    A_num = calculate_activity(lambda_decay, N_num)
    A_ana = calculate_activity(lambda_decay, N_ana)
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].plot(t, N_num, label='Numerical (Euler)', lw=1.5)
    axs[0].plot(t, N_ana, '--', label='Analytical', lw=1.5)
    axs[0].set_xlabel('Time (years)')
    axs[0].set_ylabel('Number of Nuclides')
    axs[0].set_title(f'Decay of {nuclide_name}')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(t, A_num, label='Numerical Activity', lw=1.5)
    axs[1].plot(t, A_ana, '--', label='Analytical Activity', lw=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title('Radioactive Activity Over Time')
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()
    plt.show()

# =============================================================================
# Decay chains: exact (Bateman) solutions
# =============================================================================
def bateman_daughter(N0, lambda_parent, lambda_daughter, t):
    """
    Number of daughter nuclei for Parent -> Daughter, starting with N0 parents
    and no daughters:
        N_D(t) = N0 λP / (λD - λP) * (exp(-λP t) - exp(-λD t))
    When λD = λP the limit is N0 λ t exp(-λt).

    The exact solution is used instead of Euler's method because the half-lives
    of parent and daughter are often very different (4.5 billion years for U-238,
    24 days for Th-234). Euler's method is unstable as soon as λD·dt > 2 and
    would give negative or huge numbers.
    """
    if math.isclose(lambda_parent, lambda_daughter, rel_tol=1e-12):
        return N0 * lambda_parent * t * np.exp(-lambda_parent * t)
    return (N0 * lambda_parent / (lambda_daughter - lambda_parent)
            * (np.exp(-lambda_parent * t) - np.exp(-lambda_daughter * t)))

# =============================================================================
# Simple Decay Chain (Parent -> Daughter)
# =============================================================================
def simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time):
    """
    Simple two-step chain: Parent -> Daughter (exact solution).
    Returns: t, N_parent, N_daughter, lambda1, lambda2
    """
    lambda1 = math.log(2) / half_life_parent
    lambda2 = math.log(2) / half_life_daughter
    t = np.linspace(0, total_time, L + 1)
    N_parent = N0 * np.exp(-lambda1 * t)
    N_daughter = bateman_daughter(N0, lambda1, lambda2, t)
    if lambda2 > 10 * lambda1:
        print(f"The daughter decays much faster than the parent: after a few daughter half-lives, "
              f"N_daughter/N_parent = λ1/(λ2-λ1) = {lambda1 / (lambda2 - lambda1):.3e} "
              "(secular equilibrium, equal activities).")
    return t, N_parent, N_daughter, lambda1, lambda2

def plot_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name):
    """Parent-daughter chain: populations and activities (logarithmic scales) and their ratio."""
    A_parent = calculate_activity(lambda1, N_parent)
    A_daughter = calculate_activity(lambda2, N_daughter)
    ratio_population = np.divide(N_daughter, N_parent, out=np.zeros_like(N_daughter), where=(N_parent != 0))

    fig, axs = plt.subplots(3, 1, figsize=(10, 15))

    # Subplot 1: Populations (log)
    axs[0].plot(t, N_parent, label=f'{parent_name} (Parent)', lw=1.5)
    axs[0].plot(t, N_daughter, label=f'{daughter_name} (Daughter)', lw=1.5)
    axs[0].set_xlabel('Time (years)')
    axs[0].set_ylabel('Number of Nuclides')
    axs[0].set_title('Decay Chain Population')
    axs[0].grid(True)
    axs[0].legend()
    axs[0].set_yscale('log')

    # Subplot 2: Activities (log)
    axs[1].plot(t, A_parent, label=f'{parent_name} Activity', lw=1.5)
    axs[1].plot(t, A_daughter, '--', label=f'{daughter_name} Activity', lw=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title('Decay Chain Activity')
    axs[1].grid(True)
    axs[1].legend()
    axs[1].set_yscale('log')

    # Subplot 3: Population Ratio
    axs[2].plot(t, ratio_population, label='Daughter/Parent Ratio', lw=1.5, color='purple')
    axs[2].set_xlabel('Time (years)')
    axs[2].set_ylabel('Population Ratio')
    axs[2].set_title('Daughter-to-Parent Population Ratio')
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.show()

def plot_3d_decay_chain(t, N_parent, N_daughter, parent_name, daughter_name):
    """Parent-daughter chain as a curve in the (t, N_parent, N_daughter) space."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(t, N_parent, N_daughter, label=f"{parent_name} -> {daughter_name}")
    ax.set_xlabel('Time (years)')
    ax.set_ylabel(f'{parent_name} Population')
    ax.set_zlabel(f'{daughter_name} Population')
    ax.set_title('3D Decay Chain Simulation')
    ax.legend()
    plt.show()

# =============================================================================
# Complex Decay Chain with Branching
# =============================================================================
def simulate_complex_decay_chain(N0, half_life_parent, half_life_A, half_life_B, BR_A, BR_B, L, total_time):
    """
    Complex decay chain with branching (exact solution):
      Parent --> Daughter A (branching ratio BR_A)
             --> Daughter B (branching ratio BR_B)
    Each branch is a Bateman solution multiplied by its branching ratio.
    """
    lambda_P = math.log(2) / half_life_parent
    lambda_A = math.log(2) / half_life_A
    lambda_B = math.log(2) / half_life_B
    t = np.linspace(0, total_time, L + 1)

    N_parent = N0 * np.exp(-lambda_P * t)
    N_A = BR_A * bateman_daughter(N0, lambda_P, lambda_A, t)
    N_B = BR_B * bateman_daughter(N0, lambda_P, lambda_B, t)

    return t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B

def plot_complex_decay_chain(t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B,
                             parent_name, daughter_A, daughter_B):
    """
    Plot the complex decay chain populations and activities with branching.
    """
    A_parent = calculate_activity(lambda_P, N_parent)
    A_A = calculate_activity(lambda_A, N_A)
    A_B = calculate_activity(lambda_B, N_B)

    ratio_A = np.divide(N_A, N_parent, out=np.zeros_like(N_A), where=(N_parent != 0))
    ratio_B = np.divide(N_B, N_parent, out=np.zeros_like(N_B), where=(N_parent != 0))

    fig, axs = plt.subplots(3, 1, figsize=(10, 15))

    # Subplot 1: Populations (log scale)
    axs[0].plot(t, N_parent, label=f'{parent_name} (Parent)', lw=1.5)
    axs[0].plot(t, N_A, label=f'{daughter_A} (Daughter A)', lw=1.5)
    axs[0].plot(t, N_B, label=f'{daughter_B} (Daughter B)', lw=1.5)
    axs[0].set_xlabel('Time (years)')
    axs[0].set_ylabel('Number of Nuclides')
    axs[0].set_title('Complex Decay Chain Population (Branching)')
    axs[0].grid(True)
    axs[0].legend()
    axs[0].set_yscale('log')

    # Subplot 2: Activities (log scale)
    axs[1].plot(t, A_parent, label=f'{parent_name} Activity', lw=1.5)
    axs[1].plot(t, A_A, label=f'{daughter_A} Activity', lw=1.5)
    axs[1].plot(t, A_B, label=f'{daughter_B} Activity', lw=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title('Complex Decay Chain Activity (Branching)')
    axs[1].grid(True)
    axs[1].legend()
    axs[1].set_yscale('log')

    # Subplot 3: Population Ratios
    axs[2].plot(t, ratio_A, label=f'{daughter_A} / {parent_name}', lw=1.5, color='green')
    axs[2].plot(t, ratio_B, label=f'{daughter_B} / {parent_name}', lw=1.5, color='red')
    axs[2].set_xlabel('Time (years)')
    axs[2].set_ylabel('Population Ratio')
    axs[2].set_title('Daughter-to-Parent Ratios (Branching)')
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.show()

# =============================================================================
# Export
# =============================================================================
def export_simulation_data(filename, t, *columns, headers=None):
    """Write the time array and any number of data columns to a CSV file."""
    if headers is None:
        headers = [f"Data_{k + 1}" for k in range(len(columns))]
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time (years)"] + list(headers))
        for row in zip(t, *columns):
            writer.writerow(row)
    print(f"Data exported to {filename}")

# =============================================================================
# Interactive Single Decay Simulation with ipywidgets (Jupyter only)
# =============================================================================
def interactive_decay_simulation(N0=1e6, half_life=1e4, steps=5000, time_multiplier=5):
    total_time = time_multiplier * half_life
    t, N_num, N_ana, lambda_decay = simulate_decay(int(N0), half_life, int(steps), total_time)
    plot_results(t, N_num, N_ana, lambda_decay, "Interactive Nuclide")

def interactive_widget():
    """
    Sliders for the single decay simulation. Use it in a Jupyter notebook:
        from nuclear import interactive_widget
        interactive_widget()
    Requires the ipywidgets package.
    """
    import ipywidgets as widgets
    from IPython.display import display

    interactive_sim = widgets.interactive(
        interactive_decay_simulation,
        N0=widgets.IntSlider(value=1000000, min=100000, max=10000000, step=100000, description='N0:'),
        half_life=widgets.FloatLogSlider(value=1e4, base=10, min=3, max=7, step=0.1, description='Half-life:'),
        steps=widgets.IntSlider(value=5000, min=1000, max=20000, step=500, description='Steps:'),
        time_multiplier=widgets.FloatSlider(value=5, min=1, max=10, step=0.5, description='Time mult.:')
    )
    display(interactive_sim)

# =============================================================================
# Interactive program
# =============================================================================
def run_simulation(df):
    """Run one simulation. Returns False if the user chose to exit."""
    print("Suggested parameters:")
    print(" - Single Decay: N0 = 1e6, Steps = 5000, Total Time = 5 × (half-life)")
    print(" - Simple Decay Chain: N0 = 1e6, Steps = 5000, Total Time = 5 × (parent's half-life)")
    print(" - Complex Decay Chain: same as above, but you can set branching ratios manually.")
    print(" - Monte Carlo Decay: N0 = 100 to 10000, Steps = 200, Total Time = 5 × (half-life)")
    print("\nPress Enter to accept suggestions or enter your own values.\n")

    print("Choose simulation type:\n1 - Single Nuclide Decay\n2 - Simple Decay Chain (Parent -> Daughter)\n3 - Complex Decay Chain (Branching)\n4 - Monte Carlo Decay (random decays)")
    simulation_type = input("Enter 1, 2, 3 or 4: ").strip()

    if simulation_type == '1':
        # Single decay
        choice, nuclide_symbol, nuclide_name, half_life = get_nuclide_choice(df)
        if choice is None:
            return False
        N0, L, total_time = get_simulation_parameters(1_000_000, 5000, 5 * half_life)
        t, N_num, N_ana, lambda_decay = simulate_decay(N0, half_life, L, total_time)
        plot_results(t, N_num, N_ana, lambda_decay, nuclide_name)
        filename, columns, headers = ("single_decay_simulation.csv", (N_num, N_ana),
                                      ("N numerical (Euler)", "N analytical"))

    elif simulation_type == '2':
        # Simple Decay Chain
        print("\n--- Parent Nuclide Selection ---")
        choice, parent_symbol, parent_name, half_life_parent = get_nuclide_choice(df)
        if choice is None:
            return False

        # a daughter and its half-life are suggested for some parents; otherwise they are entered
        suggestion = nuclides.get(parent_symbol, {})
        daughter_default = suggestion.get("daughter")
        half_life_daughter_default = suggestion.get("daughter_half_life")
        if daughter_default:
            print(f"\nSuggested daughter of {parent_name} ({parent_symbol}):")
            print(f"  Daughter Nuclide: {daughter_default}")
            print(f"  Daughter Half-life: {half_life_daughter_default:.3e} years")

        daughter_name_in = input(f"Enter the name of the Daughter nuclide [default: {daughter_default or 'custom'}]: ").strip()
        daughter_name = daughter_name_in or daughter_default or "CustomDaughter"

        # the suggested half-life only applies to the suggested daughter
        if daughter_name != daughter_default:
            half_life_daughter_default = None
        prompt = f"Enter the half-life (in years) for {daughter_name}"
        if half_life_daughter_default:
            prompt += f" [default: {half_life_daughter_default}]"
        half_life_daughter = ask_number(prompt + ": ", default=half_life_daughter_default,
                                        cast=float, minimum=1e-300)

        N0, L, total_time = get_simulation_parameters(1_000_000, 5000, 5 * half_life_parent)
        t, N_parent, N_daughter, lambda1, lambda2 = simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time)
        plot_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name)
        plot_3d_decay_chain(t, N_parent, N_daughter, parent_name, daughter_name)
        filename, columns, headers = ("simple_decay_chain.csv", (N_parent, N_daughter),
                                      (f"N {parent_name}", f"N {daughter_name}"))

    elif simulation_type == '3':
        # Complex Decay Chain (Branching)
        print("\n--- Parent Nuclide Selection for Complex Chain ---")
        choice, parent_symbol, parent_name, half_life_parent = get_nuclide_choice(df)
        if choice is None:
            return False

        # the two daughters, their half-lives and the branching ratio are entered by the user
        print("\nEnter Daughter A details:")
        daughterA = input("Name of Daughter A: ").strip() or "Daughter A"
        half_life_A = ask_number(f"Enter the half-life (in years) for {daughterA}: ", cast=float, minimum=1e-300)

        print("\nEnter Daughter B details:")
        daughterB = input("Name of Daughter B: ").strip() or "Daughter B"
        half_life_B = ask_number(f"Enter the half-life (in years) for {daughterB}: ", cast=float, minimum=1e-300)

        BR_A = ask_number("Enter the branching ratio for Daughter A (0 <= BR_A <= 1): ",
                          cast=float, minimum=0.0, maximum=1.0)
        BR_B = 1 - BR_A

        N0, L, total_time = get_simulation_parameters(1_000_000, 5000, 5 * half_life_parent)
        t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B = simulate_complex_decay_chain(
            N0, half_life_parent, half_life_A, half_life_B, BR_A, BR_B, L, total_time
        )
        plot_complex_decay_chain(t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B,
                                 parent_name, daughterA, daughterB)
        filename, columns, headers = ("complex_decay_chain.csv", (N_parent, N_A, N_B),
                                      (f"N {parent_name}", f"N {daughterA}", f"N {daughterB}"))

    elif simulation_type == '4':
        # Monte Carlo decay
        choice, nuclide_symbol, nuclide_name, half_life = get_nuclide_choice(df)
        if choice is None:
            return False
        N0, L, total_time = get_simulation_parameters(1000, 200, 5 * half_life)
        n_runs = ask_number("How many independent simulations? [default: 5]: ",
                            default=5, cast=int, minimum=1, maximum=100)
        runs = []
        for k in range(n_runs):
            t, N, N_expected, sigma, lambda_decay = simulate_decay_monte_carlo(N0, half_life, L, total_time)
            runs.append(N)
        plot_monte_carlo(t, np.array(runs), N_expected, sigma, nuclide_name, n_runs)
        filename, columns, headers = ("monte_carlo_decay.csv", (*runs, N_expected),
                                      tuple(f"N run {k + 1}" for k in range(n_runs)) + ("N expected",))

    else:
        print("Invalid simulation type selected.")
        return True

    if ask_yes_no("\nExport simulation data for external visualization? [yes/no]: "):
        export_simulation_data(filename, t, *columns, headers=headers)
    return True

def main():
    display_header()
    df = load_nuclides_dict()
    display_nuclides(df)
    display_decay_chain_candidates()

    while run_simulation(df):
        if not ask_yes_no("Do you want to run another simulation? [yes/no]: "):
            break
    print("Exiting simulation.")

if __name__ == "__main__":
    main()
