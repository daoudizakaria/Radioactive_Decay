import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
import ipywidgets as widgets
from IPython.display import display

# Import the enriched nuclide dictionary from nuclides_data.py
from nuclides_data import nuclides

# =============================================================================
# Header display with nuclear physics context
# =============================================================================
def display_header():
    header = (
        " " * 40 +   "-----------------------------------------\n"
        + " " * 40 + "|             Zakaria Daoudi            |\n"
        + " " * 40 + "|   Radioactivity, Decay & Decay Chains |\n"
        + " " * 40 + "-----------------------------------------\n"
    )
    print(header)
    print("This simulation demonstrates key nuclear physics concepts:")
    print("- Radioactive decay, half-life, and decay constant")
    print("- Activity (A = λN)")
    print("- Comparison of numerical (Euler) and analytical solutions")
    print("- Decay chains (parent -> daughter) and branching ratios\n")

# =============================================================================
# Data loading function (load nuclide data from dictionary)
# =============================================================================
def load_nuclides_dict():
    """
    Convert the nuclide dictionary from nuclides_data.py into a pandas DataFrame.
    """
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

def display_ideal_daughter_candidates():
    """Display simple decay chain candidates."""
    print("Ideal simple decay chain candidates (with ideal daughter suggestions):")
    found = False
    for symbol, data in nuclides.items():
        if data.get("ideal_daughter") is not None:
            print(f" - {data['name']} ({symbol}) -> Ideal daughter: {data['ideal_daughter']} "
                  f"(Suggested daughter half-life: {data['daughter_half_life']:.3e} years)")
            found = True
    if not found:
        print("   None available.")
    print("\n")

def display_complex_chain_candidates():
    """Display complex decay chain candidates (branching)."""
    print("Complex decay chain candidates (with branching):")
    found = False
    for symbol, data in nuclides.items():
        # If a nuclide is defined with branching info, e.g. "daughter_A", "daughter_B", "BR_A", "BR_B"
        # you can check them here. For now, we assume the dictionary might have these fields.
        if data.get("daughter_A") and data.get("daughter_B"):
            print(f" - {data['name']} ({symbol}) -> Daughter A: {data['daughter_A']} / Daughter B: {data['daughter_B']}")
            found = True
    if not found:
        print("   None available.")
    print("\n")

def get_nuclide_choice(df):
    """
    Prompt the user to choose a nuclide from the DataFrame.
    Returns: index, nuclide_symbol, nuclide_name, half_life
    """
    while True:
        try:
            choice = int(input("Choose a key (row index) of a nuclide (or type -1 to exit): "))
            if choice == -1:
                return None, None, None, None
            if choice < 0 or choice >= len(df):
                print("Invalid choice. Please select a valid index.")
            else:
                nuclide_symbol = df.loc[choice, "Nuclide Symbol"]
                nuclide_name = df.loc[choice, "Nuclide Name"]
                half_life = float(df.loc[choice, "Half-life (years)"])
                print("\n" + "=" * 30)
                print(f"You chose: {nuclide_name} ({nuclide_symbol})")
                print(f"Half-life: {half_life:.3e} years")
                return choice, nuclide_symbol, nuclide_name, half_life
        except ValueError:
            print("Invalid input. Please enter a numeric value.")

def get_simulation_parameters(default_N0, default_L, default_total_time):
    """
    Ask user for simulation parameters with suggested default values.
    """
    print(f"Suggested parameters: N0 = {default_N0}, Steps = {default_L}, Total Time = {default_total_time} years")
    N0_in = input(f"Enter the initial number of nuclides (N0) [default: {default_N0}]: ")
    N0 = int(N0_in) if N0_in.strip() != "" else default_N0
    
    L_in = input(f"Enter the number of simulation steps [default: {default_L}]: ")
    L = int(L_in) if L_in.strip() != "" else default_L
    
    total_time_in = input(f"Enter the total simulation time (in years) [default: {default_total_time}]: ")
    total_time = float(total_time_in) if total_time_in.strip() != "" else default_total_time
    
    return N0, L, total_time

# =============================================================================
# Single Decay
# =============================================================================
def simulate_decay(N0, half_life, L, total_time):
    """
    Single nuclide decay using Euler's method.
    Returns: t, N_numerical, N_analytical, lambda_decay
    """
    tau = half_life / math.log(2)
    dt = total_time / L
    t = np.linspace(0, total_time, L + 1)
    N_numerical = np.zeros(L + 1)
    N_numerical[0] = N0
    for i in range(L):
        N_numerical[i+1] = N_numerical[i] - (dt * N_numerical[i]) / tau
    N_analytical = N0 * np.exp(-t / tau)
    lambda_decay = math.log(2) / half_life
    return t, N_numerical, N_analytical, lambda_decay

def calculate_activity(lambda_decay, N):
    return lambda_decay * N

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
    
    A_title = 'Radioactive Activity Over Time'
    axs[1].plot(t, A_num, label='Numerical Activity', lw=1.5)
    axs[1].plot(t, A_ana, '--', label='Analytical Activity', lw=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title(A_title)
    axs[1].grid(True)
    axs[1].legend()
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# Simple Decay Chain (Parent -> Daughter)
# =============================================================================
def simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time):
    """
    Simple two-step chain: Parent -> Daughter
    Returns: t, N_parent, N_daughter, lambda1, lambda2
    """
    lambda1 = math.log(2) / half_life_parent
    lambda2 = math.log(2) / half_life_daughter
    dt = total_time / L
    t = np.linspace(0, total_time, L + 1)
    N_parent = np.zeros(L + 1)
    N_daughter = np.zeros(L + 1)
    N_parent[0] = N0
    for i in range(L):
        N_parent[i+1] = N_parent[i] - lambda1 * N_parent[i] * dt
        N_daughter[i+1] = N_daughter[i] + (lambda1 * N_parent[i] - lambda2 * N_daughter[i]) * dt
    return t, N_parent, N_daughter, lambda1, lambda2

def plot_enhanced_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name):
    """
    Enhanced 2D plotting for a simple two-step chain:
      - log scale populations
      - log scale activities
      - ratio
    """
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
    axs[1].plot(t, A_daughter, label=f'{daughter_name} Activity', lw=1.5)
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
    """
    3D plot for a simple two-step chain.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(t, N_parent, N_daughter, label=f"{parent_name} -> {daughter_name}", marker='o')
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
    Simulate a complex decay chain with branching:
      Parent --> Daughter A (branching ratio BR_A)
             --> Daughter B (branching ratio BR_B)
    """
    lambda_P = math.log(2) / half_life_parent
    lambda_A = math.log(2) / half_life_A
    lambda_B = math.log(2) / half_life_B
    dt = total_time / L
    t = np.linspace(0, total_time, L + 1)
    
    N_parent = np.zeros(L+1)
    N_A = np.zeros(L+1)
    N_B = np.zeros(L+1)
    
    N_parent[0] = N0
    for i in range(L):
        dN_parent = -lambda_P * N_parent[i] * dt
        N_parent[i+1] = N_parent[i] + dN_parent
        
        dN_A = (BR_A * lambda_P * N_parent[i] - lambda_A * N_A[i]) * dt
        N_A[i+1] = N_A[i] + dN_A
        
        dN_B = (BR_B * lambda_P * N_parent[i] - lambda_B * N_B[i]) * dt
        N_B[i+1] = N_B[i] + dN_B
        
    return t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B

def plot_complex_decay_chain(t, N_parent, N_A, N_B, parent_name, daughter_A, daughter_B):
    """
    Plot the complex decay chain populations and activities with branching.
    """
    # We do a naive approach for the lambdas if you don't pass them as parameters:
    # For a more accurate approach, pass lambda_P, lambda_A, lambda_B from simulate_complex_decay_chain.
    lambda_P = math.log(2) / (t[-1] / np.log(N_parent[0]/(N_parent[-1]+1e-10))) if N_parent[-1] != N_parent[0] else 0
    lambda_A = math.log(2) / (t[-1] / np.log((N_A[-1]+1e-10)/(1e-10))) if N_A[-1] > 1e-10 else 1
    lambda_B = math.log(2) / (t[-1] / np.log((N_B[-1]+1e-10)/(1e-10))) if N_B[-1] > 1e-10 else 1
    
    A_parent = calculate_activity(lambda_P, N_parent)
    A_A = calculate_activity(lambda_A, N_A)
    A_B = calculate_activity(lambda_B, N_B)
    
    ratio_A = np.divide(N_A, N_parent, out=np.zeros_like(N_A), where=(N_parent!=0))
    ratio_B = np.divide(N_B, N_parent, out=np.zeros_like(N_B), where=(N_parent!=0))
    
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
# Interactive Single Decay Simulation with ipywidgets
# =============================================================================
def interactive_decay_simulation(N0=1e6, half_life=1e4, steps=5000, time_multiplier=5):
    total_time = time_multiplier * half_life
    t, N_num, N_ana, lambda_decay = simulate_decay(int(N0), half_life, int(steps), total_time)
    plot_results(t, N_num, N_ana, lambda_decay, "Interactive Nuclide")

interactive_sim = widgets.interactive(
    interactive_decay_simulation,
    N0=widgets.IntSlider(value=1000000, min=100000, max=10000000, step=100000, description='N0:'),
    half_life=widgets.FloatLogSlider(value=1e4, base=10, min=3, max=7, step=0.1, description='Half-life:'),
    steps=widgets.IntSlider(value=5000, min=1000, max=20000, step=500, description='Steps:'),
    time_multiplier=widgets.FloatSlider(value=5, min=1, max=10, step=0.5, description='Time mult.:')
)
display(interactive_sim)

# =============================================================================
# Main function with 3 simulation types:
# 1) Single Decay
# 2) Simple Decay Chain (Parent -> Daughter)
# 3) Complex Decay Chain with Branching
# =============================================================================
def main():
    display_header()
    df = load_nuclides_dict()
    display_nuclides(df)
    display_ideal_daughter_candidates()
    display_complex_chain_candidates()
    
    print("Ideal parameter suggestions:")
    print(" - Single Decay: N0 = 1e6, Steps = 5000, Total Time = 5 × (half-life)")
    print(" - Simple Decay Chain: N0 = 1e6, Steps = 5000, Total Time = 5 × (parent's half-life)")
    print(" - Complex Decay Chain: same as above, but you can set branching ratios manually.")
    print("\nPress Enter to accept suggestions or enter your own values.\n")
    
    print("Choose simulation type:\n1 - Single Nuclide Decay\n2 - Simple Decay Chain (Parent -> Daughter)\n3 - Complex Decay Chain (Branching)")
    simulation_type = input("Enter 1, 2, or 3: ").strip()
    
    if simulation_type == '1':
        # Single decay
        choice, nuclide_symbol, nuclide_name, half_life = get_nuclide_choice(df)
        if choice is None:
            print("Exiting simulation.")
            return
        default_N0 = 1_000_000
        default_L = 5000
        default_total_time = 5 * half_life
        N0, L, total_time = get_simulation_parameters(default_N0, default_L, default_total_time)
        t, N_num, N_ana, lambda_decay = simulate_decay(N0, half_life, L, total_time)
        plot_results(t, N_num, N_ana, lambda_decay, nuclide_name)
    
    elif simulation_type == '2':
        # Simple Decay Chain
        print("\n--- Parent Nuclide Selection ---")
        choice, parent_symbol, parent_name, half_life_parent = get_nuclide_choice(df)
        if choice is None:
            print("Exiting simulation.")
            return
        
        # If there's an ideal daughter, use it, else user sets it manually
        if parent_symbol in nuclides and nuclides[parent_symbol].get("ideal_daughter"):
            suggestion = nuclides[parent_symbol]
            print(f"\nIdeal daughter suggestion for {parent_name} ({parent_symbol}):")
            print(f"  Daughter Nuclide: {suggestion['ideal_daughter']}")
            print(f"  Daughter Half-life: {suggestion['daughter_half_life']:.3e} years")
            daughter_default = suggestion["ideal_daughter"]
            half_life_daughter_default = suggestion["daughter_half_life"]
        else:
            daughter_default = None
            half_life_daughter_default = None
        
        daughter_name_in = input(f"Enter the name of the Daughter nuclide [default: {daughter_default if daughter_default else 'custom'}]: ")
        daughter_name = daughter_name_in.strip() if daughter_name_in.strip() != "" else (daughter_default if daughter_default else "CustomDaughter")
        
        if half_life_daughter_default:
            half_life_daughter_in = input(f"Enter the half-life (in years) for {daughter_name} [default: {half_life_daughter_default}]: ")
            half_life_daughter = float(half_life_daughter_in) if half_life_daughter_in.strip() != "" else half_life_daughter_default
        else:
            while True:
                try:
                    half_life_daughter = float(input(f"Enter the half-life (in years) for {daughter_name}: "))
                    break
                except ValueError:
                    print("Invalid input. Please enter a numeric value for the half-life.")
        
        default_N0 = 1_000_000
        default_L = 5000
        default_total_time = 5 * half_life_parent
        N0, L, total_time = get_simulation_parameters(default_N0, default_L, default_total_time)
        t, N_parent, N_daughter, lambda1, lambda2 = simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time)
        plot_enhanced_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name)
        plot_3d_decay_chain(t, N_parent, N_daughter, parent_name, daughter_name)
    
    elif simulation_type == '3':
        # Complex Decay Chain (Branching)
        print("\n--- Parent Nuclide Selection for Complex Chain ---")
        choice, parent_symbol, parent_name, half_life_parent = get_nuclide_choice(df)
        if choice is None:
            print("Exiting simulation.")
            return
        
        # For a complex chain, user must supply half-lives of Daughter A and B, and branching ratios.
        print("\nEnter Daughter A details:")
        daughterA = input("Name of Daughter A: ")
        half_life_A = float(input(f"Enter the half-life (in years) for {daughterA}: "))
        
        print("\nEnter Daughter B details:")
        daughterB = input("Name of Daughter B: ")
        half_life_B = float(input(f"Enter the half-life (in years) for {daughterB}: "))
        
        while True:
            try:
                BR_A_in = float(input("Enter the branching ratio for Daughter A (0 <= BR_A <= 1): "))
                if 0 <= BR_A_in <= 1:
                    BR_A = BR_A_in
                    break
                else:
                    print("Branching ratio must be between 0 and 1.")
            except ValueError:
                print("Invalid input. Please enter a numeric value for branching ratio.")
        BR_B = 1 - BR_A
        
        default_N0 = 1_000_000
        default_L = 5000
        default_total_time = 5 * half_life_parent
        N0, L, total_time = get_simulation_parameters(default_N0, default_L, default_total_time)
        
        t, N_parent, N_A, N_B, lambda_P, lambda_A, lambda_B = simulate_complex_decay_chain(
            N0, half_life_parent, half_life_A, half_life_B, BR_A, BR_B, L, total_time
        )
        plot_complex_decay_chain(t, N_parent, N_A, N_B, parent_name, daughterA, daughterB)
    
    else:
        print("Invalid simulation type selected. Exiting.")
        return

    # Optionally export results
    export_choice = input("\nExport simulation data for external visualization? [yes/no]: ").strip().lower()
    if export_choice == 'yes':
        if simulation_type == '1':
            # single decay
            export_simulation_data("single_decay_simulation.csv", t, N_num, N_ana)
        elif simulation_type == '2':
            # simple chain
            export_simulation_data("simple_decay_chain.csv", t, N_parent, N_daughter)
        elif simulation_type == '3':
            # complex chain
            export_simulation_data("complex_decay_chain.csv", t, N_parent, N_A, N_B)
    
    cont = input("Do you want to run another simulation? [yes/no]: ").strip().lower()
    if cont == 'yes':
        main()
    else:
        print("Exiting simulation.")

if __name__ == "__main__":
    main()


