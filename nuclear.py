import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting

# Import the nuclide dictionary from the separate file
from nuclides_data import nuclides

# Ideal decay chain suggestions for specific parent nuclides.
# For each parent (key = parent's nuclide symbol), we suggest:
# - daughter: ideal daughter nuclide symbol (or name)
# - daughter_half_life: ideal half-life (years) for the daughter nuclide
# - suggested_N0: suggested initial number of nuclides
# - suggested_steps: suggested simulation steps
# - suggested_time_multiplier: multiplier for parent's half-life to define total simulation time
ideal_decay_chain_suggestions = {
    "238U": {"daughter": "234Th", "daughter_half_life": 0.066, "suggested_N0": 1_000_000, "suggested_steps": 5000, "suggested_time_multiplier": 5},
    "235U": {"daughter": "231Pa", "daughter_half_life": 0.0013, "suggested_N0": 1_000_000, "suggested_steps": 5000, "suggested_time_multiplier": 5},
    # Add more suggestions as needed.
}

# =============================================================================
# Header display with nuclear physics context
# =============================================================================
def display_header():
    header = (
        " " * 40 + "----------------------------------------\n" +
        " " * 40 + "|        University Nuclear Physics    |\n" +
        " " * 40 + "|   Radioactivity, Decay & Decay Chains  |\n" +
        " " * 40 + "----------------------------------------\n"
    )
    print(header)
    print("This simulation demonstrates key nuclear physics concepts:")
    print("- Radioactive decay, half-life, and decay constant")
    print("- Activity (A = λN)")
    print("- Comparison of numerical (Euler) and analytical solutions")
    print("- Decay chains (parent -> daughter) with 3D visualization option\n")

# =============================================================================
# Data loading function (load nuclide data from dictionary)
# =============================================================================
def load_nuclides_dict():
    """
    Convert the nuclide dictionary from nuclides_data.py into a pandas DataFrame.
    The resulting DataFrame will have columns:
    'Nuclide Symbol', 'Nuclide Name', 'Half-life (years)'
    """
    df = pd.DataFrame.from_dict(nuclides, orient='index')
    df = df.reset_index().rename(columns={'index': 'Nuclide Symbol', 'name': 'Nuclide Name', 'half_life': 'Half-life (years)'})
    return df

def display_nuclides(df):
    """Display the available nuclides with an index for selection."""
    print("List of Radioactive Nuclides:")
    print(df[['Nuclide Symbol', 'Nuclide Name', 'Half-life (years)']])
    print("\n")

# =============================================================================
# User input function for nuclide selection (returns symbol, name, and half-life)
# =============================================================================
def get_nuclide_choice(df):
    """
    Prompt the user to choose a nuclide.
    Returns:
        index (int): index of the chosen nuclide,
        nuclide_symbol (str): symbol of the chosen nuclide,
        nuclide_name (str): name of the chosen nuclide,
        half_life (float) in years.
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
                print("You chose:", nuclide_name, f"({nuclide_symbol})")
                print(f"Half-life: {half_life:.3e} years")
                return choice, nuclide_symbol, nuclide_name, half_life
        except ValueError:
            print("Invalid input. Please enter a numeric value.")

# =============================================================================
# Function to get simulation parameters with default suggestions
# =============================================================================
def get_simulation_parameters(default_N0, default_L, default_total_time):
    """
    Ask user for simulation parameters with suggested default values.
    
    Parameters:
        default_N0: Suggested initial number of nuclides.
        default_L: Suggested number of simulation steps.
        default_total_time: Suggested total simulation time in years.
    
    Returns:
        N0: initial number of nuclides (int)
        L: number of simulation steps (int)
        total_time: total simulation time in years (float)
    """
    print(f"Suggested parameters: N0 = {default_N0}, Simulation Steps = {default_L}, Total Time = {default_total_time} years")
    N0_in = input(f"Enter the initial number of nuclides (N0) [default: {default_N0}]: ")
    N0 = int(N0_in) if N0_in.strip() != "" else default_N0
    
    L_in = input(f"Enter the number of simulation steps [default: {default_L}]: ")
    L = int(L_in) if L_in.strip() != "" else default_L
    
    total_time_in = input(f"Enter the total simulation time (in years) [default: {default_total_time}]: ")
    total_time = float(total_time_in) if total_time_in.strip() != "" else default_total_time
    
    return N0, L, total_time

# =============================================================================
# Simulation functions: Single Decay (Numerical and Analytical)
# =============================================================================
def simulate_decay(N0, half_life, L, total_time):
    """
    Simulate radioactive decay using Euler's method and compute the analytical solution.
    
    Returns:
        t: time array,
        N_numerical: numerical solution for N(t),
        N_analytical: analytical solution N(t) = N0 * exp(-t/tau),
        lambda_decay: decay constant (ln2 / half_life)
    """
    tau = half_life / math.log(2)
    dt = total_time / L
    t = np.linspace(0, total_time, L + 1)
    
    # Euler's method for numerical simulation
    N_numerical = np.zeros(L + 1)
    N_numerical[0] = N0
    for i in range(L):
        N_numerical[i+1] = N_numerical[i] - (dt * N_numerical[i]) / tau
    
    # Analytical solution: N(t) = N0 * exp(-t/tau)
    N_analytical = N0 * np.exp(-t / tau)
    
    lambda_decay = math.log(2) / half_life
    return t, N_numerical, N_analytical, lambda_decay

def calculate_activity(lambda_decay, N):
    """Calculate activity A(t) = λN(t)."""
    return lambda_decay * N

# =============================================================================
# Simulation function: Decay Chain (Two-Step Decay)
# =============================================================================
def simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time):
    """
    Simulate a two-step decay chain using Euler's method.
    Parent decays to daughter.
    
    Differential equations:
        dN1/dt = -λ1 * N1
        dN2/dt = λ1 * N1 - λ2 * N2
    
    Returns:
        t: time array,
        N1: parent nuclide population,
        N2: daughter nuclide population,
        lambda1: decay constant for parent,
        lambda2: decay constant for daughter.
    """
    lambda1 = math.log(2) / half_life_parent
    lambda2 = math.log(2) / half_life_daughter
    dt = total_time / L
    t = np.linspace(0, total_time, L + 1)
    
    N1 = np.zeros(L + 1)
    N2 = np.zeros(L + 1)
    N1[0] = N0
    N2[0] = 0  # assume no daughter at t=0
    
    for i in range(L):
        N1[i+1] = N1[i] - lambda1 * N1[i] * dt
        N2[i+1] = N2[i] + (lambda1 * N1[i] - lambda2 * N2[i]) * dt
        
    return t, N1, N2, lambda1, lambda2

# =============================================================================
# 3D Plotting Function for Decay Chain
# =============================================================================
def plot_3d_decay_chain(t, N_parent, N_daughter, parent_name, daughter_name):
    """
    Create a 3D plot of the decay chain using matplotlib's mplot3d.
    The x-axis is time, y-axis is the parent nuclide population, and
    z-axis is the daughter nuclide population.
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
# Plotting functions: 2D for Single Decay and Decay Chain
# =============================================================================
def plot_results(t, N_num, N_ana, lambda_decay, nuclide_name):
    """Plot number of nuclides and activity for single decay simulation."""
    A_num = calculate_activity(lambda_decay, N_num)
    A_ana = calculate_activity(lambda_decay, N_ana)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    axs[0].plot(t, N_num, label='Numerical (Euler)', linewidth=1.5)
    axs[0].plot(t, N_ana, '--', label='Analytical', linewidth=1.5)
    axs[0].set_xlabel('Time (years)')
    axs[0].set_ylabel('Number of Nuclides')
    axs[0].set_title(f'Decay of {nuclide_name}')
    axs[0].grid(True)
    axs[0].legend()
    
    axs[1].plot(t, calculate_activity(lambda_decay, N_num), label='Numerical Activity', linewidth=1.5)
    axs[1].plot(t, calculate_activity(lambda_decay, N_ana), '--', label='Analytical Activity', linewidth=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title('Radioactive Activity Over Time')
    axs[1].grid(True)
    axs[1].legend()
    
    plt.tight_layout()
    plt.show()

def plot_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name):
    """Plot the decay chain populations and their activities (2D plots)."""
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot populations
    axs[0].plot(t, N_parent, label=f'{parent_name} (Parent)', linewidth=1.5)
    axs[0].plot(t, N_daughter, label=f'{daughter_name} (Daughter)', linewidth=1.5)
    axs[0].set_xlabel('Time (years)')
    axs[0].set_ylabel('Number of Nuclides')
    axs[0].set_title('Decay Chain Population')
    axs[0].grid(True)
    axs[0].legend()
    
    # Plot activities
    axs[1].plot(t, calculate_activity(lambda1, N_parent), label=f'{parent_name} Activity', linewidth=1.5)
    axs[1].plot(t, calculate_activity(lambda2, N_daughter), label=f'{daughter_name} Activity', linewidth=1.5)
    axs[1].set_xlabel('Time (years)')
    axs[1].set_ylabel('Activity (decays/year)')
    axs[1].set_title('Decay Chain Activity')
    axs[1].grid(True)
    axs[1].legend()
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# (Optional) Export simulation data for advanced visualization in ParaView
# =============================================================================
def export_simulation_data(filename, t, *arrays):
    """
    Export simulation data to a CSV file.
    Each additional array will be saved as a separate column.
    This file can then be imported into ParaView for 3D visualization.
    """
    data = {"Time (years)": t}
    for i, arr in enumerate(arrays, start=1):
        data[f"Data_{i}"] = arr
    df_export = pd.DataFrame(data)
    df_export.to_csv(filename, index=False)
    print(f"Simulation data exported to {filename}")

# =============================================================================
# Main function to tie everything together
# =============================================================================
def main():
    display_header()
    # Load nuclide data from the dictionary (via our DataFrame conversion)
    df = load_nuclides_dict()
    display_nuclides(df)
    
    # Display ideal parameter suggestions to the user
    print("Ideal parameter suggestions:")
    print(" - Single Decay: N0 = 1e6, Simulation Steps = 5000, Total Time = 5 × (half-life)")
    print(" - Decay Chain: N0 = 1e6, Simulation Steps = 5000, Total Time = 5 × (parent's half-life)")
    print("   Additionally, if the decay chain plots appear jagged, consider increasing the simulation steps.")
    print("\nYou may press Enter to accept these suggestions or enter your own values.\n")
    
    simulation_type = input("Choose simulation type:\n1 - Single Nuclide Decay\n2 - Decay Chain (Parent -> Daughter)\nEnter 1 or 2: ").strip()
    
    if simulation_type == '1':
        # Single decay simulation
        choice, nuclide_symbol, nuclide_name, half_life = get_nuclide_choice(df)
        if choice is None:
            print("Exiting simulation.")
            return
        
        # Use suggested ideal parameters for single decay:
        default_N0 = 1_000_000
        default_L = 5000
        default_total_time = 0.01 * half_life
        
        N0, L, total_time = get_simulation_parameters(default_N0, default_L, default_total_time)
        t, N_num, N_ana, lambda_decay = simulate_decay(N0, half_life, L, total_time)
        plot_results(t, N_num, N_ana, lambda_decay, nuclide_name)
        
        # (Optional) Export data for external visualization
        export_choice = input("Export simulation data for external 3D visualization? [yes/no]: ").strip().lower()
        if export_choice == 'yes':
            export_simulation_data("single_decay_simulation.csv", t, N_num, N_ana)
    
    elif simulation_type == '2':
        # Decay chain simulation
        print("\n--- Parent Nuclide Selection ---")
        choice, parent_symbol, parent_name, half_life_parent = get_nuclide_choice(df)
        if choice is None:
            print("Exiting simulation.")
            return
        
        # Check for ideal daughter suggestion for this parent
        if parent_symbol in ideal_decay_chain_suggestions:
            suggestion = ideal_decay_chain_suggestions[parent_symbol]
            print(f"\nIdeal daughter suggestion for {parent_name} ({parent_symbol}):")
            print(f"  Daughter Nuclide: {suggestion['daughter']}")
            print(f"  Suggested Daughter Half-life: {suggestion['daughter_half_life']:.3e} years")
            print(f"  Suggested Parameters: N0 = {suggestion['suggested_N0']}, Steps = {suggestion['suggested_steps']}, Total Time = {suggestion['suggested_time_multiplier']} × (parent's half-life)")
            daughter_default = suggestion["daughter"]
            daughter_half_life_default = suggestion["daughter_half_life"]
            default_N0 = suggestion["suggested_N0"]
            default_L = suggestion["suggested_steps"]
            default_total_time = suggestion["suggested_time_multiplier"] * half_life_parent
        else:
            print("\nNo ideal daughter suggestion for this parent. Please enter the daughter details manually.")
            daughter_default = None
            daughter_half_life_default = None
            default_N0 = 1_000_000
            default_L = 5000
            default_total_time = 5 * half_life_parent
        
        # Prompt for daughter nuclide details:
        daughter_name_in = input(f"Enter the name of the Daughter nuclide [default: {daughter_default if daughter_default is not None else 'custom'}]: ")
        daughter_name = daughter_name_in.strip() if daughter_name_in.strip() != "" else (daughter_default if daughter_default is not None else input("Enter the Daughter nuclide name: "))
        
        if daughter_half_life_default is not None:
            daughter_half_life_in = input(f"Enter the half-life (in years) for {daughter_name} [default: {daughter_half_life_default}]: ")
            half_life_daughter = float(daughter_half_life_in) if daughter_half_life_in.strip() != "" else daughter_half_life_default
        else:
            while True:
                try:
                    half_life_daughter = float(input(f"Enter the half-life (in years) for {daughter_name}: "))
                    break
                except ValueError:
                    print("Invalid input. Please enter a numeric value for the half-life.")
        
        # Use suggested ideal parameters for decay chain:
        N0, L, total_time = get_simulation_parameters(default_N0, default_L, default_total_time)
        t, N_parent, N_daughter, lambda1, lambda2 = simulate_decay_chain(N0, half_life_parent, half_life_daughter, L, total_time)
        plot_decay_chain(t, N_parent, N_daughter, lambda1, lambda2, parent_name, daughter_name)
        
        # 3D Plot for the decay chain:
        plot_3d_decay_chain(t, N_parent, N_daughter, parent_name, daughter_name)
        
        # (Optional) Export data for external visualization
        export_choice = input("Export simulation data for external 3D visualization? [yes/no]: ").strip().lower()
        if export_choice == 'yes':
            export_simulation_data("decay_chain_simulation.csv", t, N_parent, N_daughter)
    
    else:
        print("Invalid simulation type selected. Exiting.")
        return
    
    cont = input("Do you want to run another simulation? [yes/no]: ").strip().lower()
    if cont == 'yes':
        main()
    else:
        print("Exiting simulation.")

if __name__ == "__main__":
    main()

