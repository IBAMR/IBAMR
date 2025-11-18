#!/usr/bin/env python3
"""
Analysis script for IBAMR simulation output
Extracts and plots key quantities from log files
"""

import sys
import re
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def parse_log_file(log_file):
    """Parse IBAMR log file and extract timestep information"""

    data = {
        'iteration': [],
        'time': [],
        'dt': [],
        'cfl': []
    }

    with open(log_file, 'r') as f:
        for line in f:
            # Extract iteration and time
            if 'At beginning of timestep #' in line:
                match = re.search(r'timestep # (\d+)', line)
                if match:
                    data['iteration'].append(int(match.group(1)))

            if 'Simulation time is' in line:
                match = re.search(r'Simulation time is ([\d.eE+-]+)', line)
                if match:
                    data['time'].append(float(match.group(1)))

            # Extract dt
            if 'dt =' in line:
                match = re.search(r'dt = ([\d.eE+-]+)', line)
                if match:
                    data['dt'].append(float(match.group(1)))

            # Extract CFL
            if 'CFL number' in line:
                match = re.search(r'CFL number = ([\d.eE+-]+)', line)
                if match:
                    data['cfl'].append(float(match.group(1)))

    return data

def plot_timestep_data(data, output_dir='./'):
    """Create plots of timestep data"""

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Plot timestep size
    if data['dt']:
        axes[0].plot(data['iteration'][:len(data['dt'])], data['dt'], 'b-', linewidth=2)
        axes[0].set_xlabel('Iteration')
        axes[0].set_ylabel('Timestep Size (dt)')
        axes[0].set_title('Timestep Size vs Iteration')
        axes[0].grid(True, alpha=0.3)

    # Plot CFL number
    if data['cfl']:
        axes[1].plot(data['iteration'][:len(data['cfl'])], data['cfl'], 'r-', linewidth=2)
        axes[1].set_xlabel('Iteration')
        axes[1].set_ylabel('CFL Number')
        axes[1].set_title('CFL Number vs Iteration')
        axes[1].grid(True, alpha=0.3)
        axes[1].axhline(y=0.3, color='k', linestyle='--', label='Typical CFL limit')
        axes[1].legend()

    plt.tight_layout()
    output_file = Path(output_dir) / 'timestep_analysis.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved plot: {output_file}")
    plt.close()

def print_statistics(data):
    """Print summary statistics"""

    print("\n" + "="*50)
    print("SIMULATION STATISTICS")
    print("="*50)

    if data['iteration']:
        print(f"Total iterations: {max(data['iteration'])}")

    if data['time']:
        print(f"Start time: {min(data['time']):.6e}")
        print(f"End time: {max(data['time']):.6e}")
        print(f"Total simulation time: {max(data['time']) - min(data['time']):.6e}")

    if data['dt']:
        dt_array = np.array(data['dt'])
        print(f"\nTimestep size (dt):")
        print(f"  Min: {np.min(dt_array):.6e}")
        print(f"  Max: {np.max(dt_array):.6e}")
        print(f"  Mean: {np.mean(dt_array):.6e}")
        print(f"  Std: {np.std(dt_array):.6e}")

    if data['cfl']:
        cfl_array = np.array(data['cfl'])
        print(f"\nCFL number:")
        print(f"  Min: {np.min(cfl_array):.6e}")
        print(f"  Max: {np.max(cfl_array):.6e}")
        print(f"  Mean: {np.mean(cfl_array):.6e}")
        print(f"  Std: {np.std(cfl_array):.6e}")

    print("="*50 + "\n")

def main():
    """Main function"""

    if len(sys.argv) < 2:
        print("Usage: python analyze.py <log_file>")
        print("Example: python analyze.py IB2d.log")
        sys.exit(1)

    log_file = sys.argv[1]

    if not Path(log_file).exists():
        print(f"Error: Log file not found: {log_file}")
        sys.exit(1)

    print(f"Analyzing log file: {log_file}")

    # Parse log file
    data = parse_log_file(log_file)

    # Print statistics
    print_statistics(data)

    # Create plots
    if data['dt'] or data['cfl']:
        output_dir = Path(log_file).parent
        plot_timestep_data(data, output_dir)
    else:
        print("Warning: No data found to plot")

if __name__ == "__main__":
    main()
