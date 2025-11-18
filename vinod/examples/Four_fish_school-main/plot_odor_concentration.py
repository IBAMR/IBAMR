#!/usr/bin/env python3
"""
Visualize odor concentration field with fish and fluid velocity.

This script creates publication-quality visualizations of:
- Odor concentration contours (normalized C*)
- Fluid velocity vectors/streamlines
- Fish body positions
- Vorticity fields

Based on the AIAA paper visualization methodology:
"How does vortex dynamics help undulating bodies spread odor"
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import h5py
import sys
import os
from pathlib import Path

# Configuration
VIZ_DIR = "viz_eel2d_Str"
OUTPUT_DIR = "odor_figures"
FISH_FILES = ["eel2d_1.vertex", "eel2d_2.vertex", "eel2d_3.vertex", "eel2d_4.vertex"]

# Physical parameters (from input2d)
Re = 5609.0
KAPPA = 1.0e-3
C_LOW = 0.0   # Background concentration Cl
C_HIGH = 10.0 # Source concentration Ch (from OdorSourceTerm)

def normalize_concentration(C):
    """Normalize concentration: C* = (C - Cl) / (Ch - Cl)"""
    return (C - C_LOW) / (C_HIGH - C_LOW)

def load_fish_vertices(vertex_file):
    """Load fish vertex positions from .vertex file"""
    vertices = []
    with open(vertex_file, 'r') as f:
        npts = int(f.readline().strip())
        for _ in range(npts):
            line = f.readline().strip().split()
            x, y = float(line[0]), float(line[1])
            vertices.append([x, y])
    return np.array(vertices)

def plot_odor_field(iteration, show_velocity=True, show_vorticity=False, save=True):
    """
    Create comprehensive odor visualization for a given iteration.

    Parameters:
    -----------
    iteration : int
        Iteration number to visualize
    show_velocity : bool
        Show velocity field as quiver plot
    show_vorticity : bool
        Show vorticity contours
    save : bool
        Save figure to file
    """
    # Create output directory
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    # Find visit dump file
    visit_file = f"{VIZ_DIR}/dumps.visit.{iteration:05d}.silo"

    if not os.path.exists(visit_file):
        print(f"Warning: File {visit_file} not found. Skipping iteration {iteration}")
        return

    try:
        # Load data from SILO file
        with h5py.File(visit_file, 'r') as f:
            # Load odor concentration (variable name may vary)
            if 'C' in f.keys():
                C_data = f['C'][:]
            elif 'concentration' in f.keys():
                C_data = f['concentration'][:]
            else:
                print(f"Available fields: {list(f.keys())}")
                print("Odor concentration field not found in output")
                return

            # Load velocity components
            if 'U' in f.keys():
                U_data = f['U'][:]
                V_data = f['V'][:]
            elif 'velocity_0' in f.keys():
                U_data = f['velocity_0'][:]
                V_data = f['velocity_1'][:]
            else:
                U_data = None
                V_data = None

            # Load vorticity
            if 'Omega' in f.keys():
                omega_data = f['Omega'][:]
            else:
                omega_data = None

            # Load grid coordinates
            x = f['x'][:]
            y = f['y'][:]

    except Exception as e:
        print(f"Error reading SILO file: {e}")
        print("Note: This script is designed for post-processing. Run simulation first.")
        return

    # Normalize concentration
    C_normalized = normalize_concentration(C_data)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot normalized odor concentration contours
    levels = np.linspace(0, 1, 21)
    contourf = ax.contourf(x, y, C_normalized, levels=levels, cmap='YlOrRd', alpha=0.8)
    contour = ax.contour(x, y, C_normalized, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
                         colors='black', linewidths=0.5, alpha=0.6)
    ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax, label=r'Normalized Concentration $C^*$')

    # Plot vorticity if requested
    if show_vorticity and omega_data is not None:
        vort_levels = np.linspace(-2, 2, 20)
        ax.contour(x, y, omega_data, levels=vort_levels, cmap='RdBu',
                   linewidths=0.3, alpha=0.3)

    # Plot velocity field if requested
    if show_velocity and U_data is not None and V_data is not None:
        # Subsample for clarity
        skip = max(1, len(x) // 30)
        ax.quiver(x[::skip, ::skip], y[::skip, ::skip],
                 U_data[::skip, ::skip], V_data[::skip, ::skip],
                 alpha=0.5, scale=10, width=0.003, color='blue')

    # Plot fish bodies
    for fish_file in FISH_FILES:
        if os.path.exists(fish_file):
            vertices = load_fish_vertices(fish_file)
            ax.plot(vertices[:, 0], vertices[:, 1], 'k-', linewidth=2, label='Fish')
            ax.fill(vertices[:, 0], vertices[:, 1], color='gray', alpha=0.7)

    # Mark odor source location
    source_x, source_y = -2.0, 0.0
    source_radius = 0.2
    circle = plt.Circle((source_x, source_y), source_radius,
                        color='red', fill=False, linewidth=2, linestyle='--',
                        label='Odor Source')
    ax.add_patch(circle)
    ax.plot(source_x, source_y, 'r*', markersize=15, label='Source Center')

    # Formatting
    ax.set_xlabel('x (fish lengths)', fontsize=12)
    ax.set_ylabel('y (fish lengths)', fontsize=12)
    ax.set_title(f'Odor Concentration Field - Iteration {iteration}\n' +
                 r'$Re = 5609$, $Sc = \nu/\kappa$, 4-Fish School',
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Add legend (avoid duplicates)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    plt.tight_layout()

    if save:
        output_file = f"{OUTPUT_DIR}/odor_field_iter_{iteration:05d}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")

    plt.show()

def create_animation_frames(iterations=None, max_frames=50):
    """
    Create multiple frames for animation.

    Parameters:
    -----------
    iterations : list or None
        List of iterations to process. If None, auto-detect from viz directory
    max_frames : int
        Maximum number of frames to generate
    """
    if iterations is None:
        # Auto-detect iterations from viz directory
        silo_files = sorted(Path(VIZ_DIR).glob("dumps.visit.*.silo"))
        iterations = [int(f.stem.split('.')[-1]) for f in silo_files[:max_frames]]

    print(f"Creating {len(iterations)} animation frames...")
    for i, iteration in enumerate(iterations):
        print(f"Processing frame {i+1}/{len(iterations)}: iteration {iteration}")
        plot_odor_field(iteration, show_velocity=True, show_vorticity=False, save=True)
        plt.close('all')  # Free memory

def analyze_odor_profiles(iteration, y_slice=0.0):
    """
    Analyze odor concentration profiles along a horizontal line.
    Similar to Figure 9 in the AIAA paper.

    Parameters:
    -----------
    iteration : int
        Iteration to analyze
    y_slice : float
        Y-coordinate of the slice
    """
    visit_file = f"{VIZ_DIR}/dumps.visit.{iteration:05d}.silo"

    if not os.path.exists(visit_file):
        print(f"File not found: {visit_file}")
        return

    with h5py.File(visit_file, 'r') as f:
        C_data = f['C'][:]
        x = f['x'][:]
        y = f['y'][:]

    # Find closest y index
    y_idx = np.argmin(np.abs(y[:, 0] - y_slice))

    # Extract profile
    C_profile = normalize_concentration(C_data[y_idx, :])
    x_profile = x[y_idx, :]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_profile, C_profile, 'b-', linewidth=2, label=f'$C^*$ at $y={y_slice}$')
    ax.set_xlabel('x (fish lengths)', fontsize=12)
    ax.set_ylabel(r'Normalized Concentration $C^*$', fontsize=12)
    ax.set_title(f'Odor Concentration Profile - Iteration {iteration}', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()

    output_file = f"{OUTPUT_DIR}/odor_profile_iter_{iteration:05d}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.show()

if __name__ == "__main__":
    print("=" * 70)
    print("Odor Concentration Visualization Tool")
    print("For 4-Fish School IBAMR Simulation")
    print("=" * 70)

    if len(sys.argv) > 1:
        # User specified iteration
        iteration = int(sys.argv[1])
        print(f"\nVisualizing iteration {iteration}...")
        plot_odor_field(iteration, show_velocity=True, show_vorticity=True)
        analyze_odor_profiles(iteration, y_slice=0.0)
    else:
        print("\nUsage:")
        print("  Single frame:  python plot_odor_concentration.py <iteration>")
        print("  Animation:     python plot_odor_concentration.py --animate")
        print("\nExample:")
        print("  python plot_odor_concentration.py 100")
        print("\nNote: Run IBAMR simulation first to generate viz_eel2d_Str/ directory")
