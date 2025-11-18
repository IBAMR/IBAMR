#!/usr/bin/env python3
"""
Advanced odor plume analysis for 4-fish school simulation.

This script provides detailed analysis following the AIAA paper methodology:
1. Odor concentration iso-contours (C* levels)
2. Odor flux analysis
3. Mixing efficiency quantification
4. Vortex-odor interaction metrics

References:
- "How does vortex dynamics help undulating bodies spread odor.pdf"
- "Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf"
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
import h5py
import sys
from pathlib import Path

# Configuration
VIZ_DIR = "viz_eel2d_Str"
OUTPUT_DIR = "odor_analysis"
FISH_FILES = ["eel2d_1.vertex", "eel2d_2.vertex", "eel2d_3.vertex", "eel2d_4.vertex"]

# Physical parameters
Re = 5609.0
MU = 0.785 / Re
RHO = 1.0
KAPPA = 1.0e-3
SCHMIDT = MU / (RHO * KAPPA)

# Concentration parameters
C_LOW = 0.0
C_HIGH = 10.0

def load_fish_vertices(vertex_file):
    """Load fish vertex positions"""
    vertices = []
    with open(vertex_file, 'r') as f:
        npts = int(f.readline().strip())
        for _ in range(npts):
            line = f.readline().strip().split()
            x, y = float(line[0]), float(line[1])
            vertices.append([x, y])
    return np.array(vertices)

def normalize_concentration(C):
    """Normalize: C* = (C - Cl) / (Ch - Cl)"""
    return np.clip((C - C_LOW) / (C_HIGH - C_LOW), 0, 1)

def compute_odor_gradient(C, dx, dy):
    """
    Compute odor concentration gradient ∇C.

    Returns:
    --------
    grad_C_x, grad_C_y : ndarray
        Components of concentration gradient
    grad_C_mag : ndarray
        Magnitude of gradient |∇C|
    """
    grad_C_x = np.gradient(C, dx, axis=1)
    grad_C_y = np.gradient(C, dy, axis=0)
    grad_C_mag = np.sqrt(grad_C_x**2 + grad_C_y**2)

    return grad_C_x, grad_C_y, grad_C_mag

def compute_odor_flux(C, U, V):
    """
    Compute advective odor flux: F = u*C

    Returns:
    --------
    flux_x, flux_y : ndarray
        Components of odor flux vector
    flux_mag : ndarray
        Magnitude of flux
    """
    flux_x = U * C
    flux_y = V * C
    flux_mag = np.sqrt(flux_x**2 + flux_y**2)

    return flux_x, flux_y, flux_mag

def compute_mixing_efficiency(C, omega, dx, dy):
    """
    Compute mixing efficiency metric:
    η = (Variance of C) / (Mean vorticity²)

    Higher η indicates better odor spreading by vortices.
    """
    C_var = np.var(C)
    omega_rms = np.sqrt(np.mean(omega**2))

    if omega_rms > 0:
        eta = C_var / omega_rms**2
    else:
        eta = 0.0

    return eta

def plot_odor_isocontours(iteration):
    """
    Create publication-quality iso-contour plot.
    Similar to Figures 7-10 in AIAA paper.
    """
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    visit_file = f"{VIZ_DIR}/dumps.visit.{iteration:05d}.silo"

    if not os.path.exists(visit_file):
        print(f"File not found: {visit_file}")
        return

    # Load data
    with h5py.File(visit_file, 'r') as f:
        if 'C' not in f.keys():
            print("Odor concentration not found in output")
            return

        C_data = f['C'][:]
        x = f['x'][:]
        y = f['y'][:]

        if 'Omega' in f.keys():
            omega = f['Omega'][:]
        else:
            omega = None

    # Normalize
    C_star = normalize_concentration(C_data)

    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # ======== LEFT: Filled contours ========
    ax1 = axes[0]
    levels_fill = np.linspace(0, 1, 51)
    contourf = ax1.contourf(x, y, C_star, levels=levels_fill,
                            cmap='hot', extend='both')

    # Add specific iso-contours
    levels_iso = [0.1, 0.3, 0.5, 0.7, 0.9]
    contour = ax1.contour(x, y, C_star, levels=levels_iso,
                          colors='cyan', linewidths=1.5, linestyles='solid')
    ax1.clabel(contour, inline=True, fontsize=10, fmt='%.1f')

    # Plot fish
    for fish_file in FISH_FILES:
        if Path(fish_file).exists():
            vertices = load_fish_vertices(fish_file)
            ax1.fill(vertices[:, 0], vertices[:, 1], color='white',
                    edgecolor='black', linewidth=2, alpha=0.9)

    # Mark source
    source = Circle((-2.0, 0.0), 0.2, fill=False, edgecolor='lime',
                   linewidth=2, linestyle='--', label='Odor Source')
    ax1.add_patch(source)

    ax1.set_xlabel('x (L)', fontsize=12)
    ax1.set_ylabel('y (L)', fontsize=12)
    ax1.set_title(r'Normalized Concentration $C^* = (C-C_l)/(C_h-C_l)$',
                 fontsize=13, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.2)

    cbar1 = plt.colorbar(contourf, ax=ax1, label=r'$C^*$')

    # ======== RIGHT: Vorticity + odor overlay ========
    ax2 = axes[1]

    if omega is not None:
        # Vorticity background
        vort_levels = np.linspace(-3, 3, 31)
        contourf2 = ax2.contourf(x, y, omega, levels=vort_levels,
                                cmap='RdBu_r', alpha=0.7)

        # Overlay odor contours
        contour2 = ax2.contour(x, y, C_star, levels=levels_iso,
                              colors='yellow', linewidths=2, linestyles='solid')
        ax2.clabel(contour2, inline=True, fontsize=10, fmt='%.1f',
                  colors='yellow')

        cbar2 = plt.colorbar(contourf2, ax=ax2, label=r'Vorticity $\omega$ (1/s)')
    else:
        # Just show odor with gradient magnitude
        grad_x, grad_y, grad_mag = compute_odor_gradient(C_star,
                                                         x[0, 1] - x[0, 0],
                                                         y[1, 0] - y[0, 0])
        contourf2 = ax2.contourf(x, y, grad_mag, levels=20, cmap='viridis')
        cbar2 = plt.colorbar(contourf2, ax=ax2, label=r'$|\nabla C^*|$')

    # Plot fish
    for fish_file in FISH_FILES:
        if Path(fish_file).exists():
            vertices = load_fish_vertices(fish_file)
            ax2.fill(vertices[:, 0], vertices[:, 1], color='white',
                    edgecolor='black', linewidth=2, alpha=0.9)

    # Mark source
    source2 = Circle((-2.0, 0.0), 0.2, fill=False, edgecolor='lime',
                    linewidth=2, linestyle='--')
    ax2.add_patch(source2)

    ax2.set_xlabel('x (L)', fontsize=12)
    ax2.set_ylabel('y (L)', fontsize=12)
    ax2.set_title('Vorticity-Odor Interaction', fontsize=13, fontweight='bold')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.2)

    plt.suptitle(f'Odor Plume Analysis - Iteration {iteration}\n' +
                 r'$Re = 5609$, $Sc = \nu/\kappa$, 4-Fish Rectangular School',
                 fontsize=15, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_file = f"{OUTPUT_DIR}/odor_isocontours_iter_{iteration:05d}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.show()

def analyze_odor_statistics(iterations):
    """
    Compute statistical metrics over time:
    - Mean concentration
    - Variance (spreading)
    - Mixing efficiency
    - Odor coverage area
    """
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    times = []
    mean_conc = []
    variance_conc = []
    max_conc = []
    coverage_area = []
    mixing_eff = []

    for iteration in iterations:
        visit_file = f"{VIZ_DIR}/dumps.visit.{iteration:05d}.silo"

        if not Path(visit_file).exists():
            continue

        try:
            with h5py.File(visit_file, 'r') as f:
                C_data = f['C'][:]
                x = f['x'][:]
                y = f['y'][:]

                if 'time' in f.attrs:
                    t = f.attrs['time']
                else:
                    t = iteration * 0.0001  # Approximate from dt

                omega = f['Omega'][:] if 'Omega' in f.keys() else None

            C_star = normalize_concentration(C_data)

            # Compute statistics
            times.append(t)
            mean_conc.append(np.mean(C_star))
            variance_conc.append(np.var(C_star))
            max_conc.append(np.max(C_star))

            # Coverage area (fraction of domain with C* > 0.1)
            dx = x[0, 1] - x[0, 0]
            dy = y[1, 0] - y[0, 0]
            total_area = (x.max() - x.min()) * (y.max() - y.min())
            coverage_cells = np.sum(C_star > 0.1)
            coverage = coverage_cells * dx * dy / total_area
            coverage_area.append(coverage)

            # Mixing efficiency
            if omega is not None:
                eta = compute_mixing_efficiency(C_star, omega, dx, dy)
                mixing_eff.append(eta)

        except Exception as e:
            print(f"Error processing iteration {iteration}: {e}")
            continue

    # Plot time series
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Mean concentration
    axes[0, 0].plot(times, mean_conc, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time (s)', fontsize=11)
    axes[0, 0].set_ylabel(r'Mean $C^*$', fontsize=11)
    axes[0, 0].set_title('Mean Odor Concentration', fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)

    # Variance
    axes[0, 1].plot(times, variance_conc, 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Time (s)', fontsize=11)
    axes[0, 1].set_ylabel(r'Var($C^*$)', fontsize=11)
    axes[0, 1].set_title('Odor Spreading (Variance)', fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)

    # Coverage area
    axes[1, 0].plot(times, coverage_area, 'g-', linewidth=2)
    axes[1, 0].set_xlabel('Time (s)', fontsize=11)
    axes[1, 0].set_ylabel('Coverage Fraction', fontsize=11)
    axes[1, 0].set_title(r'Odor Coverage Area ($C^* > 0.1$)', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)

    # Mixing efficiency
    if mixing_eff:
        axes[1, 1].plot(times, mixing_eff, 'm-', linewidth=2)
        axes[1, 1].set_xlabel('Time (s)', fontsize=11)
        axes[1, 1].set_ylabel(r'$\eta$ (Mixing Efficiency)', fontsize=11)
        axes[1, 1].set_title('Vortex-Enhanced Mixing', fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)
    else:
        axes[1, 1].text(0.5, 0.5, 'No vorticity data available',
                       ha='center', va='center', transform=axes[1, 1].transAxes)

    plt.suptitle('Odor Transport Statistics - 4-Fish School',
                fontsize=15, fontweight='bold')
    plt.tight_layout()

    output_file = f"{OUTPUT_DIR}/odor_statistics_timeseries.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.show()

    # Save data to CSV
    csv_file = f"{OUTPUT_DIR}/odor_statistics.csv"
    with open(csv_file, 'w') as f:
        f.write("time,mean_C,variance_C,max_C,coverage_area,mixing_efficiency\n")
        for i in range(len(times)):
            eta = mixing_eff[i] if i < len(mixing_eff) else 0.0
            f.write(f"{times[i]},{mean_conc[i]},{variance_conc[i]}," +
                   f"{max_conc[i]},{coverage_area[i]},{eta}\n")
    print(f"Saved data: {csv_file}")

def main():
    print("=" * 70)
    print("Advanced Odor Plume Analysis")
    print("4-Fish School with Odor Transport")
    print("=" * 70)

    if len(sys.argv) > 1:
        if sys.argv[1] == "--stats":
            # Analyze statistics over multiple iterations
            print("\nComputing odor transport statistics...")
            iterations = range(0, 1000, 40)  # Every viz dump
            analyze_odor_statistics(iterations)
        else:
            # Single iteration analysis
            iteration = int(sys.argv[1])
            print(f"\nAnalyzing iteration {iteration}...")
            plot_odor_isocontours(iteration)
    else:
        print("\nUsage:")
        print("  Single frame: python analyze_odor_plumes.py <iteration>")
        print("  Statistics:   python analyze_odor_plumes.py --stats")
        print("\nExamples:")
        print("  python analyze_odor_plumes.py 200")
        print("  python analyze_odor_plumes.py --stats")
        print("\nNote: Requires IBAMR simulation output in viz_eel2d_Str/")

if __name__ == "__main__":
    main()
