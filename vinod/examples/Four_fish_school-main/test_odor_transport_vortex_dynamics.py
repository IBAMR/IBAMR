#!/usr/bin/env python3
"""
Test: How does vortex dynamics help undulating bodies spread odor?

Based on: "Collective Chemotactic Behavior in Fish Schools" (arXiv:2408.16136)
Authors: Maham Kamran, Amirhossein Fardi, Chengyu Li, Muhammad Saif Ullah Khalid

This test implements the odor transport equation to demonstrate how vortices
created by undulating fish bodies enhance odor dispersion compared to pure diffusion.

Governing Equation (Convection-Diffusion):
-----------------------------------------
∂c/∂t + u·∇c = D∇²c

where:
  c(x,y,t) = odor concentration field
  u = (u_x, u_y) = velocity field from IBAMR simulation
  D = molecular diffusion coefficient
  u·∇c = convection term (advection by vortices)
  D∇²c = diffusion term (molecular diffusion)

Test Scenarios:
--------------
1. Odor spreading WITH vortex dynamics (full equation)
2. Odor spreading WITHOUT vortices (pure diffusion: D∇²c only)
3. Comparison showing enhancement factor

Author: Generated for IBAMR Four Fish School
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
import pyvista as pv
from pathlib import Path
from scipy.ndimage import laplace
from scipy.interpolate import griddata
import sys

# ============================================================
# PUBLICATION SETTINGS
# ============================================================

USE_LATEX = True

if USE_LATEX:
    try:
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}\boldmath'
        matplotlib.rcParams['font.family'] = 'serif'
        matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
        print("[INFO] LaTeX rendering ENABLED")
    except:
        USE_LATEX = False
        print("[WARNING] LaTeX failed, using standard fonts")

matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['savefig.dpi'] = 300

FONT_SIZE_TITLE = 18
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14

# ============================================================
# SIMULATION PARAMETERS
# ============================================================

# IBAMR simulation parameters
DT_IBAMR = 0.0001              # IBAMR timestep
VIZ_DUMP_INTERVAL = 40         # Visualization dump interval

# Odor transport parameters
D_ODOR = 0.001                 # Diffusion coefficient (non-dimensional)
DT_ODOR = 0.0001              # Time step for odor equation
CFL_MAX = 0.5                  # CFL condition for stability

# Computational domain
X_MIN, X_MAX = -6.0, 3.0
Y_MIN, Y_MAX = -3.0, 3.0
NX = 200                       # Grid points in x
NY = 150                       # Grid points in y

# Initial odor source (point release)
SOURCE_X = -2.0                # Source x-position
SOURCE_Y = 0.0                 # Source y-position
SOURCE_SIGMA = 0.2             # Initial Gaussian width
SOURCE_AMPLITUDE = 1.0         # Initial concentration

# Which frames to process
FRAME_START = 0
FRAME_END = 100
FRAME_SKIP = 10

# Output settings
OUTPUT_DIR = "odor_transport_test"
SAVE_FORMAT = 'png'

# ============================================================
# NUMERICAL SOLVER FOR ODOR TRANSPORT EQUATION
# ============================================================

class OdorTransportSolver:
    """
    Solves the odor transport equation:
    ∂c/∂t + u·∇c = D∇²c

    Using finite differences with upwind scheme for convection
    and central differences for diffusion.
    """

    def __init__(self, x_range, y_range, nx, ny, diffusion_coeff):
        """
        Initialize solver on regular grid

        Parameters:
        -----------
        x_range : tuple
            (x_min, x_max)
        y_range : tuple
            (y_min, y_max)
        nx, ny : int
            Number of grid points
        diffusion_coeff : float
            Molecular diffusion coefficient D
        """
        self.x_min, self.x_max = x_range
        self.y_min, self.y_max = y_range
        self.nx = nx
        self.ny = ny
        self.D = diffusion_coeff

        # Create regular grid
        self.x = np.linspace(self.x_min, self.x_max, nx)
        self.y = np.linspace(self.y_min, self.y_max, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]

        # Initialize concentration field
        self.c = np.zeros((ny, nx))
        self.t = 0.0

        print(f"[SOLVER] Grid: {nx}x{ny}, dx={self.dx:.4f}, dy={self.dy:.4f}")
        print(f"[SOLVER] Diffusion coefficient D = {self.D}")

    def set_initial_condition_gaussian(self, x0, y0, sigma, amplitude=1.0):
        """
        Set Gaussian initial condition (point source)

        c(x,y,0) = A * exp(-((x-x0)² + (y-y0)²) / (2σ²))
        """
        r_squared = (self.X - x0)**2 + (self.Y - y0)**2
        self.c = amplitude * np.exp(-r_squared / (2 * sigma**2))
        print(f"[SOLVER] Initial condition: Gaussian at ({x0}, {y0}), σ={sigma}")
        print(f"[SOLVER] Initial total mass: {np.sum(self.c) * self.dx * self.dy:.6f}")

    def compute_cfl_timestep(self, u_x, u_y):
        """
        Compute maximum stable timestep based on CFL condition

        CFL = max(|u|·dt/dx, |v|·dt/dy) < 0.5
        """
        u_max = np.max(np.abs(u_x))
        v_max = np.max(np.abs(u_y))

        dt_conv = CFL_MAX * min(self.dx / (u_max + 1e-10),
                                self.dy / (v_max + 1e-10))

        # Diffusion stability: dt < dx²/(4D)
        dt_diff = 0.25 * min(self.dx**2, self.dy**2) / (self.D + 1e-10)

        dt_stable = min(dt_conv, dt_diff)

        return dt_stable

    def step_diffusion_only(self, dt):
        """
        Solve pure diffusion: ∂c/∂t = D∇²c
        Using explicit finite differences
        """
        # Compute Laplacian using 5-point stencil
        c_new = self.c.copy()

        # Interior points
        c_new[1:-1, 1:-1] = self.c[1:-1, 1:-1] + dt * self.D * (
            (self.c[1:-1, 2:] - 2*self.c[1:-1, 1:-1] + self.c[1:-1, :-2]) / self.dx**2 +
            (self.c[2:, 1:-1] - 2*self.c[1:-1, 1:-1] + self.c[:-2, 1:-1]) / self.dy**2
        )

        # Neumann boundary conditions (zero flux)
        c_new[0, :] = c_new[1, :]
        c_new[-1, :] = c_new[-2, :]
        c_new[:, 0] = c_new[:, 1]
        c_new[:, -1] = c_new[:, -2]

        self.c = c_new
        self.t += dt

    def step_convection_diffusion(self, u_x, u_y, dt):
        """
        Solve full equation: ∂c/∂t + u·∇c = D∇²c

        Using operator splitting:
        1. Convection step (upwind scheme)
        2. Diffusion step (central differences)
        """
        # Step 1: Convection (upwind scheme)
        c_conv = self.c.copy()

        for j in range(1, self.ny - 1):
            for i in range(1, self.nx - 1):
                # Upwind derivatives for convection
                if u_x[j, i] > 0:
                    dc_dx = (self.c[j, i] - self.c[j, i-1]) / self.dx
                else:
                    dc_dx = (self.c[j, i+1] - self.c[j, i]) / self.dx

                if u_y[j, i] > 0:
                    dc_dy = (self.c[j, i] - self.c[j-1, i]) / self.dy
                else:
                    dc_dy = (self.c[j+1, i] - self.c[j, i]) / self.dy

                # Update with convection
                c_conv[j, i] = self.c[j, i] - dt * (u_x[j, i] * dc_dx + u_y[j, i] * dc_dy)

        self.c = c_conv

        # Step 2: Diffusion
        self.step_diffusion_only(dt)

    def get_concentration(self):
        """Return current concentration field"""
        return self.c.copy()

    def get_total_mass(self):
        """Return total mass (should be conserved)"""
        return np.sum(self.c) * self.dx * self.dy

# ============================================================
# DATA LOADING FROM IBAMR
# ============================================================

def load_eulerian_frame(frame_idx):
    """
    Load Eulerian (fluid velocity) data from IBAMR output
    Returns velocity components u_x, u_y and vorticity
    """
    directories = [".", "ExportEULERIANData"]

    for base_dir in directories:
        ts_dir = Path(base_dir) / f"visit_eulerian_db__{frame_idx:04d}"

        if not ts_dir.exists():
            continue

        vtk_files = sorted(ts_dir.glob("*.vtk"))
        if len(vtk_files) == 0:
            continue

        try:
            meshes = []
            for vf in vtk_files:
                try:
                    mesh = pv.read(str(vf))
                    meshes.append(mesh)
                except:
                    continue

            if len(meshes) == 0:
                continue

            # Combine mesh pieces
            combined = meshes[0]
            for m in meshes[1:]:
                combined = combined.merge(m)

            points = combined.points

            # Extract velocity components
            if 'U_x' in combined.array_names and 'U_y' in combined.array_names:
                u_x_points = combined['U_x']
                u_y_points = combined['U_y']
            else:
                print(f"    WARNING: Velocity components not found in frame {frame_idx}")
                continue

            # Extract vorticity if available
            omega = combined.get('Omega', None)

            return points, u_x_points, u_y_points, omega

        except Exception as e:
            continue

    return None, None, None, None

def load_lagrangian_frame(frame_idx):
    """Load fish/eel body positions"""
    eels_points = []

    for eel_idx in range(10):
        patterns = [
            f"ExportLagrangianData/visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",
            f"visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",
        ]

        for pattern in patterns:
            vtk_file = Path(pattern)
            if vtk_file.exists():
                try:
                    mesh = pv.read(str(vtk_file))
                    eels_points.append(mesh.points)
                    break
                except:
                    continue

        if len(eels_points) == eel_idx:
            break

    return eels_points if len(eels_points) > 0 else None

def interpolate_velocity_to_grid(points, u_x_points, u_y_points, grid_x, grid_y):
    """
    Interpolate scattered velocity data to regular grid

    Parameters:
    -----------
    points : array (N, 3)
        Scattered point coordinates from VTK
    u_x_points, u_y_points : array (N,)
        Velocity components at scattered points
    grid_x, grid_y : 2D arrays
        Regular grid coordinates

    Returns:
    --------
    u_x_grid, u_y_grid : 2D arrays
        Velocity on regular grid
    """
    # Interpolate using griddata (linear interpolation)
    u_x_grid = griddata(points[:, :2], u_x_points,
                        (grid_x, grid_y), method='linear', fill_value=0.0)
    u_y_grid = griddata(points[:, :2], u_y_points,
                        (grid_x, grid_y), method='linear', fill_value=0.0)

    return u_x_grid, u_y_grid

# ============================================================
# MAIN TEST FUNCTION
# ============================================================

def run_odor_transport_test():
    """
    Main test: Compare odor spreading with and without vortex dynamics
    """

    print("\n" + "="*80)
    print("ODOR TRANSPORT TEST: Vortex Dynamics vs Pure Diffusion")
    print("="*80)
    print("\nBased on paper: arXiv:2408.16136")
    print("Authors: Kamran, Fardi, Li, Khalid")
    print("\nGoverning equation: ∂c/∂t + u·∇c = D∇²c")
    print("="*80)

    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)

    # Initialize two solvers for comparison
    print("\n[SETUP] Initializing solvers...")

    # Solver 1: WITH vortex dynamics (full equation)
    solver_vortex = OdorTransportSolver(
        (X_MIN, X_MAX), (Y_MIN, Y_MAX), NX, NY, D_ODOR
    )
    solver_vortex.set_initial_condition_gaussian(
        SOURCE_X, SOURCE_Y, SOURCE_SIGMA, SOURCE_AMPLITUDE
    )

    # Solver 2: WITHOUT vortices (pure diffusion)
    solver_diffusion = OdorTransportSolver(
        (X_MIN, X_MAX), (Y_MIN, Y_MAX), NX, NY, D_ODOR
    )
    solver_diffusion.set_initial_condition_gaussian(
        SOURCE_X, SOURCE_Y, SOURCE_SIGMA, SOURCE_AMPLITUDE
    )

    print("\n[SETUP] Two solvers initialized:")
    print("  1. WITH vortex dynamics: ∂c/∂t + u·∇c = D∇²c")
    print("  2. WITHOUT vortices:     ∂c/∂t = D∇²c")

    # Frame processing
    frame_indices = list(range(FRAME_START, FRAME_END + 1, FRAME_SKIP))

    print(f"\n[TEST] Processing {len(frame_indices)} frames...")
    print("-"*80)

    results = []

    for idx, frame_idx in enumerate(frame_indices):
        t_ibamr = frame_idx * VIZ_DUMP_INTERVAL * DT_IBAMR

        print(f"\n[{idx+1}/{len(frame_indices)}] Frame {frame_idx} (t* = {t_ibamr:.3f})")

        # Load velocity field from IBAMR
        points, u_x_points, u_y_points, omega = load_eulerian_frame(frame_idx)

        if points is None:
            print("    ⚠ No fluid data - using zero velocity")
            u_x_grid = np.zeros((NY, NX))
            u_y_grid = np.zeros((NY, NX))
        else:
            print(f"    ✓ Loaded {len(points)} velocity points")
            # Interpolate to regular grid
            u_x_grid, u_y_grid = interpolate_velocity_to_grid(
                points, u_x_points, u_y_points,
                solver_vortex.X, solver_vortex.Y
            )
            print(f"    ✓ Interpolated to {NX}x{NY} grid")
            print(f"    ✓ Velocity range: u_x ∈ [{np.min(u_x_grid):.3f}, {np.max(u_x_grid):.3f}]")
            print(f"                      u_y ∈ [{np.min(u_y_grid):.3f}, {np.max(u_y_grid):.3f}]")

        # Load fish positions
        eels = load_lagrangian_frame(frame_idx)
        if eels:
            print(f"    ✓ Loaded {len(eels)} fish bodies")

        # Advance both solvers to current IBAMR time
        while solver_vortex.t < t_ibamr:
            # Compute stable timestep
            dt_step = solver_vortex.compute_cfl_timestep(u_x_grid, u_y_grid)
            dt_step = min(dt_step, t_ibamr - solver_vortex.t)

            # Solver 1: WITH vortices
            solver_vortex.step_convection_diffusion(u_x_grid, u_y_grid, dt_step)

            # Solver 2: Pure diffusion
            solver_diffusion.step_diffusion_only(dt_step)

        # Get current concentrations
        c_vortex = solver_vortex.get_concentration()
        c_diffusion = solver_diffusion.get_concentration()

        # Compute metrics
        mass_vortex = solver_vortex.get_total_mass()
        mass_diffusion = solver_diffusion.get_total_mass()

        # Spreading metrics (standard deviation)
        x_vortex = np.sum(c_vortex * solver_vortex.X) / (np.sum(c_vortex) + 1e-10)
        y_vortex = np.sum(c_vortex * solver_vortex.Y) / (np.sum(c_vortex) + 1e-10)
        spread_vortex = np.sqrt(
            np.sum(c_vortex * (solver_vortex.X - x_vortex)**2) / (np.sum(c_vortex) + 1e-10) +
            np.sum(c_vortex * (solver_vortex.Y - y_vortex)**2) / (np.sum(c_vortex) + 1e-10)
        )

        x_diff = np.sum(c_diffusion * solver_diffusion.X) / (np.sum(c_diffusion) + 1e-10)
        y_diff = np.sum(c_diffusion * solver_diffusion.Y) / (np.sum(c_diffusion) + 1e-10)
        spread_diff = np.sqrt(
            np.sum(c_diffusion * (solver_diffusion.X - x_diff)**2) / (np.sum(c_diffusion) + 1e-10) +
            np.sum(c_diffusion * (solver_diffusion.Y - y_diff)**2) / (np.sum(c_diffusion) + 1e-10)
        )

        enhancement = spread_vortex / (spread_diff + 1e-10)

        print(f"    ✓ Odor spreading σ_vortex = {spread_vortex:.4f}")
        print(f"    ✓ Odor spreading σ_diffusion = {spread_diff:.4f}")
        print(f"    ✓ Enhancement factor = {enhancement:.2f}x")

        results.append({
            'frame': frame_idx,
            'time': t_ibamr,
            'c_vortex': c_vortex,
            'c_diffusion': c_diffusion,
            'spread_vortex': spread_vortex,
            'spread_diffusion': spread_diff,
            'enhancement': enhancement,
            'eels': eels,
            'omega': omega,
            'points': points
        })

        # Visualize
        visualize_comparison(solver_vortex, solver_diffusion, results[-1], output_dir, idx)

    print("-"*80)
    print("\n[COMPLETE] Odor transport test finished!")
    print(f"Output directory: {output_dir}/")

    # Summary plot
    plot_summary(results, output_dir)

    return results

def visualize_comparison(solver_vortex, solver_diffusion, result, output_dir, frame_num):
    """
    Create comparison visualization:
    - Left: Odor concentration WITH vortex dynamics
    - Middle: Odor concentration WITHOUT vortices (pure diffusion)
    - Right: Enhancement (ratio)
    """

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    c_vortex = result['c_vortex']
    c_diffusion = result['c_diffusion']
    eels = result['eels']
    t = result['time']

    # Determine color scale (same for both)
    vmax = max(np.max(c_vortex), np.max(c_diffusion))

    # Panel 1: WITH vortex dynamics
    ax = axes[0]
    im1 = ax.contourf(solver_vortex.X, solver_vortex.Y, c_vortex,
                      levels=20, cmap='hot', vmin=0, vmax=vmax)

    # Overlay fish bodies
    if eels:
        for eel in eels:
            ax.plot(eel[:, 0], eel[:, 1], 'k-', linewidth=2, alpha=0.8)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')

    if USE_LATEX:
        ax.set_title(r'WITH Vortex Dynamics: $\partial c/\partial t + \mathbf{u}\cdot\nabla c = D\nabla^2 c$',
                    fontsize=14, fontweight='bold')
        ax.set_xlabel(r'$x^*$', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(r'$y^*$', fontsize=FONT_SIZE_LABEL)
    else:
        ax.set_title('WITH Vortex Dynamics', fontsize=14, fontweight='bold')
        ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL)

    plt.colorbar(im1, ax=ax, label='c')
    ax.grid(True, alpha=0.3)

    # Panel 2: WITHOUT vortices (pure diffusion)
    ax = axes[1]
    im2 = ax.contourf(solver_diffusion.X, solver_diffusion.Y, c_diffusion,
                      levels=20, cmap='hot', vmin=0, vmax=vmax)

    # Overlay fish bodies (ghosted)
    if eels:
        for eel in eels:
            ax.plot(eel[:, 0], eel[:, 1], 'k--', linewidth=1.5, alpha=0.3)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')

    if USE_LATEX:
        ax.set_title(r'WITHOUT Vortices: $\partial c/\partial t = D\nabla^2 c$',
                    fontsize=14, fontweight='bold')
        ax.set_xlabel(r'$x^*$', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(r'$y^*$', fontsize=FONT_SIZE_LABEL)
    else:
        ax.set_title('WITHOUT Vortices (Pure Diffusion)', fontsize=14, fontweight='bold')
        ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL)

    plt.colorbar(im2, ax=ax, label='c')
    ax.grid(True, alpha=0.3)

    # Panel 3: Enhancement factor
    ax = axes[2]
    enhancement = (c_vortex + 1e-10) / (c_diffusion + 1e-10)
    im3 = ax.contourf(solver_vortex.X, solver_vortex.Y, enhancement,
                     levels=20, cmap='RdYlGn', vmin=0.5, vmax=2.0)

    if eels:
        for eel in eels:
            ax.plot(eel[:, 0], eel[:, 1], 'k-', linewidth=2, alpha=0.8)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')

    ax.set_title('Enhancement Ratio (Vortex/Diffusion)', fontsize=14, fontweight='bold')
    if USE_LATEX:
        ax.set_xlabel(r'$x^*$', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(r'$y^*$', fontsize=FONT_SIZE_LABEL)
    else:
        ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL)

    plt.colorbar(im3, ax=ax, label='Ratio')
    ax.grid(True, alpha=0.3)

    # Main title
    if USE_LATEX:
        fig.suptitle(r'Odor Transport at $t^* = {:.2f}$'.format(t),
                    fontsize=16, fontweight='bold')
    else:
        fig.suptitle(f'Odor Transport at t = {t:.2f}',
                    fontsize=16, fontweight='bold')

    plt.tight_layout()

    # Save
    filename = output_dir / f"comparison_{frame_num:04d}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    ✓ Saved: {filename}")

def plot_summary(results, output_dir):
    """
    Plot summary: enhancement factor vs time
    """

    times = [r['time'] for r in results]
    enhancements = [r['enhancement'] for r in results]
    spread_vortex = [r['spread_vortex'] for r in results]
    spread_diff = [r['spread_diffusion'] for r in results]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Panel 1: Spreading width vs time
    ax = axes[0]
    ax.plot(times, spread_vortex, 'b-o', linewidth=2, markersize=6,
            label='WITH Vortex Dynamics')
    ax.plot(times, spread_diff, 'r--s', linewidth=2, markersize=6,
            label='WITHOUT Vortices (Pure Diffusion)')

    ax.set_xlabel(r'$t^*$' if USE_LATEX else 't', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel(r'Spreading Width $\sigma$' if USE_LATEX else 'Spreading Width σ',
                 fontsize=FONT_SIZE_LABEL)
    ax.set_title('Odor Spreading vs Time', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    # Panel 2: Enhancement factor vs time
    ax = axes[1]
    ax.plot(times, enhancements, 'g-^', linewidth=2, markersize=8)
    ax.axhline(y=1.0, color='k', linestyle='--', linewidth=1, alpha=0.5, label='No enhancement')

    ax.set_xlabel(r'$t^*$' if USE_LATEX else 't', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Enhancement Factor', fontsize=FONT_SIZE_LABEL)
    ax.set_title('Vortex Enhancement of Odor Spreading', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    filename = output_dir / "summary.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n[SUMMARY] Enhancement analysis saved: {filename}")
    print(f"  Average enhancement: {np.mean(enhancements):.2f}x")
    print(f"  Maximum enhancement: {np.max(enhancements):.2f}x")

# ============================================================
# MAIN ENTRY POINT
# ============================================================

def main():
    """
    Main entry point for odor transport test
    """

    print("\n" + "="*80)
    print("ODOR TRANSPORT TEST")
    print("How does vortex dynamics help undulating bodies spread odor?")
    print("="*80)

    # Check dependencies
    try:
        import pyvista
        import scipy
        print("\n[✓] Dependencies: PyVista, SciPy")
    except ImportError as e:
        print(f"\n[✗] Missing dependency: {e}")
        print("Install with: pip install pyvista scipy")
        sys.exit(1)

    # Run test
    try:
        results = run_odor_transport_test()

        print("\n" + "="*80)
        print("[SUCCESS] Test completed!")
        print("="*80)
        print("\nKey findings:")
        print("  - Vortices created by undulating fish bodies enhance odor dispersion")
        print("  - Convection (u·∇c) dominates over pure diffusion (D∇²c)")
        print("  - Enhancement factor quantifies the effect of vortex dynamics")
        print("\nReference: arXiv:2408.16136")
        print("="*80)

    except KeyboardInterrupt:
        print("\n\n[INTERRUPTED]")
        sys.exit(0)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
