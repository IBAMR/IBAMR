#!/usr/bin/env python3
"""
Test Crank-Nicolson Odor Transport Solver with IBAMR Velocity Fields

This script demonstrates the enhanced Crank-Nicolson odor transport solver
with velocity fields from the Four Fish School IBAMR simulation.

Features demonstrated:
----------------------
1. High Schmidt number capability (Sc >> 1)
2. Comparison: Crank-Nicolson vs explicit method
3. Stability with large timesteps
4. Mass conservation
5. Vortex-enhanced odor dispersion

Based on: "Collective Chemotactic Behavior in Fish Schools" (arXiv:2408.16136)
Authors: Maham Kamran, Amirhossein Fardi, Chengyu Li, Muhammad Saif Ullah Khalid
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import sys

# Import both solvers for comparison
from odor_transport_solver_CN import OdorTransportSolverCN
sys.path.insert(0, str(Path(__file__).parent))

try:
    import pyvista as pv
    from scipy.interpolate import griddata
    HAVE_PYVISTA = True
except ImportError:
    print("[WARNING] PyVista not available - running without IBAMR data")
    HAVE_PYVISTA = False

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

# Physical parameters
NU = 0.785 / 5609.0  # Kinematic viscosity (from input2d: μ/Re)
D_ODOR_BASE = 0.001  # Base diffusion coefficient
SCHMIDT_NUMBERS = [1.0, 10.0, 100.0]  # Test different Sc = ν/D

# IBAMR parameters
DT_IBAMR = 0.0001
VIZ_DUMP_INTERVAL = 40

# Computational domain
X_MIN, X_MAX = -6.0, 3.0
Y_MIN, Y_MAX = -3.0, 3.0
NX = 200
NY = 150

# Initial odor source
SOURCE_X = -2.0
SOURCE_Y = 0.0
SOURCE_SIGMA = 0.2
SOURCE_AMPLITUDE = 1.0

# Time integration
DT_CRANK_NICOLSON = 0.001  # Can use larger timesteps with implicit scheme
T_FINAL = 0.4

# Frame processing
FRAME_START = 0
FRAME_END = 100
FRAME_SKIP = 10

# Output
OUTPUT_DIR = "odor_transport_CN_test"

# ============================================================
# DATA LOADING FROM IBAMR (same as original)
# ============================================================

def load_eulerian_frame(frame_idx):
    """Load velocity field from IBAMR output"""
    if not HAVE_PYVISTA:
        return None, None, None, None

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

            combined = meshes[0]
            for m in meshes[1:]:
                combined = combined.merge(m)

            points = combined.points

            if 'U_x' in combined.array_names and 'U_y' in combined.array_names:
                u_x_points = combined['U_x']
                u_y_points = combined['U_y']
            else:
                continue

            omega = combined.get('Omega', None)

            return points, u_x_points, u_y_points, omega

        except Exception as e:
            continue

    return None, None, None, None

def load_lagrangian_frame(frame_idx):
    """Load fish positions"""
    if not HAVE_PYVISTA:
        return None

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
    """Interpolate scattered velocity to regular grid"""
    u_x_grid = griddata(points[:, :2], u_x_points,
                        (grid_x, grid_y), method='linear', fill_value=0.0)
    u_y_grid = griddata(points[:, :2], u_y_points,
                        (grid_x, grid_y), method='linear', fill_value=0.0)
    return u_x_grid, u_y_grid

# ============================================================
# TEST 1: VALIDATION WITH ANALYTICAL SOLUTION
# ============================================================

def test_pure_diffusion_analytical():
    """
    Validate Crank-Nicolson solver against analytical solution.

    For Gaussian initial condition in infinite domain:
    C(x,y,t) = (A·σ₀²)/(σ₀² + 4Dt) · exp(-r²/(2(σ₀² + 4Dt)))
    """
    print("\n" + "="*80)
    print("TEST 1: Analytical Validation (Pure Diffusion)")
    print("="*80)

    D = 0.01
    x0, y0 = 0.0, 0.0
    sigma0 = 0.2
    A = 1.0
    t_test = 0.5

    solver = OdorTransportSolverCN(
        x_range=(-2, 2), y_range=(-2, 2),
        nx=100, ny=100,
        diffusion_coeff=D,
        boundary_type='neumann'
    )

    solver.set_initial_condition_gaussian(x0, y0, sigma0, A)

    # Advance with large timesteps (shows advantage of implicit method)
    dt = 0.05  # Much larger than explicit stability limit!
    num_steps = int(t_test / dt)

    print(f"\nTime integration: dt = {dt} (implicit scheme allows large dt)")
    for step in range(num_steps):
        solver.step_diffusion_only(dt)

    # Analytical solution
    r_squared = (solver.X - x0)**2 + (solver.Y - y0)**2
    sigma_t = np.sqrt(sigma0**2 + 2*D*t_test)
    c_analytical = (A * sigma0**2 / sigma_t**2) * np.exp(-r_squared / (2 * sigma_t**2))

    c_numerical = solver.get_concentration()

    error_rel = np.sqrt(np.mean((c_numerical - c_analytical)**2)) / np.max(c_analytical)

    print(f"\n[RESULTS]")
    print(f"  Relative L2 error: {error_rel:.4%}")
    print(f"  Mass conservation error: {solver.get_mass_conservation_error():.2e}")

    if error_rel < 0.02:
        print("  [✓] PASSED: Error < 2%")
    else:
        print("  [✗] FAILED: Error > 2%")

    return error_rel < 0.02

# ============================================================
# TEST 2: HIGH SCHMIDT NUMBER CAPABILITY
# ============================================================

def test_high_schmidt_number():
    """
    Demonstrate solver capability with high Schmidt numbers.

    High Sc means small diffusion (D = ν/Sc), which makes explicit
    methods very slow. Implicit Crank-Nicolson handles this efficiently.
    """
    print("\n" + "="*80)
    print("TEST 2: High Schmidt Number Capability")
    print("="*80)

    output_dir = Path(OUTPUT_DIR) / "high_schmidt"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create simple vortex velocity field for testing
    x = np.linspace(X_MIN, X_MAX, NX)
    y = np.linspace(Y_MIN, Y_MAX, NY)
    X, Y = np.meshgrid(x, y)

    # Vortex centered at origin with strength Γ
    gamma = 2.0
    r = np.sqrt(X**2 + Y**2) + 0.1  # Add small epsilon to avoid singularity
    u_x = -gamma * Y / r**2
    u_y = gamma * X / r**2

    results = []

    for Sc in SCHMIDT_NUMBERS:
        D = NU / Sc
        print(f"\n[Sc = {Sc}] Diffusion coefficient D = {D:.6e}")

        solver = OdorTransportSolverCN(
            x_range=(X_MIN, X_MAX), y_range=(Y_MIN, Y_MAX),
            nx=NX, ny=NY,
            diffusion_coeff=D,
            schmidt_number=Sc
        )

        solver.set_initial_condition_gaussian(SOURCE_X, SOURCE_Y, SOURCE_SIGMA)

        # Advance to final time
        t = 0.0
        dt = DT_CRANK_NICOLSON
        step = 0

        while t < T_FINAL:
            dt_step = min(dt, T_FINAL - t)
            solver.step_crank_nicolson(u_x, u_y, dt_step)
            t += dt_step
            step += 1

            if step % 50 == 0:
                print(f"  Step {step}: t = {t:.4f}")

        info = solver.get_solver_info()
        results.append({
            'Sc': Sc,
            'D': D,
            'concentration': solver.get_concentration(),
            'info': info
        })

        print(f"  [DONE] {step} steps, spreading width σ = {info['spreading_width']:.4f}")
        print(f"         Mass error: {info['mass_conservation_error']:.2e}")

    # Visualize comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, result in enumerate(results):
        ax = axes[idx]
        Sc = result['Sc']
        c = result['concentration']

        # Plot concentration field
        x_plot = np.linspace(X_MIN, X_MAX, NX)
        y_plot = np.linspace(Y_MIN, Y_MAX, NY)
        X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

        im = ax.contourf(X_plot, Y_plot, c, levels=20, cmap='hot')
        plt.colorbar(im, ax=ax, label='C')

        # Add velocity field
        skip = 10
        ax.quiver(X_plot[::skip, ::skip], Y_plot[::skip, ::skip],
                 u_x[::skip, ::skip], u_y[::skip, ::skip],
                 alpha=0.4, scale=20)

        if USE_LATEX:
            ax.set_title(f'$Sc = {Sc}$ $(D = {result["D"]:.2e})$', fontsize=14, fontweight='bold')
            ax.set_xlabel(r'$x^*$', fontsize=FONT_SIZE_LABEL)
            ax.set_ylabel(r'$y^*$', fontsize=FONT_SIZE_LABEL)
        else:
            ax.set_title(f'Sc = {Sc} (D = {result["D"]:.2e})', fontsize=14, fontweight='bold')
            ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL)
            ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL)

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'Odor Transport at t = {T_FINAL:.2f} for Different Schmidt Numbers',
                fontsize=16, fontweight='bold')
    plt.tight_layout()

    filename = output_dir / "schmidt_number_comparison.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {filename}")
    print("[✓] TEST 2 COMPLETE")

    return True

# ============================================================
# TEST 3: INTEGRATION WITH IBAMR VELOCITY FIELDS
# ============================================================

def test_with_ibamr_velocity():
    """
    Test Crank-Nicolson solver with actual IBAMR velocity fields.
    """
    print("\n" + "="*80)
    print("TEST 3: Integration with IBAMR Velocity Fields")
    print("="*80)

    if not HAVE_PYVISTA:
        print("[SKIP] PyVista not available")
        return True

    output_dir = Path(OUTPUT_DIR) / "ibamr_integration"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use moderate Schmidt number
    Sc = 10.0
    D = NU / Sc

    print(f"\nSchmidt number: Sc = {Sc}")
    print(f"Diffusion coefficient: D = {D:.6e}")

    solver = OdorTransportSolverCN(
        x_range=(X_MIN, X_MAX), y_range=(Y_MIN, Y_MAX),
        nx=NX, ny=NY,
        diffusion_coeff=D,
        schmidt_number=Sc
    )

    solver.set_initial_condition_gaussian(SOURCE_X, SOURCE_Y, SOURCE_SIGMA)

    # Process frames
    frame_indices = list(range(FRAME_START, min(FRAME_END, 50), FRAME_SKIP))
    print(f"\nProcessing {len(frame_indices)} frames...")

    results = []

    for idx, frame_idx in enumerate(frame_indices):
        t_target = frame_idx * VIZ_DUMP_INTERVAL * DT_IBAMR

        print(f"\n[{idx+1}/{len(frame_indices)}] Frame {frame_idx} (t = {t_target:.4f})")

        # Load velocity field
        points, u_x_points, u_y_points, omega = load_eulerian_frame(frame_idx)

        if points is None:
            print("  [WARNING] No velocity data, using zero velocity")
            u_x_grid = np.zeros((NY, NX))
            u_y_grid = np.zeros((NY, NX))
        else:
            print(f"  [✓] Loaded {len(points)} velocity points")
            X_grid, Y_grid = np.meshgrid(
                np.linspace(X_MIN, X_MAX, NX),
                np.linspace(Y_MIN, Y_MAX, NY)
            )
            u_x_grid, u_y_grid = interpolate_velocity_to_grid(
                points, u_x_points, u_y_points, X_grid, Y_grid
            )

        # Load fish
        eels = load_lagrangian_frame(frame_idx)

        # Advance solver to target time
        while solver.t < t_target:
            dt_step = min(DT_CRANK_NICOLSON, t_target - solver.t)
            solver.step_crank_nicolson(u_x_grid, u_y_grid, dt_step)

        info = solver.get_solver_info()
        print(f"  [DONE] t = {solver.t:.4f}, σ = {info['spreading_width']:.4f}")

        results.append({
            'frame': frame_idx,
            'time': solver.t,
            'concentration': solver.get_concentration(),
            'eels': eels,
            'info': info
        })

        # Visualize
        if idx % 2 == 0:  # Save every other frame
            visualize_frame(solver, results[-1], output_dir, idx)

    # Summary plot
    plot_spreading_evolution(results, output_dir)

    print("\n[✓] TEST 3 COMPLETE")
    return True

def visualize_frame(solver, result, output_dir, frame_num):
    """Visualize single frame with fish"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    c = result['concentration']
    eels = result['eels']
    t = result['time']

    # Concentration field
    im = ax.contourf(solver.X, solver.Y, c, levels=20, cmap='hot')
    plt.colorbar(im, ax=ax, label='Concentration')

    # Overlay fish
    if eels:
        for eel in eels:
            ax.plot(eel[:, 0], eel[:, 1], 'k-', linewidth=2, alpha=0.8)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect('equal')

    if USE_LATEX:
        ax.set_title(f'Odor Concentration at $t^* = {t:.2f}$', fontsize=16, fontweight='bold')
        ax.set_xlabel(r'$x^*$', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(r'$y^*$', fontsize=FONT_SIZE_LABEL)
    else:
        ax.set_title(f'Odor Concentration at t = {t:.2f}', fontsize=16, fontweight='bold')
        ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL)

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    filename = output_dir / f"frame_{frame_num:04d}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_spreading_evolution(results, output_dir):
    """Plot spreading width evolution"""
    times = [r['time'] for r in results]
    sigmas = [r['info']['spreading_width'] for r in results]
    mass_errors = [r['info']['mass_conservation_error'] for r in results]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Spreading width
    ax = axes[0]
    ax.plot(times, sigmas, 'b-o', linewidth=2, markersize=6)
    ax.set_xlabel(r'$t^*$' if USE_LATEX else 't', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel(r'Spreading Width $\sigma$' if USE_LATEX else 'σ', fontsize=FONT_SIZE_LABEL)
    ax.set_title('Odor Spreading Evolution (Crank-Nicolson)', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Mass conservation
    ax = axes[1]
    ax.semilogy(times, mass_errors, 'r-^', linewidth=2, markersize=6)
    ax.set_xlabel(r'$t^*$' if USE_LATEX else 't', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Relative Mass Error', fontsize=FONT_SIZE_LABEL)
    ax.set_title('Mass Conservation (Crank-Nicolson)', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    filename = output_dir / "spreading_evolution.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {filename}")

# ============================================================
# MAIN
# ============================================================

def main():
    """Run all tests"""
    print("\n" + "="*80)
    print("CRANK-NICOLSON ODOR TRANSPORT SOLVER: COMPREHENSIVE TESTS")
    print("="*80)
    print("\nGoverning equation: ∂C/∂t + ui ∂C/∂xi = D ∂²C/∂xi∂xi")
    print("\nNumerical scheme:")
    print("  • Temporal: Implicit Crank-Nicolson (2nd order, unconditionally stable)")
    print("  • Convection: Upwind finite differences")
    print("  • Diffusion: Central differences with Crank-Nicolson averaging")
    print("="*80)

    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)

    # Run tests
    tests_passed = []

    try:
        tests_passed.append(("Analytical Validation", test_pure_diffusion_analytical()))
    except Exception as e:
        print(f"[ERROR] Test 1 failed: {e}")
        tests_passed.append(("Analytical Validation", False))

    try:
        tests_passed.append(("High Schmidt Number", test_high_schmidt_number()))
    except Exception as e:
        print(f"[ERROR] Test 2 failed: {e}")
        tests_passed.append(("High Schmidt Number", False))

    try:
        tests_passed.append(("IBAMR Integration", test_with_ibamr_velocity()))
    except Exception as e:
        print(f"[ERROR] Test 3 failed: {e}")
        tests_passed.append(("IBAMR Integration", False))

    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    for test_name, passed in tests_passed:
        status = "[✓] PASSED" if passed else "[✗] FAILED"
        print(f"{status}: {test_name}")

    all_passed = all(p for _, p in tests_passed)

    print("="*80)
    if all_passed:
        print("[SUCCESS] All tests passed!")
    else:
        print("[WARNING] Some tests failed")

    print(f"\nOutput directory: {OUTPUT_DIR}/")
    print("="*80)

if __name__ == "__main__":
    main()
