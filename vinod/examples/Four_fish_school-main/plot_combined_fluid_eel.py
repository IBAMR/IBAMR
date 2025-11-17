#!/usr/bin/env python3
"""
Create PUBLICATION-QUALITY animation: Vorticity + 4 Dark Gray Eels
WITH LaTeX rendering for professional output

Files:
- ExportEULERIANData/visit_eulerian_db__0000/*.vtk
- ExportLagrangianData/visit_lagrangian_db__00__0000.vtk

Author: Vinod
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pyvista as pv
from pathlib import Path
import sys

# ============================================================
# PUBLICATION SETTINGS WITH LaTeX
# ============================================================

USE_LATEX = True  # Set to False if LaTeX not installed

if USE_LATEX:
    try:
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}\boldmath'
        matplotlib.rcParams['font.family'] = 'serif'
        matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
        print("[INFO] LaTeX rendering ENABLED for publication quality")
    except Exception as e:
        print(f"[WARNING] LaTeX failed to initialize: {e}")
        print("[INFO] Falling back to non-LaTeX rendering")
        USE_LATEX = False
        matplotlib.rcParams['font.family'] = 'sans-serif'
else:
    matplotlib.rcParams['font.family'] = 'sans-serif'
    print("[INFO] LaTeX rendering DISABLED")

matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1

# Publication-quality font sizes
FONT_SIZE_TITLE = 20
FONT_SIZE_LABEL = 18
FONT_SIZE_TICK = 16
FONT_SIZE_COLORBAR = 16
LINEWIDTH_AXES = 2.0
LINEWIDTH_GRID = 1.0

# ============================================================
# USER CONFIGURATION
# ============================================================

# Simulation parameters
DT = 0.0001
VIZ_DUMP_INTERVAL = 40

# Which frames to process (START SMALL FOR TESTING!)
FRAME_START = 0
FRAME_END = 2500        # Change to 100 or 2500 after testing
FRAME_SKIP = 100

# FLUID VARIABLE
FLUID_VARIABLE = 'U_x'      # Options: 'Omega', 'P', 'U_x', 'U_y'
                               # Omega = vorticity (most common for eel swimming)
                               # P = pressure field
                               # U_x = x-velocity component
                               # U_y = y-velocity component

# Color scale limits for fluid field Omega
#FLUID_MIN = -10.0
#FLUID_MAX = 10.0

# Color scale limits for fluid field U_x
FLUID_MIN = -0.40
FLUID_MAX = 0.40

# Output settings
OUTPUT_DIR = "combined_frames_U_x"
SAVE_FORMAT = 'png'
FRAME_DPI = 300
FIGURE_SIZE = (14, 7)

# Domain
X_MIN = -6
X_MAX = 3
Y_MIN = -3
Y_MAX = 3

# Eel visualization - DARK GRAY (publication quality)
EEL_COLOR = '#2d2d2d'         # Dark gray (hex color)
EEL_LINEWIDTH = 4.5           # Thick lines for visibility
EEL_ALPHA = 1.0               # Fully opaque

# Position-based labels for eels (for legend)
EEL_POSITION_LABELS = {
    0: 'BL',  # Bottom-Left (rear-bottom)
    1: 'TL',  # Top-Left (rear-top)
    2: 'BR',  # Bottom-Right (front-bottom)
    3: 'TR'   # Top-Right (front-top)
}
# Adjust these based on your actual eel positions!

SHOW_EEL_LEGEND = True        # Show legend with eel labels

# Fluid field visualization
#FLUID_COLORMAP = 'RdBu_r'     # Red-White-Blue (standard in CFD) Omega  Colormap for Omega
FLUID_COLORMAP = 'seismic'     # Colormap for U_x

FLUID_ALPHA = 1.0
FLUID_POINT_SIZE = 1

# ============================================================

def get_variable_latex_name(var_name):
    """
    Get LaTeX symbol for variable name (publication standard)
    """
    names = {
        'Omega': r'$\omega^{\ast}_z$',
        'P': r'$p^{\ast}$',
        'U_x': r'$u^{\ast}_x$',
        'U_y': r'$u^{\ast}_y$'
    }
    return names.get(var_name, var_name)

def get_variable_plain_name(var_name):
    """
    Get plain text name for variable (fallback)
    """
    names = {
        'Omega': 'Vorticity',
        'P': 'Pressure',
        'U_x': 'X-Velocity',
        'U_y': 'Y-Velocity'
    }
    return names.get(var_name, var_name)

def load_eulerian_frame(frame_idx):
    """
    Load Eulerian (fluid) data for a single frame
    YOUR STRUCTURE: ExportEULERIANData/visit_eulerian_db_XXXX/*.vtk
    """
    
    # Try current directory and ExportEULERIANData
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
            
            # Combine all mesh pieces
            combined = meshes[0]
            for m in meshes[1:]:
                combined = combined.merge(m)
            
            # Check if variable exists
            if FLUID_VARIABLE not in combined.array_names:
                print(f"\n    WARNING: '{FLUID_VARIABLE}' not in frame {frame_idx}")
                print(f"    Available: {combined.array_names}")
                continue
            
            points = combined.points
            values = combined[FLUID_VARIABLE]
            
            return points, values
            
        except Exception as e:
            continue
    
    return None, None

def load_lagrangian_frame(frame_idx):
    """
    Load Lagrangian (eel) data for a single frame
    YOUR STRUCTURE: ExportLagrangianData/visit_lagrangian_db__XX__YYYY.vtk
    """
    
    eels_points = []
    
    # Try to load up to 10 eels (00, 01, 02, 03, ...)
    for eel_idx in range(10):
        # Your exact file pattern with DOUBLE underscores
        patterns = [
            f"ExportLagrangianData/visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",
            f"visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",
            # Also try single underscore in case
            f"ExportLagrangianData/visit_lagrangian_db_{eel_idx:02d}_{frame_idx:04d}.vtk",
            f"visit_lagrangian_db_{eel_idx:02d}_{frame_idx:04d}.vtk",
        ]
        
        for pattern in patterns:
            vtk_file = Path(pattern)
            if vtk_file.exists():
                try:
                    mesh = pv.read(str(vtk_file))
                    eels_points.append(mesh.points)
                    break  # Found this eel, move to next
                except Exception as e:
                    continue
        
        # If we didn't find this eel index, stop looking for higher indices
        if len(eels_points) == eel_idx:
            break
    
    return eels_points if len(eels_points) > 0 else None

def determine_fluid_range(frame_indices):
    """
    Auto-detect or use user-defined fluid variable range
    """
    
    if FLUID_MIN is not None and FLUID_MAX is not None:
        print(f"\n[INFO] Using USER-DEFINED range for {FLUID_VARIABLE}:")
        print(f"  Min: {FLUID_MIN:.2f}")
        print(f"  Max: {FLUID_MAX:.2f}")
        return FLUID_MIN, FLUID_MAX
    
    print(f"\n[INFO] Auto-detecting range for {FLUID_VARIABLE}...")
    
    global_min = 0
    global_max = 0
    
    # Sample a few frames
    sample_indices = [0, len(frame_indices)//4, len(frame_indices)//2,
                     3*len(frame_indices)//4, -1]
    sample_frames = [frame_indices[i] for i in sample_indices if i < len(frame_indices)]
    
    for frame_idx in sample_frames:
        points, values = load_eulerian_frame(frame_idx)
        if values is not None:
            global_min = min(global_min, values.min())
            global_max = max(global_max, values.max())
            print(f"  Frame {frame_idx}: min={values.min():.4f}, max={values.max():.4f}")
    
    vmax = max(abs(global_min), abs(global_max))
    print(f"\n[INFO] AUTO-DETECTED range: [{-vmax:.2f}, {vmax:.2f}]")
    
    return -vmax, vmax

def create_animation():
    """
    Main function to create combined animation
    """
    
    print("\n" + "="*80)
    print("CREATING PUBLICATION-QUALITY VORTICITY + DARK GRAY EELS ANIMATION")
    print("="*80)
    
    frame_indices = list(range(FRAME_START, FRAME_END + 1, FRAME_SKIP))
    
    actual_timestep_start = FRAME_START * VIZ_DUMP_INTERVAL
    actual_timestep_end = FRAME_END * VIZ_DUMP_INTERVAL
    t_start = actual_timestep_start * DT
    t_end = actual_timestep_end * DT
    
    print(f"\n[CONFIG] Settings:")
    print(f"  LaTeX rendering: {'ENABLED ✓' if USE_LATEX else 'DISABLED'}")
    print(f"  Fluid variable: {FLUID_VARIABLE} (Options: Omega, P, U_x, U_y)")
    print(f"  Frames: {FRAME_START} to {FRAME_END} (every {FRAME_SKIP})")
    print(f"  Total frames: {len(frame_indices)}")
    print(f"  Time range: t* = {t_start:.3f} to {t_end:.3f}")
    print(f"  Domain: X=[{X_MIN}, {X_MAX}], Y=[{Y_MIN}, {Y_MAX}]")
    print(f"  Eel color: {EEL_COLOR} (dark gray)")
    print(f"  Eel line width: {EEL_LINEWIDTH}")
    print(f"  Eel labels: {EEL_POSITION_LABELS}")
    print(f"  Show legend: {SHOW_EEL_LEGEND}")
    
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)
    print(f"  Output: {output_dir}")
    
    # Determine fluid color scale
    fluid_min, fluid_max = determine_fluid_range(frame_indices)
    
    print("\n[INFO] Creating frames...")
    print("-"*80)
    
    successful_frames = 0
    failed_frames = 0
    
    for i, frame_idx in enumerate(frame_indices):
        actual_timestep = frame_idx * VIZ_DUMP_INTERVAL
        t_nondim = actual_timestep * DT
        
        print(f"  [{i+1:3d}/{len(frame_indices):3d}] Frame {frame_idx:4d} (t* = {t_nondim:.3f})", end='')
        
        # Load both Eulerian and Lagrangian data
        fluid_points, fluid_values = load_eulerian_frame(frame_idx)
        eels_points = load_lagrangian_frame(frame_idx)
        
        if fluid_points is None:
            print(" - SKIP (no fluid data)")
            failed_frames += 1
            continue
        
        if eels_points is None:
            print(" - WARNING (no eel data, fluid only)")
        
        try:
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            
            # Plot fluid field (vorticity) as background
            scatter = ax.scatter(
                fluid_points[:, 0],
                fluid_points[:, 1],
                c=fluid_values,
                cmap=FLUID_COLORMAP,
                s=FLUID_POINT_SIZE,
                alpha=FLUID_ALPHA,
                vmin=fluid_min,
                vmax=fluid_max,
                zorder=1,
                rasterized=True  # Better for large datasets
            )
            
            # Overlay each eel body in DARK GRAY
            if eels_points is not None:
                for eel_idx, points in enumerate(eels_points):
                    # Use position label from dictionary
                    label = EEL_POSITION_LABELS.get(eel_idx, f'Eel {eel_idx + 1}')
                    
                    ax.plot(points[:, 0], points[:, 1],
                           color=EEL_COLOR,
                           linewidth=EEL_LINEWIDTH,
                           alpha=EEL_ALPHA,
                           label=label if SHOW_EEL_LEGEND else None,
                           zorder=10,  # Foreground
                           solid_capstyle='round')
                
                # Add legend if enabled and multiple eels
                if SHOW_EEL_LEGEND and len(eels_points) > 1:
                    ax.legend(loc='upper left', fontsize=12, framealpha=0.95, 
                             edgecolor='black', fancybox=False)
            
            # Set domain
            ax.set_xlim(X_MIN, X_MAX)
            ax.set_ylim(Y_MIN, Y_MAX)
            
            # Labels and title (LaTeX or plain)
            if USE_LATEX:
                ax.set_xlabel(r'$x^{\ast}$', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_ylabel(r'$y^{\ast}$', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                title = r'{} at $t^{{\ast}} = {:.2f}$'.format(
                    get_variable_latex_name(FLUID_VARIABLE), t_nondim)
                ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
            else:
                ax.set_xlabel('x*', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_ylabel('y*', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                title = '{} at t* = {:.2f}'.format(
                    get_variable_plain_name(FLUID_VARIABLE), t_nondim)
                ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            ax.set_aspect('equal', adjustable='box')
            ax.grid(True, alpha=0.3, linewidth=LINEWIDTH_GRID, linestyle='--')
            
            ax.tick_params(axis='both', which='major',
                          labelsize=FONT_SIZE_TICK,
                          width=LINEWIDTH_AXES, length=6)
            
            for spine in ax.spines.values():
                spine.set_linewidth(LINEWIDTH_AXES)
            
            # Colorbar
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
            
            if USE_LATEX:
                cbar.set_label(get_variable_latex_name(FLUID_VARIABLE),
                             fontsize=FONT_SIZE_COLORBAR,
                             rotation=0, labelpad=20)
            else:
                cbar.set_label(get_variable_plain_name(FLUID_VARIABLE),
                             fontsize=FONT_SIZE_COLORBAR)
            
            cbar.ax.tick_params(labelsize=FONT_SIZE_COLORBAR-2,
                              width=LINEWIDTH_AXES)
            
            # Add text showing number of eels (if no legend)
            if eels_points is not None and not SHOW_EEL_LEGEND:
                ax.text(0.02, 0.98, f'{len(eels_points)} Eels', 
                       transform=ax.transAxes,
                       fontsize=12, fontweight='bold',
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Save with high quality
            frame_file = output_dir / f"frame_{i:04d}.{SAVE_FORMAT}"
            plt.savefig(frame_file, dpi=FRAME_DPI,
                       bbox_inches='tight', format=SAVE_FORMAT)
            plt.close(fig)
            
            successful_frames += 1
            print(" - OK")
            
        except Exception as e:
            print(f" - ERROR: {e}")
            failed_frames += 1
            plt.close('all')
            continue
    
    print("-"*80)
    
    print("\n" + "="*80)
    print("[COMPLETE]")
    print("="*80)
    print(f"  Successful: {successful_frames}")
    print(f"  Failed: {failed_frames}")
    
    if successful_frames > 0:
        print(f"\n[SUCCESS] Publication-quality frames created!")
        print(f"  Output: {output_dir}/")
        print(f"\n[NEXT STEP] Create video:")
        print(f"  cd {output_dir}")
        print(f"  ffmpeg -framerate 10 -i frame_%04d.png \\")
        print(f"         -c:v libx264 -crf 18 -pix_fmt yuv420p ../animation.mp4")
        print("="*80)
        return True
    else:
        print("\n[ERROR] No frames created!")
        return False

def main():
    print("\n" + "="*80)
    print("PUBLICATION-QUALITY VORTICITY + DARK GRAY EELS VISUALIZATION")
    print("="*80)
    
    # Check dependencies
    try:
        import pyvista
        print("\n[OK] PyVista is installed")
    except ImportError:
        print("\n[ERROR] PyVista not installed!")
        print("Install with: pip install pyvista")
        sys.exit(1)
    
    # Check LaTeX installation
    if USE_LATEX:
        print("[INFO] LaTeX rendering enabled - checking installation...")
        print("       If you get errors, install LaTeX:")
        print("       Ubuntu/Debian: sudo apt-get install texlive-latex-extra dvipng")
        print("       Or set USE_LATEX = False in line 23")
    
    # Check directories
    eulerian_exists = Path("ExportEULERIANData").exists()
    lagrangian_exists = Path("ExportLagrangianData").exists()
    
    print(f"\n[CHECK] Directories:")
    print(f"  ExportEULERIANData:   {'EXISTS ✓' if eulerian_exists else 'NOT FOUND ✗'}")
    print(f"  ExportLagrangianData: {'EXISTS ✓' if lagrangian_exists else 'NOT FOUND ✗'}")
    
    if not eulerian_exists:
        print("\n[ERROR] ExportEULERIANData/ not found in current directory!")
        print("Please run from: /mnt/d/September2025/29september/XJWang/Fourfish/T1/")
        sys.exit(1)
    
    try:
        success = create_animation()
        sys.exit(0 if success else 1)
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