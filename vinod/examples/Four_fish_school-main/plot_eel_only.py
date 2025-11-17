#!/usr/bin/env python3
"""
Create publication-quality animation of ONLY the eel body (Lagrangian structure)
No fluid field - just the eel swimming motion

Author: Generated for IBAMR visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pyvista as pv
from pathlib import Path
import sys

# ============================================================
# PUBLICATION SETTINGS
# ============================================================

USE_LATEX = True

if USE_LATEX:
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}\boldmath'
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']

matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1

FONT_SIZE_TITLE = 20
FONT_SIZE_LABEL = 18
FONT_SIZE_TICK = 16
LINEWIDTH_AXES = 2.0
LINEWIDTH_GRID = 1.0

# ============================================================
# USER CONFIGURATION
# ============================================================

# Simulation parameters
DT = 0.0001                   # Timestep size
VIZ_DUMP_INTERVAL = 40        # Viz dump interval

# Which frames to process
FRAME_START = 0               # First frame
FRAME_END = 500              # Last frame
FRAME_SKIP = 1                # Process every Nth frame

# Output settings
OUTPUT_DIR = "eel_only_frames"
SAVE_FORMAT = 'png'
FRAME_DPI = 300
FIGURE_SIZE = (12, 6)

# Plot domain (adjust to zoom in/out)
X_MIN = -6
X_MAX = 3
Y_MIN = -3
Y_MAX = 3

# Eel visualization settings (optimized for 4 eels)
EEL_COLORS = ['black', 'red', 'blue', 'green']  # Distinct colors for 4 eels
EEL_LINEWIDTH = 2.5           # Slightly thinner so 4 eels don't overlap visually
EEL_MARKER_SIZE = 2           # Size of dots at Lagrangian points
SHOW_POINTS = False           # Show individual Lagrangian points?

# Background
BACKGROUND_COLOR = 'white'    # 'white' or 'lightgray'

# ============================================================

def load_lagrangian_frame(frame_idx):
    """
    Load Lagrangian (eel body) data for a single frame
    Handles multiple eels (files ending with _00_, _01_, etc.)
    
    Parameters:
    -----------
    frame_idx : int
        Frame number (e.g., 2250 for files like visit_lagrangian_db_00_2250.vtk)
    
    Returns:
    --------
    eels_points : list of numpy arrays
        List of eel body coordinates, one array per eel
    """
    
    eels_points = []
    
    # Try to load multiple eels (00, 01, 02, etc.)
    for eel_idx in range(10):  # Support up to 10 eels
        # Try naming patterns (both single _ and double __ underscore)
        # Also check in ExportLagrangianData subdirectory
        patterns = [
            f"visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",  # Current dir, double __
            f"visit_lagrangian_db_{eel_idx:02d}_{frame_idx:04d}.vtk",   # Current dir, single _
            f"ExportLagrangianData/visit_lagrangian_db__{eel_idx:02d}__{frame_idx:04d}.vtk",  # Subfolder, double __
            f"ExportLagrangianData/visit_lagrangian_db_{eel_idx:02d}_{frame_idx:04d}.vtk",   # Subfolder, single _
            f"visit_lagrangian_db_{frame_idx:04d}.vtk" if eel_idx == 0 else None  # Single eel fallback
        ]
        
        for pattern in patterns:
            if pattern is None:
                continue
                
            vtk_file = Path(pattern)
            
            if vtk_file.exists():
                try:
                    mesh = pv.read(str(vtk_file))
                    points = mesh.points
                    eels_points.append(points)
                    break  # Found this eel, move to next
                except Exception as e:
                    continue
        
        # If we didn't find this eel index, assume no more eels
        if eel_idx > 0 and len(eels_points) == eel_idx:
            break
    
    return eels_points if len(eels_points) > 0 else None

def create_animation():
    """
    Main function to create eel-only animation
    """
    
    print("\n" + "="*70)
    print("CREATING EEL BODY ANIMATION")
    print("="*70)
    
    # Create list of frames to process
    frame_indices = list(range(FRAME_START, FRAME_END + 1, FRAME_SKIP))
    
    # Calculate actual timesteps
    actual_timestep_start = FRAME_START * VIZ_DUMP_INTERVAL
    actual_timestep_end = FRAME_END * VIZ_DUMP_INTERVAL
    t_start = actual_timestep_start * DT
    t_end = actual_timestep_end * DT
    
    print("\n[CONFIG] Settings:")
    print("  Frames: {} to {} (every {})".format(FRAME_START, FRAME_END, FRAME_SKIP))
    print("  Total frames: {}".format(len(frame_indices)))
    print("  Time range: t* = {:.3f} to {:.3f}".format(t_start, t_end))
    print("  Domain: X=[{}, {}], Y=[{}, {}]".format(X_MIN, X_MAX, Y_MIN, Y_MAX))
    print("  Output: {} at {} DPI".format(SAVE_FORMAT.upper(), FRAME_DPI))
    
    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)
    print("  Output directory: {}".format(output_dir))
    
    print("\n[INFO] Creating frames...")
    print("-"*70)
    
    # Detect number of eels from first frame
    test_eels = load_lagrangian_frame(FRAME_START)
    if test_eels:
        print(f"  🐟 Detected {len(test_eels)} eels in the simulation!")
        for idx, points in enumerate(test_eels):
            print(f"     Eel {idx+1}: {len(points)} Lagrangian points")
        print("-"*70)
    
    successful_frames = 0
    failed_frames = 0
    
    for i, frame_idx in enumerate(frame_indices):
        # Calculate time
        actual_timestep = frame_idx * VIZ_DUMP_INTERVAL
        t_nondim = actual_timestep * DT
        
        print("  [{:3d}/{:3d}] Frame {:4d} (t* = {:.3f})".format(
            i+1, len(frame_indices), frame_idx, t_nondim), end='')
        
        # Load Lagrangian data (may contain multiple eels)
        eels_points = load_lagrangian_frame(frame_idx)
        
        if eels_points is None:
            print(" - SKIP")
            failed_frames += 1
            continue
        
        try:
            # Create figure
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            
            # Set background color
            ax.set_facecolor(BACKGROUND_COLOR)
            
            # Plot each eel
            for eel_idx, points in enumerate(eels_points):
                color = EEL_COLORS[eel_idx % len(EEL_COLORS)]
                label = f'Eel {eel_idx + 1}' if len(eels_points) > 1 else 'Eel body'
                
                ax.plot(points[:, 0], points[:, 1], 
                       color=color, 
                       linewidth=EEL_LINEWIDTH,
                       linestyle='-',
                       marker='o' if SHOW_POINTS else None,
                       markersize=EEL_MARKER_SIZE if SHOW_POINTS else 0,
                       label=label)
            
            # Add legend if multiple eels
            if len(eels_points) > 1:
                ax.legend(loc='upper right', fontsize=11, framealpha=0.9, ncol=2 if len(eels_points) > 2 else 1)
            
            # Set domain limits
            ax.set_xlim(X_MIN, X_MAX)
            ax.set_ylim(Y_MIN, Y_MAX)
            
            # Labels and title
            if USE_LATEX:
                ax.set_xlabel(r'$x^{\ast}$', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_ylabel(r'$y^{\ast}$', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_title(r'Eel Body at $t^{\ast} = %.2f$' % t_nondim, 
                           fontsize=FONT_SIZE_TITLE, fontweight='bold')
            else:
                ax.set_xlabel('x', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_ylabel('y', fontsize=FONT_SIZE_LABEL, fontweight='bold')
                ax.set_title('Eel Body at t = {:.2f}'.format(t_nondim), 
                           fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            # Equal aspect ratio
            ax.set_aspect('equal', adjustable='box')
            
            # Grid
            ax.grid(True, alpha=0.3, linewidth=LINEWIDTH_GRID, linestyle='--')
            
            # Tick parameters
            ax.tick_params(axis='both', which='major', 
                          labelsize=FONT_SIZE_TICK,
                          width=LINEWIDTH_AXES, length=6)
            
            # Spine thickness
            for spine in ax.spines.values():
                spine.set_linewidth(LINEWIDTH_AXES)
            
            # Save figure
            frame_file = output_dir / "frame_{:04d}.{}".format(i, SAVE_FORMAT)
            plt.savefig(frame_file, dpi=FRAME_DPI, bbox_inches='tight', format=SAVE_FORMAT)
            plt.close(fig)
            
            successful_frames += 1
            
            # Count total points across all eels
            total_points = sum(len(points) for points in eels_points)
            print(" - OK ({} eels, {} total points)".format(len(eels_points), total_points))
            
        except Exception as e:
            print(" - ERROR: {}".format(e))
            failed_frames += 1
            plt.close('all')
            continue
    
    print("-"*70)
    
    # Summary
    print("\n" + "="*70)
    print("[COMPLETE]")
    print("="*70)
    print("  Successful: {}".format(successful_frames))
    print("  Failed: {}".format(failed_frames))
    print("  Output: {}".format(output_dir))
    
    if successful_frames > 0:
        print("\n[NEXT STEP] Create video:")
        print("  cd {}".format(output_dir))
        if SAVE_FORMAT == 'png':
            print("  ffmpeg -framerate 10 -i frame_%04d.png \\")
            print("         -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" \\")
            print("         -c:v libx264 -pix_fmt yuv420p ../eel_only.mp4")
        print("="*70)
        return True
    else:
        print("\n[ERROR] No frames created!")
        return False

def main():
    """
    Main entry point
    """
    
    print("\n" + "="*70)
    print("EEL BODY VISUALIZATION")
    print("="*70)
    
    print("\n[CONFIG] Visualization settings:")
    print("  Eel colors: {}".format(EEL_COLORS))
    print("  Line width: {}".format(EEL_LINEWIDTH))
    print("  Show points: {}".format(SHOW_POINTS))
    print("  Background: {}".format(BACKGROUND_COLOR))
    
    try:
        success = create_animation()
        sys.exit(0 if success else 1)
        
    except KeyboardInterrupt:
        print("\n\n[INTERRUPTED]")
        sys.exit(0)
        
    except Exception as e:
        print("\n[ERROR] {}".format(e))
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()