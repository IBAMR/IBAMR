#!/usr/bin/env python3
"""
Performance Analysis Tool for Undulatory Foil Propulsion
=========================================================

Analyzes the effects of Reynolds number and thickness on swimming performance
of an undulatory self-propelled foil with adaptive kinematics.

Usage:
    python analyze_performance.py [performance_files...]

If no files are specified, it will search for all performance_*.dat files
in the current directory.

Author: IBAMR Implementation
Date: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from pathlib import Path

class PerformanceAnalyzer:
    """Analyzes performance metrics from undulatory foil simulations"""

    def __init__(self, filename):
        """Initialize analyzer with a performance data file"""
        self.filename = filename
        self.data = None
        self.reynolds_number = None
        self.thickness_ratio = None
        self.swimming_mode = None
        self.load_data()

    def load_data(self):
        """Load performance data from file"""
        # Read header to extract parameters
        with open(self.filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'Reynolds number:' in line:
                    self.reynolds_number = float(line.split(':')[1].strip())
                elif 'Thickness ratio:' in line:
                    self.thickness_ratio = float(line.split(':')[1].strip())
                elif 'Swimming mode:' in line:
                    self.swimming_mode = float(line.split(':')[1].strip())

        # Load numerical data
        self.data = np.loadtxt(self.filename, comments='#')

        # Extract columns
        self.time = self.data[:, 0]
        self.adapted_amplitude = self.data[:, 1]
        self.adapted_frequency = self.data[:, 2]
        self.swimming_speed = self.data[:, 3]
        self.thrust = self.data[:, 4]
        self.power = self.data[:, 5]
        self.efficiency = self.data[:, 6]

    def compute_statistics(self):
        """Compute time-averaged statistics after initial transient"""
        # Skip first 20% as transient
        start_idx = int(0.2 * len(self.time))

        stats = {
            'Re': self.reynolds_number,
            'h/L': self.thickness_ratio,
            'mode': self.swimming_mode,
            'avg_speed': np.mean(self.swimming_speed[start_idx:]),
            'avg_thrust': np.mean(self.thrust[start_idx:]),
            'avg_power': np.mean(self.power[start_idx:]),
            'avg_efficiency': np.mean(self.efficiency[start_idx:]),
            'avg_amplitude': np.mean(self.adapted_amplitude[start_idx:]),
            'avg_frequency': np.mean(self.adapted_frequency[start_idx:]),
            'std_speed': np.std(self.swimming_speed[start_idx:]),
            'std_efficiency': np.std(self.efficiency[start_idx:]),
        }

        # Compute Strouhal number: St = f * A / U
        if stats['avg_speed'] > 0:
            stats['strouhal'] = (stats['avg_frequency'] * stats['avg_amplitude']) / stats['avg_speed']
        else:
            stats['strouhal'] = 0.0

        return stats

    def plot_time_series(self, ax=None):
        """Plot time series of key performance metrics"""
        if ax is None:
            fig, ax = plt.subplots(3, 1, figsize=(10, 10))

        # Swimming speed
        ax[0].plot(self.time, self.swimming_speed, 'b-', linewidth=1.5)
        ax[0].set_ylabel('Swimming Speed', fontsize=12)
        ax[0].grid(True, alpha=0.3)
        ax[0].set_title(f'Re={self.reynolds_number:.0f}, h/L={self.thickness_ratio:.3f}',
                       fontsize=14, fontweight='bold')

        # Thrust and Power
        ax1_twin = ax[1].twinx()
        ax[1].plot(self.time, self.thrust, 'r-', linewidth=1.5, label='Thrust')
        ax1_twin.plot(self.time, self.power, 'g-', linewidth=1.5, label='Power')
        ax[1].set_ylabel('Thrust', fontsize=12, color='r')
        ax1_twin.set_ylabel('Power', fontsize=12, color='g')
        ax[1].tick_params(axis='y', labelcolor='r')
        ax1_twin.tick_params(axis='y', labelcolor='g')
        ax[1].grid(True, alpha=0.3)

        # Efficiency
        ax[2].plot(self.time, self.efficiency, 'purple', linewidth=1.5)
        ax[2].set_ylabel('Propulsive Efficiency', fontsize=12)
        ax[2].set_xlabel('Time', fontsize=12)
        ax[2].grid(True, alpha=0.3)
        ax[2].axhline(y=np.mean(self.efficiency[int(0.2*len(self.time)):]),
                     color='k', linestyle='--', alpha=0.5, label='Mean')
        ax[2].legend()

        plt.tight_layout()
        return ax

    def plot_kinematics_adaptation(self, ax=None):
        """Plot adapted kinematics parameters"""
        if ax is None:
            fig, ax = plt.subplots(2, 1, figsize=(10, 6))

        # Adapted amplitude
        ax[0].plot(self.time, self.adapted_amplitude, 'b-', linewidth=1.5)
        ax[0].set_ylabel('Adapted Amplitude', fontsize=12)
        ax[0].grid(True, alpha=0.3)
        ax[0].set_title('Kinematic Adaptation', fontsize=14, fontweight='bold')

        # Adapted frequency
        ax[1].plot(self.time, self.adapted_frequency, 'r-', linewidth=1.5)
        ax[1].set_ylabel('Adapted Frequency', fontsize=12)
        ax[1].set_xlabel('Time', fontsize=12)
        ax[1].grid(True, alpha=0.3)

        plt.tight_layout()
        return ax


def compare_reynolds_effects(analyzers):
    """Compare performance across different Reynolds numbers"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Effects of Reynolds Number and Thickness on Foil Propulsion',
                 fontsize=16, fontweight='bold')

    stats_list = [a.compute_statistics() for a in analyzers]

    # Sort by Reynolds number
    stats_list.sort(key=lambda x: x['Re'])

    Re_values = [s['Re'] for s in stats_list]
    thickness_values = [s['h/L'] for s in stats_list]

    # Plot 1: Swimming speed vs Re
    colors = plt.cm.viridis(np.linspace(0, 1, len(stats_list)))
    for i, stats in enumerate(stats_list):
        axes[0, 0].scatter(stats['Re'], stats['avg_speed'],
                          s=200*stats['h/L'], c=[colors[i]],
                          label=f"h/L={stats['h/L']:.3f}", alpha=0.7, edgecolors='k')
    axes[0, 0].set_xlabel('Reynolds Number', fontsize=12)
    axes[0, 0].set_ylabel('Average Swimming Speed', fontsize=12)
    axes[0, 0].set_xscale('log')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=9)

    # Plot 2: Efficiency vs Re
    for i, stats in enumerate(stats_list):
        axes[0, 1].scatter(stats['Re'], stats['avg_efficiency'],
                          s=200*stats['h/L'], c=[colors[i]], alpha=0.7, edgecolors='k')
    axes[0, 1].set_xlabel('Reynolds Number', fontsize=12)
    axes[0, 1].set_ylabel('Average Propulsive Efficiency', fontsize=12)
    axes[0, 1].set_xscale('log')
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Strouhal number vs Re
    for i, stats in enumerate(stats_list):
        axes[1, 0].scatter(stats['Re'], stats['strouhal'],
                          s=200*stats['h/L'], c=[colors[i]], alpha=0.7, edgecolors='k')
    axes[1, 0].set_xlabel('Reynolds Number', fontsize=12)
    axes[1, 0].set_ylabel('Strouhal Number', fontsize=12)
    axes[1, 0].set_xscale('log')
    axes[1, 0].axhline(y=0.3, color='r', linestyle='--', alpha=0.5, label='Optimal St~0.3')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()

    # Plot 4: Power vs Speed (power curves)
    for i, stats in enumerate(stats_list):
        axes[1, 1].scatter(stats['avg_speed'], stats['avg_power'],
                          s=200*stats['h/L'], c=[colors[i]],
                          label=f"Re={stats['Re']:.0f}", alpha=0.7, edgecolors='k')
    axes[1, 1].set_xlabel('Swimming Speed', fontsize=12)
    axes[1, 1].set_ylabel('Power Required', fontsize=12)
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=9)

    plt.tight_layout()
    return fig


def generate_summary_table(analyzers):
    """Generate a summary table of all simulations"""
    print("\n" + "="*100)
    print("PERFORMANCE SUMMARY TABLE")
    print("="*100)
    print(f"{'Case':<25} {'Re':<10} {'h/L':<8} {'Speed':<10} {'Efficiency':<12} {'Strouhal':<10} {'Mode':<15}")
    print("-"*100)

    stats_list = [a.compute_statistics() for a in analyzers]
    stats_list.sort(key=lambda x: (x['Re'], x['h/L']))

    for i, stats in enumerate(stats_list):
        mode_str = "Anguilliform" if stats['mode'] < 0.5 else "Carangiform" if stats['mode'] > 0.7 else "Mixed"
        case_name = Path(analyzers[i].filename).stem
        print(f"{case_name:<25} {stats['Re']:<10.0f} {stats['h/L']:<8.3f} "
              f"{stats['avg_speed']:<10.4f} {stats['avg_efficiency']:<12.4f} "
              f"{stats['strouhal']:<10.4f} {mode_str:<15}")

    print("="*100 + "\n")


def main():
    """Main analysis routine"""
    # Find all performance data files
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        files = glob.glob('performance*.dat')

    if not files:
        print("No performance data files found!")
        print("Usage: python analyze_performance.py [performance_files...]")
        return

    print(f"Found {len(files)} performance data file(s)")

    # Create analyzers
    analyzers = []
    for f in files:
        try:
            analyzer = PerformanceAnalyzer(f)
            analyzers.append(analyzer)
            print(f"  - Loaded: {f} (Re={analyzer.reynolds_number}, h/L={analyzer.thickness_ratio})")
        except Exception as e:
            print(f"  - Error loading {f}: {e}")

    if not analyzers:
        print("No valid data files could be loaded!")
        return

    # Generate summary table
    generate_summary_table(analyzers)

    # Plot individual time series
    for analyzer in analyzers:
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))
        analyzer.plot_time_series(axes)
        output_name = f"timeseries_{Path(analyzer.filename).stem}.png"
        plt.savefig(output_name, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_name}")
        plt.close()

        # Plot kinematics adaptation
        fig, axes = plt.subplots(2, 1, figsize=(10, 6))
        analyzer.plot_kinematics_adaptation(axes)
        output_name = f"kinematics_{Path(analyzer.filename).stem}.png"
        plt.savefig(output_name, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_name}")
        plt.close()

    # Comparative plots
    if len(analyzers) > 1:
        fig = compare_reynolds_effects(analyzers)
        plt.savefig('reynolds_thickness_comparison.png', dpi=150, bbox_inches='tight')
        print("Saved: reynolds_thickness_comparison.png")
        plt.close()

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
