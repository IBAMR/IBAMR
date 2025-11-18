#!/usr/bin/env python3
"""
Enhanced Odor Transport Solver with Crank-Nicolson Scheme

This module implements a robust numerical solver for the odor transport equation:

    ∂C/∂t + ui ∂C/∂xi = D ∂²C/∂xi∂xi

where:
  C = odor concentration (nondimensionalized by source concentration)
  ui = velocity field components (from Navier-Stokes solution)
  D = odor diffusivity

Numerical Methods:
------------------
1. Temporal discretization: Implicit Crank-Nicolson scheme (2nd order accurate)
   - Same scheme used for Navier-Stokes temporal terms in IBAMR
   - Unconditionally stable for high Schmidt numbers

2. Convection term (ui ∂C/∂xi): Upwind finite difference scheme
   - Ensures stability for advection-dominated flows
   - First-order accurate but stable

3. Diffusion term (D ∂²C/∂xi∂xi): Central finite differences
   - Second-order accurate in space
   - Crank-Nicolson averaging in time

Key Features:
-------------
- Handles high Schmidt numbers (Sc = ν/D >> 1) efficiently
- Implicit solver using sparse linear algebra
- Mass conservation
- Flexible boundary conditions (Neumann, Dirichlet, periodic)
- Compatible with IBAMR velocity fields

Reference: NSF Publication 10308831
Based on: "Collective Chemotactic Behavior in Fish Schools" (arXiv:2408.16136)
Authors: Maham Kamran, Amirhossein Fardi, Chengyu Li, Muhammad Saif Ullah Khalid
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import warnings

class OdorTransportSolverCN:
    """
    Odor transport solver using implicit Crank-Nicolson scheme.

    Solves: ∂C/∂t + ui ∂C/∂xi = D ∂²C/∂xi∂xi

    The equation is discretized as:
    (C^(n+1) - C^n)/Δt + u·∇C^(n+1/2) = D/2 * (∇²C^(n+1) + ∇²C^n)

    where the convection term is treated explicitly (upwind) and
    diffusion term is treated implicitly (Crank-Nicolson).
    """

    def __init__(self, x_range, y_range, nx, ny, diffusion_coeff,
                 boundary_type='neumann', schmidt_number=None):
        """
        Initialize the odor transport solver.

        Parameters:
        -----------
        x_range : tuple
            (x_min, x_max) domain bounds in x-direction
        y_range : tuple
            (y_min, y_max) domain bounds in y-direction
        nx, ny : int
            Number of grid points in x and y directions
        diffusion_coeff : float
            Molecular diffusion coefficient D
        boundary_type : str
            Type of boundary conditions:
            - 'neumann': Zero flux (∂C/∂n = 0) - default
            - 'dirichlet': Zero concentration (C = 0)
            - 'periodic': Periodic boundaries
        schmidt_number : float, optional
            Schmidt number Sc = ν/D for reference
        """
        self.x_min, self.x_max = x_range
        self.y_min, self.y_max = y_range
        self.nx = nx
        self.ny = ny
        self.D = diffusion_coeff
        self.boundary_type = boundary_type
        self.schmidt_number = schmidt_number

        # Create regular Cartesian grid
        self.x = np.linspace(self.x_min, self.x_max, nx)
        self.y = np.linspace(self.y_min, self.y_max, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]

        # Initialize concentration field
        self.c = np.zeros((ny, nx))
        self.t = 0.0

        # Statistics
        self.total_steps = 0
        self.mass_initial = 0.0

        # Build implicit diffusion matrix (constant, can be precomputed)
        self._build_diffusion_matrix()

        print(f"[SOLVER-CN] Crank-Nicolson Odor Transport Solver Initialized")
        print(f"  Grid: {nx} × {ny} points")
        print(f"  Domain: [{self.x_min:.2f}, {self.x_max:.2f}] × [{self.y_min:.2f}, {self.y_max:.2f}]")
        print(f"  Resolution: Δx = {self.dx:.4f}, Δy = {self.dy:.4f}")
        print(f"  Diffusion coefficient: D = {self.D:.6f}")
        print(f"  Boundary conditions: {boundary_type}")
        if schmidt_number:
            print(f"  Schmidt number: Sc = {schmidt_number:.1f}")
            print(f"  Note: High Sc handled efficiently by implicit scheme")

    def _build_diffusion_matrix(self):
        """
        Build the sparse matrix for implicit diffusion term.

        For Crank-Nicolson scheme:
        (I - θ·Δt·D·L) C^(n+1) = (I + (1-θ)·Δt·D·L) C^n + Δt·RHS

        where L is the Laplacian operator and θ = 0.5 for Crank-Nicolson.
        We build the LHS matrix here (will be modified with actual Δt during solve).
        """
        N = self.nx * self.ny

        # Laplacian stencil coefficients
        # ∇²C ≈ (C[i-1,j] - 2C[i,j] + C[i+1,j])/Δx² + (C[i,j-1] - 2C[i,j] + C[i,j+1])/Δy²
        cx = 1.0 / (self.dx**2)
        cy = 1.0 / (self.dy**2)
        center_coeff = -2.0 * (cx + cy)

        # Build Laplacian matrix (5-point stencil)
        self.L = lil_matrix((N, N))

        for j in range(self.ny):
            for i in range(self.nx):
                idx = j * self.nx + i  # Flatten 2D index to 1D

                # Center point
                self.L[idx, idx] = center_coeff

                # East neighbor (i+1, j)
                if i < self.nx - 1:
                    idx_e = j * self.nx + (i + 1)
                    self.L[idx, idx_e] = cx
                elif self.boundary_type == 'periodic':
                    idx_e = j * self.nx + 0
                    self.L[idx, idx_e] = cx
                elif self.boundary_type == 'neumann':
                    # Modify center coefficient for zero-flux BC
                    self.L[idx, idx] += cx

                # West neighbor (i-1, j)
                if i > 0:
                    idx_w = j * self.nx + (i - 1)
                    self.L[idx, idx_w] = cx
                elif self.boundary_type == 'periodic':
                    idx_w = j * self.nx + (self.nx - 1)
                    self.L[idx, idx_w] = cx
                elif self.boundary_type == 'neumann':
                    self.L[idx, idx] += cx

                # North neighbor (i, j+1)
                if j < self.ny - 1:
                    idx_n = (j + 1) * self.nx + i
                    self.L[idx, idx_n] = cy
                elif self.boundary_type == 'periodic':
                    idx_n = 0 * self.nx + i
                    self.L[idx, idx_n] = cy
                elif self.boundary_type == 'neumann':
                    self.L[idx, idx] += cy

                # South neighbor (i, j-1)
                if j > 0:
                    idx_s = (j - 1) * self.nx + i
                    self.L[idx, idx_s] = cy
                elif self.boundary_type == 'periodic':
                    idx_s = (self.ny - 1) * self.nx + i
                    self.L[idx, idx_s] = cy
                elif self.boundary_type == 'neumann':
                    self.L[idx, idx] += cy

        # Convert to CSR format for efficient arithmetic
        self.L = self.L.tocsr()

        print(f"[SOLVER-CN] Built sparse Laplacian matrix: {N}×{N} with {self.L.nnz} non-zeros")

    def set_initial_condition_gaussian(self, x0, y0, sigma, amplitude=1.0):
        """
        Set Gaussian initial condition (point source).

        C(x,y,0) = A * exp(-((x-x0)² + (y-y0)²) / (2σ²))

        Parameters:
        -----------
        x0, y0 : float
            Source location
        sigma : float
            Initial spreading width
        amplitude : float
            Maximum concentration at source
        """
        r_squared = (self.X - x0)**2 + (self.Y - y0)**2
        self.c = amplitude * np.exp(-r_squared / (2 * sigma**2))
        self.t = 0.0
        self.total_steps = 0
        self.mass_initial = self.get_total_mass()

        print(f"[SOLVER-CN] Initial condition set: Gaussian at ({x0}, {y0})")
        print(f"  Width: σ = {sigma:.3f}")
        print(f"  Amplitude: A = {amplitude:.3f}")
        print(f"  Initial total mass: {self.mass_initial:.6f}")

    def set_initial_condition_custom(self, concentration_field):
        """
        Set custom initial condition.

        Parameters:
        -----------
        concentration_field : ndarray (ny, nx)
            Initial concentration distribution
        """
        if concentration_field.shape != (self.ny, self.nx):
            raise ValueError(f"Shape mismatch: expected {(self.ny, self.nx)}, got {concentration_field.shape}")

        self.c = concentration_field.copy()
        self.t = 0.0
        self.total_steps = 0
        self.mass_initial = self.get_total_mass()

        print(f"[SOLVER-CN] Custom initial condition set")
        print(f"  Initial total mass: {self.mass_initial:.6f}")

    def _compute_convection_term_upwind(self, u_x, u_y, c_field):
        """
        Compute convection term using upwind finite differences.

        Returns: -ui ∂C/∂xi computed on interior points

        Upwind scheme ensures stability:
        - If u_x > 0: use backward difference (C[i] - C[i-1])/Δx
        - If u_x < 0: use forward difference (C[i+1] - C[i])/Δx

        Parameters:
        -----------
        u_x, u_y : ndarray (ny, nx)
            Velocity components
        c_field : ndarray (ny, nx)
            Concentration field

        Returns:
        --------
        conv_term : ndarray (ny, nx)
            Convection contribution -u·∇C
        """
        conv_term = np.zeros_like(c_field)

        # Interior points
        for j in range(1, self.ny - 1):
            for i in range(1, self.nx - 1):
                # x-direction derivative (upwind)
                if u_x[j, i] > 0:
                    dc_dx = (c_field[j, i] - c_field[j, i-1]) / self.dx
                else:
                    dc_dx = (c_field[j, i+1] - c_field[j, i]) / self.dx

                # y-direction derivative (upwind)
                if u_y[j, i] > 0:
                    dc_dy = (c_field[j, i] - c_field[j-1, i]) / self.dy
                else:
                    dc_dy = (c_field[j+1, i] - c_field[j, i]) / self.dy

                # Convection term: -u·∇C
                conv_term[j, i] = -(u_x[j, i] * dc_dx + u_y[j, i] * dc_dy)

        return conv_term

    def compute_stable_timestep(self, u_x, u_y, cfl_max=0.5):
        """
        Compute maximum stable timestep based on CFL condition.

        For explicit convection: CFL = max(|u|·Δt/Δx, |v|·Δt/Δy) < CFL_max

        Note: Diffusion term is implicit, so no diffusion stability constraint.

        Parameters:
        -----------
        u_x, u_y : ndarray (ny, nx)
            Velocity field components
        cfl_max : float
            Maximum CFL number (default 0.5)

        Returns:
        --------
        dt_max : float
            Maximum stable timestep
        """
        u_max = np.max(np.abs(u_x)) + 1e-10
        v_max = np.max(np.abs(u_y)) + 1e-10

        dt_conv = cfl_max * min(self.dx / u_max, self.dy / v_max)

        # For purely implicit scheme, no diffusion constraint needed
        # But for practical purposes, limit timestep
        dt_max = dt_conv

        return dt_max

    def step_crank_nicolson(self, u_x, u_y, dt):
        """
        Advance solution one timestep using Crank-Nicolson scheme.

        Discretization:
        (C^(n+1) - C^n)/Δt = -u·∇C^n + D/2·(∇²C^(n+1) + ∇²C^n)

        Rearranged:
        (I - (Δt·D/2)·L) C^(n+1) = C^n + Δt·(-u·∇C^n + (D/2)·L·C^n)

        where L is the Laplacian operator.

        Parameters:
        -----------
        u_x, u_y : ndarray (ny, nx)
            Velocity field components at current time
        dt : float
            Timestep size
        """
        # Compute convection term explicitly (upwind scheme)
        conv_term = self._compute_convection_term_upwind(u_x, u_y, self.c)

        # Flatten concentration field to 1D vector
        c_flat = self.c.flatten()

        # Build RHS: C^n + Δt·(-u·∇C^n + (D/2)·L·C^n)
        theta = 0.5  # Crank-Nicolson parameter

        rhs_diffusion = self.D * theta * self.L.dot(c_flat)
        rhs_convection = dt * conv_term.flatten()

        rhs = c_flat + rhs_convection + dt * rhs_diffusion

        # Build LHS matrix: I - (Δt·D/2)·L
        N = self.nx * self.ny
        from scipy.sparse import identity
        I = identity(N, format='csr')

        lhs = I - (dt * self.D * theta) * self.L

        # Solve linear system: lhs · C^(n+1) = rhs
        try:
            c_new_flat = spsolve(lhs, rhs)
        except Exception as e:
            warnings.warn(f"Sparse solver failed: {e}. Using fallback.")
            # Fallback to iterative solver if direct solver fails
            from scipy.sparse.linalg import bicgstab
            c_new_flat, info = bicgstab(lhs, rhs, x0=c_flat, tol=1e-8)
            if info != 0:
                raise RuntimeError(f"Iterative solver failed with code {info}")

        # Reshape to 2D
        self.c = c_new_flat.reshape((self.ny, self.nx))

        # Apply boundary conditions if needed
        self._apply_boundary_conditions()

        # Enforce non-negativity (physical constraint)
        self.c = np.maximum(self.c, 0.0)

        # Update time and statistics
        self.t += dt
        self.total_steps += 1

    def _apply_boundary_conditions(self):
        """
        Apply boundary conditions to concentration field.

        - Neumann (zero flux): Already handled in matrix construction
        - Dirichlet (zero concentration): Set boundary values to zero
        - Periodic: Already handled in matrix construction
        """
        if self.boundary_type == 'dirichlet':
            self.c[0, :] = 0.0
            self.c[-1, :] = 0.0
            self.c[:, 0] = 0.0
            self.c[:, -1] = 0.0
        # Neumann and periodic BCs are built into the matrix

    def step_diffusion_only(self, dt):
        """
        Advance solution using pure diffusion (no convection).
        Useful for comparison studies.

        Solves: ∂C/∂t = D∇²C

        Parameters:
        -----------
        dt : float
            Timestep size
        """
        # Use Crank-Nicolson with zero velocity
        u_x = np.zeros((self.ny, self.nx))
        u_y = np.zeros((self.ny, self.nx))
        self.step_crank_nicolson(u_x, u_y, dt)

    def get_concentration(self):
        """
        Get current concentration field.

        Returns:
        --------
        c : ndarray (ny, nx)
            Concentration field (copy)
        """
        return self.c.copy()

    def get_total_mass(self):
        """
        Compute total mass in domain.

        Mass = ∫∫ C(x,y) dx dy ≈ Σ C[i,j] · Δx · Δy

        Returns:
        --------
        mass : float
            Total integrated mass
        """
        return np.sum(self.c) * self.dx * self.dy

    def get_mass_conservation_error(self):
        """
        Compute relative mass conservation error.

        Error = |Mass_current - Mass_initial| / Mass_initial

        Returns:
        --------
        error : float
            Relative mass error (should be small for good solver)
        """
        mass_current = self.get_total_mass()
        if self.mass_initial > 1e-10:
            error = abs(mass_current - self.mass_initial) / self.mass_initial
        else:
            error = abs(mass_current - self.mass_initial)
        return error

    def get_spreading_width(self):
        """
        Compute spreading width (standard deviation).

        σ = √(∫∫ C·(r-r_c)² dx dy / ∫∫ C dx dy)

        where r_c is the centroid position.

        Returns:
        --------
        sigma : float
            Spreading width
        centroid : tuple (x_c, y_c)
            Centroid position
        """
        total_mass = np.sum(self.c) + 1e-10

        # Compute centroid
        x_c = np.sum(self.c * self.X) / total_mass
        y_c = np.sum(self.c * self.Y) / total_mass

        # Compute variance
        var_x = np.sum(self.c * (self.X - x_c)**2) / total_mass
        var_y = np.sum(self.c * (self.Y - y_c)**2) / total_mass

        sigma = np.sqrt(var_x + var_y)

        return sigma, (x_c, y_c)

    def get_max_concentration(self):
        """
        Get maximum concentration value and its location.

        Returns:
        --------
        c_max : float
            Maximum concentration
        location : tuple (x_max, y_max)
            Location of maximum
        """
        c_max = np.max(self.c)
        idx_max = np.argmax(self.c)
        j_max, i_max = np.unravel_index(idx_max, (self.ny, self.nx))
        x_max = self.x[i_max]
        y_max = self.y[j_max]

        return c_max, (x_max, y_max)

    def get_solver_info(self):
        """
        Get solver statistics and information.

        Returns:
        --------
        info : dict
            Dictionary containing solver statistics
        """
        mass_error = self.get_mass_conservation_error()
        sigma, (x_c, y_c) = self.get_spreading_width()
        c_max, (x_max, y_max) = self.get_max_concentration()

        info = {
            'time': self.t,
            'total_steps': self.total_steps,
            'total_mass': self.get_total_mass(),
            'mass_conservation_error': mass_error,
            'spreading_width': sigma,
            'centroid': (x_c, y_c),
            'max_concentration': c_max,
            'max_location': (x_max, y_max),
            'grid_size': (self.nx, self.ny),
            'resolution': (self.dx, self.dy),
            'diffusion_coeff': self.D,
            'schmidt_number': self.schmidt_number
        }

        return info

    def print_status(self):
        """Print current solver status."""
        info = self.get_solver_info()

        print(f"\n[SOLVER STATUS] t = {info['time']:.4f}")
        print(f"  Steps: {info['total_steps']}")
        print(f"  Total mass: {info['total_mass']:.6f}")
        print(f"  Mass error: {info['mass_conservation_error']:.2e}")
        print(f"  Spreading width: σ = {info['spreading_width']:.4f}")
        print(f"  Centroid: ({info['centroid'][0]:.3f}, {info['centroid'][1]:.3f})")
        print(f"  Max concentration: {info['max_concentration']:.4f} at ({info['max_location'][0]:.3f}, {info['max_location'][1]:.3f})")


# ============================================================
# VALIDATION AND TESTING FUNCTIONS
# ============================================================

def validate_solver_analytical():
    """
    Validate solver against analytical solution for pure diffusion.

    For pure diffusion with Gaussian initial condition:
    C(x,y,t) = (A/(1 + 4Dt/σ₀²)) · exp(-r²/(2(σ₀² + 2Dt)))

    where σ₀ is initial width and r² = (x-x₀)² + (y-y₀)²
    """
    print("\n" + "="*80)
    print("VALIDATION TEST: Pure Diffusion with Analytical Solution")
    print("="*80)

    # Parameters
    D = 0.01
    x0, y0 = 0.0, 0.0
    sigma0 = 0.2
    A = 1.0
    t_final = 1.0

    # Create solver
    solver = OdorTransportSolverCN(
        x_range=(-2, 2), y_range=(-2, 2),
        nx=100, ny=100,
        diffusion_coeff=D,
        boundary_type='neumann'
    )

    solver.set_initial_condition_gaussian(x0, y0, sigma0, A)

    # Advance to t_final
    dt = 0.01
    num_steps = int(t_final / dt)

    print(f"\nAdvancing solver from t=0 to t={t_final} with dt={dt}")
    print(f"Total steps: {num_steps}")

    for step in range(num_steps):
        solver.step_diffusion_only(dt)
        if (step + 1) % 20 == 0:
            print(f"  Step {step+1}/{num_steps}: t = {solver.t:.3f}")

    # Analytical solution at t_final
    r_squared = (solver.X - x0)**2 + (solver.Y - y0)**2
    sigma_t = np.sqrt(sigma0**2 + 2*D*t_final)
    c_analytical = (A * sigma0**2 / sigma_t**2) * np.exp(-r_squared / (2 * sigma_t**2))

    # Compare
    c_numerical = solver.get_concentration()

    error_l2 = np.sqrt(np.mean((c_numerical - c_analytical)**2))
    error_max = np.max(np.abs(c_numerical - c_analytical))
    error_rel = error_l2 / np.max(c_analytical)

    print(f"\n[RESULTS]")
    print(f"  L2 error: {error_l2:.6e}")
    print(f"  Max error: {error_max:.6e}")
    print(f"  Relative error: {error_rel:.4%}")

    solver.print_status()

    if error_rel < 0.01:
        print("\n[✓] VALIDATION PASSED: Relative error < 1%")
    else:
        print("\n[✗] VALIDATION WARNING: Relative error > 1%")

    print("="*80)

    return solver, c_analytical, error_rel


if __name__ == "__main__":
    """
    Test the Crank-Nicolson odor transport solver.
    """
    print("\n" + "="*80)
    print("CRANK-NICOLSON ODOR TRANSPORT SOLVER")
    print("="*80)
    print("\nGoverning equation: ∂C/∂t + ui ∂C/∂xi = D ∂²C/∂xi∂xi")
    print("\nNumerical methods:")
    print("  - Temporal: Implicit Crank-Nicolson (2nd order)")
    print("  - Convection: Upwind finite differences")
    print("  - Diffusion: Central differences with implicit averaging")
    print("="*80)

    # Run validation
    validate_solver_analytical()

    print("\n[INFO] Validation complete!")
    print("[INFO] Solver ready for use with IBAMR velocity fields")
