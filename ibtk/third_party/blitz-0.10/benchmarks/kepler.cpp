#include <blitz/vector.h>
#include <blitz/tinyvec.h>

BZ_USING_NAMESPACE(blitz)

/*
 * Test a 12th order symmetric multistep method for solving the equations
 * of motion of a single planet circling the Sun.  The Sun is fixed in
 * space.
 *
 * Original F77 version written by John K. Prentice, Quetzal Computational
 * Associates, 21 Decmber 1992
 * 
 * Blitz++ version by Todd Veldhuizen, 17 August 1997
 * The C++ version is a faithful translation of the Fortran 90 version,
 * so apologies for the "C++Tran" style.
 */

inline double relativeError(double a, double b)
{
    if (b != 0.0)
        return (a - b) / b;
    else
        return a;
}

int main()
{
    Vector<double> x_position_numerical(13), y_position_numerical(13),
        alpha(13), beta(13), gamma(13), x_acceleration(13), y_acceleration(13);
   
    /*
     * 12th order symmetric method coefficients
     *
     * Reference: "Symmetric Multistep Methods for the Numerical
     * Integration of Planetary Orbits", G. D. Quinlan and
     * S. Tremaine, The Astronomical Journal, 100 (1990), page 1695.
     *
     * Note!! The beta below are actually 53,222,400 times the
     * real beta.  This common factor is divided out in the
     * symmetric multistep calculation itself, in order to minimize
     * round-off
     */

    const double beta_factor = 53222400.0;
    alpha = 1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, -2.0, 1.0;
    beta  = 0.0, 90987349.0, -229596838.0, 812627169.0, -1628539944.0, 
            2714971338.0, -3041896548.0, 2714971338.0, -1628539944.0,
            812627169.0, -229596838.0, 90987349.0, 0.0;

    /*
     * 12th order Cowell predictor coefficients
     *
     * Reference:  "Astronomical Papers Prepared for the Use of the
     * American Ephemeris and Nautical Almanac", C. J. Cohen, E. C.
     * Hubbard, and C. Oesterwinter, 22 (1973), page 20-21.
     *
     * Note!!  The gamma below are actually 1,743,565,824,000 times
     * the real gamma.  This common factor is divided out in the
     * Cowell predictor calculation itself, in order to minimize
     * round-off
     */

    const double gamma_factor = 1743565824000.0;
    gamma = 9072652009253.0, -39726106418680.0, 140544566352762.0, 
       -344579280210129.0, 613137294629235.0, -811345852376496.0,
       807012356281740.0, -602852367932304.0, 333888089374395.0, 
       -133228219027160.0, 36262456774618.0, -6033724094760.0,
       463483373517.0;

    // Initialize variables

    const double time_step = 0.25,
                 stop_time = 365000.0,
                 radius    = 1.0;
    double time = - time_step;

    cout << " Position solution via 12th order symmetric multistep method\n"
         << " Velocity solution via 12th order Cowell predictor method\n"
         << "     radius = " << radius << ", time step = " << time_step
         << endl;

    // Define a constant which is needed later by the exact solution
    const double gaussian_constant = 0.01720209895;
    const double gravitational_constant = pow(gaussian_constant,2);
    const double constant = sqrt(gravitational_constant/pow(radius,3));

    // Initialize the first 12 numerical values using the exact values

    double x_position_exact, y_position_exact;

    for (int j=-1; j <= 11; ++j)
    {
        if (j >= 0)
            time += time_step;

        x_position_exact = radius * cos(constant * time);
        y_position_exact = radius * sin(constant * time);

        if (j >= 0)
        {
            x_position_numerical(j) = x_position_exact;
            y_position_numerical(j) = y_position_exact;
        }

        x_acceleration(j+1) = -gravitational_constant/pow(radius,3) 
            * x_position_exact;
        y_acceleration(j+1) = -gravitational_constant/pow(radius,3)
            * y_position_exact;
    }

    /*
     * Compute exact kinetic and potential energies, and the
     * angular momentum.  These values are all divided by the mass
     * of the object.  Since they are conserved, they will never change
     * and hence do not have to be recalculated later.
     */

    double x_dot_exact = -radius * constant * sin(constant*time),
         y_dot_exact =  radius * constant * cos(constant*time),
         exact_velocity_squared = pow(x_dot_exact,2) + pow(y_dot_exact,2),
         exact_kinetic_energy = 0.5 * exact_velocity_squared,
         exact_potential_energy = -gravitational_constant / radius,
         exact_total_energy = exact_potential_energy + exact_kinetic_energy,
         exact_angular_momentum = x_position_exact * y_dot_exact
             - y_position_exact * x_dot_exact;

    double x_dot_numerical, y_dot_numerical;

    // Perform loop over time

    while (time <= stop_time)                           
    {
        // Advance time step (eek!)
        time += time_step;        

        // Calculate new acceleration of body at time=time-time_step
        double numerical_radius_squared = pow(x_position_numerical(11),2)
            + pow(y_position_numerical(11),2);
        x_acceleration(12) = -gravitational_constant
            / pow(numerical_radius_squared, 1.5) * x_position_numerical(11);
        y_acceleration(12) = -gravitational_constant
            / pow(numerical_radius_squared, 1.5) * y_position_numerical(11);

        // Numerically solve for the new positions using a 12th order
        // symmetric multistep method.

        // First sum the first and second terms

        double x_alpha_sum = dot(alpha(Range(0,11)), 
            x_position_numerical(Range(0,11)));
        double y_alpha_sum = dot(alpha(Range(0,11)), 
            y_position_numerical(Range(0,11)));

        double x_beta_sum = dot(beta(Range(0,11)), x_acceleration(Range(1,12)));
        double y_beta_sum = dot(beta(Range(0,11)), y_acceleration(Range(1,12)));
        x_position_numerical(12) = (-x_alpha_sum) + pow(time_step,2) 
            * (x_beta_sum / beta_factor);
        y_position_numerical(12) = (-y_alpha_sum) + pow(time_step,2)
            * (y_beta_sum / beta_factor);

        // Numerically solve for the new velocities using a 12th order
        // Cowell predictor method.

        // First sum the gamma terms

        double x_gamma_sum = dot(gamma, x_acceleration.reverse()),
               y_gamma_sum = dot(gamma, y_acceleration.reverse());

        x_dot_numerical = (x_position_numerical(11)
            - x_position_numerical(10)) / time_step + time_step 
            * (x_gamma_sum / gamma_factor);
        y_dot_numerical = (y_position_numerical(11)
            - y_position_numerical(10)) / time_step + time_step
            * (y_gamma_sum / gamma_factor);

        // Push the stack down one

        for (int j=0; j <= 11; ++j)
        {
            x_position_numerical(j) = x_position_numerical(j+1);
            y_position_numerical(j) = y_position_numerical(j+1);
            x_acceleration(j) = x_acceleration(j+1);
            y_acceleration(j) = y_acceleration(j+1);
        }
    }

    // Print results

    // First compute energies and angular momenta (add divided by the mass
    // of the object)

    double numerical_velocity_squared = pow(x_dot_numerical,2) +
               pow(y_dot_numerical,2),
           numerical_radius = sqrt(pow(x_position_numerical(12),2)
               + pow(y_position_numerical(12),2)),
           numerical_kinetic_energy = 0.5 * numerical_velocity_squared,
           numerical_potential_energy = -gravitational_constant 
               / numerical_radius,
           numerical_total_energy = numerical_potential_energy
               + numerical_kinetic_energy,
           numerical_angular_momentum = x_position_numerical(12)
               * y_dot_numerical - y_position_numerical(12) * x_dot_numerical;

    // Compute exact results for comparison to the numerical results

   x_position_exact = radius * cos(constant * time);
   y_position_exact = radius * sin(constant * time);
   x_dot_exact = -radius * constant * sin(constant * time);
   y_dot_exact =  radius * constant * cos(constant * time);

    // Next compute relative errors

    double radius_error = relativeError(numerical_radius, radius),
           x_error = relativeError(x_position_numerical(12), x_position_exact),
           y_error = relativeError(y_position_numerical(12), y_position_exact),
           x_dot_error = relativeError(x_dot_numerical, x_dot_exact),
           y_dot_error = relativeError(y_dot_numerical, y_dot_exact);

    double kinetic_energy_error = relativeError(numerical_kinetic_energy,
               exact_kinetic_energy),
           potential_energy_error = relativeError(numerical_potential_energy,
               exact_potential_energy),
           total_energy_error = relativeError(numerical_total_energy,
               exact_total_energy),
           angular_momentum_error = relativeError(numerical_angular_momentum,
               exact_angular_momentum);

    cout << " Time = " << time << endl
         << "    x rel error  = " << x_error << " y rel error  = " << y_error
         << endl
         << "    vx rel error = " << x_dot_error << " vy rel error = " 
         << y_dot_error << endl
         << "    KE rel error = " << kinetic_energy_error 
         << " PE rel error = " << potential_energy_error << endl
         << "    TE rel error = " << total_energy_error << " AM rel error = "
         << angular_momentum_error << endl
         << "    numerical radius = " << numerical_radius 
         << " radius rel error = " << radius_error << endl;

    return 0;
}

