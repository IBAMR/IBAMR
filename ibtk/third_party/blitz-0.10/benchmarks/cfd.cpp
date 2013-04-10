#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

BZ_USING_NAMESPACE(blitz)

/*
 * The current implementation of stencil objects forces these variables
 * to be placed in global scope.  Ugh.  This restriction will be removed
 * eventually.
 */
double rho;                        // Density of fluid
double recip_rho;                  // 1/rho
double eta;                        // Kinematic viscosity
double time_now;                   // Elapsed seconds
double delta_t;                    // Time step
double volume;                     // Volume of a cell
double airPressure;                // Air pressure (Pa)
double spatialStep;                // Grid element size
double gravity;                    // Acceleration due to gravity
double gravityPressureGradient;    // Pressure gradient due to gravity
/*
 * The "geometry" object specifies how an array is mapped into real-world
 * space.  In this case, "UniformCubicGeometry" is used, which means that
 * the real-world grid is orthogonal, regularly spaced, with the same spatial
 * step in each dimension.
 */

UniformCubicGeometry<3> geom;      // Geometry
/*
 * Some typedefs to make life easier.
 */

typedef TinyVector<double,3> vector3d;
typedef Array<vector3d,3> vectorField;
typedef Array<double,3>   scalarField;

/***********          Timestep the velocity field           ************
 * This is a 63-point stencil.  For example, Laplacian3DVec4 turns into
 * a 45-point stencil: each 2nd derivative is a 5-point stencil, and
 * there are 9 of these derivatives to take the Laplacian of a 3D vector
 * field.
 */

BZ_DECLARE_STENCIL5(timestep, V, nextV, P, advect, force)

    nextV = *V + delta_t * ( recip_rho * (
      eta * Laplacian3DVec4(V,geom) - grad3D4(P, geom) + *force) - *advect);

BZ_END_STENCIL
/*
 * Allocate arrays and set their initial state
 */
void setup(const int N, vectorField& V, vectorField& nextV, scalarField& P,
    scalarField& P_rhs, vectorField& advect, vectorField& force)
{
    // A 1m x 1m x 1m domain
    spatialStep = 1.0 / (N - 1);
    geom = UniformCubicGeometry<3>(spatialStep);

    // Allocate arrays
    allocateArrays(shape(N,N,N), advect, V, nextV, force);  // vector fields
    allocateArrays(shape(N,N,N), P, P_rhs);                 // scalar fields

    // Since incompressibility is assumed, pressure only shows up as
    // derivative terms in the equations.  We choose airPressure = 0
    // as an arbitrary datum.

    airPressure = 0;             // Pa
    rho = 1000;                  // density of fluid, kg/m^3
    recip_rho = 1.0 / rho;       // inverse of density
    eta = 1.0e-6;                // kinematic viscosity of fluid, m^2/s
    gravity = 9.81;              // m/s^2
    delta_t = 0.001;             // initial time step, in seconds
    volume = pow3(spatialStep);  // cubic volume associated with grid point

    // Kludge: Set eta high, so that the flow will spread faster.
    // This means the cube is filled with molasses, rather than water.
    eta *= 1000;

    // Initial conditions: quiescent
    V = 0.0;
    P_rhs = 0.0;
    advect = 0.0;
    nextV = 0.0;
    P = 0.0;
    force = 0.0;
}

// Calculate a simple check on a vector field
void record(vectorField& V)
{
    // Calculate the magnitude of a field
    const int x=0, y=1, z=2;
    double magx = sum(pow2(V[x])) / V.numElements();
    double magy = sum(pow2(V[y])) / V.numElements();
    double magz = sum(pow2(V[z])) / V.numElements();

    cout << "norm = [" << magx
        << " " << magy << " " << magz << " ]" << endl;
}

void iterate(vectorField& V, vectorField& nextV, scalarField& P,
    scalarField& P_rhs, vectorField& advect, vectorField& force)
{
    // Time step
    applyStencil(timestep(), V, nextV, P, advect, force);
}

int main()
{
    vectorField V, nextV;        // Velocity fields
    scalarField P, P_rhs;        // Pressure fields
    vectorField advect;          // Advection field
    vectorField force;           // Forcing function

    const int N = 50;            // Arrays are NxNxN

    setup(N, V, nextV, P, P_rhs, advect, force);

    const int nIters = 10;

    for (int i=0; i < nIters; ++i)
    {
        iterate(V, nextV, P, P_rhs, advect, force);
    }

    return 0;
}

