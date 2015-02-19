#ifndef included_DenseMobilityMatrix
#define included_DenseMobilityMatrix

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
  
enum DENSE_MATRIX_TYPE
{
    FILE=0,
    RPY,
    EMPIRICAL,
    HODLR,
    M_UNKNOWN_TYPE = -1
};

enum DENSE_INVERSE_METHOD
{
    LAPACK_CHOLESKY=0,
    LAPACK_LU,
    LAPACK_SVD,
    INV_UNKNOWN_TYPE = -1
};

class DenseMobilityMatrix
    :public SAMRAI::tbox::DescribedClass
{

/////////////////////////////// PUBLIC //////////////////////////////////////

public:
    /*!
     * \brief The only constructor of this class.
     */
    DenseMobilityMatrix(  
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	const int struct_ID);
  
    /*!
     * \brief Destructor for this class.
     */
    ~DenseMobilityMatrix();

    /*!
     * \brief Initialize DenseMobilityMatrix.
     */
    void 
    initializeDenseMobilityMatrix(
    const double Dt, 
    PetscScalar *X, 
    PetscScalar *W);

    /*!
     * \brief Deallocate storage for the mobility(friction) matrix
     *   
     */
    void
    deallocateDenseMobilityMatrix();

    /*!
     * \brief Computes the Mobility matrix
     */
    void 
    generateDenseMobilityMatrix();

    /*!
     * \brief Computes the resistance matrix by inverting the mobility matrix
     *
     */
    void 
    generateDenseFrictionMatrix();
    
    /*!
     * \brief Computes the solution on RHS.
     */
    void
    computeSolution(
    PetscScalar *rhs_b)

     /*!
     * \brief Returns id of the structore accosiated with current mobility matrix.
     */
    int
    getStructID();

private:
    /*!
     * \brief Scales mobility  matrix MM=scale_factor*[MM]+delta/weight
     *
     */
    void 
    scaleMobilityMatrix();

    /*!
     * \brief Reads stored Mobility matrix(es) from file(s).
     */
    void
    readMobilityMatrixfromFile(const double Dt);

   // Solver stuff
    std::string d_object_name;
    bool d_is_MobilityMatrixInitialized;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_input_db;
    SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator;

    //positions and weigths.
    Vec* d_X;
    Vec* d_W;		

    // Array storing friction matrix.
    double *d_friction_matrix;
    
    int d_struct_ID;           //structure number 
    int d_structure_no_nodes;  //number of nodes for current structure
    int d_total_no_nodes;      //total number of nodes of the system 
 
    // System physical parameters.
    double d_mu;       // Fluid viscosity  
    double d_rho;      // Fluid density  
    double d_DX;       // grid spacing
    double d_L_DOMAIN; //Length of domain used in 2D steady stokes empirical fitting

    
    //parameters used in this class.
    char* d_IBKernelName;
    double d_F_PERIODIC_CORRECTION;
    DENSE_INVERSE_METHOD d_INVERSE_METHOD;
    DENSE_MATRIX_TYPE d_MATRIX_TYPE;
    DENSE_MATRIX_TYPE d_TRUE_SUBMATRIX_TYPE;
    // Scaling factors for moblity matrix 
    // ScaledMM = d_mob_scale_1*[MM] + d_mob_scale_2*[I]
    double d_mob_scale_1, d_mob_scale_2; 

    // For LAPACK call 
    int* d_IPIV;
    //SVD method parameters
    double svd_inverse_tolerance, svd_inverse_epsilon;
};// DenseMobilityMatrix

}// IBAMR

#endif // #ifndef included_DenseMobilityMatrix

