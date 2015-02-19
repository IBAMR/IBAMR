/////////////////////////////// INCLUDES /////////////////////////////////////

#include <zlib.h>
#include <math.h>
#include <string>
#include <algorithm>

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h"
#include "petsc-private/petscimpl.h"
#include "DenseMobilityMatrix.h"

namespace IBAMR
{

#ifndef KRON
#define KRON(i, j) ((i==j)?1:0) //Kronecker symbol
#endif

namespace 
{
    
//**********************************************
// LAPACK functions declaration for dense matrix direct solver 
//**********************************************
// LAPACK function for LU factorization
extern "C" int 
dgetrf_(
    const int& N1, 
    const int& N2, 
    double* A, 
    const int& LDA, 
    int *IPIV,
    int& INFO);

// LAPACK function for that uses LU factorization to find solution
extern "C" int 
dgetrs_(
    const char* TRANS, 
    const int& N, 
    const int& NRHS, 
    const double* A, 
    const int& LDA, 
    const int * IPIV,
    double *B, 
    const int& LDB, 
    int& INFO);  

// LAPACK function for Cholesky factorization
extern "C" int 
dpotrf_(
    const char* UPLO, 
    const int& N, 
    double* A, 
    const int& LDA, 
    int& INFO);

// LAPACK function for that uses Cholesky factorization to find solution
extern "C" int 
dpotrs_(
    const char* UPLO, 
    const int& N, 
    const int& NRHS, 
    const double* A, 
    const int& LDA, 
    double *B, 
    const int& LDB, 
    int& INFO);  
//LAPACK function for SVD factorization
extern "C" void 
dsyevr_(
    const char* jobz, 
    const char* range, 
    const char* uplo, 
    const int& n, 
    double* a,
    const int& lda, 
    const double& vl, 
    const double& vu, 
    const int& il, 
    const int& iu, 
    const double& abstol,
    int& m, 
    double* w, 
    double* z, 
    const int& ldz, 
    int* isuppz, 
    double* work,
    const int& lwork, 
    int* iwork, 
    const int& liwork, 
    int& info );

//**********************************************
// End of LAPACK functions declaration 


//**********************************************
// MobilityFunctions declaration  
//**********************************************

//Empirical formula f(r) and g(r) coefficients
extern "C" void
getEmpiricalMobilityComponents(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double rr,
    const double DX,
    const bool resetAllConstants,
    const double L_domain,
    double *F_MobilityValue,
    double *G_Mobilityvalue);

//Empirical formula Mobility matrix generator
extern "C" void
getEmpiricalMobilityMatrix(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double DX,
    const double *X, 
    const int N,
    const bool resetAllConstants,
    const double PERIODIC_CORRECTION,
    const double L_domain,
    double *MM);

//RPY Mobility matrix generator
extern "C" void
getRPYMobilityMatrix(
    const char* IBKernelName,
    const double MU,
    const double DX,
    const double *X, 
    const int N,
    const double PERIODIC_CORRECTION,
    double *MM);

//Hydrodynamic Radius value
extern "C" double
getHydroRadius(const char* IBKernelName);

/*!
 * returns squared norm of the vector
 */
inline double
compute_sqnorm(
    const double *a_vec)
{
#if (NDIM==3)
    return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1]+a_vec[2]*a_vec[2];
#elif(NDIM==2)
  return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1];
#endif
}
//**********************************************
// End of MobilityFunctions functions declaration 
}// anonymous

DenseMobilityMatrix::DenseMobilityMatrix(  
    Pointer<Database> input_db,
    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
    const int struct_ID):
    d_input_db(input_db),
    d_ins_integrator(navier_stokes_integrator),
    d_friction_matrix(NULL),
    d_struct_ID(struct_ID),
    d_mu(1.0),
    d_rho(1.0),
    d_DX(0.0),
    d_L_DOMAIN(0.0),
    d_mob_scale_1(1.0), 
    d_mob_scale_2(0.0), 
    d_IBKernelName(NULL),
    d_F_PERIODIC_CORRECTION(0.0),
    d_INVERSE_METHOD(INV_UNKNOWN_TYPE),
    d_MATRIX_TYPE(M_UNKNOWN_TYPE),
    d_TRUE_SUBMATRIX_TYPE(M_UNKNOWN_TYPE),
    d_IPIV(NULL),
    d_structure_no_nodes(0),
    d_total_no_nodes(0) 

{
    // Get from input
    const std::string matrix_type =  input_db->getString("dense_mobility_type");
    if      (matrix_type == "FILE")            d_MATRIX_TYPE = FILE;
    else if (matrix_type == "RPY")             d_MATRIX_TYPE = RPY;
    else if (matrix_type == "EMPIRICAL")       d_MATRIX_TYPE = EMPIRICAL;
    else if (matrix_type == "HODLR")           d_MATRIX_TYPE = HODLR;
    else    TBOX_ERROR("DenseMobilityMatrix::DenseMobilityMatrix:  Unknown mobility matrix type is provided." << std::endl);

    if (d_MATRIX_TYPE == HODLR)
    {
	svd_inverse_tolerance = input_db->getDatabase("HODLR")->getDouble("hodlr_tolerance");
    }
    else
    {        
	const std::string inverse_name =  input_db->getString("inverse_method");
	if      (inverse_name == "LAPACK_CHOLESKY") d_INVERSE_METHOD = LAPACK_CHOLESKY;
	else if (inverse_name == "LAPACK_LU")       d_INVERSE_METHOD = LAPACK_LU;
	else if (inverse_name == "LAPACK_SVD")    
	{
	    d_INVERSE_METHOD    = LAPACK_SVD;
	    svd_inverse_tolerance = input_db->getDatabase("LAPACK_SVD")->getDouble("eigenvalue_replace_value");
	    svd_inverse_epsilon   = input_db->getDatabase("LAPACK_SVD")->getDouble("min_eigenvalue_threshold");
	}
	
	
	
	// Check if moblity matrix is provided in the file
	if (d_MATRIX_TYPE == FILE)
	{
	    tbox::Array<std::string> files = input_db->getDatabase("FILE")->getStringArray("mob_matrix_filenames");
	    if (files.size() < 1) TBOX_ERROR("DenseMobilityMatrix::DenseMobilityMatrix  File name for mobility matrix is missing." << std::endl);
	    else
	    {
		for (int cnt=0;cnt<files.size();cnt++) d_mobfile_names.push_back(files[cnt]);
	    }
	    if (files.size()>1)
	    {
		const std::string submatrix_type = input_db->getDatabase("FILE")->getString("mobility_type_subblock");
		if      (submatrix_type == "RPY")             d_TRUE_SUBMATRIX_TYPE = RPY;
		else if (submatrix_type == "EMPIRICAL")       d_TRUE_SUBMATRIX_TYPE = EMPIRICAL;
		else    TBOX_ERROR("DenseMobilityMatrix::DenseMobilityMatrix: Unknown mobility type for true moblity submatrix is provided." << std::endl);
	    }
	} else 
	{
	    // Get the correction term in f function due to periodic BC
	    d_F_PERIODIC_CORRECTION = input_db->getDouble("f_periodic_correction");
	    //scalling(regularization) factors
	    d_mob_scale_1 = input_db->getDoubleWithDefault("mob_scale_1",1.0);
	    d_mob_scale_2 = input_db->getDoubleWithDefault("mob_scale_2",0.0);
	}//end else file
    } //end else hodlr

    return;
}
  
DenseMobilityMatrix::~DenseMobilityMatrix()
{
    if (d_friction_matrix) delete [] d_friction_matrix;
    if (d_IPIV)            delete [] d_IPIV;
    if (d_X)               delete [] d_X;
    if (d_W)               delete [] d_W;
    d_is_MobilityMatrixInitialized=false;
}


void
DenseMobilityMatrix::deallocateDenseMobilityMatrix()
{
    if (d_friction_matrix) delete [] d_friction_matrix;
    if (d_IPIV)            delete [] d_IPIV;
    if (d_X)               delete [] d_X;
    if (d_W)               delete [] d_W;
    d_is_MobilityMatrixInitialized=false;

    return;
}// deallocateSolverState

void 
DenseMobilityMatrix::initializeDenseMobilityMatrix(
    const double Dt, 
    PetscScalar *X, 
    PetscScalar *W)
{
    d_rho = d_ins_integrator->getStokesSpecifications()->getRho(); 
    d_mu  = d_ins_integrator->getStokesSpecifications()->getMu();

    const int size = d_structure_no_nodes*NDIM;
    d_X=new double[size];
    d_W=new double[size];

    for (int i = 0; i < size; i++) 
    {
	d_X[i] = X[i];
	d_W[i] = W[i];
    }

    Pointer<PatchLevel<NDIM> > struct_patch_level        = d_ins_integrator->getPatchHierarchy()->getPatchLevel(getFinestLevelNumber());
    const IntVector<NDIM>& ratio                         = struct_patch_level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom      = d_ins_integrator->getPatchHierarchy()->getGridGeometry();
    const double* dx                                     = grid_geom->getDx();
    d_DX=  (dx[0]/ratio(0));

    //set kernel name
    d_IBKernelName=(char*) d_cib_method->getLDataManager()->getDefaultInterpKernelFunction().c_str();
    
#if (NDIM==2)
    d_L_DOMAIN = input_db->getDoubleWithDefault("Domain_length",0.0);
    if (!d_L_DOMAIN)
    {
	for(int d = 0; d < NDIM; ++d)  L_DOMAIN = MAX(L_DOMAIN, *(grid_geom->getXUpper()+d)-*(grid_geom->getXLower()+d));
    }
#endif

    //set to initialized
    d_is_MobilityMatrixInitialized=true;
    return;
}

void 
DenseMobilityMatrix::generateMobilityMatrix()
{ 
    if (!d_is_MobilityMatrixInitialized) TBOX_ERROR("DenseMobilityMatrix::generateMobilityMatrix(): The DenseMobilityMatrix must be initialized before matrix generation call." << std::endl);

    const int size = d_structure_no_nodes*NDIM;
    if (!d_friction_matrix) d_friction_matrix = new double[size*size]; 

    getLocalPositionWeightArrays(X,W);

    if (d_MATRIX_TYPE == FILE)
    {
	getTrueMobilityMatrix(Dt);
    }
    else if (d_MATRIX_TYPE == RPY)
    {
	getRPYMobilityMatrix(d_IBKernelName, d_mu, d_DX, X, d_structure_no_blobs, d_F_PERIODIC_CORRECTION, d_friction_matrix);
    }
    else if (d_MATRIX_TYPE == EMPIRICAL)
    {
	getEmpiricalMobilityMatrix(d_IBKernelName, d_mu, d_rho, Dt, d_DX, X, d_structure_no_blobs, 0, d_F_PERIODIC_CORRECTION, d_L_DOMAIN, d_friction_matrix);
    }
    else if (d_MATRIX_TYPE == HODLR)
    {
	TBOX_ERROR("DenseMobilityMatrix::generateMobilityMatrix(): HODLR is not implemented yet." <<std::endl);
    }
    else 
	TBOX_ERROR("DenseMobilityMatrix::generateMobilityMatrix(): Invalid type of a mobility matrix." <<std::endl);
    
    //Scale MM
    scaleMobilityMatrix();

    if (d_X) delete [] d_X;
    if (d_W) delete [] d_W;
    return;
}// generateMobilityMatrix

void 
DenseMobilityMatrix::generateFrictionMatrix()
{
    int info = 0;
    const int size = d_structure_no_blobs*NDIM;
    
    if (d_INVERSE_METHOD == LAPACK_CHOLESKY)
    {
	dpotrf_((char*) "L", size, d_friction_matrix, size, info);
    }
    else if (d_INVERSE_METHOD == LAPACK_LU) 
    {
	if (!d_IPIV) d_IPIV = new int[size];
	dgetrf_(size, size, d_friction_matrix, size, d_IPIV, info);
    }
    else if (d_INVERSE_METHOD == LAPACK_SVD) 
    {
	
	/* Locals */
	int il, iu, m, lwork, liwork;
	double abstol, vl, vu;
	int iwkopt;
	double wkopt;
	/* Local arrays */
	int *isuppz = new int[2*size];
	double *w = new double [size];
	double *z = new double [size*size];
	
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Query and allocate the optimal workspace */
	lwork = -1;
	liwork = -1;
	dsyevr_((char*)"V", (char*)"A",(char*)"L", size, d_friction_matrix, size, 
		vl, vu, il, iu, abstol, 
		m, w, z, size, isuppz, &wkopt, lwork, &iwkopt, liwork, info );
	
	lwork = (int) wkopt;
	double* work = new double[lwork];
	liwork = iwkopt;
	int* iwork = new int[liwork];
	/* Solve eigenproblem */
	dsyevr_((char*)"V", (char*)"A",(char*)"L", size, d_friction_matrix, size,
		vl, vu, il, iu,abstol, 
		m, w, z, size, isuppz, work, lwork, iwork, liwork, info );
	int counter=0, counter_zero=0;
	//make negative eigenvalues be equal min_value from input
	for (int ipart=0; ipart<size; ++ipart) 
	{
	    if (w[ipart]< svd_inverse_epsilon)
	    {
		w[ipart] = svd_inverse_tolerance;
		counter++;
	    }
	    
	}
	for (int ipart=0; ipart<size; ++ipart) 
	{
	    if (w[ipart]==0.) counter_zero++;
	    
	    for (int jpart=0; jpart<size; ++jpart) 
	    {
		if (w[jpart]==0.)
		{
		    d_friction_matrix[jpart*size+ipart] = 0.;
		}
		else
		{
		    if (ipart==jpart) 
			d_friction_matrix[jpart*size+ipart] = z[jpart*size+ipart]/sqrt(w[jpart]);
		    else
			d_friction_matrix[jpart*size+ipart] = z[jpart*size+ipart]/sqrt(w[jpart]);
		}
	    }
	}
	
	std::cout<<"DenseMobilityMatrix::generateFrictionMatrix(): \n For structure "<<d_struct_ID<<",  "<<counter
	    <<" eigenvalues have been changed, number of pseudoinverse "<<counter_zero<< std::endl;

	delete[] isuppz;
	delete[] w;
	delete[] z;
	delete[] work;
	delete[] iwork;
    }        
    else if (d_MATRIX_TYPE == HODLR) 
    {
	//hodlr_mobility_matrix->compute_Factor();
    }
    else 
    {
	TBOX_ERROR("DenseMobilityMatrix::generateFrictionMatrix(): Invalid method for inverting mobility matrix." <<std::endl);
    }
    
    if (info)
    {
	TBOX_ERROR("DenseMobilityMatrix::generateFrictionMatrix(): Inverting mobility matrix failed for structure "<<d_struct_ID<<",  with err " << info<< std::endl);
    }
    return;
}// generateFrictionMatrix


void
DenseMobilityMatrix::computeSolution(
    PetscScalar *rhs_b)
{
    int info=0;
    const int size = d_structure_no_blobs*NDIM;
    
    d_b=new double[size];
    for (int i = 0; i < size; i++) 
    {
	d_b[i] = rhs_b[i];
    }
    
    if (d_MATRIX_TYPE == HODLR) 
    {

    }else if (d_INVERSE_METHOD == LAPACK_CHOLESKY)
    {
	dpotrs_((char*)"L", size, 1, d_friction_matrix, size, d_b, size, info);
    }
    else if (d_INVERSE_METHOD == LAPACK_LU)
    {
	dgetrs_((char*)"N", size, 1, d_friction_matrix, size, d_IPIV, d_b, size, info);
    }
    else if (d_INVERSE_METHOD == LAPACK_SVD) 
    {
	double *d_b2 = new double[size];
	for (int ipart=0; ipart<size; ++ipart) 
	{	
	    d_b2[ipart]=0.;
	    for (int jpart=0; jpart<size; ++jpart) 
	    {
		d_b2[ipart] +=d_friction_matrix[ipart*size+jpart]*d_b[jpart];
	    }
	}
	
	for (int ipart=0; ipart<size; ++ipart) 
	{	
	    d_b[ipart]=0.;
	    for (int jpart=0; jpart<size; ++jpart) 
	    {
		d_b[ipart] +=d_friction_matrix[jpart*size+ipart]*d_b2[jpart];
	    }
	}
	delete[] d_b2;
    }
    if (info)
    {
	TBOX_ERROR("DenseMobilityMatrix::computeSolution() Solution failed, err:" << info << std::endl);
    }
    
    for (int i = 0; i < size; i++) 
    {
	rhs_b[i] = d_b[i];
    }

    delete [] d_b;
    return;
}// computeSolution 

int 
DenseMobilityMatrix::getStructID()
{
    return d_struct_ID;
}

/////////////////////////////// PRIVATE //////////////////////////////////////
void 
DenseMobilityMatrix::scaleMobilityMatrix()
{
    const int size = d_structure_no_blobs*NDIM;
    for (int ipart = 0; ipart < size; ++ipart) 
	for (int jpart = 0; jpart < size; ++jpart) 
	{
	    d_friction_matrix[ipart*size+jpart] *= d_mob_scale_1;
	    
	    if (ipart==jpart)
	    {
		d_friction_matrix[ipart*size+jpart] += d_mob_scale_2*d_W[ipart];
		
	    }
	}	
}
void
DenseMobilityMatrix::readMobilityMatrixFromFile(
    const double DDt)
{

//later

}//end readTrueMobilityMatrix

}// IBAMR


