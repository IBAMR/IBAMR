srcs-ibtk.cpp += $(call thisdir, \
	BGaussSeidelPreconditioner.cpp \
	BJacobiPreconditioner.cpp \
	CCLaplaceOperator.cpp \
	CCPoissonBoxRelaxationFACOperator.cpp \
	CCPoissonLevelRelaxationFACOperator.cpp \
	CCPoissonPETScLevelSolver.cpp \
	CCPoissonPointRelaxationFACOperator.cpp \
	CCPoissonSolverManager.cpp \
	FACPreconditioner.cpp \
	KrylovLinearSolverManager.cpp \
	KrylovLinearSolverPoissonSolverInterface.cpp \
	LaplaceOperator.cpp \
	NewtonKrylovSolverManager.cpp \
	PETScKrylovLinearSolver.cpp \
	PETScKrylovPoissonSolver.cpp \
	PETScLevelSolver.cpp \
	PETScMFFDJacobianOperator.cpp \
	PETScNewtonKrylovSolver.cpp \
	PoissonFACPreconditioner.cpp \
	PoissonFACPreconditionerStrategy.cpp \
	PoissonSolver.cpp \
	SCLaplaceOperator.cpp \
	SCPoissonPETScLevelSolver.cpp \
	SCPoissonPointRelaxationFACOperator.cpp \
	SCPoissonSolverManager.cpp \
	VCSCViscousOpPointRelaxationFACOperator.cpp \
	VCSCViscousOperator.cpp \
	VCSCViscousPETScLevelSolver.cpp \
	)

ifneq ($(HYPRE_LIB),)
srcs-ibtk.cpp += $(call thisdir, \
	CCPoissonHypreLevelSolver.cpp \
	SCPoissonHypreLevelSolver.cpp \
	)
endif

include $(call incsubdirs,fortran)
