#pragma once

#define OPEN_MP
#ifdef OPEN_MP
#  include <omp.h>
#  ifndef _MSC_VER
#    define OMP_FOR _Pragma("omp parallel for")
#  else
#    define OMP_FOR __pragma(omp parallel for)
#  endif // vs compiler
#else // no omp
#  define OMP_FOR
#endif // #ifndef OPEN_MP

#include "dynamics_solver_interface.h"

#include <Eigen/Dense>

#include "sparseMatrix.h"
#include "cloth_data.h"
#include "geom_funcs.h"

#define ADDCOMPRESSIONRESISTENCE 1

class ImplicitHyperElasticFEMSolver :
	public Dynamics_Solver_Interface
	{
	public:
		ImplicitHyperElasticFEMSolver(void);
		~ImplicitHyperElasticFEMSolver(void);

		void initialize(Cloth_Data *cloth);
		//void initializeSparseMatrixFromOutline(Cloth_Data *cloth);
		
		void resetParameters();

		virtual void advance_time_step(Cloth_Data* cloth);

		void addShearComponents(Cloth_Data* cloth);
		void addGravityComponents(Cloth_Data* cloth);
		void addLinearDamping(Cloth_Data* cloth);
		//void addBendingComponents(Cloth_Data* cloth);
		void finalizeAndSolve(Cloth_Data* cloth);


	private:
		Eigen::VectorXd force_;

		//Later on we will need this 
		Eigen::VectorXd RHSVector_;
		SparseMatrix *LHSMatrix_;

		int lastFrameId_;
		double *massMatrixDiagonal_;

		//Constrain information 
		//@TODO Later on shift to cloth_data 
		int numConstrainedVerts_;
		int *constrainedVerts_;
	};

