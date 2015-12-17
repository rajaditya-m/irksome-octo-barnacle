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
		void computedFdU(Cloth_Data* cloth);
		~ImplicitHyperElasticFEMSolver(void);

		void initialize(Cloth_Data *cloth);
		void initializeSparseMatrixFromOutline(Cloth_Data *cloth);
		
		void resetParameters();

		virtual void advance_time_step(Cloth_Data* cloth);
		
		double computeGammaValues(int i,int j,std::vector<double> &sigma,double IIIC, std::vector<double> &gradient,std::vector<double> &hessian);
		int compute4x4TensorIndex(int i, int j, int m, int n);
		int compute4x9TensorIndex(int i, int j, int m, int n);
	  int compute6x9TensorIndex(int i, int j, int m, int n);
	  void computeDPDF_Hat(std::vector<double>& sigma, std::vector<double>& gradients, std::vector<double>& hessians, double IIIC, std::vector<double>& DPDF_Hat);
		Eigen::MatrixXd convert6VectorToEigen3x2Matrix(std::vector<double> &vec);
		Eigen::MatrixXd convert4VectorToEigen2x2Matrix(std::vector<double> &vec);
		void computeDPDF(std::vector<double> &DPDF_Hat, Eigen::MatrixXd &U, Eigen::MatrixXd &V,std::vector<double> &DPDF);
		void computeDGDF(std::vector<double> &DPDF, Eigen::Matrix2d bVec, std::vector<double> &DGDF);
		void computeElementStiffnessMatrix(int el, std::vector<double>& DGDF, std::vector<double>& KELEM);
		void addElementStiffnessMatrixToGlobalStiffnessMatrix(Cloth_Data* cloth, int el, std::vector<double>& KELEM);
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
		SparseMatrix* tangentStiffnessMatrix_;

		int lastFrameId_;
		double *massMatrixDiagonal_;

		//Constrain information 
		//@TODO Later on shift to cloth_data 
		int numConstrainedVerts_;
		int *constrainedVerts_;

		//That weird ordering scheme for some reason 
		std::vector<int> rowMajorMatrixToTeran;
		std::vector<int> teranToRowMajorMatrix;
		
		std::vector<double> dDSdU;
		std::vector<double> dFdUs;
	};

