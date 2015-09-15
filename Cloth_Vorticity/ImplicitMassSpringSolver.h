#pragma once

#include "dynamics_solver_interface.h"
#include "cloth_data.h"
#include <Eigen\Dense>
#include "sparseMatrix.h"
#include "insertRows.h"
#include "CGSolver.h"
#include "PardisoSolver.h"

class ImplicitMassSpringSolver :
	public Dynamics_Solver_Interface
	{
	public:
		ImplicitMassSpringSolver(void);
		~ImplicitMassSpringSolver(void);

		void initialize(Cloth_Data *cloth);
		void generateAllSprings(Cloth_Data *cloth);
		void initializeSparseMatrixFromOutline(Cloth_Data *cloth);
		
		void resetParameters();

		virtual void advance_time_step(Cloth_Data* cloth);

		void addShearSprings(Cloth_Data* cloth);
		void addBendingSprings(Cloth_Data* cloth);
		void addAllSprings(Cloth_Data* cloth); //Use this NOT THE OTHER TWO BECAUSE YOUR FRIGGIN LIFE IS ALREADY COMPLEX ENOUGH WITHOUT TWO FUNCTIONS
		void addGravityComponents(Cloth_Data* cloth);
		//void addDampingComponents(Cloth_Data* cloth);
		//void addBendingComponents(Cloth_Data* cloth);
		void finalizeAndSolve(Cloth_Data* cloth);

	private:
		Eigen::VectorXd RHSVector_;
		SparseMatrix *LHSMatrix_;
		int lastFrameId_;
		double *massMatrixDiagonal_;

		//Constrain information 
		//@TODO Later on shift to cloth_data 
		int numConstrainedVerts_;
		int *constrainedVerts_;

		//All the spring information for faster processing 
		std::vector<std::pair<int,int> > stretchSprings_;
		std::vector<double> restLenStretchSprings_;
		std::vector<std::pair<int,int> > bendingSprings_;
		std::vector<double> restLenbendingSprings_;

		//Combine them into one comprehensuive stuff so that I dont have to differential processing 
		std::vector<std::pair<int,int> > allSprings_;
		std::vector<double> restLenAllSprings_;

	};

