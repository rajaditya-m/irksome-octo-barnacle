#include "ImplicitHyperElasticFEMSolver.h"

#ifndef ELT
  #define ELT(numRows,i,j) (((long)j)*((long)numRows)+((long)i))
#endif

ImplicitHyperElasticFEMSolver::ImplicitHyperElasticFEMSolver(void)
{
  rowMajorMatrixToTeran.resize(4);
  rowMajorMatrixToTeran[0] = 0;
	rowMajorMatrixToTeran[1] = 2;
	rowMajorMatrixToTeran[2] = 3;
	rowMajorMatrixToTeran[3] = 1;

	teranToRowMajorMatrix.resize(4);
	for(int i=0;i<4;i++)
	{
		teranToRowMajorMatrix[rowMajorMatrixToTeran[i]] = i;
	}

	dDSdU.resize(54,0);
		
	dDSdU[compute6x9TensorIndex(0, 0, 0, 0)] = -1.0;
	dDSdU[compute6x9TensorIndex(1, 0, 0, 1)] = -1.0;
	dDSdU[compute6x9TensorIndex(2, 0, 0, 2)] = -1.0;
	dDSdU[compute6x9TensorIndex(0, 1, 1, 0)] = -1.0;
	dDSdU[compute6x9TensorIndex(1, 1, 1, 1)] = -1.0;
	dDSdU[compute6x9TensorIndex(2, 1, 1, 2)] = -1.0;

	dDSdU[compute6x9TensorIndex(0, 0, 2, 0)] = 1.0;
	dDSdU[compute6x9TensorIndex(0, 1, 2, 0)] = 1.0;
	dDSdU[compute6x9TensorIndex(1, 0, 2, 1)] = 1.0;
	dDSdU[compute6x9TensorIndex(1, 1, 2, 1)] = 1.0;
	dDSdU[compute6x9TensorIndex(2, 0, 2, 2)] = 1.0;
	dDSdU[compute6x9TensorIndex(2, 1, 2, 2)] = 1.0;

}

/*void ImplicitHyperElasticFEMSolver::computedFdU(Cloth_Data* cloth) {
	int numElements = cloth->getMesh()->get_number_triangles();
  for (int el = 0; el < numElements; el++) {
    double * dFdU = &dFdUs[36 * el];
		Eigen::Matrix2d dmInv = cloth->getDmInv(el);
    for (int index = 0; index < 36; index++) {
      int n = index % 3;
      int m = (int)(index / 3) % 3;
      int j = (int)(index / 9) % 2;
      int i = (int)(index / 18) % 2;
      double result = 0.0;
      for (int k = 0; k < 2; k++)
        result += dDSdU[compute4x9TensorIndex(i, k, m, n)] * dmInv(k,j);
      dFdU[compute4x9TensorIndex(i, j, m, n)] = result;
    }
  }
}*/

void ImplicitHyperElasticFEMSolver::initializeSparseMatrixFromOutline(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
  int numElementVertices = 3;
	std::vector<int> vertices(numElementVertices);

  // build the non-zero locations of the tangent stiffness matrix
  SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);
	int numElements = cloth->getMesh()->get_number_triangles();

  for (int el = 0; el < numElements; el++) {

		Triangles tri = cloth->getMesh()->get_triangle(el);
    vertices[0] = tri.a;
    vertices[1] = tri.b;
    vertices[2] = tri.c;

    for (int i = 0; i < numElementVertices; i++)
      for (int j = 0; j < numElementVertices; j++) {
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++) {
            // only add the entry if both vertices are free (non-fixed)
            // the corresponding elt is in row 3*i+k, column 3*j+l
            emptyMatrix->AddEntry( 3 * vertices[i] + k, 3 * vertices[j] + l, 0.0 );
          }
      }
  }

  tangentStiffnessMatrix_ = new SparseMatrix(emptyMatrix);
  delete(emptyMatrix);

	rayleighDampingMatrix_ = new SparseMatrix(*tangentStiffnessMatrix_);
}

void ImplicitHyperElasticFEMSolver::initializeMassMatrixFromOutline(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);
	for(int p=0;p<numVertices;p++)
	{
		float vertexMass = cloth->get_vertex_mass(p);
		for(int d=0;d<3;d++)
		{
			emptyMatrix->AddEntry(3*p+d,3*p+d,vertexMass);
		}
	}
	massMatrix_ = new SparseMatrix(emptyMatrix);
	delete emptyMatrix;
}

void ImplicitHyperElasticFEMSolver::computedFdU(Cloth_Data* cloth) {
	int numElements = cloth->getMesh()->get_number_triangles();
  for (int el = 0; el < numElements; el++) {
    double * dFdU = &dFdUs[54 * el];
		Eigen::Matrix2d dmInv = cloth->getDmInv(el);
    for (int index = 0; index < 54; index++) {
      int n = index % 3;
      int m = (int)(index / 3) % 3;
      int j = (int)(index / 9) % 2;
      int i = (int)(index / 18) % 3;
      double result = 0.0;
      for (int k = 0; k < 2; k++)
        result += dDSdU[compute6x9TensorIndex(i, k, m, n)] * dmInv(k,j);
      dFdU[compute6x9TensorIndex(i, j, m, n)] = result;
    }
  }
}

ImplicitHyperElasticFEMSolver::~ImplicitHyperElasticFEMSolver(void)
{
  delete tangentStiffnessMatrix_;
	delete massMatrix_;
	delete rayleighDampingMatrix_;
}

void ImplicitHyperElasticFEMSolver::initialize(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();

	massMatrixDiagonal_ = new double[3*numVertices];
	//OMP_FOR
	for(int p=0; p<numVertices;p++) {
		float vertex_mass = cloth->get_vertex_mass(p);
		massMatrixDiagonal_[3*p+0] = vertex_mass;
		massMatrixDiagonal_[3*p+1] = vertex_mass;
		massMatrixDiagonal_[3*p+2] = vertex_mass;
	}

	force_ = Eigen::VectorXd::Zero(numVertices*3);

	//Constrained verts
	numConstrainedVerts_ = 2;
	constrainedVerts_ = new int[6];
	constrainedVerts_[3] = 3*CORNERV2+0;
	constrainedVerts_[4] = 3*CORNERV2+1;
	constrainedVerts_[5] = 3*CORNERV2+2;
	constrainedVerts_[0] = 3*CORNERV1+0;
	constrainedVerts_[1] = 3*CORNERV1+1;
	constrainedVerts_[2] = 3*CORNERV1+2;

	//Add the complemented stuff 
	dFdUs.resize(54*numTriangles,0);
	computedFdU(cloth);

	//Compute the tangent Stiffness Matrix Structure 
	initializeSparseMatrixFromOutline(cloth);

	//Compute the mass matrix structure 
	initializeMassMatrixFromOutline(cloth);

	//Establish substructures
	rayleighDampingMatrix_->BuildSubMatrixIndices(*massMatrix_);
  tangentStiffnessMatrix_->BuildSubMatrixIndices(*massMatrix_);

	//Build the correct internal and external forces vector 
	internalForce_ =  Eigen::VectorXd::Zero(numVertices*3);
	externalForce_ =  Eigen::VectorXd::Zero(numVertices*3);
	residuals_ = Eigen::VectorXd::Zero(numVertices*3);
	velocities_ = Eigen::VectorXd::Zero(numVertices*3);

}

void ImplicitHyperElasticFEMSolver::advance_time_step(Cloth_Data* cloth)
{
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;
	
	//Preapre the components as usual 
	tangentStiffnessMatrix_->ResetToZero();
	rayleighDampingMatrix_->ResetToZero();

	//Prepare the RHS Vector for usage
	force_.setZero();
	internalForce_.setZero();
	externalForce_.setZero();
	residuals_.setZero();
	velocities_.setZero();

	//Add the physics components
	addShearComponents(cloth);
	addGravityComponents(cloth);

	//Solve and report
	//finalizeAndSolve(cloth);
	finalizeAndSolve_Explicit(cloth);
}

double ImplicitHyperElasticFEMSolver::computeGammaValues(int i,int j,std::vector<double> &sigma,double IIIC, std::vector<double> &gradient,std::vector<double> &hessian)
{
	double tempGammaVec1[3];
  tempGammaVec1[0] = 2.0 * sigma[i];
  tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
  //tempGammaVec1[2] = 2.0 * IIIC / sigma[i];
  tempGammaVec1[2] = 2.0 * ((sigma[i]*sigma[i]*sigma[j]*sigma[j])/sigma[i]);
  double tempGammaVec2[3];
  tempGammaVec2[0] = 2.0 * sigma[j];
  tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
  //tempGammaVec2[2] = 2.0 * IIIC / sigma[j];
  tempGammaVec2[2] = 2.0 * ((sigma[i]*sigma[i]*sigma[j]*sigma[j])/sigma[j]);
  double productResult[3];
  productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] +
                      tempGammaVec2[2] * hessian[2]);
  productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] +
                      tempGammaVec2[2] * hessian[4]);
  productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] +
                      tempGammaVec2[2] * hessian[5]);
  //return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
  //       tempGammaVec1[2] * productResult[2] + 4.0 * IIIC * gradient[2] / (sigma[i] * sigma[j]));
	return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
          tempGammaVec1[2] * productResult[2] + 4.0 * gradient[2] * sigma[i] * sigma[j]);
}

int ImplicitHyperElasticFEMSolver::compute4x4TensorIndex(int i, int j, int m, int n)
{
	int rowIndex_in9x9Matrix = rowMajorMatrixToTeran[2 * i + j];
  int columnIndex_in9x9Matrix = rowMajorMatrixToTeran[2 * m + n];
  return (4 * rowIndex_in9x9Matrix + columnIndex_in9x9Matrix);
}

int ImplicitHyperElasticFEMSolver::compute6x9TensorIndex(int i, int j, int m, int n) {
  int rowIndex = 2 * i + j;
  int columnIndex = 3 * m + n;
	if(9*rowIndex+columnIndex>=54)
	{
		std::cout << i << " " << j << " " << m << " " << n << " " << (9*rowIndex+columnIndex) << "\n";
		return 0;
	}
  return (9 * rowIndex + columnIndex);
}

int ImplicitHyperElasticFEMSolver::compute6x6TensorIndex(int i, int j, int m, int n)
{
	int rowIndex = i * 2 + j;
	int colIndex = m * 2 + n;
	if(6*rowIndex+colIndex>=36)
	{
		std::cout << i << " " << j << " " << m << " " << n << " " << (6*rowIndex+colIndex) << "\n";
		return 0;
	}
	return (6*rowIndex+colIndex);
}

/*void ImplicitHyperElasticFEMSolver::computeDPDF_Hat(std::vector<double> &sigma, std::vector<double> &gradients, std::vector<double> &hessians, double IIIC,std::vector<double> &DPDH_Hat)
{
  double l1_Sq = sigma[0]*sigma[0];
  double l2_Sq = sigma[1]*sigma[1];

	double alpha11 = 2.0 * gradients[0] + 8.0 * l1_Sq * gradients[1];
	double alpha22 = 2.0 * gradients[0] + 8.0 * l2_Sq * gradients[1];

	double alpha12 = 2.0 * gradients[0] + 4.0 * (l1_Sq + l2_Sq) * gradients[1];

	double beta11 = 4.0 * l1_Sq * gradients[1] - (2.0 * IIIC * gradients[2]) / l1_Sq;
	double beta22 = 4.0 * l2_Sq * gradients[1] - (2.0 * IIIC * gradients[2]) / l2_Sq;

	double beta12 = 4.0 * sigma[0] * sigma[1] * gradients[1] - (2.0 * IIIC *gradients[2]) / (sigma[0] * sigma[1]);

	double gamma11 = computeGammaValues(0, 0, sigma, IIIC, gradients, hessians);
	double gamma22 = computeGammaValues(1, 1, sigma, IIIC, gradients, hessians);
	double gamma12 = computeGammaValues(0, 1, sigma, IIIC, gradients, hessians);

	double x1111, x2222;
	double x2211;
	double x2121;
	double x2112;

	x1111 = alpha11 + beta11 + gamma11;
	x2222 = alpha22 + beta22 + gamma22;

	x2211 = gamma12;
	x2121 = beta12;
	x2112 = alpha12;

	//x2211 = beta12;
	//x2121 = gamma12;
	//x2112 = alpha12;
	
	DPDH_Hat.resize(16,0);

#ifdef ENSURE_SEMI_POSITIVEDEFINITENESS
	Eigen::Matrix2d A;
	A(0,0) = x1111; A(0,1) = x2211;
	A(1,0) = x2211; A(1,1) = x2222;
	Eigen::EigenSolver<Eigen::Matrix2d> eigsolver1(A);
	Eigen::Vector2d d = eigsolver1.eigenvalues().real();
	Eigen::Matrix2d v = eigsolver1.eigenvectors().real();
	bool corrected = false;
	for (int i = 0; i < 2 && d[i] < 0; ++i) {
		d[i] = 0;
		corrected = true;
	}
	Eigen::Matrix2d tmp;
	if (corrected) {
		tmp(0,0) = d(0)*v(0,0); tmp(0,1)  = d(1)*v(0,1);
		tmp(1,0) = d(0)*v(1,0); tmp(1,1)  = d(1)*v(1,1);

		Eigen::Matrix2d A = v.transpose()*tmp;

		double x1111 = A(0,0);
		double x2222 = A(1,1);
		double x2211 = A(0,1);

		DPDH_Hat[compute4x4TensorIndex(0, 0, 0, 0)] = x1111;
		DPDH_Hat[compute4x4TensorIndex(0, 0, 1, 1)] = x2211;

		DPDH_Hat[compute4x4TensorIndex(1, 1, 0, 0)] = x2211;
		DPDH_Hat[compute4x4TensorIndex(1, 1, 1, 1)] = x2222;
	} 
	else 
	{
		DPDH_Hat[compute4x4TensorIndex(0, 0, 0, 0)] = x1111;
		DPDH_Hat[compute4x4TensorIndex(0, 0, 1, 1)] = x2211;

		DPDH_Hat[compute4x4TensorIndex(1, 1, 0, 0)] = x2211;
		DPDH_Hat[compute4x4TensorIndex(1, 1, 1, 1)] = x2222;
	}

	A(0,0) = x2121; A(0,1) = x2112;
	A(1,0) = x2112; A(1,1) = x2121;
	Eigen::EigenSolver<Eigen::Matrix2d> eigsolver2(A);
	d = eigsolver2.eigenvalues().real();
	v = eigsolver2.eigenvectors().real();
	corrected = false;
	for (int i = 0; i < 2 && d[i] < 0; ++i) {
		d[i] = 0;
		corrected = true;
	}
	if (corrected) 
	{
		tmp(0,0) = d(0)*v(0,0); tmp(0,1)  = d(1)*v(0,1);
		tmp(1,0) = d(0)*v(1,0); tmp(1,1)  = d(1)*v(1,1);

		Eigen::Matrix2d A = v.transpose()*tmp;

		double x2121 = A(0,0);
		double x2112 = A(0,1);

		DPDH_Hat[compute4x4TensorIndex(0, 1, 0, 1)] = x2121;
		DPDH_Hat[compute4x4TensorIndex(0, 1, 1, 0)] = x2112;

		DPDH_Hat[compute4x4TensorIndex(1, 0, 0, 1)] = x2112;
		DPDH_Hat[compute4x4TensorIndex(1, 0, 1, 0)] = x2121;
	} 
	else
	{
		DPDH_Hat[compute4x4TensorIndex(0, 1, 0, 1)] = x2121;
		DPDH_Hat[compute4x4TensorIndex(0, 1, 1, 0)] = x2112;

		DPDH_Hat[compute4x4TensorIndex(1, 0, 0, 1)] = x2112;
		DPDH_Hat[compute4x4TensorIndex(1, 0, 1, 0)] = x2121;
	}
#else
	DPDH_Hat[compute4x4TensorIndex(0, 0, 0, 0)] = x1111;
	DPDH_Hat[compute4x4TensorIndex(1, 1, 0, 0)] = x2211;

	DPDH_Hat[compute4x4TensorIndex(0, 1, 0, 1)] = x2121;
	DPDH_Hat[compute4x4TensorIndex(1, 0, 0, 1)] = x2112;

	DPDH_Hat[compute4x4TensorIndex(0, 1, 1, 0)] = x2112;
	DPDH_Hat[compute4x4TensorIndex(1, 0, 1, 0)] = x2121;

	DPDH_Hat[compute4x4TensorIndex(0, 0, 1, 1)] = x2211;
	DPDH_Hat[compute4x4TensorIndex(1, 1, 1, 1)] = x2222;

#endif
}*/

void ImplicitHyperElasticFEMSolver::computeDPDF_Hat(std::vector<double> &sigma, std::vector<double> &gradients, std::vector<double> &hessians, double IIIC,std::vector<double> &DPDH_Hat)
{
  double l1_Sq = sigma[0]*sigma[0];
  double l2_Sq = sigma[1]*sigma[1];

	double l1_Cube = sigma[0]*sigma[0]*sigma[0];
	double l2_Cube = sigma[1]*sigma[1]*sigma[1];

	DPDH_Hat.resize(36,0);

	double x0000 = 12.0*l1_Sq*gradients[1] + 2.0*l2_Sq*gradients[2] + 2.0*gradients[0] + 4.0*l1_Cube*(2.0*sigma[0]*l2_Sq*hessians[4] + 2.0*sigma[0]*hessians[1] + 4.0*l1_Cube*hessians[3]) 
                                                                              + 2.0*sigma[0]*l2_Sq*(2.0*sigma[0]*l2_Sq*hessians[5] + 2.0*sigma[0]*hessians[2] + 4.0*l1_Cube*hessians[4])
																																							      + 2.0*sigma[0]*(2.0*sigma[0]*l2_Sq*hessians[2] + 2.0*sigma[0]*hessians[0] + 4.0*l1_Cube*hessians[1]) ;


	double x1100 = 4.0*sigma[0]*sigma[1]*gradients[2] + 4.0*l1_Cube*(2.0*sigma[1]*l1_Sq*hessians[4] + 2.0*sigma[1]*hessians[1] + 4.0*l2_Cube*hessians[3]) 
                                             + 2.0*sigma[0]*l2_Sq*(2.0*sigma[1]*l1_Sq*hessians[5] + 2.0*sigma[1]*hessians[2] + 4.0*l2_Cube*hessians[4])
																							     + 2.0*sigma[0]*(2.0*sigma[1]*l1_Sq*hessians[2] + 2.0*sigma[1]*hessians[0] + 4.0*l2_Cube*hessians[1]) ;

	double x0101 = (4.0*l1_Sq + 4.0*l2_Sq)*gradients[1] + 2.0*gradients[0];
	double x1001 = 2.0*sigma[1]*sigma[2]*(2.0*gradients[1] - gradients[2]);

	double x0011 = 4.0*sigma[0]*sigma[1]*gradients[2] + 4.0*l2_Cube*(2.0*sigma[0]*l2_Sq*hessians[4] + 2.0*sigma[0]*hessians[1] + 4.0*l1_Cube*hessians[3]) 
                                             + 2.0*sigma[1]*l1_Sq*(2.0*sigma[0]*l2_Sq*hessians[5] + 2.0*sigma[0]*hessians[2] + 4.0*l1_Cube*hessians[4])
																							     + 2.0*sigma[1]*(2.0*sigma[0]*l2_Sq*hessians[2] + 2.0*sigma[0]*hessians[0] + 4.0*l1_Cube*hessians[1]) ;


	double x1111 = 12.0*l2_Sq*gradients[1] + 2.0*l1_Sq*gradients[2] + 2.0*gradients[0] + 4.0*l2_Cube*(2.0*sigma[1]*l1_Sq*hessians[4] + 2.0*sigma[1]*hessians[1] + 4.0*l2_Cube*hessians[3]) 
                                                                              + 2.0*sigma[1]*l1_Sq*(2.0*sigma[1]*l1_Sq*hessians[5] + 2.0*sigma[1]*hessians[2] + 4.0*l2_Cube*hessians[4])
																																							      + 2.0*sigma[1]*(2.0*sigma[1]*l1_Sq*hessians[2] + 2.0*sigma[1]*hessians[0] + 4.0*l2_Cube*hessians[1]) ;

	double x2020 = 4.0*l1_Sq*gradients[1] + 2.0*l2_Sq*gradients[2] + 2.0*gradients[0];
	double x2121 = 4.0*l2_Sq*gradients[1] + 2.0*l1_Sq*gradients[2] + 2.0*gradients[0];

	DPDH_Hat[compute6x6TensorIndex(0,0,0,0)] = x0000;
	DPDH_Hat[compute6x6TensorIndex(1,1,0,0)] = x1100;

	DPDH_Hat[compute6x6TensorIndex(0,1,0,1)] = x0101;
	DPDH_Hat[compute6x6TensorIndex(1,0,0,1)] = x1001;
	
	DPDH_Hat[compute6x6TensorIndex(0,1,1,0)] = x1001;
	DPDH_Hat[compute6x6TensorIndex(1,0,1,0)] = x0101;

	DPDH_Hat[compute6x6TensorIndex(0,0,1,1)] = x0011;
	DPDH_Hat[compute6x6TensorIndex(1,1,1,1)] = x1111;

	DPDH_Hat[compute6x6TensorIndex(2,0,2,0)] = x2020;

	DPDH_Hat[compute6x6TensorIndex(2,1,2,1)] = x2121;

	
}

Eigen::MatrixXd ImplicitHyperElasticFEMSolver::convert6VectorToEigen3x2Matrix(std::vector<double> &vec)
{
	Eigen::MatrixXd res(3,2);
	res(0,0) = vec[0]; res(0,1) = vec[1];
	res(1,0) = vec[2]; res(1,1) = vec[3];
	res(2,0) = vec[4]; res(2,1) = vec[5];
	return res;
}

Eigen::MatrixXd ImplicitHyperElasticFEMSolver::convert4VectorToEigen2x2Matrix(std::vector<double> &vec)
{
	Eigen::MatrixXd res(2,2);
	res(0,0) = vec[0]; res(0,1) = vec[1];
	res(1,0) = vec[2]; res(1,1) = vec[3];
	return res;
}

/*
void ImplicitHyperElasticFEMSolver::computeDPDFFrustrating(Eigen::MatrixXd &F,std::vector<double> &gradients, std::vector<double> &hessians,std::vector<double> &DPDF)
{
	//First some convinienet shorthands because whynot 
  double x1 = F(0,0); double x2 = F(0,1);
  double y1 = F(1,0); double y2 = F(1,1);
  double z1 = F(2,0); double z2 = F(2,1);

	double x1_2 = x1*x1; double x2_2 = x2*x2;
	double y1_2 = y1*y1; double y2_2 = y2*y2;
	double z1_2 = z1*z1; double z2_2 = z2*z2;

	DPDF.resize(36,0);

	//A total of 36 numbers 
	DPDF[compute6x6TensorIndex(0,0,0,0)] = (2.0*y2_2 + 2.0*z2_2)*gradients[2] + (8.0*x1_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4*x1*(x1_2 + y1_2 + z1_2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*x1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,0)] = (-2.0*y1*y2 -2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*x1)* ((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1] + 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*x1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,0)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*x1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,0)] = -2.0*x2*z2*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*x1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,0)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] +2.0*z2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] +2.0*z2*hessians[1])
	+ (2.0*x1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] +2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2*y1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,0)] = (4.0*x2*y1 -2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*y1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,0)] = (2.0*x2_2 + 2.0*z2_2)*gradients[2] + (8.0*y1_2 + 4.0*y2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*y1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,0)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*y1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +2.0*z1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +2.0*z1*hessians[1])
	+ (2.0*y1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,0)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*y1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);
	
	DPDF[compute6x6TensorIndex(0,0,2,0)] = (-2.0*x2*z2)*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*z1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,0)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*z1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*z1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,0)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*z1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,0)] = (2.0*x2_2 + 2.0*y2_2)*gradients[2] + (8.0*z1_2 + 4.0*(x1_2 + y1_2 + z1_2) + 4.0*z2_2)*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*z1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,0)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*z1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,0,1)] = (-2.0*y1*y2 - 2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*x2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,1)] = (2.0*y1_2 + 2.0*z1_2)*gradients[2] + (4.0*x1_2 + 8.0*x2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*x2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,1)] = (4.0*x2*y1 - 2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*x2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*x2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,1)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +  2.0*z1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +  2.0*z1*hessians[1])
	+ (2.0*x2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +  2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*x2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,1)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*y2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1] 
	+  (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*y2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,1)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*y2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,1)] = (2.0*x1_2 + 2.0*z1_2)*gradients[2] + (4.0*y1_2 + 8.0*y2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.*gradients[0]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*y2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,1)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] + 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*y2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*y2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,2,1)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,1)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,1)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,1)] = (2.0*x1_2 + 2.0*y1_2)*gradients[2] + (4.0*z1_2 + 8.0*z2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
  + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	//DBG__Matrix_Dump(&DPDF[0],6,6);
}
*/

void ImplicitHyperElasticFEMSolver::computeDPDFFrustrating(Eigen::MatrixXd &F,std::vector<double> &gradients, std::vector<double> &hessians,std::vector<double> &DPDF)
{
	//First some convinienet shorthands because whynot 
  double x1 = F(0,0); double x2 = F(0,1);
  double y1 = F(1,0); double y2 = F(1,1);
  double z1 = F(2,0); double z2 = F(2,1);

	double x1_2 = x1*x1; double x2_2 = x2*x2;
	double y1_2 = y1*y1; double y2_2 = y2*y2;
	double z1_2 = z1*z1; double z2_2 = z2*z2;

	DPDF.resize(36,0);

	//A total of 36 numbers 
	DPDF[compute6x6TensorIndex(0,0,0,0)] = (2.0*y2_2 + 2.0*z2_2)*gradients[2] + (8.0*x1_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4*x1*(x1_2 + y1_2 + z1_2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*x1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,0)] = (-2.0*y1*y2 -2.0*z1*z2)*gradients[2]
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*x1)* ((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1] + 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*x1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,0)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*x1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,0)] = -2.0*x2*z2*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*x1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,0)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1] 
	+ (-2.0*x2*y1*y2 + 2.0*x1*y2_2 -2.0*x2*z1*z2 + 2.0*x1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] +2.0*z2*hessians[2])
	+ (4.0*x1*(x1_2 + y1_2 + z1_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] +2.0*z2*hessians[1])
	+ (2.0*x1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] +2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,0)] = -2.0*x2*y2*gradients[2] + (8.0*x1*y1 + 4.0*x2*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2*y1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,0)] = (4.0*x2*y1 -2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*y1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,0)] = (2.0*x2_2 + 2.0*z2_2)*gradients[2] + (8.0*y1_2 + 4.0*y2_2 + 4.0*(x1_2 + y1_2 + z1_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*y1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,0)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*y1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +2.0*z1*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +2.0*z1*hessians[1])
	+ (2.0*y1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,0)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*y1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);
	
	DPDF[compute6x6TensorIndex(0,0,2,0)] = (-2.0*x2*z2)*gradients[2] + (8.0*x1*z1 + 4.0*x2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*z1)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,0)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*z1)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,0)] = (-2.0*y2*z2)*gradients[2] + (8.0*y1*z1 + 4.0*y2*z2)*gradients[1]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*z1)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,0)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*z1)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,0)] = (2.0*x2_2 + 2.0*y2_2)*gradients[2] + (8.0*z1_2 + 4.0*(x1_2 + y1_2 + z1_2) + 4.0*z2_2)*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*z1)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,0)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*z1)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,0,1)] = (-2.0*y1*y2 - 2.0*z1*z2)*gradients[2] + (4.0*x1*x2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*x2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,0,1)] = (2.0*y1_2 + 2.0*z1_2)*gradients[2] + (4.0*x1_2 + 8.0*x2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*x2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,0,1)] = (4.0*x2*y1 - 2.0*x1*y2)*gradients[2] + (4.0*x1*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*x2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,0,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*x2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,0,1)] = (4.0*x2*z1 - 2.0*x1*z2)*gradients[2] + (4.0*x1*z2)*gradients[1] 
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] +  2.0*z1*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] +  2.0*z1*hessians[1])
	+ (2.0*x2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] +  2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,0,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*x2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,1,1)] = (-2.0*x2*y1 + 4.0*x1*y2)*gradients[2] + (4.0*x2*y1)*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*y2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1* z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,1,1)] = (-2.0*x1*y1)*gradients[2] + (4.0*x1*y1 + 8.0*x2*y2)*gradients[1] 
	+  (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*y2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,1,1)] = (-2.0*x1*x2 - 2.0*z1*z2)*gradients[2] + (4.0*y1*y2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*y2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,1,1)] = (2.0*x1_2 + 2.0*z1_2)*gradients[2] + (4.0*y1_2 + 8.0*y2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.*gradients[0]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*y2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,1,1)] = (4.0*y2*z1 - 2.0*y1*z2)*gradients[2] + (4.0*y1*z2)*gradients[1] + 
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*y2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,1,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
	+ (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*y2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	DPDF[compute6x6TensorIndex(0,0,2,1)] = (-2.0*x2*z1 + 4.0*x1*z2)*gradients[2] + (4.0*x2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[5] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*x1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[4] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*x1*hessians[1])
	+ (2.0*z2)*((-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*hessians[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*x1*hessians[0]);

	DPDF[compute6x6TensorIndex(0,1,2,1)] = (-2.0*x1*z1)*gradients[2] + (4.0*x1*z1 + 8.0*x2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[5] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*x2*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[4] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*x2*hessians[1])
	+ (2.0*z2)*((2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*hessians[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*x2*hessians[0]);

	DPDF[compute6x6TensorIndex(1,0,2,1)] = (-2.0*y2*z1 + 4.0*y1*z2)*gradients[2] + (4.0*y2*z1)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[5] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*y1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[4] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*y1*hessians[1])
	+ (2.0*z2)*((2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*hessians[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*y1*hessians[0]);

	DPDF[compute6x6TensorIndex(1,1,2,1)] = (-2.0*y1*z1)*gradients[2] + (4.0*y1*z1 + 8.0*y2*z2)*gradients[1]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[5] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*y2*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[4] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*y2*hessians[1])
	+ (2.0*z2)*((-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*hessians[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*y2*hessians[0]);

	DPDF[compute6x6TensorIndex(2,0,2,1)] = (-2.0*x1*x2 - 2.0*y1*y2)*gradients[2] + (4.0*z1*z2 + 4.0*(x1*x2 + y1*y2 + z1*z2))*gradients[1] 
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[5] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[4] + 2.0*z1*hessians[2])
	+ (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[4] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[3] + 2.0*z1*hessians[1])
	+ (2.0*z2)*((2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*hessians[2] + (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*hessians[1] + 2.0*z1*hessians[0]);

	DPDF[compute6x6TensorIndex(2,1,2,1)] = (2.0*x1_2 + 2.0*y1_2)*gradients[2] + (4.0*z1_2 + 8.0*z2_2 + 4.0*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*gradients[0]
	+ (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[5] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[4] + 2.0*z2*hessians[2])
  + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[4] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[3] + 2.0*z2*hessians[1])
	+ (2.0*z2)*((-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*hessians[2] + (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*hessians[1] + 2.0*z2*hessians[0]);

	//DBG__Matrix_Dump(&DPDF[0],6,6);
}

void ImplicitHyperElasticFEMSolver::computeFirstPiolaStress(Eigen::MatrixXd &F, std::vector<double> &gradients, Eigen::MatrixXd &P)
{
	double x1 = F(0,0); double x2 = F(0,1);
  double y1 = F(1,0); double y2 = F(1,1);
  double z1 = F(2,0); double z2 = F(2,1);

	double x1_2 = x1*x1; double x2_2 = x2*x2;
	double y1_2 = y1*y1; double y2_2 = y2*y2;
	double z1_2 = z1*z1; double z2_2 = z2*z2;

	P = Eigen::MatrixXd(3,2);

	P(0,0) = (-2.0*x2*y1*y2 + 2.0*x1*y2_2 - 2.0*x2*z1*z2 + 2.0*x1*z2_2)*gradients[2] + (4.0*x1*(x1_2 + y1_2 + z1_2) + 4.0*x2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] +  2.0*x1*gradients[0];
	P(0,1) = (2.0*x2*y1_2 - 2.0*x1*y1*y2 + 2.0*x2*z1_2 - 2.0*x1*z1*z2)*gradients[2] + (4.0*x1*(x1*x2 + y1*y2 + z1*z2) + 4.0*x2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*x2*gradients[0];

	P(1,0) = (2.0*x2_2*y1 - 2.0*x1*x2*y2 - 2.0*y2*z1*z2 + 2.0*y1*z2_2)*gradients[2] + (4.0*y1*(x1_2 + y1_2 + z1_2) + 4.0*y2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] + 2.0*y1*gradients[0];
	P(1,1) = (-2.0*x1*x2*y1 + 2.0*x1_2*y2 + 2.0*y2*z1_2 - 2.0*y1*z1*z2)*gradients[2] + (4.0*y1*(x1*x2 + y1*y2 + z1*z2) + 4.0*y2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*y2*gradients[0];

	P(2,0) = (2.0*x2_2*z1 + 2.0*y2_2*z1 - 2.0*x1*x2*z2 - 2.0*y1*y2*z2)*gradients[2] +  (4.0*z1*(x1_2 + y1_2 + z1_2) + 4.0*z2*(x1*x2 + y1*y2 + z1*z2))*gradients[1] + 2.0*z1*gradients[0];
	P(2,1) = (-2.0*x1*x2*z1 - 2.0*y1*y2*z1 + 2.0*x1_2*z2 + 2.0*y1_2*z2)*gradients[2] +  (4.0*z1*(x1*x2 + y1*y2 + z1*z2) + 4.0*z2*(x2_2 + y2_2 + z2_2))*gradients[1] + 2.0*z2*gradients[0];
}


void ImplicitHyperElasticFEMSolver::computeDPDF(std::vector<double> &DPDF_Hat, Eigen::MatrixXd &U, Eigen::MatrixXd &V,std::vector<double> &DPDF)
{
	Eigen::MatrixXd UT = U.transpose();
	Eigen::MatrixXd VT = V.transpose();

	std::vector<double> eiejVector(6,0);
	DPDF.resize(36,0);

	for (int column = 0; column < 6; column++) {
    eiejVector[column] = 1.0;
		Eigen::MatrixXd ei_ej = convert6VectorToEigen3x2Matrix(eiejVector);
    Eigen::MatrixXd ut_eiej_v = UT * ei_ej * V;

    /*std::vector<double> ut_eiej_v_TeranVector(4); //in Teran order
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v(0,0);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v(0,1);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v(1,0);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v(1,1);*/

		/*std::vector<double> ut_eiej_v_Vector(4);
		ut_eiej_v_Vector[0] = ut_eiej_v(0,0);
    ut_eiej_v_Vector[1] = ut_eiej_v(0,1);
    ut_eiej_v_Vector[2] = ut_eiej_v(1,0);
    ut_eiej_v_Vector[3] = ut_eiej_v(1,1);
    
		std::vector<double> dPdF_resultVector(6); // not in Teran order
    for (int innerRow = 0; innerRow < 6; innerRow++) {
      double tempResult = 0.0;
      for (int innerColumn = 0; innerColumn < 4; innerColumn++) {
        tempResult += DPDF_Hat[innerRow * 6 + innerColumn] *
                      ut_eiej_v_Vector[innerColumn];
      }
      //dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
      dPdF_resultVector[innerRow] = tempResult;
    }*/

		Eigen::MatrixXd tmpMat(3,2);
		int m = column/2;
		int n = column%2;
		tmpMat(0,0) = DPDF_Hat[compute6x6TensorIndex(0,0,m,n)];
		tmpMat(0,1) = DPDF_Hat[compute6x6TensorIndex(0,1,m,n)];
		tmpMat(1,0) = DPDF_Hat[compute6x6TensorIndex(1,0,m,n)];
		tmpMat(1,1) = DPDF_Hat[compute6x6TensorIndex(1,1,m,n)];
		tmpMat(2,0) = DPDF_Hat[compute6x6TensorIndex(2,0,m,n)];
		tmpMat(2,1) = DPDF_Hat[compute6x6TensorIndex(2,1,m,n)];

		Eigen::MatrixXd dPdF_resultMatrix_6 = tmpMat * ut_eiej_v;
		Eigen::MatrixXd dPdF_resultMatrix(2,2);
		dPdF_resultMatrix(0,0) = dPdF_resultMatrix_6(0,0);
		dPdF_resultMatrix(0,1) = dPdF_resultMatrix_6(0,1);
		dPdF_resultMatrix(1,0) = dPdF_resultMatrix_6(1,0);
		dPdF_resultMatrix(1,1) = dPdF_resultMatrix_6(1,1);
		
		//Eigen::MatrixXd dPdF_resultMatrix = convert4VectorToEigen2x2Matrix(dPdF_resultVector);
    Eigen::MatrixXd u_dpdf_vt = U * dPdF_resultMatrix * VT;
    /*DPDF[column +  0] = u_dpdf_vt(0,0);
    DPDF[column +  6] = u_dpdf_vt(0,1);
    DPDF[column + 12] = u_dpdf_vt(1,0);
    DPDF[column + 18] = u_dpdf_vt(1,1);
    DPDF[column + 24] = u_dpdf_vt(2,0);
    DPDF[column + 30] = u_dpdf_vt(2,1);*/

		DPDF[compute6x6TensorIndex(0,0,m,n)] = u_dpdf_vt(0,0);
    DPDF[compute6x6TensorIndex(0,1,m,n)] = u_dpdf_vt(0,1);
    DPDF[compute6x6TensorIndex(1,0,m,n)] = u_dpdf_vt(1,0);
    DPDF[compute6x6TensorIndex(1,1,m,n)] = u_dpdf_vt(1,1);
    DPDF[compute6x6TensorIndex(2,0,m,n)] = u_dpdf_vt(2,0);
    DPDF[compute6x6TensorIndex(2,1,m,n)] = u_dpdf_vt(2,1);

    // reset
    eiejVector[column] = 0.0;
  }

}

void ImplicitHyperElasticFEMSolver::computeDGDF(std::vector<double> &DPDF, Eigen::Matrix2d bVec, std::vector<double> &DGDF)
{
	DGDF.resize(36,0);
	for (int abc = 0; abc < 2; abc++)
		for (int i = 0; i < 3; i++)
			for (int column = 0; column < 6; column++)
				for (int k = 0; k < 2; k++)
				{ 
					DGDF[18 * abc + 6 * i + column] += DPDF[(2 * i + k) * 6 + column] * (bVec(abc,k));
				}
				//std::cout << "\n";
	//DBG__Matrix_Dump(&DGDF[0],6,6);
}

void ImplicitHyperElasticFEMSolver::computeElementStiffnessMatrix(int el,std::vector<double> &DGDF,std::vector<double> &KELEM)
{
	double * dFdU = &dFdUs[54 * el];
	//DBG__Matrix_Dump(dFdU,6,9);
		// K is stored column-major (however, it doesn't matter because K is symmetric)
		// Now its time to initalize element K which is a 9x9 matrix  arranged as a double vector 
		KELEM.resize(81,0);
		for (int row = 0; row < 6; row++) {
			for (int column = 0; column < 9; column++) {
				double result = 0;
				for (int inner = 0; inner < 6; inner++) {
					//dGdF is 6x6, and dFdU is 6x9
					result += DGDF[6 * row + inner] * dFdU[9 * inner + column];
				}
				KELEM[9 * column + row] = result;
			}
		}

		for (int row = 0; row < 9; row++) {
    //17th column
    KELEM[9 * row + 6] = -KELEM[9 * row + 0] - KELEM[9 * row + 3] ;
    //8th column
    KELEM[9 * row + 7] = -KELEM[9 * row + 1] - KELEM[9 * row + 4] ;
    //9th column
    KELEM[9 * row + 8] = -KELEM[9 * row + 2] - KELEM[9 * row + 5] ;
  }

	/*for(int i=0;i<9;i++)
	{
		for(int j=0;j<9;j++)
		{
			std::cout << std::setw(4) << KELEM[9*j+i] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "------------\n";*/
}

void ImplicitHyperElasticFEMSolver::addElementStiffnessMatrixToGlobalStiffnessMatrix(Cloth_Data* cloth,int el,std::vector<double> &KELEM)
{
  Triangles tri = cloth->getMesh()->get_triangle(el);
  std::vector<int> vertices(3);
	vertices[0] = tri.a;
	vertices[1] = tri.b;
	vertices[2] = tri.c;
	
	for (int vtxIndexA = 0; vtxIndexA < 3; vtxIndexA++)
	{
		for (int vtxIndexB = 0; vtxIndexB < 3; vtxIndexB++) {
			int vtxA = vertices[vtxIndexA];
			int vtxB = vertices[vtxIndexB];
			for (int i = 0; i < 3; i++) 
			{
				for (int j = 0; j < 3; j++) 
				{
					int row = 3 * vtxA + i;
					double * value = &KELEM[ELT(9, 3 * vtxIndexA + i, 3 * vtxIndexB + j)];
					int columnIndex = tangentStiffnessMatrix_->GetInverseIndex(row, 3 * vtxB + j);
					tangentStiffnessMatrix_->AddEntry(row, columnIndex, *value);
				}
			}
		}
	}
}

void ImplicitHyperElasticFEMSolver::addShearComponents(Cloth_Data *cloth)
{
	float young_modulus = cloth->get_property_obj()->get_youngs_modulus();
	float poisson = cloth->get_property_obj()->get_poisson_ratio();

	float lame_lambda = (young_modulus*poisson)/((1.0+poisson)*(1.0-(2.0*poisson)));
	float lame_mu = young_modulus/(2.0*(1.0+poisson));

	float h = cloth->get_property_obj()->get_timestep();

	int num_vertices = cloth->getMesh()->get_number_vertices();
	int num_triangles = cloth->getMesh()->get_number_triangles();
	
	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = cloth->getMesh()->get_triangle(t);

		Eigen::Vector3d pa = cloth->getMesh()->get_point_data( tri.a,lastFrameId_);
		Eigen::Vector3d pb = cloth->getMesh()->get_point_data( tri.b,lastFrameId_);
		Eigen::Vector3d pc = cloth->getMesh()->get_point_data( tri.c,lastFrameId_);

		Eigen::Vector3d va = cloth->get_velocity(tri.a);
		Eigen::Vector3d vb = cloth->get_velocity(tri.b);
		Eigen::Vector3d vc = cloth->get_velocity(tri.c);

		Eigen::Matrix2d Bm = cloth->getDmInv(t);
		Eigen::Matrix2d BmT = Bm.transpose();
		/*double W = cloth->getW(t);

		float r_ua = (cloth->get_vertex_distribution(t))(0);
		float r_va = (cloth->get_vertex_distribution(t))(1);
		float r_ub = (cloth->get_vertex_distribution(t))(2);
		float r_vb = (cloth->get_vertex_distribution(t))(3);
		float r_uc = (cloth->get_vertex_distribution(t))(4);
		float r_vc = (cloth->get_vertex_distribution(t))(5);
		float d = (cloth->get_vertex_distribution(t))(6);

		Eigen::Vector3d U = (r_ua * pa) + (r_ub * pb) + (r_uc * pc);
		Eigen::Vector3d V = (r_va * pa) + (r_vb * pb) + (r_vc * pc);

		Eigen::MatrixXd F(3,2);
		F(0,0) = U(0); F(0,1) = V(0);
		F(1,0) = U(1); F(1,1) = V(1);
		F(2,0) = U(2); F(2,1) = V(2);*/

		Eigen::Vector3d diff1 = pc-pa;
		Eigen::Vector3d diff2 = pc-pb;

		Eigen::MatrixXd tmp(3,2);
		tmp(0,0) = diff1.x(); tmp(0,1) = diff2.x();
		tmp(1,0) = diff1.y(); tmp(1,1) = diff2.y();
		tmp(2,0) = diff1.z(); tmp(2,1) = diff2.z();

		Eigen::MatrixXd F = tmp*Bm;

		//Assume here that we are using NeoHookean Materials 
		std::vector<double> gradients(3);
		std::vector<double> hessians(6);

		double IIIC = (F(0,1)*F(0,1)*F(1,0)*F(1,0)) - (2.0*F(0,0)*F(0,1)*F(1,0)*F(1,1)) + (F(0,0)*F(0,0)*F(1,1)*F(1,1))+ (F(0,1)*F(0,1)*F(2,0)*F(2,0)) + (F(1,1)*F(1,1)*F(2,0)*F(2,0)) - (2.0*F(0,0)*F(0,1)*F(2,0)*F(2,1)) - (2.0*F(1,0)*F(1,1)*F(2,0)*F(2,1)) + (F(0,0)*F(0,0)*F(2,1)*F(2,1)) + (F(1,0)*F(1,0)*F(2,1)*F(2,1));


#ifdef NEOHOOKEAN_MATERIAL 
		gradients[0] = 0.5 * lame_mu;
		gradients[1] = 0.0;
		gradients[2] = (-0.5 * lame_mu + 0.25 * lame_lambda * log(IIIC)) / IIIC;

		hessians[0] = 0.0;
		hessians[1] = 0.0;
		hessians[2] = 0.0;
		hessians[3] = 0.0;
		hessians[4] = 0.0;
		hessians[5] = (0.25 * lame_lambda + 0.5 * lame_mu - 0.25 * lame_lambda * log(IIIC)) / (IIIC * IIIC);

#else
		gradients[0] = 0.25 * lame_lambda * (IC - 3.0) - 0.5 * lame_mu;
		gradients[1] = 0.25 * lame_mu;
		gradients[2] = 0.0;

		hessians[0] = 0.25 * lame_lambda;
		hessians[1] = 0.0;
		hessians[2] = 0.0;
		hessians[3] = 0.0;
		hessians[4] = 0.0;
		hessians[5] = 0.0;
#endif		
		

//Add compression resistance gradient if needed 
#ifdef ADDCOMPRESSIONRESISTENCE
		double J = sqrt(IIIC);
		if (J < 1.0)
		{
			double compressionResistanceFactor =1.0;
			gradients[2] += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
			hessians[5] += compressionResistanceFactor * (1.0 - J) * (1.0 + J) / (3456.0 * J * J * J);
		}
#endif

		//Directly compute the piola kirchoff stress tensor 
		Eigen::MatrixXd P_1;
		computeFirstPiolaStress(F,gradients,P_1);

		Eigen::Vector2d ewan0 = cloth->getEdgeWeightedTriangleNormals(t,0);
		Eigen::Vector2d ewan1 = cloth->getEdgeWeightedTriangleNormals(t,1);
		Eigen::Vector2d ewan2 = cloth->getEdgeWeightedTriangleNormals(t,2);

		Eigen::MatrixXd internalForce_0 = P_1 * ewan0;
		Eigen::MatrixXd internalForce_1 = P_1 * ewan1;
		Eigen::MatrixXd internalForce_2 = P_1 * ewan2;

		force_[3*tri.a + 0] -= internalForce_0(0,0);
		force_[3*tri.a + 1] -= internalForce_0(1,0);
		force_[3*tri.a + 2] -= internalForce_0(2,0);

		force_[3*tri.b + 0] -= internalForce_1(0,0);
		force_[3*tri.b + 1] -= internalForce_1(1,0);
		force_[3*tri.b + 2] -= internalForce_1(2,0);

		force_[3*tri.c + 0] -= internalForce_2(0,0);
		force_[3*tri.c + 1] -= internalForce_2(1,0);
		force_[3*tri.c + 2] -= internalForce_2(2,0);

		internalForce_[3*tri.a + 0] += internalForce_0(0,0);
		internalForce_[3*tri.a + 1] += internalForce_0(1,0);
		internalForce_[3*tri.a + 2] += internalForce_0(2,0);

		internalForce_[3*tri.b + 0] += internalForce_1(0,0);
		internalForce_[3*tri.b + 1] += internalForce_1(1,0);
		internalForce_[3*tri.b + 2] += internalForce_1(2,0);

		internalForce_[3*tri.c + 0] += internalForce_2(0,0);
		internalForce_[3*tri.c + 1] += internalForce_2(1,0);
		internalForce_[3*tri.c + 2] += internalForce_2(2,0);

		//Now we will compute the stiffness matrices one by one 
		
		//First we need to compute the value of dPdF_atFHat which is a 4x4 tensor
		std::vector<double> dPdF;
		std::vector<double> dGdF;
		std::vector<double> KElem;

		//Next we compute teh dPdF which is a 6x6 tensor
		computeDPDFFrustrating(F,gradients,hessians,dPdF);
		
		//Next we compute the dGdF which is a 6x6 tensor
		Eigen::Matrix2d bVec;
		bVec(0,0) = ewan0(0);
		bVec(0,1) = ewan0(1);
		bVec(1,0) = ewan1(0);
		bVec(1,1) = ewan1(1);
		computeDGDF(dPdF,bVec,dGdF);
		//Finally using DGDF and precomputede dFdUs we can compute the elementary stiffness matrix 
		computeElementStiffnessMatrix(t,dGdF,KElem);
		//Insert it in the correct place in the global matrix
		addElementStiffnessMatrixToGlobalStiffnessMatrix(cloth,t,KElem);

		if(t==1000)
		{
			std::cout << F << "\n";

			std::cout << "--------------\n";

			std::cout << lame_lambda << "  " << lame_mu << "\n";

				std::cout << "--------------\n";

				for(int m=0;m<3;m++)
				{
					for(int n=0;n<2;n++)
					{
						for(int i=0;i<3;i++)
						{
							for(int j=0;j<2;j++)
							{
								std::cout << dPdF[compute6x6TensorIndex(i,j,m,n)] << " ";
							}
							std::cout << "\n";
						}
					}
				}

			std::cout << "++++++++++++++\n";
		}
	}
}

void ImplicitHyperElasticFEMSolver::addGravityComponents(Cloth_Data *cloth) {
	int num_vertices = cloth->getMesh()->get_number_vertices();
	//OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float gravForce = cloth->get_vertex_mass(p) * 9.8f; 
		externalForce_[3*p+1] -= gravForce;
		force_[3*p+1] -= gravForce;
	}

}

void ImplicitHyperElasticFEMSolver::addLinearDamping(Cloth_Data* cloth)
{
	float kd = 0.001;
	int num_vertices = cloth->getMesh()->get_number_vertices();
	//OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
	
	  Eigen::Vector3d velocity = cloth->get_velocity(p);
		Eigen::Vector3d dampingForce = -kd*velocity;
		force_[3*p+0] += dampingForce[0]; 
		force_[3*p+1] += dampingForce[1]; 
		force_[3*p+2] += dampingForce[2]; 
	}
}

void ImplicitHyperElasticFEMSolver::importPreviousStepVelocities(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	//OMP_FOR
	for(int i=0;i<numVertices;i++)
	{
	  Eigen::Vector3d v = cloth->get_velocity(i);
		velocities_[3*i+0] = v[0];
		velocities_[3*i+1] = v[1];
		velocities_[3*i+2] = v[2];
	}
}

void ImplicitHyperElasticFEMSolver::computeFinalLinearSystem(Cloth_Data* cloth)
{
	float timeStep = cloth->get_property_obj()->get_timestep();
	int num_vertices = cloth->getMesh()->get_number_vertices();

	double *tempMatVecPdk = new double[3*num_vertices];

#ifdef ENABLE_RAYLEIGH_DAMPING
	double dampingStiffness = 0.01;
	double dampingMass = 1.0;

	tangentStiffnessMatrix_->ScalarMultiply(dampingStiffness, rayleighDampingMatrix_);
  rayleighDampingMatrix_->AddSubMatrix(dampingMass,*massMatrix_);
#endif

	*tangentStiffnessMatrix_ *= timeStep;
  tangentStiffnessMatrix_->MultiplyVectorAdd(velocities_.data(), tempMatVecPdk);

#ifdef ENABLE_RAYLEIGH_DAMPING
	*tangentStiffnessMatrix_ += *rayleighDampingMatrix_;
#endif

  *tangentStiffnessMatrix_ *= timeStep;
  tangentStiffnessMatrix_->AddSubMatrix(1.0, *massMatrix_);

	// add externalForces, internalForces
  //OMP_FOR 
	for(int i=0; i<(3*num_vertices); i++)
  {
	  residuals_[i] += (externalForce_[i]-internalForce_[i] - tempMatVecPdk[i])*timeStep;
  }
	delete[] tempMatVecPdk;

}

void ImplicitHyperElasticFEMSolver::finalizeAndSolve(Cloth_Data *cloth) {
	float timeStep = cloth->get_property_obj()->get_timestep();
	int num_vertices = cloth->getMesh()->get_number_vertices();
	
	//Add the UI Force if any 
	int lastVClicked = cloth->getLastClickedVertex();
	Eigen::Vector3d uiForce = cloth->getCurrentUIForce();
	if(lastVClicked!=-1) {
		externalForce_[3*lastVClicked+0] += uiForce[0];
		externalForce_[3*lastVClicked+1] += uiForce[1];
		externalForce_[3*lastVClicked+2] += uiForce[2];
	}
	
	//This is the stuff for theimplicit method (these are the finalization steps)
	importPreviousStepVelocities(cloth);
	computeFinalLinearSystem(cloth);
	
	

	//Setup the constrained version of the whole thingy 
	int n3 = num_vertices*3;
	int numConstrainedDOFs = n3 - (3*numConstrainedVerts_);
	//Make a copy of the LHSMatrix_ and remove rows and columns from it
	SparseMatrix *tempSpMatCopy = new SparseMatrix(*tangentStiffnessMatrix_);
	tempSpMatCopy->RemoveRowsColumns(numConstrainedVerts_*3,constrainedVerts_);
	double sum = tempSpMatCopy->SumEntries();
	std::cout << "Checksum:" << sum << "\n";
	//Make a copy of RHS Vector and remove rows and columns from it 
	double* rhsConstrained = new double[numConstrainedDOFs];
	RemoveRows(num_vertices*3,rhsConstrained,residuals_.data(),numConstrainedVerts_*3,constrainedVerts_);
	//Make a correct size row vector to hold the reult
	double *resultConstrained = new double[numConstrainedDOFs];
	memset(resultConstrained,0,sizeof(double)*numConstrainedDOFs);

#if defined USE_PARDISO_SOLVER
	PardisoSolver solver(tempSpMatCopy,7,0,0,0);
	solver.ComputeCholeskyDecomposition(tempSpMatCopy);
	int retVal = solver.SolveLinearSystem(resultConstrained,rhsConstrained);
#elif defined USE_CG_SOLVER
	CGSolver cgsolver(tempSpMatCopy);
	cgsolver.SolveLinearSystemWithJacobiPreconditioner(resultConstrained, rhsConstrained, 1e-6, 10000);
#else
	int numGaussSeidelIterations = 8000;
		for(int iter=0; iter<numGaussSeidelIterations; iter++)
			tempSpMatCopy->DoOneGaussSeidelIteration(resultConstrained, rhsConstrained);
#endif
	//Expand the rows now to fill it 
	double *delV = new double[num_vertices*3];
	InsertRows(num_vertices*3,resultConstrained,delV,numConstrainedVerts_*3,constrainedVerts_);
	//Free the stuff we allocated
	delete(tempSpMatCopy);
	delete[] rhsConstrained;
	delete[] resultConstrained;

	//For now just leapfrog this one 
	//OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		Eigen::Vector3d oldPos = cloth->getMesh()->get_point_data(p,lastFrameId_);
		Eigen::Vector3d oldVel(velocities_[3*p+0],velocities_[3*p+1],velocities_[3*p+2]);;

		//Eigen::Vector3d elemDelQ(delV[3*p+0],delV[3*p+1],delV[3*p+2]);
		Eigen::Vector3d elemDelV(delV[3*p+0],delV[3*p+1],delV[3*p+2]);
		Eigen::Vector3d newVel = oldVel + elemDelV;
		Eigen::Vector3d newPos = oldPos + (timeStep*elemDelV);
		cloth->set_next_step_pos(newPos,p);
		cloth->set_next_step_velocity(newVel,p);
	}

	delete[] delV;
}

void ImplicitHyperElasticFEMSolver::finalizeAndSolve_Explicit(Cloth_Data* cloth)
{
	float time_step = cloth->get_property_obj()->get_timestep();
  int num_vertices = cloth->getMesh()->get_number_vertices();
 	
 	//Add the UI Force if any 
 	int lastVClicked = cloth->getLastClickedVertex();
 	Eigen::Vector3d uiForce = cloth->getCurrentUIForce();
 	if(lastVClicked!=-1) {
 		force_[3*lastVClicked+0] += uiForce[0];
 		force_[3*lastVClicked+1] += uiForce[1];
 		force_[3*lastVClicked+2] += uiForce[2];
 	}
 
 	
 	OMP_FOR
 	for(int p=0;p<num_vertices;p++)
 	{
 		Eigen::Vector3d old_pos = cloth->getMesh()->get_point_data(p,lastFrameId_);
 		Eigen::Vector3d old_vel = cloth->get_velocity(p);
 
 		Eigen::Vector3d a(force_[3*p+0]/massMatrixDiagonal_[3*p+0],force_[3*p+1]/massMatrixDiagonal_[3*p+1],force_[3*p+2]/massMatrixDiagonal_[3*p+2]);
 
 		Eigen::Vector3d elemDelV = a*time_step;
		Eigen::Vector3d new_vel =old_vel +  elemDelV;
		Eigen::Vector3d new_pos = old_pos + (time_step * new_vel);
 
 		if(p==CORNERV1 || p==CORNERV2) {
 			new_pos = old_pos;
 			new_vel = old_vel;
 		}
 		cloth->set_next_step_pos(new_pos,p);
 		cloth->set_next_step_velocity(new_vel,p);
 	}
}

void ImplicitHyperElasticFEMSolver::resetParameters() {
	force_.setZero();
}