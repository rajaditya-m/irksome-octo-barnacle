#include "ImplicitHyperElasticFEMSolver.h"

#ifndef ELT
  #define ELT(numRows,i,j) (((long)j)*((long)numRows)+((long)i))
#endif

ImplicitHyperElasticFEMSolver::ImplicitHyperElasticFEMSolver(void)
{
  rowMajorMatrixToTeran.resize(4);
  rowMajorMatrixToTeran[0] = 0;
	rowMajorMatrixToTeran[1] = 3;
	rowMajorMatrixToTeran[2] = 1;
	rowMajorMatrixToTeran[3] = 2;

	teranToRowMajorMatrix.resize(4);
  rowMajorMatrixToTeran[0] = 0;
	rowMajorMatrixToTeran[1] = 2;
	rowMajorMatrixToTeran[2] = 3;
	rowMajorMatrixToTeran[3] = 1;

	 /*for (int abc = 0; abc < 3; abc++)
    for (int i = 0; i < 3; i++)
      for (int column = 0; column < 9; column++)
        for (int k = 0; k < 3; k++)
         std::cout << 27 * abc + 9 * i + column << "\n";//] += dPdF[(3 * i + k) * 9 + column] * (*(bVec[abc]))[k];*/

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
}


void ImplicitHyperElasticFEMSolver::computedFdU(Cloth_Data* cloth) {
	int numElements = cloth->getMesh()->get_number_triangles();
  for (int el = 0; el < numElements; el++) {
    double * dFdU = &dFdUs[54 * el];
		Eigen::Matrix2d dmInv = cloth->getDmInv(el);
    for (int index = 0; index < 36; index++) {
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
}

void ImplicitHyperElasticFEMSolver::initialize(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
	int numTriangles = cloth->getMesh()->get_number_triangles();

	massMatrixDiagonal_ = new double[3*numVertices];
	OMP_FOR
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
	constrainedVerts_[3] = 7800;
	constrainedVerts_[4] = 7801;
	constrainedVerts_[5] = 7802;
	constrainedVerts_[0] = 150;
	constrainedVerts_[1] = 151;
	constrainedVerts_[2] = 152;

	//Add the complemented stuff 
	dFdUs.resize(54*numTriangles);
	computedFdU(cloth);

	//Compute the tangent Stiffness Matrix Structure 
	initializeSparseMatrixFromOutline(cloth);

}

void ImplicitHyperElasticFEMSolver::advance_time_step(Cloth_Data* cloth)
{
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;
	
	//Preapre the components as usual 
	tangentStiffnessMatrix_->ResetToZero();

	//Prepare the RHS Vector for usage
	force_.setZero();

	//Add the physics components
	addShearComponents(cloth);
	//addLinearDamping(cloth);
	//addShearComponentsCorotational(cloth);
	addGravityComponents(cloth);

	//Solve and report
	finalizeAndSolve(cloth);
}

double ImplicitHyperElasticFEMSolver::computeGammaValues(int i,int j,std::vector<double> &sigma,double IIIC, std::vector<double> &gradient,std::vector<double> &hessian)
{
	double tempGammaVec1[3];
  tempGammaVec1[0] = 2.0 * sigma[i];
  tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
  tempGammaVec1[2] = 2.0 * IIIC / sigma[i];
  double tempGammaVec2[3];
  tempGammaVec2[0] = 2.0 * sigma[j];
  tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
  tempGammaVec2[2] = 2.0 * IIIC / sigma[j];
  double productResult[3];
  productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] +
                      tempGammaVec2[2] * hessian[2]);
  productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] +
                      tempGammaVec2[2] * hessian[4]);
  productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] +
                      tempGammaVec2[2] * hessian[5]);
  return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
          tempGammaVec1[2] * productResult[2] + 4.0 * IIIC * gradient[2] / (sigma[i] * sigma[j]));
}

int ImplicitHyperElasticFEMSolver::compute4x4TensorIndex(int i, int j, int m, int n)
{
	int rowIndex_in9x9Matrix = rowMajorMatrixToTeran[2 * i + j];
  int columnIndex_in9x9Matrix = rowMajorMatrixToTeran[2 * m + n];
  return (4 * rowIndex_in9x9Matrix + columnIndex_in9x9Matrix);
}

int ImplicitHyperElasticFEMSolver::compute4x9TensorIndex(int i, int j, int m, int n) {
  int rowIndex = 3 * i + j;
  int columnIndex = 2 * m + n;
  return (9 * rowIndex + columnIndex);
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

void ImplicitHyperElasticFEMSolver::computeDPDF_Hat(std::vector<double> &sigma, std::vector<double> &gradients, std::vector<double> &hessians, double IIIC,std::vector<double> &DPDH_Hat)
{
  double l1_Sq = sigma[0]*sigma[0];
  double l2_Sq = sigma[0]*sigma[0];

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
	x2121 = alpha12;
	x2112 = beta12;

	DPDH_Hat.resize(16,0);

	DPDH_Hat[compute4x4TensorIndex(0, 0, 0, 0)] = x1111;
	DPDH_Hat[compute4x4TensorIndex(0, 0, 1, 1)] = x2211;

	DPDH_Hat[compute4x4TensorIndex(1, 1, 0, 0)] = x2211;
	DPDH_Hat[compute4x4TensorIndex(1, 1, 1, 1)] = x2222;

	DPDH_Hat[compute4x4TensorIndex(0, 1, 0, 1)] = x2121;
	DPDH_Hat[compute4x4TensorIndex(0, 1, 1, 0)] = x2112;

	DPDH_Hat[compute4x4TensorIndex(1, 0, 0, 1)] = x2112;
	DPDH_Hat[compute4x4TensorIndex(1, 0, 1, 0)] = x2121;
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

    std::vector<double> ut_eiej_v_TeranVector(4); //in Teran order
    
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v(0,0);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v(0,1);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v(1,0);
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v(1,1);
    
		std::vector<double> dPdF_resultVector(4); // not in Teran order
    for (int innerRow = 0; innerRow < 4; innerRow++) {
      double tempResult = 0.0;
      for (int innerColumn = 0; innerColumn < 4; innerColumn++) {
        tempResult += DPDF_Hat[innerRow * 4 + innerColumn] *
                      ut_eiej_v_TeranVector[innerColumn];
      }
      dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
    }

		Eigen::MatrixXd dPdF_resultMatrix = convert4VectorToEigen2x2Matrix(dPdF_resultVector);
    Eigen::MatrixXd u_dpdf_vt = U * dPdF_resultMatrix * VT;
    DPDF[column +  0] = u_dpdf_vt(0,0);
    DPDF[column +  6] = u_dpdf_vt(0,1);
    DPDF[column + 12] = u_dpdf_vt(1,0);
    DPDF[column + 18] = u_dpdf_vt(1,1);
    DPDF[column + 24] = u_dpdf_vt(2,0);
    DPDF[column + 30] = u_dpdf_vt(2,1);
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
					  DGDF[18 * abc + 6 * i + column] += DPDF[(2 * i + k) * 6 + column] * (bVec(abc,k));
}

void ImplicitHyperElasticFEMSolver::computeElementStiffnessMatrix(int el,std::vector<double> &DGDF,std::vector<double> &KELEM)
{
	double * dFdU = &dFdUs[54 * el];

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
	
	OMP_FOR
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
		double W = cloth->getW(t);

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
		F(2,0) = U(2); F(2,1) = V(2);

		
		//Compute the singular value decompositions now
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(F,Eigen::ComputeThinU|Eigen::ComputeThinV);
		Eigen::VectorXd fHat = svd.singularValues();
		Eigen::MatrixXd UMat = svd.matrixU();
		Eigen::MatrixXd VMat = svd.matrixV();

		//@TODO:Compute the inversion threshold factor here

		//Now compute the value of the first piola stress from the force vectors (rather only the principal stretches)
		std::vector<double> sigma(2);
		sigma[0] = fHat(0);
		sigma[1] = fHat(1);
		double l1_Sq = fHat(0)*fHat(0);
		double l2_Sq = fHat(1)*fHat(1);

		double IC = l1_Sq + l2_Sq;
		double IIC = l1_Sq*l1_Sq + l2_Sq*l2_Sq;
		double IIIC = l1_Sq*l2_Sq;

		//Assume here that we are using NeoHookean Materials 
		std::vector<double> gradients(3);
		std::vector<double> hessians(6);

		//@TODO:Make this from the APIs later on

		gradients[0] = 0.5 * lame_mu;
		gradients[1] = 0.0;
		gradients[2] = (-0.5 * lame_mu + 0.25 * lame_lambda * log(IIIC)) / IIIC;

		hessians[0] = 0.0;
		hessians[1] = 0.0;
		hessians[2] = 0.0;
		hessians[3] = 0.0;
		hessians[4] = 0.0;
		hessians[5] = (0.25 * lame_lambda + 0.5 * lame_mu - 0.25 * lame_lambda * log(IIIC)) / (IIIC * IIIC);

		//Add compression resistance gradient if needed 
		if (ADDCOMPRESSIONRESISTENCE)
		{
			double J = sqrt(IIIC);
			if (J < 1.0)
			{
				double compressionResistanceFactor =1.0;
				gradients[2] += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
			}
		}

		//Now compute the values of P_Hat
		double P_Hat_0 = gradients[0]*2.0*fHat(0) + gradients[1]*4.0*l1_Sq*fHat(0) + gradients[2]*2.0*IIIC*(1.0/fHat(0));
		double P_Hat_1 = gradients[0]*2.0*fHat(1) + gradients[1]*4.0*l2_Sq*fHat(1) + gradients[2]*2.0*IIIC*(1.0/fHat(1));

		Eigen::Matrix2d P_Hat;
		P_Hat.setZero();
		P_Hat(0,0) = P_Hat_0;
		P_Hat(1,1) = P_Hat_1;
		Eigen::MatrixXd P_1 = UMat*P_Hat*(VMat.transpose());

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

		//Now we will compute the stiffness matrices one by one 
		
		//First we need to compute the value of dPdF_atFHat which is a 4x4 tensor
		std::vector<double> dPdF_atFHat;
		std::vector<double> dPdF;
		std::vector<double> dGdF;
		std::vector<double> KElem;

		computeDPDF_Hat(sigma,gradients,hessians,IIIC,dPdF_atFHat);
		//Next we compute teh dPdF which is a 6x6 tensor
		computeDPDF(dPdF_atFHat,UMat,VMat,dPdF);
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



	}
}

void ImplicitHyperElasticFEMSolver::addGravityComponents(Cloth_Data *cloth) {
	int num_vertices = cloth->getMesh()->get_number_vertices();
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
		float gravForce = cloth->get_vertex_mass(p) * 9.8f;
		force_[3*p+1] -= gravForce; 
	}

}

void ImplicitHyperElasticFEMSolver::addLinearDamping(Cloth_Data* cloth)
{
	float kd = 0.001;
	int num_vertices = cloth->getMesh()->get_number_vertices();
	OMP_FOR
	for(int p=0;p<num_vertices;p++)
	{
	
	  Eigen::Vector3d velocity = cloth->get_velocity(p);
		Eigen::Vector3d dampingForce = -kd*velocity;
		force_[3*p+0] += dampingForce[0]; 
		force_[3*p+1] += dampingForce[1]; 
		force_[3*p+2] += dampingForce[2]; 
	}
}

void ImplicitHyperElasticFEMSolver::finalizeAndSolve(Cloth_Data *cloth) {
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
		//DBG__Vector_Dump(elemDelV,p,"DElv");
		Eigen::Vector3d new_vel =old_vel +  elemDelV;
		Eigen::Vector3d new_pos = old_pos + (time_step * new_vel);

		if(p==2600 || p==50) {
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