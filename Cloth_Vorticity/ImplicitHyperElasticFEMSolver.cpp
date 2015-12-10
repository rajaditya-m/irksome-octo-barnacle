#include "ImplicitHyperElasticFEMSolver.h"


ImplicitHyperElasticFEMSolver::ImplicitHyperElasticFEMSolver(void)
	{
	}


ImplicitHyperElasticFEMSolver::~ImplicitHyperElasticFEMSolver(void)
	{
	}

void ImplicitHyperElasticFEMSolver::initialize(Cloth_Data* cloth)
{
	int numVertices = cloth->getMesh()->get_number_vertices();
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
}

void ImplicitHyperElasticFEMSolver::advance_time_step(Cloth_Data* cloth)
{
	//Set the last frame information
	lastFrameId_ = (cloth->getMesh()->get_number_frames()) - 1;
	
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
		double l1_Sq = fHat(0)*fHat(0);
		double l2_Sq = fHat(1)*fHat(1);
		double IC = l1_Sq + l2_Sq;
		double IIC = l1_Sq*l1_Sq + l2_Sq*l2_Sq;
		double IIIC = l1_Sq*l2_Sq;

		//Assume here that we are using NeoHookean Materials 
		double dPsidI1 = 0.5 * lame_mu;
		double dPsidI2 = 0.0;
		double dPsidI3 = (-0.5 * lame_mu + 0.25 * lame_lambda * log(IIIC)) / IIIC;

		//Add compression resistance gradient if needed 
		if (ADDCOMPRESSIONRESISTENCE)
		{
			double J = sqrt(IIIC);
			if (J < 1.0)
			{
				double compressionResistanceFactor =1.0;
				dPsidI3 += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
			}
		}

		//Now compute the values of P_Hat
		double P_Hat_0 = dPsidI1*2.0*fHat(0) + dPsidI2*4.0*l1_Sq*fHat(0) + dPsidI3*2.0*IIIC*(1.0/fHat(0));
		double P_Hat_1 = dPsidI1*2.0*fHat(1) + dPsidI2*4.0*l2_Sq*fHat(1) + dPsidI3*2.0*IIIC*(1.0/fHat(1));

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