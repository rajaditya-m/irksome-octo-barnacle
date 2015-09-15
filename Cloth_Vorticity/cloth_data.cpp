#include "cloth_data.h"


Cloth_Data::Cloth_Data( const char* name,
				const char* file_location_vertex,
				const char* file_location_mesh,
				const char* file_location_texture,
				const char* property_xml,
				const char* material_xml):
	RenderObject(CLOTH,false)
{
	//Initialize the mesh now 
	mesh_ = new TriMesh(name,file_location_vertex,file_location_mesh,file_location_texture,material_xml,1.0);

	//Set Cloth property
	cloth_prop_obj = new ClothProperties(property_xml);

	int num_vertices = mesh_->get_number_vertices();
	int num_triangles = mesh_->get_number_triangles();

	float density = cloth_prop_obj->get_density();
	float cloth_thickness = cloth_prop_obj->get_thickness();

	//Resize the velocity vectors and mass vector
	velocity.resize(num_vertices,Eigen::Vector3d::Zero());
	perVertexVectorBuffer_.resize(num_vertices,Eigen::Vector3d::Zero());
	vertex_mass.resize(num_vertices,0.0f);

	//Reset the buffers to contain all zeros
	next_step_position.resize(num_vertices,Eigen::Vector3d::Zero());
	next_step_velocity.resize(num_vertices,Eigen::Vector3d::Zero());

	for(int t=0;t<num_triangles;t++)
	{
		Triangles tri = mesh_->get_triangle(t);
		//Compute the mass of each triangle 
		Eigen::Vector3d p1 = mesh_->get_point_data(tri.a);
		Eigen::Vector3d p2 = mesh_->get_point_data(tri.b);
		Eigen::Vector3d p3 = mesh_->get_point_data(tri.c);
		float tri_area = triangle_area(p1,p2,p3);
		float mass = density*tri_area*cloth_thickness;
		vertex_mass[tri.a] += 0.33f * mass;
		vertex_mass[tri.b] += 0.33f * mass;
		vertex_mass[tri.c] += 0.33f * mass;

		//Precompute the UV Parameterization also 
		Eigen::Vector2d i0_t = mesh_->get_uv_data(tri.at);
		Eigen::Vector2d i1_t = mesh_->get_uv_data(tri.bt);
		Eigen::Vector2d i2_t = mesh_->get_uv_data(tri.ct);
	
		//Get the Isometric UV Parameterization computerd by rotation 
		Eigen::Vector3d triNorm = tri.getTriangleNormal();
		Eigen::Vector3d xzAxisNorm = Eigen::Vector3d(0.0,1.0,0.0);
		//Eigen::Matrix4f rotMatrix = compute_rotation_matrix_a_to_b(triNorm,xzAxisNorm);

		float ua = p1.x();
		float va = p1.z();
		float ub = p2.x();
		float vb = p2.z();
		float uc = p3.x();
		float vc = p3.z();


		//std::cout << "Triangle#" << t << " " << ua << "," << va << "," << ub << "," << vb << "," << uc << "," << vc << "\n";

		Vector7d temp_vertex_distribution;

		double det = ua*(vb - vc)+ ub*(vc - va)+ uc*(va - vb);
		double inv_det = 1.0f/det;
		temp_vertex_distribution(0) = (vb - vc) * inv_det;
		temp_vertex_distribution(1) = (uc - ub) * inv_det;
		temp_vertex_distribution(2) = (vc - va) * inv_det;
		temp_vertex_distribution(3) = (ua - uc) * inv_det;
		temp_vertex_distribution(4) = (va - vb) * inv_det;
		temp_vertex_distribution(5) = (ub - ua) * inv_det;
		temp_vertex_distribution(6) = det;
		vertex_distribution.push_back(temp_vertex_distribution);

		if(t==34)
			std::cout << "INIT# " << t << " "<< temp_vertex_distribution(0) << " " << temp_vertex_distribution(1) << " " << temp_vertex_distribution(2) << " " << temp_vertex_distribution(3) << " " << temp_vertex_distribution(4) << " " << temp_vertex_distribution(5) << "\n"; 

	}
}

	Cloth_Data::Cloth_Data( const char* name,
				const char* file_location,
				const char* property_xml,
				const char* material_xml):
	RenderObject(CLOTH,false)
{
	//Initialize the mesh now 
	mesh_ = new TriMesh(name,file_location,material_xml,1.0);

	//Set Cloth property
	cloth_prop_obj = new ClothProperties(property_xml);

	int num_vertices = mesh_->get_number_vertices();
	int num_triangles = mesh_->get_number_triangles();

	float density = cloth_prop_obj->get_density();
	float cloth_thickness = cloth_prop_obj->get_thickness();

	//Resize the velocity vectors and mass vector
	velocity.resize(num_vertices,Eigen::Vector3d::Zero());
	perVertexVectorBuffer_.resize(num_vertices,Eigen::Vector3d::Zero());
	vertex_mass.resize(num_vertices,0.0f);

	//Reset the buffers to contain all zeros
	next_step_position.resize(num_vertices,Eigen::Vector3d::Zero());
	next_step_velocity.resize(num_vertices,Eigen::Vector3d::Zero());

	for(int t=0;t<num_triangles;t++)
	{

		Triangles tri = mesh_->get_triangle(t);
		//Compute the mass of each triangle 
		Eigen::Vector3d p1 = mesh_->get_point_data(tri.a);
		Eigen::Vector3d p2 = mesh_->get_point_data(tri.b);
		Eigen::Vector3d p3 = mesh_->get_point_data(tri.c);


		float tri_area = triangle_area(p1,p2,p3);
		float mass = density*tri_area*cloth_thickness;
		vertex_mass[tri.a] += 0.33f * mass;
		vertex_mass[tri.b] += 0.33f * mass;
		vertex_mass[tri.c] += 0.33f * mass;

		//Precompute the UV Parameterization also 
		Eigen::Vector2d i0_t = mesh_->get_uv_data(tri.at);
		Eigen::Vector2d i1_t = mesh_->get_uv_data(tri.bt);
		Eigen::Vector2d i2_t = mesh_->get_uv_data(tri.ct);
	
		//Get the Isometric UV Parameterization computerd by rotation 
		Eigen::Vector3d triNorm = tri.getTriangleNormal();
		Eigen::Vector3d xzAxisNorm = Eigen::Vector3d(0.0,1.0,0.0);
		//Eigen::Matrix4f rotMatrix = compute_rotation_matrix_a_to_b(triNorm,xzAxisNorm);

		/*float ua1 = p1.x();
		float va1 = p1.z();
		float ub1 = p2.x();
		float vb1 = p2.z();
		float uc1 = p3.x();
		float vc1 = p3.z();*/

		/*float ua = p1.x();
		float va = p1.z();
		float ub = p2.x();
		float vb = p2.z();
		float uc = p3.x();
		float vc = p3.z();*/


		/*float ua = i0_t.x();
		float va = i0_t.y();
		float ub = i1_t.x();
		float vb = i1_t.y();
		float uc = i2_t.x();
		float vc = i2_t.y();*/

		Eigen::Vector3d x0 = p1-p2;
		double lenAB = x0.norm();

		Eigen::Vector3d xn0 = x0.normalized();
		Eigen::Vector3d x1 = p1-p3;
		double lenAC = x1.norm();


		float ua = 0.0;
		float va = 0.0;
		float ub = lenAB;
		float vb = 0.0;
		float uc = xn0.dot(x1);
		float vc = sqrt(lenAC * lenAC - uc * uc);


		/*if(t%75==0) {
			std::cout << "Triangle#" << t << " " << ua << "," << va << "," << ub << "," << vb << "," << uc << "," << vc << "\n";
			std::cout << "Triangle#" << t << " " << ua1 << "," << va1 << "," << ub1 << "," << vb1 << "," << uc1 << "," << vc1 << "\n";
		}*/

		Vector7d temp_vertex_distribution;

		double det = ua*(vb - vc)+ ub*(vc - va)+ uc*(va - vb);
		//std::cout << "INIT#" << t << " " <<  det << "\n";
		double inv_det = 1.0f/det;
		temp_vertex_distribution(0) = (vb - vc) * inv_det;
		temp_vertex_distribution(1) = (uc - ub) * inv_det;
		temp_vertex_distribution(2) = (vc - va) * inv_det;
		temp_vertex_distribution(3) = (ua - uc) * inv_det;
		temp_vertex_distribution(4) = (va - vb) * inv_det;
		temp_vertex_distribution(5) = (ub - ua) * inv_det;
		temp_vertex_distribution(6) = det;
		vertex_distribution.push_back(temp_vertex_distribution);

		//if(t==34)
			//std::cout << "INIT# " << t << " "<< temp_vertex_distribution(0) << " " << temp_vertex_distribution(1) << " " << temp_vertex_distribution(2) << " " << temp_vertex_distribution(3) << " " << temp_vertex_distribution(4) << " " << temp_vertex_distribution(5) << "\n"; 

	}

	//std::cout << "rendermode value:" << renderMode_ << "\n";
}

Cloth_Data::~Cloth_Data(void)
{
}

void Cloth_Data::finalize_velocity_position()
{
	//Transfer the buffer
	velocity = next_step_velocity;
	mesh_->add_frame_data(next_step_position);

	//Clear the next step buffers
	std::fill(next_step_position.begin(),next_step_position.end(),Eigen::Vector3d::Zero());
	std::fill(next_step_velocity.begin(),next_step_velocity.end(),Eigen::Vector3d::Zero());
}

//@TODO : Implement a nice way to do vector valued heatmaps rendering 
void Cloth_Data::render(int frameId)
{
	if(render_)
	{
		switch(renderMode_)
		{
			case(SHADING) : mesh_->draw_mesh(frameId); break;
			case(WIREFRAME) : mesh_->draw_wireframe(frameId); break;
			case(HMAP_VELOCITY) : // break;
			case(HMAP_ACCLRN) : 
			case(HMAP_SHEAR_FORCE) :
			case(HMAP_BENDING_FORCE) :
			case(HMAP_DAMPING_FORCE) :  break;
		}
	}
}

void Cloth_Data::resetParameters()
{
	std::fill(velocity.begin(),velocity.end(),Eigen::Vector3d::Zero());
	std::fill(next_step_position.begin(),next_step_position.end(),Eigen::Vector3d::Zero());
	std::fill(next_step_velocity.begin(),next_step_velocity.end(),Eigen::Vector3d::Zero());
	mesh_->resetMeshDefaults();
}