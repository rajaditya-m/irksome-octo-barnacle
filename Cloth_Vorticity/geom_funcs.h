#pragma once

#include <Eigen\Dense>

inline float triangle_area(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c)
{
	Eigen::Vector3d vBA = b-a;
	Eigen::Vector3d vCA = c-a;
	Eigen::Vector3d crossPdk = vBA.cross(vCA);
	float area = crossPdk.norm()*0.5f;
	return area;
}

inline Eigen::Vector2f compute_outward_normal(Eigen::Vector2f a, Eigen::Vector2f b, Eigen::Vector2f c)
{
	Eigen::Vector2f ab = b-a;
	Eigen::Vector2f ac = c-a;
	Eigen::Vector2f N(ab.y(),-ab.x());

	float d = N.dot(ac);

	if(d>0.0)
		N = -N;

	N.normalized();
	
	return N;
}

inline Eigen::Vector3d compute_triangle_normal(Eigen::Vector3d a,Eigen::Vector3d b,Eigen::Vector3d c)
{
	Eigen::Vector3d ab = b-a;
	Eigen::Vector3d ac = c-a;
	Eigen::Vector3d normal = ab.cross(ac);
	normal = normal.normalized();
	return normal;
}


inline Eigen::Matrix3d compute_rotation_matrix_a_to_b(Eigen::Vector3d norm_in_a,Eigen::Vector3d norm_in_b)
{
	float cos_theta = norm_in_a.dot(norm_in_b);
	
	Eigen::Vector3d unit_axis = norm_in_a.cross(norm_in_b);
	unit_axis.normalized();

	float x = unit_axis.x();
	float y = unit_axis.y();
	float z = unit_axis.z();

	float c = cos_theta;
	float term_under_sqrt = 1.0f-c*c;
	if(term_under_sqrt<0.0)
		std::cout << "[ERROR] Error while computing the rotation matrix\n";
	float s = sqrt(term_under_sqrt);
	float C = 1-c;

	Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Zero();

	rot_mat(0,0) = x*x*C+c;
	rot_mat(0,1) = x*y*C-z*s;
	rot_mat(0,2) = x*z*C+y*s;
	rot_mat(1,0) = y*x*C+z*s;
	rot_mat(1,1) = y*y*C+c;
	rot_mat(1,2) = y*z*C-x*s;
	rot_mat(2,0) = z*x*C-y*s;
	rot_mat(2,1) = z*y*C+x*s;
	rot_mat(2,2) = z*z*C+c;

	return rot_mat;
}

inline bool compare_point_leq(Eigen::Vector3d a, Eigen::Vector3d b)
{
	for(int i=0;i<3;i++)
	{
		if(a(i)>b(i))
			return false;
	}
	return true;
}

inline bool compare_point_geq(Eigen::Vector3d a, Eigen::Vector3d b)
{
	for(int i=0;i<3;i++)
	{
		if(a(i)<b(i))
			return false;
	}
	return true;
}

/*
inline Eigen::Matrix3d compute_rotation_matrix_a_to_b(Eigen::Vector3d norm_in_a,Eigen::Vector3d norm_in_b)
{
	norm_in_a.normalize();
	norm_in_b.normalize();

	Eigen::Vector3d unit_axis = norm_in_a.cross(norm_in_b);
	float len = unit_axis.norm();
	if(fabs(len) < 1.0e-14)
	{
		Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
		return identity;
	}
	unit_axis.normalized();

	Vector3d v2 = cross(v0,v1);
			if (v2.length()<1e-20)	// TODO
			{
				//std::cout<<"("<<v2.x<<","<<v2.y<<","<<v2.z<<") ";
				Matrix3d id;
				id.identity();
				return id;				
			}
			v2.normalize();



			Matrix3d U(v0,v2,cross(v0,v2));
			Matrix3d V(v1,v2,cross(v1,v2));
			V.transpose();
			return V*U;
}
*/