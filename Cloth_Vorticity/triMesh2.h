#pragma once

#include <vector>
#include <set>
#include <fstream>
#include <iostream>

#include "Eigen\Dense"

#include <QGLWidget>


#include "triangles.h"
#include "edge.h"
#include "global_typedefs.h"
#include "struct_decls.h"
#include "mesh_material.h"
#include "grid.h"


class TriMesh
{
public:
	TriMesh(const char* meshName,
				 const char* file_location_vertex,
				 const char* file_location_mesh,
				 const char* file_location_texture,
				 const char* material_xml,
				 float pdfs);

	TriMesh(const char* meshName,
				 const char* file_location_vertex,
				 const char* file_location_mesh,
				 const char* material_xml,
				 float pdfs);

	TriMesh(const char* meshName,
				 const char* file_location_gen,
				 const char* material_xml,
				 int num_frames,
				 float pdfs);

	~TriMesh(void);

	//More important functions 
	int get_number_frames()	const										{ return point_data.size();				}
	int get_number_vertices() const 									{ return point_data[0].size() ;			}
	int get_number_triangles() const									{ return mesh_data.size();				}
	int get_number_edges() const										{ return edge_list.size();				}
	float get_average_edge_length()	const								{ return average_edge_length_;			}
	float get_longest_edge_length()	const								{ return longest_edge_length_;			}
	Pair_Vec3f& get_super_bounding_box()								{ return super_bounding_box;			}

	//Accessor functions for this class 
	Triangles get_triangle(int tri_idx)	const							{ return mesh_data[tri_idx];			}
	bool is_part_of_mesh(int vert_idx)	const							{ return part_of_mesh[vert_idx];		}
	Eigen::Vector3d get_point_data(int idx, int frame = 0)	const		{ return point_data[frame][idx];		}
	Eigen::Vector3d get_normal_data(int idx, int frame = 0)	const		{ return normal_data[frame][idx];		}
	Eigen::Vector2f get_uv_data(int idx) const							{ return texture_data[idx];				}
	Eigen::Vector2f get_texture_data(int idx) const						{ return texture_data[idx];				}
	std::set<Edge>& get_edge_list()										{ return edge_list;						}
	Grid* get_grid_obj() const											{ return voxel_grid;					}

	//Mutator functions for this class 
	void set_point_data(const Eigen::Vector3d p,int idx,int frame=0)	{ point_data[frame][idx] = p;			}
	void add_point_data_frame(std::vector<Eigen::Vector3d> &vec)		{ point_data.push_back(vec);			}

	//Other functions
	void generate_frame_normals(int frame_idx) ;
	void add_frame_data(std::vector<Eigen::Vector3d> &vec);

	//Other functions 
	void setStatic(bool res)     { isStatic_ = res;}
	void resetMeshDefaults();

	//Drawing Functions 
	void draw_mesh(int idx);
	void draw_wireframe(int idx);
	void draw_bounding_box(int idx);
	void draw_super_bounding_box();
	
	
	void update_voxelization();					//This takes every frame data that is present and updates the voxelization
	void voxelize_frame(int frame_id);			//This takes one frame id and voxelizes it 

private:
	//Data containign frame by frame stuff 
	Vector2d_Vec3f point_data;
	Vector2d_Vec3f normal_data;
	Vector_PairVec3f bounding_box;

	//Mesh specific Data
	std::vector<Triangles> mesh_data;
	std::vector<bool> part_of_mesh;
	std::set<Edge> edge_list;
	Vector_Vec2f texture_data;

	const float CONST_SCALING;

	//Useful Statistics will be updated as and when needed
	float average_edge_length_;
	float longest_edge_length_;

	//Super bounding box measurements
	Pair_Vec3f super_bounding_box;
	
	//Other data
	bool has_uv_coods_;

	//This should have some material properties 
	Mesh_Material* material;
	//@TODO : Support for texture rendering 
	//@TODO : Material file reading support 

	//Voxelization wrt some grid boundaries
	bool voxelized_ ;
	Grid* voxel_grid;
	//Voxel to primitive DS
	std::vector< std::vector< std::vector <int> > > voxel_to_triangle;
	std::vector< std::vector< std::vector <int> > > voxel_to_point;
	std::vector< std::vector< std::vector <int> > > voxel_to_edge;
	//Primitive to voxel DS
	std::vector< std::vector< std::vector <Voxel> > > triangle_to_voxel;
	std::vector< std::vector< std::vector <Voxel> > > edge_to_voxel;
	std::vector< std::vector< Voxel > > point_to_voxel;

	//static mesh
	bool isStatic_;

};