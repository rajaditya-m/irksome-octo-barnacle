#include <iostream>
#include <Eigen\Dense>

inline void DBG__Vector_Dump(Eigen::Vector3d vec)
{
	std::cout << "[" << vec[0] << "," << vec[1] << "," << vec[2] << "]\n";
}

inline void DBG__Vector_Dump(Eigen::Vector3d vec,int seq_num,const char* header="")
{
	std::cout << header << " Seq#" << seq_num << " [" << vec[0] << "," << vec[1] << "," << vec[2] << "]\n";
}
