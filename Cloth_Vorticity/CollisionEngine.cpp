#include "CollisionEngine.h"


CollisionEngine::CollisionEngine(void)
{
}


CollisionEngine::~CollisionEngine(void)
{
}

void CollisionEngine::resolveClothBodyCollision(Cloth_Data* cloth,Body_Data* body)
{
	int lastFrameId = cloth->getMesh()->get_number_frames();
	//Eigen::Vector3d clothVertexPositions = cloth->getMesh()->
	//int lastFrameId = cloth->getMesh()->get_number_frames();

	//We have 
}
