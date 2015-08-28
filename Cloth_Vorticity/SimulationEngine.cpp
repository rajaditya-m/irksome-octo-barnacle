#include "SimulationEngine.h"


SimulationEngine::SimulationEngine(Cloth_Data* cloth_data,FEM_Solver* solver,Body_Data* body_data)
{
	clothData_ = cloth_data;
	solver_ = solver;
	solver_->initialize(clothData_);
	bodyData_ = body_data;
}


SimulationEngine::~SimulationEngine(void)
{
}

void SimulationEngine::generate_next_frame()
{
	//Advance the time step to generate the next set of collision as per internal stuff
	solver_->advance_time_step(clothData_);
	
	//Then apply the collision engine
	//collision_engine->resolve_cloth_body_collisions(cloth_data,body_data);
	
	//Then resolve self collisions 
	//collision_engine->resolve_cloth_cloth_collisions(cloth_data);
	
	//Finalize the velocity and positions
	clothData_->finalize_velocity_position();

	//Control the flow of parameter to the cloth engine depending upon the type of rendering 
	populatePerVertexBuffer();
}

void SimulationEngine::populatePerVertexBuffer()
{
	switch(clothData_->getRenderMode())
	{
		case(HMAP_ACCLRN) : clothData_->setPerVertexVectorBuffer(solver_->getAcceleration()); break;
		case(HMAP_BENDING_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getBendingForce()); break;
		case(HMAP_SHEAR_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getShearForce()); break;
		case(HMAP_DAMPING_FORCE) : clothData_->setPerVertexVectorBuffer(solver_->getDampingForce()); break;
	}
}
