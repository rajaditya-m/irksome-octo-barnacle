#pragma once

#include "cloth_data.h"
//#include "ImplicitFEMSolver.h"
#include "ImplicitMassSpringSolver.h"
#include "body_data.h"


class SimulationEngine
{
public:
	SimulationEngine(Cloth_Data* cloth_data,ImplicitMassSpringSolver* solver,Body_Data* body_data);
	~SimulationEngine(void);

	void generate_next_frame();
	void populatePerVertexBuffer();

public:
	

	Cloth_Data* clothData_;
	ImplicitMassSpringSolver* solver_;
	Body_Data* bodyData_;
};

