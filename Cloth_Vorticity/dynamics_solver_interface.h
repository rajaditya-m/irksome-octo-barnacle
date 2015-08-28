#pragma once

//#include "cloth_data.h"
class Cloth_Data;

class Dynamics_Solver_Interface
{
public:
	virtual void advance_time_step(Cloth_Data* cloth) = 0;
};