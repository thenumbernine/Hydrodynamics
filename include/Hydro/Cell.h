#pragma once

#include <vector>

//cell-centered values
struct Cell {
	//static values
	double x;	//position

	//dynamic values
	std::vector<double> state;	//state variables

	//primitives ... only used in getPrimitives ... by draw ... hmm ...
	std::vector<double> primitives;

	//used for Burgers
	//aux variables:
	double pressure;
	
	//used for Godunov
	std::vector<double> stateLeft, stateRight;

	//used for integration
	std::vector<double> dq_dt;
	std::vector<double> tmpState0;
	std::vector<double> tmpState1;
	std::vector<double> tmpState2;
	std::vector<double> tmpState3;
	std::vector<double> tmpState4;
};


