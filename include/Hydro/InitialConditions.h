#pragma once

class IHydro;

class InitialConditions {
public:
	virtual void operator()(IHydro *hydro) = 0;
};

