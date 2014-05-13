#pragma once

struct IHydro;

struct InitialConditions {
	virtual void operator()(IHydro *hydro, double noise) = 0;
};

