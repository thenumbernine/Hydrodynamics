#pragma once

struct IHydro;

struct BoundaryMethod {
	virtual void operator()(IHydro *ihydro) = 0;
};

