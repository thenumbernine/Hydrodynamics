#pragma once

class IHydro;

class BoundaryMethod {
public:
	virtual void operator()(IHydro *ihydro) = 0;
};

