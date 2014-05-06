#pragma once

#include "Hydro/ExplicitMethod.h"

#include <vector>
#include <functional>

class Hydro;

class ExplicitMethod {
protected:
	void copyState(Hydro *hydro, std::vector<double> Cell::*dst, std::vector<double> Cell::*src);
	void addMulState(Hydro *hydro, std::vector<double> Cell::*dst, std::vector<double> Cell::*src, double dt);
public:
	virtual void operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv) = 0;
};

class ForwardEulerExplicitMethod : public ExplicitMethod {
public:
	virtual void operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv);
};

class RK2ExplicitMethod : public ExplicitMethod {
public:
	virtual void operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv);
};

class RK4ExplicitMethod : public ExplicitMethod {
public:
	virtual void operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv);
};

class ICN3ExplicitMethod : public ExplicitMethod {
public:
	virtual void operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv);
};

