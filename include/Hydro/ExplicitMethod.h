#pragma once

#include <functional>

template<typename Hydro>
class ExplicitMethod {
protected:
	enum { numberOfStates = Hydro::numberOfStates };
	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Cell::StateVector StateVector;
	
	void copyState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src);
	void addMulState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src, Real dt);
public:
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv) = 0;
};

template<typename Hydro>
void ExplicitMethod<Hydro>::copyState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src) {
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
		for (int state = 0; state < numberOfStates; ++state) {
			(cell.*dst)(state) = (cell.*src)(state);
		}
	});
}

template<typename Hydro>
void ExplicitMethod<Hydro>::addMulState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src, Real dt) {
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
		for (int state = 0; state < numberOfStates; ++state) {
			(cell.*dst)(state) += (cell.*src)(state) * dt;
		}
	});
}

