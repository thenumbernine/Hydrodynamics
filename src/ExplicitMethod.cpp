#include "Hydro/Hydro.h"
#include "Hydro/ExplicitMethod.h"
#include "Hydro/EquationOfState.h"

void ExplicitMethod::copyState(Hydro *hydro, std::vector<double> Cell::*dst, std::vector<double> Cell::*src) {
	for (int i = 0; i < hydro->size; ++i) {
		for (int j = 0; j < hydro->equationOfState->numberOfStates(); ++j) {
			(hydro->cells[i].*dst)[j] = (hydro->cells[i].*src)[j];
		}
	}
}

void ExplicitMethod::addMulState(Hydro *hydro, std::vector<double> Cell::*dst, std::vector<double> Cell::*src, double dt) {
	for (int i = 0; i < hydro->size; ++i) {
		for (int j = 0; j < hydro->equationOfState->numberOfStates(); ++j) {
			(hydro->cells[i].*dst)[j] += (hydro->cells[i].*src)[j] * dt;
		}
	}
}

void ForwardEulerExplicitMethod::operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv) {
	std::vector<double> Cell::*dq_dt = &Cell::tmpState0;
	copyState(hydro, dq_dt, &Cell::state);
	deriv(hydro, dt, dq_dt);
	addMulState(hydro, &Cell::state, dq_dt, dt);
}

void RK2ExplicitMethod::operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv) {
	std::vector<double> Cell::*src = &Cell::tmpState0;
	copyState(hydro, src, &Cell::state);
	std::vector<double> Cell::*k1 = &Cell::tmpState1;
	copyState(hydro, k1, &Cell::state);
	deriv(hydro, dt, k1);
	addMulState(hydro, &Cell::state, k1, .5 * dt);
	std::vector<double> Cell::*k2 = &Cell::tmpState2;
	deriv(hydro, dt, k2);
	copyState(hydro, &Cell::state, src);
	addMulState(hydro, &Cell::state, k2, dt);
}

void RK4ExplicitMethod::operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv) {
	std::vector<double> Cell::*src = &Cell::tmpState0;
	copyState(hydro, src, &Cell::state);
	std::vector<double> Cell::*k1 = &Cell::tmpState1;
	copyState(hydro, k1, &Cell::state);
	deriv(hydro, dt, k1);
	addMulState(hydro, &Cell::state, k1, .5 * dt);
	std::vector<double> Cell::*k2 = &Cell::tmpState2;
	copyState(hydro, k2, &Cell::state);
	deriv(hydro, dt, k2);
	copyState(hydro, &Cell::state, src);
	addMulState(hydro, &Cell::state, k2, .5 * dt);
	std::vector<double> Cell::*k3 = &Cell::tmpState3;
	copyState(hydro, k3, &Cell::state);
	deriv(hydro, dt, k3);
	copyState(hydro, &Cell::state, src);
	addMulState(hydro, &Cell::state, k3, dt);
	std::vector<double> Cell::*k4 = &Cell::tmpState4;
	copyState(hydro, k4, &Cell::state);
	deriv(hydro, dt, k4);
	copyState(hydro, &Cell::state, src);
	addMulState(hydro, &Cell::state, k1, dt / 6.);
	addMulState(hydro, &Cell::state, k2, dt / 3.);
	addMulState(hydro, &Cell::state, k3, dt / 3.);
	addMulState(hydro, &Cell::state, k4, dt / 6.);
}
void ICN3ExplicitMethod::operator()(Hydro *hydro, double dt, std::function<void(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt)> deriv) {
	//first iteration
	std::vector<double> Cell::*srcQ = &Cell::tmpState0;
	copyState(hydro, srcQ, &Cell::state);
	std::vector<double> Cell::*firstK = &Cell::tmpState1;
	copyState(hydro, firstK, &Cell::state);
	deriv(hydro, dt, firstK);
	addMulState(hydro, &Cell::state, firstK, dt);
	std::vector<double> Cell::*k = &Cell::tmpState2;
	copyState(hydro, k, &Cell::state);

	//second and so on
	for (int i = 1; i < 3; ++i) {
		deriv(hydro, dt, k);
		copyState(hydro, &Cell::state, srcQ);
		addMulState(hydro, &Cell::state, k, .5 * dt);
		addMulState(hydro, &Cell::state, firstK, .5 * dt);
	}
}
