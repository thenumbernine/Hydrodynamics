#include "Hydro/BoundaryMethod.h"
#include "Hydro/Hydro.h"

//putting them here for now

void PeriodicBoundaryMethod::operator()(Hydro *hydro) {
	hydro->cells[0].state[0] = hydro->cells[hydro->size-4].state[0];
	hydro->cells[1].state[0] = hydro->cells[hydro->size-3].state[0];
	hydro->cells[hydro->size-2].state[0] = hydro->cells[2].state[0];
	hydro->cells[hydro->size-1].state[0] = hydro->cells[3].state[0];
	hydro->cells[0].state[1] = hydro->cells[hydro->size-4].state[1];
	hydro->cells[1].state[1] = hydro->cells[hydro->size-3].state[1];
	hydro->cells[hydro->size-2].state[1] = hydro->cells[2].state[1];
	hydro->cells[hydro->size-1].state[1] = hydro->cells[3].state[1];
	hydro->cells[0].state[2] = hydro->cells[hydro->size-4].state[2];
	hydro->cells[1].state[2] = hydro->cells[hydro->size-3].state[2];
	hydro->cells[hydro->size-2].state[2] = hydro->cells[2].state[2];
	hydro->cells[hydro->size-1].state[2] = hydro->cells[3].state[2];
}

void MirrorBoundaryMethod::operator()(Hydro *hydro) {
	hydro->cells[0].state[0] = hydro->cells[3].state[0];
	hydro->cells[1].state[0] = hydro->cells[2].state[0];
	hydro->cells[hydro->size-2].state[0] = hydro->cells[hydro->size-3].state[0];
	hydro->cells[hydro->size-1].state[0] = hydro->cells[hydro->size-4].state[0];
	hydro->cells[0].state[1] = -hydro->cells[3].state[1];
	hydro->cells[1].state[1] = -hydro->cells[2].state[1];
	hydro->cells[hydro->size-2].state[1] = -hydro->cells[hydro->size-3].state[1];
	hydro->cells[hydro->size-1].state[1] = -hydro->cells[hydro->size-4].state[1];
	hydro->cells[0].state[2] = hydro->cells[3].state[2];
	hydro->cells[1].state[2] = hydro->cells[2].state[2];
	hydro->cells[hydro->size-2].state[2] = hydro->cells[hydro->size-3].state[2];
	hydro->cells[hydro->size-1].state[2] = hydro->cells[hydro->size-4].state[2];
}

void ConstantBoundaryMethod::operator()(Hydro *hydro) {
	hydro->cells[0].state[0] = 0;
	hydro->cells[1].state[0] = 0;
	hydro->cells[hydro->size-2].state[0] = 0;
	hydro->cells[hydro->size-1].state[0] = 0;
	hydro->cells[0].state[1] = 0;
	hydro->cells[1].state[1] = 0;
	hydro->cells[hydro->size-2].state[1] = 0;
	hydro->cells[hydro->size-1].state[1] = 0;
	hydro->cells[0].state[2] = 0;
	hydro->cells[1].state[2] = 0;
	hydro->cells[hydro->size-2].state[2] = 0;
	hydro->cells[hydro->size-1].state[2] = 0;
}

void FreeFlowBoundaryMethod::operator()(Hydro *hydro) {
	hydro->cells[0].state[0] = hydro->cells[1].state[0] = hydro->cells[2].state[0];
	hydro->cells[hydro->size-1].state[0] = hydro->cells[hydro->size-2].state[0] = hydro->cells[hydro->size-3].state[0];
	hydro->cells[0].state[1] = hydro->cells[1].state[1] = hydro->cells[2].state[1];
	hydro->cells[hydro->size-1].state[1] = hydro->cells[hydro->size-2].state[1] = hydro->cells[hydro->size-3].state[1];
	hydro->cells[0].state[2] = hydro->cells[1].state[2] = hydro->cells[2].state[2];
	hydro->cells[hydro->size-1].state[2] = hydro->cells[hydro->size-2].state[2] = hydro->cells[hydro->size-3].state[2];
}

