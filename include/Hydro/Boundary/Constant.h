#pragma once

#include "Hydro/Boundary/Boundary.h"

namespace Boundary {

template<typename Hydro>
class Constant : public Boundary {
public:
	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void Constant<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->cells(0).second.state(0) = 0;
	hydro->cells(1).second.state(0) = 0;
	hydro->cells(hydro->size(0)-2).second.state(0) = 0;
	hydro->cells(hydro->size(0)-1).second.state(0) = 0;
	hydro->cells(0).second.state(1) = 0;
	hydro->cells(1).second.state(1) = 0;
	hydro->cells(hydro->size(0)-2).second.state(1) = 0;
	hydro->cells(hydro->size(0)-1).second.state(1) = 0;
	hydro->cells(0).second.state(2) = 0;
	hydro->cells(1).second.state(2) = 0;
	hydro->cells(hydro->size(0)-2).second.state(2) = 0;
	hydro->cells(hydro->size(0)-1).second.state(2) = 0;
}

};

