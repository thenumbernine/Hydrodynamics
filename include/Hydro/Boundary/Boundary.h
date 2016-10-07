#pragma once

#include "Hydro/IHydro.h"

namespace Boundary {

struct Boundary {
	virtual void operator()(IHydro *ihydro) = 0;
};

};
