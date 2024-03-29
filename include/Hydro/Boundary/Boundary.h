#pragma once

#include "Hydro/IHydro.h"

namespace Hydrodynamics {
namespace Boundary {

struct Boundary {
	virtual ~Boundary() {}
	virtual void operator()(IHydro *ihydro) = 0;
};

}
}
