#pragma once

#include "Tensor/Grid.h"	//RangeObj
#include "Hydro/Boundary/Boundary.h"

namespace Hydrodynamics {
namespace Boundary {

template<typename Hydro>
class Periodic : public Boundary {
public:
	static constexpr auto rank = Hydro::rank;
	static constexpr auto numberOfStates = Hydro::numberOfStates;

	using IVector = typename Hydro::IVector;

	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void Periodic<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	for (int side = 0; side < rank; ++side) {
		//volume range over all other dimensions ...
		using BoundaryIntVector = Tensor::intN<MinOneMinusOne<rank>::value>; 
		using BoundaryRangeObj = Tensor::RangeObj<MinOneMinusOne<rank>::value>;
		BoundaryIntVector min, max;
		for (int k = 0; k < rank - 1; ++k) {
			if (k < side) {
				max(k) = hydro->size(k);
			} else {
				max(k) = hydro->size(k+1);
			}
		}
		BoundaryRangeObj boundary(min, max);
		
		auto callback = [&](BoundaryIntVector boundaryIndex) {
			IVector index;
			for (int k = 0; k < rank-1; ++k) {
				if (k < side) {
					index(k) = boundaryIndex(k);
				} else {
					index(k+1) = boundaryIndex(k);
				}
			}
		
			IVector di;
			di(side) = 1;

			IVector leftIndex = index;
			leftIndex(side) = 0;
			
			IVector rightIndex = index;
			rightIndex(side) = hydro->size(side)-1;
			
			//lhs
			hydro->cells(leftIndex).second.state = hydro->cells(rightIndex - di*3).second.state;
			hydro->cells(leftIndex + di).second.state = hydro->cells(rightIndex - di*2).second.state;

			//rhs
			hydro->cells(rightIndex).second.state = hydro->cells(leftIndex + di*3).second.state;
			hydro->cells(rightIndex - di).second.state = hydro->cells(leftIndex + di*2).second.state;
		};

//if rank is one then our range object is zero, even though we want to run the lambda once ...
		if (rank == 1) {
			callback(BoundaryIntVector());
		} else {
#if PLATFORM_MSVC	//Microsoft is why we can't have nice things.
			for (typename BoundaryRangeObj::iterator i = boundary.begin(); i != boundary.end(); ++i) {
				callback(*i);
			}
#else
			std::for_each(boundary.begin(), boundary.end(), callback);
#endif
		}
	}
}

}
}
