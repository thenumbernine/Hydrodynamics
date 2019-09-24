#pragma once

#include "Tensor/Grid.h"	//RangeObj
#include "Hydro/Boundary/Boundary.h"

namespace Hydrodynamics {
namespace Boundary {

template<typename Hydro>
class Mirror : public Boundary {
public:
	static constexpr auto rank = Hydro::rank;
	static constexpr auto numberOfStates = Hydro::numberOfStates;

	using IVector = typename Hydro::IVector;
	
	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void Mirror<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	for (int side = 0; side < rank; ++side) {
		//volume range over all other dimensions ...
		using BoundaryIntVector = Tensor::Vector<int, MinOneMinusOne<rank>::value>; 
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

			//lhs
			index(side) = 0;
			hydro->cells(index).second.state(0) = hydro->cells(index + di*3).second.state(0);
			hydro->cells(index + di).second.state(0) = hydro->cells(index + di*2).second.state(0);
			for (int k = 0; k < rank; ++k) {
				hydro->cells(index).second.state(k+1) = -hydro->cells(index + di*3).second.state(k+1);
				hydro->cells(index + di).second.state(k+1) = -hydro->cells(index + di*2).second.state(k+1);
			}
			hydro->cells(index).second.state(rank+1) = hydro->cells(index + di*3).second.state(rank+1);
			hydro->cells(index + di).second.state(rank+1) = hydro->cells(index + di*2).second.state(rank+1);

			//rhs
			index(side) = hydro->size(side)-1;
			hydro->cells(index).second.state(0) = hydro->cells(index - di*3).second.state(0);
			hydro->cells(index - di).second.state(0) = hydro->cells(index - di*2).second.state(0);
			for (int k = 0; k < rank; ++k) {
				hydro->cells(index).second.state(k+1) = -hydro->cells(index - di*3).second.state(k+1);
				hydro->cells(index - di).second.state(k+1) = -hydro->cells(index - di*2).second.state(k+1);
			}
			hydro->cells(index).second.state(rank+1) = hydro->cells(index - di*3).second.state(rank+1);
			hydro->cells(index - di).second.state(rank+1) = hydro->cells(index - di*2).second.state(rank+1);
		};

//if rank is one then our range object is zero, even though we want to run the lambda once ...
		if (rank == 1) {
			callback(BoundaryIntVector());
		} else {
#ifdef PLATFORM_msvc	//Microsoft is why we can't have nice things.
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
