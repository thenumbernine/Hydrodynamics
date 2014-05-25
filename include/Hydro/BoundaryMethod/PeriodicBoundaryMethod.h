#pragma once

#include "TensorMath/Grid.h"	//RangeObj
#include "Hydro/BoundaryMethod.h"

template<typename Hydro>
class PeriodicBoundaryMethod : public BoundaryMethod {
public:
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };

	typedef typename Hydro::IVector IVector;

	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void PeriodicBoundaryMethod<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	for (int side = 0; side < rank; ++side) {
		//volume range over all other dimensions ...
		typedef ::Vector<int, rank-1> BoundaryIntVector;
		BoundaryIntVector min, max;
		for (int k = 0; k < rank - 1; ++k) {
			if (k < side) {
				max(k) = hydro->size(k);
			} else {
				max(k) = hydro->size(k+1);
			}
		}
		RangeObj<rank-1> boundary(min, max);
		
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
			std::for_each(boundary.begin(), boundary.end(), callback);
		}
	}
}

