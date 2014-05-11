#pragma once

#include "TensorMath/Grid.h"	//RangeObj
#include "Hydro/BoundaryMethod.h"

template<typename Hydro>
class MirrorBoundaryMethod : public BoundaryMethod {
public:
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };

	typedef typename Hydro::IVector IVector;
	
	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void MirrorBoundaryMethod<Hydro>::operator()(IHydro *ihydro) {
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

			//lhs
			index(side) = 0;
			hydro->cells(index).state(0) = hydro->cells(index + di*3).state(0);
			hydro->cells(index + di).state(0) = hydro->cells(index + di*2).state(0);
			for (int k = 0; k < rank; ++k) {
				hydro->cells(index).state(k+1) = -hydro->cells(index + di*3).state(k+1);
				hydro->cells(index + di).state(k+1) = -hydro->cells(index + di*2).state(k+1);
			}
			hydro->cells(index).state(rank+1) = hydro->cells(index + di*3).state(rank+1);
			hydro->cells(index + di).state(rank+1) = hydro->cells(index + di*2).state(rank+1);

			//rhs
			index(side) = hydro->size(side)-1;
			hydro->cells(index).state(0) = hydro->cells(index - di*3).state(0);
			hydro->cells(index - di).state(0) = hydro->cells(index - di*2).state(0);
			for (int k = 0; k < rank; ++k) {
				hydro->cells(index).state(k+1) = -hydro->cells(index - di*3).state(k+1);
				hydro->cells(index - di).state(k+1) = -hydro->cells(index - di*2).state(k+1);
			}
			hydro->cells(index).state(rank+1) = hydro->cells(index - di*3).state(rank+1);
			hydro->cells(index - di).state(rank+1) = hydro->cells(index - di*2).state(rank+1);
		};

//if rank is one then our range object is zero, even though we want to run the lambda once ...
		if (rank == 1) {
			callback(BoundaryIntVector());
		} else {
			std::for_each(boundary.begin(), boundary.end(), callback);
		}
	}
	
#if 0	
	hydro->cells(0).state(0) = hydro->cells(3).state(0);
	hydro->cells(1).state(0) = hydro->cells(2).state(0);
	hydro->cells(hydro->size(0)-2).state(0) = hydro->cells(hydro->size(0)-3).state(0);
	hydro->cells(hydro->size(0)-1).state(0) = hydro->cells(hydro->size(0)-4).state(0);
	hydro->cells(0).state(1) = -hydro->cells(3).state(1);
	hydro->cells(1).state(1) = -hydro->cells(2).state(1);
	hydro->cells(hydro->size(0)-2).state(1) = -hydro->cells(hydro->size(0)-3).state(1);
	hydro->cells(hydro->size(0)-1).state(1) = -hydro->cells(hydro->size(0)-4).state(1);
	hydro->cells(0).state(2) = hydro->cells(3).state(2);
	hydro->cells(1).state(2) = hydro->cells(2).state(2);
	hydro->cells(hydro->size(0)-2).state(2) = hydro->cells(hydro->size(0)-3).state(2);
	hydro->cells(hydro->size(0)-1).state(2) = hydro->cells(hydro->size(0)-4).state(2);
#endif
}
