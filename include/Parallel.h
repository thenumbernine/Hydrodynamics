#pragma once

#include "TensorMath/Vector.h"

#include <thread>

template<int rank, typename Callback>
void RangeParallelFor(Vector<int, rank> min, Vector<int, rank> max, Callback callback) {
	typedef Vector<int, rank> IVector;
	IVector mid = (min + max) / 2;
	
	RangeObj<rank> ranges[1 << rank];

	int i = 0;
	RangeObj<rank> splits(IVector(0), IVector(2));
	std::for_each(splits.begin(), splits.end(), [&](IVector split) {
		RangeObj<rank> &r = ranges[i++];
		for (int k = 0; k < rank; ++k) {
			r.min(k) = split(k) == 0 ? min(k) : mid(k);
			r.max(k) = split(k) == 0 ? mid(k) : max(k);
		}
	});
	
	std::thread threads[1 << rank];
	for (int i = 0; i < 1 << rank; ++i) {
		RangeObj<rank> &r = ranges[i];
		threads[i] = std::thread([&](){
			std::for_each(r.begin(), r.end(), [&](IVector index) {
				callback(index);
			});
		});
	}

	std::for_each(begin(threads), end(threads), [](std::thread &t){
		t.join();
	});
}

