#pragma once

#include "TensorMath/Vector.h"

#include <thread>

template<int rank, typename Callback>
void RangeParallelFor(Vector<int, rank> min, Vector<int, rank> max, Callback callback) {
	typedef Vector<int, rank> IVector;
	IVector mid = (min + max) / 2;
	std::vector<RangeObj<rank>> ranges;
	RangeObj<rank> splits(IVector(0), IVector(2));
	std::for_each(splits.begin(), splits.end(), [&](IVector split) {
		RangeObj<rank> r;
		for (int k = 0; k < rank; ++k) {
			r.min(k) = split(k) == 0 ? min(k) : mid(k);
			r.max(k) = split(k) == 0 ? mid(k) : max(k);
		}
		ranges.push_back(r);
	});
	
	std::vector<std::thread> threads;
	for (int i = 0; i < (int)ranges.size(); ++i) {
		RangeObj<rank> &r = ranges[i];
		threads.push_back(std::thread([&](){
			std::for_each(r.begin(), r.end(), [&](IVector index) {
				callback(index);
			});
		}));
	}

	std::for_each(threads.begin(), threads.end(), [](std::thread &t){
		t.join();
	});
}

