#pragma once

#include "TensorMath/Vector.h"

#include <thread>

namespace Parallel {
#if 0
#define ParallelFor	std::for_each
#else
template<typename ObjectType, typename Callback, int numThreads = 4>
void For(ObjectType *begin, ObjectType *end, Callback callback) {
	int totalRange = end - begin;
	std::thread threads[numThreads];
	for (int i = 0; i < numThreads; ++i) {
		int beginIndex = i * totalRange / numThreads;
		int endIndex = (i + 1) * totalRange / numThreads;
		threads[i] = std::thread([&,beginIndex,endIndex](){
			std::for_each(begin + beginIndex, begin + endIndex, callback);
		});
	}
	for (int i = 0; i < numThreads; ++i) {
		threads[i].join();
	}
};
#endif
}

