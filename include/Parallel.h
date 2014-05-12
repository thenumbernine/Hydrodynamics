#pragma once

#include "TensorMath/Vector.h"

#include <thread>

namespace Parallel {

#if 0	//single-threaded
#define ParallelFor	std::for_each
#else	//multi-threaded
template<typename Iterator, typename Callback, int numThreads = 4>
void For(Iterator begin, Iterator end, Callback callback) {
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

template<
	typename Iterator, 
	typename Result, 
	int numThreads = 4,
	typename Callback = std::function<Result (typename std::iterator_traits<Iterator>::value_type &)>,
	typename Combiner = std::function<Result (Result, Result)>
>
Result Reduce(
	Iterator begin, 
	Iterator end, 
	Callback callback,
	Result initialValue = Result(),
	Combiner combiner = [&](Result a, Result b) -> Result { return a + b; })
{
	int totalRange = end - begin;
	std::thread threads[numThreads];
	Result results[numThreads];
	for (int i = 0; i < numThreads; ++i) {
		int beginIndex = i * totalRange / numThreads;
		int endIndex = (i + 1) * totalRange / numThreads;
		threads[i] = std::thread([&,beginIndex,endIndex,i]() {
			results[i] = initialValue;
			std::for_each(begin + beginIndex, begin + endIndex, [&](typename std::iterator_traits<Iterator >::value_type &value) {
				results[i] = combiner(results[i], callback(value));
			});
		});
	}
	for (int i = 0; i < numThreads; ++i) {
		threads[i].join();
		initialValue = combiner(initialValue, results[i]);
	}
	
	return initialValue;
}

}

