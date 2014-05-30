#pragma once

#include <thread>
#include <functional>
#include <algorithm>

template<int numThreads = 4>
struct ParallelCount {
	template<typename Iterator, typename Callback>
	void foreach(Iterator begin, Iterator end, Callback callback) {
		int totalRange = end - begin;
		std::thread threads[numThreads];
		for (int i = 0; i < numThreads; ++i) {
			int beginIndex = i * totalRange / numThreads;
			int endIndex = (i + 1) * totalRange / numThreads;
			threads[i] = std::thread([&,beginIndex,endIndex]() {
				std::for_each(begin + beginIndex, begin + endIndex, callback);
			});
		}
		for (int i = 0; i < numThreads; ++i) {
			threads[i].join();
		}
	};

	// a shy step away from std::accumulate
	// in that the values in the iterator are mapped first (via callback)
	// before they are accumulated.
	// I could make a new structure for buffering my intermediate values, but I don't really want to.
	template<
		typename Iterator, 
		typename Result, 
		typename Callback = std::function<Result (typename std::iterator_traits<Iterator>::value_type &)>,
		typename Combiner = std::function<Result (Result, Result)>
	>
	Result reduce(
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
};

typedef ParallelCount<4> Parallel;

#include "Common/Singleton.h"
extern Singleton<Parallel> parallel;

