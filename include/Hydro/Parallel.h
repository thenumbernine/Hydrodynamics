#pragma once

#include <memory>

#ifdef DISABLE_MULTITHREAD

#include <functional>
#include <algorithm>

namespace Hydrodynamics {

//TODO this behavior in the Parallel library
// but putting it as a condition for the Parallel class would mean a runtime check

struct Parallel {
	Parallel(int) {}

	static size_t getNumThreads() { return 1; }

	template<typename Iterator, typename Callback>
	void foreach(Iterator begin, Iterator end, Callback callback) {
		std::for_each(begin, end, callback);
	}

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
		Combiner combiner = std::plus<Result>())
	{
		for (Iterator i = begin; i != end; ++i) {
			initialValue = combiner(initialValue, callback(*i));
		}
		return initialValue;
	}
};

}

#else	//DISABLE_MULTITHREAD

#include "Parallel/Parallel.h"

namespace Hydrodynamics {
using Parallel = ::Parallel::Parallel;
}

#endif	//DISABLE_MULTITHREAD

namespace Hydrodynamics {
extern std::shared_ptr<Parallel> parallel;
}
