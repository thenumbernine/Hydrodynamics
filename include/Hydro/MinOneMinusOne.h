#pragma once
//this is used often enough, and MSVC can't understand a 0-sized array ...

namespace Hydrodynamics {

template<int n>
struct MinOneMinusOne {
	enum { value = n - 1 };
};

template<>
struct MinOneMinusOne<1> {
	enum { value = 1 };
};

}
