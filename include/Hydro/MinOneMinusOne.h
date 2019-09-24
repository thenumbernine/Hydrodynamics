#pragma once
//this is used often enough, and MSVC can't understand a 0-sized array ...

namespace Hydrodynamics {

template<int n>
struct MinOneMinusOne {
	static constexpr auto value = n - 1;
};

template<>
struct MinOneMinusOne<1> {
	static constexpr auto value = 1;
};

}
