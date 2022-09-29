#pragma once

#include "Tensor/Inverse.h"
#include "Hydro/MinOneMinusOne.h"

namespace Hydrodynamics {

template<int rank>
struct BuildPerpendicularBasis {
	template<typename Real>
	static void go(
		Tensor::_tensor<Real, rank> normal, 
		Tensor::_vec<Tensor::_tensor<Real, rank>, MinOneMinusOne<rank>::value> &tangents) 
	{
		//1) pick normal's max abs component
		//2) fill in all axii but that component
		//3) apply Graham-Schmidt
		// what about coordinate system handedness?

		int maxAxis = -1;
		Real maxValue = -HUGE_VAL;
		for (int k = 0; k < rank; ++k) {
			Real absNormal = (Real)fabs(normal(k));
			if (absNormal > maxValue) {
				maxValue = absNormal;
				maxAxis = k;
			}
		}

		for (int j = 0; j < rank-1; ++j) {
			if (j < maxAxis) {
				tangents(j)(j) = 1;
			} else {
				tangents(j)(j+1) = 1;
			}
		
			for (int k = j-1; k >= 0; --k) {
				Real num = Real(0), denom = Real(0);
				for (int i = 0; i < rank; ++i) {
					num += tangents(j)(i) * tangents(k)(i);
					denom += tangents(j)(i) * tangents(j)(i);
				}
				tangents(j) -= tangents(j) * (num / denom);
			}
			{
				Real num = Real(0), denom = Real(0);
				for (int i = 0; i < rank; ++i) {
					num += tangents(j)(i) * normal(i);
					denom += tangents(j)(i) * tangents(j)(i);
				}
				tangents(j) -= tangents(j) * (num / denom);
			}
			{
				Real len = Real(0);
				for (int k = 0; k < rank; ++k) {
					len += sqrt(tangents(j)(k) * tangents(j)(k));
				}
				tangents(j) *= Real(1) / len;
			}
		}
	}
};

template<>
struct BuildPerpendicularBasis<1> {
	template<typename Real>
	static void go(
		Tensor::_tensor<Real, 1> normal, 
		//the 2nd param should be ,0, but MSVC complains if I allocate a size-zero array
		// so the array has to be size 1, so the pattern has to be 1 here as well ...
		Tensor::_vec<Tensor::_tensor<Real, 1>, 1> &tangents)
	{
	}
};

template<>
struct BuildPerpendicularBasis<2> {
	template<typename Real>
	static void go(
		Tensor::_tensor<Real, 2> normal, 
		Tensor::_vec<Tensor::_tensor<Real, 2>, 1> &tangents) 
	{
		tangents(0)(0) = -normal(1);
		tangents(0)(1) = normal(0);
	}
};

template<typename OutputType, typename InputType>
struct InverseCramer {
	static OutputType go(InputType input) {
		return Tensor::inverse<InputType>(input);
	}
};

//from Numerical Reicipes 2nd ed:
template<typename OutputType, typename InputType>
struct InverseGaussJordan {
	static OutputType go(InputType input) {
		static_assert(InputType::rank == 2, "input must be a matrix!");
		static_assert(OutputType::rank == 2, "output must be a matrix!");
		static_assert(InputType::template dim<0> == InputType::template dim<1>, "input must be square!"); 
		static_assert(OutputType::template dim<0> == OutputType::template dim<1>, "output must be square!"); 
		static_assert(InputType::template dim<0> == OutputType::template dim<0>, "input and output dimensions must match!");
		
		static constexpr auto dim = InputType::template dim<0>;
		using Real = typename InputType::ScalarType;

		int indxc[dim], indxr[dim], ipiv[dim];
		int i, icol = 0, irow = 0, j, k, l, ll;
		Real big, dum, pivinv, temp;
		for (j = 0; j < dim; ++j) {
			ipiv[j] = 0;
		}
		for (i = 0; i < dim; ++i) {
			big = Real(0);
			for ( j = 0; j < dim; ++j) {
				if (ipiv[j] != 1) {
					for (k = 0; k < dim; ++k) {
						if (ipiv[k] == 0) {
							if ((Real)fabs(input(j, k)) >= big) {
								big = (Real)fabs(input(j, k));
								irow = j;
								icol = k;
							}
						} else {
							if(ipiv[k] > 1) {
								throw Common::Exception() << "singular!";
							}
						}
					}
				}
			}
			++ipiv[icol];
			if(irow != icol) {
				for (int k = 0; k < dim; ++k) {
					temp = input(irow, k);
					input(irow, k) = input(icol, k);
					input(icol, k) = temp;
				}
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if(input(icol, icol) == Real(0)) {
				throw Common::Exception() << "singular!";
			}
			pivinv = Real(1) / input(icol, icol);
			input(icol, icol) = Real(1);
			for (l = 0; l < dim; ++l) {
				input(icol, l) *= pivinv;
			}
			for (ll = 0; ll < dim; ++ll) {
				if (ll != icol) {
					dum = input(ll, icol);
					input(ll, icol) = Real(0);
					for (l = 0; l < dim; ++l) {
						input(ll, l) -= input(icol, l) * dum;
					}
				}
			}
		}
		
		for (l = dim-1; l >= 0; --l) {
			if (indxr[l] != indxc[l]) {
				for (k = 0; k < dim; ++k) {
					temp = input(k, indxr[l]);
					input(k, indxr[l]) = input(k, indxc[l]);
					input(k, indxc[l]) = temp;
				}
			}
		}

		OutputType output;
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				output(i,j) = input(i,j);
			}
		}
		return output;
	}
};

}
