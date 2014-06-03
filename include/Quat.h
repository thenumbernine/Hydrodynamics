#pragma once

#include "Tensor/Vector.h"

struct Quat : public Tensor::Vector<float, 4> {
	
	Quat() {
		(*this)(3) = 1.;
	}

	Quat(float x, float y, float z, float w) {
		(*this)(0) = x;
		(*this)(1) = y;
		(*this)(2) = z;
		(*this)(3) = w;
	}

	Quat fromAngleAxis() {
		float c = cos((*this)(3) * .5);
		float n = sqrt((*this)(0) * (*this)(0) + (*this)(1) * (*this)(1) + (*this)(2) * (*this)(2));
		float sn = sin((*this)(3) * .5) / n;
		return Quat(sn * (*this)(0), sn * (*this)(1), sn * (*this)(2), c);
	}
	
	Quat toAngleAxis() {
		float cosHalfAngle = (*this)(3);
		if ((*this)(3) < -1.) (*this)(3) = -1.;
		if ((*this)(3) > 1.) (*this)(3) = 1.;
		float halfAngle = acos(cosHalfAngle);
		float scale = sin(halfAngle);

		if (scale >= -1e-4 && scale <= 1e-4) {
			return Quat(0,0,1,0);
		}

		scale = 1. / scale;
		return Quat(
			(*this)(0) * scale,
			(*this)(1) * scale,
			(*this)(2) * scale,
			2. * halfAngle);
	}

	static Quat mul(Quat &q, Quat &r) {
		float a = (q(3) + q(0)) * (r(3) + r(0));
		float b = (q(2) - q(1)) * (r(1) - r(2));
		float c = (q(0) - q(3)) * (r(1) + r(2));
		float d = (q(1) + q(2)) * (r(0) - r(3));
		float e = (q(0) + q(2)) * (r(0) + r(1));
		float f = (q(0) - q(2)) * (r(0) - r(1));
		float g = (q(3) + q(1)) * (r(3) - r(2));
		float h = (q(3) - q(1)) * (r(3) + r(2));

		return Quat(
			a - .5 * ( e + f + g + h),
			-c + .5 * ( e - f + g - h),
			-d + .5 * ( e - f - g + h),
			b + .5 * (-e - f + g + h));
	}
};

Quat operator*(Quat a, Quat b) {
	return Quat::mul(a,b);
}

