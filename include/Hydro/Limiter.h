#pragma once

#include <cmath>
#include <algorithm>

namespace Hydrodynamics {
namespace Limiter {

template<typename Real>
struct Limiter {
	virtual ~Limiter() {}
	virtual Real operator()(Real) = 0;
};

template<typename Real>
struct DonorCell : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct LaxWendroff : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct BeamWarming : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Fromm : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct CHARM : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct HCUS : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct HQUICK : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Koren : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct MinMod : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Oshker : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Ospre : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Smart : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Sweby : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct UMIST : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct VanAlbada1 : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct VanAlbada2 : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct VanLeer : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct MonotonizedCentral : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct Superbee : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real>
struct BarthJespersen : public Limiter<Real> {
	virtual Real operator()(Real);
};

template<typename Real> Real DonorCell<Real>::operator()(Real r) { return Real(0.); }
template<typename Real> Real LaxWendroff<Real>::operator()(Real r) { return Real(1.); }
template<typename Real> Real BeamWarming<Real>::operator()(Real r) { return r; }
template<typename Real> Real Fromm<Real>::operator()(Real r) { return Real(.5) * (Real(1.) + r); }
template<typename Real> Real CHARM<Real>::operator()(Real r) {if (r < Real(0.)) return Real(0.); return r*(Real(3.)*r+Real(1.))/((r+Real(1.))*(r+Real(1.))); }
template<typename Real> Real HCUS<Real>::operator()(Real r) { return std::max<Real>(Real(0.), Real(1.5) * (r + fabs(r)) / (r + 2.) ); }
template<typename Real> Real HQUICK<Real>::operator()(Real r) { return std::max<Real>(Real(0.), Real(2.) * (r + fabs(r)) / (r + Real(3.)) ); }
template<typename Real> Real Koren<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.)*r, std::min<Real>((Real(1.) + Real(2.)*r)/Real(3.), Real(2.)))); }
template<typename Real> Real MinMod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(r,Real(1.)) ); }
template<typename Real> Real Oshker<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(r,Real(1.5)) ); }	//replace 1.5 with 1 <= beta <= 2	
template<typename Real> Real Ospre<Real>::operator()(Real r) { return Real(.5) * (r*r + r) / (r*r + r + Real(.1)); }
template<typename Real> Real Smart<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.) * r, std::min<Real>(Real(.25) + Real(.75) * r, Real(4.)))); }
template<typename Real> Real Sweby<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::max<Real>(std::min<Real>(Real(1.5) * r, Real(1.)), std::min<Real>(r, Real(1.5)))); }	//replace 1.5 with 1 <= beta <= 2
template<typename Real> Real UMIST<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(std::min<Real>(Real(2.)*r, Real(.75) + Real(.25)*r), std::min<Real>(Real(.25) + Real(.75)*r, Real(2.)))); }	
template<typename Real> Real VanAlbada1<Real>::operator()(Real r) { return (r * r + r) / (r * r + Real(1.)); }
template<typename Real> Real VanAlbada2<Real>::operator()(Real r) { return Real(2.) * r / (r * r + Real(1.)); }
template<typename Real> Real VanLeer<Real>::operator()(Real r) { return (r + fabs(r)) / (Real(1.) + fabs(r)); }
template<typename Real> Real MonotonizedCentral<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.), std::min<Real>(Real(.5) * (Real(1.) + r), Real(2.) * r))); }
template<typename Real> Real Superbee<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::max<Real>(std::min<Real>(Real(1.), Real(2.)*r), std::min<Real>(Real(2.), r))); }
template<typename Real> Real BarthJespersen<Real>::operator()(Real r) { return Real(.5) * (r + Real(1.)) * std::min<Real>(Real(1.), std::min<Real>(Real(4.)*r/(r+Real(1.)), Real(4.)/(r+Real(1.)))); }

}
}
