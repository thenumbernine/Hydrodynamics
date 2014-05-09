#pragma once

#include <math.h>
#include <algorithm>

template<typename Real>
class FluxMethod {
public:
	virtual Real operator()(Real) = 0;
};

template<typename Real>
class DonorCellFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class LaxWendroffFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class BeamWarmingFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class FrommFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class CHARMFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class HCUSFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class HQUICKFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class KorenFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class MinModFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class OshkerFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class OspreFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class SmartFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class SwebyFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class UMISTFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class VanAlbada1FluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class VanAlbada2FluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class VanLeerFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class MonotonizedCentralFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class SuperbeeFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real>
class BarthJespersenFluxMethod : public FluxMethod<Real> {
public:
	virtual Real operator()(Real);
};

template<typename Real> Real DonorCellFluxMethod<Real>::operator()(Real r) { return Real(0.); }
template<typename Real> Real LaxWendroffFluxMethod<Real>::operator()(Real r) { return Real(1.); }
template<typename Real> Real BeamWarmingFluxMethod<Real>::operator()(Real r) { return r; }
template<typename Real> Real FrommFluxMethod<Real>::operator()(Real r) { return Real(.5) * (Real(1.) + r); }
template<typename Real> Real CHARMFluxMethod<Real>::operator()(Real r) {if (r < Real(0.)) return Real(0.); return r*(Real(3.)*r+Real(1.))/((r+Real(1.))*(r+Real(1.))); }
template<typename Real> Real HCUSFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), Real(1.5) * (r + fabs(r)) / (r + 2.) ); }
template<typename Real> Real HQUICKFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), Real(2.) * (r + fabs(r)) / (r + Real(3.)) ); }
template<typename Real> Real KorenFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.)*r, std::min<Real>((Real(1.) + Real(2.)*r)/Real(3.), Real(2.)))); }
template<typename Real> Real MinModFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(r,Real(1.)) ); }
template<typename Real> Real OshkerFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(r,Real(1.5)) ); }	//replace 1.5 with 1 <= beta <= 2	
template<typename Real> Real OspreFluxMethod<Real>::operator()(Real r) { return Real(.5) * (r*r + r) / (r*r + r + Real(.1)); }
template<typename Real> Real SmartFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.) * r, std::min<Real>(Real(.25) + Real(.75) * r, Real(4.)))); }
template<typename Real> Real SwebyFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::max<Real>(std::min<Real>(Real(1.5) * r, Real(1.)), std::min<Real>(r, Real(1.5)))); }	//replace 1.5 with 1 <= beta <= 2
template<typename Real> Real UMISTFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(std::min<Real>(Real(2.)*r, Real(.75) + Real(.25)*r), std::min<Real>(Real(.25) + Real(.75)*r, Real(2.)))); }	
template<typename Real> Real VanAlbada1FluxMethod<Real>::operator()(Real r) { return (r * r + r) / (r * r + Real(1.)); }
template<typename Real> Real VanAlbada2FluxMethod<Real>::operator()(Real r) { return Real(2.) * r / (r * r + Real(1.)); }
template<typename Real> Real VanLeerFluxMethod<Real>::operator()(Real r) { return (r + fabs(r)) / (Real(1.) + fabs(r)); }
template<typename Real> Real MonotonizedCentralFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::min<Real>(Real(2.), std::min<Real>(Real(.5) * (Real(1.) + r), Real(2.) * r))); }
template<typename Real> Real SuperbeeFluxMethod<Real>::operator()(Real r) { return std::max<Real>(Real(0.), std::max<Real>(std::min<Real>(Real(1.), Real(2.)*r), std::min<Real>(Real(2.), r))); }
template<typename Real> Real BarthJespersenFluxMethod<Real>::operator()(Real r) { return Real(.5) * (r + Real(1.)) * std::min<Real>(Real(1.), std::min<Real>(Real(4.)*r/(r+Real(1.)), Real(4.)/(r+Real(1.)))); }
