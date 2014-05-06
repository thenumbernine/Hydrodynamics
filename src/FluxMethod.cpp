#include "Hydro/FluxMethod.h"

#include <math.h>
#include <algorithm>

double DonorCellFluxMethod::operator()(double r) { return 0.; }
double LaxWendroffFluxMethod::operator()(double r) { return 1.; }
double BeamWarmingFluxMethod::operator()(double r) { return r; }
double FrommFluxMethod::operator()(double r) { return .5 * (1. + r); }
double CHARMFluxMethod::operator()(double r) {if (r < 0.) return 0.; return r*(3.*r+1.)/((r+1.)*(r+1.)); }
double HCUSFluxMethod::operator()(double r) { return std::max<double>(0., 1.5 * (r + fabs(r)) / (r + 2.) ); }
double HQUICKFluxMethod::operator()(double r) { return std::max<double>(0., 2. * (r + fabs(r)) / (r + 3.) ); }
double KorenFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(2.*r, std::min<double>((1. + 2.*r)/3., 2.))); }
double MinModFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(r,1.) ); }
double OshkerFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(r,1.5) ); }	//replace 1.5 with 1 <= beta <= 2	
double OspreFluxMethod::operator()(double r) { return .5 * (r*r + r) / (r*r + r + .1); }
double SmartFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(2. * r, std::min<double>(.25 + .75 * r, 4.))); }
double SwebyFluxMethod::operator()(double r) { return std::max<double>(0., std::max<double>(std::min<double>(1.5 * r, 1.), std::min<double>(r, 1.5))); }	//replace 1.5 with 1 <= beta <= 2
double UMISTFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(std::min<double>(2.*r, .75 + .25*r), std::min<double>(.25 + .75*r, 2.))); }	
double VanAlbada1FluxMethod::operator()(double r) { return (r * r + r) / (r * r + 1.); }
double VanAlbada2FluxMethod::operator()(double r) { return 2. * r / (r * r + 1.); }
double VanLeerFluxMethod::operator()(double r) { return (r + fabs(r)) / (1. + fabs(r)); }
double MonotonizedCentralFluxMethod::operator()(double r) { return std::max<double>(0., std::min<double>(2., std::min<double>(.5 * (1. + r), 2. * r))); }
double SuperbeeFluxMethod::operator()(double r) { return std::max<double>(0., std::max<double>(std::min<double>(1., 2.*r), std::min<double>(2., r))); }
double BarthJespersenFluxMethod::operator()(double r) { return .5 * (r + 1.) * std::min<double>(1., std::min<double>(4.*r/(r+1.), 4./(r+1.))); }
