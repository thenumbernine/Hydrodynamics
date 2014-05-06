#pragma once

class FluxMethod {
public:
	virtual double operator()(double) = 0;
};

class DonorCellFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class LaxWendroffFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class BeamWarmingFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class FrommFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class CHARMFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class HCUSFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class HQUICKFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class KorenFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class MinModFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class OshkerFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class OspreFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class SmartFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class SwebyFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class UMISTFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class VanAlbada1FluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class VanAlbada2FluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class VanLeerFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class MonotonizedCentralFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class SuperbeeFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

class BarthJespersenFluxMethod : public FluxMethod {
public:
	virtual double operator()(double);
};

