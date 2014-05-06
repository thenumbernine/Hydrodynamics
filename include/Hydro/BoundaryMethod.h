#pragma once

class Hydro;
class BoundaryMethod {
public:
	virtual void operator()(Hydro *hydro) = 0;
};

class PeriodicBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(Hydro *hydro);
};

class MirrorBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(Hydro *hydro);
};

class ConstantBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(Hydro *hydro);
};

class FreeFlowBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(Hydro *hydro);
};

