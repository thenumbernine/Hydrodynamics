#pragma once

#include "Hydro/Cell.h"
#include "Hydro/Interface.h"
#include "TensorMath/Grid.h"
#include "TensorMath/Vector.h"

class BoundaryMethod;

template<typename Real>
class Solver;

template<typename Hydro>
class ExplicitMethod;

template<typename Real>
class FluxMethod;

class InitialConditions;

class IHydro {
public:
	virtual void update() = 0;
	virtual void draw() = 0;
};

template<typename Real_, int rank_, typename EquationOfState_>
class Hydro : public IHydro {
public:
	typedef Real_ Real;
	enum { rank = rank_ };

	typedef EquationOfState_ EquationOfState;
	typedef Solver<Real> Solver;
	typedef ExplicitMethod<Hydro> ExplicitMethod;
	typedef FluxMethod<Real> FluxMethod;
	
	enum { numberOfStates = EquationOfState::numberOfStates };
	typedef ::Cell<Real, rank, numberOfStates> Cell;
	typedef ::Interface<Real, rank, numberOfStates> Interface;
	
	typedef ::Vector<int, rank> IVector;
	typedef typename Cell::Vector Vector;
	typedef typename Cell::StateVector StateVector;
public:	//hydro args
	IVector size;
	bool useCFL;
	Real cfl;
	Real fixedDT;
	Real gamma;
	Vector externalForce;
	BoundaryMethod *boundaryMethod;
	EquationOfState *equationOfState;
	Solver *solver;
	ExplicitMethod *explicitMethod;
	FluxMethod *fluxMethod;
	InitialConditions *initialConditions;

public:	//'til I can work out access
	int nghost;
	Vector xmin, xmax;

	typedef ::Grid<Cell, rank> CellGrid;
	typedef ::Vector<Interface, rank> InterfaceVector;
	typedef ::Grid<InterfaceVector, rank> InterfaceGrid;
	CellGrid cells;	//size^n
	InterfaceGrid interfaces; //(size+1)^n * n
public:
	Hydro(IVector size_,
		bool useCFL_,
		Real cfl_,
		Real fixedDT_,
		Real gamma_,
		BoundaryMethod *boundaryMethod_,
		EquationOfState *equationOfState_,
		Solver *solver_,
		ExplicitMethod *explicitMethod_,
		FluxMethod *fluxMethod_,
		InitialConditions *initialConditions_);
	
	void resetCoordinates(Vector xmin_, Vector xmax_);
	void step(Real dt);
	void boundary();
	void getPrimitives();
public:
	virtual void update();
	virtual void draw();
};

#include "Hydro/Solver.h"
#include <OpenGL/gl.h>

#define numberof(x) (sizeof(x)/sizeof((x)[0]))
#define endof(x) ((x)+numberof(x))

template<typename Real, int rank, typename EquationOfState>
Hydro<Real, rank, EquationOfState>::Hydro(IVector size_,
	bool useCFL_,
	Real cfl_,
	Real fixedDT_,
	Real gamma_,
	BoundaryMethod *boundaryMethod_,
	EquationOfState *equationOfState_,
	Solver *solver_,
	ExplicitMethod *explicitMethod_,
	FluxMethod *fluxMethod_,
	InitialConditions *initialConditions_)
: nghost(2) 
, cells(size_)
, interfaces(size_+1)
{
	size = size_;
	useCFL = useCFL_;
	cfl = cfl_;
	fixedDT = fixedDT_;
	gamma = gamma_;
	boundaryMethod = boundaryMethod_;
	equationOfState = equationOfState_;
	solver = solver_;
	explicitMethod = explicitMethod_;
	fluxMethod = fluxMethod_;
	initialConditions = initialConditions_;

	resetCoordinates(xmin, xmax);

	(*initialConditions)(this);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::resetCoordinates(Vector xmin_, Vector xmax_) {
	xmin = xmin_;
	xmax = xmax_;

	for (typename CellGrid::iterator i = cells.begin(); i != cells.end(); ++i) {
		for (int j = 0; j < rank; ++j) {
			i->x(j) = xmin(j) + (xmax(j) - xmin(j)) * (Real)i.index(j) / (Real)(size(j)-1);
		}
	}

	for (typename InterfaceGrid::iterator i = interfaces.begin(); i != interfaces.end(); ++i) {
		InterfaceVector &interface = *i;
		for (int k = 0; k < rank; ++k) {	//which side
			if (!(i.index(k) == 0 || i.index(k) == size(k))) {
				IVector nextIndex = i.index;
				IVector prevIndex = nextIndex;
				--prevIndex(k);
				
				for (int j = 0; j < 3; ++j) {	//which element of the position vector
					interface(k).x(j) = .5 * (cells(nextIndex).x(j) + cells(prevIndex).x(j));
				}
			}
		}
	}
	for (typename InterfaceGrid::iterator i = interfaces.begin(); i != interfaces.end(); ++i) {
		InterfaceVector &interface = *i;
		for (int k = 0; k < rank; ++k) {
			//extrapolate based on which edge it is
			if (i.index(k) == size(k)) {
				IVector prevIndex = i.index;
				--prevIndex(k);
				IVector prev2Index = prevIndex;
				--prev2Index(k);
				for (int j = 0; j < 3; ++j) {
					interface(k).x(j) = 2. * interfaces(prevIndex)(k).x(j) - interfaces(prev2Index)(k).x(j);
				}
			} else if (i.index(k) == 0) {
				IVector nextIndex = i.index;
				++nextIndex(k);
				IVector next2Index = nextIndex;
				++next2Index(k);
				for (int j = 0; j < 3; ++j) {
					interface(k).x(j) = 2. * interfaces(nextIndex)(k).x(j) - interfaces(next2Index)(k).x(j);
				}
			}
		}
	}
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::boundary() {
	(*boundaryMethod)(this);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::step(Real dt) {
	boundary();
	solver->step(this, dt);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::update() {
	solver->initStep(this);

	Real dt = useCFL 
		? solver->calcCFLTimestep(this)
		: fixedDT;
	
	step(dt);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::getPrimitives() {
	std::for_each(cells.begin(), cells.end(), [&](Cell &cell) {
		ICell *icell = dynamic_cast<ICell*>(&cell);
		equationOfState->getPrimitives(icell);
	});
}

template<typename Real, int rank> void plotVertex(::Tensor<Real, Upper<rank> > x, Real value);

template<> void plotVertex<float, 1>(::Tensor<float, Upper<1> > x, float value) { glVertex2f(x(0), value); }
template<> void plotVertex<double, 1>(::Tensor<double, Upper<1> > x, double value) { glVertex2d(x(0), value); }
template<> void plotVertex<float, 2>(::Tensor<float, Upper<2> > x, float value) { glVertex3f(x(0), x(1), value); }
template<> void plotVertex<double, 2>(::Tensor<double, Upper<2> > x, double value) { glVertex3d(x(0), x(1), value); }
template<> void plotVertex<float, 3>(::Tensor<float, Upper<3> > x, float value) { glVertex4f(x(0), x(1), value, x(2)); }
template<> void plotVertex<double, 3>(::Tensor<double, Upper<3> > x, double value) { glVertex4d(x(0), x(1), value, x(2)); }

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::draw() {
	getPrimitives();

	static ::Vector<::Vector<float, 3>, numberOfStates> colors;
	static bool initColors = false;
	if (!initColors) {
		initColors = true;
		for (int i = 0; i < numberOfStates; ++i) {
			float lenSq = 0.;
			for (int j = 0; j < 3; ++j) {
				colors(i)(j) = (float)rand() / (float)RAND_MAX;
				lenSq += colors(i)(j) * colors(i)(j);
			}
			float len = sqrt(lenSq);
			for (int j = 0; j < 3; ++j) {
				colors(i)(j) /= len;
			}
		}
	}

	const Real plotScalar = .1;
	glBegin(GL_LINES);
	for (typename CellGrid::iterator i = cells.begin(); i != cells.end(); ++i) {
		Cell &cell = *i;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (i.index(side) < nghost-1 || i.index(side) >= size(side) - nghost) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				Vector dx(0.);
				dx(side) = (xmax(side) - xmin(side)) / (Real)size(side);
				IVector indexL = i.index;
				--indexL(side);
				
				for (int state = 0; state < numberOfStates; ++state) {
					glColor3fv(colors(state).v);
					plotVertex<Real, rank>(cell.x - dx, plotScalar * cells(indexL).primitives(state));
					plotVertex<Real, rank>(cell.x, plotScalar * cell.primitives(state));
				}
			}
		}
	}
	glEnd();
}

