#pragma once

#include "Hydro/Cell.h"
#include "Hydro/Interface.h"
#include "TensorMath/Grid.h"
#include "TensorMath/Vector.h"
#include "Parallel.h"

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
	virtual void resize(int width, int height) = 0;
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
	typedef typename Interface::StateMatrix StateMatrix;
	typedef typename Interface::StateInverseMatrix StateInverseMatrix;
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
	CellGrid cells;
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
	virtual void resize(int width, int height);
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
, cells(size_+1)		
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

	Parallel::For(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		for (int j = 0; j < rank; ++j) {
			cell.x(j) = xmin(j) + (xmax(j) - xmin(j)) * (Real)index(j) / (Real)(size(j)-1);
		}
	});

	Parallel::For(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		InterfaceVector &interface = v.second.interfaces;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {	//which side
			if (index(side) < 1 || index(side) >= size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexR = index;
				IVector indexL = index;
				--indexL(side);
				for (int k = 0; k < rank; ++k) {
					interface(side).x(k) = (cells(indexR).x(k) + cells(indexL).x(k)) * Real(.5);
				}
			}
		}
	});

	Parallel::For(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		InterfaceVector &interface = v.second.interfaces;
		for (int k = 0; k < rank; ++k) {
			//extrapolate based on which edge it is
			if (index(k) == size(k)) {
				IVector indexL = index;
				--indexL(k);
				IVector indexL2 = indexL;
				--indexL2(k);
				for (int j = 0; j < 3; ++j) {
					interface(k).x(j) = 2. * cells(indexL).interfaces(k).x(j) - cells(indexL2).interfaces(k).x(j);
				}
			} else if (index(k) == 0) {
				IVector indexR = index;
				++indexR(k);
				IVector indexR2 = indexR;
				++indexR2(k);			
				for (int j = 0; j < 3; ++j) {
					interface(k).x(j) = 2. * cells(indexR).interfaces(k).x(j) - cells(indexR2).interfaces(k).x(j);
				}
			}
		}
	});
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::boundary() {
	PROFILE()
	(*boundaryMethod)(this);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::step(Real dt) {
	PROFILE()
	boundary();
	solver->step(this, dt);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::update() {
	PROFILE()
	
	solver->initStep(this);

	Real dt = useCFL 
		? solver->calcCFLTimestep(this)
		: fixedDT;

	step(dt);
}

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::getPrimitives() {
	PROFILE()
	Parallel::For(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		cell.primitives = equationOfState->getPrimitives(cell.state);
	});
}

template<int rank>
struct HydroResize {
	static void resize(int width, int height) {
		const float zNear = 1;
		const float zFar = 1000;
		float aspectRatio = (float)width / (float)height;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-aspectRatio * zNear, aspectRatio * zNear, -zNear, zNear, zNear, zFar);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0,0,-3);
	}
};

template<>
struct HydroResize<1> {
	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glOrtho(-aspectRatio, aspectRatio, -.25, .5, -1., 1.);
	}
};

template<>
struct HydroResize<2> {
	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glOrtho(-aspectRatio, aspectRatio, -1., 1., -1., 1.);
	}
};

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::resize(int width, int height) {
	HydroResize<rank>::resize(width, height);
}

template<typename Real, int rank> void plotVertex(::Tensor<Real, Upper<rank> > x, Real value);

template<> void plotVertex<float, 1>(::Tensor<float, Upper<1> > x, float value) { glVertex2f(x(0), value); }
template<> void plotVertex<double, 1>(::Tensor<double, Upper<1> > x, double value) { glVertex2d(x(0), value); }
template<> void plotVertex<float, 2>(::Tensor<float, Upper<2> > x, float value) { glVertex3f(x(0), x(1), value); }
template<> void plotVertex<double, 2>(::Tensor<double, Upper<2> > x, double value) { glVertex3d(x(0), x(1), value); }
template<> void plotVertex<float, 3>(::Tensor<float, Upper<3> > x, float value) { glVertex4f(x(0), x(1), value, x(2)); }
template<> void plotVertex<double, 3>(::Tensor<double, Upper<3> > x, double value) { glVertex4d(x(0), x(1), value, x(2)); }

template<int rank>
struct HydroPlot {
	template<typename Hydro>
	static void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		glBegin(GL_TRIANGLE_STRIP);
		std::for_each(hydro.cells.begin(), hydro.cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			Cell &cell = v.second;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < hydro.nghost-1 || index(side) >= hydro.size(side) - hydro.nghost) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				const float plotScalar = .1;

				//color by state value, neglect height or use it for coordinates
				glColor3f(0, cell.state(0), 0);	//color by density
				plotVertex<Real, rank>(cell.x, plotScalar * cell.primitives(0));		
			}
		});
		glEnd();
	}
};

//1D case
template<>
struct HydroPlot<1> {
	template<typename Hydro>
	static void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::StateVector StateVector;
		typedef typename Hydro::IVector IVector;
		const double plotScalar = .1;
		for (int state = 0; state < 3; ++state) {
			::Vector<float,3> color;
			color(state) = 1;
			glColor3fv(color.v);
			for (int i = 0; i < hydro.size(0); ++i) {
				Cell &cell = hydro.cells(IVector(i));
				StateVector leftPrims = hydro.equationOfState->getPrimitives(cell.stateLeft(0));
				StateVector rightPrims = hydro.equationOfState->getPrimitives(cell.stateRight(0));
				glBegin(GL_LINE_STRIP);
				glVertex2d(cell.interfaces(0).x(0), plotScalar * leftPrims(state));
				glVertex2d(cell.x(0), plotScalar * cell.primitives(state));
				glVertex2d(hydro.cells(i+1).interfaces(0).x(0), plotScalar * rightPrims(state));
				glEnd();
			}
		}
	}
};

//2D case
template<>
struct HydroPlot<2> {
	template<typename Hydro>
	static void draw(Hydro &hydro) {
		static GLuint texID = 0;
		static bool init = false;
		if (!init) {
			init = true;
			glGenTextures(1, &texID);
			glBindTexture(GL_TEXTURE_1D, texID);

			float colors[][3] = {
				{.5, .5, .5},
				{1, 1, 0},
				{1, .5, 0},
				{1, 0, 0}
			};

			const int width = 256;
			float data[width*3];
			for (int i = 0; i < width; ++i) {
				::Vector<float,3> c;
				float f = (float)i / (float)width * (float)numberof(colors);
				int ci = (int)f;
				float s = f - (float)ci;

				for (int j = 0; j < 3; ++j) {
					data[3 * i + j] = colors[ci][j] * (1.f - s) + colors[ci+1][j] * s;
				}
			}

			glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, width, 0, GL_RGB, GL_FLOAT, data);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
		
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		glEnable(GL_TEXTURE_1D);
		for (int y = hydro.nghost-1; y < hydro.size(1)-1; ++y) {
			glBegin(GL_TRIANGLE_STRIP);
			for (int x = hydro.nghost-1; x < hydro.size(0); ++x) {
			
				for (int offset = 0; offset < 2; ++offset) {
					int index = x + hydro.cells.size(0) * (y + offset);
					Cell &cell = hydro.cells.v[index].second;

					//color by state value, neglect height or use it for coordinates
					glTexCoord1f(1.f - cell.state(0));
					glVertex2d(cell.x(0), cell.x(1));
				}
			}
			glEnd();
		}
		glDisable(GL_TEXTURE_1D);
	}
};

template<typename Real, int rank, typename EquationOfState>
void Hydro<Real, rank, EquationOfState>::draw() {
	PROFILE()
	getPrimitives();
	HydroPlot<rank>::template draw<Hydro>(*this);
}

