#pragma once

#include "Hydro/Cell.h"
#include "Hydro/Interface.h"
#include "Hydro/ISolver.h"
#include "TensorMath/Grid.h"
#include "TensorMath/Vector.h"
#include "Quat.h"
#include "Parallel.h"

#include <OpenGL/gl.h>


template<int rank>
struct HydroPlot {
	Quat viewAngle;
	float dist;

	HydroPlot() : dist(2.) {}

	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::StateVector StateVector;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0,0,-dist);
		Quat angleAxis = viewAngle.toAngleAxis();
		glRotatef(angleAxis(3) * 180. / M_PI, angleAxis(0), angleAxis(1), angleAxis(2));
		
		glEnable(GL_TEXTURE_1D);
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE);
		glBegin(GL_POINTS);
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
				//color by state value, neglect height or use it for coordinates
				glTexCoord1f(2. * cell.state(0));	//color by density
				glVertex3d(cell.x(0), cell.x(1), cell.x(2));
			}
		});
		glEnd();
		glDisable(GL_TEXTURE_1D);
	}

	static void resize(int width, int height) {
		const float zNear = .01;
		const float zFar = 10;
		float aspectRatio = (float)width / (float)height;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-aspectRatio * zNear, aspectRatio * zNear, -zNear, zNear, zNear, zFar);
	}

	void pan(int dx, int dy) {
		float magn = sqrt(dx * dx + dy * dy);
		float fdx = (float)dx / magn;
		float fdy = (float)dy / magn;
		Quat rotation = Quat(fdy, fdx, 0, magn * M_PI / 180.).fromAngleAxis();
		viewAngle = rotation * viewAngle;
		viewAngle /= Quat::length(viewAngle);
	}

	void zoom(int dz) {
		dist *= (float)exp((float)dz * -.03f);
	}
};

//1D case
template<>
struct HydroPlot<1> {
	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::StateVector StateVector;
		typedef typename Hydro::IVector IVector;
#if 0	//show piecewise step functions - in anticipation of getting PPM method working	
		for (int state = 0; state < 3; ++state) {
			::Vector<float,3> color;
			color(state) = 1;
			glColor3fv(color.v);
			for (int i = 0; i < hydro.size(0); ++i) {
				Cell &cell = hydro.cells(IVector(i));
				StateVector primitivesLeft = hydro.equationOfState->getPrimitives(cell.stateLeft(0));
				StateVector primitives = hydro.equationOfState->getPrimitives(cell.state);
				StateVector primitivesRight = hydro.equationOfState->getPrimitives(cell.stateRight(0));
				glBegin(GL_LINE_STRIP);
				glVertex2d(cell.interfaces(0).x(0), primitivesLeft(state));
				glVertex2d(cell.x(0), primitives(state));
				glVertex2d(hydro.cells(i+1).second.interfaces(0).x(0), primitivesRight(state));
				glEnd();
			}
		}
#endif
#if 1	//good ol fashioned graph
		for (int state = 0; state < 3; ++state) {
			::Vector<float,3> color;
			color(state) = 1;
			glColor3fv(color.v);
			glBegin(GL_LINE_STRIP);
			std::for_each(hydro.cells.begin(), hydro.cells.end(), [&](typename CellGrid::value_type &v) {
				Cell &cell = v.second;
				StateVector primitives = hydro.equationOfState->getPrimitives(cell.state);
				glVertex2d(cell.x(0), primitives(state));
			});
			glEnd();
		}

#endif
	}

	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glOrtho(-aspectRatio, aspectRatio, -1., 3., -1., 1.);
	}

	void pan(int dx, int dy) {
	}
	
	void zoom(int dz) {}
};

//2D case
template<>
struct HydroPlot<2> {
	::Vector<float, 2> viewPos;
	float viewZoom;

	HydroPlot() : viewZoom(1.) {}

	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		glPushMatrix();
		glTranslatef(-viewPos(0), -viewPos(1), 0);
		glScalef(viewZoom, viewZoom, viewZoom);
		glEnable(GL_TEXTURE_1D);
		for (int y = hydro.nghost-1; y < hydro.size(1)-1; ++y) {
			glBegin(GL_TRIANGLE_STRIP);
			for (int x = hydro.nghost-1; x < hydro.size(0); ++x) {
			
				for (int offset = 0; offset < 2; ++offset) {
					int index = x + hydro.cells.size(0) * (y + offset);
					Cell &cell = hydro.cells.v[index].second;

					//color by state value, neglect height or use it for coordinates
					glTexCoord1f(2. * cell.state(0));
					glVertex2d(cell.x(0), cell.x(1));
				}
			}
			glEnd();
		}
		glDisable(GL_TEXTURE_1D);
		glPopMatrix();
	}

	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glOrtho(-aspectRatio, aspectRatio, -1., 1., -1., 1.);
	}

	void pan(int dx, int dy) {
		viewPos(0) -= (float)dx * 0.01f;
		viewPos(1) += (float)dy * 0.01f;
	}
	
	void zoom(int dz) {
		viewZoom *= exp((float)dz * -.03f); 
	}
};



struct BoundaryMethod;

template<typename Real>
struct ISolver;

template<typename Hydro>
struct ExplicitMethod;

template<typename Real>
struct FluxMethod;

struct IHydro {
	virtual void update() = 0;
	virtual void draw() = 0;
	virtual void resize(int width, int height) = 0;
	virtual void pan(int dx, int dy) = 0;
	virtual void zoom(int dz) = 0;
};

template<typename EquationOfState_>
struct Hydro : public IHydro {
	typedef EquationOfState_ EquationOfState;
	
	typedef typename EquationOfState::Real Real;
	enum { rank = EquationOfState::rank };

	typedef ISolver<Real> ISolver;
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
	Real minPotentialEnergy;
	BoundaryMethod *boundaryMethod;
	EquationOfState *equationOfState;
	ISolver *solver;
	ExplicitMethod *explicitMethod;
	FluxMethod *fluxMethod;

public:	//'til I can work out access
	int nghost;
	Vector xmin, xmax;

	typedef ::Grid<std::pair<IVector, Cell>, rank> CellGrid;
	typedef ::Vector<Interface, rank> InterfaceVector;
	CellGrid cells;

	HydroPlot<rank> plot;
public:
	Hydro(IVector size_,
		bool useCFL_,
		Real cfl_,
		Real fixedDT_,
		Real gamma_,
		BoundaryMethod *boundaryMethod_,
		EquationOfState *equationOfState_,
		ISolver *solver_,
		ExplicitMethod *explicitMethod_,
		FluxMethod *fluxMethod_);
	
	void resetCoordinates(Vector xmin_, Vector xmax_);
	void step(Real dt);
	void boundary();
public:
	virtual void update();
	virtual void draw();
	virtual void resize(int width, int height);
	virtual void pan(int dx, int dy);
	virtual void zoom(int dz);
};

template<typename EquationOfState>
Hydro<EquationOfState>::Hydro(IVector size_,
	bool useCFL_,
	Real cfl_,
	Real fixedDT_,
	Real gamma_,
	BoundaryMethod *boundaryMethod_,
	EquationOfState *equationOfState_,
	ISolver *solver_,
	ExplicitMethod *explicitMethod_,
	FluxMethod *fluxMethod_)
: minPotentialEnergy(0)
, nghost(2) 
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

	resetCoordinates(xmin, xmax);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::resetCoordinates(Vector xmin_, Vector xmax_) {
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
					interface(side).x(k) = (cells(indexR).second.x(k) + cells(indexL).second.x(k)) * Real(.5);
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
					interface(k).x(j) = 2. * cells(indexL).second.interfaces(k).x(j) - cells(indexL2).second.interfaces(k).x(j);
				}
			} else if (index(k) == 0) {
				IVector indexR = index;
				++indexR(k);
				IVector indexR2 = indexR;
				++indexR2(k);			
				for (int j = 0; j < 3; ++j) {
					interface(k).x(j) = 2. * cells(indexR).second.interfaces(k).x(j) - cells(indexR2).second.interfaces(k).x(j);
				}
			}
		}
	});
}

template<typename EquationOfState>
void Hydro<EquationOfState>::boundary() {
	PROFILE()
	(*boundaryMethod)(this);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::step(Real dt) {
	PROFILE()
	boundary();
	solver->step(this, dt);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::update() {
	PROFILE()
	
	solver->initStep(this);

	Real dt = useCFL 
		? solver->calcCFLTimestep(this)
		: fixedDT;

	step(dt);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::resize(int width, int height) {
	plot.resize(width, height);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::pan(int dx, int dy) {
	plot.pan(dx, dy);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::zoom(int dz) {
	plot.zoom(dz);
}

template<typename EquationOfState>
void Hydro<EquationOfState>::draw() {
	PROFILE()
	plot.template draw<Hydro>(*this);
}

