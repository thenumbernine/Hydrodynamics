#pragma once

#include "Hydro/IHydro.h"
#include "Hydro/Plot.h"
#include "Hydro/Cell.h"
#include "Hydro/Interface.h"
#include "Hydro/Solver/ISolver.h"
#include "Hydro/Explicit/Explicit.h"
#include "Hydro/Boundary/Boundary.h"
#include "Hydro/Limiter.h"
#include "Tensor/Grid.h"
#include "Tensor/Vector.h"
#include "Parallel/Parallel.h"

template<typename Equation_>
struct Hydro : public IHydro {
	typedef Equation_ Equation;
	
	typedef typename Equation::Real Real;
	enum { rank = Equation::rank };

	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::Explicit::Explicit<Hydro> Explicit;
	typedef ::Limiter::Limiter<Real> Limiter;
	typedef typename ::DisplayMethod<Hydro> DisplayMethod;
	
	enum { numberOfStates = Equation::numberOfStates };
	typedef ::Cell<Real, rank, numberOfStates> Cell;
	typedef ::Interface<Real, rank, numberOfStates> Interface;
	
	typedef Tensor::Vector<int, rank> IVector;
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
	std::shared_ptr<::Boundary::Boundary> boundaryMethod;
	std::shared_ptr<Equation> equation;
	std::shared_ptr<ISolver> solver;
	std::shared_ptr<Explicit> explicitMethod;
	std::shared_ptr<Limiter> limiter;
	std::shared_ptr<DisplayMethod> displayMethod;

public:	//'til I can work out access
	int nghost;
	Vector xmin, xmax;

	typedef Tensor::Grid<std::pair<IVector, Cell>, rank> CellGrid;
	typedef Tensor::Vector<Interface, rank> InterfaceVector;
	CellGrid cells;

	Plot<rank> plot;
public:
	Hydro(IVector size_,
		bool useCFL_,
		Real cfl_,
		Real fixedDT_,
		Real gamma_,
		std::shared_ptr<::Boundary::Boundary> boundaryMethod_,
		std::shared_ptr<Equation> equation_,
		std::shared_ptr<ISolver> solver_,
		std::shared_ptr<Explicit> explicitMethod_,
		std::shared_ptr<Limiter> limiter_,
		std::shared_ptr<DisplayMethod> displayMethod_);
	
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

template<typename Equation>
Hydro<Equation>::Hydro(IVector size_,
	bool useCFL_,
	Real cfl_,
	Real fixedDT_,
	Real gamma_,
	std::shared_ptr<::Boundary::Boundary> boundaryMethod_,
	std::shared_ptr<Equation> equation_,
	std::shared_ptr<ISolver> solver_,
	std::shared_ptr<Explicit> explicitMethod_,
	std::shared_ptr<Limiter> limiter_,
	std::shared_ptr<DisplayMethod> displayMethod_)
: size(size_)
, useCFL(useCFL_)
, cfl(cfl_)
, fixedDT(fixedDT_)
, gamma(gamma_)
, minPotentialEnergy(0)
, boundaryMethod(boundaryMethod_)
, equation(equation_)
, solver(solver_)
, explicitMethod(explicitMethod_)
, limiter(limiter_)
, displayMethod(displayMethod_)
, nghost(2) 
, cells(size_+1)
{
	Tensor::RangeObj<rank> range(IVector(), cells.size);
	for (IVector index : range) {
		typename CellGrid::Type &v = cells(index);
		v.first = index;
	}

	resetCoordinates(xmin, xmax);
}

template<typename Equation>
void Hydro<Equation>::resetCoordinates(Vector xmin_, Vector xmax_) {
	xmin = xmin_;
	xmax = xmax_;

	Parallel::parallel->foreach(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		for (int j = 0; j < rank; ++j) {
			cell.x(j) = xmin(j) + (xmax(j) - xmin(j)) * (Real)index(j) / (Real)(size(j)-1);
		}
	});

	Parallel::parallel->foreach(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
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

	Parallel::parallel->foreach(cells.begin(), cells.end(), [&](typename CellGrid::value_type &v) {
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

template<typename Equation>
void Hydro<Equation>::boundary() {
	PROFILE()
	(*boundaryMethod)(this);
}

template<typename Equation>
void Hydro<Equation>::step(Real dt) {
	PROFILE()
	boundary();
	solver->step(this, dt);
}

template<typename Equation>
void Hydro<Equation>::update() {
	PROFILE()

	solver->initStep(this);

	Real dt = useCFL 
		? solver->calcCFLTimestep(this)
		: fixedDT;

	step(dt);
}

template<typename Equation>
void Hydro<Equation>::resize(int width, int height) {
	plot.resize(width, height);
}

template<typename Equation>
void Hydro<Equation>::pan(int dx, int dy) {
	plot.pan(dx, dy);
}

template<typename Equation>
void Hydro<Equation>::zoom(int dz) {
	plot.zoom(dz);
}

template<typename Equation>
void Hydro<Equation>::draw() {
	PROFILE()
	plot.template draw<Hydro>(*this, displayMethod);
}

