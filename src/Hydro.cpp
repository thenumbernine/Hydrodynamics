#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"
#include "Hydro/BoundaryMethod.h"
#include "Hydro/EquationOfState.h"
#include "Hydro/Solver.h"
#include <OpenGL/gl.h>

#define numberof(x) (sizeof(x)/sizeof((x)[0]))
#define endof(x) ((x)+numberof(x))

Hydro::Hydro(const HydroArgs &args) 
: nghost(2) {
	size = args.size;
	useCFL = args.useCFL;
	cfl = args.cfl;
	fixedDT = args.fixedDT;
	gamma = args.gamma;
	boundaryMethod = args.boundaryMethod;
	equationOfState = args.equationOfState;
	solver = args.solver;
	explicitMethod = args.explicitMethod;
	fluxMethod = args.fluxMethod;
	initialConditions = args.initialConditions;

	std::vector<double> Cell::*cellVectorFields[] = {
		&Cell::state,
		&Cell::primitives,
		&Cell::stateLeft,
		&Cell::stateRight,
		&Cell::dq_dt,
		&Cell::tmpState0,
		&Cell::tmpState1,
		&Cell::tmpState2,
		&Cell::tmpState3,
		&Cell::tmpState4,
	};
	cells.resize(size);
	std::for_each(cells.begin(), cells.end(), [&](Cell &cell) {
		std::for_each(cellVectorFields, endof(cellVectorFields), [&](std::vector<double> Cell::*field) {
			(cell.*field).resize(equationOfState->numberOfStates());
		});
	});
	
	std::vector<double> Interface::*interfaceVectorFields[] = {
		&Interface::r,
		&Interface::flux,
		&Interface::eigenvalues,
		&Interface::stateMid,
		&Interface::rTilde,
		&Interface::deltaQTilde,
	};
	std::vector<std::vector<double> > Interface::*interfaceMatrixFields[] = {
		&Interface::jacobian,
		&Interface::eigenvectors,
		&Interface::eigenvectorsInverse,
	};
	interfaces.resize(size+1);
	std::for_each(interfaces.begin(), interfaces.end(), [&](Interface &interface){
		std::for_each(interfaceVectorFields, endof(interfaceVectorFields), [&](std::vector<double> Interface::*field) {
			(interface.*field).resize(equationOfState->numberOfStates());
		});
		std::for_each(interfaceMatrixFields, endof(interfaceMatrixFields), [&](std::vector<std::vector<double> > Interface::*field) {
			(interface.*field).resize(equationOfState->numberOfStates());
			std::for_each((interface.*field).begin(), (interface.*field).end(), [&](std::vector<double> &col) {
				col.resize(equationOfState->numberOfStates());
			});
		});
	});
	resetCoordinates(xmin, xmax);

	(*args.initialConditions)(this);
}

void Hydro::resetCoordinates(double xmin_, double xmax_) {
	xmin = xmin_;
	xmax = xmax_;

	for (int i = 0; i < size; ++i) {
		cells[i].x = xmin + (xmax - xmin) * (double)i / (double)(size-1);
	}

	for (int i = 1; i <= size; ++i) {
		interfaces[i].x = .5 * (cells[i].x + cells[i-1].x);
	}
	interfaces[0].x = 2. * interfaces[1].x - interfaces[2].x;
	interfaces[size].x = 2. * interfaces[size-1].x - interfaces[size-2].x;
}

void Hydro::boundary() {
	(*boundaryMethod)(this);
}

void Hydro::step(double dt) {
	boundary();
	solver->step(this, dt);
}

void Hydro::update() {
	solver->initStep(this);

	double dt = useCFL 
		? solver->calcCFLTimestep(this)
		: fixedDT;
	
	step(dt);
}

void Hydro::getPrimitives() {
	equationOfState->getPrimitives(this);
}

void Hydro::draw() {
	getPrimitives();
	for (int j = 0; j < equationOfState->numberOfStates(); ++j) {
		float color[3] = {0,0,0};
		color[j] = 1;
		glColor3fv(color);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < size; ++i) {
			glVertex2f( (float)i / (float)(size-1) * 2. - 1., cells[i].primitives[j] * .2);
		}
		glEnd();
	}
}


