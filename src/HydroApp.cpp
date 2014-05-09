#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <OpenGL/gl.h>
#include <SDL/SDL.h>

#include "GLApp/GLApp.h" 

#include "Common/Exception.h"

#include "Hydro/Hydro.h"

#include "Hydro/EquationOfState/EulerEquationOfState.h"

#include "Hydro/InitialConditions/SodInitialConditions.h"
#include "Hydro/InitialConditions/SedovInitialConditions.h"
#include "Hydro/InitialConditions/AdvectInitialConditions.h"
#include "Hydro/InitialConditions/WaveInitialConditions.h"

#include "Hydro/BoundaryMethod/PeriodicBoundaryMethod.h"
#include "Hydro/BoundaryMethod/MirrorBoundaryMethod.h"
#include "Hydro/BoundaryMethod/ConstantBoundaryMethod.h"
#include "Hydro/BoundaryMethod/FreeFlowBoundaryMethod.h"

#include "Hydro/Solver/EulerEquationBurgersSolverExplicit.h"
#include "Hydro/Solver/EulerEquationGodunovSolverExplicit.h"
//#include "Hydro/Solver/EulerEquationRoeSolverExplicit.h"

#include "Hydro/ExplicitMethod/ForwardEulerExplicitMethod.h"
#include "Hydro/ExplicitMethod/RungeKutta2ExplicitMethod.h"
#include "Hydro/ExplicitMethod/RungeKutta4ExplicitMethod.h"
#include "Hydro/ExplicitMethod/IterativeCrankNicolson3ExplicitMethod.h"

#include "Hydro/FluxMethod.h"

class HydroArgs {
public:
	int dim;
	int size;
	bool useCFL;
	double cfl;
	double fixedDT;
	double gamma;
	std::string precision;
	std::string boundaryMethodName;
	std::string equationOfStateName;
	std::string solverName;
	std::string explicitMethodName;
	std::string fluxMethodName;
	std::string initialConditionsName;
	HydroArgs() 
	: dim(1)
	, size(256)
	, useCFL(true)
	, cfl(.5)
	, fixedDT(.1)
	, gamma(1.4)
	, precision("double")
	, boundaryMethodName("Mirror")
	, equationOfStateName("Euler")
	, solverName("EulerEquationGodunovSolverExplicit")//("EulerEquationRoeSolverExplicit")
	, explicitMethodName("ForwardEuler")
	, fluxMethodName("Superbee")
	, initialConditionsName("Sod")
	{}
};


class HydroApp : public GLApp {
	IHydro *hydro;
	HydroArgs args;
	int dim;

public:
	HydroApp()
	: hydro(NULL)
	, dim(1)
	{}

	virtual int main(int argc, char **argv) {
		for (int i = 0; i < argc; ++i) {
			if (i < argc-1) {
				if (!strcmp(argv[i], "--initialConditions")) {
					args.initialConditionsName = argv[++i];
				} else if (!strcmp(argv[i], "--boundaryMethod")) {
					args.boundaryMethodName = argv[++i];
				} else if (!strcmp(argv[i], "--equationOfState")) {
					args.equationOfStateName = argv[++i];
				} else if (!strcmp(argv[i], "--solver")) {
					args.solverName = argv[++i];
				} else if (!strcmp(argv[i], "--explicitMethod")) {
					args.explicitMethodName = argv[++i];
				} else if (!strcmp(argv[i], "--fluxMethod")) {
					args.fluxMethodName = argv[++i];
				} else if (!strcmp(argv[i], "--size")) {
					args.size = atoi(argv[++i]);
				} else if (!strcmp(argv[i], "--useCFL")) {
					args.useCFL = !strcmp(argv[++i], "true") ? true : false;
				} else if (!strcmp(argv[i], "--cfl")) {
					args.cfl = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--fixedDT")) {
					args.fixedDT = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--gamma")) {
					args.gamma = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--dim")) {
					dim = atoi(argv[++i]);
				} else if (!strcmp(argv[i], "--precision")) {
					args.precision = argv[++i];
				}
			}
		}

		return GLApp::main(argc, argv);
	}

	template<typename Real, int rank, typename EquationOfState>
	void initType() {
		typedef ::Hydro<Real, rank, EquationOfState> Hydro;
		typedef ::Solver<Real> Solver;
		typedef ::ExplicitMethod<Hydro> ExplicitMethod;
		typedef ::FluxMethod<Real> FluxMethod;

		EquationOfState *equationOfState = new EquationOfState();

		InitialConditions *initialConditions = NULL;
		if (args.initialConditionsName == "Sod") {
			initialConditions = new SodInitialConditions<Real, rank, EquationOfState>();
		} else if (args.initialConditionsName == "Sedov") {
			initialConditions = new SedovInitialConditions<Real, rank, EquationOfState>();
		} else if (args.initialConditionsName == "Advect") {
			initialConditions = new AdvectInitialConditions<Real, rank, EquationOfState>();
		} else if (args.initialConditionsName == "Wave") {
			initialConditions = new WaveInitialConditions<Real, rank, EquationOfState>();
		} else {
			throw Exception() << "unknown initial conditions " << args.initialConditionsName;
		}

		BoundaryMethod *boundaryMethod = NULL;
		if (args.boundaryMethodName == "Periodic") {
			boundaryMethod = new PeriodicBoundaryMethod<Hydro>();
		} else if (args.boundaryMethodName =="Mirror") {
			boundaryMethod = new MirrorBoundaryMethod<Hydro>();
		} else if (args.boundaryMethodName =="Constant") {
			boundaryMethod = new ConstantBoundaryMethod<Hydro>();
		} else if (args.boundaryMethodName =="FreeFlow") {
			boundaryMethod = new FreeFlowBoundaryMethod<Hydro>();
		} else {
			throw Exception() << "unknown boundary method " << args.boundaryMethodName;
		}

		Solver *solver = NULL;
		if (args.solverName == "EulerEquationBurgersSolverExplicit") {
			solver = new EulerEquationBurgersSolverExplicit<Hydro>();
		} else if (args.solverName == "EulerEquationGodunovSolverExplicit") {
			solver = new EulerEquationGodunovSolverExplicit<Hydro>();
		//} else if (args.solverName == "EulerEquationRoeSolverExplicit") {
		//	solver = new EulerEquationRoeSolverExplicit<Hydro>();
		} else {
			throw Exception() << "unknown solver " << args.solverName;
		}

		ExplicitMethod *explicitMethod = NULL;
		if (args.explicitMethodName == "ForwardEuler") {
			explicitMethod = new ForwardEulerExplicitMethod<Hydro>();
		} else if (args.explicitMethodName == "RungeKutta2") {
			explicitMethod = new RungeKutta2ExplicitMethod<Hydro>();
		} else if (args.explicitMethodName == "RungeKutta4") {
			explicitMethod = new RungeKutta4ExplicitMethod<Hydro>();
		} else if (args.explicitMethodName == "IterativeCrankNicolson3") {
			explicitMethod = new IterativeCrankNicolson3ExplicitMethod<Hydro>();
		} else {
			throw Exception() << "unknown explicit method " << args.explicitMethodName;
		}

		FluxMethod *fluxMethod = NULL;
		if (args.fluxMethodName == "DonorCell") {
			fluxMethod = new DonorCellFluxMethod<Real>();
		} else if (args.fluxMethodName == "LaxWendroff") {
			fluxMethod = new LaxWendroffFluxMethod<Real>();
		} else if (args.fluxMethodName == "BeamWarming") {
			fluxMethod = new BeamWarmingFluxMethod<Real>();
		} else if (args.fluxMethodName == "Fromm") {
			fluxMethod = new FrommFluxMethod<Real>();
		} else if (args.fluxMethodName == "CHARM") {
			fluxMethod = new CHARMFluxMethod<Real>();
		} else if (args.fluxMethodName == "HCUS") {
			fluxMethod = new HCUSFluxMethod<Real>();
		} else if (args.fluxMethodName == "HQUICK") {
			fluxMethod = new HQUICKFluxMethod<Real>();
		} else if (args.fluxMethodName == "Koren") {
			fluxMethod = new KorenFluxMethod<Real>();
		} else if (args.fluxMethodName == "MinMod") {
			fluxMethod = new MinModFluxMethod<Real>();
		} else if (args.fluxMethodName == "Oshker") {
			fluxMethod = new OshkerFluxMethod<Real>();
		} else if (args.fluxMethodName == "Ospre") {
			fluxMethod = new OspreFluxMethod<Real>();
		} else if (args.fluxMethodName == "Smart") {
			fluxMethod = new SmartFluxMethod<Real>();
		} else if (args.fluxMethodName == "Sweby") {
			fluxMethod = new SwebyFluxMethod<Real>();
		} else if (args.fluxMethodName == "UMIST") {
			fluxMethod = new UMISTFluxMethod<Real>();
		} else if (args.fluxMethodName == "VanAlbada1") {
			fluxMethod = new VanAlbada1FluxMethod<Real>();
		} else if (args.fluxMethodName == "VanAlbada2") {
			fluxMethod = new VanAlbada2FluxMethod<Real>();
		} else if (args.fluxMethodName == "VanLeer") {
			fluxMethod = new VanLeerFluxMethod<Real>();
		} else if (args.fluxMethodName == "MonotonizedCentral") {
			fluxMethod = new MonotonizedCentralFluxMethod<Real>();
		} else if (args.fluxMethodName == "Superbee") {
			fluxMethod = new SuperbeeFluxMethod<Real>();
		} else if (args.fluxMethodName == "BarthJespersen") {
			fluxMethod = new BarthJespersenFluxMethod<Real>();
		} else {
			throw Exception() << "unknown flux method " << args.fluxMethodName;
		}

		typedef ::Vector<int, rank> IVector;
		hydro = new Hydro(
			IVector(args.size),
			args.useCFL,
			args.cfl,
			args.fixedDT,
			args.gamma,
			boundaryMethod,
			equationOfState,
			solver,
			explicitMethod,
			fluxMethod,
			initialConditions);
	}

	template<typename Real, int rank>
	void initPrecision() {
		if (args.equationOfStateName == "Euler") {
			initType<Real, rank, EulerEquationOfState<Real, rank> >();
		} else {
			throw Exception() << "unknown equation of state " << args.equationOfStateName;
		}
	}

	template<int rank>
	void initSize() {
		if (args.precision == "single") {
			initPrecision<float,rank>();
		} else if (args.precision == "double") {
			initPrecision<double,rank>();
		} else {
			throw Exception() << "unknown precision " << args.precision;
		}
	}

	virtual void init() {
		GLApp::init();
		switch (args.dim) {
		case 1:
			initSize<1>();
			break;
#if 0	// fix the boundary condition static asserts before enabling these
		case 2:
			initSize<2>();
			break;
		case 3:
			initSize<3>();
			break;
#endif
		default:
			throw Exception() << "unknown dim " << args.dim;
		}
	}

	virtual void resize(int width, int height) {
		GLApp::resize(width, height);
		float aspectRatio = (float)width / (float)height;	
		glOrtho(-aspectRatio, aspectRatio, -1, 1, -1, 1);
	}
	
	virtual void update() {
		GLApp::update();
		hydro->update();
		hydro->draw();
	}
};
GLAPP_MAIN(HydroApp)

