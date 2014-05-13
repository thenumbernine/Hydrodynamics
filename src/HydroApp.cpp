#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>

#include <OpenGL/gl.h>
#include <SDL/SDL.h>

#define numberof(x) (sizeof(x)/sizeof((x)[0]))
#define endof(x) ((x)+numberof(x))
#define crand()	((double)rand()/(double)RAND_MAX)
#include "GLApp/GLApp.h" 

#include "Common/Exception.h"

#include "Profile.h"

#include "Hydro/Hydro.h"

#include "Hydro/EquationOfState/EulerEquationOfState.h"

#include "Hydro/BoundaryMethod/MirrorBoundaryMethod.h"
#include "Hydro/BoundaryMethod/PeriodicBoundaryMethod.h"
//TODO get these working for all dimensions
//#include "Hydro/BoundaryMethod/ConstantBoundaryMethod.h"
//#include "Hydro/BoundaryMethod/FreeFlowBoundaryMethod.h"

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
	double noise;
	double fixedDT;
	double gamma;
	std::vector<double> externalForce;	//relies on dim ... hmm ... reason to use config file over arguments?
	std::string precision;
	std::string boundaryMethodName;
	std::string equationOfStateName;
	std::string solverName;
	std::string explicitMethodName;
	std::string fluxMethodName;
	std::string initialConditionsName;
	HydroArgs() 
	: dim(2)
	, size(256)
	, useCFL(true)
	, cfl(.5)
	, noise(.01)
	, fixedDT(.1)
	, gamma(1.4)
	, precision("double")
	, boundaryMethodName("Mirror")
	, equationOfStateName("Euler")
	, solverName("Roe")
	, explicitMethodName("ForwardEuler")
	, fluxMethodName("Superbee")
	, initialConditionsName("Sod")
	{}
};


class HydroApp : public GLApp {
	IHydro *ihydro;
	HydroArgs args;

public:
	HydroApp()
	: ihydro(NULL)
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
				} else if (!strcmp(argv[i], "--noise")) {
					args.noise = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--fixedDT")) {
					args.fixedDT = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--gamma")) {
					args.gamma = atof(argv[++i]);
				} else if (!strcmp(argv[i], "--dim")) {
					if (args.externalForce.size()) throw Exception() << "you must set dim before you set externalForce";
					args.dim = atoi(argv[++i]);
				} else if (!strcmp(argv[i], "--precision")) {
					args.precision = argv[++i];
				}
			}
			if (i < argc-args.dim) {	//dim must be set first
				if (!strcmp(argv[i], "--externalForce")) {
					args.externalForce.resize(args.dim);
					for (int k = 0; k < args.dim; ++k) {
						args.externalForce[k] = atof(argv[++i]);
					}
				}
			}
		}

		return GLApp::main(argc, argv);
	}

	template<typename Real, int rank, typename EquationOfState>
	void initType() {
		typedef ::Hydro<Real, rank, EquationOfState> Hydro;
		typedef ::ISolver<Real> ISolver;
		typedef ::ExplicitMethod<Hydro> ExplicitMethod;
		typedef ::FluxMethod<Real> FluxMethod;

		EquationOfState *equationOfState = new EquationOfState();

		//these are eos-specific
		InitialConditions *initialConditions = equationOfState->initialConditions.create(args.initialConditionsName);
		ISolver *solver = equationOfState->solvers.create(args.solverName);

		BoundaryMethod *boundaryMethod = NULL;
		if (args.boundaryMethodName =="Mirror") {
			boundaryMethod = new MirrorBoundaryMethod<Hydro>();
		} else if (args.boundaryMethodName == "Periodic") {
			boundaryMethod = new PeriodicBoundaryMethod<Hydro>();
#if 0
		} else if (args.boundaryMethodName =="Constant") {
			boundaryMethod = new ConstantBoundaryMethod<Hydro>();
		} else if (args.boundaryMethodName =="FreeFlow") {
			boundaryMethod = new FreeFlowBoundaryMethod<Hydro>();
#endif
		} else {
			throw Exception() << "unknown boundary method " << args.boundaryMethodName;
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
		Hydro *hydro = new Hydro(
			IVector(args.size),
			args.useCFL,
			args.cfl,
			args.fixedDT,
			args.gamma,
			boundaryMethod,
			equationOfState,
			solver,
			explicitMethod,
			fluxMethod);
		for (int i = 0; i < rank && i < args.externalForce.size(); ++i) {
			hydro->externalForce(i) = args.externalForce[i];
		}
		ihydro = hydro;
		(*initialConditions)(ihydro, args.noise);
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
		
		//not sure where to put this yet:
		{
			static GLuint texID = 0;

			glGenTextures(1, &texID);
			glBindTexture(GL_TEXTURE_1D, texID);

			float colors[][3] = {
				{0, 0, 0},
				{0, 0, .5},
				{1, .5, 0},
				{1, 0, 0}
			};

			const int width = 256;
			unsigned char data[width*3];
			for (int i = 0; i < width; ++i) {
				::Vector<float,3> c;
				float f = (float)i / (float)width * (float)numberof(colors);
				int ci = (int)f;
				float s = f - (float)ci;
				if (ci >= numberof(colors)) {
					ci = numberof(colors)-1;
					s = 0;
				}

				for (int j = 0; j < 3; ++j) {
					data[3 * i + j] = (unsigned char)(255. * (colors[ci][j] * (1.f - s) + colors[ci+1][j] * s));
				}
			}

			glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, width, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		}
		
		switch (args.dim) {
		case 1:
			initSize<1>();
			break;
		case 2:
			initSize<2>();
			break;
		case 3:
			initSize<3>();
			break;
		default:
			throw Exception() << "unknown dim " << args.dim;
		}
	}

	virtual void resize(int width, int height) {
		GLApp::resize(width, height);
		ihydro->resize(width, height);
	}
	
	virtual void update() {
		PROFILE_BEGIN_FRAME()	
		{
			PROFILE()

			GLApp::update();
			ihydro->update();
			ihydro->draw();
		}	
		PROFILE_END_FRAME()
	}

	virtual void shutdown() {
		PROFILE_DONE()
	}
};
GLAPP_MAIN(HydroApp)

