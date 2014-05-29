#include "Profiler/Profiler.h"	//placed at the top so everyone can use it without including it because I am lazy like that
#include "Common/Macros.h"		//similar laziness
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
#include "GLApp/GLApp.h" 
#include "Common/Exception.h"
#include <OpenGL/gl.h>
#include <SDL2/SDL.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

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
	: GLApp()
	, ihydro(NULL)
	{}

	virtual int main(int argc, char **argv) {
		for (int i = 0; i < argc; ++i) {
			if (!strcmp(argv[i], "--help")) {
				std::cout << "usage: hydro <args>" << std::endl;
				std::cout << "args:" << std::endl;
				std::cout << "  --equationOfState <equationOfState>" << std::endl;
				std::cout << "    can be one of the following: (Euler)" << std::endl;
				std::cout << "  --initialConditions <initialConditions>" << std::endl;
				std::cout << "    can be one of the following:" << std::endl;
				std::cout << "      for Euler equation of state:" << std::endl;
				std::cout << "        (Sod) Sedov Advect Wave KelvinHemholtz RayleighTaylor" << std::endl;
				std::cout << "  --boundaryMethod <boundaryMethod>" << std::endl;
				std::cout << "    can be one of the following: (Mirror) Peroiod" << std::endl;
				std::cout << "  --solver <solver>" << std::endl;
				std::cout << "    can be one of the following:" << std::endl;
				std::cout << "    for Euler equation of state: Burgers Godunov (Roe)" << std::endl;
				std::cout << "  --explicitMethod <explicitMethod>" << std::endl;
				std::cout << "    can be one of the following:" << std::endl;
				std::cout << "      ForwardEuler RungeKutta2 (RungeKutta4) IterativeCrankNicolson3" << std::endl;
				std::cout << "  --fluxMethod <fluxMethod>" << std::endl;
				std::cout << "    can be one of the following:" << std::endl;
				std::cout << "      DonorCell LaxWendroff BeamWarming Fromm CHARM HCUS HQUICK Koren MinMod" << std::endl;
				std::cout << "      Oshker Ospre Smart Sweby UMIST VanAlbada1 VanAlbada2 VanLeer" << std::endl;
				std::cout << "      MonotonizedCentral (Superbee) BarthJespersen" << std::endl;
				std::cout << "  --dim <dim>" << std::endl;
				std::cout << "    the dimension, can be 1, 2, or 3.  default " << args.dim << std::endl;
				std::cout << "  --size <size1> <size2> ... <sizeN>" << std::endl;
				std::cout << "    the grid size, for N = the dimension of the grid.  default " << args.size << std::endl;
				std::cout << "  --useCFL <true|false> = whether to use CFL or a fixed timestep.  default " << args.useCFL << std::endl;
				std::cout << "  --cfl <CFL> = the CFL number.  default" << args.cfl << std::endl;
				std::cout << "  --fixedDT <dt> = the fixed timestep. default " << args.fixedDT << std::endl;
				std::cout << "  --noise <noise> = noise amplitude to apply to initial velocity.  default " << args.noise << std::endl;
				std::cout << "  --precision <precision> = precision: single double.  default: " << args.precision << std::endl;
				return 1;
			}
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
		typedef ::Hydro<EquationOfState> Hydro;
		typedef ::ISolver<Real> ISolver;
		typedef ::ExplicitMethod<Hydro> ExplicitMethod;
		typedef ::FluxMethod<Real> FluxMethod;

		EquationOfState *equationOfState = new EquationOfState();

		//these are eos-specific
		InitialConditions<Real, rank> *initialConditions = equationOfState->initialConditions.create(args.initialConditionsName);
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

		hydro->resetCoordinates(initialConditions->xmin, initialConditions->xmax);
		
		//after external force is determined, recalibrate potential energy
		//check the minimum of the potential energy at all corners
		typedef typename Hydro::IVector IVector;
		RangeObj<rank> range(IVector(), IVector(2));
		hydro->minPotentialEnergy = HUGE_VAL;
		std::for_each(
			range.begin(), 
			range.end(),
			[&](IVector index)
		{
			Real energyPotential = 0.; 
			::Vector<Real, rank> x;
			for (int k = 0; k < rank; ++k) {
				x(k) = index(k) ? hydro->xmax(k) : hydro->xmin(k);
				energyPotential += x(k) * hydro->externalForce(k);
			}
			std::cout << " corner " << index << " potential " << energyPotential << " corner " << x << " extrnal force " << hydro->externalForce << std::endl;
			hydro->minPotentialEnergy = std::min<Real>(hydro->minPotentialEnergy, energyPotential);
		});
		//add its negative to all potential energy calculations
		// to keep potential energy from being a negative value
	std::cout << " min potential " << hydro->minPotentialEnergy << std::endl;
		hydro->minPotentialEnergy = 0;//-hydro->minPotentialEnergy;
		
		//once min potential energy is determined, set up initial conditions
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

	virtual void sdlEvent(SDL_Event &event) {
		static bool leftButtonDown = false;
		static bool leftShiftDown = false;
		static bool rightShiftDown = false;
		bool shiftDown = leftShiftDown | rightShiftDown;

		switch (event.type) {
		case SDL_MOUSEMOTION:
			{	
				int idx = event.motion.xrel;
				int idy = event.motion.yrel;
				if (leftButtonDown) {
					if (shiftDown) {
						if (idy) {
							ihydro->zoom(idy);
						} 
					} else {
						if (idx || idy) {
							ihydro->pan(idx, idy);
						}
					}
				}
			}
			break;
		case SDL_MOUSEBUTTONDOWN:
			if (event.button.button == SDL_BUTTON_LEFT) {
				leftButtonDown = true;
			}
			break;
		case SDL_MOUSEBUTTONUP:
			if (event.button.button == SDL_BUTTON_LEFT) {
				leftButtonDown = false;
			}
			break;
		case SDL_KEYDOWN:
			if (event.key.keysym.sym == SDLK_LSHIFT) {
				leftShiftDown = true;
			} else if (event.key.keysym.sym == SDLK_RSHIFT) {
				rightShiftDown = true;
			}
			break;
		case SDL_KEYUP:
			if (event.key.keysym.sym == SDLK_LSHIFT) {
				leftShiftDown = false;
			} else if (event.key.keysym.sym == SDLK_RSHIFT) {
				rightShiftDown = false;
			}
			break;
		}
	}
};
GLAPP_MAIN(HydroApp)

