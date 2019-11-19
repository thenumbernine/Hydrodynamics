#include "Profiler/Profiler.h"	//placed at the top so everyone can use it without including it because I am lazy like that
#include "Common/Macros.h"		//similar laziness
#include "Hydro/Hydro.h"
#include "Hydro/Equation/Euler.h"
//#include "Hydro/Equation/MHD.h"
#include "Hydro/Equation/SRHD.h"
#include "Hydro/Boundary/Mirror.h"
#include "Hydro/Boundary/Periodic.h"
//TODO get these working for all dimensions
//#include "Hydro/Boundary/Constant.h"
//#include "Hydro/Boundary/FreeFlow.h"
#include "Hydro/Explicit/ForwardEuler.h"
#include "Hydro/Explicit/RungeKutta2.h"
#include "Hydro/Explicit/RungeKutta4.h"
#include "Hydro/Explicit/IterativeCrankNicolson3.h"
#include "Hydro/Limiter.h"
#include "GLApp/gl.h"
#include "GLApp/GLApp.h"
#include "Common/Exception.h"
#include "SDL.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

class HydroArgs {
public:
	int dim;
	std::vector<int> size;
	int numThreads;
	bool useCFL;
	double cfl;
	double noise;
	double fixedDT;
	double gamma;

	std::vector<double> externalForce;	//relies on dim ... hmm ... reason to use config file over arguments?
	std::string precision;
	std::string boundaryName;
	std::string equationName;
	std::string solverName;
	std::string explicitName;
	std::string limiterName;
	std::string initialConditionsName;
	std::string displayName;

	HydroArgs()
	: dim(2)
	, size({256,256})
	, numThreads(std::thread::hardware_concurrency())
	, useCFL(true)
	, cfl(.5)
	, noise(0.)
	, fixedDT(.1)
	, gamma(1.4)
	, precision("double")
	, boundaryName("Mirror")
	, equationName("Euler")
	, solverName("Roe")
	, explicitName("ForwardEuler")
	, limiterName("Superbee")
	, initialConditionsName("Sod")
	, displayName("density")
	{}
};

namespace Hydrodynamics {

struct HydroApp : public GLApp::GLApp {
	using Super = ::GLApp::GLApp;

	std::shared_ptr<IHydro> ihydro;
	HydroArgs hydroArgs;

	virtual const char* getTitle() { return "Hydrodynamics"; }

	virtual void init(const Init& args);

	template<typename Real, int rank, typename Equation>
	void initType();

	template<typename Real, int rank>
	void initPrecision();
	
	template<int rank>
	void initSize();
	
	virtual void onResize();
	virtual void onUpdate();
	virtual void onSDLEvent(SDL_Event &event);
};

template<typename Real, int rank, typename Equation>
void HydroApp::initType() {
	using Hydro = Hydrodynamics::Hydro<Equation>;
	using ISolver = Hydrodynamics::Solver::ISolver<Real>;
	using Explicit = Hydrodynamics::Explicit::Explicit<Hydro>;
	using Limiter = Hydrodynamics::Limiter::Limiter<Real>;

	std::shared_ptr<Equation> equation = std::make_shared<Equation>();

	//these are eos-specific
	std::shared_ptr<InitialConditions::InitialConditions<Real, rank>> initialConditions = equation->initialConditions(hydroArgs.initialConditionsName);
	std::shared_ptr<ISolver> solver = equation->solvers(hydroArgs.solverName);

	AllocatorMap<Boundary::Boundary> boundaryMethods;
	boundaryMethods.template add<Boundary::Mirror<Hydro>>("Mirror");
	boundaryMethods.template add<Boundary::Periodic<Hydro>>("Periodic");
	//boundaryMethods.template add<Boundary::Constant<Hydro>>("Constant");
	//boundaryMethods.template add<Boundary::FreeFlow<Hydro>>("FreeFlow");
	std::shared_ptr<Boundary::Boundary> boundaryMethod = boundaryMethods(hydroArgs.boundaryName);

	AllocatorMap<Explicit> explicitMethods;
	explicitMethods.template add<Hydrodynamics::Explicit::ForwardEuler<Hydro>>("ForwardEuler");
	explicitMethods.template add<Hydrodynamics::Explicit::RungeKutta2<Hydro>>("RungeKutta2");
	explicitMethods.template add<Hydrodynamics::Explicit::RungeKutta4<Hydro>>("RungeKutta4");
	explicitMethods.template add<Hydrodynamics::Explicit::IterativeCrankNicolson3<Hydro>>("IterativeCrankNicolson3");
	std::shared_ptr<Explicit> explicitMethod = explicitMethods(hydroArgs.explicitName);

	AllocatorMap<Limiter> limiters;
	limiters.template add<Hydrodynamics::Limiter::DonorCell<Real>>("DonorCell");
	limiters.template add<Hydrodynamics::Limiter::LaxWendroff<Real>>("LaxWendroff");
	limiters.template add<Hydrodynamics::Limiter::BeamWarming<Real>>("BeamWarming");
	limiters.template add<Hydrodynamics::Limiter::Fromm<Real>>("Fromm");
	limiters.template add<Hydrodynamics::Limiter::CHARM<Real>>("CHARM");
	limiters.template add<Hydrodynamics::Limiter::HCUS<Real>>("HCUS");
	limiters.template add<Hydrodynamics::Limiter::HQUICK<Real>>("HQUICK");
	limiters.template add<Hydrodynamics::Limiter::Koren<Real>>("Koren");
	limiters.template add<Hydrodynamics::Limiter::MinMod<Real>>("MinMod");
	limiters.template add<Hydrodynamics::Limiter::Oshker<Real>>("Oshker");
	limiters.template add<Hydrodynamics::Limiter::Ospre<Real>>("Ospre");
	limiters.template add<Hydrodynamics::Limiter::Smart<Real>>("Smart");
	limiters.template add<Hydrodynamics::Limiter::Sweby<Real>>("Sweby");
	limiters.template add<Hydrodynamics::Limiter::UMIST<Real>>("UMIST");
	limiters.template add<Hydrodynamics::Limiter::VanAlbada1<Real>>("VanAlbada1");
	limiters.template add<Hydrodynamics::Limiter::VanAlbada2<Real>>("VanAlbada2");
	limiters.template add<Hydrodynamics::Limiter::VanLeer<Real>>("VanLeer");
	limiters.template add<Hydrodynamics::Limiter::MonotonizedCentral<Real>>("MonotonizedCentral");
	limiters.template add<Hydrodynamics::Limiter::Superbee<Real>>("Superbee");
	limiters.template add<Hydrodynamics::Limiter::BarthJespersen<Real>>("BarthJespersen");
	std::shared_ptr<Limiter> limiter = limiters(hydroArgs.limiterName);

	AllocatorMap<DisplayMethod<Hydro>> displayMethods;
	displayMethods.template add<DensityColoring<Hydro>>("density");
	displayMethods.template add<VelocityColoring<Hydro>>("velocity");
	displayMethods.template add<PressureColoring<Hydro>>("pressure");
	std::shared_ptr<DisplayMethod<Hydro>> displayMethod = displayMethods(hydroArgs.displayName);

	using IVector = typename Hydro::IVector;
	IVector sizev;
	for (int i = 0; i < rank; ++i) {
		sizev(i) = hydroArgs.size[i];
	}

	parallel = std::make_shared<::Parallel::Parallel>(hydroArgs.numThreads);

	std::shared_ptr<Hydro> hydro = std::make_shared<Hydro>(
		sizev,
		hydroArgs.useCFL,
		(Real)hydroArgs.cfl,
		(Real)hydroArgs.fixedDT,
		(Real)hydroArgs.gamma,
		boundaryMethod,
		equation,
		solver,
		explicitMethod,
		limiter,
		displayMethod);
	
	for (int i = 0; i < rank && i < (int)hydroArgs.externalForce.size(); ++i) {
		hydro->externalForce(i) = hydroArgs.externalForce[i];
	}
	ihydro = hydro;

	hydro->resetCoordinates(initialConditions->xmin, initialConditions->xmax);
	
	//after external force is determined, recalibrate potential energy
	//check the minimum of the potential energy at all corners
	Tensor::RangeObj<rank> range(IVector(), IVector(2));
	hydro->minPotentialEnergy = HUGE_VAL;
	for (IVector index : range) {
		Real potentialSpecificEnergy = 0.;
		Tensor::Vector<Real, rank> x;
		for (int k = 0; k < rank; ++k) {
			x(k) = index(k) ? hydro->xmax(k) : hydro->xmin(k);
			potentialSpecificEnergy += x(k) * hydro->externalForce(k);
		}
		std::cout << " corner " << index << " potential " << potentialSpecificEnergy << " corner " << x << " external force " << hydro->externalForce << std::endl;
		hydro->minPotentialEnergy = std::min<Real>(hydro->minPotentialEnergy, potentialSpecificEnergy);
	}
	//add its negative to all potential energy calculations
	// to keep potential energy from being a negative value
std::cout << " min potential " << hydro->minPotentialEnergy << std::endl;
	hydro->minPotentialEnergy = 0;//-hydro->minPotentialEnergy;
	
	//once min potential energy is determined, set up initial conditions
	(*initialConditions)(&*ihydro, hydroArgs.noise);

	//manually call, since the initial onResize() was called before ihydro was assigned
	ihydro->resize(screenSize(0), screenSize(1));
}

template<typename Real, int rank>
void HydroApp::initPrecision() {
	if (hydroArgs.equationName == "Euler") {
		initType<Real, rank, Hydrodynamics::Equation::Euler<Real, rank>>();
	//} else if (hydroArgs.equationName == "MHD") {
	//	initType<Real, rank, Hydrodynamics::Equation::MHD<Real, rank>>();
	} else if (hydroArgs.equationName == "SRHD") {
		initType<Real, rank, Hydrodynamics::Equation::SRHD<Real, rank>>();
	} else {
		throw Common::Exception() << "unknown equation of state " << hydroArgs.equationName;
	}
}

template<int rank>
void HydroApp::initSize() {
	if (hydroArgs.precision == "single") {
		initPrecision<float,rank>();
	} else if (hydroArgs.precision == "double") {
		initPrecision<double,rank>();
	} else {
		throw Common::Exception() << "unknown precision " << hydroArgs.precision;
	}
}

void HydroApp::init(const Init& args) {
	bool setSize = false;
	bool setExternalForce = false;
	for (int i = 1; i < (int)args.size(); ++i) {
		if (args[i] ==  "--help") {
			std::cout << "usage: hydro <args>" << std::endl;
			std::cout << "args:" << std::endl;
			std::cout << "  --equation <equation>" << std::endl;
			std::cout << "    can be one of the following: (Euler)" << std::endl;
			std::cout << "  --initialConditions <initialConditions>" << std::endl;
			std::cout << "    can be one of the following:" << std::endl;
			std::cout << "      for Euler equation of state:" << std::endl;
			std::cout << "        (Sod) Sedov Advect Wave KelvinHelmholtz RayleighTaylor" << std::endl;
			std::cout << "  --boundary <boundary>" << std::endl;
			std::cout << "    can be one of the following: (Mirror) Peroiod" << std::endl;
			std::cout << "  --solver <solver>" << std::endl;
			std::cout << "    can be one of the following:" << std::endl;
			std::cout << "    for Euler equation of state: Burgers Godunov (Roe)" << std::endl;
			std::cout << "  --explicit <explicit>" << std::endl;
			std::cout << "    can be one of the following:" << std::endl;
			std::cout << "      (ForwardEuler) RungeKutta2 RungeKutta4 IterativeCrankNicolson3" << std::endl;
			std::cout << "  --limiter <limiter>" << std::endl;
			std::cout << "    can be one of the following:" << std::endl;
			std::cout << "      DonorCell LaxWendroff BeamWarming Fromm CHARM HCUS HQUICK Koren MinMod" << std::endl;
			std::cout << "      Oshker Ospre Smart Sweby UMIST VanAlbada1 VanAlbada2 VanLeer" << std::endl;
			std::cout << "      MonotonizedCentral (Superbee) BarthJespersen" << std::endl;
			std::cout << "  --display <display>" << std::endl;
			std::cout << "    can be one of the following: (density) velocity pressure" << std::endl;
			std::cout << "  --dim <dim>" << std::endl;
			std::cout << "    the dimension, can be 1, 2, or 3.  default " << hydroArgs.dim << std::endl;
			std::cout << "  --size <size1> <size2> ... <sizeN>" << std::endl;
			std::cout << "    the grid size, for N = the dimension of the grid.  default ";
			const char *comma = "";
			for (int i : hydroArgs.size) { std::cout << comma << i; comma = ", "; }
			std::cout << std::endl;
			std::cout << "  --numThreads <numThreads> = the number of threads to use.  default " << hydroArgs.numThreads << std::endl;
			std::cout << "  --useCFL <true|false> = whether to use CFL or a fixed timestep.  default " << (hydroArgs.useCFL ? "true" : "false") << std::endl;
			std::cout << "  --cfl <CFL> = the CFL number.  default " << hydroArgs.cfl << std::endl;
			std::cout << "  --fixedDT <dt> = the fixed timestep. default " << hydroArgs.fixedDT << std::endl;
			std::cout << "  --noise <noise> = noise amplitude to apply to initial velocity.  default " << hydroArgs.noise << std::endl;
			std::cout << "  --precision <precision> = precision: single double.  default: " << hydroArgs.precision << std::endl;
			requestExit(1);
			return;
		}
		if (i < (int)args.size()-1) {
			if (args[i] == "--initialConditions") {
				hydroArgs.initialConditionsName = args[++i];
				continue;
			} else if (args[i] == "--boundary") {
				hydroArgs.boundaryName = args[++i];
				continue;
			} else if (args[i] == "--equation") {
				hydroArgs.equationName = args[++i];
				continue;
			} else if (args[i] == "--solver") {
				hydroArgs.solverName = args[++i];
				continue;
			} else if (args[i] == "--explicit") {
				hydroArgs.explicitName = args[++i];
				continue;
			} else if (args[i] == "--limiter") {
				hydroArgs.limiterName = args[++i];
				continue;
			} else if (args[i] == "--display") {
				hydroArgs.displayName = args[++i];
				continue;
			} else if (args[i] == "--size") {
				setSize = true;
				hydroArgs.size.resize(hydroArgs.dim);
				for (int k = 0; k < hydroArgs.dim; ++k) {
					hydroArgs.size[k] = std::stoi(args[++i]);
				}
				continue;
			} else if (args[i] == "--numThreads") {
				hydroArgs.numThreads = std::stoi(args[++i]);
				continue;
			} else if (args[i] == "--useCFL") {
				hydroArgs.useCFL = args[++i] == "true" ? true : false;
				continue;
			} else if (args[i] == "--cfl") {
				hydroArgs.cfl = std::stof(args[++i]);
				continue;
			} else if (args[i] == "--noise") {
				hydroArgs.noise = std::stof(args[++i]);
				continue;
			} else if (args[i] == "--fixedDT") {
				hydroArgs.fixedDT = std::stof(args[++i]);
				continue;
			} else if (args[i] == "--gamma") {
				hydroArgs.gamma = std::stof(args[++i]);
				continue;
			} else if (args[i] == "--dim") {
				if (setSize) throw Common::Exception() << "you must set dim before you set size";
				if (setExternalForce) throw Common::Exception() << "you must set dim before you set externalForce";
				hydroArgs.dim = std::stoi(args[++i]);
				continue;
			} else if (args[i] == "--precision") {
				hydroArgs.precision = args[++i];
				continue;
			}
		}
		if (i < (int)args.size()-hydroArgs.dim) {	//dim must be set first
			if (args[i] == "--externalForce") {
				setExternalForce = true;
				hydroArgs.externalForce.resize(hydroArgs.dim);
				for (int k = 0; k < hydroArgs.dim; ++k) {
					hydroArgs.externalForce[k] = std::stof(args[++i]);
				}
				continue;
			}
		}
		throw Common::Exception() << "got unknown cmdline argument: " << args[i];
	}
	
	Super::init(args);

	//not sure where to put this yet:
	{
		static GLuint texID = 0;

		glGenTextures(1, &texID);
		glBindTexture(GL_TEXTURE_1D, texID);

		float colors[][3] = {
			{0,0,0},	// black ... ?
			{0,0,1},	// blue
			{0,1,1},	// cyan
			{0,1,0},	// green
			{1,1,0},	// yellow
			{1,.5,0},	// orange
			{1,0,0},	// red
			{1,1,1},	// white
		};

		const int width = 256;
		unsigned char data[width*3];
		for (int i = 0; i < width; ++i) {
			float f = (float)i / (float)width * (float)numberof(colors);
			int ci = (int)f;
			float s = f - (float)ci;
			if (ci >= (int)numberof(colors)) {
				ci = (int)numberof(colors)-1;
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
	
	switch (hydroArgs.dim) {
	case 1:	//SRHD only works for 3D
		initSize<1>();
		break;
	case 2:
		initSize<2>();
		break;
	case 3:
		initSize<3>();
		break;
	default:
		throw Common::Exception() << "unknown dim " << hydroArgs.dim;
	}
}

void HydroApp::onResize() {
	Super::onResize();
	if (ihydro) ihydro->resize(screenSize(0), screenSize(1));
}
	
void HydroApp::onUpdate() {
	PROFILE_BEGIN_FRAME()
	{
		PROFILE()
		Super::onUpdate();
		ihydro->update();
		ihydro->draw();
	}
	PROFILE_END_FRAME()
}

void HydroApp::onSDLEvent(SDL_Event &event) {
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

}

GLAPP_MAIN(Hydrodynamics::HydroApp)
