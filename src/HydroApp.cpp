#ifdef DISABLE_PROFILER
#define PROFILE()
#define PROFILE_BEGIN_FRAME()
#define PROFILE_END_FRAME()
#else
#include "Profiler/Profiler.h"	//placed at the top so everyone can use it without including it because I am lazy like that
#endif

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
#include "GLApp/GLApp.h"
#include "GLCxx/Texture.h"
#include "GLCxx/gl.h"
#include "Common/Exception.h"
#include "SDL.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <thread>

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

std::shared_ptr<Parallel> parallel;

struct HydroApp : public GLApp::GLApp {
	using Super = ::GLApp::GLApp;

	std::shared_ptr<IHydro> ihydro;
	HydroArgs hydroArgs;

	virtual std::string getTitle() { return "Hydrodynamics"; }

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
//printf("here\n"); exit(0); // works

	//these are eos-specific
	std::cout << "initialConditions = " << hydroArgs.initialConditionsName << std::endl;	
	std::shared_ptr<InitialConditions::InitialConditions<Real, rank>> initialConditions = equation->initialConditions[hydroArgs.initialConditionsName]();

	std::cout << "solver = " << hydroArgs.solverName << std::endl;	
	std::shared_ptr<ISolver> solver = equation->solvers[hydroArgs.solverName]();

#define ALLOCATOR(BaseType, Name, Type) {#Name, []() -> std::shared_ptr<BaseType> { return std::make_shared<Type>(); }}

	std::shared_ptr<Boundary::Boundary> boundaryMethod;
	{
		AllocatorMap<Boundary::Boundary> boundaryMethods = {
			ALLOCATOR(Boundary::Boundary, Mirror, Boundary::Mirror<Hydro>),
			ALLOCATOR(Boundary::Boundary, Periodic, Boundary::Periodic<Hydro>),
			//ALLOCATOR(Constant, Boundary::Boundary, Boundary::Constant<Hydro>),
			//ALLOCATOR(FreeFlow, Boundary::Boundary, Boundary::FreeFlow<Hydro>),
		};
		std::cout << "boundary = " << hydroArgs.boundaryName << std::endl;	
		boundaryMethod = boundaryMethods[hydroArgs.boundaryName]();
	}

	std::shared_ptr<Explicit> explicitMethod;
	{
		AllocatorMap<Explicit> explicitMethods = {
			ALLOCATOR(Explicit, ForwardEuler, Hydrodynamics::Explicit::ForwardEuler<Hydro>),
			ALLOCATOR(Explicit, RungeKutta2, Hydrodynamics::Explicit::RungeKutta2<Hydro>),
			ALLOCATOR(Explicit, RungeKutta4, Hydrodynamics::Explicit::RungeKutta4<Hydro>),
			ALLOCATOR(Explicit, IterativeCrankNicolson3, Hydrodynamics::Explicit::IterativeCrankNicolson3<Hydro>),
		};
		std::cout << "explicit = " << hydroArgs.explicitName << std::endl;	
		explicitMethod = explicitMethods[hydroArgs.explicitName]();
	}

	std::shared_ptr<Limiter> limiter; 
	{
		AllocatorMap<Limiter> limiters = {
			ALLOCATOR(Limiter, DonorCell, Hydrodynamics::Limiter::DonorCell<Real>),
			ALLOCATOR(Limiter, LaxWendroff, Hydrodynamics::Limiter::LaxWendroff<Real>),
			ALLOCATOR(Limiter, BeamWarming, Hydrodynamics::Limiter::BeamWarming<Real>),
			ALLOCATOR(Limiter, Fromm, Hydrodynamics::Limiter::Fromm<Real>),
			ALLOCATOR(Limiter, CHARM, Hydrodynamics::Limiter::CHARM<Real>),
			ALLOCATOR(Limiter, HCUS, Hydrodynamics::Limiter::HCUS<Real>),
			ALLOCATOR(Limiter, HQUICK, Hydrodynamics::Limiter::HQUICK<Real>),
			ALLOCATOR(Limiter, Koren, Hydrodynamics::Limiter::Koren<Real>),
			ALLOCATOR(Limiter, MinMod, Hydrodynamics::Limiter::MinMod<Real>),
			ALLOCATOR(Limiter, Oshker, Hydrodynamics::Limiter::Oshker<Real>),
			ALLOCATOR(Limiter, Ospre, Hydrodynamics::Limiter::Ospre<Real>),
			ALLOCATOR(Limiter, Smart, Hydrodynamics::Limiter::Smart<Real>),
			ALLOCATOR(Limiter, Sweby, Hydrodynamics::Limiter::Sweby<Real>),
			ALLOCATOR(Limiter, UMIST, Hydrodynamics::Limiter::UMIST<Real>),
			ALLOCATOR(Limiter, VanAlbada1, Hydrodynamics::Limiter::VanAlbada1<Real>),
			ALLOCATOR(Limiter, VanAlbada2, Hydrodynamics::Limiter::VanAlbada2<Real>),
			ALLOCATOR(Limiter, VanLeer, Hydrodynamics::Limiter::VanLeer<Real>),
			ALLOCATOR(Limiter, MonotonizedCentral, Hydrodynamics::Limiter::MonotonizedCentral<Real>),
			ALLOCATOR(Limiter, Superbee, Hydrodynamics::Limiter::Superbee<Real>),
			ALLOCATOR(Limiter, BarthJespersen, Hydrodynamics::Limiter::BarthJespersen<Real>),
		};
		std::cout << "limiter = " << hydroArgs.limiterName << std::endl;	
		limiter = limiters[hydroArgs.limiterName]();
	}

	std::shared_ptr<DisplayMethod<Hydro>> displayMethod;
	{
		AllocatorMap<DisplayMethod<Hydro>> displayMethods = {
			ALLOCATOR(DisplayMethod<Hydro>, density, DensityColoring<Hydro>),
			ALLOCATOR(DisplayMethod<Hydro>, velocity, VelocityColoring<Hydro>),
			ALLOCATOR(DisplayMethod<Hydro>, pressure, PressureColoring<Hydro>),
		};
		std::cout << "display = " << hydroArgs.displayName << std::endl;	
		displayMethod = displayMethods[hydroArgs.displayName]();
	}

#undef ALLOCATOR

	using IVector = typename Hydro::IVector;
	IVector sizev;
	for (int i = 0; i < rank; ++i) {
		std::cout << "size("<<i<<") = " << hydroArgs.size[i] << std::endl;	
		sizev(i) = hydroArgs.size[i];
	}

	std::cout << "numThreads = " << hydroArgs.numThreads << std::endl;	
	parallel = std::make_shared<Parallel>(hydroArgs.numThreads);

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
		Tensor::vec<Real, rank> x;
		for (int k = 0; k < rank; ++k) {
			x(k) = index(k) ? hydro->xmax(k) : hydro->xmin(k);
			potentialSpecificEnergy += x(k) * hydro->externalForce(k);
		}
// TODO i've broken operator<< for tensors
//		std::cout << " corner " << index << " potential " << potentialSpecificEnergy << " corner " << x << " external force " << hydro->externalForce << std::endl;
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
			std::cout << "    can be one of the following: (Mirror) Peroiodic" << std::endl;
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

		static auto tex = GLCxx::Texture1D()
			.bind()
			.create1D(width, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, data)
			.setParam<GL_TEXTURE_MAG_FILTER>(GL_LINEAR)
			.setParam<GL_TEXTURE_MIN_FILTER>(GL_NEAREST)
			.setParam<GL_TEXTURE_WRAP_S>(GL_REPEAT);
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
