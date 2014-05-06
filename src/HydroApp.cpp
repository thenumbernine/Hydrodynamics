#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <map>
#include <string>

#include <OpenGL/gl.h>
#include <SDL/SDL.h>

#include "GLApp/GLApp.h" 
#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"
#include "Hydro/BoundaryMethod.h"
#include "Hydro/EquationOfState.h"
#include "Hydro/Solver.h"
#include "Hydro/ExplicitMethod.h"
#include "Hydro/FluxMethod.h"
#include "Common/Exception.h"

class Allocator {
public:
	virtual void *create() = 0;
};

template<typename SubType>
class TypeAllocator : public Allocator {
public:
	virtual void *create() { return new SubType(); }
};

template<typename Type> struct AllocatorMap {
	static std::map<std::string, Allocator*> allocators;
};
template<> std::map<std::string, Allocator*> AllocatorMap<InitialConditions>::allocators = std::map<std::string, Allocator*>();
template<> std::map<std::string, Allocator*> AllocatorMap<BoundaryMethod>::allocators = std::map<std::string, Allocator*>();
template<> std::map<std::string, Allocator*> AllocatorMap<EquationOfState>::allocators = std::map<std::string, Allocator*>();
template<> std::map<std::string, Allocator*> AllocatorMap<Solver>::allocators = std::map<std::string, Allocator*>();
template<> std::map<std::string, Allocator*> AllocatorMap<ExplicitMethod>::allocators = std::map<std::string, Allocator*>();
template<> std::map<std::string, Allocator*> AllocatorMap<FluxMethod>::allocators = std::map<std::string, Allocator*>();

template<typename Type>
class Create {
public:
	Type *operator()(int &i, char **argv) {
		std::string name = argv[++i];
		std::map<std::string, Allocator*>::iterator iter = AllocatorMap<Type>::allocators.find(name);
		if (iter == AllocatorMap<Type>::allocators.end()) throw Exception() << "failed to find " << name;
		return (Type*)(iter->second->create());
	}
};

class HydroApp : public GLApp {
	Hydro *hydro;
	HydroArgs args;

public:
	HydroApp()
	: hydro(NULL)
	{}

	virtual int main(int argc, char **argv) {

		AllocatorMap<InitialConditions>::allocators["Sod"] = new TypeAllocator<SodInitialConditions>();
		AllocatorMap<InitialConditions>::allocators["Sedov"] = new TypeAllocator<SedovInitialConditions>();
		AllocatorMap<InitialConditions>::allocators["Advect"] = new TypeAllocator<AdvectInitialConditions>();
		AllocatorMap<InitialConditions>::allocators["Wave"] = new TypeAllocator<WaveInitialConditions>();

		AllocatorMap<BoundaryMethod>::allocators["Periodic"] = new TypeAllocator<PeriodicBoundaryMethod>();
		AllocatorMap<BoundaryMethod>::allocators["Mirror"] = new TypeAllocator<MirrorBoundaryMethod>();
		AllocatorMap<BoundaryMethod>::allocators["Constant"] = new TypeAllocator<ConstantBoundaryMethod>();
		AllocatorMap<BoundaryMethod>::allocators["FreeFlow"] = new TypeAllocator<FreeFlowBoundaryMethod>();

		AllocatorMap<EquationOfState>::allocators["Euler"] = new TypeAllocator<EulerEquationOfState>();

		AllocatorMap<Solver>::allocators["EulerEquationBurgersSolverExplicit"] = new TypeAllocator<EulerEquationBurgersSolverExplicit>();
		AllocatorMap<Solver>::allocators["EulerEquationGodunovSolverExplicit"] = new TypeAllocator<EulerEquationGodunovSolverExplicit>();
		AllocatorMap<Solver>::allocators["EulerEquationRoeSolverExplicit"] = new TypeAllocator<EulerEquationRoeSolverExplicit>();

		AllocatorMap<ExplicitMethod>::allocators["ForwardEuler"] = new TypeAllocator<ForwardEulerExplicitMethod>();
		AllocatorMap<ExplicitMethod>::allocators["RK2"] = new TypeAllocator<RK2ExplicitMethod>();
		AllocatorMap<ExplicitMethod>::allocators["RK4"] = new TypeAllocator<RK4ExplicitMethod>();
		AllocatorMap<ExplicitMethod>::allocators["ICN3"] = new TypeAllocator<ICN3ExplicitMethod>();

		AllocatorMap<FluxMethod>::allocators["DonorCell"] = new TypeAllocator<DonorCellFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["LaxWendroff"] = new TypeAllocator<LaxWendroffFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["BeamWarming"] = new TypeAllocator<BeamWarmingFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Fromm"] = new TypeAllocator<FrommFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["CHARM"] = new TypeAllocator<CHARMFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["HCUS"] = new TypeAllocator<HCUSFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["HQUICK"] = new TypeAllocator<HQUICKFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Koren"] = new TypeAllocator<KorenFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["MinMod"] = new TypeAllocator<MinModFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Oshker"] = new TypeAllocator<OshkerFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Ospre"] = new TypeAllocator<OspreFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Smart"] = new TypeAllocator<SmartFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Sweby"] = new TypeAllocator<SwebyFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["UMIST"] = new TypeAllocator<UMISTFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["VanAlbada1"] = new TypeAllocator<VanAlbada1FluxMethod>();
		AllocatorMap<FluxMethod>::allocators["VanAlbada2"] = new TypeAllocator<VanAlbada2FluxMethod>();
		AllocatorMap<FluxMethod>::allocators["VanLeer"] = new TypeAllocator<VanLeerFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["MonotonizedCentral"] = new TypeAllocator<MonotonizedCentralFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["Superbee"] = new TypeAllocator<SuperbeeFluxMethod>();
		AllocatorMap<FluxMethod>::allocators["BarthJespersen"] = new TypeAllocator<BarthJespersenFluxMethod>();
	
		for (int i = 0; i < argc; ++i) {
			if (!strcmp(argv[i], "--help")) {
				std::cout << "args:" << std::endl;
				std::cout << "--initialConditions" << std::endl;
				std::for_each(AllocatorMap<InitialConditions>::allocators.begin(), AllocatorMap<InitialConditions>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--boundaryMethod" << std::endl;
				std::for_each(AllocatorMap<BoundaryMethod>::allocators.begin(), AllocatorMap<BoundaryMethod>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--equationOfState" << std::endl;
				std::for_each(AllocatorMap<EquationOfState>::allocators.begin(), AllocatorMap<EquationOfState>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--solver" << std::endl;
				std::for_each(AllocatorMap<Solver>::allocators.begin(), AllocatorMap<Solver>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--explicitMethod" << std::endl;
				std::for_each(AllocatorMap<ExplicitMethod>::allocators.begin(), AllocatorMap<ExplicitMethod>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--fluxMethod" << std::endl;
				std::for_each(AllocatorMap<FluxMethod>::allocators.begin(), AllocatorMap<FluxMethod>::allocators.end(), [&](std::map<std::string, Allocator*>::value_type &value) {
					std::cout << "\t" << value.first << std::endl;
				});
				std::cout << "--size" << std::endl;
				std::cout << "--useCFL" << std::endl;
				std::cout << "--cfl" << std::endl;
				std::cout << "--fixedDT" << std::endl;
				std::cout << "--gamma" << std::endl;
				exit(0);
			}
			if (i < argc-1) {
				if (!strcmp(argv[i], "--initialConditions")) {
					args.initialConditions = Create<InitialConditions>()(i, argv);
				} else if (!strcmp(argv[i], "--boundaryMethod")) {
					args.boundaryMethod = Create<BoundaryMethod>()(i, argv);
				} else if (!strcmp(argv[i], "--equationOfState")) {
					args.equationOfState = Create<EquationOfState>()(i, argv);
				} else if (!strcmp(argv[i], "--solver")) {
					args.solver = Create<Solver>()(i, argv);
				} else if (!strcmp(argv[i], "--explicitMethod")) {
					args.explicitMethod = Create<ExplicitMethod>()(i, argv);
				} else if (!strcmp(argv[i], "--fluxMethod")) {
					args.fluxMethod = Create<FluxMethod>()(i, argv);
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
				}
			}
		}

		if (!args.initialConditions) args.initialConditions = new SodInitialConditions();
		if (!args.boundaryMethod) args.boundaryMethod = new MirrorBoundaryMethod();
		if (!args.equationOfState) args.equationOfState = new EulerEquationOfState();
		if (!args.solver) args.solver = new EulerEquationRoeSolverExplicit();
		if (!args.explicitMethod) args.explicitMethod = new ForwardEulerExplicitMethod();
		if (!args.fluxMethod) args.fluxMethod = new SuperbeeFluxMethod();
		
		return GLApp::main(argc, argv);
	}

	virtual void init() {
		GLApp::init();
	
		hydro = new Hydro(args);
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

