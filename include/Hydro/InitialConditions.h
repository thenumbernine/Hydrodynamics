#pragma once

class Hydro;

class InitialConditions {
public:
	virtual void operator()(Hydro *hydro) = 0;
};

class SodInitialConditions : public InitialConditions {
public:
	virtual void operator()(Hydro *hydro); 
};

class SedovInitialConditions : public InitialConditions {
public:
	virtual void operator()(Hydro *hydro); 
};

class AdvectInitialConditions : public InitialConditions {
public:
	virtual void operator()(Hydro *hydro); 
};

class WaveInitialConditions : public InitialConditions {
public:
	virtual void operator()(Hydro *hydro); 
};
