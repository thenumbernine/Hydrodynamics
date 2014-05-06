#pragma once

//cell interface values
struct Interface {
	//static values
	double x;	//position

	//dynamic values

	//used by Burgers
	std::vector<double> r;
	std::vector<double> flux;
	double velocity;

	//used for Godunov
	std::vector<double> stateMid;
	std::vector<std::vector<double> > jacobian;
	std::vector<double> eigenvalues;
	std::vector<std::vector<double> > eigenvectors;
	std::vector<std::vector<double> > eigenvectorsInverse;
	std::vector<double> rTilde;	//r projected into the eigenvector basis
	std::vector<double> deltaQTilde;	//dq projected into eigenvector basis
};


