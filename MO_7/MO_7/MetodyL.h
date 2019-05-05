#pragma once
#include <vector>

class MetodyL
{
	std::vector<std::vector<double>> matrix;
	std::vector<double> b;
	std::vector<double> x0, x1;
	void initVec();
	void matrixSet();
	void bSet();
	double omega, eps;
	int size;

public:
	MetodyL();
	~MetodyL();

	void jacobiego();
	void gaussa_seidl();
	void sor();
	bool checkSol(std::vector<double> m);
	void copyMatrix(std::vector<std::vector<double> > &mat);

};