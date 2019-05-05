#include "pch.h"
#include "MetodyL.h"


void MetodyL::initVec()
{
	x0.resize(size);
	x1.resize(size);
	for (int i = 0; i < size; i++) {
		x0.at(i) = 1;
		x1.at(i) = 0;
	}
}

MetodyL::MetodyL()
{
	
	size =4; //define size of matrix
	eps = 1e-7; // precision of calculations
	omega = 1.0 / 2.0; // param 
	matrix.resize(4);
	matrixSet();
	bSet();	
}

void MetodyL::bSet() {
	b.push_back(395); b.push_back(603); b.push_back(-415); b.push_back(-606);
}

void MetodyL::matrixSet() {
	matrix[0].push_back(100); matrix[0].push_back(1); matrix[0].push_back(-2); matrix[0].push_back(3);
	matrix[1].push_back(4); matrix[1].push_back(300); matrix[1].push_back(-5); matrix[1].push_back(6);
	matrix[2].push_back(7); matrix[2].push_back(-8); matrix[2].push_back(400); matrix[2].push_back(9);
	matrix[3].push_back(-10); matrix[3].push_back(11); matrix[3].push_back(-12); matrix[3].push_back(200);
}

MetodyL::~MetodyL()
{
}

void MetodyL::jacobiego()
{
	initVec();
	std::vector<std::vector<double>> m;
	m.resize(size);
	copyMatrix(m);
	int iter = 0;
	std::cout << "----------METODA JACOBIEGO----------------\n";
	while ((iter < 100) && fabs(x0[0] - x1[0]) > eps && fabs(x0[1] - x1[1]) > eps && fabs(x0[2] - x1[2]) > eps && fabs(x0[3] - x1[3]) > eps) {
		for (int i = 0; i < size; i++) {
			x0.at(i) = x1.at(i);
		}
		for (int i = 0; i < size; i++) {
			x1[i]=0;
			for (int j = 0; j < size; j++) {
				if (i != j) {
					x1[i] += (-1.0) *  (m[i][j] / m.at(i).at(i)) * x0.at(j);
				}
			}
			x1.at(i) += b.at(i)/m.at(i).at(i);

		}
				
		iter++;
		printf("%3d   %5.5lf	 %5.5lf	 %5.5lf   %5.5lf \n ", iter, x0[0], x0[1], x0[2], x0[3]);
		if (checkSol(x0)) break;
	}

	
}

void MetodyL::gaussa_seidl()
{
	//initiates vector by default values
	initVec();
	std::cout << "----------------Gaussa-Seidl------------------\n";
	std::vector<std::vector<double>> m;
	m.resize(size);
	copyMatrix(m);
	int iter = 0;
	while ((iter < 100) && fabs(x0[0] - x1[0]) > eps && fabs(x0[1] - x1[1]) > eps && fabs(x0[2] - x1[2]) > eps && fabs(x0[3] - x1[3]) > eps) {
		for (int i = 0; i < size; i++) {
			x0.at(i) = x1.at(i);
			x1[i] = 0;
			for (int j = 0; j < size; j++) {
				if (i != j) {
					x1[i] += (-1.0) *  (m[i][j] / m.at(i).at(i)) * x0.at(j);
				}
			}
			x1.at(i) += b.at(i) / m.at(i).at(i);
			

		}

		iter++;
		printf("%3d   %5.5lf	 %5.5lf	 %5.5lf   %5.5lf \n ", iter, x0[0], x0[1], x0[2], x0[3]);
		if (checkSol(x0)) break;
	}
	
}

//solving equalations by sor method
void MetodyL::sor()
{
	initVec();
	std::cout << "----------------SOR------------------\n";
	std::vector<std::vector<double>> m;
	m.resize(size);
	copyMatrix(m);
	double s1 = 0, s2 = 0;
	int iter = 0;
	while ((iter < 100) && fabs(x0[0] - x1[0]) > eps && fabs(x0[1] - x1[1]) > eps && fabs(x0[2] - x1[2]) > eps && fabs(x0[3] - x1[3]) > eps) {
		for (int i = 0; i < size; i++) {
			x0.at(i) = x1.at(i);
			x1[i] = 0; s1 = 0; s2 = 0;
			for (int j = 0; j < size; j++) {
				if (i < j) {
					s1 += m[i][j] * x0[j];
				}
				if (i >= j) {
					s2 += m[i][j] * x0[j];
				}
			}
			x1.at(i) += x0[i] - ((omega / m[i][i])*(s1 + s2 - b[i]));


		}
		
		iter++;
		printf("%3d   %5.5lf	 %5.5lf	 %5.5lf   %5.5lf \n ", iter, x0[0], x0[1], x0[2], x0[3]);
		if (checkSol(x0)) break;
	}
}

//test third criterion of convergence
bool MetodyL::checkSol(std::vector<double> m)
{
	double k = 0;
	std::vector<double> v;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			v.push_back(matrix[i][j] * m[j]);
		}
		if (fabs(v[i] - b[i]) < eps)
			k++;
	}
	if (k == 4)
		return true;
	else
		return false;
}

//copy matrix to matrix given in param
void MetodyL::copyMatrix(std::vector<std::vector<double>> &mat)
{

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			mat[i].push_back(matrix[i][j]);
		}
	}
}
