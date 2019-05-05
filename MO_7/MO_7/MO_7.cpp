#include "pch.h"
#include "MetodyL.h"
#include <iostream>

int main()
{

	std::shared_ptr<MetodyL> jac = std::make_shared<MetodyL>();
	jac->jacobiego();
	jac->gaussa_seidl();
	jac->sor();

	return 0;
}
