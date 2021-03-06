// MO8.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include "pch.h"
#include "derivative_calc.h"
#include <iostream>
#include <fstream>
int main()
{
	const long double exact_value_central = 0.35355339059327376220042218105242451964241796884423701829;
	const long double exact_value_left_bound = 0;
	const long double exact_value_right_bound = 0.5;
	const long double exact_s_c = 0.17677669529663688110021109052621225982120898442211850914;
	const long double exact_s_r = 0;
	const long double exact_s_l = 0.25;
	double leftLimit = 0;
	double rightLimit = M_PI;
	std::shared_ptr<DerivativeCalc<float>> floatCalc = std::make_shared<DerivativeCalc<float>>();
	std::shared_ptr<DerivativeCalc<double>> doubleCalc = std::make_shared<DerivativeCalc<double>>();
	//file to save data
	FILE *file;
	#pragma warning(disable:4996)
	file = fopen("wykres.csv", "w");
	fprintf(file, "h;f_backward;f_central;f_forward;d_backward;d_central;d_forward;f_s;d_s\n");
	for (double h = 1e-10; h < 1 - 2*2e-4; h += 2e-4) {
		fprintf(file, "%1.5lf ;", h);
		fprintf(file, "%.15lf; ", fabs(fabs(floatCalc->first_deriv_backward(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(floatCalc->first_deriv_central(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(floatCalc->first_deriv_forward(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(doubleCalc->first_deriv_backward(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(doubleCalc->first_deriv_central(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(doubleCalc->first_deriv_forward(M_PI / 2, 0, M_PI, h)) - exact_value_central));
		fprintf(file, "%.15lf; ", fabs(fabs(floatCalc->second_deriv(M_PI / 2, 0, M_PI, h)) - exact_s_c));
		fprintf(file, "%.15lf; ", fabs(fabs(doubleCalc->second_deriv(M_PI / 2, 0, M_PI, h)) - exact_s_c));

		/*
		file << fabs(fabs(floatCalc->first_deriv_backward(M_PI / 2, 0, M_PI, h)) - exact_value_central) << ";";
		file << fabs(fabs(floatCalc->first_deriv_central(M_PI / 2, 0, M_PI, h)) - exact_value_central) << ";";
		file << fabs(fabs(floatCalc->first_deriv_forward(M_PI / 2, 0, M_PI, h)) - exact_value_central)<< ";";
		file << fabs(fabs(doubleCalc->first_deriv_backward(M_PI / 2, 0, M_PI, h)) - exact_value_central) << ";";
		file << fabs(fabs(doubleCalc->first_deriv_central(M_PI / 2, 0, M_PI, h)) - exact_value_central) << ";";
		file << fabs(fabs(doubleCalc->first_deriv_forward(M_PI / 2, 0, M_PI, h)) - exact_value_central) << ";";
		file << fabs(doubleCalc->second_deriv(M_PI / 2, 0, M_PI, h)) << ";";
		*/
		fprintf(file, "\n");
	}
	fclose(file);
	
	return 0;
}
