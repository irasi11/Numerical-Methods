#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

template <typename T> 
class DerivativeCalc {
public:
	T fun(T x) {
		return cos(x / 2);
	}
	T first_deriv_backward(T x, T begin, T end, T h) {
		if (h != 0 && x >= begin && x < end) return (fun(x) - fun(x - h)) / h;;
		if (h != 0 && x == begin) return (-3.0/2 * fun(x) + 2* fun(x+h) - 0.5*fun(x+h+h)) /h;
		if (h != 0 && x == end) return (0.5*fun(x-h-h)-2*fun(x+h) + 3.0/2 *fun(x))/h;
		else return -1;
	}

	T first_deriv_central(T x, T begin, T end, T h) {
		if(h!=0) return (fun(x + h) - fun(x - h)) / (2 * h);
		else return -1;
	}

	T first_deriv_forward(T x, T begin, T end, T h) {
		if (h != 0 && x >= begin && x < end) return (fun(x + h) - fun(x)) / h;
		if (h != 0 && x == begin) return (-3.0 / 2 * fun(x) + 2 * fun(x + h) - 0.5*fun(x + h + h)) / h;
		if (h != 0 && x == end) return (0.5*fun(x - h - h) - 2 * fun(x + h) + 3.0 / 2 * fun(x)) / h;
		else return -1;
	}

	T second_deriv(T x, T begin, T end, T h) {
		if (h != 0 && x >= begin && x < end) return (fun(x) - 2 * fun(x+h) + fun(x + h +h)) / (h*h);
		if (h != 0 && x == begin) return (fun(x) - 2 * fun(x+h) + fun(x + h + h)) / (h * h);
		if (h != 0 && x == end) return (fun(x - h -h) - 2 * fun(x - h) + fun(h)) / (h * h);
		else return -1;
	}
};