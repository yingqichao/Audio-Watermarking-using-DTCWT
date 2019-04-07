#pragma once

#include <cstring>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
using namespace std;

//简易Map.Entity
//Author:Qichao Ying

class entity
{

public:
	entity();
	entity(int ind, double val);
	~entity();
	int index;
	double value;


	/*** overloaded > operator  重载大于号***/
	bool operator > (const entity& d) {
		return value > d.value;
	}
	/*** overloaded < operator 重载小于号 ***/
	bool operator <(const entity& d) {
		return value < d.value;
	}




};