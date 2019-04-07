#pragma once
#include <cstring>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
#include "AudioFile.h"
using namespace std;

class DTCWT_BANDS
{

public:
	DTCWT_BANDS();
	~DTCWT_BANDS();
	vector<vector<complex<double>>> BAND;
	vector<double> Lo;

};