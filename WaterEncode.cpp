#include "WaterEncode.h"
#include <stdio.h>
#include "pch.h"

using namespace std;

WaterEncode::WaterEncode()
{
	sg0o = 19;
	sg1o = 13;
	sh0o = 13;
	sh1o = 19;

}

WaterEncode::~WaterEncode()
{
	delete g0o;
	delete g1o;
	delete h0o;
	delete h1o;
	delete g0a;
	delete g0b;
	delete g1a;
	delete g1b;
	delete h0a;
	delete h0b;
	delete h1a;
	delete h1b;
}

WaterEncode::WaterEncode(int len, int w_width, double alpha, int h, int w)
{
	m_len = len;

	sg0o = 19;
	sg1o = 13;
	sh0o = 13;
	sh1o = 19;

}

void WaterEncode::coldfilt(vector<double>& X, int H, string filt1, string filt2,
						   vector<double>& out)
{
	int W = 1;
	int filtersize = 14;
	int filtersize2 = filtersize / 2;
	int xe_len = H + 2 * filtersize;
	int *xe = new int[xe_len];

	double* filter1 = h0b;
	if (filt1 == "h1b")
		filter1 = h1b;
	double* filter2 = h0a;
	if (filt2 == "h1a")
		filter2 = h1a;


	for (int i = 0; i < H; i++)
		xe[filtersize + i] = i + 1;

	for (int j = 0; j < filtersize; j++)		//小于min的序号 负数
		xe[filtersize - 1 - j] = xe[filtersize + j];

	for (int j = 0; j < filtersize; j++)		//超过max的序号
		xe[filtersize + H + j] = xe[filtersize + H - 1 - j];

	double *filter1o = new double[filtersize2];
	double *filter1e = new double[filtersize2];
	double *filter2o = new double[filtersize2];
	double *filter2e = new double[filtersize2];

	for (int i = 0; i < filtersize2; i++)
	{
		filter1o[i] = filter1[i * 2];
		filter1e[i] = filter1[i * 2 + 1];
		filter2o[i] = filter2[i * 2];
		filter2e[i] = filter2[i * 2 + 1];
	}

	int length = (H + 2 * filtersize) / 4 - 1;
	int *t = new int[length];

	for (int i = 0; i < length; i++)
	{
		t[i] = i * 4 + 6;
	}

	double* conv1 = new double[length * W];
	double* conv2 = new double[length * W];
	double* conv3 = new double[length * W];
	double* conv4 = new double[length * W];

	for (int i = 0; i < length; i++)
	{
		int index1 = t[i] - 1 - 1;
		int index2 = t[i] - 3 - 1;
		int index3 = t[i] - 1;
		int index4 = t[i] - 2 - 1;
		index1 = xe[index1] - 1;
		index2 = xe[index2] - 1;
		index3 = xe[index3] - 1;
		index4 = xe[index4] - 1;
		for (int j = 0; j < W; j++)
		{
			conv1[i*W + j] = X[index1*W + j];
			conv2[i*W + j] = X[index2*W + j];
			conv3[i*W + j] = X[index3*W + j];
			conv4[i*W + j] = X[index4*W + j];
		}
	}

	double sum = 0;
	for (int i = 0; i < filtersize; i++)
		sum += filter1[i] * filter2[i];

	double con1, con2, con3, con4;
	int Hout = length - filtersize2 + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < W; j++)
		{
			con1 = 0;
			con2 = 0;
			con3 = 0;
			con4 = 0;
			for (int m = 0; m < filtersize2; m++)
			{
				con1 += conv1[(i + m) * W + j] * filter1o[filtersize2 - 1 - m];
				con2 += conv2[(i + m) * W + j] * filter1e[filtersize2 - 1 - m];
				con3 += conv3[(i + m) * W + j] * filter2o[filtersize2 - 1 - m];
				con4 += conv4[(i + m) * W + j] * filter2e[filtersize2 - 1 - m];
			}
			if (sum > 0)
			{
				out[j * Hout * 2 + 2 * i] = con1 + con2;		//转置!!!!!这里判断
				out[j * Hout * 2 + 2 * i + 1] = con3 + con4;
			}else{
				out[j * Hout * 2 + 2 * i + 1] = con1 + con2;	//转置!!!!!这里判断
				out[j * Hout * 2 + 2 * i] = con3 + con4;
			}
		}
	}

	delete filter1o;
	delete filter1e;
	delete filter2o;
	delete filter2e;

	delete conv1;
	delete conv2;
	delete conv3;
	delete conv4;

	delete xe;
	delete t;
}

void WaterEncode::colfilter(vector<double>& X, int index, int XHeight, string filt, int filtersize,
							vector<double>& out, bool addflag)
{
	int XWidth = 1;
	int size2 = filtersize / 2;
	int length = XHeight + 2 * size2;
	double* filter = h0o;
	if (filt == "h1o")
		filter = h1o;
	else if(filt=="g1o")
		filter = g1o;
	else if (filt == "g0o")
		filter = g0o;
	//Width always equals to 1 for 1-dimen
	int xe_len = XHeight + 2 * size2;//XHeight是长度，XWeight = 1;
	int* xe = new int[xe_len];

	for (int i = 0; i < XHeight; i++)
		xe[size2 + i] = i + 1;

	for (int j = 0; j < size2; j++)		//小于min的序号 负数
		xe[size2 - 1 - j] = xe[size2 + j];

	for (int j = 0; j < size2; j++)		//超过max的序号
		xe[size2 + XHeight + j] = xe[size2 + XHeight - 1 - j];

	double* conv = new double[length * XWidth];
	for (int i = 0; i < length; i++)
	{
		int ind = xe[i] - 1;
		for (int j = 0; j < XWidth; j++)
			conv[i * XWidth + j] = X[index+ind * XWidth + j];
	}

	int Hout = xe_len - filtersize + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < XWidth; j++)
		{
			double con = 0;
			for (int m = 0; m < filtersize; m++)
				con += conv[(i + m) * XWidth + j] * filter[filtersize - 1 - m];

			if (addflag)						//转置
				out[j * Hout + i] += con;
			else {
				out[j * Hout + i] = con;
				//cout << j * Hout + i << " : " << con << endl;
			}
		}
	}

	delete xe;
	delete conv;
}

void WaterEncode::colifilt(vector<double>& X, int H,  string filt1, string filt2,
						   vector<double>& out){
	int W = 1;
	int filtersize = 14;
	int filtersize2 = filtersize / 2;
	int xe_len = H + 2 * filtersize2;

	double* filter1 = g0b;
	if (filt1 == "g1b")
		filter1 = g1b;
	double* filter2 = g0a;
	if (filt2 == "g1a")
		filter2 = g1a;

	int* xe = new int[xe_len];

	double* filter1o = new double[filtersize2];
	double* filter1e = new double[filtersize2];
	double* filter2o = new double[filtersize2];
	double* filter2e = new double[filtersize2];

	int length = (H + filtersize) / 2 - 1;
	int* t = new int[length];
	double* conv1 = new double[length * W];
	double* conv2 = new double[length * W];

	for (int i = 0; i < H; i++)
		xe[filtersize2 + i] = i + 1;
	for (int j = 0; j < filtersize2; j++)		//小于min的序号
		xe[filtersize2 - 1 - j] = xe[filtersize2 + j];
	for (int j = 0; j < filtersize2; j++)		//超过max的序号
		xe[filtersize2 + H + j] = xe[filtersize2 + H - 1 - j];

	for (int i = 0; i < filtersize2; i++)
	{
		filter1o[i] = filter1[i * 2];
		filter1e[i] = filter1[i * 2 + 1];
		filter2o[i] = filter2[i * 2];
		filter2e[i] = filter2[i * 2 + 1];
	}

	for (int i = 0; i < length; i++)
		t[i] = i * 2 + 3;

	double sum = 0;
	for (int i = 0; i < filtersize; i++)
		sum += filter1[i] * filter2[i];


	int index1, index2;
	for (int i = 0; i < length; i++)
	{
		if (sum > 0)	//ta tb
		{
			index1 = t[i] - 1;
			index2 = t[i] - 2;
		}
		else
		{
			index1 = t[i] - 2;
			index2 = t[i] - 1;
		}
		index1 = xe[index1] - 1;
		index2 = xe[index2] - 1;

		for (int j = 0; j < W; j++)
		{
			conv1[i * W + j] = X[index2 * W + j];
			conv2[i * W + j] = X[index1 * W + j];
		}
	}

	double con1, con2, con3, con4;
	int Hout = length - filtersize2 + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < W; j++)
		{
			con1 = 0;
			con2 = 0;
			con3 = 0;
			con4 = 0;
			for (int m = 0; m < filtersize2; m++)
			{
				con1 += conv1[(i + m) * W + j] * filter1o[filtersize2 - 1 - m];
				con2 += conv2[(i + m) * W + j] * filter2o[filtersize2 - 1 - m];
				con3 += conv1[(i + m) * W + j] * filter1e[filtersize2 - 1 - m];
				con4 += conv2[(i + m) * W + j] * filter2e[filtersize2 - 1 - m];
			}
			out[j * Hout * 4 + 4 * i] += con1;		//转置
			out[j * Hout * 4 + 4 * i + 1] += con2;
			out[j * Hout * 4 + 4 * i + 2] += con3;
			out[j * Hout * 4 + 4 * i + 3] += con4;
		}
	}

	delete conv1;
	delete conv2;

	delete filter1o;
	delete filter1e;
	delete filter2o;
	delete filter2e;

	delete xe;
	delete t;
}


DTCWT_BANDS* WaterEncode::dtwavexfm_1d(vector<double>& X, int index, int len, int level)
{
	DTCWT_BANDS* res = new DTCWT_BANDS();
	vector<vector<complex<double>>> bands;
	vector<complex<double>> Yh1;
	vector<double> Hi; vector<double> Lo;
	Hi.resize(len); Lo.resize(len);

	//Level 1
	
	colfilter(X, index, len, "h0o", sh0o, Lo, false);
	//test
	//for (int i = 0; i < len; i++) {
	//	cout << Lo[0][i] << endl;
	//}
	colfilter(X, index, len, "h1o", sh1o, Hi, false);

	for (int i = 0; i < len/2; i ++) {
		complex<double> tmp(Hi[2 * i],Hi[2 * i +1]);
		Yh1.push_back(tmp);
	}


	//Level 2
	vector<complex<double>> Yh2;
	vector<double> Hi1; vector<double> Lo1;
	Hi1.resize(len/2); Lo1.resize(len/2);
	coldfilt(Lo, len, "h1b", "h1a", Hi1);
	//test
	//for (int i = 0; i < len / 2; i++) {
	//	cout << Hi1[0][i] << endl;
	//}


	coldfilt(Lo, len, "h0b", "h0a", Lo1);

	for (int i = 0; i < len / 4; i++) {
		complex<double> tmp(Hi1[2*i], Hi1[2*i + 1]);
		Yh2.push_back(tmp);
	}

	bands.push_back(Yh1); bands.push_back(Yh2);

	res->Lo = Lo1; res->BAND = bands;

	//cout << sizeof(complex<double>) << endl;
	//cout << res->BAND[1].size() << endl;
	//cout << sizeof(res->BAND[1].size()) / sizeof(complex<double>) << endl;
	
	//cout << "Finished DTCWT at Frame " << (int)(index / len) << endl;

	return res;

}


void WaterEncode::dtwaveifm_1d(DTCWT_BANDS* bands, int H, int level,
								vector<double>& out, int index)
{
		int W = 1;
		//Level 2
		int len = bands->Lo.size();
		vector<double> Lo1; vector<double> Lo2;
		Lo1.resize(H); Lo2.resize(H);

		vector<double> Hi;
		int len1 = bands->BAND[1].size();//
		/*cout << sizeof(complex<double>) << endl;
		cout << sizeof(bands->BAND[1].size()) << endl;
		cout << sizeof(bands->BAND[1].size()) / sizeof(complex<double>) << endl;*/
		Hi.resize(len1*2);
			
		c2q1d(bands->BAND[1],Hi);

		colifilt(bands->Lo, len, "g0b", "g0a", Lo1);
		colifilt(Hi,len,"g1b","g1a",Lo2);

		for(int i = 0; i < H; i++)
			Lo1[i] += Lo2[i];

		//Level 1
		int len2 = bands->BAND[0].size();
		vector<double> Hi1;
		Hi1.resize(len2*2);

		c2q1d(bands->BAND[0],Hi1);

		vector<double> Lo3; vector<double> Lo4;
		Lo3.resize(H); Lo4.resize(H);

		colfilter(Lo1, 0, H, "g0o", sg0o, Lo3, false);
		colfilter(Hi1, 0, H, "g1o", sg1o, Lo4, false);

		for (int i = 0; i < len * 2; i++)
			out[i+index] = (Lo3[i]+ Lo4[i]);

}

void WaterEncode::c2q1d(vector<complex<double>>& comp, vector<double>& res) {
	int len = comp.size();//-1 -> '\0'
	/*AudioFile<double>::AudioBuffer res; 
	res.resize(1);
	res[0].resize(len*2); */
	for (int i = 0; i < len; i++) {
		res[2 * i] = comp[i].real();
		res[2 * i+1] = comp[i].imag();
	}
	//return res;
}
