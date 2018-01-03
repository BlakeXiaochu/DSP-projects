#include <iostream>
#include <complex>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

#define PI 3.141592653

//global varibles
string path1 = "D:\\application\\Courses and Assignments\\大三上\\数字信号处理\\第二次课程设计\\附件1\\";
string path2 = "D:\\application\\Courses and Assignments\\大三上\\数字信号处理\\第二次课程设计\\附件2\\data.wav";
string file_list[10] = {"data1081.wav", "data1107.wav" , "data1140.wav" , "data1219.wav" , "data1234.wav",
			"data1489.wav", "data1507.wav", "data1611.wav", "data1942.wav", "data1944.wav" };
int wav_size[10] = {1081, 1107, 1140, 1219, 1234, 1489, 1507, 1611, 1942, 1944};
int offset = 44;
double fs = 8000.0;

double row_freq[4] = {697.0, 770.0, 852.0, 941.0};
double col_freq[4] = {1209.0, 1336.0, 1477.0, 1633.0};

//************
//*    FFT   *
//************
void invert(double* xn, int n, complex<double>* inverted_xn) {
	int bit_num = int(round(log2(n)));
	int half_width = int(floor(bit_num / 2.0));

	for (int i = 0; i < n; i++) {

		int t = 0;
		int bit_pos1 = 1;
		int bit_pos2 = n >> 1;
		for (int j = 0; j < half_width; j++) {
			t = (t | ((i & bit_pos1) << (bit_num - j * 2 - 1)));
			t = (t | ((i & bit_pos2) >> (bit_num - j * 2 - 1)));
			bit_pos1 = bit_pos1 << 1;
			bit_pos2 = bit_pos2 >> 1;
		}

		if ((bit_num & 0x01) == 1) {
			t = t | (i & bit_pos1);
		}

		inverted_xn[i] = xn[t];
	}

}

void dit_butterfly(complex<double> &x1, complex<double> &x2, complex<double> &rotation_factor, complex<double> &y1, complex<double> &y2) {
	complex<double> t = x2 * rotation_factor;
	y2 = x1 - t;
	y1 = x1 + t;
}

void dit_fft(double* xn, int n, complex<double>* Xk) {
	int m = int(round(log2(n)));

	invert(xn, n, Xk);

	for (int i = 0; i < m; i++) {
		int N = 1 << (i + 1);
		complex<double> Wn(cos(-2 * PI / N), sin(-2 * PI / N));

		for (int j = 0; j < n / N; j++) {
			complex<double> Wr = 1.0;

			for (int k = 0; k < N / 2; k++) {
				dit_butterfly(Xk[j*N + k], Xk[j*N + k + N / 2], Wr, Xk[j*N + k], Xk[j*N + k + N / 2]);
				Wr *= Wn;
			}
		}
	}
}

void fft(double* xn, int n, complex<double>* Xk) {
	dit_fft(xn, n, Xk);
}

//************
//* Goertzel *
//************
//v: the hidden states of filter, length 2. v[0], v[1] represent vk[n - 1], vk[n - 2]
//every time the func is called, v is updated
//a1, b1: cos(omega_k), Wk
void Goertzel_filter(double &x, double *v, complex<double> &y, double &a1, complex<double> &b1) {
	double vn = 0.0;

	vn = x + 2.0 * a1 * v[0] - v[1];
	y = vn - b1 * v[0];

	//update state
	v[1] = v[0];
	v[0] = vn;

	return;
}

void Goertzel(int N, int k, double *x, complex<double> &Xk) {
	//compute filter coeff
	double a1 = cos(2.0 * PI * k / N);
	complex<double> b1(cos(-2.0 * PI * k / N), sin(-2.0 * PI * k / N));

	double v[2] = { 0.0, 0.0 };
	complex<double> y;

	for (int i = 0; i < N; i++) {
		Goertzel_filter(x[i], v, y, a1, b1);
	}

	double xN = 0.0;
	Goertzel_filter(xN, v, y, a1, b1);

	Xk = y;
}

//look up encoding table
void lookup(double row_H[4], double col_H[4], char &c) {
	double t1 = 0, t2 = 0;
	int index1 = -1, index2 = -1;

	for (int j = 0; j < 4; j++) {
		if (row_H[j] > t1) { t1 = row_H[j]; index1 = j; }
		if (col_H[j] > t2) { t2 = col_H[j]; index2 = j; }
	}

	switch (index1)
	{
	case 0:
	{
		c = (index2 == 0) ? '1' :
			((index2 == 1) ? '2' :
			((index2 == 2) ? '3' :
				((index2 == 3) ? 'A' : 'E')));
		break;
	}
	case 1:
	{
		c = (index2 == 0) ? '4' :
			((index2 == 1) ? '5' :
			((index2 == 2) ? '6' :
				((index2 == 3) ? 'B' : 'E')));
		break;
	}
	case 2:
	{
		c = (index2 == 0) ? '7' :
			((index2 == 1) ? '8' :
			((index2 == 2) ? '9' :
				((index2 == 3) ? 'C' : 'E')));
		break;
	}
	case 3:
	{
		c = (index2 == 0) ? '*' :
			((index2 == 1) ? '0' :
			((index2 == 2) ? '#' :
				((index2 == 3) ? 'D' : 'E')));
		break;
	}
	default:
		c = 'E';
		break;
	}
}

//***********
//* (1) FFT *
//***********
void FFT_anaysis(double* xn, int n, int file_id) {
	int N = 2 << (int)(ceil(log2(n)));
	double *ext_xn = new double[N];
	complex<double> *Xk = new complex<double>[N];
	for (int j = 0; j < n; j++) {
		ext_xn[j] = xn[j];
	}
	for (int j = n; j < N; j++) {
		ext_xn[j] = 0;
	}

	//compute freq index
	double delta_freq = (2 * PI) / N;
	int row_freq_index[4];
	int col_freq_index[4];
	for (int j = 0; j < 4; j++) {
		row_freq_index[j] = int(row_freq[j] / delta_freq);
		col_freq_index[j] = int(col_freq[j] / delta_freq);
	}

	fft(ext_xn, N, Xk);
	double row_H[4];
	double col_H[4];
	for (int j = 0; j < 4; j++) {
		//choose the max value of two side of index
		double H1, H2;
		H1 = abs(Xk[row_freq_index[j]]);
		H2 = abs(Xk[row_freq_index[j] + 1]);
		row_H[j] = (H1 > H2) ? H1 : H2;

		H1 = abs(Xk[col_freq_index[j]]);
		H2 = abs(Xk[col_freq_index[j] + 1]);
		col_H[j] = (H1 > H2) ? H1 : H2;
	}

	char c;
	lookup(row_H, col_H, c);
	cout << "FFT, " << file_list[file_id] << ": " << c;
	cout << endl;
}

//****************
//* (2) Goertzel *
//****************
void Goertzel_anaysis(double* xn, int n, int file_id) {
	int N = n;
	//compute freq index
	double delta_freq = (2 * PI) / N;
	int row_freq_index[4];
	int col_freq_index[4];
	for (int j = 0; j < 4; j++) {
		row_freq_index[j] = int(row_freq[j] / delta_freq);
		col_freq_index[j] = int(col_freq[j] / delta_freq);
	}

	double row_H[4];
	double col_H[4];
	//compute Amp response of specific index
	for (int j = 0; j < 4; j++) {
		//choose the max value of two side of index
		double H1, H2;
		complex<double> K1, K2;
		Goertzel(N, row_freq_index[j], xn, K1);
		H1 = abs(K1);
		Goertzel(N, row_freq_index[j] + 1, xn, K2);
		H2 = abs(K2);
		row_H[j] = (H1 > H2) ? H1 : H2;

		Goertzel(N, col_freq_index[j], xn, K1);
		H1 = abs(K1);
		Goertzel(N, col_freq_index[j] + 1, xn, K2);
		H2 = abs(K2);
		col_H[j] = (H1 > H2) ? H1 : H2;
	}

	char c;
	lookup(row_H, col_H, c);
	cout << "Goertzel, " << file_list[file_id] << ": " << c;
	cout << endl;
}

//****************
//* (3) long wav *
//****************
void signal_judge(double row_H[4], double col_H[4], char &c) {
	double sum = 0.0;
	double max = 0.0;

	for (int i = 0; i < 4; i++) {
		sum += row_H[i];
		max = (row_H[i] > max) ? row_H[i] : max;
	}
	sum -= max;
	sum /= 3.0;
	if (max / sum <= 15.0) {
		c = '_';
		return;
	}

	sum = 0.0;
	max = 0.0;
	for (int i = 0; i < 4; i++) {
		sum += col_H[i];
		max = (col_H[i] > max) ? col_H[i] : max;
	}
	sum -= max;
	sum /= 3.0;
	//行频率中的最大值需要大于其他三个频率点的平均值的5倍，才被认为是非噪声信号
	if (max / sum <= 5.0) {
		c = '_';
		return;
	}

	lookup(row_H, col_H, c);
	return;
}

void long_wav_anaysis(bool judge = true) {
	//short-time anaysis window length and step
	int win_len = 300;
	int step = 100;

	//hamming window function
	double *win = new double[win_len];
	for (int j = 0; j < win_len; j++) {
		win[j] = 0.54 - (0.46 * cos(2 * PI * j / win_len));
	}

	double delta_freq = (2 * PI) / win_len;
	int row_freq_index[4];
	int col_freq_index[4];
	for (int j = 0; j < 4; j++) {
		row_freq_index[j] = int(row_freq[j] / delta_freq);
		col_freq_index[j] = int(col_freq[j] / delta_freq);
	}

	string wav_path = path2;
	int N = 53082;
	short* data = new short[N];
	ifstream ifs(wav_path, ifstream::binary);
	if (!ifs)
		cout << "reading wav file data failed..." << endl;

	ifs.seekg(offset);
	ifs.read((char*)data, N * 2);
	ifs.close();

	// convert into double
	double *xn = new double[N];
	for (int j = 0; j < N; j++) {
		//normalized to [-1, 1]
		xn[j] = double(data[j] / 32767.0);
	}

	cout << endl << "long wav anaysis: " << endl;

	int start = 0;
	double *region = new double[win_len];
	while (true) {
		if ((start + win_len) >= N) break;

		//hamming
		for (int j = 0; j < win_len; j++) {
			region[j] = xn[start + j] * win[j];
		}

		double row_H[4];
		double col_H[4];
		//compute Amp response of specific index
		for (int j = 0; j < 4; j++) {
			//choose the max value of two side of index
			double H1, H2;
			complex<double> K1, K2;
			Goertzel(win_len, row_freq_index[j], region, K1);
			H1 = abs(K1);
			Goertzel(win_len, row_freq_index[j] + 1, region, K2);
			H2 = abs(K2);
			row_H[j] = (H1 > H2) ? H1 : H2;

			Goertzel(win_len, col_freq_index[j], region, K1);
			H1 = abs(K1);
			Goertzel(win_len, col_freq_index[j] + 1, region, K2);
			H2 = abs(K2);
			col_H[j] = (H1 > H2) ? H1 : H2;
		}

		char c;
		//判断是否为噪声
		if (judge) signal_judge(row_H, col_H, c);
		else lookup(row_H, col_H, c);
		cout << c;

		start += step;
	}
	cout << endl;
}

int main() {
	//please change the path before using
	string wav_path = path1;

	//freq normalization(to [0, 2*pi])
	for (int i = 0; i < 4; i++) {
		row_freq[i] /= fs;
		row_freq[i] *= (2 * PI);
		col_freq[i] /= fs;
		col_freq[i] *= (2 * PI);
	}

	for (int i = 0; i < 10; i++) {
		// read wav files
		int N = wav_size[i];
		short* data = new short[N];

		ifstream ifs(wav_path + file_list[i], ifstream::binary);
		if (!ifs)
			cout << "reading wav file data failed..." << endl;

		ifs.seekg(offset);
		ifs.read((char*)data, N * 2);
		ifs.close();
;
		// convert into double
		double *xn = new double[N];
		for (int j = 0; j < N; j++) {
			//normalized to [-1, 1]
			xn[j] = double(data[j] / 32767.0);
		}
		
		//***********
		//* (1) FFT *
		//***********
		FFT_anaysis(xn, N, i);


		//****************
		//* (2) Goertzel *
		//****************
		Goertzel_anaysis(xn, N, i);
	}

	//****************
	//* (3) long wav *
	//****************
	long_wav_anaysis();

	return 0;
}
