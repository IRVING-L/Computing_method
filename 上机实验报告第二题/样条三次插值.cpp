#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
using namespace std;
/*
该程序为C++版本的样条三次插值函数
1.需要预先定义插值数据的个数
2.数据定义好后，将109行的for循环的起点、终点、精度设置为自己的插值数据范围的正确数字
3.该程序正常设置参数并运行后，将会在D盘中保存一个Excel文件，其内容是插值点的数值
*/
//定义插值数据的个数
const int arr_size = 27;

void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n, int width);
inline float get_y(const float x0[], const float y0[], const float mat_h[], const float mat_M[], const float x);
int main()
{
	
	//自己按照理解写、
	//样条插值大块上分为两部分
	//第一：A*M=D，先利用插值数据点和边界条件三，算出带状矩阵A，矩阵D，然后利用高斯消元反求M
	//第二：算出矩阵M后，再回到书上公式，求出x对应的函数值


	//第一步，先获取插值数据点
	float x0[arr_size], y0[arr_size];
	cout << "请输入插值数据点" << endl;
	for (int i = 0; i < arr_size; ++i)
	{
		cin >> x0[i] >> y0[i];
	}
	cout << "读取插值数据点完成" << endl;
	//第二步，算出D
	float mat_D[arr_size], x0_temp[3], y0_temp[3];
	for (int i = 1; i < arr_size - 1; ++i)
	{
		for (int j = -1; j <= 1; ++j)
		{
			x0_temp[1 + j] = x0[i + j];
			y0_temp[1 + j] = y0[i + j];
		}
		for (int k = 1; k < 3; ++k)
		{
			for (int q = 3 - 1; q >= k; --q)
			{
				y0_temp[q] = (y0_temp[q] - y0_temp[q - 1]) / (x0_temp[q] - x0_temp[q - k]);
			}
		}
		mat_D[i] = 6 * y0_temp[2];
	}
	//mat_D的首尾元素是边界条件d0,dn决定的,因此要算出d0,dn的值

	float mat_h[arr_size] = {};
	//先计算矩阵h的值
	for (int i = 1; i < arr_size; ++i)
	{
		mat_h[i] = x0[i] - x0[i - 1];
	}
	//算d0，他是四阶差商
	float temp[arr_size];
	for (int i = 0; i < 4; ++i) temp[i] = y0[i];
	for (int k = 1; k < 4; ++k)
	{
		for (int i = 3; i >= k; --i)
		{
			temp[i] = (temp[i]-temp[i-1]) / (x0[i]-x0[i-k]);
		}
	}
	float d0 = -12*mat_h[1]*temp[3];
	//算dn，他是四阶差商
	for (int i = 0; i < 4; ++i) temp[i] = y0[arr_size - 4 + i];
	for (int k = 1; k < 4; ++k)
	{
		for (int i = 3; i >= k; --i)
		{
			temp[i] = (temp[i] - temp[i - 1]) / (x0[arr_size - 4 + i] - x0[arr_size - 4 + i - k]);//x0的角标要注意
		}
	}
	float dn = 12*mat_h[arr_size-1]*temp[3];
	mat_D[0] = d0;
	mat_D[arr_size - 1] = dn;
	
	//第三步，算矩阵A
	float mat_A[arr_size][arr_size] = {};
	float mat_p1[arr_size-1], mat_p2[arr_size-1];
	
	for (int i = 0; i < arr_size - 1; ++i)
	{
		mat_p1[i] = mat_h[i + 1] / (mat_h[i] + mat_h[i + 1]); //λi
		mat_p2[i] = 1 - mat_p1[i];//μi
	}
	mat_p1[0] = -2;
	mat_p2[arr_size - 2] = -2;
	for (int i = 0; i < arr_size; ++i)
	{
		mat_A[i][i] = 2;
		if (i - 1 >= 0) mat_A[i][i - 1] = mat_p2[i - 1];
		if (i + 1 < arr_size) mat_A[i][i + 1] = mat_p1[i];
	}//矩阵A计算完成

	
	//第四步，使用高斯消元求解矩阵M
	float mat_M[arr_size];
	gauss_solve(mat_A, mat_D, mat_M, arr_size, 1);
	cout << "写入文件中..." << endl;
	ofstream fout("D:\\data.csv", ofstream::out | ofstream::trunc);
	if (!fout.is_open()) cerr << "打开文件失败" << endl;
	for (float i = 0; i < 52; i += 0.1)
	{
		fout <<i<<","<< -get_y(x0, y0, mat_h, mat_M, i) << endl;
	}
	fout.close();
	cout << "插值计算完成" << endl;
	return 0;
}


void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n,int width)
{
	for (int k = 0; k < n - 1; k++)
	{
		if (arr_A[k][k] == 0) cerr << "主元为0，程序失败" << endl;
		for (int i = k + 1; /*i <=k+width &&*/ i<n; i++)
		{
			float temp = arr_A[i][k] / arr_A[k][k];
			for (int j = k; /*j < k + 2 * width + 1 &&*/ j < n; j++)
			{
				arr_A[i][j] = arr_A[i][j] - temp * arr_A[k][j];
			}
			arr_B[i] = arr_B[i] - temp * arr_B[k];
		}
	}
	arrx[n - 1] = arr_B[n - 1] / arr_A[n - 1][n - 1];
	for (int k = n - 1; k >= 0; k--)
	{
		float temp = arr_B[k];
		for (int j = k + 1; /*j <= k+width &&*/ j < n; j++)
		{
			temp = temp - (double)arr_A[k][j] * arrx[j];
		}
		arrx[k] = temp / arr_A[k][k];
	}
}

inline float get_y(const float x0[],const float y0[], const float mat_h[], const float mat_M[], const float x)
{
	int k = 0;
	for (int i = 1; i < arr_size; ++i)
	{
		if (x0[i - 1] <= x && x < x0[i])
		{
			k = i;
			break;
		}
	}
	if (k > 0)
	{
		float a1 = (x0[k] - x);
		float a2 = (x - x0[k - 1]);
		float h = mat_h[k];
		float y = ( mat_M[k - 1] * a1 * a1 * a1 / 6
			+ mat_M[k] * a2 * a2 * a2 / 6
			+ (y0[k - 1] - mat_M[k - 1] * h * h / 6) * a1
			+ (y0[k] - mat_M[k] * h * h / 6) * a2 ) / h;
		return y;
	}
	else
	{
		cerr << "输入的x值超出了能够插值的范围" << endl
			<< "你输入的x值：" << x << endl;
		return -1;
	}
}