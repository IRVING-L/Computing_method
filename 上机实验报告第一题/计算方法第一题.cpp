#include<iostream>
#include<string>
#include<fstream>
#include<time.h>

using namespace std;
//提前定义读取的矩阵阶数以及带宽
/***********************
每一个文件里的带状矩阵的阶数和带宽都不一样：
data20191.dat:阶数 10， 带宽 3
data20192.dat:阶数 20， 带宽 4
data20193.dat:阶数 3120 带宽 6
data20194.dat:阶数 52100带宽 8
************************/
const int arr_size = 52100;
const int width = 8;

//函数的声明
void arr_exchange(float arrtemp[arr_size][width * 2 + 1], float arra[arr_size][arr_size], int arr_size, int width);
void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n, int width);
struct FileInfo {
	int id;			// 数据文件标示
	int ver;		// 数据文件版本号(非压缩为258，压缩为514）
	//long int id1;		// 备用标志
};

struct HeadInfo {
	int n;		  	// 方程组的阶数
	int q;          // 带状矩阵上带宽
	int p;          // 带状矩阵下带宽
};

int main()
{
	//结构体的实例化
	FileInfo f1;
	HeadInfo h1;

	//打开数据文件
	cout << "请输入文件名\t(提示：data20191.dat)" << endl;
	string filename;
	cin >> filename;
	ifstream fin(filename, ifstream::binary);
	if (!fin.is_open()) cerr << "error" << endl;

	//使用read函数直接读取相应长度的数据
	fin.read((char*)&f1, sizeof(f1));
	fin.read((char*)&h1, sizeof(h1));
	//打印矩阵的阶数、带宽信息
	cout << "矩阵的文件版本号：" << f1.ver << endl
		<< "矩阵阶数：" << h1.n << endl
		<< "矩阵的带宽：" << h1.p << endl;
	//判断全局变量定义的数组长度和带宽与读取的矩阵信息是否一致
	if (arr_size != h1.n || width != h1.q)
	{
		cerr << "数组长度定义不正确，请退出程序检查数据！" << endl;
		return -1;
	}
	cout << "是否继续？（Y 0r N）" << endl;
	char confirm;
	cin >> confirm;

	//在堆区创建矩阵
	clock_t t_beg = clock();//时间记录
	float(*arrtemp)[width * 2 + 1] = new float[arr_size][width * 2 + 1]{};
	float(*arra)[arr_size] = new float[arr_size][arr_size]{};
	float* arrb = new float[arr_size];
	float* arrx = new float[arr_size];
	clock_t t_mid1 = clock();
	cout << "在堆区创建数组的时间：" << t_mid1 - t_beg <<"ms"<< endl;
	if (confirm == 'Y' || confirm == 'y')
	{
		
		switch (f1.ver)//f1.ver用于判断是否为压缩矩阵
		{
		case 258://非压缩矩阵
			fin.read((char*)arra, 4 * arr_size * arr_size);
			fin.read((char*)arrb, 4 * arr_size);
			cout << "非压缩矩阵读取完成" << "\t";
			break;
		case 514://压缩矩阵
			fin.read((char*)arrtemp, 4 * arr_size * (2 * width + 1));
			fin.read((char*)arrb, 4 * arr_size);
			arr_exchange(arrtemp, arra, arr_size, width);
			cout << "压缩矩阵读取完成" << "\t";
			break;
		}
		fin.close();
		clock_t t_mid2 = clock();
		cout << "读取矩阵的时间为：" << t_mid2 - t_mid1 << "ms" << endl;
		cout << "正在计算，请等待..." << "\t";
		//矩阵读取完毕，调用高斯函数进行计算
		gauss_solve(arra, arrb, arrx, arr_size, width);

		clock_t t_end = clock();
		cout << "计算矩阵的时间为：" << t_end - t_mid2 << "ms" << endl;
		cout << "计算完毕，打印结果" << endl;
		for (int i = 0; i < 10; ++i)
		{
			cout << "第" << i + 1 << "个解：" << arrx[i] << endl;
		}
		cout << "..." << endl;
		for (int i = 0; i < 10; ++i)
		{
			cout << "第" << arr_size - 9 + i << "个解：" << arrx[arr_size - 10 + i] << endl;
		}
	}
	else cout << "已退出程序" << endl;
	delete[]arrtemp;
	delete[]arra;
	delete[]arrb;
	delete[]arrx;
	return 0;
}

//压缩矩阵转换为非压缩矩阵
void arr_exchange(float arrtemp[arr_size][width * 2 + 1], float arra[arr_size][arr_size], int arr_size, int width)
{
	for (int i = 0; i < arr_size; ++i)
	{
		arra[i][i] = arrtemp[i][width];
		for (int j = 1; j <= width; ++j)
		{
			if (arrtemp[i][width + j] != 0)
			{
				arra[i][i + j] = arrtemp[i][width + j];
			}
			if (arrtemp[i][width - j] != 0)
			{
				arra[i][i - j] = arrtemp[i][width - j];
			}
		}
	}
}

//高斯消元求解函数
void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n, int width)
{
	for (int k = 0; k < n - 1; k++)
	{
		if (arr_A[k][k] == 0) cerr << "主元为0，程序失败" << endl;
		for (int i = k + 1; i <= k + width && i < n; i++)
		{
			float temp = arr_A[i][k] / arr_A[k][k];
			for (int j = k; j <= k + 2 * width + 1 && j < n; j++)
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
		for (int j = k + 1; j < n && j <= k + width; j++)
		{
			temp = temp - (double)arr_A[k][j] * arrx[j];
		}
		arrx[k] = temp / arr_A[k][k];
	}
}
