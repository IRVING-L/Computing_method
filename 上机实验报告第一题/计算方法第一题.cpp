#include<iostream>
#include<string>
#include<fstream>
#include<time.h>

using namespace std;
//��ǰ�����ȡ�ľ�������Լ�����
const int arr_size = 52100;
const int width = 8;

//����������
void arr_exchange(float arrtemp[arr_size][width * 2 + 1], float arra[arr_size][arr_size], int arr_size, int width);
void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n, int width);
struct FileInfo {
	int id;			// �����ļ���ʾ
	int ver;		// �����ļ��汾��(��ѹ��Ϊ258��ѹ��Ϊ514��
	//long int id1;		// ���ñ�־
};

struct HeadInfo {
	int n;		  	// ������Ľ���
	int q;          // ��״�����ϴ���
	int p;          // ��״�����´���
};

int main()
{
	//�ṹ���ʵ����
	FileInfo f1;
	HeadInfo h1;

	//�������ļ�
	cout << "�������ļ���\t(��ʾ��data20191.dat)" << endl;
	string filename;
	cin >> filename;
	ifstream fin(filename, ifstream::binary);
	if (!fin.is_open()) cerr << "error" << endl;

	//ʹ��read����ֱ�Ӷ�ȡ��Ӧ���ȵ�����
	fin.read((char*)&f1, sizeof(f1));
	fin.read((char*)&h1, sizeof(h1));
	//��ӡ����Ľ�����������Ϣ
	cout << "������ļ��汾�ţ�" << f1.ver << endl
		<< "���������" << h1.n << endl
		<< "����Ĵ���" << h1.p << endl;
	//�ж�ȫ�ֱ�����������鳤�Ⱥʹ������ȡ�ľ�����Ϣ�Ƿ�һ��
	if (arr_size != h1.n || width != h1.q)
	{
		cerr << "���鳤�ȶ��岻��ȷ�����˳����������ݣ�" << endl;
		return -1;
	}
	cout << "�Ƿ��������Y 0r N��" << endl;
	char confirm;
	cin >> confirm;

	//�ڶ�����������
	clock_t t_beg = clock();//ʱ���¼
	float(*arrtemp)[width * 2 + 1] = new float[arr_size][width * 2 + 1]{};
	float(*arra)[arr_size] = new float[arr_size][arr_size]{};
	float* arrb = new float[arr_size];
	float* arrx = new float[arr_size];
	clock_t t_mid1 = clock();
	cout << "�ڶ������������ʱ�䣺" << t_mid1 - t_beg <<"ms"<< endl;
	if (confirm == 'Y' || confirm == 'y')
	{
		
		switch (f1.ver)//f1.ver�����ж��Ƿ�Ϊѹ������
		{
		case 258://��ѹ������
			fin.read((char*)arra, 4 * arr_size * arr_size);
			fin.read((char*)arrb, 4 * arr_size);
			cout << "��ѹ�������ȡ���" << "\t";
			break;
		case 514://ѹ������
			fin.read((char*)arrtemp, 4 * arr_size * (2 * width + 1));
			fin.read((char*)arrb, 4 * arr_size);
			arr_exchange(arrtemp, arra, arr_size, width);
			cout << "ѹ�������ȡ���" << "\t";
			break;
		}
		fin.close();
		clock_t t_mid2 = clock();
		cout << "��ȡ�����ʱ��Ϊ��" << t_mid2 - t_mid1 << "ms" << endl;
		cout << "���ڼ��㣬��ȴ�..." << "\t";
		//�����ȡ��ϣ����ø�˹�������м���
		gauss_solve(arra, arrb, arrx, arr_size, width);

		clock_t t_end = clock();
		cout << "��������ʱ��Ϊ��" << t_end - t_mid2 << "ms" << endl;
		cout << "������ϣ���ӡ���" << endl;
		for (int i = 0; i < 10; ++i)
		{
			cout << "��" << i + 1 << "���⣺" << arrx[i] << endl;
		}
		cout << "..." << endl;
		for (int i = 0; i < 10; ++i)
		{
			cout << "��" << arr_size - 9 + i << "���⣺" << arrx[arr_size - 10 + i] << endl;
		}
	}
	else cout << "���˳�����" << endl;
	delete[]arrtemp;
	delete[]arra;
	delete[]arrb;
	delete[]arrx;
	return 0;
}

//ѹ������ת��Ϊ��ѹ������
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

//��˹��Ԫ��⺯��
void gauss_solve(float arr_A[arr_size][arr_size], float arr_B[arr_size], float arrx[arr_size], int n, int width)
{
	for (int k = 0; k < n - 1; k++)
	{
		if (arr_A[k][k] == 0) cerr << "��ԪΪ0������ʧ��" << endl;
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
