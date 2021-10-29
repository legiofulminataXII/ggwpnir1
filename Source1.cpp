#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <random>
#include <cmath>
using namespace std;

int main()
{
	ifstream fin("input.txt");
	ofstream fout("output.txt");
	if (!fin.is_open()) // Если файл не был открыт
		cout << "File can/t be opened.\n";
	else
	{
		int N, T;
		double Lx, Ly, Lz;
		double m0 = 6.63e-26;
		double k = 1.38e-23;
		fin >> N;
		fin >> Lx;
		fin >> T;
		fin.close();
		double sigma = sqrt(k * T / m0);
		std::random_device rd{};
		std::mt19937 gen{ rd() };
		std::normal_distribution<> d{ 0, sigma };
		vector<double>x;
		vector<double>y;
		vector<double>z;
		vector<double>Vx;
		vector<double>Vy;
		vector<double>Vz;
		x.resize(N);
		y.resize(N);
		z.resize(N);
		Vx.resize(N);
		Vy.resize(N);
		Vz.resize(N);
		int ID = 0;
		int n = pow(N, 1.0 / 3.0);
		Ly = Lz = Lx;
		double Z = 0.0;
		for (int k = 0; k < n; k++) {
			double Y = 0.0;
			for (int j = 0; j < n; j++) {
				double X = 0.0;
				for (int i = 0; i < n; i++) {
					x[ID] = X;
					y[ID] = Y;
					z[ID] = Z;
					Vx[ID] = d(gen);
					Vy[ID] = d(gen);
					Vz[ID] = d(gen);
					fout << '{' << z[ID] << ", ";
					fout << y[ID] << ", ";
					fout << x[ID] << "}" << '\t';
					fout << '{' << Vz[ID] << ", ";
					fout << Vy[ID] << ", ";
					fout << Vx[ID] << "}" << '\n';
					ID++;
					X = X + Lx / (n - 1);
				}
				Y = Y + Ly / (n - 1);
			}
			Z = Z + Lz / (n - 1);
		}
		double sumx = 0;
		double sumy = 0;
		double sumz = 0;
		double sumx2 = 0;
		double sumy2 = 0;
		double sumz2 = 0;
		for (ID = 0; ID < N; ID++) {
			sumx += Vx[ID];
			sumy += Vy[ID];
			sumz += Vz[ID];
			sumx2 += Vx[ID] * Vx[ID];
			sumy2 += Vy[ID] * Vy[ID];
			sumz2 += Vz[ID] * Vz[ID];
		}
		fout << sumx / N << endl;
		fout << sumy / N << endl;
		fout << sumz / N << endl;
		fout << sumx2 / N << endl;
		fout << sumy2 / N << endl;
		fout << sumz2 / N << endl;
		fout.close();
	}
	return 0;
}