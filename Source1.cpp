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
	if (!fin.is_open()) // Esli file ne bil otkrit
		cout << "File can/t be opened.\n";
	else
	{
		int N;
		double Lx, Ly, Lz, T0, Tx, Ty, Tz, T, c, Txd, Tyd, Tzd, Td, Txm, Tym, Tzm, Tm;
		double m0 = 6.63e-26;
		double k = 1.38e-23;
		fin >> N;
		fin >> Lx;
		fin >> T0;
		fin.close();
		double sigma = sqrt(k * T0 / m0);
		fout << sigma*sigma << endl;
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
		int n = ceil(pow(N, 1.0 / 3.0));
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
		Tx = sumx2 * m0 / (N * k);
		Ty = sumy2 * m0 / (N * k);
		Tz = sumz2 * m0 / (N * k);
		T = (Tx + Ty + Tz) / 3; // samaia pervaia temperatura, bez dreifa i mashtabirovania
		fout << "Vx средняя равна " << sumx / N << endl;
		fout << "Vy средняя равна " << sumy / N << endl;
		fout << "Vz средняя равна " << sumz / N << endl;
		fout << "V^2 средняя равна " << ((sumx2 / N) + (sumy2 / N) + (sumz2 / N)) / 3 << endl;
		fout << "T0 была равна " << T0 << "K, T теперь равнa " << T << "K" << endl;
		vector<double>Vxd;
		vector<double>Vyd;
		vector<double>Vzd;
		Vxd.resize(N);
		Vyd.resize(N);
		Vzd.resize(N);
		for (ID = 0; ID < N; ID++) {
			Vxd[ID] = Vx[ID] - (sumx / N);
			Vyd[ID] = Vy[ID] - (sumy / N);
			Vzd[ID] = Vz[ID] - (sumz / N); // ustranili dreif
		}
		double sumxd = 0;
		double sumyd = 0;
		double sumzd = 0;
		double sumxd2 = 0;
		double sumyd2 = 0;
		double sumzd2 = 0;
		for (int ID = 0; ID < N; ID++) {
			sumxd += Vxd[ID];
			sumyd += Vyd[ID];
			sumzd += Vzd[ID];
			sumxd2 += Vxd[ID] * Vxd[ID];
			sumyd2 += Vyd[ID] * Vyd[ID];
			sumzd2 += Vzd[ID] * Vzd[ID];
		}
		Txd = sumxd2 * m0 / (N * k);
		Tyd = sumyd2 * m0 / (N * k);
		Tzd = sumzd2 * m0 / (N * k);
		Td = (Txd + Tyd + Tzd) / 3; // temperatura posle dreifa, no do mashtabirovania
		fout << "После устранения дрейфа проекции скоростей поменялись:" << endl;
		fout << "Vx средняя равна " << sumxd / N << endl;
		fout << "Vy средняя равна " << sumyd / N << endl;
		fout << "Vz средняя равна " << sumzd / N << endl;
		fout << "Температура после устранения дрейфа равна " << Td << "K." << endl;
		vector<double>Vxm;
		vector<double>Vym;
		vector<double>Vzm;
		Vxm.resize(N);
		Vym.resize(N);
		Vzm.resize(N);
		double sumxm = 0;
		double sumym = 0;
		double sumzm = 0;
		double sumxm2 = 0;
		double sumym2 = 0;
		double sumzm2 = 0;
		c = T0 / Td; // vveli coefficient
		for (int ID = 0; ID < N; ID++) {
			Vxm[ID] = Vxd[ID] * c;
			Vym[ID] = Vyd[ID] * c;
			Vzm[ID] = Vzd[ID] * c; // domnozhili vse skorosti
			sumxm += Vxm[ID];
			sumym += Vym[ID];
			sumzm += Vzm[ID];
			sumxm2 += Vxm[ID] * Vxm[ID];
			sumym2 += Vym[ID] * Vym[ID];
			sumzm2 += Vzm[ID] * Vzm[ID];
		}
		Txm = sumxm2 * m0 / (N * k);
		Tym = sumym2 * m0 / (N * k);
		Tzm = sumzm2 * m0 / (N * k);
		Tm = (Txm + Tym + Tzm) / 3; // temperatura posle vseh manipuliatsiy
		fout << "Температура после масштабирования равна " << Tm << "K." << endl;
		fout.close();
	}
	return 0;
}
