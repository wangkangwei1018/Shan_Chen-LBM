#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
using namespace std;

const int imax = 200;
const int jmax = 200;
const int Q = 9;
const double w[Q] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
const int cx[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int cy[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
const int opp[] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

double fin[imax + 1][jmax + 1][Q]{}, fout[imax + 1][jmax + 1][Q]{}, P[imax + 1][jmax + 1]{};
double rho[imax + 1][jmax + 1]{}, psi[imax + 1][jmax + 1]{}, U[imax + 1][jmax + 1]{}, V[imax + 1][jmax + 1]{};

double Feq(double rho, double U, double V, int k);
void init();
void Evolution();
void Output_plt(int iter);

double Gc = -6.5;
double tau = 1.0;
double rhoavg = 0.693;
int iter = 0;

int main()
{
	init();

	for (; iter <= 10000; iter++)
	{
		Evolution();

		if (iter % 1000 == 0)
			Output_plt(iter);
	}

	return 0;
}

double Feq(double rho, double U, double V, int k)
{
	double t1, t2;

	t1 = cx[k] * U + cy[k] * V;
	t2 = U * U + V * V;

	return w[k] * rho * (1 + 3 * t1 + 4.5 * t1 * t1 - 1.5 * t2);
}

void init()
{
	for(int i = 0; i <= imax; i++)
		for (int j = 0; j <= jmax; j++)
		{
			rho[i][j] = rhoavg - 0.5 * 0.01 * rhoavg + 0.01 * rhoavg * rand() / RAND_MAX;

			for (int k = 0; k < Q; k++)
			{
				fin[i][j][k] = Feq(rho[i][j], U[i][j], V[i][j], k);
			}
		}
}

void Evolution()
{
	//Collision
	for(int i = 0; i <= imax; i++)
		for (int j = 0; j <= jmax; j++)
		{
			for (int k = 0; k < Q; k++)
			{
				fout[i][j][k] = fin[i][j][k] - (1. / tau) * (fin[i][j][k] - Feq(rho[i][j], U[i][j], V[i][j], k));
			}
		}

	//Streaming
	for(int i = 0; i <= imax; i++)
		for (int j = 0; j <= jmax; j++)
		{
			for (int k = 0; k < Q; k++)
			{
				int ip = (i - cx[k] + imax + 1) % (imax + 1);
				int jp = (j - cy[k] + jmax + 1) % (jmax + 1);

				fin[i][j][k] = fout[ip][jp][k];
			}
		}

	
	//Macro-parameters
	for(int i = 0; i <= imax; i++)
		for (int j = 0; j <= jmax; j++)
		{
			double sum = 0.0;
			for (int k = 0; k < Q; k++)
			{
				sum += fin[i][j][k];
			}

			rho[i][j] = sum;
			psi[i][j] = 1.0 - exp(-rho[i][j]);
		}

	for(int i = 0; i <= imax; i++)
		for (int j = 0; j <= jmax; j++)
		{
			double Fcx = 0.0, Fcy = 0.0;
			for (int k = 0; k < Q; k++)
			{
				int ip = (i + cx[k] + imax + 1) % (imax + 1);
				int jp = (j + cy[k] + jmax + 1) % (jmax + 1);

				Fcx += w[k] * psi[ip][jp] * cx[k];
				Fcy += w[k] * psi[ip][jp] * cy[k];
			}
			Fcx *= -Gc * psi[i][j];
			Fcy *= -Gc * psi[i][j];

			double sumx = 0.0, sumy = 0.0;
			for (int k = 0; k < Q; k++)
			{
				sumx += cx[k] * fin[i][j][k];
				sumy += cy[k] * fin[i][j][k];
			}

			U[i][j] = (sumx + tau * Fcx) / rho[i][j];
			V[i][j] = (sumy + tau * Fcy) / rho[i][j];
		}
}

void Output_plt(int iter)
{
	ostringstream f1;
	f1 << "fpcTecplot_" << iter << ".plt";
	f1 << scientific; //使用科学计数
	f1.precision(4);  //设置精度为4
	ofstream out(f1.str().c_str());

	out << "VARIABLES = \"X\"\n";
	out << "\"Y\"" << endl;
	out << "\"Rho\"" << endl;
	out << "\"U\"" << endl;
	out << "\"V\"" << endl;
	out << "\"Vel\"" << endl;
	out << "\"P\"" << endl;
	out << "ZONE T=\"G1\" " << "I=" << imax + 1 << ',' << ' ' << "J=" << jmax + 1 << endl;
	out << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)" << endl;

	for (int j = 0; j <= jmax; j++)
		for (int i = 0; i <= imax; i++)
		{
			P[i][j] = rho[i][j] / 3.0 + Gc * psi[i][j] * psi[i][j] / 6.0;
			out << i << '\t' << j << '\t' << rho[i][j] << '\t' << U[i][j] << '\t' << V[i][j] << '\t'
				<< sqrt(U[i][j] * U[i][j] + V[i][j] * V[i][j]) << '\t' << P[i][j] << endl;
		}
}
