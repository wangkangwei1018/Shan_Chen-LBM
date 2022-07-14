#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;

const int Q = 9;
const int NX = 100;
const int NY = 100;

const double w[Q] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
const int cx[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int cy[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

double arho[NX + 1][NY + 1], brho[NX + 1][NY + 1], au[NX + 1][NY + 1], av[NX + 1][NY + 1], bu[NX + 1][NY + 1], bv[NX + 1][NY + 1], u[NX + 1][NY + 1], v[NX + 1][NY + 1];
double afin[NX + 1][NY + 1][Q], afout[NX + 1][NY + 1][Q], bfin[NX + 1][NY + 1][Q], bfout[NX + 1][NY + 1][Q];

double atau = 1.0, btau = 1.0;
double Gc = 1.0;
double Feq(int k, double rho, double u, double v);
double xForce(int i, int j, int component);
double yForce(int i, int j, int component);
double aphi(int i, int j);
double bphi(int i, int j);

void Init();
void Evolution();
void Output(int iter);

int main()
{
	Init();
	for (int iter = 0; iter <= 80000; iter++)
	{
		Evolution();

		if (iter % 2000 == 0)
			Output(iter);
	}

	return 0;
}


void Init()
{
	for(int i = 0; i <= NX; i++)
		for (int j = 0; j <= NY; j++)
		{
			arho[i][j] = 0;
			brho[i][j] = 1;
			au[i][j] = av[i][j] = 0.0;
			bu[i][j] = bv[i][j] = 0.0;
			u[i][j] = v[i][j] = 0.0;

			if (j >= 0 && hypot(i - NX / 2, j) < 10)
			{
				arho[i][j] = 1;
				brho[i][j] = 0;
			}
		}


	for(int i = 0; i <= NX; i++)
		for (int j = 0; j <= NY; j++)
		{
			for (int k = 0; k < Q; k++)
			{
				afin[i][j][k] = Feq(k, arho[i][j], au[i][j], av[i][j]);
				bfin[i][j][k] = Feq(k, brho[i][j], bu[i][j], bv[i][j]);
			}
		}
}

void Evolution()
{
	int i, j, k;
	double deltfs = 0.1;

	for(i = 0; i <= NX; i++)
		for(j = 0; j <= NY; j++)
			for (k = 0; k < Q; k++)
			{
				afout[i][j][k] = afin[i][j][k] - (afin[i][j][k] - Feq(k, arho[i][j], au[i][j], av[i][j])) / atau;
				bfout[i][j][k] = bfin[i][j][k] - (bfin[i][j][k] - Feq(k, brho[i][j], bu[i][j], bv[i][j])) / btau;
			}

	for(i = 0; i <= NX; i++)
		for(j = 0; j <= NY; j++)
			for (k = 0; k < Q; k++)
			{
				int ip = (i - cx[k] + NX + 1) % (NX + 1);
				int jp = (j - cy[k] + NY + 1) % (NY + 1);

				afin[i][j][k] = afout[ip][jp][k];
				bfin[i][j][k] = bfout[ip][jp][k];
			}

	for(i = 0; i <= NX; i++)
		for (j = 0; j <= NY; j++)
		{
			arho[i][j] = 0;
			brho[i][j] = 0;
			au[i][j] = av[i][j] = 0;
			bu[i][j] = bv[i][j] = 0;

			for (k = 0; k < Q; k++)
			{
				arho[i][j] += afin[i][j][k];
				brho[i][j] += bfin[i][j][k];

				au[i][j] += afin[i][j][k] * cx[k];
				av[i][j] += afin[i][j][k] * cy[k];
				bu[i][j] += bfin[i][j][k] * cx[k];
				bv[i][j] += bfin[i][j][k] * cy[k];
			}

			u[i][j] = (au[i][j] / atau + bu[i][j] / btau) / (arho[i][j] / atau + brho[i][j] / btau);
			v[i][j] = (av[i][j] / atau + bv[i][j] / btau) / (arho[i][j] / atau + brho[i][j] / btau);

			au[i][j] = u[i][j] + atau * xForce(i, j, 1) / (arho[i][j] + 1.0e-32);
			av[i][j] = v[i][j] + atau * yForce(i, j, 1) / (arho[i][j] + 1.0e-32);
			bu[i][j] = u[i][j] + btau * xForce(i, j, 2) / (brho[i][j] + 1.0e-32);
			bv[i][j] = v[i][j] + btau * yForce(i, j, 2) / (brho[i][j] + 1.0e-32);
		}
}

void Output(int iter)
{
	ostringstream f1;
	f1 << "fpcTecplot_" << iter << ".plt";
	f1 << scientific; //使用科学计数
	f1.precision(4);  //设置精度为4
	ofstream out(f1.str().c_str());

	out << "VARIABLES = \"X\"\n";
	out << "\"Y\"" << endl;
	out << "\"arho\"" << endl;
	out << "\"brho\"" << endl;
	out << "\"au\"" << endl;
	out << "\"av\"" << endl;
	out << "\"bu\"" << endl;
	out << "\"bv\"" << endl;
	out << "\"U\"" << endl;
	out << "\"V\"" << endl;
	out << "\"P\"" << endl;
	out << "ZONE T=\"G1\" " << "I=" << NX + 1 << ',' << ' ' << "J=" << NY + 1 << endl;
	out << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)" << endl;

	for (int j = 0; j <= NY; j++)
		for (int i = 0; i <= NX; i++)
		{
			out << i << '\t' << j << '\t' << arho[i][j] << '\t' << brho[i][j] << '\t' << au[i][j] << '\t' << av[i][j] << '\t'
				<< bu[i][j] << '\t' << bv[i][j] << '\t' << u[i][j] << '\t' << v[i][j] << '\t'
				<< (arho[i][j] + brho[i][j]) / 3. + Gc / 3. * arho[i][j] * brho[i][j] << endl;
		}
}

double Feq(int k, double rho, double u, double v)
{
	double t1, t2;
	t1 = u * cx[k] + v * cy[k];
	t2 = u * u + v * v;

	return w[k] * rho * (1.0 + 3.0 * t1 + 4.5 * t1 * t1 - 1.5 * t2);
}

double xForce(int i, int j, int component)
{
	double xforce = 0.0;

	if (component == 1)
	{
		xforce = -Gc * aphi(i, j) * ((bphi(i + 1, j) - bphi(i - 1, j)) / 3. + (bphi(i + 1, j + 1) - bphi(i - 1, j + 1) - bphi(i - 1, j - 1) + bphi(i + 1, j - 1)) / 12.);
	}
	if (component == 2)
	{
		xforce = -Gc * bphi(i, j) * ((aphi(i + 1, j) - aphi(i - 1, j)) / 3. + (aphi(i + 1, j + 1) - aphi(i - 1, j + 1) - aphi(i - 1, j - 1) + aphi(i + 1, j - 1)) / 12.);
	}

	return xforce;
}

double yForce(int i, int j, int component)
{
	double yforce = 0.0;

	if (component == 1)
	{	
		yforce = -Gc * aphi(i, j) * ((bphi(i, j + 1) - bphi(i, j - 1)) / 3. + (bphi(i + 1, j + 1) + bphi(i - 1, j + 1) - bphi(i - 1, j - 1) - bphi(i + 1, j - 1)) / 12.);
	}
	if (component == 2)
	{
		yforce = -Gc * bphi(i, j) * ((aphi(i, j + 1) - aphi(i, j - 1)) / 3. + (aphi(i + 1, j + 1) + aphi(i - 1, j + 1) - aphi(i - 1, j - 1) - aphi(i + 1, j - 1)) / 12.);
	}

	return yforce;
}

double aphi(int i, int j)
{
	double aphi;
	int ip = (i + NX + 1) % (NX + 1);
	int jp = (j + NY + 1) % (NY + 1);

	aphi = arho[ip][jp];

	return aphi;
}

double bphi(int i, int j)
{
	double bphi;

	int ip = (i + NX + 1) % (NX + 1);
	int jp = (j + NY + 1) % (NY + 1);

	bphi = brho[ip][jp];

	return bphi;
}
