#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void main(int argc, char **argv)
{
	int i;
	double dr, r, rF, f1, f2, k, T_M, T_W;
	double mu, gamma, psi, q, G0, EG, p;
	double B0, Ediff, n, m, Edisp, A0;
	double T, Tref;
	double eta_diff, eta_disp, GI, GG, A, B, fI, tau;
	double temp;
	double dt, t, vis;
	double y2s = 3.1536e7;
	double visc_A();
	double visc_B();
	double visc_GG();
	double visc_fI();

	mu = 0.72;
	gamma = 1.0;

	q = 4.0;
	p = 2.0;

	T_M = 1500.0;
	T_W = 285.0 ;
	Tref = 2500.0;

	f1 = 1e-12;
	f2 = 1e-2;
	k = 2;

	n = 3.0;
	A0 = 1.1e5 / pow(1e6, n);
	Edisp = 530.0 * 1e3;

	m = 3.0;
	B0 = 13.6 * pow(1e-6, m) / 1e6;
	Ediff = 300 * 1e3;

	G0 = 2e4 * pow(1e-6, p);
	EG = 300 * 1e3;

	dt = 1.0e1;

	r = 1e-2;
	tau = 10e6;

	t = 0.0;

	i = 0;
	while (t < 100e6)
	{
		T = 1100.0;
		A = visc_A(A0, Edisp, T);
		B = visc_B(B0, Ediff, T, m);
		GG = visc_GG(G0, EG, T);
		GI = q/p * pow(1e-6, q-p) * GG / 250.0;
		fI = visc_fI(f1, f2, T_M, T_W, T, k);
		fI = 5e-3;

		temp = B / (A * pow(tau, n - 1));
		rF = pow(temp, 1.0 / m);

		psi = 2 * A * pow(tau, n+1) * (1.0 + pow(rF/r, m));

		dr = mu * GI / (q * pow(r, q-1)) - fI * r * r * psi / (gamma * mu);

		r += dr * dt * y2s;

		vis = 0.5 * tau / (A*pow(tau,n) + B*tau/pow(r, m));
		t += dt;
		i++;
		if(i % 50 == 0)
			printf("%e %e %e\n",t, r, vis);
	}

	return;
}

double visc_A(A0, Edisp, T)
	double A0, Edisp, T;
{
	return A0*exp(0.0-Edisp/(8.3145*T));
}

double visc_B(B0, Ediff, T, m)
	double B0, Ediff, T, m;
{
	return pow(M_PI*0.5, 0.0-m)*B0*exp(0.0-Ediff/(8.3145*T));
}

double visc_GG(G0, EG, T)
	double G0, EG, T;
{
	return G0 * exp(0.0 - EG / (8.3145 * T));
}

double visc_fI(f1, f2, T_M, T_W, T, k)
	double f1, f2, T_M, T_W, T, k;
{
	double value, temp;

	temp = (pow(T_M, k) - pow(T, k)) / (pow(T_M, k) - pow(T_W, k));
	value = f1 * pow(f2/f1, temp);

	return value;
}
