#include "global_variables.h"

void visc_from_grain_size(lev)
	int lev;
{
	static int been_here=0;
	int i;
	double dr, r, rF, f1, f2, k, T_M, T_W;
	double mu, gamma, psi, q, G0, EG, p;
	double B0, Ediff, n, m, Edisp, A0;
	double T, Tref;
	double eta_diff, eta_disp, GI, GG, A, B, fI, tau;
	double temp;
	double dt, dt_ad;
	static double t = 0.0;
	double y2s = 3.1536e7;

	double stress_scale;
	double time_scale;
	double vis_ref, kappa, D, alpha, g, rho;

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


	alpha = 1e-5;
	kappa = 1e-6;
	D = 2890e3;
	rho = 3300.0;
	g = 9.8;

	vis_ref =  (rho * g * alpha * Tref * D * D * D) / (G.Ra * kappa);
	stress_scale = (vis_ref * kappa) / (D * D);

	time_scale = D * D / kappa;

	G.gs[lev] = realloc(G.gs[lev], G.nno[lev] * sizeof(double));
	if(been_here == 0)
	{
		for(i=0; i<G.nno[lev]; i++)	
		{
			G.gs[lev][i] = 1e-3;
		}
		been_here = 1;
		return;
	}

	for(i=0; i<G.nno[lev]; i++)
		G.stress[lev][i] *= stress_scale;

	dt = 1e2 * y2s;//need to be small enough to converge

	while (t < G.model_time * time_scale)
	{
		for(i=0; i<G.nno[lev]; i++)	
		{
			if(2890 * (1.0 - G.X[lev][1][i]) > 200.0)
			{
				G.stress[lev][i] = 0.0;
				continue;
			}

			T = G.T[lev][i] * Tref + 273.5 + 2890.0*(1.0-G.X[lev][1][i])*0.3;

			r = G.gs[lev][i];

			A = visc_A(A0, Edisp, T);
			B = visc_B(B0, Ediff, T, m);

			tau = G.stress[lev][i];

			GG = visc_GG(G0, EG, T);
			GI = q/p * pow(1e-6, q-p) * GG / 250.0;
			fI = visc_fI(f1, f2, T_M, T_W, T, k);

			temp = B / (A * pow(tau, n - 1));
			rF = pow(temp, 1.0 / m);

			psi = 2 * A * pow(tau, n+1) * (1.0 + pow(rF/r, m));

			dr = mu * GI / (q * pow(r, q-1)) - fI * r * r * psi / (gamma * mu);

			G.gs[lev][i] += dr * dt;

			G.vis[lev][i] = 0.5 * tau / (A*pow(tau,n) + B*tau/pow(G.gs[lev][i], m));

			G.vis[lev][i] /= vis_ref;

			if(G.vis[lev][i] < G.viscosity.lowest) G.vis[lev][i] = G.viscosity.lowest;
			if(G.vis[lev][i] > G.viscosity.highest) G.vis[lev][i] = G.viscosity.highest;
		}
		t += dt;
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

void compute_stress()
{
	int lev;
	int i;
        double eta, sr;
        double u1_1,u2_2,u1_2,u2_1;
        double dNa_dx,dNa_dz;
	double *stress;


	int e,a,d,node;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.stress[lev] = realloc(G.stress[lev], G.nno[lev] * sizeof(double));

		stress = malloc(G.nel[lev] * sizeof(double));
		for(e=0; e<G.nel[lev]; e++)
		{
			eta=0.0;
			for(a=0; a<NEN; a++)
			{
				i = node_index(G.ien[lev][a][e], lev);
				eta += G.N[a].center * G.vis[lev][i];
			}

			u1_1=u1_2=u2_2=0.0;
			for(a=0; a<NEN; a++)
			{
				node=G.ien[lev][a][e];
				i = node_index(node, lev);

				dNa_dx=G.dN_dxi[ a].center*G.dz_deta[lev][e]
				      -G.dN_deta[a].center*G.dz_dxi;
				dNa_dz=G.dN_deta[a].center*G.dx_dxi[lev][e]
				      -G.dN_dxi[ a].center*G.dx_deta;
				u1_1+=dNa_dx*G.V[lev][0][i];
				u2_2+=dNa_dz*G.V[lev][1][i];
				u1_2+=dNa_dz*G.V[lev][0][i];
			}

			u1_1 /= G.jacobi[lev][e];
			u2_2 /= G.jacobi[lev][e];
			u1_2 /= G.jacobi[lev][e];
			u2_1 = u1_2;

			sr = u1_1*u1_1 + u2_2*u2_2 + u1_2*u2_1*2.0;
			sr = sqrt(0.5*sr);

			stress[e] = 2.0 * eta * sr;
		}


		p_to_nodes(stress, G.stress[lev], lev);

		free(stress);
	}

	return;
}
