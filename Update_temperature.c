#include "global_variables.h"

void update_temperature()
{
	static int been_here = 0;
	int lev = G.max_level - 1;
	int i,j,neq;
	double *B,*d0;
	double *CT, *CTdot;
	double val;

	double **Arow_val;
	int *Ancol, **Arow_col;

	neq = G.heat[lev].neq;

	calculate_artificial_diffusivity(lev);
	K_matrix_heat(lev);
	F_matrix_heat(lev);
	C_matrix_heat(lev);

	other_KFC_heat_contribution(lev);

	if(G.lg_mul)
		build_lagrange_heat_KC_matrix(lev);

	if(G.control.AMR == 1 || G.solution_cycle == 0)
		setup_MPI_comm_Tframe(lev);

	if(been_here == 0)
	{
		initial_T_dot(lev);
		been_here=1;
	}


	G.heat[lev].dt = advection_dt(lev) * 0.75;
	//if(G.ip == 0)
	//	printf("%e\n",G.heat[lev].dt);
	//G.heat[lev].dt = 1e-4;//advection_dt(lev);

	G.model_time += G.heat[lev].dt;

	B = malloc(neq * sizeof(double));
	d0 = malloc(neq * sizeof(double));

	for(i=0; i<neq; i++)
		d0[i] = G.heat[lev].d[i];

	CT = malloc(neq * sizeof(double));
	CTdot = malloc(neq * sizeof(double));

	K_prod_d_heat(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, G.heat[lev].d, CT, neq, neq, lev, 'T');
	K_prod_d_heat(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, G.heat[lev].ddot, CTdot, neq, neq, lev, 'T');

	/*
	for(i=0; i<neq; i++)
		printf("here: %3d %3d %3d %3d %+.10e %+.10e %+.10e %+.10e\n",G.ip, i, G.heat[lev].eq_node[i], G.smpi[lev].Teq_nsh[i], G.heat[lev].d[i], CT[i], G.heat[lev].ddot[i], CTdot[i]);
	terminate();
	*/

	for(i=0; i<neq; i++)
	{
		B[i] = G.heat[lev].F[i] + CT[i] / (G.control.theta * G.heat[lev].dt);
		    			+ CTdot[i] * (1.0 - G.control.theta) / G.control.theta;
	}

        Ancol = malloc(neq * sizeof(int));
        Arow_col = malloc(neq * sizeof(int *)); 
        Arow_val = malloc(neq * sizeof(double *));
        for(i=0; i<neq; i++)
	{
                Ancol[i] = G.heat[lev].Kncol[i];
                Arow_col[i] = malloc(Ancol[i] * sizeof(int));
                Arow_val[i] = malloc(Ancol[i] * sizeof(double));
		for(j=0; j<Ancol[i]; j++)
		{
			val = G.heat[lev].Crow_val[i][j] / (G.control.theta * G.heat[lev].dt)
			    + G.heat[lev].Krow_val[i][j];

			Arow_val[i][j] = val;
			Arow_col[i][j] = G.heat[lev].Krow_col[i][j];
		}
	}


	/*
	if(G.np == 1)
		gauss_seidel(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, G.control.accuracy, 400);
	else
	{
		//Conj_Grad_heat(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, G.control.accuracy, 400, 'T');
		Conj_Grad_heat(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, 1e-3, G.control.max_citcom_cycle, 'T');
//		jacobi_heat(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, G.control.accuracy, 400, 'T');
	}
	*/

	//jacobi_heat(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, 1e-3, 400, 'T');
	Conj_Grad_heat(Arow_val, Arow_col, Ancol, B, G.heat[lev].d, neq, lev, 1e-3, 1000, 'T');

        for(i=0; i<G.nno[lev]; i++)
        {
                if(G.heat[lev].id[i] >= 0)
                {
                        G.T[lev][i] = G.heat[lev].d[G.heat[lev].id[i]];
                }
                else
                {
                        G.T[lev][i]=G.bc.T_bc_val[lev][i];
                }
        }

        for(i=0; i<neq; i++)
        {
                G.heat[lev].ddot[i]=1.0/(G.control.theta*G.heat[lev].dt)*(G.heat[lev].d[i]-d0[i])
                        -(1.0-G.control.theta)/G.control.theta*G.heat[lev].ddot[i];
        }

        free(B);
        free(d0);
	free(Ancol);
        for(i=0; i<neq; i++)
	{
		free(Arow_col[i]);
		free(Arow_val[i]);
	}
	free(Arow_col);
	free(Arow_val);
	free(CT);
	free(CTdot);

	return;
}

double advection_dt(lev)
	int lev;
{
        double dt,ts;
        double u1,u2,u,step;
        double diff_timestep;
	double L[NSD];
	double global_dt;
        int e,a,d,node,level;

        diff_timestep=1e10;
        for(e=0; e<G.nel[lev]; e++)
        {
		level = Mcode_lev(G.Mcode[lev][e]);
		L[0] = G.dx[level];
		L[1] = G.dz[level];

                for(d=0; d<NDF; d++)
                {
                        ts = L[d] * L[d];
                        if(diff_timestep > ts)
                                diff_timestep = ts;
                }
        }

        dt=1e10;
        for(e=0; e<G.nel[lev]; e++)
        {
                u = u1 = u2 = 0.0;
                for(a=0; a<NEN; a++)
                {
                        node = node_index(G.ien[lev][a][e], lev);
                        u1 += G.N[a].center * G.V[lev][0][node];
                        u2 += G.N[a].center * G.V[lev][1][node];
                }
                u = fabs(u1) / L[0] + fabs(u2) / L[1];
                step = 0.5 / u;
                dt = MIN(step, dt);
        }
        dt = MIN(dt, diff_timestep);
	MPI_Allreduce(&dt,&global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        return global_dt;
}

void initial_T_dot(lev)
	int lev;
{
	double *B;
	int i,j,neq,neq_m;

	neq = G.heat[lev].neq;
	neq_m = G.heat[lev].neq_m;
	G.heat[lev].d = realloc(G.heat[lev].d, neq * sizeof(double));
	G.heat[lev].ddot = realloc(G.heat[lev].ddot, neq * sizeof(double));

	j = 0;
	for(i=0; i<G.nno[lev]; i++)
	{
		if(G.heat[lev].id[i] >= 0)
		{
			G.heat[lev].d[j++] = G.T[lev][i];
		}
	}

	for(i=neq_m; i<neq; i++)
		G.heat[lev].d[i] = 0.0;

	B = malloc(neq * sizeof(double));

	K_prod_d_heat(G.heat[lev].Krow_val, G.heat[lev].Krow_col, G.heat[lev].Kncol, G.heat[lev].d, B, neq, neq, lev, 'T');

	/*
	for(i=0; i<neq; i++)
		printf("%4d %4d %4d %4d %+16e\n",G.ip, i, G.heat[lev].eq_node[i], G.smpi[lev].Teq_nsh[i], B[i]);
	terminate();
	*/



	for(i=0; i<neq; i++)
	{
		B[i] = G.heat[lev].F[i] - B[i];
	}

	for(i=0; i<neq; i++)
		G.heat[lev].ddot[i] = 0.0;

	/*
	if(G.np == 1)
	{
		gauss_seidel(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, B, G.heat[lev].ddot, neq, lev, G.control.accuracy, 400);
	}
	else
	{
//		jacobi_heat(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, B, G.heat[lev].ddot, neq, lev, G.control.accuracy, 400, 'T');
		Conj_Grad_heat(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, B, G.heat[lev].ddot, neq, lev, G.control.accuracy, 400, 'T');
	}
	*/


	Conj_Grad_heat(G.heat[lev].Crow_val, G.heat[lev].Crow_col, G.heat[lev].Cncol, B, G.heat[lev].ddot, neq, lev, 1e-3, 400, 'T');

	free(B);
	return;
}

void calculate_artificial_diffusivity(lev)
	int lev;
{
        int a,e,d,j,level;
        double u[NSD],u_scalar;
        double alpha[NSD];
        double xi[NSD];
        double k_hat;
	double L[NSD];
        double diffusivity[2][2]={{1.0,0.0},
                                  {0.0,1.0}};
        for(d=0; d<NSD; d++)
	{
		G.heat[lev].tau[d] = realloc(G.heat[lev].tau[d], G.nel[lev] * sizeof(double));
	}

        for(e=0; e<G.nel[lev]; e++)
        {
		level = Mcode_lev(G.Mcode[lev][e]);
		L[0] = G.dx[level];
		L[1] = G.dz[level];

                for(d=0; d<NSD; d++)
                {
                        u[d]=0.0;
                        for(a=0; a<NEN; a++)
                        {
				j = node_index(G.ien[lev][a][e], lev);
                                u[d] += G.N[a].center * G.V[lev][d][j];
                        }

                        alpha[d] = 0.5 * u[d] * L[d] / diffusivity[d][d];

                        if(fabs(alpha[d]) < G.control.accuracy)
                        {
                                xi[d]=0.0;
                        }
                        else
                        {
                                xi[d]=cosh(alpha[d])/sinh(alpha[d]) - 1.0/alpha[d];
                        }
                }
                u_scalar=(u[0] * u[0] + u[1] * u[1]);
                k_hat=(xi[0] * u[0] * L[0] + xi[1] * u[1] * L[1]) * 0.5;
                for(d=0; d<NSD; d++)
                {
                        G.heat[lev].tau[d][e] = k_hat * u[d] / u_scalar;
                }
        }
	return;
}

void K_matrix_heat(lev)
	int lev;
{
	int i,e,a,b,P,Q,nel,neq,j,c,kk;
	int index_arr[VPTS];
	double K[VPTS][VPTS];

	neq = G.heat[lev].neq;
	nel = G.nel[lev];

        for(i=0; i<neq; i++)
        {
                G.heat[lev].Kncol[i] = 0;
                for(j=0; j<MAX_T_NCOL; j++)
                        G.heat[lev].Krow_val[i][j] = 0.0;
        }

	G.heat[lev].F = realloc(G.heat[lev].F, neq * sizeof(double));
	for(i=0; i<neq; i++)
		G.heat[lev].F[i] = 0.0;

        for(e=0; e<nel; e++)
        {
                for(b=0; b<VPTS; b++)
                        index_arr[b] = node_index(G.ien[lev][b][e],lev);

                get_element_K_heat(e, K, lev);
                for(a=0; a<NEN; a++)
                {
                        P = G.heat[lev].lm[a][e];
                        for(b=0; b<NEN; b++)
                        {
                                Q = G.heat[lev].lm[b][e];

				if(P >= 0 && Q >= 0)
				{
					c = value_in_array(Q, G.heat[lev].Krow_col[P], G.heat[lev].Kncol[P]);
					if(c >= 0)
					{
						G.heat[lev].Krow_val[P][c] += K[a][b];
					}
					else
					{
						G.heat[lev].Krow_col[P][ G.heat[lev].Kncol[P] ] = Q;
						G.heat[lev].Krow_val[P][ G.heat[lev].Kncol[P] ] += K[a][b];
						G.heat[lev].Kncol[P]++;
					}
				}

				if(P >= 0)
				{
					j = index_arr[b]; 
					if(G.bc.T_bc[lev][j] == 'T')
					{
						G.heat[lev].F[P] -= K[a][b] * G.bc.T_bc_val[lev][j];
					}
				}
                        }
                }
        }

	/*
        for(i=0; i<neq; i++)
        {
                G.heat[lev].Krow_col[i] = realloc(G.heat[lev].Krow_col[i], G.heat[lev].Kncol[i] * sizeof(int));
                G.heat[lev].Krow_val[i] = realloc(G.heat[lev].Krow_val[i], G.heat[lev].Kncol[i] * sizeof(double));
        }
	*/

	return;
}

void get_element_K_heat(e,k,lev)
	int e,lev;
	double k[VPTS][VPTS];
{
        int a,b,node,g;
        double dNa_dx,dNa_dz;
        double dNb_dx,dNb_dz;
        double temp1,temp2,temp3;
        double Na_hat;
        double u1,u2;
        double diffusivity[2][2]={{1.0,0.0},
                                  {0.0,1.0}};

        u1 = u2 = 0.0;
        for(a=0; a<NEN; a++)
        {
		node = node_index(G.ien[lev][a][e],lev);
                u1 += G.V[lev][0][node] * G.N[a].center;
                u2 += G.V[lev][1][node] * G.N[a].center;
        }

        for(a=0; a<NEN; a++)
        {
                for(b=0; b<NEN; b++)
                {
                        k[a][b] = 0.0;
                        for(g=0; g<NEG; g++)
                        {

                                dNa_dx=G.dN_dxi[ a].gauss[g]*G.dz_deta[lev][e]
                                      -G.dN_deta[a].gauss[g]*G.dz_dxi;
                                dNa_dz=G.dN_deta[a].gauss[g]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ a].gauss[g]*G.dx_deta;

                                dNb_dx=G.dN_dxi[ b].gauss[g]*G.dz_deta[lev][e]
                                      -G.dN_deta[b].gauss[g]*G.dz_dxi;
                                dNb_dz=G.dN_deta[b].gauss[g]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ b].gauss[g]*G.dx_deta;

                                Na_hat=G.heat[lev].tau[0][e] * dNa_dx
                                      +G.heat[lev].tau[1][e] * dNa_dz;

                                temp1 = u1 * G.N[a].gauss[g] * dNb_dx
                                      + u2 * G.N[a].gauss[g] * dNb_dz;

                                temp2 = u1 * Na_hat * dNb_dx
                                      + u2 * Na_hat * dNb_dz;

                                temp3=dNa_dx*diffusivity[0][0]*dNb_dx
                                     +dNa_dx*diffusivity[0][1]*dNb_dz
                                     +dNa_dz*diffusivity[1][0]*dNb_dx
                                     +dNa_dz*diffusivity[1][1]*dNb_dz;
                                k[a][b] += temp1+(temp2+temp3)/G.jacobi[lev][e];
			}
		}
	}

	return;
}

void F_matrix_heat(lev)
	int lev;
{
        int a,e,P;
        double f[VPTS];

        for(e=0; e<G.nel[lev]; e++)
        {
                get_element_F_heat(e,f,lev);
                for(a=0; a<NEN; a++)
                {
                        P = G.heat[lev].lm[a][e];
                        if(P >= 0)
			{
                                 G.heat[lev].F[P] += f[a];
			}
                }
        }

	return;
}

void get_element_F_heat(e,f,lev)
        int e,lev;
        double f[VPTS];
{
        int a,d,node;
        double Na_hat;
	double dNa_dx,dNa_dz;
        int natural_bc,level;
        double L;


	level = Mcode_lev(G.Mcode[lev][e]);

        natural_bc = 0;
        for(a=0; a<NEN; a++)
        {
                node = node_index(G.ien[lev][a][e],lev);
                if(G.bc.T_bc[lev][node] == 'F')
                {
                        natural_bc += a;
                }
        }
        if(natural_bc==3)//horizontal bc
        {
                natural_bc=1;
        }
        else if(natural_bc==1||natural_bc==5)//vertical bc
        {
                natural_bc=2;
        }
        else if(natural_bc==0 || natural_bc==1
             || natural_bc==2 || natural_bc==3)//ignore corner of domain
        {
                natural_bc=0;
        }
        for(a=0; a<NEN; a++)
        {
                f[a] = 0.0;
                for(d=0; d<NEG; d++)
                {
                                dNa_dx=G.dN_dxi[ a].gauss[d]*G.dz_deta[lev][e]
                                      -G.dN_deta[a].gauss[d]*G.dz_dxi;
                                dNa_dz=G.dN_deta[a].gauss[d]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ a].gauss[d]*G.dx_deta;

                        Na_hat=G.heat[lev].tau[0][e]*dNa_dx
                              +G.heat[lev].tau[1][e]*dNa_dz;

                        f[a]+=(G.N[a].gauss[d]+Na_hat) * G.control.Q;
                }
                f[a] *= G.jacobi[lev][e];

                if(natural_bc!=0)
                {
                        node=node_index(G.ien[lev][a][e],lev);
                        L=G.dz[level];
                        f[a] -= 0.5 * L * G.bc.T_bc_val[lev][node];
                }

        }
        return;
}

void C_matrix_heat(lev)
	int lev;
{
        double c[VPTS][VPTS];
        int e,a,P,Q,b,i,j,neq,nel,k,kk,*size;

	neq = G.heat[lev].neq;
	nel = G.nel[lev];

        for(i=0; i<neq; i++)
        {
                G.heat[lev].Cncol[i] = 0;
                for(j=0; j<MAX_T_NCOL; j++)
                        G.heat[lev].Crow_val[i][j] = 0.0;
        }

        for(e=0; e<nel; e++)
        {
                get_element_C_heat(e,c,lev);
                for(a=0; a<NEN; a++)
                {
                        P = G.heat[lev].lm[a][e];
                        for(b=0; b<NEN; b++)
                        {
                                Q = G.heat[lev].lm[b][e];

				if(P >= 0 && Q >= 0)
				{
					k = value_in_array(Q, G.heat[lev].Crow_col[P], G.heat[lev].Cncol[P]);
					if(k >= 0)
					{
						G.heat[lev].Crow_val[P][k] += c[a][b];
					}
					else
					{
						G.heat[lev].Crow_col[P][ G.heat[lev].Cncol[P] ] = Q;
						G.heat[lev].Crow_val[P][ G.heat[lev].Cncol[P] ] += c[a][b];
						G.heat[lev].Cncol[P]++;
					}
				}
                        }
                }
        }

	/*
        for(i=0; i<neq; i++)
        {
                G.heat[lev].Crow_col[i] = realloc(G.heat[lev].Crow_col[i], G.heat[lev].Cncol[i] * sizeof(int));
                G.heat[lev].Crow_val[i] = realloc(G.heat[lev].Crow_val[i], G.heat[lev].Cncol[i] * sizeof(double));
        }
	*/

	return;
}


void get_element_C_heat(e,c,lev)
        int e,lev;
        double c[VPTS][VPTS];
{
        int a,b,d;
        double Na_hat,dNa_dx,dNa_dz;
        for(a=0; a<NEN; a++)
        {
                for(b=0; b<NEN; b++)
                {
                        c[a][b]=0.0;
                        for(d=0; d<NEG; d++)
                        {
                                dNa_dx=G.dN_dxi[ a].gauss[d]*G.dz_deta[lev][e]
                                      -G.dN_deta[a].gauss[d]*G.dz_dxi;
                                dNa_dz=G.dN_deta[a].gauss[d]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ a].gauss[d]*G.dx_deta;

                                Na_hat=G.heat[lev].tau[0][e]*dNa_dx
                                      +G.heat[lev].tau[1][e]*dNa_dz;

                                c[a][b]+=(G.N[a].gauss[d]+Na_hat)*G.N[b].gauss[d];
                        }
                        c[a][b] *= G.jacobi[lev][e];
                }
        }

        return;
}


void other_KFC_heat_contribution(lev)
	int lev;
{
	int i,j,k,p,q,n,c,size;
	double *send_K,*recv_K;
	double *send_C,*recv_C;
	double *send_F,*recv_F;
	MPI_Status stats;

	//for F
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		send_F = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));

		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			p = G.smpi[lev].Teq[i][j];
			send_F[j] = G.heat[lev].F[p];
		}
		MPI_Send(send_F,G.smpi[lev].n_Teq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);

		free(send_F);
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		recv_F = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));

		MPI_Recv(recv_F,G.smpi[lev].n_Teq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			p = G.smpi[lev].Teq[i][j];
			G.heat[lev].F[p] += recv_F[j];
		}

		free(recv_F);
	}


	//for K
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			p = G.smpi[lev].Teq[i][j];
			send_K = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));
			n = 0;
			for(k=0; k<G.smpi[lev].n_Teq[i]; k++)
			{
				q = G.smpi[lev].Teq[i][k];

				c = value_in_array(q, G.heat[lev].Krow_col[p], G.heat[lev].Kncol[p]);

				if(c >= 0)
					send_K[n] = G.heat[lev].Krow_val[p][c];
				else
					send_K[n] = 0.0;
				n++;
			}
			MPI_Send(send_K,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
			free(send_K);
		}
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			recv_K = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));

			MPI_Recv(recv_K,G.smpi[lev].n_Teq[i]*G.smpi[lev].n_Teq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
			p = G.smpi[lev].Teq[i][j];
			size = G.heat[lev].Kncol[p];

			n=0;
			for(k=0; k<G.smpi[lev].n_Teq[i]; k++)
			{
				q = G.smpi[lev].Teq[i][k];

				c = value_in_array(q, G.heat[lev].Krow_col[p], G.heat[lev].Kncol[p]);

				if (c >= 0)
				{
					G.heat[lev].Krow_val[p][c] += recv_K[n];
				}
				else
				{
					if(fabs(recv_K[n])>0.0)
					{
						c = G.heat[lev].Kncol[p]++;

						/*
						if(G.heat[lev].Kncol[p] >= size)
						{
							size *= 2;
							G.heat[lev].Krow_val[p] = realloc(G.heat[lev].Krow_val[p], size * sizeof(double));
							G.heat[lev].Krow_col[p] = realloc(G.heat[lev].Krow_col[p], size * sizeof(int));
						}
						*/

						G.heat[lev].Krow_val[p][c] = recv_K[n];
						G.heat[lev].Krow_col[p][c] = q;
					}
				}
				n++;
			}
			/*
			G.heat[lev].Krow_val[p] = realloc(G.heat[lev].Krow_val[p], G.heat[lev].Kncol[p] * sizeof(double));
			G.heat[lev].Krow_col[p] = realloc(G.heat[lev].Krow_col[p], G.heat[lev].Kncol[p] * sizeof(int));
			*/
			free(recv_K);
		}

	}

	// for C
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			p = G.smpi[lev].Teq[i][j];
			send_C = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));
			n = 0;
			for(k=0; k<G.smpi[lev].n_Teq[i]; k++)
			{
				q = G.smpi[lev].Teq[i][k];

				c = value_in_array(q, G.heat[lev].Crow_col[p], G.heat[lev].Cncol[p]);

				if(c >= 0)
					send_C[n] = G.heat[lev].Crow_val[p][c];
				else
					send_C[n] = 0.0;
				n++;
			}
			MPI_Send(send_C,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
			free(send_C);
		}
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<G.smpi[lev].n_Teq[i]; j++)
		{
			recv_C = malloc(G.smpi[lev].n_Teq[i] * sizeof(double));

			MPI_Recv(recv_C,G.smpi[lev].n_Teq[i]*G.smpi[lev].n_Teq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
			p = G.smpi[lev].Teq[i][j];
			size = G.heat[lev].Cncol[p];
			n=0;
			for(k=0; k<G.smpi[lev].n_Teq[i]; k++)
			{
				q = G.smpi[lev].Teq[i][k];

				c = value_in_array(q, G.heat[lev].Crow_col[p], G.heat[lev].Cncol[p]);

				if (c >= 0)
				{
					G.heat[lev].Crow_val[p][c] += recv_C[n];
				}
				else
				{
					if(fabs(recv_C[n])>0.0)
					{
						c = G.heat[lev].Cncol[p]++;

						/*
						if(G.heat[lev].Cncol[p] >= size)
						{
							size *= 2;
							G.heat[lev].Crow_val[p] = realloc(G.heat[lev].Crow_val[p], size * sizeof(double));
							G.heat[lev].Crow_col[p] = realloc(G.heat[lev].Crow_col[p], size * sizeof(int));
						}
						*/

						G.heat[lev].Crow_val[p][c] = recv_C[n];
						G.heat[lev].Crow_col[p][c] = q;
					}
				}
				n++;
			}
			/*
			G.heat[lev].Crow_val[p] = realloc(G.heat[lev].Crow_val[p], G.heat[lev].Cncol[p] * sizeof(double));
			G.heat[lev].Crow_col[p] = realloc(G.heat[lev].Crow_col[p], G.heat[lev].Cncol[p] * sizeof(int));
			*/
			free(recv_C);
		}

	}

        //sort col for each row, this is a must
        int log;
        for(i=0; i<G.heat[lev].neq; i++)
        {
                log = 1;
                while(log)
                {
                        log = 0;
                        for(j=0; j<G.heat[lev].Kncol[i]-1; j++)
                        {
                                if(G.heat[lev].Krow_col[i][j] > G.heat[lev].Krow_col[i][j+1])
                                {
                                        swap(&G.heat[lev].Krow_col[i][j], &G.heat[lev].Krow_col[i][j+1]);
                                        swap_d(&G.heat[lev].Krow_val[i][j], &G.heat[lev].Krow_val[i][j+1]);
                                        log = 1;
                                }
                        }
                }
        }

        for(i=0; i<G.heat[lev].neq; i++)
        {
                log = 1;
                while(log)
                {
                        log = 0;
                        for(j=0; j<G.heat[lev].Cncol[i]-1; j++)
                        {
                                if(G.heat[lev].Crow_col[i][j] > G.heat[lev].Crow_col[i][j+1])
                                {
                                        swap(&G.heat[lev].Crow_col[i][j], &G.heat[lev].Crow_col[i][j+1]);
                                        swap_d(&G.heat[lev].Crow_val[i][j], &G.heat[lev].Crow_val[i][j+1]);
                                        log = 1;
                                }
                        }
                }
        }

	return;
}

void enforce_constrain_Teq(x, lev, data_type)
	int lev;
	double *x;
	char data_type;
{
	int i,eq,j,c,node,k,eqc,inode;
	double temp;
	double *eqval;
	int *hnode;
	int number,size;

	int *irecv;
	double *drecv;
	int iamount,damount;

	MPI_Status stats;

	number = 0;
	size = 1;
	hnode = malloc(size * sizeof(int));
	eqval = malloc(size * sizeof(double));


	for(i=0; i<G.hangle[lev].nh; i++)
	{
		node = G.hangle[lev].node[i];

		inode = G.hangle[lev].ni[i];

		hnode[number] = node;

		temp = 0.0;
		for(c=0; c<2; c++)
		{
			k = G.hangle[lev].ci[c][i];
			eqc = G.heat[lev].id[k];
			if(eqc != -1)
				temp += x[eqc];
			else
				temp += G.bc.T_bc_val[lev][k];
		}
		temp *= 0.5;

		eqval[number] = temp;

		number++;

		if(number >= size)
		{
			size *= 2;
			hnode = realloc(hnode, size * sizeof(int));
			eqval = realloc(eqval, size * sizeof(double));
		}

		if(inode >= 0)
		{
			eq = G.heat[lev].id[inode];
			x[eq] = temp;
		}
	}


	hnode = realloc(hnode, number * sizeof(int));
	eqval = realloc(eqval, number * sizeof(double));

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(hnode, number, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(eqval, number, MPI_DOUBLE, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}

	free(hnode);
	free(eqval);


	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &iamount);
                irecv = malloc(iamount * sizeof(int));
                MPI_Recv(irecv,iamount,MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		damount = iamount;
                drecv = malloc(damount * sizeof(double));
                MPI_Recv(drecv,damount,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		for(j=0; j<iamount; j++)
		{
			node = irecv[j];
			k = node_index(node,lev);
			if(k >= 0)
			{
				eq = G.heat[lev].id[k];
				if(eq >= 0)
				{
					x[eq] = drecv[j];
				}
			}
		}
		free(irecv);
		free(drecv);
	}

	return;
}

void K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev, data_type)
        double **matr_row_val;
        int **matr_row_col, *matr_ncol;
        double *arr,*prod;
        int nrow , ncol, lev;
        char data_type;
{
        int i, j, c;
	void other_Teq_contribution_slow();


        for(i=0; i<nrow; i++)
                prod[i] = 0.0;//do not forget to initialize, otherwiseprod[c] below is not be defined yet

        for(i=0; i<nrow; i++)
        {
                for(j=0; j<matr_ncol[i]; j++)
                {
                        c = matr_row_col[i][j];
                        prod[i] += matr_row_val[i][j] * arr[c];
                }
        }

        other_Teq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);
        other_lg_Teq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);

        return;
}


void jacobi_heat(matr_row_val, matr_row_col, matr_ncol, arr, x, m, lev, accu, max_cycle, data_type)
        double **matr_row_val;
        int **matr_row_col, *matr_ncol;
        double *arr, *x;
        double accu;
        int m,lev,max_cycle;
	int data_type;
{
        int i, j, c, cycle, idiag;
        double *x_tmp, *prod, diag;
        double diff, temp, diff_global;

        cycle = 0;
        x_tmp = malloc(m * sizeof(double));
        prod = malloc(m * sizeof(double));
        while(1)
        {
                for(i=0; i<m; i++)
                {
                        x_tmp[i] = x[i];
                        prod[i] = 0.0;
                }

                for(i=0; i<m; i++)
                {
			idiag = 0; diag = 1.0;
                        for(j=0; j<matr_ncol[i]; j++)
                        {
                                c = matr_row_col[i][j];
				if(c != i)
				{
	                                prod[i] += matr_row_val[i][j] * x_tmp[c];
				}
				else
				{
					idiag = 1;
					diag = matr_row_val[i][j];
				}
                        }
			if(idiag == 0)
			{
				fprintf(stderr,"wrong here\n");
			}
                        x[i] = (arr[i] - prod[i]) / diag;
                }


                other_Teq_contribution(matr_row_val,matr_row_col,matr_ncol,x,x_tmp,m,lev,0,data_type);

                enforce_constrain_Teq(x, lev);

                cycle++;

                diff = 0.0;
                for(i=0; i<m; i++)
                {
                        temp = fabs(x[i] - x_tmp[i]);
                        if(temp > diff)
                                diff = temp;
                }

                MPI_Allreduce(&diff,&diff_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

                if(diff_global < accu || cycle > max_cycle)
                        break;
        }


        free(prod);
        free(x_tmp);
        return;
}

void Conj_Grad_heat_new(matr_row_val, matr_row_col, matr_ncol, b, x, n, lev, accu, max_cycle, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *b, *x, accu;
	int n, lev, max_cycle;
	char data_type;
{
        int count, i, c;
        double *r0, *r1, *r2, *z0, *z1, *p1, *p2, *Ap, *shuffle;
        double alpha, beta, dotprod, dotr1z1, dotr0z0;
        double temp, residual;
	void Tproduct();


        r0 = malloc(n * sizeof(double));
        r1 = malloc(n * sizeof(double));
        r2 = malloc(n * sizeof(double));
        z0 = malloc(n * sizeof(double));
        z1 = malloc(n * sizeof(double));
        p1 = malloc(n * sizeof(double));
        p2 = malloc(n * sizeof(double));
        Ap = malloc(n * sizeof(double));

        for(i=0; i<n; i++)
        {
                r1[i] = b[i];
                x[i] = 0.0;
        }


	residual = global_T_prod(r1, r1, lev);


        assert(residual != 0.0  /* initial residual for CG = 0.0 */);
        count = 0;
        while(((residual > accu) && (count < max_cycle)) || count == 0)
        {
                for(i=0; i<n; i++)
			z1[i] = r1[i];

		dotr1z1 = global_T_prod(r1, z1, lev);

                if(count == 0)
                {
                        for(i=0; i<n; i++)
                                p2[i] = z1[i];
                }
                else
                {
                        assert(dotr0z0 != 0.0 /* in head of conj_grad */);
                        beta = dotr1z1 / dotr0z0;
                        for(i=0; i<n; i++)
                                p2[i] = z1[i] + beta * p1[i];
                }

                dotr0z0 = dotr1z1;

                K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, p2, Ap, n, n, lev, data_type);

		/*
		for(i=0; i<n; i++)
			printf("%d %d %d %e\n",G.ip, i, G.heat[lev].eq_node[i], Ap[i]);
		terminate();
		*/


		dotprod = global_T_prod(p2, Ap, lev);


                if(dotprod == 0.0)
                        alpha = 1e-3;
                else
                        alpha = dotr1z1 / dotprod;

                for(i=0; i<n; i++)
                {
                        x[i] += alpha * p2[i];
                        r2[i] = r1[i] - alpha * Ap[i];
                }
//		if(G.lg_mul == 0)
			enforce_constrain_Teq(x, lev);

		residual = global_T_prod(r2, r2, lev);

                shuffle = r0; r0 = r1; r1 = r2; r2 = shuffle;
                shuffle = z0; z0 = z1; z1 = shuffle;
                shuffle = p1; p1 = p2; p2 = shuffle;
                count++;
        }
	if(G.ip == 0)
		printf("Temperature CG cycle = %d\n", count);

        free(r0);
        free(r1);
        free(r2);
        free(z0);
        free(z1);
        free(p1);
        free(p2);
        free(Ap);
}

void Conj_Grad_heat(matr_row_val, matr_row_col, matr_ncol, b, x, n, lev, accu, max_cycle, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *b, *x, accu;
	int n, lev, max_cycle;
	char data_type;
{
	int i;
	double *r,*d,*s,*q;
	double residual,global_residual;
        double delta_new,delta_old,delta0;
        double alpha,beta;
        double temp1;
        double temp;
        int cycle;


	if(G.control.cg_precondition)
	{
		int c;
		r=malloc(n*sizeof(double));
		s=malloc(n*sizeof(double));
		d=malloc(n*sizeof(double));
		q=malloc(n*sizeof(double));


		K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
		for(i=0; i<n; i++)
		{
			r[i] = b[i] - r[i];
			c = value_in_sorted_array(i, matr_row_col[i], matr_ncol[i]);
			d[i] = r[i] / matr_row_val[i][c];
		}

		/*
		temp=0.0;
		for(i=0;i<n;i++)
			temp += r[i] * d[i];
		MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		*/

		delta_new = global_T_prod(r, d, lev);

		delta0 = delta_new;

		cycle = 0;
		while(cycle < max_cycle && delta_new > accu * delta0)
		{
			K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, d, q, n, n, lev, data_type);

			/*
			temp=0.0;
			for(i=0; i<n; i++)
				temp += d[i] * q[i];
			MPI_Allreduce(&temp, &temp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			temp1 = global_T_prod(d, q, lev);

			alpha = delta_new / temp1;

			for(i=0; i<n; i++)
				x[i] += alpha * d[i];

			enforce_constrain_Teq(x, lev);

			if(cycle % 50 == 0)
			{
				K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
				for(i=0; i<n; i++)
					r[i] = b[i] - r[i];
			}
			else
			{
				for(i=0; i<n; i++)
					r[i] -= alpha * q[i];
			}

			for(i=0; i<n; i++)
			{
				c = value_in_sorted_array(i, matr_row_col[i], matr_ncol[i]);
				s[i] = r[i] / matr_row_val[i][c];
			}

			delta_old = delta_new;

			/*
			temp = 0.0;
			for(i=0; i<n; i++)
				temp += r[i] * s[i];
			MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			delta_new = global_T_prod(r, s, lev);


			beta = delta_new / delta_old;
			for(i=0; i<n; i++)
				d[i] = s[i] + beta * d[i];

			cycle++;

			residual = 0.0;
			for(i=0; i<n; i++)
			{
				if(fabs(alpha * r[i]) > residual)
					residual = fabs(alpha * r[i]);
			}
			MPI_Allreduce(&residual,&global_residual,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

			if(global_residual < accu)
				break;

		}
		if(G.ip == 0)
			printf("Temperature CG cycle = %d\n",cycle);
		free(r);
		free(s);
		free(d);
		free(q);
	}
	else
	{
		r=malloc(n*sizeof(double));
		d=malloc(n*sizeof(double));
		q=malloc(n*sizeof(double));


		K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
		for(i=0; i<n; i++)
		{
			r[i] = b[i] - r[i];
			d[i] = r[i];
		}

		/*
		temp=0.0;
		for(i=0;i<n;i++)
			temp += r[i] * d[i];
		MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		*/

		delta_new = global_T_prod(r, d, lev);

		delta0 = delta_new;

		cycle = 0;
		while(cycle < max_cycle && delta_new > accu * delta0)
		{
			K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, d, q, n, n, lev, data_type);

			/*
			temp=0.0;
			for(i=0; i<n; i++)
				temp += d[i] * q[i];
			MPI_Allreduce(&temp, &temp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			temp1 = global_T_prod(d, q, lev);

			alpha = delta_new / temp1;

			for(i=0; i<n; i++)
				x[i] += alpha * d[i];
			enforce_constrain_Teq(x, lev);

			if(cycle % 50 == 0)
			{
				K_prod_d_heat(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
				for(i=0; i<n; i++)
					r[i] = b[i] - r[i];
			}
			else
			{
				for(i=0; i<n; i++)
					r[i] -= alpha * q[i];
			}

			delta_old = delta_new;

			/*
			temp = 0.0;
			for(i=0; i<n; i++)
				temp += r[i] * r[i];
			MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			delta_new = global_T_prod(r, r, lev);

			beta = delta_new / delta_old;
			for(i=0; i<n; i++)
				d[i] = r[i] + beta * d[i];

			cycle++;

			residual = 0.0;
			for(i=0; i<n; i++)
			{
				if(fabs(alpha * r[i]) > residual)
					residual = fabs(alpha * r[i]);
			}
			MPI_Allreduce(&residual,&global_residual,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			if(global_residual < accu && cycle > 10)
				break;
		}
		if(G.ip == 0)
			printf("Temperature CG cycle = %d\n",cycle);

		free(r);
		free(d);
		free(q);
	}

	return;
}

void other_Teq_contribution(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr_out, *arr_in;
	int m,lev,method;
	char data_type;
{
        int i, j, r, c, l, a;
        double *send, **recv;
        int *nrecv;
        MPI_Status stats;

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                send = malloc(G.smpi[lev].Tsend_number[i] * sizeof(double));
                for(j=0; j<G.smpi[lev].Tsend_number[i]; j++)
                {
                        l = G.smpi[lev].Tsend_row[i][j];
                        r = G.smpi[lev].Tsend_matr_row[i][j];
                        c = G.smpi[lev].Tsend_matr_col[i][j];
                        send[j] = matr_row_val[r][c] * arr_in[l];
                }
                MPI_Send(send, G.smpi[lev].Tsend_number[i], MPI_DOUBLE, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
                free(send);
        }

        recv = malloc(G.smpi[lev].nshareip * sizeof(double *));
        nrecv = malloc(G.smpi[lev].nshareip * sizeof(int));

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_DOUBLE, &nrecv[i]);
                recv[i] = malloc(nrecv[i] * sizeof(double));
                MPI_Recv(recv[i],nrecv[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
        }

        for(i=0; i<G.smpi[lev].nshareip; i++)
        for(j=0; j<G.smpi[lev].Trecv_number[i]; j++)
        {
                l = G.smpi[lev].Trecv_row[i][j];
        //      ip = G.smpi[lev].recv_matr_ip[i][j]; ip should be equal to i
                c = G.smpi[lev].Trecv_matr_col[i][j];

                if(method == 0)
                {
                        a = value_in_sorted_array(l, matr_row_col[l], matr_ncol[l]);
                        arr_out[l] -= recv[i][c]/matr_row_val[l][a];
                }
                else if(method == 1)
                {
                        arr_out[l] += recv[i][c];
                }
                else if(method == 2)
                {
                        arr_out[l] -= recv[i][c];
                }
        }


        for(i=0; i<G.smpi[lev].nshareip; i++)
                free(recv[i]);
        free(recv);
        free(nrecv);
        return;
}




void setup_MPI_comm_Tframe(lev)
	int lev;
{
	int i, j, k, c, l, m;
	int *nseq, nsh, *ip, **seq, *eq_node, *index, *checked, ncheck;
	int *send_row_node, *send_col_node, nsend, size, *send_lg_info;
	int **recv_row_node, **recv_col_node, *nrecv, **recv_lg_info;
	int *tmp_size, *number, inode, node;

	int **matr_row_col, *matr_ncol;
	MPI_Status stats;

	nseq = G.smpi[lev].n_Teq;
	seq = G.smpi[lev].Teq;

	nsh = G.smpi[lev].nshareip;
	ip = G.smpi[lev].ip;

	eq_node = G.heat[lev].eq_node;

	send_row_node = NULL;
	send_col_node = NULL;
	send_lg_info = NULL;
	m = G.heat[lev].neq;

        matr_row_col = G.heat[lev].Krow_col;
        matr_ncol = G.heat[lev].Kncol;

	index = malloc(m * sizeof(int));

        tmp_size = malloc(nsh * sizeof(int));
        number = malloc(nsh * sizeof(int));

        for(i=0; i<nsh; i++)
        {
                number[i] = 0;
                tmp_size[i] = 1;
                G.smpi[lev].Tsend_row[i] = realloc(G.smpi[lev].Tsend_row[i], tmp_size[i] * sizeof(int));
                G.smpi[lev].Tsend_matr_row[i] = realloc(G.smpi[lev].Tsend_matr_row[i], tmp_size[i] * sizeof(int));
                G.smpi[lev].Tsend_matr_col[i] = realloc(G.smpi[lev].Tsend_matr_col[i], tmp_size[i] * sizeof(int));
        }

	for(i=0; i<nsh; i++)
	{
		nsend = 0;
		size = 1;
		send_row_node = malloc(size * sizeof(int));
		send_col_node = malloc(size * sizeof(int));
		send_lg_info = malloc(size * sizeof(int));

		for(j=0; j<m; j++)
			index[j] = value_in_sorted_array(j, seq[i], nseq[i]);

		for(j=0; j<nseq[i]; j++)
		{
			k = seq[i][j];
			for(l=0; l<matr_ncol[k]; l++)
			{
				c = matr_row_col[k][l];
				if(index[c] == -1)
				{
					send_row_node[nsend] = eq_node[k];
					send_col_node[nsend] = eq_node[c];

                                        G.smpi[lev].Tsend_row[i][ number[i] ] = c;
                                        G.smpi[lev].Tsend_matr_row[i][ number[i] ] = k;
                                        G.smpi[lev].Tsend_matr_col[i][ number[i] ] = l;
                                        number[i]++;
                                        if(number[i] >= tmp_size[i])
                                        {
                                                tmp_size[i] *= 2;
                                                G.smpi[lev].Tsend_row[i] = realloc(G.smpi[lev].Tsend_row[i], tmp_size[i] * sizeof(int));
                                                G.smpi[lev].Tsend_matr_row[i] = realloc(G.smpi[lev].Tsend_matr_row[i], tmp_size[i] * sizeof(int));
                                                G.smpi[lev].Tsend_matr_col[i] = realloc(G.smpi[lev].Tsend_matr_col[i], tmp_size[i] * sizeof(int));
                                        }
					if(value_in_array(eq_node[c], G.lg_hangle[lev].node, G.lg_hangle[lev].nh) >= 0)
                                                send_lg_info[nsend] = 1;
                                        else
                                                send_lg_info[nsend] = 0;

					nsend++;
					if(nsend >= size)
					{
						size *= 2;
						send_row_node = realloc(send_row_node, size * sizeof(int));
						send_col_node = realloc(send_col_node, size * sizeof(int));
						send_lg_info = realloc(send_lg_info, size * sizeof(int));
					}
				}
			}
		}

		send_row_node = realloc(send_row_node, nsend * sizeof(int));
		send_col_node = realloc(send_col_node, nsend * sizeof(int));
		send_lg_info = realloc(send_lg_info, nsend * sizeof(int));

		MPI_Send(send_row_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(send_col_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(send_lg_info, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		free(send_row_node);
		free(send_col_node);
		free(send_lg_info);
	}
        for(i=0; i<nsh; i++)
        {
                G.smpi[lev].Tsend_row[i] = realloc(G.smpi[lev].Tsend_row[i], number[i] * sizeof(int));
                G.smpi[lev].Tsend_matr_row[i] = realloc(G.smpi[lev].Tsend_matr_row[i], number[i] * sizeof(int));
                G.smpi[lev].Tsend_matr_col[i] = realloc(G.smpi[lev].Tsend_matr_col[i], number[i] * sizeof(int));
                G.smpi[lev].Tsend_number[i] = number[i];
        }

	recv_row_node = malloc(nsh * sizeof(int *));
	recv_col_node = malloc(nsh * sizeof(int *));
	recv_lg_info = malloc(nsh * sizeof(int *));
	nrecv = malloc(nsh * sizeof(int));
	for(i=0; i<nsh; i++)
	{
                MPI_Probe(ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                recv_row_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_lg_info[i] = malloc(nrecv[i] * sizeof(int));

                MPI_Recv(recv_row_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_lg_info[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
	}

        for(i=0; i<nsh; i++)
        {
                number[i] = 0;
                tmp_size[i] = 1;
                G.smpi[lev].Trecv_row[i] = realloc(G.smpi[lev].Trecv_row[i], tmp_size[i] * sizeof(int));
                G.smpi[lev].Trecv_matr_ip[i] = realloc(G.smpi[lev].Trecv_matr_ip[i], tmp_size[i] * sizeof(int));
                G.smpi[lev].Trecv_matr_col[i] = realloc(G.smpi[lev].Trecv_matr_col[i], tmp_size[i] * sizeof(int));
        }


	/*



	int *share_eq, nshare_eq;

	nshare_eq = 0;
	size = 1;
	share_eq = malloc(size * sizeof(int));

	for(i=0; i<nsh; i++)
	{
		for(j=0; j<nseq[i]; j++)
		{
			share_eq[ nshare_eq++ ] = seq[i][j];
			if(nshare_eq >= size)
			{
				size *= 2;
				share_eq = realloc(share_eq, size * sizeof(int));
			}
		}
	}
	share_eq = realloc(share_eq, nshare_eq * sizeof(int));
	sort(share_eq, nshare_eq);
	nshare_eq = remove_dup(share_eq, nshare_eq);

	for(i=0; i<nshare_eq; i++)
	{
		ncheck = 0;
		size = 1;
		checked = malloc(size * sizeof(int));

		l = share_eq[i];
		for(j=0; j<nsh; j++)
		{
			for(k=0; k<nrecv[j]; k++)
			{
				if(eq_node[l] == recv_row_node[j][k])
				{
					if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1 || recv_lg_info[j][k] == 1)
					{
						G.smpi[lev].Trecv_row[j][ number[j] ] = l;
						G.smpi[lev].Trecv_matr_ip[j][ number[j] ] = j;
						G.smpi[lev].Trecv_matr_col[j][ number[j] ] = k;
						number[j]++;
						if(number[j] >= tmp_size[j])
						{
							tmp_size[j] *= 2;
							G.smpi[lev].Trecv_row[j] = realloc(G.smpi[lev].Trecv_row[j], tmp_size[j] * sizeof(int));
							G.smpi[lev].Trecv_matr_ip[j] = realloc(G.smpi[lev].Trecv_matr_ip[j], tmp_size[j] * sizeof(int));
							G.smpi[lev].Trecv_matr_col[j] = realloc(G.smpi[lev].Trecv_matr_col[j], tmp_size[j] * sizeof(int));
						}

						checked[ ncheck++ ] = recv_col_node[j][k];
						if(ncheck >= size)
						{
							size *= 2;
							checked = realloc(checked, size * sizeof(int));
						}
					}
				}
			}
		}
		free(checked);
	}
	*/
	int *share_node, nshare_node;

	nshare_node = G.smpi[lev].nshare_node;
	share_node = G.smpi[lev].share_node;
	for(i=0; i<nshare_node; i++)
	{
		node = share_node[i];
		inode = node_index(node, lev);
		l = G.heat[lev].id[inode];
		if(l >= 0)
		{
			 ncheck = 0;
                         size = 1;
                         checked = malloc(size * sizeof(int));
			 for(j=0; j<nsh; j++)
			 {
				for(k=0; k<nrecv[j]; k++)
				{
					if(eq_node[l] == recv_row_node[j][k])
					{
						if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1 || recv_lg_info[j][k] == 1)
						{
							G.smpi[lev].Trecv_row[j][ number[j] ] = l;
							G.smpi[lev].Trecv_matr_ip[j][ number[j] ] = j;
							G.smpi[lev].Trecv_matr_col[j][ number[j] ] = k;
							number[j]++;
							if(number[j] >= tmp_size[j])
							{
								tmp_size[j] *= 2;
								G.smpi[lev].Trecv_row[j] = realloc(G.smpi[lev].Trecv_row[j], tmp_size[j] * sizeof(int));
								G.smpi[lev].Trecv_matr_ip[j] = realloc(G.smpi[lev].Trecv_matr_ip[j], tmp_size[j] * sizeof(int));
								G.smpi[lev].Trecv_matr_col[j] = realloc(G.smpi[lev].Trecv_matr_col[j], tmp_size[j] * sizeof(int));
							}

							checked[ ncheck++ ] = recv_col_node[j][k];
							if(ncheck >= size)
							{
								size *= 2;
								checked = realloc(checked, size * sizeof(int));
							}
						}
					}
				}
			 }
			 free(checked);
		}
	}

        for(i=0; i<nsh; i++)
        {
                G.smpi[lev].Trecv_row[i] = realloc(G.smpi[lev].Trecv_row[i], number[i] * sizeof(int));
                G.smpi[lev].Trecv_matr_ip[i] = realloc(G.smpi[lev].Trecv_matr_ip[i], number[i] * sizeof(int));
                G.smpi[lev].Trecv_matr_col[i] = realloc(G.smpi[lev].Trecv_matr_col[i], number[i] * sizeof(int));
                G.smpi[lev].Trecv_number[i] = number[i];
        }
        free(tmp_size);
        free(number);

	for(i=0; i<nsh; i++)
	{
		free(recv_row_node[i]);
		free(recv_col_node[i]);
		free(recv_lg_info[i]);
	}
	free(recv_row_node);
	free(recv_col_node);
	free(recv_lg_info);
	free(index);

	return;

}



void other_lg_Teq_contribution(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
        double **matr_row_val;
        int **matr_row_col, *matr_ncol;
        double *arr_out, *arr_in;
        int m,lev,method;
        char data_type;
{
        int i, j, c, d, k, hnode, inode, eq, lg_eq;
        int *isend, number, **irecv, *nrecv;
        double *dsend, **drecv;
        MPI_Status stats;

        //each pair of isend contains 5 components, including 2 cnode, 1 hnode, 1 hnode identification, and 1 dof info
        isend = malloc(G.heat[lev].neq_lr * 5 * sizeof(int));
        dsend = malloc(G.heat[lev].neq_lr * 5 * sizeof(double));

        for(i=0; i<G.heat[lev].neq_lr * 5; i++)
        {
                isend[i] = -1;
                dsend[i] = 0.0;
        }

        number = 0;
        for(i=G.heat[lev].neq_m; i<m; i++)
        {
                for(j=0; j<matr_ncol[i]; j++)
                {
                        c = matr_row_col[i][j];
                        isend[number * 5 + j] = G.heat[lev].eq_node[c];
                        dsend[number * 5 + j] = matr_row_val[i][j] * arr_in[c];
                }
                isend[number*5 + 3] = G.heat[lev].eq_node[i];
                isend[number*5 + 4] = 0;//only 1 dimension for T eqs,
                number++;
        }

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                MPI_Send(isend, number*5, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(dsend, number*5, MPI_DOUBLE, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
        }

        free(isend);
        free(dsend);


        nrecv = malloc(G.smpi[lev].nshareip * sizeof(int));
        irecv = malloc(G.smpi[lev].nshareip * sizeof(int *));
        drecv = malloc(G.smpi[lev].nshareip * sizeof(double *));

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                irecv[i] = malloc(nrecv[i] * sizeof(int));
                drecv[i] = malloc(nrecv[i] * sizeof(double));

                MPI_Recv(irecv[i],nrecv[i],MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(drecv[i],nrecv[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
        }

        for(i=G.heat[lev].neq_m; i<m; i++)
        {
                hnode = G.heat[lev].eq_node[i];

                if(matr_ncol[i] == 3)//this eq already has all 3 columns, not need to do anything
                        continue;

                for(j=0; j<G.smpi[lev].nshareip; j++)
                {
                        for(k=0; k<nrecv[j]/5; k++)
                        {
                                if(hnode == irecv[j][k*5+3])
                                {
                                        for(c=0; c<3; c++)
                                        {
                                                if(node_index(irecv[j][k*5+c],lev) == -1)
                                                {
                                                        arr_out[i] += drecv[j][k*5+c];
                                                }
                                        }
                                }
                        }
                }
        }

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                free(irecv[i]);
                free(drecv[i]);
        }
        free(nrecv);
        free(irecv);
        free(drecv);
        return;
}



void other_Teq_contribution_slow(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
	double **matr_row_val;
        int **matr_row_col, *matr_ncol;
        double *arr_out, *arr_in;
        int m,lev,method;
        char data_type;
{
	int i, j, k, c, l;
	int *nseq, nsh, *ip, **seq, *eq_node, *index, *checked, ncheck;
	int *send_row_node, *send_col_node, nsend, size, *send_lg_info;
	int **recv_row_node, **recv_col_node, *nrecv, **recv_lg_info;
	int node, inode;
	double *send_prod;
	double **recv_prod;

	MPI_Status stats;

	nseq = G.smpi[lev].n_Teq;
	seq = G.smpi[lev].Teq;

	nsh = G.smpi[lev].nshareip;
	ip = G.smpi[lev].ip;

	eq_node = G.heat[lev].eq_node;

	send_row_node = NULL;
	send_col_node = NULL;
	send_lg_info = NULL;
	send_prod = NULL;

	index = malloc(m * sizeof(int));

	for(i=0; i<nsh; i++)
	{
		nsend = 0;
		size = 1;
		send_row_node = malloc(size * sizeof(int));
		send_col_node = malloc(size * sizeof(int));
		send_lg_info = malloc(size * sizeof(int));
		send_prod = malloc(size * sizeof(double));

		for(j=0; j<m; j++)
			index[j] = value_in_sorted_array(j, seq[i], nseq[i]);

		for(j=0; j<nseq[i]; j++)
		{
			k = seq[i][j];
			for(l=0; l<matr_ncol[k]; l++)
			{
				c = matr_row_col[k][l];

				if(index[c] == -1)
				{

					send_row_node[nsend] = eq_node[k];
					send_col_node[nsend] = eq_node[c];
                                      	send_prod[nsend] = matr_row_val[k][l] * arr_in[c]; 
					if(value_in_array(eq_node[c], G.lg_hangle[lev].node, G.lg_hangle[lev].nh) >= 0)
                                                send_lg_info[nsend] = 1;
                                        else
                                                send_lg_info[nsend] = 0;

					nsend++;
					if(nsend >= size)
					{
						size *= 2;
						send_row_node = realloc(send_row_node, size * sizeof(int));
						send_col_node = realloc(send_col_node, size * sizeof(int));
						send_lg_info = realloc(send_lg_info, size * sizeof(int));
						send_prod = realloc(send_prod, size * sizeof(double));
					}
				}
			}
		}

		send_row_node = realloc(send_row_node, nsend * sizeof(int));
		send_col_node = realloc(send_col_node, nsend * sizeof(int));
		send_lg_info = realloc(send_lg_info, nsend * sizeof(int));
		send_prod = realloc(send_prod, nsend * sizeof(double));

		MPI_Send(send_row_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(send_col_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(send_lg_info, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(send_prod, nsend, MPI_DOUBLE, ip[i], 0, MPI_COMM_WORLD);

		free(send_row_node);
		free(send_col_node);
		free(send_lg_info);
		free(send_prod);
	}

	recv_row_node = malloc(nsh * sizeof(int *));
	recv_col_node = malloc(nsh * sizeof(int *));
	recv_lg_info = malloc(nsh * sizeof(int *));
	recv_prod = malloc(nsh * sizeof(double *));
	nrecv = malloc(nsh * sizeof(int));
	for(i=0; i<nsh; i++)
	{
                MPI_Probe(ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);

                recv_row_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_lg_info[i] = malloc(nrecv[i] * sizeof(int));
                recv_prod[i] = malloc(nrecv[i] * sizeof(double));

                MPI_Recv(recv_row_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_lg_info[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_prod[i],nrecv[i],MPI_DOUBLE,ip[i],0,MPI_COMM_WORLD,&stats);
	}

	int *share_node, nshare_node;

	nshare_node = G.smpi[lev].nshare_node;
	share_node = G.smpi[lev].share_node;
	for(i=0; i<nshare_node; i++)
	{
		node = share_node[i];
		inode = node_index(node, lev);
		l = G.heat[lev].id[inode];
		if(l >= 0)
		{
			 ncheck = 0;
                         size = 1;
                         checked = malloc(size * sizeof(int));
			 for(j=0; j<nsh; j++)
			 {
				for(k=0; k<nrecv[j]; k++)
				{
					if(eq_node[l] == recv_row_node[j][k])
					{
						if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1 || recv_lg_info[j][k] == 1)
						{
							if(method == 0)
							{
								c = value_in_sorted_array(l, matr_row_col[l], matr_ncol[l]);
                                                                arr_out[l] -= recv_prod[j][k]/matr_row_val[l][c];
							}
							else if(method == 1)
							{
								arr_out[l] += recv_prod[j][k];
							}

							checked[ ncheck++ ] = recv_col_node[j][k];
                                                        if(ncheck >= size)
                                                        {
                                                                size *= 2;
                                                                checked = realloc(checked, size * sizeof(int));
                                                        }
						}
					}
				}
			 }
			 free(checked);
		}
	}


	/*
	int *share_eq, nshare_eq;

	nshare_eq = 0;
	size = 1;
	share_eq = malloc(size * sizeof(int));

	for(i=0; i<nsh; i++)
	{
		for(j=0; j<nseq[i]; j++)
		{
			share_eq[ nshare_eq++ ] = seq[i][j];
			if(nshare_eq >= size)
			{
				size *= 2;
				share_eq = realloc(share_eq, size * sizeof(int));
			}
		}
	}
	share_eq = realloc(share_eq, nshare_eq * sizeof(int));
	sort(share_eq, nshare_eq);
	nshare_eq = remove_dup(share_eq, nshare_eq);

	for(i=0; i<nshare_eq; i++)
	{
		ncheck = 0;
		size = 1;
		checked = malloc(size * sizeof(int));

		l = share_eq[i];
		for(j=0; j<nsh; j++)
		{
			for(k=0; k<nrecv[j]; k++)
			{
				if(eq_node[l] == recv_row_node[j][k])
				{
					if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1 || recv_lg_info[j][k] == 1)
					{
						arr_out[l] += recv_prod[j][k];
						checked[ ncheck++ ] = recv_col_node[j][k];
						if(ncheck >= size)
						{
							size *= 2;
							checked = realloc(checked, size * sizeof(int));
						}
					}
				}
			}
		}
		free(checked);
	}
	*/

	for(i=0; i<nsh; i++)
	{
		free(recv_row_node[i]);
		free(recv_col_node[i]);
		free(recv_lg_info[i]);
		free(recv_prod[i]);
	}
	free(recv_row_node);
	free(recv_col_node);
	free(recv_lg_info);
	free(recv_prod);
	free(index);

	return;

}
