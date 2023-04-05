#include "global_variables.h"

void initial_data_structure()
{
	int lev, d, a;

	G.max_elx = malloc(G.max_level*sizeof(int));
	G.max_elz = malloc(G.max_level*sizeof(int));
	G.max_nel = malloc(G.max_level*sizeof(int));
	G.max_nox = malloc(G.max_level*sizeof(int));
	G.max_noz = malloc(G.max_level*sizeof(int));
	G.max_nno = malloc(G.max_level*sizeof(int));

	G.nel = malloc(G.max_level*sizeof(int));
	G.nno = malloc(G.max_level*sizeof(int));
	G.nno_prime = malloc(G.max_level*sizeof(int));

	G.node = malloc(G.max_level*sizeof(int *));
	G.hash_node = malloc(G.max_level*sizeof(struct HASH_ITEM *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.node[lev] = NULL;
		G.hash_node[lev] = malloc(sizeof(struct HASH_ITEM *));
	}

	G.Mcode = malloc(G.max_level*sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.Mcode[lev] = NULL;

	G.ien = malloc(G.max_level*sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.ien[lev] = malloc(NEN * sizeof(int *));
		for(d=0; d<NEN; d++)
			G.ien[lev][d] = NULL;
	}

	G.global.min_Mcode = malloc(G.max_level*sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.global.min_Mcode[lev] = malloc(G.np * sizeof(int));

	G.global.max_Mcode = malloc(G.max_level*sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.global.max_Mcode[lev] = malloc(G.np * sizeof(int));

	G.smpi = malloc(G.max_level*sizeof(struct MPI_SHARE));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.smpi[lev].veq_nsh = NULL;
		G.smpi[lev].Teq_nsh = NULL;
		for(a=0; a<MAX_SHARE_IP; a++)
		{
			G.smpi[lev].veq[a] = NULL;
			G.smpi[lev].Teq[a] = NULL;

			G.smpi[lev].node[a] = NULL;
			G.smpi[lev].nel[a] = NULL;
			G.smpi[lev].el[a] = malloc(VPTS * sizeof(int *));
			for(d=0; d<VPTS; d++)
				G.smpi[lev].el[a][d] = NULL;


			G.smpi[lev].send_row[a] = NULL;
			G.smpi[lev].send_matr_row[a] = NULL;
			G.smpi[lev].send_matr_col[a] = NULL;

			G.smpi[lev].recv_row[a] = NULL;
			G.smpi[lev].recv_matr_ip[a] = NULL;
			G.smpi[lev].recv_matr_col[a] = NULL;

			G.smpi[lev].Tsend_row[a] = NULL;
			G.smpi[lev].Tsend_matr_row[a] = NULL;
			G.smpi[lev].Tsend_matr_col[a] = NULL;

			G.smpi[lev].Trecv_row[a] = NULL;
			G.smpi[lev].Trecv_matr_ip[a] = NULL;
			G.smpi[lev].Trecv_matr_col[a] = NULL;
		}
		G.smpi[lev].share_node = NULL;
		G.smpi[lev].Tshare_node = NULL;
	}


	G.bc.v_bc = malloc(G.max_level * sizeof(char *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.bc.v_bc[lev] = malloc(NSD * sizeof(char *));
		for(d=0; d<NSD; d++)
			G.bc.v_bc[lev][d] = NULL;
	}

	G.p2n = malloc(G.max_level * sizeof(struct P2N));
	G.eco = malloc(G.max_level * sizeof(struct ECO));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.p2n[lev].nel = NULL;
		for(a=0; a<NEN; a++)
			G.p2n[lev].el[a] = NULL;

		G.eco[lev].area = NULL;
	}


	G.bc.v_bc_val = malloc(G.max_level * sizeof(double *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.bc.v_bc_val[lev] = malloc(NSD * sizeof(double *));
		for(d=0; d<NSD; d++)
			G.bc.v_bc_val[lev][d] = NULL;
	}

	G.id = malloc(G.max_level * sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.id[lev] = malloc(NSD * sizeof(int *));
		for(d=0; d<NSD; d++)
			G.id[lev][d] = NULL;
	}

	G.veq_node = malloc(G.max_level * sizeof(int *));
	G.veq_dof = malloc(G.max_level * sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.veq_node[lev] = NULL;
		G.veq_dof[lev] = NULL;
	}

	G.bc.T_bc = malloc(G.max_level * sizeof(char *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.bc.T_bc[lev] = NULL;

	G.bc.T_bc_val = malloc(G.max_level * sizeof(double *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.bc.T_bc_val[lev] = NULL;

	G.neq = malloc(G.max_level * sizeof(int));
	G.neq_f = malloc(G.max_level * sizeof(int));
	G.neq_c = malloc(G.max_level * sizeof(int));
	G.neq_m = malloc(G.max_level * sizeof(int));
	G.neq_lr = malloc(G.max_level * sizeof(int));

	G.lm = malloc(G.max_level * sizeof(int *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.lm[lev] = malloc(NSD * sizeof(int *));
		for(d=0; d<NSD; d++)
		{
			G.lm[lev][d] = malloc(NEN * sizeof(int *));
			for(a=0; a<NEN; a++)
			{
				G.lm[lev][d][a] = NULL;
			}
		}
	}

	G.dx_dxi = malloc(G.max_level * sizeof(double *));
	G.dz_deta = malloc(G.max_level * sizeof(double *));
	G.jacobi = malloc(G.max_level * sizeof(double *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.dx_dxi[lev] = NULL;
		G.dz_deta[lev] = NULL;
		G.jacobi[lev] = NULL;
	}

	G.dx = malloc(G.max_level * sizeof(double));
	G.dz = malloc(G.max_level * sizeof(double));

	G.X = malloc(G.max_level * sizeof(double *));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.X[lev] = malloc(NSD * sizeof(double *));
		for(d=0; d<NSD; d++)
			G.X[lev][d] = NULL;
	}

	G.T = malloc(G.max_level * sizeof(double *));
	G.buoyancy = malloc(G.max_level * sizeof(double *));
	G.vis = malloc(G.max_level * sizeof(double *));
	G.V = malloc(G.max_level * sizeof(double *));

	G.d = malloc(G.max_level * sizeof(double *));
	G.p = malloc(G.max_level * sizeof(double *));
	G.err = malloc(G.max_level * sizeof(double *));


	G.stress = malloc(G.max_level * sizeof(double *));
	G.gs = malloc(G.max_level * sizeof(double *));

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.T[lev] = NULL;
		G.buoyancy[lev] = NULL;
		G.vis[lev] = NULL;
		G.V[lev] = malloc(NSD * sizeof(double *));
		for(d=0; d<NSD; d++)
			G.V[lev][d] = NULL;
		G.d[lev] = NULL;
		G.p[lev] = NULL;
		G.err[lev] = NULL;
		G.stress[lev] = NULL;
		G.gs[lev] = NULL;
	}

	G.Kncol = malloc(G.max_level * sizeof(int *));
	G.Krow_col = malloc(G.max_level * sizeof(int *));
	G.Krow_val = malloc(G.max_level * sizeof(double *));
	G.Kdiag = malloc(G.max_level * sizeof(double *));

	G.Gncol = malloc(G.max_level * sizeof(int *));
	G.Grow_col = malloc(G.max_level * sizeof(int *));
	G.Grow_val = malloc(G.max_level * sizeof(double *));

	G.F = malloc(G.max_level * sizeof(double *));
	for(lev=G.min_level; lev<G.max_level; lev++)
		G.F[lev] = NULL;

	G.hangle = malloc(G.max_level * sizeof(struct HANGLE));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.hangle[lev].node = NULL;
		G.hangle[lev].ni = NULL;

		for(d=0; d<NSD; d++)
		{
			G.hangle[lev].cnode[d] = NULL;
			G.hangle[lev].ci[d] = NULL;
		}
	}

	G.lg_hangle = malloc(G.max_level * sizeof(struct HANGLE));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.lg_hangle[lev].node = NULL;
		G.lg_hangle[lev].ni = NULL;

		for(d=0; d<NSD; d++)
		{
			G.lg_hangle[lev].cnode[d] = NULL;
			G.lg_hangle[lev].ci[d] = NULL;
		}
	}

	G.interp = malloc(G.max_level * sizeof(struct INTERPOLATE));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.interp[lev].local_eq_n = NULL;
		G.interp[lev].node_amount = NULL;
		G.interp[lev].eq_amount = NULL;
		for(d=0; d<VPTS; d++)
			G.interp[lev].local_eq[d] = NULL;

		G.interp[lev].remote_node = malloc(G.np * sizeof(int *));
		G.interp[lev].remote_node_n = malloc(G.np * sizeof(int *));
		G.interp[lev].remote_eq = malloc(G.np * sizeof(int *));

		for(a=0; a<G.np; a++)
		{
			G.interp[lev].remote_node[a] = NULL;
			G.interp[lev].remote_node_n[a] = NULL;
			G.interp[lev].remote_eq[a] = NULL;
		}
		G.interp[lev].ip_send = NULL;
		G.interp[lev].ip_recv = NULL;

		for(d=0; d<VPTS; d++)
			G.interp[lev].local_lv[d] = NULL;

		G.interp[lev].local_u = NULL;
		G.interp[lev].local_nlv = NULL;

		for(a=0; a<MAX_SHARE_IP; a++)
		{
			for(d=0; d<VPTS; d++)
				G.interp[lev].send_lv[a][d] = NULL;

			G.interp[lev].send_u[a] = NULL;
			G.interp[lev].send_nlv[a] = NULL;

			G.interp[lev].recv_l[a] = NULL;
			G.interp[lev].recv_u[a] = NULL;
		}
	}

	G.project = malloc(G.max_level * sizeof(struct INTERPOLATE));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.project[lev].local_eq_n = NULL;
		G.project[lev].node_amount = NULL;
		G.project[lev].eq_amount = NULL;
		for(d=0; d<VPTS; d++)
			G.project[lev].local_eq[d] = NULL;

		G.project[lev].remote_node = malloc(G.np * sizeof(int *));
		G.project[lev].remote_node_n = malloc(G.np * sizeof(int *));
		G.project[lev].remote_eq = malloc(G.np * sizeof(int *));

		for(a=0; a<G.np; a++)
		{
			G.project[lev].remote_node[a] = NULL;
			G.project[lev].remote_node_n[a] = NULL;
			G.project[lev].remote_eq[a] = NULL;
		}
		G.project[lev].ip_send = NULL;
		G.project[lev].ip_recv = NULL;


		G.project[lev].local_l = NULL;
		G.project[lev].local_u = NULL;

		for(a=0; a<MAX_SHARE_IP; a++)
		{
			G.project[lev].send_l[a] = NULL;
			G.project[lev].send_u[a] = NULL;

			G.project[lev].recv_l[a] = NULL;
			G.project[lev].recv_u[a] = NULL;
		}
	}


	G.heat = malloc(G.max_level * sizeof(struct HEAT));
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.heat[lev].id = NULL;
		G.heat[lev].eq_node = NULL;
		G.heat[lev].lm = malloc(NEN * sizeof(int *));
		for(a=0; a<NEN; a++)
			G.heat[lev].lm[a] = NULL;

		for(d=0; d<NSD; d++)
			G.heat[lev].tau[d] = NULL;

		G.heat[lev].d = NULL;
		G.heat[lev].ddot = NULL;

		G.heat[lev].F = NULL;
	}

	G.head = NULL;
	return;
}

void release_memory()
{
        int lev, i;
	//dynamic memory control is more complex for stiffness matrix, 
	//as they are multidimension whose size are often not well defined at the beginning

        for(lev=G.min_level; lev<G.max_level; lev++)
        {
                free(G.Kncol[lev]);
                free(G.Kdiag[lev]);
                for(i=0; i<G.neq[lev]; i++)
                {
                        free(G.Krow_col[lev][i]);
                        free(G.Krow_val[lev][i]);
                }
                free(G.Krow_col[lev]);
                free(G.Krow_val[lev]);

                free(G.Gncol[lev]);
                for(i=0; i<G.neq[lev]; i++)
                {
                        free(G.Grow_col[lev][i]);
                        free(G.Grow_val[lev][i]);
                }
                free(G.Grow_col[lev]);
                free(G.Grow_val[lev]);
        }


        lev = G.max_level - 1;//heat matrix is only assigned on levmax, no mg solution
        for(i=0; i<G.heat[lev].neq; i++)
        {
                free(G.heat[lev].Krow_col[i]);
                free(G.heat[lev].Krow_val[i]);
        }
        free(G.heat[lev].Krow_col);
        free(G.heat[lev].Krow_val);
        free(G.heat[lev].Kncol);


        for(i=0; i<G.heat[lev].neq; i++)
        {
                free(G.heat[lev].Crow_col[i]);
                free(G.heat[lev].Crow_val[i]);
        }
        free(G.heat[lev].Crow_col);
        free(G.heat[lev].Crow_val);
        free(G.heat[lev].Cncol);

        return;
}

void get_time_independent_variables()
{
        int lev,mul;

        G.np = G.npx * G.npz;
        G.ipz = G.ip % G.npz;
        G.ipx = (G.ip - G.ipz) / G.npz;

        G.length = G.global.length / G.npx;
        G.height = G.global.height / G.npz;

        for(lev=0; lev<G.max_level; lev++)
        {
                mul = pow(2, lev);

                G.max_elx[lev] =  mul;
                G.max_elz[lev] =  mul;
                G.max_nel[lev] = mul * mul;

                G.max_nox[lev] = mul + 1;
                G.max_noz[lev] = mul + 1;
                G.max_nno[lev] = G.max_nox[lev] * G.max_noz[lev];

                G.dx[lev] = G.length / G.max_elx[lev];
                G.dz[lev] = G.height / G.max_elz[lev];
        }

        G.gnox = G.max_elx[G.max_level-1] * G.npx + 1;
        G.gnoz = G.max_elz[G.max_level-1] * G.npz + 1;

        return;
}


void node_position()
{
        int lev,d,i,ix,iz;

        for(lev=G.min_level; lev<G.max_level; lev++)
        for(d=0; d<NSD; d++)
        {
                G.X[lev][d] = realloc(G.X[lev][d], G.nno[lev] * sizeof(double));
        }

        for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(i=0; i<G.nno[lev]; i++)
		{
			ix = G.node[lev][i] / G.gnoz;
			iz = G.node[lev][i] - ix*G.gnoz;

			G.X[lev][0][i] = ix * G.dx[G.max_level-1];
			G.X[lev][1][i] = iz * G.dz[G.max_level-1];
		}
	}
        return;
}


void assign_matrix_memory()
{
	int lev, i, neq;
        for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.Kncol[lev] = malloc(G.neq[lev] * sizeof(int));
		G.Kdiag[lev] = malloc(G.neq[lev] * sizeof(double));
		G.Krow_col[lev] = malloc(G.neq[lev] * sizeof(int *));
		G.Krow_val[lev] = malloc(G.neq[lev] * sizeof(double *));

		for(i=0; i<G.neq[lev]; i++)
		{
			G.Krow_col[lev][i] = malloc(MAX_V_NCOL * sizeof(int));
			G.Krow_val[lev][i] = malloc(MAX_V_NCOL * sizeof(double));
		}


		G.Gncol[lev] = malloc(G.neq[lev] * sizeof(int));
		G.Grow_col[lev] = malloc(G.neq[lev] * sizeof(int *));
		G.Grow_val[lev] = malloc(G.neq[lev] * sizeof(double *));

		for(i=0; i<G.neq[lev]; i++)
		{
			G.Grow_col[lev][i] = malloc(MAX_V_NCOL * sizeof(int));
			G.Grow_val[lev][i] = malloc(MAX_V_NCOL * sizeof(double));
		}

	}


	lev = G.max_level - 1;
	neq = G.heat[lev].neq;

        G.heat[lev].Kncol = malloc(neq * sizeof(int));
        G.heat[lev].Krow_col = malloc(neq * sizeof(int *));
        G.heat[lev].Krow_val = malloc(neq * sizeof(double *));

        for(i=0; i<neq; i++)
        {
                G.heat[lev].Krow_col[i] = malloc(MAX_T_NCOL * sizeof(int));
                G.heat[lev].Krow_val[i] = malloc(MAX_T_NCOL * sizeof(double));
        }

        G.heat[lev].Cncol = malloc(neq * sizeof(int));
        G.heat[lev].Crow_col = malloc(neq * sizeof(int *));
        G.heat[lev].Crow_val = malloc(neq * sizeof(double *));

        for(i=0; i<neq; i++)
        {
                G.heat[lev].Crow_col[i] = malloc(MAX_T_NCOL * sizeof(int));
                G.heat[lev].Crow_val[i] = malloc(MAX_T_NCOL * sizeof(double));
        }



	return;
}
