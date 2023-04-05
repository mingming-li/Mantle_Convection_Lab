#include "global_variables.h"


void velocity_pressure_KGF_matrix()
{
	int lev, i;
	
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.F[lev] = realloc(G.F[lev], G.neq[lev] * sizeof(double));
		for(i=0; i<G.neq_m[lev]; i++)
			G.F[lev][i]=0.0;

		K_matrix(lev);
		F_matrix(lev);
		G_matrix(lev);
		remote_KF_contribution(lev);

		if(G.lg_mul)
			build_lagrange_K_matrix(lev);

		if(G.control.AMR == 1 || G.solution_cycle == 0)
		{
			setup_MPI_comm_vframe(lev);
		}
	}
	/*
        lev = G.max_level - 1;
        for(i=0; i<G.neq[lev]; i++)
        for(int j=0; j<G.Kncol[lev][i]; j++)
	{
		int k = G.Krow_col[lev][i][j];
                printf("%d %d %d %e\n",G.ip, i,k, G.Krow_val[lev][i][j]);
	}
	*/



	return;
}

void K_matrix(lev)
	int lev;
{
	int e,i,a,j,b,P,p,q,Q,node,k, kk;
	int index_arr[VPTS];
	double K[NEE][NEE];

	for(i=0; i<G.neq_m[lev]; i++)
	{
		G.Kncol[lev][i] = 0;
		for(j=0; j<MAX_V_NCOL; j++)
			G.Krow_val[lev][i][j] = 0.0;
	}

	for(e=0; e<G.nel[lev]; e++)
	{
		for(b=0; b<VPTS; b++)
               		index_arr[b] = node_index(G.ien[lev][b][e],lev);

		get_element_K(e,K,lev);

                for(i=0; i<NDF; i++)
                for(a=0; a<NEN; a++)
                {
                        P = G.lm[lev][i][a][e];
                        p = NED * a + i;
                        for(j=0; j<NDF; j++)
                        for(b=0; b<NEN; b++)
                        {
                                Q = G.lm[lev][j][b][e];
				q = NED * b + j;

				//if(P <= Q && P >= 0 && Q >= 0) disregard K sym for simplicity
				if(P >= 0 && Q >= 0)
				{
					k = value_in_array(Q, G.Krow_col[lev][P], G.Kncol[lev][P]);
					if(k >= 0)
					{
						G.Krow_val[lev][P][k] += K[p][q];
					}
					else
					{
						G.Krow_col[lev][P][ G.Kncol[lev][P] ] = Q;
						G.Krow_val[lev][P][ G.Kncol[lev][P] ] += K[p][q];
						G.Kncol[lev][P]++;
					}
				}
					
				if(P >= 0)
				{
					node = index_arr[b];
					if(G.bc.v_bc[lev][j][node] == 'V')
					{
						G.F[lev][P] -= K[p][q] * G.bc.v_bc_val[lev][j][node];
					}
				}
                        }
                }
	}




	/*
	for(i=0; i<G.neq[lev]; i++)
	{
		G.Krow_col[lev][i] = realloc(G.Krow_col[lev][i], G.Kncol[lev][i] * sizeof(int));
		G.Krow_val[lev][i] = realloc(G.Krow_val[lev][i], G.Kncol[lev][i] * sizeof(double));
	}
	*/


	return;
}

void F_matrix(lev)
	int lev;
{
        int a,e,P,i,p;
        double F[NEE];

        for(e=0; e<G.nel[lev]; e++)
        {
                get_element_F(e,F,lev);
                for(i=0; i<NDF; i++)
                for(a=0; a<NEN; a++)
                {
                        P = G.lm[lev][i][a][e];
                        if(P >= 0)
			{
				p = NED * a + i;
				G.F[lev][P] += F[p];
			}
                }
        }


	return;
}

void G_matrix(lev)
	int lev;
{
        int e,a,i,P,p,k,j,kk;
        double g[NEE];

        for(i=0; i<G.neq_m[lev]; i++)
        {
                G.Gncol[lev][i] = 0;
                for(j=0; j<MAX_V_NCOL; j++)
                        G.Grow_val[lev][i][j] = 0.0;
        }


        for(e=0; e<G.nel[lev]; e++)
        {
                get_element_G(e,g,lev);
                for(i=0; i<NDF; i++)
                for(a=0; a<NEN; a++)
                {
                        P = G.lm[lev][i][a][e];
                        p = NED * a + i;
			if(P == -1)
				continue;

			k = value_in_array(e, G.Grow_col[lev][P], G.Gncol[lev][P]);
			if(k >= 0)
			{
				G.Grow_val[lev][P][k] += g[p];
			}
			else
			{
				G.Grow_col[lev][P][ G.Gncol[lev][P] ] = e;
				G.Grow_val[lev][P][ G.Gncol[lev][P] ] += g[p];
				G.Gncol[lev][P]++;
			}
                }
        }

	/*
	for(i=0; i<G.neq[lev]; i++)
	{
		G.Grow_col[lev][i] = realloc(G.Grow_col[lev][i], G.Gncol[lev][i] * sizeof(int));
		G.Grow_val[lev][i] = realloc(G.Grow_val[lev][i], G.Gncol[lev][i] * sizeof(double));
	}
	*/

	return;
}

void get_element_K(e,K,lev)
	int e,lev;
	double K[NEE][NEE];
{
        int i,j,a,b,p,q,g,d,node[VPTS];
        double eta;
        double D[3][3]={{0.0,0.0,0.0},
                        {0.0,0.0,0.0},
                        {0.0,0.0,0.0}};
        double dNa_dx,dNa_dz;
        double dNb_dx,dNb_dz;


	for(d=0; d<NEN; d++)
	{
		node[d] = node_index(G.ien[lev][d][e],lev);
	}

        for(a=0; a<NEN; a++)
        for(i=0; i<NED; i++)
	{
		p = NED * a + i;
		for(b=a; b<NEN; b++)
		for(j=0; j<NED; j++)
		{
			q = NED * b + j;
			K[p][q] = 0.0;
			for(g=0; g<NEG; g++)
			{
				eta = 0.0;
				for(d=0; d<NEN; d++)
				{
					eta += G.N[d].gauss[g] * G.vis[lev][node[d]];
				}
			        D[0][0]=D[1][1]=2.0*eta;
			        D[2][2]=eta;


                                dNa_dx=G.dN_dxi[ a].gauss[g]*G.dz_deta[lev][e]
                                      -G.dN_deta[a].gauss[g]*G.dz_dxi;             
                                dNa_dz=G.dN_deta[a].gauss[g]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ a].gauss[g]*G.dx_deta;

                                dNb_dx=G.dN_dxi[ b].gauss[g]*G.dz_deta[lev][e]
                                      -G.dN_deta[b].gauss[g]*G.dz_dxi;             
                                dNb_dz=G.dN_deta[b].gauss[g]*G.dx_dxi[lev][e]
                                      -G.dN_dxi[ b].gauss[g]*G.dx_deta;

                                if(i==0 && j==0)
                                {
                                        K[p][q]+=((dNa_dx*D[0][0]+dNa_dz*D[2][0])*dNb_dx
                                                 +(dNa_dx*D[0][2]+dNa_dz*D[2][2])*dNb_dz);
                                }
                                else if(i==1 && j==0)
                                {
                                        K[p][q]+=((dNa_dz*D[1][0]+dNa_dx*D[2][0])*dNb_dx
                                                 +(dNa_dz*D[1][2]+dNa_dx*D[2][2])*dNb_dz);
                                }
                                else if(i==0 && j==1)
                                {
                                        K[p][q]+=((dNa_dx*D[0][1]+dNa_dz*D[2][1])*dNb_dz
                                                 +(dNa_dx*D[0][2]+dNa_dz*D[2][2])*dNb_dx);
                                }
                                else if(i==1 && j==1)
                                {
                                        K[p][q]+=((dNa_dz*D[1][1]+dNa_dx*D[2][1])*dNb_dz
                                                 +(dNa_dz*D[1][2]+dNa_dx*D[2][2])*dNb_dx);
                                }

			}
			K[p][q]/=G.jacobi[lev][e];
		}
	}
        for(p=0; p<NEE; p++)
        for(q=0; q<NEE; q++)
        {
                K[q][p] = K[p][q];
        }
	return;
}

void get_element_F(e,F,lev)
	int e,lev;
	double F[NEE];
{
        int a,d,node,i,p,index_arr[VPTS];
        double force_at_gs[4];
        int natural_bc[2];
        double L;


        for(a=0; a<NEN; a++)
		index_arr[a] = node_index(G.ien[lev][a][e],lev);

        for(d=0; d<NEG; d++)
        {
                force_at_gs[d] = 0.0;
                for(a=0; a<NEN; a++)
                {
                        force_at_gs[d] += G.N[a].gauss[d] * G.buoyancy[lev][index_arr[a]];
                }
        }

        natural_bc[0] = natural_bc[1] = 0;
        for(i=0; i<NED; i++)
        {
                for(a=0; a<NEN; a++)
                {
                        node = index_arr[a];
                        if(G.bc.v_bc[lev][i][node] == 'S')
                        {
                                natural_bc[i] += a;
                        }
                }
                if(natural_bc[i]==3)//horizontal bc
                {
                        natural_bc[i]=1;
                }
                else if(natural_bc[i]==1||natural_bc[i]==5)//vertical bc
                {
                        natural_bc[i]=2;
                }
                else if(natural_bc[i]==0 || natural_bc[i]==1
                     || natural_bc[i]==2 || natural_bc[i]==3)//ignore corner of domain
                {
                        natural_bc[i]=0;
                }
        }

        for(i=0; i<NED; i++)
        for(a=0; a<NEN; a++)
        {
                p = NED * a + i;
                F[p]=0.0;
                if(i == 1)
                {
                        for(d=0; d<NEG; d++)
                        {
                                F[p] += G.N[a].gauss[d] * force_at_gs[d];
                        }
                        F[p] *= G.jacobi[lev][e];
                }

                if(natural_bc[i] != 0)
                {
                        node = index_arr[a];
                        L = G.X[lev][i][index_arr[2]] - G.X[lev][i][index_arr[0]];
                        F[p] += 0.5 * L * G.bc.v_bc_val[lev][i][node];
                }
        }

	return;
}

void get_element_G(e,g,lev)
	int e,lev;
	double g[NEE];
{
        int i,a,p,gs;
	double dNa_dx,dNa_dz;


        i=0;
        for(a=0; a<NEN; a++)
        {
                p = NED * a + i;
                g[p] = 0.0;
                for(gs=0; gs<NEG; gs++)
                {
                        dNa_dx=G.dN_dxi[ a].gauss[gs]*G.dz_deta[lev][e]
                              -G.dN_deta[a].gauss[gs]*G.dz_dxi;             
                        g[p] -= dNa_dx;
                }
        }

        i=1;
        for(a=0; a<NEN; a++)
        {
                p = NED * a + i;
                g[p] = 0.0;
                for(gs=0; gs<NEG; gs++)
                {
                        dNa_dz=G.dN_deta[a].gauss[gs]*G.dx_dxi[lev][e]
                              -G.dN_dxi[ a].gauss[gs]*G.dx_deta;
                        g[p] -= dNa_dz;
                }
        }

	return;
}

void remote_KF_contribution(lev)
	int lev;
{
        int i,j,k,p,n,q,c;
        MPI_Status stats;

	int size;

	if(G.control.mpi_KF_method == 1)//memory good, time bad
	{
	        double *send_K,*recv_K;
       	 	double *send_F,*recv_F;
		//for F
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			send_F = malloc(G.smpi[lev].n_veq[i] * sizeof(double));

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				send_F[j]=G.F[lev][p];
			}
			MPI_Send(send_F,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
			free(send_F);
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			recv_F = malloc(G.smpi[lev].n_veq[i] * sizeof(double));
			MPI_Recv(recv_F,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				G.F[lev][p] += recv_F[j];
			}
			free(recv_F);
		}

		//send K
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				send_K = malloc(G.smpi[lev].n_veq[i] * sizeof(double));
				n = 0;
				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);

					if(c >= 0)
						send_K[n] = G.Krow_val[lev][p][c];
					else
						send_K[n] = 0.0;
					n++;
				}
				//send for each row, this requires much less memory, but more MPI communication
				MPI_Send(send_K,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
				free(send_K);
			}
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				recv_K = malloc(G.smpi[lev].n_veq[i] * sizeof(double));
				
				//recv for each row, this requires much less memory, but more MPI communication
				MPI_Recv(recv_K,G.smpi[lev].n_veq[i]*G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

				p = G.smpi[lev].veq[i][j];

				size = G.Kncol[lev][p];
				 
				n=0;
				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);

					if (c >= 0)
					{
						G.Krow_val[lev][p][c] += recv_K[n];
					}
					else
					{
						if(fabs(recv_K[n])>0.0)
						{
							c = G.Kncol[lev][p]++;

							/*
							if(G.Kncol[lev][p] >= size)
							{
								size *= 2;
								G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], size * sizeof(double));
								G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], size * sizeof(int));
							}
							*/

							G.Krow_val[lev][p][c] = recv_K[n];
							G.Krow_col[lev][p][c] = q;
						}
					}
					n++;
				}

				/*
				G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], G.Kncol[lev][p] * sizeof(double));
				G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], G.Kncol[lev][p] * sizeof(int));
				*/

				free(recv_K);
			}
		}
	}
	else if(G.control.mpi_KF_method == 2)//time good, memory bad, default
	{
		/* working and only using 1 MPI, but require large memory */
		double *send, *recv;
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			n = G.smpi[lev].n_veq[i] + G.smpi[lev].n_veq[i] * G.smpi[lev].n_veq[i];// the first half is for F, and the other half after + is for K

			send = malloc(n * sizeof(double));

			n = 0;

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				send[n] = G.F[lev][p];
				n++;
			}

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);
					// c can not -1, meaning that this value does not exist in the matrix. 
					// however, for simplicty of matching index between send and recv, i still send values when they are zero

					if(c >= 0)
						send[n] = G.Krow_val[lev][p][c];
					else
						send[n] = 0.0;
					n++;
				}
			}

			MPI_Send(send,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);

			free(send);
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			n = G.smpi[lev].n_veq[i] + G.smpi[lev].n_veq[i] * G.smpi[lev].n_veq[i];

			recv = malloc(n * sizeof(double));

			MPI_Recv(recv,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

			n = 0;
			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				G.F[lev][p] += recv[n];
				n++;
			}

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];

				size = G.Kncol[lev][p];

				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);

					if (c >= 0)//the index already exist, no need to increase the array size
					{
						G.Krow_val[lev][p][c] += recv[n];
					}
					else//if the index does not exist locally, but is calcualted remotely, need to increase the size of Kncol array
					{
						if(fabs(recv[n])>0.0)// only add new number when its mag is larger than zero
						{
							
							c = G.Kncol[lev][p]++;
							/*

							if(G.Kncol[lev][p] >= size)
							{
								size *= 2;
								G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], size * sizeof(double));
								G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], size * sizeof(int));
							}
							*/

							G.Krow_val[lev][p][c] = recv[n];
							G.Krow_col[lev][p][c] = q;
						}
					}

					n++;
				}

				/*
				G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], G.Kncol[lev][p] * sizeof(double));
				G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], G.Kncol[lev][p] * sizeof(int));
				*/
			}

			free(recv);
		}
	}
	else if(G.control.mpi_KF_method == 3)//time and memory between 1 and 2, memory slightly better than 2, but still fast
	{
		double *send_K,*recv_K;
		double *send_F,*recv_F;
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			send_F = malloc(G.smpi[lev].n_veq[i] * sizeof(double));

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				send_F[j]=G.F[lev][p];
			}
			MPI_Send(send_F,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
			free(send_F);
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			recv_F = malloc(G.smpi[lev].n_veq[i] * sizeof(double));
			MPI_Recv(recv_F,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				G.F[lev][p] += recv_F[j];
			}
			free(recv_F);
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			send_K = malloc(G.smpi[lev].n_veq[i] * G.smpi[lev].n_veq[i] * sizeof(double));

			n = 0;
			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);
					if(c >= 0)
						send_K[n] = G.Krow_val[lev][p][c];
					else
						send_K[n] = 0.0;
					n++;
				}
			}
			MPI_Send(send_K,n,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
			free(send_K);
		}

		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			recv_K = malloc(G.smpi[lev].n_veq[i] * G.smpi[lev].n_veq[i] * sizeof(double));

			MPI_Recv(recv_K,G.smpi[lev].n_veq[i]*G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

			n=0;
			for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			{
				p = G.smpi[lev].veq[i][j];
				size = G.Kncol[lev][p];
				for(k=0; k<G.smpi[lev].n_veq[i]; k++)
				{
					q = G.smpi[lev].veq[i][k];

					c = value_in_array(q, G.Krow_col[lev][p], G.Kncol[lev][p]);

					if (c >= 0)
					{
						G.Krow_val[lev][p][c] += recv_K[n];
					}
					else
					{
						if(fabs(recv_K[n])>0.0)
						{
							c = G.Kncol[lev][p]++;

							/*
							if(G.Kncol[lev][p] >= size)
							{
								size *= 2;
								G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], size * sizeof(double));
								G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], size * sizeof(int));
							}
							*/

							G.Krow_val[lev][p][c] = recv_K[n];
							G.Krow_col[lev][p][c] = q;
						}
					}


					n++;
				}
				/*
				G.Krow_val[lev][p] = realloc(G.Krow_val[lev][p], G.Kncol[lev][p] * sizeof(double));
				G.Krow_col[lev][p] = realloc(G.Krow_col[lev][p], G.Kncol[lev][p] * sizeof(int));
				*/
			}

			free(recv_K);
		}
	
	}

	//sort col for each row, this is a must
	int log;
	for(i=0; i<G.neq[lev]; i++)
	{
		log = 1;
		while(log)
		{
			log = 0;
			for(j=0; j<G.Kncol[lev][i]-1; j++)
			{
				if(G.Krow_col[lev][i][j] > G.Krow_col[lev][i][j+1])
				{
					swap(&G.Krow_col[lev][i][j], &G.Krow_col[lev][i][j+1]);
					swap_d(&G.Krow_val[lev][i][j], &G.Krow_val[lev][i][j+1]);
					log = 1;
				}
			}
		}
	}

	for(i=0; i<G.neq[lev]; i++)
	{
		k = value_in_array(i, G.Krow_col[lev][i], G.Kncol[lev][i]);
		G.Kdiag[lev][i] = G.Krow_val[lev][i][k];
	}

	return;
}

