#include "global_variables.h"

void solve_velocity_pressure()
{
	int lev;
	initial_d_p();
	citcom_solver();
	v_to_nodes();
	compute_stress();
//	p_to_nodes();
	return;
}

void citcom_solver()
{
	double *r0,*r1,*r2,*z0,*z1,*s1,*s2,*Ah,*u1;
	double dvelocity,dpressure;
        double r1dotz1,r0dotr0,delta,s2dotAhat,alpha;
	double temp1,temp2,temp3,temp4,temp;
	double *shuffle;
	double *BAB;
	int levmax = G.max_level - 1;
	int i,cycle,c,j;
	int nel = G.nel[levmax];
	int neq = G.neq[levmax];
	int neq_m = G.neq_m[levmax];

	r0 = malloc(nel*sizeof(double));
	r1 = malloc(nel*sizeof(double));
	r2 = malloc(nel*sizeof(double));
	z0 = malloc(nel*sizeof(double));
	z1 = malloc(nel*sizeof(double));
	s1 = malloc(nel*sizeof(double));
	s2 = malloc(nel*sizeof(double));
	Ah = malloc(neq*sizeof(double));
	u1 = malloc(neq*sizeof(double));

	BAB = malloc(nel * sizeof(double));


	if(G.control.auto_accuracy)
	{
		double v_res;
		int gneq;

		MPI_Allreduce(&neq, &gneq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		v_res = sqrt(global_v_prod(G.F[levmax], G.F[levmax], levmax)/gneq);
		G.control.accuracy = G.control.citcom_accu * v_res;
		if(G.ip == 0)
			fprintf(stderr,"accuracy = %e %e\n", v_res, G.control.accuracy);
	}

	//prepare data for precondition
	for(i=0; i<nel; i++)
		BAB[i] = 0.0;

	for(i=0; i<neq_m; i++)
	{
		for(j=0; j<G.Gncol[levmax][i]; j++)
		{
			c = G.Grow_col[levmax][i][j];
			BAB[c] += G.Grow_val[levmax][i][j] * G.Grow_val[levmax][i][j] / G.Kdiag[levmax][i];
		}
	}


	K_prod_d(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], G.d[levmax], u1, neq, neq, levmax, 'v');

	for(i=0; i<neq; i++)
		Ah[i] = 0.0;

	G_prod_p(G.Grow_val[levmax], G.Grow_col[levmax], G.Gncol[levmax], G.p[levmax], Ah, neq_m, nel, levmax);

	for(i=0; i<neq; i++)
	{
		Ah[i] = G.F[levmax][i] - Ah[i] - u1[i];
	}

	if(G.ip == 0)
		fprintf(stderr,"citcom solver\n");

        if(strcmp(G.solver, "jc") == 0)
        	jacobi(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy, 1000);
	else if(strcmp(G.solver, "gs") == 0)
        	gauss_seidel(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy, 1000);
	else if(strcmp(G.solver, "cg") == 0)
                Conj_Grad(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy, G.control.max_cg_cycle, 'v');
	else if(strcmp(G.solver, "mg") == 0)
		multigrid(G.Krow_val, G.Krow_col, G.Kncol, Ah, u1, G.neq, G.control.accuracy);


	/*
	double **K;
	int k;

	K = malloc(neq * sizeof(double *));
	for(i=0; i<neq; i++)
		K[i] = malloc(neq * sizeof(double));

	for(i=0; i<neq; i++)
	for(j=0; j<neq; j++)
		K[i][j] = 0.0;

	for(i=0; i<neq; i++)
	for(j=0; j<G.Kncol[levmax][i]; j++)
	{
		k = G.Krow_col[levmax][i][j];
		K[i][k] = G.Krow_val[levmax][i][j];
	}

	for(i=0; i<neq; i++)
	{
		for(j=0; j<neq - 1; j++)
			printf("%e ",K[i][j]);
		printf("%e\n",K[i][j]);
	}

	terminate();
	*/

	for(i=0; i<neq; i++)
	{
		G.d[levmax][i] += u1[i];
	}

	GT_prod_d(G.Grow_val[levmax], G.Grow_col[levmax], G.Gncol[levmax], G.d[levmax], r1, neq_m, nel, levmax);

        dvelocity=1.0;
        dpressure=1.0;
        cycle=0;
        while((dvelocity > G.control.citcom_accu || (dpressure > G.control.citcom_accu && G.control.check_v_only == 0) ) && cycle < G.control.max_citcom_cycle)
        {
                for(i=0; i<nel; i++)
                        z1[i] = BAB[i] * r1[i];

		MPI_Allreduce(&temp, &r1dotz1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		r1dotz1 = global_p_prod(r1, z1, levmax);


                if(cycle == 0)
                {
                        for(i=0; i<nel; i++)
                                s2[i] = z1[i];
                }
                else
                {
			r0dotr0 = global_p_prod(r0, z0, levmax);

			assert(r0dotr0 > 0.0);

                        delta = r1dotz1 / r0dotr0;

                        for(i=0; i<nel; i++)
                                s2[i] = z1[i] + delta * s1[i];
                }

		for(i=0; i<neq_m; i++)
			Ah[i] = 0.0;

		G_prod_p(G.Grow_val[levmax], G.Grow_col[levmax], G.Gncol[levmax], s2, Ah, neq_m, nel, levmax);



		for(i=0; i<neq; i++)
			u1[i] = 0.0;

		if(strcmp(G.solver, "jc") == 0)
			jacobi(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy, 1000);
		else if(strcmp(G.solver, "gs") == 0)
			gauss_seidel(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy, 1000);
		else if(strcmp(G.solver, "cg") == 0)
			Conj_Grad(G.Krow_val[levmax], G.Krow_col[levmax], G.Kncol[levmax], Ah, u1, G.neq[levmax], levmax, G.control.accuracy,G.control.max_cg_cycle, 'v');
		else if(strcmp(G.solver, "mg") == 0)
			multigrid(G.Krow_val, G.Krow_col, G.Kncol, Ah, u1, G.neq, G.control.accuracy);

		for(i=0; i<neq; i++)
			Ah[i] = 0.0;

		GT_prod_d(G.Grow_val[levmax], G.Grow_col[levmax], G.Gncol[levmax], u1, Ah, neq_m, nel, levmax);

		s2dotAhat = global_p_prod(s2, Ah, levmax);

                alpha = r1dotz1 / s2dotAhat;
                for(i=0; i<nel; i++)
                {
                        r2[i] = r1[i] - alpha * Ah[i];
                        G.p[levmax][i] += alpha * s2[i];
                }

                for(i=0; i<neq; i++)
                        G.d[levmax][i] -= alpha * u1[i];

		temp1 = global_p_prod(s2, s2, levmax);

		temp2 = global_p_prod(G.p[levmax], G.p[levmax], levmax);

		temp3 = global_v_prod(u1, u1, levmax);

		temp4 = global_v_prod(G.d[levmax], G.d[levmax], levmax);

                dpressure = alpha * sqrt(temp1 / (1.0e-32 + temp2));

                dvelocity = alpha * sqrt(temp3 / (1.0e-32 + temp4));

		if(G.ip == 0 && G.control.show_convg)
			printf("%d dv = %.10e dp = %.10e\n",cycle,dvelocity,dpressure);

		cycle++;

                shuffle = s1;
                s1 = s2;
                s2 = shuffle;
                shuffle = r0;
                r0 = r1;
                r1 = r2;
                r2 = shuffle;
                shuffle = z0;
                z0 = z1;
                z1 = shuffle;

	}

	if(G.ip == 0)
		printf("%d dv = %.10e dp = %.10e\n",cycle,dvelocity,dpressure);

	/*
	for(i=0;i<neq;i++)
	{
		printf("%d %e\n",i,G.d[levmax][i]);
	}
	*/

//	for(lev=levmax; lev>G.min_level+1; lev--)
//	        project_res(G.d[lev],G.d[lev-1],G.neq[lev],G.neq[lev-1],lev,lev-1);


        free(r0);
        free(r1);
        free(r2);
        free(z0);
        free(z1);
        free(s1);
        free(s2);
        free(Ah);
        free(u1);
	free(BAB);
	return;
}

void initial_d_p()
{
	int lev,i;
	static int been_here=0;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.d[lev] = realloc(G.d[lev], G.neq[lev] * sizeof(double));

		for(i=0; i<G.neq[lev]; i++)
			G.d[lev][i] = 0.0;
	}


	if(been_here && G.control.AMR == 0 && G.control.zero_P == 0)
		return;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.p[lev] = realloc(G.p[lev], G.nel[lev] * sizeof(double));
		for(i=0; i<G.nel[lev]; i++)
			G.p[lev][i] = 0.0;
	}

	been_here = 1;

	return;
}
void product(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr,*prod;
	int nrow , ncol, lev;
	char data_type;
{
	int i, j, c;

	/*
	for(i=0; i<nrow; i++)
        {
		prod[i] = 0.0;
                for(j=0; j<matr_ncol[i]; j++)
                {
                        c = matr_row_col[i][j];
                        prod[i] += matr_row_val[i][j] * arr[c] / G.smpi[lev].veq_nsh[i];
                }
        }

	int l, k;
	double *send, *recv;
        MPI_Status stats;

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                send = malloc(G.smpi[lev].n_veq[i] * sizeof(double));

                for(j=0; j<G.smpi[lev].n_veq[i]; j++)
                {
                        l = G.smpi[lev].veq[i][j];
                        send[j] = 0.0;
                        for(k=0; k<matr_ncol[l]; k++)
                        {
                        	c = matr_row_col[l][k];
                        	send[j] += matr_row_val[l][k] * arr[c] / G.smpi[lev].veq_nsh[l];
                        }
                }
                MPI_Send(send,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
                free(send);
        }

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                recv = malloc(G.smpi[lev].n_veq[i] * sizeof(double));
                MPI_Recv(recv,G.smpi[lev].n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

                for(j=0; j<G.smpi[lev].n_veq[i]; j++)
                {
                        k = G.smpi[lev].veq[i][j];
                        arr[k] += recv[j];
                }

                free(recv);
        }
	*/





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

        other_eq_contribution_slow(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);

	other_lg_eq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);

	return;
}


void K_prod_d(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr,*prod;
	int nrow , ncol, lev;
	char data_type;
{
	int i, j, c;

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

        other_eq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);
	if(G.lg_mul)
		other_lg_eq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);

	return;
}

void other_lg_eq_contribution(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
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
	isend = malloc(G.neq_lr[lev] * 5 * sizeof(int));
	dsend = malloc(G.neq_lr[lev] * 5 * sizeof(double));

	for(i=0; i<G.neq_lr[lev] * 5; i++)
	{
		isend[i] = -1;
		dsend[i] = 0.0;
	}

	number = 0;
	for(i=G.neq_m[lev]; i<m; i++)
	{
		for(j=0; j<matr_ncol[i]; j++)
		{
			c = matr_row_col[i][j];
			isend[number * 5 + j] = G.veq_node[lev][c];
			dsend[number * 5 + j] = matr_row_val[i][j] * arr_in[c];
		}
		isend[number*5 + 3] = G.veq_node[lev][i];
		isend[number*5 + 4] = G.veq_dof[lev][i];
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

	for(i=G.neq_m[lev]; i<m; i++)
	{
		hnode = G.veq_node[lev][i];
		d = G.veq_dof[lev][i];


		if(matr_ncol[i] == 3)//this eq already has all 3 columns, not need to do anything
			continue;

		for(j=0; j<G.smpi[lev].nshareip; j++)
		{
			for(k=0; k<nrecv[j]/5; k++)
			{
				if(d == irecv[j][k*5+4] && hnode == irecv[j][k*5+3])
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

		/*
		if(G.ip == 3)
		for(j=0; j<number/4; j++)
		{
			printf("%d %d %d %d %d %d\n",G.smpi[lev].ip[i], j, irecv[4*j+0], irecv[4*j+1], irecv[4*j+2], irecv[4*j+3]);
			hnode = irecv[j*4+3];
			k = value_in_array(hnode, G.lg_hangle[lev].node, G.lg_hangle[lev].nh);
			if(k == -1)
				continue;
			inode = G.lg_hangle[lev].ni[k];
			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
			}

		}

		for(j=0; j<number/4; j++)
		{
			hnode = irecv[j*4+3];
			k = value_in_array(hnode, G.lg_hangle[lev].node, G.lg_hangle[lev].nh);
			if(k == -1)
				continue;

			inode = G.lg_hangle[lev].ni[k];
			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
				lg_eq = -1;//just to initialize
				for(int ii=G.neq_m[lev]; ii<m; ii++)
				{
					for(int jj=0; jj<matr_ncol[ii]; jj++)
					{
						if(matr_row_col[ii][jj] == eq)
						{
							lg_eq = ii;
							break;
						}
					}
				}

				for(k=0; k<3; k++)
				{
					if(irecv[j*4+k] != hnode && irecv[j*4+k] != -1)
					{
						arr_out[lg_eq] += drecv[j*4+k];
					}
				}
			}
		}
		free(irecv);
		free(drecv);
		*/


	return;
}

void K_prod_d_old(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr,*prod;
	int nrow , ncol, lev;
	char data_type;
{
	int i, j, c;


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

        other_eq_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1,data_type);

	return;
}

void G_prod_p(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr,*prod;
	int nrow, ncol, lev;
{

	int i, j, c;

	for(i=0; i<nrow; i++)
	{
		prod[i] = 0.0;
		for(j=0; j<matr_ncol[i]; j++)
		{
			c = matr_row_col[i][j];
			prod[i] += matr_row_val[i][j] * arr[c];
		}
	}

        other_vel_contribution(matr_row_val,matr_row_col,matr_ncol,prod,arr,nrow,lev,1);

	return;
}

void other_vel_contribution(matr_row_val, matr_row_col, matr_ncol, arr_out, arr_in, m, lev, log)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr_out, *arr_in;
	int m,lev,log;
{

	int i,j,k,l,iel,c;
	double *send;
	double *recv;
	int *n_veq;
	MPI_Status stats;

	if(G.lg_mul)
		n_veq = G.smpi[lev].n_veq_m;
	else
		n_veq = G.smpi[lev].n_veq;

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		send = malloc(n_veq[i] * sizeof(double));
		for(j=0; j<n_veq[i]; j++)
		{
			send[j] = 0.0;
			k = G.smpi[lev].veq[i][j];
			for(iel=0; iel<G.smpi[lev].nel[i][j]; iel++)
			{
				l = G.smpi[lev].el[i][iel][j];
				c = value_in_sorted_array(l, matr_row_col[k], matr_ncol[k]);
				send[j] += matr_row_val[k][c] * arr_in[l];
			}
		}
		MPI_Send(send,n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD);
		free(send);
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		recv = malloc(n_veq[i] * sizeof(double));
		MPI_Recv(recv,n_veq[i],MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
		for(j=0; j<n_veq[i]; j++)
		{
			k = G.smpi[lev].veq[i][j];
			arr_out[k] += recv[j];
		}

		free(recv);
	}

	return;
}

void GT_prod_d(matr_row_val, matr_row_col, matr_ncol, arr, prod, nrow, ncol, lev)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	int nrow, ncol, lev;
	double *arr,*prod;
{
	int i, j, c;

	for(i=0; i<ncol; i++)
		prod[i] = 0.0;

	for(i=0; i<nrow; i++)
	{
		for(j=0; j<matr_ncol[i]; j++)
		{
			c = matr_row_col[i][j];
			prod[c] += matr_row_val[i][j] * arr[i];
		}
	}

	return;
}

void multigrid(matr_row_val, matr_row_col, matr_ncol, arr, x, m, accu)
	double ***matr_row_val;
	int ***matr_row_col, **matr_ncol;
        double *arr, *x;
	double accu;
        int *m;
{
        double *x_tmp;
        int convg,cycle,levmax,lev,i;

	if(G.min_level >= G.max_level -1)
	{
		if(G.ip == 0)
			fprintf(stderr,"Error: for MG solver, min_level needs to be smaller than max_level-1\n");
		terminate();
	}


	levmax = G.max_level - 1;

        for(lev=G.min_level; lev<G.max_level; lev++)
                G.err[lev] = realloc(G.err[lev], m[lev] * sizeof(double));

	for(lev=G.min_level; lev<G.max_level; lev++)
		initial_error(G.err[lev], m[lev]);


	x_tmp = malloc(m[levmax] * sizeof(double));

        convg = 0;
        cycle = 0;
        while(convg==0 && cycle < G.control.mg_cycle)
        {
                for(i=0; i<m[levmax]; i++)
                        x_tmp[i] = x[i];

               	MG(matr_row_val, matr_row_col, matr_ncol, arr, x, m, levmax);

                convg = check_convergency(x, x_tmp, m[levmax], accu);

                cycle++;
        }

	free(x_tmp);
	return;
}

void MG(matr_row_val, matr_row_col, matr_ncol, arr, x, m, lev)
	double ***matr_row_val;
	int ***matr_row_col, **matr_ncol;
        double *arr, *x;
        int *m, lev;
{
	int i;
        double *res_u,*res_l;

        res_u = malloc(m[lev] * sizeof(double));
        res_l = malloc(m[lev-1] * sizeof(double));

	if(G.np == 1)
		gauss_seidel(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], arr, x, m[lev], lev, G.control.accuracy, G.control.mg_finest_iteration);
	else 
//		jacobi(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], arr, x, m[lev], lev, G.control.accuracy, G.control.mg_finest_iteration);

		Conj_Grad(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], arr, x, m[lev], lev, G.control.accuracy, G.control.mg_finest_iteration, 'v');


        get_residual(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], arr, x, res_u, m[lev], lev, 'v');


        project_res(res_u, res_l, m[lev], m[lev-1], lev, lev-1);


        lev--;

        if(lev==G.min_level)
	{
		if(G.np == 1)
			gauss_seidel(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], res_l, G.err[lev], m[lev], lev, G.control.accuracy, 1000);
		else 
		{
			jacobi(matr_row_val[lev], matr_row_col[lev], matr_ncol[lev], res_l, G.err[lev], m[lev], lev, G.control.accuracy, 1000);
		}
	}
        else
               	MG(matr_row_val, matr_row_col, matr_ncol, res_l, G.err[lev], m, lev);

        interp_error(G.err[lev],G.err[lev+1],m[lev],m[lev+1],lev,lev+1);

        for(i=0; i<m[lev+1]; i++)
               x[i]+=G.err[lev+1][i];

	if(G.np == 1)
		gauss_seidel(matr_row_val[lev+1], matr_row_col[lev+1], matr_ncol[lev+1], arr, x, m[lev+1], lev+1, G.control.accuracy, G.control.mg_coarse_iteration);
	else
		jacobi(matr_row_val[lev+1], matr_row_col[lev+1], matr_ncol[lev+1], arr, x, m[lev+1], lev+1, G.control.accuracy, G.control.mg_coarse_iteration);

	free(res_u);
	free(res_l);

	return;
}
void project_res(res_u, res_l, mu, ml, ulev, llev)
        double *res_u, *res_l;
        int ulev, llev, mu, ml;
{
	int i,l,u,j;
	int number;
	int ip;
	double *buff_d = NULL;
	MPI_Status stats;


	for(i=0; i<G.project[llev].nlocal; i++)
	{
		l = G.project[llev].local_l[i];
		u = G.project[llev].local_u[i];
		res_l[l] = res_u[u];
	}

//	if(G.project[llev].nlocal == ml)
//		return;//all nodes are local


	for(i=0; i<G.project[llev].nip_send; i++)
	{
		ip = G.project[llev].ip_send[i];
		number = G.project[llev].node_amount[i] * 2;
                buff_d = malloc(number * sizeof(double));
		for(j=0; j<number; j++)
			buff_d[j] = 0.0;


		for(j=0; j<G.project[llev].nsend[i]; j++)
		{
			l = G.project[llev].send_l[i][j];
			u = G.project[llev].send_u[i][j];
			if(u >= 0)
				buff_d[l] = res_u[u];
		}

		MPI_Send(buff_d,number,MPI_DOUBLE,ip,0,MPI_COMM_WORLD);
		free(buff_d);
	}

	for(i=0; i<G.project[llev].nip_recv; i++)
        {
                ip = G.project[llev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_DOUBLE, &number);
                buff_d = malloc(number * sizeof(double));
                MPI_Recv(buff_d,number,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&stats);

		for(j=0; j<G.project[llev].nrecv[i]; j++)
		{
			l = G.project[llev].recv_l[i][j];
			u = G.project[llev].recv_u[i][j];
			res_l[l] = buff_d[u];
		}
		free(buff_d);
        }

	return;
}

void project_res_slow(res_u, res_l, mu, ml, ulev, llev)
        double *res_u, *res_l;
        int ulev, llev, mu, ml;
{
	int i,l,u,d,k,j,node,eq;
	int number;
	int ip,loc;
	char *checked;
	double *buff_d = NULL;
	int *buff_i = NULL;
	MPI_Status stats;


	checked = malloc(ml * sizeof(char));
	for(l=0; l<ml; l++) 
		checked[l] = 'n';

	for(l=0; l<ml; l++)
	{
		if(G.project[llev].local_eq_n[l] == 1)
		{
			u = G.project[llev].local_eq[0][l];
			res_l[l] = res_u[u];
			checked[l] = 'y';
		}
	}

	for(i=0; i<G.project[llev].nip_send; i++)
	{
		ip = G.project[llev].ip_send[i];
		number = G.project[llev].node_amount[i];
                buff_d = malloc(number * 2 * sizeof(double));

		loc = 0; 
		for(j=0; j<G.project[llev].node_amount[i]; j++)
		{
			for(d=0; d<NSD; d++)
                        {
                                eq = j * NSD + d;
                                u = G.project[llev].remote_eq[i][loc++];
				if(u >= 0)
	                                buff_d[eq] = res_u[u];
				else 
					buff_d[eq] = 0.0;
                        }
		}
		MPI_Send(G.project[llev].remote_node[i],number,MPI_INT,ip,0,MPI_COMM_WORLD);
		MPI_Send(buff_d,number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD);
		free(buff_d);
	}


	for(i=0; i<G.project[llev].nip_recv; i++)
        {
                ip = G.project[llev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &number);
                buff_i = malloc(number * sizeof(int));
                MPI_Recv(buff_i,number,MPI_INT,ip,0,MPI_COMM_WORLD,&stats);

                buff_d = malloc(number * 2 * sizeof(double));
                MPI_Recv(buff_d,number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&stats);

                for(k=0; k<number; k++)
                {
                        node = buff_i[k];
                        j = node_index(node, llev);
                        if(j >= 0)
                        {
                                for(d=0; d<NSD; d++)
                                {
                                        l = G.id[llev][d][j];
					if(l >= 0 && checked[l] == 'n')
					{
						res_l[l] = buff_d[d + k * NSD];
						checked[l] = 'y';
					}
                                }
                        }
                }

		free(buff_i);
		free(buff_d);
        }


	free(checked);

	return;
}


void interp_error_fast(err_l,err_u,ml,mu,llev,ulev)
        double *err_l,*err_u;
        int ulev,llev,ml,mu;
{
	int i,j,k,d,eq;
	int u,l;
	int ip;
	int loc;
	int number;
	double *buff_d = NULL;
	char *checked;
	MPI_Status stats;

	checked = malloc(mu * sizeof(char));
        for(u=0; u<mu; u++)
                checked[u] = 'n';
	for(i=0; i<G.interp[ulev].nlocal; i++)
	{
		u = G.interp[ulev].local_u[i];
		err_u[u] = 0.0;
		for(d=0; d<G.interp[ulev].local_nlv[i]; d++)
		{
			l = G.interp[ulev].local_lv[d][i];
			err_u[u] += err_l[l];
		}
		err_u[u] /= G.interp[ulev].local_eq_n[u];
		checked[u] = 'y';
	}

	for(i=0; i<G.interp[ulev].nip_send; i++)
	{
		ip = G.interp[ulev].ip_send[i];

		number = G.interp[ulev].node_amount[i];
		buff_d = malloc(number * 2 * sizeof(double));

		loc = 0;
		for(j=0; j<G.interp[ulev].node_amount[i]; j++)
		{
			for(d=0; d<NSD; d++)
			{
				eq = j * NSD + d;
				buff_d[eq] = 0.0;
				for(k=0; k<G.interp[ulev].remote_node_n[i][j]; k++)
				{
					l = G.interp[ulev].remote_eq[i][loc++];
					if(l != -1)
					{
						buff_d[eq] += err_l[l];
					}
				}
				buff_d[eq] /= G.interp[ulev].remote_node_n[i][j];
			}
		}
		MPI_Send(buff_d,number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD);

		free(buff_d);
	}

        for(i=0; i<G.interp[ulev].nip_recv; i++)
        {
                ip = G.interp[ulev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_DOUBLE, &number);
                buff_d = malloc(number * sizeof(double));
                MPI_Recv(buff_d,number,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&stats);

                for(j=0; j<G.interp[llev].nrecv[i]; j++)
                {
                        l = G.interp[llev].recv_l[i][j];
                        u = G.interp[llev].recv_u[i][j];
			err_u[u] = buff_d[l];
			checked[u] = 'y';
                }
                free(buff_d);
        }

	free(checked);
	return;
}

void interp_error(err_l,err_u,ml,mu,llev,ulev)
        double *err_l,*err_u;
        int ulev,llev,ml,mu;
{
	int i,j,k,d,eq,node;
	int u,l;
	char *checked;
	int ip;
	int loc;
	int number;
	double **buff_d, *recv_d[100];
	int **buff_i;
	MPI_Status stats;

	MPI_Request irequest[100];
	MPI_Request drequest[100];
	int nir, ndr;

	checked = malloc(mu * sizeof(char));
	for(u=0; u<mu; u++)
		checked[u] = 'n';

	number = 0;
	for(u=0; u<mu; u++)
	{
		if(G.interp[ulev].local_eq_n[u] > 0)
		{
			err_u[u] = 0.0;
			for(i=0; i<G.interp[ulev].local_eq_n[u]; i++)
			{
				l = G.interp[ulev].local_eq[i][u];
				if(l != -1)//boundary nodes
				{
					err_u[u] += err_l[l];
				}
			}
			err_u[u] /= G.interp[ulev].local_eq_n[u];
			checked[u] = 'y';
			number++;
		}
	}


	buff_d = malloc(G.interp[ulev].nip_send * sizeof(double *));
	nir = 0;
	ndr = 0;
	for(i=0; i<G.interp[ulev].nip_send; i++)
	{
		ip = G.interp[ulev].ip_send[i];

		number = G.interp[ulev].node_amount[i];

		buff_d[i] = malloc(number * 2 * sizeof(double));


		loc = 0;
		for(j=0; j<G.interp[ulev].node_amount[i]; j++)
		{
			for(d=0; d<NSD; d++)
			{
				eq = j * NSD + d;
				buff_d[i][eq] = 0.0;
				for(k=0; k<G.interp[ulev].remote_node_n[i][j]; k++)
				{
					l = G.interp[ulev].remote_eq[i][loc++];
					if(l != -1)
					{
						buff_d[i][eq] += err_l[l];
					}
				}
				buff_d[i][eq] /= G.interp[ulev].remote_node_n[i][j];
			}
		}
		MPI_Isend(G.interp[ulev].remote_node[i],number,MPI_INT,ip,0,MPI_COMM_WORLD,&irequest[nir++]);
		MPI_Isend(buff_d[i],number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&drequest[ndr++]);
	}

	buff_i = malloc(G.interp[ulev].nip_recv * sizeof(int *));
        for(i=0; i<G.interp[ulev].nip_recv; i++)
        {
                ip = G.interp[ulev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &number);
                buff_i[i] = malloc(number * sizeof(int));
                MPI_Irecv(buff_i[i],number,MPI_INT,ip,0,MPI_COMM_WORLD,&irequest[nir++]);

                recv_d[i] = malloc(number * 2 * sizeof(double));
                MPI_Irecv(recv_d[i],number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&drequest[ndr++]);
	}

	MPI_Waitall(nir, irequest, &stats);
	MPI_Waitall(ndr, drequest, &stats);


        for(i=0; i<G.interp[ulev].nip_recv; i++)
	{
		for(k=0; k<number; k++)
		{
			node = buff_i[i][k];	
			j = node_index(node,ulev);
			if(j >= 0)
			{
				for(d=0; d<NSD; d++)
				{
					u = G.id[ulev][d][j];
					if(u >= 0 && checked[u] == 'n')
					{
						err_u[u] = recv_d[i][d + k * NSD];
						checked[u] = 'y';
					}
				}
			}
		}
        }

        for(i=0; i<G.interp[ulev].nip_recv; i++)
		free(recv_d[i]);

        for(i=0; i<G.interp[ulev].nip_recv; i++)
		free(buff_i[i]);
	free(buff_i);


	for(i=0; i<G.interp[ulev].nip_send; i++)
		free(buff_d[i]);
	free(buff_d);




	free(checked);

	return;
}

void interp_error_blocked_comm(err_l,err_u,ml,mu,llev,ulev)
        double *err_l,*err_u;
        int ulev,llev,ml,mu;
{
	int i,j,k,d,eq,node;
	int u,l;
	char *checked;
	int ip;
	int loc;
	int number;
	double *buff_d = NULL;
	int *buff_i = NULL;
	MPI_Status stats;

	checked = malloc(mu * sizeof(char));
	for(u=0; u<mu; u++)
		checked[u] = 'n';

	number = 0;
	for(u=0; u<mu; u++)
	{
		if(G.interp[ulev].local_eq_n[u] > 0)
		{
			err_u[u] = 0.0;
			for(i=0; i<G.interp[ulev].local_eq_n[u]; i++)
			{
				l = G.interp[ulev].local_eq[i][u];
				if(l != -1)//boundary nodes
				{
					err_u[u] += err_l[l];
				}
			}
			err_u[u] /= G.interp[ulev].local_eq_n[u];
			checked[u] = 'y';
			number++;
		}
	}


	for(i=0; i<G.interp[ulev].nip_send; i++)
	{
		ip = G.interp[ulev].ip_send[i];

		number = G.interp[ulev].node_amount[i];

		buff_d = malloc(number * 2 * sizeof(double));


		loc = 0;
		for(j=0; j<G.interp[ulev].node_amount[i]; j++)
		{
			for(d=0; d<NSD; d++)
			{
				eq = j * NSD + d;
				buff_d[eq] = 0.0;
				for(k=0; k<G.interp[ulev].remote_node_n[i][j]; k++)
				{
					l = G.interp[ulev].remote_eq[i][loc++];
					if(l != -1)
					{
						buff_d[eq] += err_l[l];
					}
				}
				buff_d[eq] /= G.interp[ulev].remote_node_n[i][j];
			}
		}
		MPI_Send(G.interp[ulev].remote_node[i],number,MPI_INT,ip,0,MPI_COMM_WORLD);
		MPI_Send(buff_d,number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD);

		free(buff_d);
	}

        for(i=0; i<G.interp[ulev].nip_recv; i++)
        {
                ip = G.interp[ulev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &number);
                buff_i = malloc(number * sizeof(int));
                MPI_Recv(buff_i,number,MPI_INT,ip,0,MPI_COMM_WORLD,&stats);

                buff_d = malloc(number * 2 * sizeof(double));
                MPI_Recv(buff_d,number*2,MPI_DOUBLE,ip,0,MPI_COMM_WORLD,&stats);

		for(k=0; k<number; k++)
		{
			node = buff_i[k];	
			j = node_index(node,ulev);
			if(j >= 0)
			{
				for(d=0; d<NSD; d++)
				{
					u = G.id[ulev][d][j];
					if(u >= 0 && checked[u] == 'n')
					{
						err_u[u] = buff_d[d + k * NSD];
						checked[u] = 'y';
					}
				}
			}
		}

		free(buff_d);
		free(buff_i);
        }

	free(checked);

	return;
}

void initial_error(arr,m)
        double *arr;
        int m;
{
        int i;
        for(i=0; i<m; i++)
                arr[i]=0.0;
        return;
}


int check_convergency(new,old,m,error)
        double *new,*old,error;
        int m;
{
        int i;
        double diff,temp,diff_global;

        diff=0.0;
        for(i=0; i<m; i++)
        {
                temp = fabs(new[i] - old[i]);
                if(temp > diff)
                        diff = temp;
        }

	MPI_Allreduce(&diff,&diff_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

        if(diff_global < error)
                return 1;
        else
                return 0;

}

void get_residual(matr_row_val, matr_row_col, matr_ncol, arr, x, res, m, lev, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
        double *arr,*x,*res;
        int m,lev;
	char data_type;
{
	int i, j, c;

	for(i=0; i<m; i++)
		res[i] = arr[i];//do not forget to initialize, otherwiseprod[c] below is not be defined yet


	for(i=0; i<m; i++)
	{
		for(j=0; j<matr_ncol[i]; j++)
		{
			c = matr_row_col[i][j];
			res[i] -= matr_row_val[i][j] * x[c];
		}
	}

        other_eq_contribution(matr_row_val, matr_row_col, matr_ncol, res, x, m, lev, 2, data_type);

        return;
}


void v_to_nodes()
{
	int i,d,j;
	int lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(d=0; d<NSD; d++)
			G.V[lev][d] = realloc(G.V[lev][d], G.nno[lev] * sizeof(double));

		for(i=0; i<G.nno[lev]; i++)
		{
			for(d=0; d<NSD; d++)
			{
				j = G.id[lev][d][i];
				if(j >= 0)
					G.V[lev][d][i] = G.d[lev][j];
				else
					G.V[lev][d][i] = G.bc.v_bc_val[lev][d][i];
			}
		}
	}

	return;
}

void p_to_nodes(P, PN, lev)
	double *P, *PN;
	int lev;
{
	int i, a, e;
	double total_area, area;

	for(i=0; i<G.nno[lev]; i++)
	{
		PN[i] = 0.0;
		total_area = 0.0;
		for(a=0; a<G.p2n[lev].nel[i]; a++)
		{
			e = G.p2n[lev].el[a][i];
			area = G.eco[lev].area[e];
			PN[i] += area * P[e];
			total_area += area;
		}
		PN[i] /= total_area;
	}
	return;
}

void enforce_constrain_veq(x, lev)
	int lev;
	double *x;
{
	int i,eq,j,c,node,k,d,eqc,inode;
	double temp[NSD];
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
	eqval = malloc(size * NSD * sizeof(double));


	for(i=0; i<G.hangle[lev].nh; i++)
	{
		node = G.hangle[lev].node[i];

		inode = G.hangle[lev].ni[i];

		hnode[number] = node;

		for(d=0; d<NSD; d++)
		{
			temp[d] = 0.0;
			for(c=0; c<NSD; c++)
			{
				k = G.hangle[lev].ci[c][i];
				eqc = G.id[lev][d][k];

				if(eqc != -1)
					temp[d] += x[eqc];
				else
					temp[d] += G.bc.v_bc_val[lev][d][k];
			}
			temp[d] *= 0.5;

			eqval[number * NSD + d] = temp[d];
		}

		number++;

		if(number >= size)
		{
			size *= 2;
			hnode = realloc(hnode, size * sizeof(int));
			eqval = realloc(eqval, size * NSD * sizeof(double));
		}

		if(inode >= 0)
		{
			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
				x[eq] = temp[d];
			}
		}
	}


	hnode = realloc(hnode, number * sizeof(int));
	eqval = realloc(eqval, number * NSD * sizeof(double));

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(hnode, number, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
		MPI_Send(eqval, number * NSD, MPI_DOUBLE, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}

	free(hnode);
	free(eqval);

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &iamount);
                irecv = malloc(iamount * sizeof(int));
                MPI_Recv(irecv,iamount,MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);



		damount = iamount * NSD;

                drecv = malloc(damount * sizeof(double));
                MPI_Recv(drecv,damount,MPI_DOUBLE,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		for(j=0; j<iamount; j++)
		{
			node = irecv[j];
			k = node_index(node,lev);
			if(k >= 0)
			{
				for(d=0; d<NSD; d++)
				{
					eq = G.id[lev][d][k];
					if(eq >= 0)
					{
						x[eq] = drecv[j * NSD + d];
					}
				}
			}
		}
		free(irecv);
		free(drecv);
	}
	return;
}

void Conj_Grad_new(matr_row_val, matr_row_col, matr_ncol, b, x, n, lev, accu, max_cycle, data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *b, *x, accu;
	int n, lev, max_cycle;
	char data_type;
{
	int count, i;
	double *r0, *r1, *r2, *z0, *z1, *p1, *p2, *Ap, *shuffle;
	double alpha, beta, dotprod, dotr1z1, dotr0z0;
	double residual;
	double gprod();
	void product();


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

	residual = global_v_prod(r1, r1, lev);

	assert(residual != 0.0  /* initial residual for CG = 0.0 */);

	count = 0;
	while(((residual > accu) && (count < max_cycle)) || count == 0)
	{
		if(G.lg_mul)
		{
			for(i=0; i<n; i++)
				z1[i] = r1[i];
		}
		else
		{
			for(i=0; i<n; i++)
				z1[i] = r1[i]/G.Kdiag[lev][i];
		}

		dotr1z1 = global_v_prod(r1, z1, lev);


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

		K_prod_d(matr_row_val, matr_row_col, matr_ncol, p2, Ap, n, n, lev, data_type);

		dotprod = global_v_prod(p2, Ap, lev);

		if(dotprod == 0.0)
			alpha = 1e-3;
		else
			alpha = dotr1z1 / dotprod;

		for(i=0; i<n; i++)
		{
			x[i] += alpha * p2[i];
		        r2[i] = r1[i] - alpha * Ap[i];
		}
		//if(G.lg_mul == 0)
			enforce_constrain_veq(x, lev);

		residual = global_v_prod(r2, r2, lev);

		shuffle = r0; r0 = r1; r1 = r2; r2 = shuffle;
		shuffle = z0; z0 = z1; z1 = shuffle;
		shuffle = p1; p1 = p2; p2 = shuffle;
		count++;
	}
	if(G.ip == 0)
		printf("CG cycle = %d\n",count);

	free(r0);
	free(r1);
	free(r2);
	free(z0);
	free(z1);
	free(p1);
	free(p2);
	free(Ap);

	return;
}


void Conj_Grad(matr_row_val, matr_row_col, matr_ncol, b, x, n, lev, accu, max_cycle, data_type)
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
		r=malloc(n*sizeof(double));
		s=malloc(n*sizeof(double));
		d=malloc(n*sizeof(double));
		q=malloc(n*sizeof(double));


		K_prod_d(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
		for(i=0; i<n; i++)
		{
			r[i] = b[i] - r[i];
			d[i] = r[i] / G.Kdiag[lev][i];
		}

		/*
		temp=0.0;
		for(i=0;i<n;i++)
			temp += r[i] * d[i];
		MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		*/

		delta_new = global_v_prod(r, d, lev);

		delta0 = delta_new;

		cycle = 0;
		while(cycle < max_cycle && delta_new > accu * delta0)
		{
			K_prod_d(matr_row_val, matr_row_col, matr_ncol, d, q, n, n, lev, data_type);

			/*
			temp=0.0;
			for(i=0; i<n; i++)
				temp += d[i] * q[i];
			MPI_Allreduce(&temp, &temp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			temp1 = global_v_prod(d, q, lev);

			alpha = delta_new / temp1;

			for(i=0; i<n; i++)
				x[i] += alpha * d[i];

			enforce_constrain_veq(x, lev);

			if(cycle % 50 == 0)
			{
				K_prod_d(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
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
				s[i] = r[i] / G.Kdiag[lev][i];
			}

			delta_old = delta_new;

			/*
			temp = 0.0;
			for(i=0; i<n; i++)
				temp += r[i] * s[i];
			MPI_Allreduce(&temp, &delta_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			delta_new = global_v_prod(r, s, lev);

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


		K_prod_d(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
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

		delta_new = global_v_prod(r, d, lev);

		delta0 = delta_new;

		cycle = 0;
		while(cycle < max_cycle && delta_new > accu * delta0)
		{
			K_prod_d(matr_row_val, matr_row_col, matr_ncol, d, q, n, n, lev, data_type);

			/*
			temp=0.0;
			for(i=0; i<n; i++)
				temp += d[i] * q[i];
			MPI_Allreduce(&temp, &temp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			*/

			temp1 = global_v_prod(d, q, lev);

			alpha = delta_new / temp1;

			for(i=0; i<n; i++)
				x[i] += alpha * d[i];

			enforce_constrain_veq(x, lev);

			if(cycle % 50 == 0)
			{
				K_prod_d(matr_row_val, matr_row_col, matr_ncol, x, r, n, n, lev, data_type);
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
			delta_new = global_v_prod(r, r, lev);


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
			if(global_residual < accu)
				break;

		}

		free(r);
		free(d);
		free(q);
	}

	return;
}


void jacobi(matr_row_val, matr_row_col, matr_ncol, arr, x, m, lev, accu, max_cycle)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
        double *arr, *x;
	double accu;
        int m,lev,max_cycle;
{
	int i, j, c, cycle, idiag;
	double *x_tmp, *prod;
	double diff, temp, diff_global,diag;

	cycle = 0;
	x_tmp = malloc(m * sizeof(double));
	prod = malloc(m * sizeof(double));
	while(cycle < max_cycle)
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
			assert(idiag != 0);//for security reason only, very likely not needed
			x[i] = (arr[i] - prod[i]) / diag;
		}

		other_eq_contribution(matr_row_val,matr_row_col,matr_ncol,x,x_tmp,m,lev,0,'v');
//		enforce_constrain_veq(x, lev);

		cycle++;

		diff = 0.0;
		for(i=0; i<m; i++)
		{
			temp = fabs(x[i] - x_tmp[i]);
			if(temp > diff)
				diff = temp;
		}

		MPI_Allreduce(&diff,&diff_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

		if(diff_global < accu)
			break;
	}

	free(prod);
	free(x_tmp);
	return;
}

void gauss_seidel(matr_row_val, matr_row_col, matr_ncol, arr, x, m, lev, accu, max_cycle)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
        double *arr, *x;
	double accu;
        int m,lev,max_cycle;
{
	int i, j, c, cycle, idiag;
	double *x_tmp, *prod;
	double diff, temp, diff_global,diag;


	cycle = 0;
	x_tmp = malloc(m * sizeof(double));
	prod = malloc(m * sizeof(double));
	while(cycle < max_cycle)
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
					prod[i] += matr_row_val[i][j] * x[c];
				}
				else
				{
					idiag = 1;
					diag = matr_row_val[i][j];
				}
			}
			assert(idiag != 0);//for security reason only, very likely not needed
			x[i] = (arr[i] - prod[i]) / diag;
		}

		cycle++;

		diff = 0.0;
		for(i=0; i<m; i++)
		{
			temp = fabs(x[i] - x_tmp[i]);
			if(temp > diff)
				diff = temp;
		}

		MPI_Allreduce(&diff,&diff_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

		if(diff_global < accu)
			break;
	}

	free(prod);
	free(x_tmp);
	return;
}


void other_eq_contribution(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr_out, *arr_in;
	int m,lev,method;
	char data_type;
{
	int i, j, r, c, l;
	double *send, **recv;
	int *nrecv;
	MPI_Status stats;

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		send = malloc(G.smpi[lev].send_number[i] * sizeof(double));
		for(j=0; j<G.smpi[lev].send_number[i]; j++)
		{
			l = G.smpi[lev].send_row[i][j];
			r = G.smpi[lev].send_matr_row[i][j];
			c = G.smpi[lev].send_matr_col[i][j];
			send[j] = matr_row_val[r][c] * arr_in[l];
		}
                MPI_Send(send, G.smpi[lev].send_number[i], MPI_DOUBLE, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
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
	for(j=0; j<G.smpi[lev].recv_number[i]; j++)
	{
		l = G.smpi[lev].recv_row[i][j];
	//	ip = G.smpi[lev].recv_matr_ip[i][j];
		c = G.smpi[lev].recv_matr_col[i][j];

		if(method == 0)
		{
	//		a = value_in_sorted_array(l, matr_row_col[l], matr_ncol[l]);
	//		arr_out[l] -= recv[i][c]/matr_row_val[l][a];
			arr_out[l] -= recv[i][c]/G.Kdiag[lev][l];//G.Kdiag is used here, meaning this function is only for velocity
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


void other_eq_contribution_slow(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
	double **matr_row_val;
	int **matr_row_col, *matr_ncol;
	double *arr_out, *arr_in;
	int m,lev,method;
	char data_type;
{
        int i, j, k, c, l, a, b, drow, dcol, node, inode, d;
        int *nseq, nsh, *ip, **seq, *eq_node, *eq_dof, *index;
        int *send_row_node, *send_col_node, *send_row_dof, *send_col_dof, nsend, size, *send_lg_info;
        double *send_prod;
        int **recv_row_node, **recv_col_node, *nrecv, **recv_lg_info;
        int **recv_row_dof, **recv_col_dof;
	int *share_node, nshare_node, *checked, ncheck;
        double **recv_prod;

	static double time_used = 0.0;
	double start_time, end_time;
        MPI_Status stats;

	start_time = MPI_Wtime();;

        nseq = G.smpi[lev].n_veq;
        seq = G.smpi[lev].veq;

        nsh = G.smpi[lev].nshareip;
        ip = G.smpi[lev].ip;

        eq_node = G.veq_node[lev];
        eq_dof = G.veq_dof[lev];

        send_row_node = NULL; //node for each eq on the row index
        send_col_node = NULL; //node for each eq on the col index
        send_row_dof = NULL;  //dof for each eq on the row index
        send_col_dof = NULL;  //dof for each eq on the col index
        send_lg_info = NULL;  
        send_prod = NULL;    //matrix production for each eq

        index = malloc(m * sizeof(int)); // will be used to check if an eq is shared. Only contribution of not-shared eq need to added because the shared components have already been considered in K_prod_d and other function that call other_eq_contribution.


	//get data for each unshared eqs and send to sharing ips for later usage
        for(i=0; i<nsh; i++)
        {
                nsend = 0;
                size = 1;
                send_row_node = malloc(size * sizeof(int));
                send_col_node = malloc(size * sizeof(int));
                send_row_dof = malloc(size * sizeof(int));
                send_col_dof = malloc(size * sizeof(int));
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

//				if(G.ip == 2 && k == 3 && G.smpi[lev].ip[i] == 1)
//					printf("%d %d %d %d\n",k, c, G.veq_node[lev][c],index[c]);

				if(index[c] == -1)
				{
					send_row_node[nsend] = eq_node[k];
					send_col_node[nsend] = eq_node[c];
					/*
					if(l >= G.Kncol_m[lev][k])
						send_col_node[nsend] = 0 - eq_node[c];
					else
						send_col_node[nsend] = eq_node[c];
						*/

					send_row_dof[nsend] = eq_dof[k];
					send_col_dof[nsend] = eq_dof[c];

					send_prod[nsend] = matr_row_val[k][l] * arr_in[c];

					//if(l >= G.Kncol_m[lev][k])
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
						send_row_dof = realloc(send_row_dof, size * sizeof(int));
						send_col_dof = realloc(send_col_dof, size * sizeof(int));
						send_lg_info = realloc(send_lg_info, size * sizeof(int));
						send_prod = realloc(send_prod, size * sizeof(double));
					}
				}
			}
		}

                send_row_node = realloc(send_row_node, nsend * sizeof(int));
                send_col_node = realloc(send_col_node, nsend * sizeof(int));
                send_row_dof = realloc(send_row_dof, nsend * sizeof(int));
                send_col_dof = realloc(send_col_dof, nsend * sizeof(int));
		send_lg_info = realloc(send_lg_info, nsend * sizeof(int));
                send_prod = realloc(send_prod, nsend * sizeof(double));

                MPI_Send(send_row_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_col_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_row_dof, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_col_dof, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_lg_info, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_prod, nsend, MPI_DOUBLE, ip[i], 0, MPI_COMM_WORLD);

                free(send_row_node);
                free(send_col_node);
                free(send_row_dof);
                free(send_col_dof);
                free(send_lg_info);
                free(send_prod);
        }

	//recv data
        recv_row_node = malloc(nsh * sizeof(int *));
        recv_col_node = malloc(nsh * sizeof(int *));
        recv_row_dof = malloc(nsh * sizeof(int *));
        recv_col_dof = malloc(nsh * sizeof(int *));
        recv_lg_info = malloc(nsh * sizeof(int *));
        recv_prod = malloc(nsh * sizeof(double *));
        nrecv = malloc(nsh * sizeof(int));
        for(i=0; i<nsh; i++)
        {
                MPI_Probe(ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                recv_row_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_row_dof[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_dof[i] = malloc(nrecv[i] * sizeof(int));
                recv_lg_info[i] = malloc(nrecv[i] * sizeof(int));
                recv_prod[i] = malloc(nrecv[i] * sizeof(double));

                MPI_Recv(recv_row_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_row_dof[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_dof[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_lg_info[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_prod[i],nrecv[i],MPI_DOUBLE,ip[i],0,MPI_COMM_WORLD,&stats);
        }

	nshare_node = G.smpi[lev].nshare_node;
	share_node = G.smpi[lev].share_node;

	for(dcol=0; dcol<NSD; dcol++)//dof for col of eq, need to use all NSD because both direction of V contribute to the matrix
	{
		for(i=0; i<nshare_node; i++)
		{
			node = share_node[i];
			inode = node_index(node, lev);
			for(drow=0; drow<NSD; drow++)//dof for row of eq, this is used to get eqs related to nshare_node
			{
				l = G.id[lev][drow][inode];

				if(l >= 0)
				{
					ncheck = 0;
					size = 1;
					checked = malloc(size * sizeof(int));

					for(j=0; j<nsh; j++)
					{
						for(k=0; k<nrecv[j]; k++)
						{
							if(eq_node[l] == recv_row_node[j][k] && recv_col_dof[j][k] == dcol && recv_row_dof[j][k] == drow)
							{
								//must void consider the contribution of one node more than once for each dof, remove duplication

								if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1
										|| recv_lg_info[j][k] == 1)

//									       	|| (value_in_array(recv_col_node[j][k], G.lg_hangle[lev].node, G.lg_hangle[lev].nh) == -1 && node_index(recv_col_node[j][k], lev) >= 0))
								{


//				if(G.ip == 1 && l == 8 && G.smpi[lev].ip[j] == 0)
//					printf("%e %d %d %d %d %e\n",arr_out[l],recv_row_node[j][k], recv_row_dof[j][k],recv_col_node[j][k], recv_col_dof[j][k],recv_prod[j][k]);

									if(method == 0)
									{
										c = value_in_sorted_array(l, matr_row_col[l], matr_ncol[l]);
										arr_out[l] -= recv_prod[j][k]/matr_row_val[l][c];
									}
									else if(method == 1)
									{
										arr_out[l] += recv_prod[j][k];
									}
									else if(method == 2)
									{
										arr_out[l] -= recv_prod[j][k];
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
		}
	}
        for(i=0; i<nsh; i++)
        {       
                free(recv_row_node[i]);
                free(recv_col_node[i]);
                free(recv_row_dof[i]);
                free(recv_col_dof[i]);
                free(recv_lg_info[i]);
                free(recv_prod[i]);
        }
        free(recv_row_node);
        free(recv_col_node);
        free(recv_row_dof);
        free(recv_col_dof);
        free(recv_lg_info);
        free(recv_prod);
        free(index);

	end_time = MPI_Wtime();

	time_used += end_time - start_time;

	G.CPU_time = time_used;

//	if(G.ip == 0)
//		fprintf(stderr, "time used in other_eq_contribution: %e\n", time_used);


	return;
}


void other_eq_contribution_old(matr_row_val,matr_row_col,matr_ncol,arr_out,arr_in,m,lev,method,data_type)
        double **matr_row_val;
        int **matr_row_col, *matr_ncol;
        double *arr_out, *arr_in;
        int m,lev,method;
        char data_type;
{
        int i, j, k, c, l;
        int *nseq, nsh, *ip, **seq;
        double *send, *recv;
        int *index;
        MPI_Status stats;


        nseq = NULL;
        seq = NULL;

        nseq = G.smpi[lev].n_veq;
        seq = G.smpi[lev].veq;

        nsh = G.smpi[lev].nshareip;
        ip = G.smpi[lev].ip;

        index = malloc(m * sizeof(int));

        for(i=0; i<nsh; i++)
        {
                send = malloc(nseq[i] * sizeof(double));

                for(j=0; j<m; j++)
                        index[j] = value_in_sorted_array(j, seq[i], nseq[i]);

                for(j=0; j<nseq[i]; j++)
                {
                        send[j] = 0.0;
                        k = seq[i][j];
                        for(l=0; l<matr_ncol[k]; l++)
                        {
                                c = matr_row_col[k][l];
                                if(index[c] == -1)
                                {
                                        send[j] += matr_row_val[k][l] * arr_in[c];
                                }
                        }
                }

                MPI_Send(send, nseq[i], MPI_DOUBLE, ip[i], 0, MPI_COMM_WORLD);
                free(send);
        }

        free(index);

        for(i=0; i<nsh; i++)
        {
                recv = malloc(nseq[i] * sizeof(double));
                MPI_Recv(recv, nseq[i], MPI_DOUBLE, ip[i], 0, MPI_COMM_WORLD, &stats);
                for(j=0; j<nseq[i]; j++)
                {
                        k = seq[i][j];
                        if(method == 0)
                        {
                                c = value_in_sorted_array(k, matr_row_col[k], matr_ncol[k]);
                                arr_out[k] -= recv[j]/matr_row_val[k][c];
                        }
                        else if(method == 1)
                                arr_out[k] += recv[j];
                        else if(method == 2)
                                arr_out[k] -= recv[j];
                }

                free(recv);
        }
        return;
}


void setup_MPI_comm_vframe(lev)
	int lev;
{
        int i, j, k, c, l, drow, dcol, node, inode, m;
        int *nseq, nsh, *ip, **seq, *eq_node, *eq_dof, *index;
        int *send_row_node, *send_col_node, *send_row_dof, *send_col_dof, nsend, size, *send_lg_info;
        int **recv_row_node, **recv_col_node, *nrecv, **recv_lg_info;
        int **recv_row_dof, **recv_col_dof;
	int *share_node, nshare_node, *checked, ncheck;
	int *matr_ncol, **matr_row_col;
	int *number, *tmp_size;

        MPI_Status stats;

        nseq = G.smpi[lev].n_veq;
        seq = G.smpi[lev].veq;

        nsh = G.smpi[lev].nshareip;
        ip = G.smpi[lev].ip;

        eq_node = G.veq_node[lev];
        eq_dof = G.veq_dof[lev];

        send_row_node = NULL; //node for each eq on the row index
        send_col_node = NULL; //node for each eq on the col index
        send_row_dof = NULL;  //dof for each eq on the row index
        send_col_dof = NULL;  //dof for each eq on the col index
        send_lg_info = NULL;  //dof for each eq on the col index
	m = G.neq[lev];//number of eqs


        index = malloc(m * sizeof(int)); // will be used to check if an eq is shared. Only contribution of not-shared eq need to added because the shared components have already been considered in K_prod_d and other function that call other_eq_contribution.

	matr_row_col = G.Krow_col[lev];
	matr_ncol = G.Kncol[lev];


	nshare_node = G.smpi[lev].nshare_node; //all shareing node without duplication
	share_node = G.smpi[lev].share_node;


	tmp_size = malloc(nsh * sizeof(int));
	number = malloc(nsh * sizeof(int));

	for(i=0; i<nsh; i++)
	{
		number[i] = 0;
		tmp_size[i] = 1;
		G.smpi[lev].send_row[i] = realloc(G.smpi[lev].send_row[i], tmp_size[i] * sizeof(int));
		G.smpi[lev].send_matr_row[i] = realloc(G.smpi[lev].send_matr_row[i], tmp_size[i] * sizeof(int));
		G.smpi[lev].send_matr_col[i] = realloc(G.smpi[lev].send_matr_col[i], tmp_size[i] * sizeof(int));
	}

	//get data for each unshared eqs and send to sharing ips for later usage
        for(i=0; i<nsh; i++)
        {
                nsend = 0;
                size = 1;
                send_row_node = malloc(size * sizeof(int));
                send_col_node = malloc(size * sizeof(int));
                send_row_dof = malloc(size * sizeof(int));
                send_col_dof = malloc(size * sizeof(int));
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

					send_row_dof[nsend] = eq_dof[k];
					send_col_dof[nsend] = eq_dof[c];
					if(value_in_array(eq_node[c], G.lg_hangle[lev].node, G.lg_hangle[lev].nh) >= 0)
                                                send_lg_info[nsend] = 1;
                                        else
                                                send_lg_info[nsend] = 0;

					G.smpi[lev].send_row[i][ number[i] ] = c;
					G.smpi[lev].send_matr_row[i][ number[i] ] = k;
					G.smpi[lev].send_matr_col[i][ number[i] ] = l;
					number[i]++;
					if(number[i] >= tmp_size[i])
					{
						tmp_size[i] *= 2;
						G.smpi[lev].send_row[i] = realloc(G.smpi[lev].send_row[i], tmp_size[i] * sizeof(int));
						G.smpi[lev].send_matr_row[i] = realloc(G.smpi[lev].send_matr_row[i], tmp_size[i] * sizeof(int));
						G.smpi[lev].send_matr_col[i] = realloc(G.smpi[lev].send_matr_col[i], tmp_size[i] * sizeof(int));
					}
					

					nsend++;
					if(nsend >= size)
					{
						size *= 2;
						send_row_node = realloc(send_row_node, size * sizeof(int));
						send_col_node = realloc(send_col_node, size * sizeof(int));
						send_row_dof = realloc(send_row_dof, size * sizeof(int));
						send_col_dof = realloc(send_col_dof, size * sizeof(int));
						send_lg_info = realloc(send_lg_info, size * sizeof(int));
					}
				}
			}
		}

                send_row_node = realloc(send_row_node, nsend * sizeof(int));
                send_col_node = realloc(send_col_node, nsend * sizeof(int));
                send_row_dof = realloc(send_row_dof, nsend * sizeof(int));
                send_col_dof = realloc(send_col_dof, nsend * sizeof(int));
                send_lg_info = realloc(send_lg_info, nsend * sizeof(int));

                MPI_Send(send_row_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_col_node, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_row_dof, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_col_dof, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);
                MPI_Send(send_lg_info, nsend, MPI_INT, ip[i], 0, MPI_COMM_WORLD);

                free(send_row_node);
                free(send_col_node);
                free(send_row_dof);
                free(send_col_dof);
                free(send_lg_info);
        }

	for(i=0; i<nsh; i++)
	{
		G.smpi[lev].send_row[i] = realloc(G.smpi[lev].send_row[i], number[i] * sizeof(int));
		G.smpi[lev].send_matr_row[i] = realloc(G.smpi[lev].send_matr_row[i], number[i] * sizeof(int));
		G.smpi[lev].send_matr_col[i] = realloc(G.smpi[lev].send_matr_col[i], number[i] * sizeof(int));
		G.smpi[lev].send_number[i] = number[i];
	}


	//recv data
        recv_row_node = malloc(nsh * sizeof(int *));
        recv_col_node = malloc(nsh * sizeof(int *));
        recv_row_dof = malloc(nsh * sizeof(int *));
        recv_col_dof = malloc(nsh * sizeof(int *));
        recv_lg_info = malloc(nsh * sizeof(int *));
        nrecv = malloc(nsh * sizeof(int));
        for(i=0; i<nsh; i++)
        {
                MPI_Probe(ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                recv_row_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_node[i] = malloc(nrecv[i] * sizeof(int));
                recv_row_dof[i] = malloc(nrecv[i] * sizeof(int));
                recv_col_dof[i] = malloc(nrecv[i] * sizeof(int));
                recv_lg_info[i] = malloc(nrecv[i] * sizeof(int));

                MPI_Recv(recv_row_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_node[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_row_dof[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_col_dof[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
                MPI_Recv(recv_lg_info[i],nrecv[i],MPI_INT,ip[i],0,MPI_COMM_WORLD,&stats);
        }



	for(i=0; i<nsh; i++)
	{
		number[i] = 0;
		tmp_size[i] = 1;
		G.smpi[lev].recv_row[i] = realloc(G.smpi[lev].recv_row[i], tmp_size[i] * sizeof(int));
		G.smpi[lev].recv_matr_ip[i] = realloc(G.smpi[lev].recv_matr_ip[i], tmp_size[i] * sizeof(int));
		G.smpi[lev].recv_matr_col[i] = realloc(G.smpi[lev].recv_matr_col[i], tmp_size[i] * sizeof(int));
	}


	for(dcol=0; dcol<NSD; dcol++)//dof for col of eq, need to use all NSD because both direction of V contribute to the matrix
	{
		for(i=0; i<nshare_node; i++)
		{
			node = share_node[i];
			inode = node_index(node, lev);
			for(drow=0; drow<NSD; drow++)//dof for row of eq, this is used to get eqs related to nshare_node
			{
				l = G.id[lev][drow][inode];

				if(l >= 0)
				{
					ncheck = 0;
					size = 1;
					checked = malloc(size * sizeof(int));

					for(j=0; j<nsh; j++)
					{
						for(k=0; k<nrecv[j]; k++)
						{
							if(eq_node[l] == recv_row_node[j][k] && recv_col_dof[j][k] == dcol && recv_row_dof[j][k] == drow)
							{
								//must void consider the contribution of one node more than once for each dof, remove duplication
								if(value_in_array(recv_col_node[j][k], checked, ncheck) == -1 || recv_lg_info[j][k] == 1)
								{

									G.smpi[lev].recv_row[j][ number[j] ] = l;
									G.smpi[lev].recv_matr_ip[j][ number[j] ] = j;
									G.smpi[lev].recv_matr_col[j][ number[j] ] = k;
									number[j]++;
									if(number[j] >= tmp_size[j])
									{
										tmp_size[j] *= 2;
										G.smpi[lev].recv_row[j] = realloc(G.smpi[lev].recv_row[j], tmp_size[j] * sizeof(int));
										G.smpi[lev].recv_matr_ip[j] = realloc(G.smpi[lev].recv_matr_ip[j], tmp_size[j] * sizeof(int));
										G.smpi[lev].recv_matr_col[j] = realloc(G.smpi[lev].recv_matr_col[j], tmp_size[j] * sizeof(int));
									}

									//arr_out[l] += recv_prod[j][k];
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
		}
	}


	for(i=0; i<nsh; i++)
	{
		G.smpi[lev].recv_row[i] = realloc(G.smpi[lev].recv_row[i], number[i] * sizeof(int));
		G.smpi[lev].recv_matr_ip[i] = realloc(G.smpi[lev].recv_matr_ip[i], number[i] * sizeof(int));
		G.smpi[lev].recv_matr_col[i] = realloc(G.smpi[lev].recv_matr_col[i], number[i] * sizeof(int));
		G.smpi[lev].recv_number[i] = number[i];
	}
	free(tmp_size);
	free(number);

        for(i=0; i<nsh; i++)
        {       
                free(recv_row_node[i]);
                free(recv_col_node[i]);
                free(recv_row_dof[i]);
                free(recv_col_dof[i]);
                free(recv_lg_info[i]);
        }
        free(recv_row_node);
        free(recv_col_node);
        free(recv_row_dof);
        free(recv_col_dof);
        free(recv_lg_info);
        free(index);
	/*
	if(G.ip == 0)
		fprintf(stderr, "time used in other_eq_contribution: lev = %d %e\n", lev, time_used);
		*/

	return;
}
