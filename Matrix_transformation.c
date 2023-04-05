#include "global_variables.h"

void add_lagrange_neq(lev)
	int lev;
{
	int i, j, d, eq, inode, nh, size, node, cnode[VPTS], c, ci[VPTS], lg_row, lg_row_T;
	int *send, number, *recv;
	MPI_Status stats;
	void reassign_lg_hangle_memory();

	//the hangle data structure include nodes that are not in the currect ip but are in other ips, it also does not consider sharing hangling nodes
	//This is becuase hangle data is used in enforce_constrain function, which does not need sharing hangling nodes but need to find avereage of solutions on a node to be shared, even the node does not exist in the current ip.
	//the lg_hangle data contain sharing nodes and remove ghost nodes. This data is required for lagrange multiplier method
		nh = 0;
		size = 1;
		reassign_lg_hangle_memory(size, lev);

		send = malloc(G.hangle[lev].nh * (1 + NSD) * sizeof(int));
		number = 0;
		//first, for hangling node existing in the current ip
		for(i=0; i<G.hangle[lev].nh; i++)
		{
			inode = G.hangle[lev].ni[i];
			node = G.hangle[lev].node[i];
                        for(c=0; c<2; c++)
			{
				cnode[c] = G.hangle[lev].cnode[c][i];
				ci[c] = G.hangle[lev].ci[c][i];
			}

			send[number++] = node;
			send[number++] = cnode[0];
			send[number++] = cnode[1];

			if(inode >= 0)
			{
				G.lg_hangle[lev].node[nh] = node;
				G.lg_hangle[lev].ni[nh] = inode;

                                for(c=0; c<2; c++)
				{
                                        G.lg_hangle[lev].cnode[c][nh] = cnode[c];
                                        G.lg_hangle[lev].ci[c][nh] = ci[c];
				}
				nh++;
				G.lg_hangle[lev].nh = nh;
				if(nh >= size)
				{
					size *= 2;
					reassign_lg_hangle_memory(size, lev);
				}
			}
		}

		//sharing hangle nodes to other ips
		send = realloc(send, number * sizeof(int));
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			MPI_Send(send, number, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
		}
		free(send);

		//append remote hangle nodes to this current ip
		for(i=0; i<G.smpi[lev].nshareip; i++)
		{
			MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
			MPI_Get_count(&stats, MPI_INT, &number);
			recv = malloc(number * sizeof(int));
			MPI_Recv(recv,number,MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

			for(j=0; j<number; j += 3)
			{
				node = recv[j];
				cnode[0] = recv[j + 1];
				cnode[1] = recv[j + 2];

				inode = node_index(node, lev);
				if(inode >= 0)
				{
					G.lg_hangle[lev].node[nh] = node;
					G.lg_hangle[lev].ni[nh] = inode;
					for(c=0; c<2; c++)
					{
						G.lg_hangle[lev].cnode[c][nh] = cnode[c];
						ci[c] = node_index(cnode[c], lev);
						G.lg_hangle[lev].ci[c][nh] = ci[c];
					}
					nh++;
					G.lg_hangle[lev].nh = nh;
					if(nh >= size)
					{
						size *= 2;
						reassign_lg_hangle_memory(size, lev);
					}
				}
			}
			free(recv);
		}
		reassign_lg_hangle_memory(nh, lev);


		// calculate neq_c, neq_lr, and update neq
		G.neq_c[lev] = 0;
		G.heat[lev].neq_c = 0;
		for(i=0; i<G.lg_hangle[lev].nh; i++)
		{
			inode = G.lg_hangle[lev].ni[i];
			if(G.heat[lev].id[inode] >= 0)
				G.heat[lev].neq_c++;
			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
				if(eq >= 0)
				{
					G.neq_c[lev]++;
				}
			}
		}

		G.neq_lr[lev] = G.neq_c[lev];
		G.neq[lev] = G.neq_m[lev] + G.neq_lr[lev];

		G.veq_node[lev] = realloc(G.veq_node[lev], G.neq[lev] * sizeof(int));
		G.veq_dof[lev] = realloc(G.veq_dof[lev], G.neq[lev] * sizeof(int));

		G.heat[lev].neq_lr = G.heat[lev].neq_c;
		G.heat[lev].neq = G.heat[lev].neq_m + G.heat[lev].neq_lr;

		G.heat[lev].eq_node = realloc(G.heat[lev].eq_node, G.heat[lev].neq * sizeof(int));

		lg_row = G.neq_m[lev];
		lg_row_T = G.heat[lev].neq_m;
		for(i=0; i<G.lg_hangle[lev].nh; i++)
		{
			node = G.lg_hangle[lev].node[i];
			inode = G.lg_hangle[lev].ni[i];

			if(G.heat[lev].id[inode] >= 0)
			{
				G.heat[lev].eq_node[lg_row_T] = node;
				lg_row_T++;
			}

			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
				if(eq >= 0)
				{
					G.veq_node[lev][lg_row] = node;
					G.veq_dof[lev][lg_row] = d;
					lg_row++;
				}
			}
		}
	return;
}

void reassign_lg_hangle_memory(size, lev)
        int size, lev;
{
        int c;
        G.lg_hangle[lev].node = realloc(G.lg_hangle[lev].node, size * sizeof(int));
        G.lg_hangle[lev].ni = realloc(G.lg_hangle[lev].ni, size * sizeof(int));
        for(c=0; c<NSD; c++)
        {
                G.lg_hangle[lev].cnode[c] = realloc(G.lg_hangle[lev].cnode[c], size * sizeof(int));
                G.lg_hangle[lev].ci[c] = realloc(G.lg_hangle[lev].ci[c], size * sizeof(int));
        }
}


void build_lagrange_K_matrix(lev)
	int lev;
{
	int i, j, inode, d, eq, eq_c, ic, dc;
	int lg_row;

	for(i=G.neq_m[lev]; i<G.neq[lev]; i++)
		G.Kncol[lev][i] = 0;

	for(i=G.neq_m[lev]; i<G.neq[lev]; i++)
		G.F[lev][i] = 0.0;


	//fill in lagrange values for K and F matrix

	lg_row = G.neq_m[lev];
	for(i=0; i<G.lg_hangle[lev].nh; i++)
	{
		inode = G.lg_hangle[lev].ni[i];
		for(d=0; d<NSD; d++)
		{
			eq = G.id[lev][d][inode];
			if(eq >= 0)
			{
				//append col
				G.Krow_val[lev][eq][ G.Kncol[lev][eq] ] = 2.0;
				G.Krow_col[lev][eq][ G.Kncol[lev][eq] ] = lg_row ;
				G.Kncol[lev][eq]++;

				//append row
				G.Krow_val[lev][lg_row][ G.Kncol[lev][lg_row] ] = 2.0;
				G.Krow_col[lev][lg_row][ G.Kncol[lev][lg_row] ] = eq;
				G.Kncol[lev][lg_row]++;


				for(dc=0; dc<2; dc++)
				{
					ic = G.lg_hangle[lev].ci[dc][i];
					if(ic >= 0)
						eq_c = G.id[lev][d][ic];//this should be d, not dc
					else
						eq_c = -1;

					if(eq_c >= 0)
					{
						//append col
						G.Krow_val[lev][eq_c][ G.Kncol[lev][eq_c] ] = -1.0;
						G.Krow_col[lev][eq_c][ G.Kncol[lev][eq_c] ] = lg_row;
						G.Kncol[lev][eq_c]++;

						//append row
						G.Krow_val[lev][lg_row][ G.Kncol[lev][lg_row] ] = -1.0;
						G.Krow_col[lev][lg_row][ G.Kncol[lev][lg_row] ] = eq_c;
						G.Kncol[lev][lg_row]++;
					}
					else
					{
						G.F[lev][lg_row] += 0.5 * G.bc.v_bc_val[lev][dc][ic];
					}
				}

				lg_row++;
			}
		}
	}

	//sort cols of lagrange rows
        int log;
        for(i=G.neq_m[lev]; i<G.neq[lev]; i++)
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

	find_nshare_for_lg_eqs(lev);
	/*
	for(i=0; i<G.neq[lev]; i++)
	{
		for(j=0; j<G.Kncol[lev][i]; j++)
		{
			printf("%d %d %e\n",i, G.Krow_col[lev][i][j],G.Krow_val[lev][i][j]);
		}
	}

	terminate();
	*/


	return;
}

void find_nshare_for_lg_eqs(lev)
	int lev;
{
	int i, j, eq, neq, neq_m, node, inode, d;
	int *send_hnode, number, lg_eq;
	int *recv_hnode;
	static int been_here = 0;
	MPI_Status stats;

	if(been_here)
		return;
	been_here++;

	neq = G.neq[lev];
	neq_m = G.neq_m[lev];

	//increase array size
	G.smpi[lev].veq_nsh = realloc(G.smpi[lev].veq_nsh, neq * sizeof(int));

	//for convinience, also append sharing lg eqs 
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].n_veq_m[i] = G.smpi[lev].n_veq[i];
		G.smpi[lev].veq[i] = realloc(G.smpi[lev].veq[i], (G.smpi[lev].n_veq[i]+2*G.lg_hangle[lev].nh) * sizeof(int));
	}

	for(i=neq_m; i<neq; i++)
	{
		G.smpi[lev].veq_nsh[i] = 1;
	}

	send_hnode = malloc(G.neq_lr[lev] * sizeof(int));

	number = 0;
	for(i=0; i<G.lg_hangle[lev].nh; i++)
	{
		node = G.lg_hangle[lev].node[i];
		send_hnode[number++] = node;
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(send_hnode, number, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
		MPI_Get_count(&stats, MPI_INT, &number);
		recv_hnode = malloc(number * sizeof(int));
		MPI_Recv(recv_hnode,number,MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		for(j=0; j<number; j++)
		{
			node = recv_hnode[j];	
			inode = node_index(node, lev);
			if(inode >= 0)
			{
				for(d=0; d<NSD; d++)
				{
					eq = G.id[lev][d][inode];
					lg_eq = -1;
					for(int ii=neq_m; ii<neq; ii++)
					{
						for(int jj=0; jj<G.Kncol[lev][ii]; jj++)
						{
							if(G.Krow_col[lev][ii][jj] == eq)
							{
								lg_eq = ii;
								break;
							}
						}
					}
					G.smpi[lev].veq_nsh[lg_eq]++;
					G.smpi[lev].veq[i][ G.smpi[lev].n_veq[i]++ ] = lg_eq;
				}
			}
		}
		free(recv_hnode);
	}

	/*
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<G.smpi[lev].n_veq[i]; j++)
		{
			printf("%d %d %d %d %d\n",G.ip, G.smpi[lev].ip[i],j,G.smpi[lev].veq[i][j],G.veq_node[lev][G.smpi[lev].veq[i][j]]);
		}
	}
	*/

	free(send_hnode);
	return;
}

void build_lagrange_heat_KC_matrix(lev)
	int lev;
{
	int i, j, inode, eq, eq_c, ic, dc;
	int lg_row;


	for(i=G.heat[lev].neq_m; i<G.heat[lev].neq; i++)
	{
		G.heat[lev].Kncol[i] = 0;
		G.heat[lev].Cncol[i] = 0;
	}

	for(i=G.heat[lev].neq_m; i<G.heat[lev].neq; i++)
		G.heat[lev].F[i] = 0.0;


	//fill in lagrange values for K and F matrix

	lg_row = G.heat[lev].neq_m;
	for(i=0; i<G.lg_hangle[lev].nh; i++)
	{
		inode = G.lg_hangle[lev].ni[i];
		eq = G.heat[lev].id[inode];

		if(eq == -1)
			continue;

		//append col
		G.heat[lev].Krow_val[eq][ G.heat[lev].Kncol[eq] ] = 2.0;
		G.heat[lev].Krow_col[eq][ G.heat[lev].Kncol[eq] ] = lg_row;
		G.heat[lev].Kncol[eq]++;

		G.heat[lev].Crow_val[eq][ G.heat[lev].Cncol[eq] ] = 2.0;
		G.heat[lev].Crow_col[eq][ G.heat[lev].Cncol[eq] ] = lg_row;
		G.heat[lev].Cncol[eq]++;

		//append row
		G.heat[lev].Krow_val[lg_row][ G.heat[lev].Kncol[lg_row] ] = 2.0;
		G.heat[lev].Krow_col[lg_row][ G.heat[lev].Kncol[lg_row] ] = eq;
		G.heat[lev].Kncol[lg_row]++;

		G.heat[lev].Crow_val[lg_row][ G.heat[lev].Cncol[lg_row] ] = 2.0;
		G.heat[lev].Crow_col[lg_row][ G.heat[lev].Cncol[lg_row] ] = eq;
		G.heat[lev].Cncol[lg_row]++;


		for(dc=0; dc<2; dc++)
		{
			ic = G.lg_hangle[lev].ci[dc][i];
			if(ic >= 0)
				eq_c = G.heat[lev].id[ic];
			else
				eq_c = -1;

			if(eq_c >= 0)
			{
				//append col
				G.heat[lev].Krow_val[eq_c][ G.heat[lev].Kncol[eq_c] ] = -1.0;
				G.heat[lev].Krow_col[eq_c][ G.heat[lev].Kncol[eq_c] ] = lg_row;
				G.heat[lev].Kncol[eq_c]++;

				G.heat[lev].Crow_val[eq_c][ G.heat[lev].Cncol[eq_c] ] = -1.0;
				G.heat[lev].Crow_col[eq_c][ G.heat[lev].Cncol[eq_c] ] = lg_row;
				G.heat[lev].Cncol[eq_c]++;

				//append row
				G.heat[lev].Krow_val[lg_row][ G.heat[lev].Kncol[lg_row] ] = -1.0;
				G.heat[lev].Krow_col[lg_row][ G.heat[lev].Kncol[lg_row] ] = eq_c;
				G.heat[lev].Kncol[lg_row]++;

				G.heat[lev].Crow_val[lg_row][ G.heat[lev].Cncol[lg_row] ] = -1.0;
				G.heat[lev].Crow_col[lg_row][ G.heat[lev].Cncol[lg_row] ] = eq_c;
				G.heat[lev].Cncol[lg_row]++;
			}
			else
			{
				G.heat[lev].F[lg_row] += 0.5 * G.bc.v_bc_val[lev][dc][ic];
			}
		}
		lg_row++;
	}

	//sort cols of lagrange rows
        int log;
        for(i=G.heat[lev].neq_m; i<G.heat[lev].neq; i++)
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

        for(i=G.heat[lev].neq_m; i<G.heat[lev].neq; i++)
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

	find_nshare_for_lg_T_eqs(lev);
	return;
}

void find_nshare_for_lg_T_eqs(lev)
	int lev;
{
	int i, j, eq, neq, neq_m, node, inode;
	int *send_hnode, number, lg_eq;
	int *recv_hnode;
	static int been_here = 0;
	MPI_Status stats;

	if(been_here)
		return;
	been_here++;

	neq = G.heat[lev].neq;
	neq_m = G.heat[lev].neq_m;

	//increase array size
	G.smpi[lev].Teq_nsh = realloc(G.smpi[lev].Teq_nsh, neq * sizeof(int));

	//for convinience, also append sharing lg eqs 
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].n_Teq_m[i] = G.smpi[lev].n_Teq[i];
		G.smpi[lev].Teq[i] = realloc(G.smpi[lev].Teq[i], (G.smpi[lev].n_Teq[i]+G.lg_hangle[lev].nh) * sizeof(int));
	}

	for(i=neq_m; i<neq; i++)
	{
		G.smpi[lev].Teq_nsh[i] = 1;
	}

	send_hnode = malloc(G.heat[lev].neq_lr * sizeof(int));

	number = 0;
	for(i=0; i<G.lg_hangle[lev].nh; i++)
	{
		node = G.lg_hangle[lev].node[i];
		send_hnode[number++] = node;
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(send_hnode, number, MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
		MPI_Get_count(&stats, MPI_INT, &number);
		recv_hnode = malloc(number * sizeof(int));
		MPI_Recv(recv_hnode,number,MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);

		for(j=0; j<number; j++)
		{
			node = recv_hnode[j];	
			inode = node_index(node, lev);
			if(inode >= 0)
			{
				eq = G.heat[lev].id[inode];
				lg_eq = -1;
				for(int ii=neq_m; ii<neq; ii++)
				{
					for(int jj=0; jj<G.heat[lev].Kncol[ii]; jj++)
					{
						if(G.heat[lev].Krow_col[ii][jj] == eq)
						{
							lg_eq = ii;
							break;
						}
					}
				}
				G.smpi[lev].Teq_nsh[lg_eq]++;
				G.smpi[lev].Teq[i][ G.smpi[lev].n_Teq[i]++ ] = lg_eq;
			}
		}
		free(recv_hnode);
	}


	/*
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		printf("%d %d %d\n",G.ip, G.smpi[lev].ip[i],G.smpi[lev].n_Teq[i]);
	}
	terminate();
	*/

	free(send_hnode);
	return;
}
