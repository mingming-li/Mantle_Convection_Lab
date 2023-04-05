#include "global_variables.h"

void construct_connection()
{
        int lev;

        for(lev=G.min_level; lev<G.max_level; lev++)
	{
                proc_conn(lev);
		find_hangling_nodes(lev);
		if(G.lg_mul)
			add_lagrange_neq(lev);
		find_sharing_nodes(lev);
		find_sharing_eqs(lev);
		elements_for_sharing_eqs(lev);//used in G_prod_d function
		find_nshare_for_eqs(lev);
		find_nshare_for_Teqs(lev);
	}

	if(strcmp(G.solver, "mg") == 0)
	{
		for(lev=G.min_level; lev<G.max_level-1; lev++)
		{
			build_interp_conn(lev);//lev is lower level, should be <G.max_leve-1
			setup_interp_comm_frame(lev);
			build_project_conn(lev);//lev is lower level, should be <G.max_leve-1
			setup_project_comm_frame(lev);
		}
	}

        return;
}

void proc_conn(lev)
	int lev;
{
        int i,node,ip[VPTS],n,j,p,nsh;

	nsh = 0;
        for(i=0; i<G.nno[lev]; i++)
        {
                node = G.node[lev][i];
                n = what_ip_node_in(node, lev, ip);
                for(j=0; j<n; j++)
                {
			p = value_in_array(ip[j], G.smpi[lev].ip, nsh);
			if(p == -1)
                        {
                                G.smpi[lev].ip[nsh] = ip[j];
                                nsh++;
                                G.smpi[lev].nshareip = nsh;
                        }
                }
	}

	return;
}

void find_hangling_nodes(lev)
	int lev;
{
	int e, level, i, j, n, nh, c, p, k;
	int node[VPTS], ip[VPTS];
	int *nsend, **send, *size;
	int *nrecv, **recv;
	MPI_Status stats;
	void reassign_hangle_memory();

	//assign memory for sending
	//nsend: number of data send
	//size: for dynamic memory control
	//send: data
	//this is sending to multiple ips, so nsend is an array
	nsend = malloc(G.smpi[lev].nshareip * sizeof(int));
	size = malloc(G.smpi[lev].nshareip * sizeof(int));
	send = malloc(G.smpi[lev].nshareip * sizeof(int *));
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		nsend[i] = 0;
		size[i] = 1 + NSD;
		send[i] = malloc(size[i] * sizeof(int));
	}


	//hangle is an struct data type
	//nh: the number of hangling nodes
	//node: global index for hangling nodes
	//ni: local index for hangling nodes
	//cnode: global index for connecting nodes, with a size of NSD
	//ci: local index for cnode

	G.hangle[lev].nh = 0;
	nh = 1;
	reassign_hangle_memory(nh, lev);
        for(e=0; e<G.nel[lev]; e++)
        {
                level = Mcode_lev(G.Mcode[lev][e]);
                if(level == lev)
			continue;//the element is already in its finest mesh at this level, no hangling nodes

		//the node[i] are for global index for potential hangling nodes, 1 element has 4 pot. hangling nodes
                node[0] = (G.ien[lev][0][e]+G.ien[lev][1][e])/2;
                node[1] = (G.ien[lev][1][e]+G.ien[lev][2][e])/2;
                node[2] = (G.ien[lev][2][e]+G.ien[lev][3][e])/2;
                node[3] = (G.ien[lev][3][e]+G.ien[lev][0][e])/2;

                for(i=0; i<VPTS; i++)
                {
                        j = node_index(node[i], lev);
			if(j >=0 )//local hangle node
			{
				//add hangling nodes
				G.hangle[lev].node[ G.hangle[lev].nh ] = node[i];
				for(c=0; c<2; c++)
					G.hangle[lev].cnode[c][ G.hangle[lev].nh ] = G.ien[lev][(c + i) % VPTS][e];
				G.hangle[lev].nh++;

				//dynamic memory control
				if(G.hangle[lev].nh >= nh)
				{
					nh *= 2;
					reassign_hangle_memory(nh, lev);
				}
			}
			else
			{
				//when the hangling node does not exist locally, they may exist in other processors
				//therefore, this hangling node is also important, whose solutions can be shared with other processors
                		n = what_ip_node_in(node[i], lev, ip);
				for(k=0; k<n; k++)
				{
					p = value_in_array(ip[k], G.smpi[lev].ip, G.smpi[lev].nshareip);
					//preparing data for sending
					//This data including the global index of the hangling nodes not existing locally,
					//and their connecting nodes
					//therefore, each pair of data contain 1+NSD integers
					send[p][ nsend[p]++ ] = node[i];
                        		for(c=0; c<NSD; c++)
					{
						send[p][ nsend[p]++ ] = G.ien[lev][(c + i) % VPTS][e];
					}

					if(nsend[p] >= size[p])
					{
						size[p] *= 2;
						send[p] = realloc(send[p], size[p] * sizeof(int));
					}
				}
			}
		}
	}

	//send potentially hangle nodes not existing locally to other processor for checking
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(send[i], nsend[i], MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
		free(send[i]);
	}
	free(send);
	free(nsend);

	//receiving
	nrecv = malloc(G.smpi[lev].nshareip * sizeof(int));
	recv = malloc(G.smpi[lev].nshareip * sizeof(int *));
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                recv[i] = malloc(nrecv[i] * sizeof(int));
                MPI_Recv(recv[i],nrecv[i],MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
	}

	//after receiving potential hangling nodes, we need to check whether the nodes are really hangling nodes, 
	//meaning that if they exist in the processor
	
	//prepare real remote hangle nodes for sending

	nsend = malloc(G.smpi[lev].nshareip * sizeof(int));
	size = malloc(G.smpi[lev].nshareip * sizeof(int));
	send = malloc(G.smpi[lev].nshareip * sizeof(int *));
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		nsend[i] = 0;
		size[i] = 1 + NSD;
		send[i] = malloc(size[i] * sizeof(int));
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		n = nrecv[i]/(1+NSD);//number of hangling nodes, because each pair contains 1+NSD numbers
		for(j=0; j<n; j++)
		{
			k = j * (1 + NSD);//index in recv for center node

			// if this node exists in the processor, this is real hangling nodes, 
			// otherwise, it is not and should not be considered
			if(node_index(recv[i][k], lev) >= 0) 
			{
				//data for real hangling remote nodes
				send[i][ nsend[i]++ ] = recv[i][k];
				for(c=0; c<NSD; c++)
					send[i][ nsend[i]++ ] = recv[i][k + c + 1];
				if(nsend[i] >= size[i])
				{
					size[i] *= 2;
					send[i] = realloc(send[i], size[i] * sizeof(int));
				}
			}
		}
	}

	//send real remote hangle nodes back 
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(send[i], nsend[i], MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}

	//receiving
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nrecv[i]);
                recv[i] = malloc(nrecv[i] * sizeof(int));
                MPI_Recv(recv[i],nrecv[i],MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
	}


	//append remote hangle nodes, and remove duplication
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		nh = G.hangle[lev].nh;//previous number local hangling nodes
		n = nrecv[i]/(1+NSD); // number of remote hangling nodes

		reassign_hangle_memory(n + nh, lev);//using it maximum memory, and will be reallocate later after removing duplications

		for(j=0; j<n; j++)
		{
			k = j * (1 + NSD);//index in recv for center node
			p = value_in_array(recv[i][k], G.hangle[lev].node, G.hangle[lev].nh);
			if(p == -1)//if it does not exist in previous array
			{
				G.hangle[lev].node[j + nh] = recv[i][k];//index for 
				for(c=0; c<NSD; c++)
				{
					// c+k is for indexes of connection nodes
					G.hangle[lev].cnode[c][j + nh] = recv[i][c + k + 1];
				}
				G.hangle[lev].nh++;
			}
		}
	}

	reassign_hangle_memory(G.hangle[lev].nh, lev);

	for(i=0; i<G.hangle[lev].nh; i++)
	{
		G.hangle[lev].ni[i] = node_index(G.hangle[lev].node[i], lev);
		for(c=0; c<NSD; c++)
			G.hangle[lev].ci[c][i] = node_index(G.hangle[lev].cnode[c][i], lev);

//		printf("%d | %d %d %d %d\n", G.ip, i, G.hangle[lev].node[i],G.hangle[lev].cnode[0][i],G.hangle[lev].cnode[1][i]);
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		free(send[i]);
		free(recv[i]);
	}
	free(send);
	free(recv);
	free(nsend);
	free(nrecv);
	free(size);

	return;
}

void reassign_hangle_memory(size, lev)
	int size, lev;
{
	int c;
	G.hangle[lev].node = realloc(G.hangle[lev].node, size * sizeof(int));
	G.hangle[lev].ni = realloc(G.hangle[lev].ni, size * sizeof(int));
	for(c=0; c<NSD; c++)
	{
		G.hangle[lev].cnode[c] = realloc(G.hangle[lev].cnode[c], size * sizeof(int));
		G.hangle[lev].ci[c] = realloc(G.hangle[lev].ci[c], size * sizeof(int));
	}
}

void find_sharing_nodes(lev)
	int lev;
{
	int i, ip[NEN], n, j, p, k, node;
	int *size;

	size = malloc(G.smpi[lev].nshareip * sizeof(int));
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		size[i] = 1;
		G.smpi[lev].nno[i] = 0;
		G.smpi[lev].node[i] = realloc(G.smpi[lev].node[i], size[i] * sizeof(int));
	}

	for(i=0; i<G.nno[lev]; i++)
	{
		node = G.node[lev][i];
		n = what_ip_node_in(node, lev, ip);
		for(j=0; j<n; j++)
		{
			 p = value_in_array(ip[j], G.smpi[lev].ip, G.smpi[lev].nshareip);
			 k = G.smpi[lev].nno[p];
			 G.smpi[lev].node[p][k] = node;
			 G.smpi[lev].nno[p]++;
			 if(G.smpi[lev].nno[p] >= size[p])
			 {
				size[p] *= 2;

				G.smpi[lev].node[p] = realloc(G.smpi[lev].node[p], size[p] * sizeof(int));
			 }
		}
	}
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].node[i] = realloc(G.smpi[lev].node[i], G.smpi[lev].nno[i] * sizeof(int));
	}

	//remove invalid share
	//this is because the function 'what_ip_node_in' sometimes return fake ips that having the node on edges
	//buff: global sharing nodes

	int *nbuff, **buff;
	MPI_Status stats;

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(G.smpi[lev].node[i], G.smpi[lev].nno[i], MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}

	nbuff = malloc(G.smpi[lev].nshareip * sizeof(int));
	buff = malloc(G.smpi[lev].nshareip * sizeof(int *));

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nbuff[i]);
                buff[i] = malloc(nbuff[i] * sizeof(int));
                MPI_Recv(buff[i],nbuff[i],MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
	}

	//if buff[i][j] does not exist, mark it to -1
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		for(j=0; j<nbuff[i]; j++)
		{
			k = node_index(buff[i][j], lev);
			if(k == -1)
			{
				buff[i][j] = -1;
			}
		}
	}


	//send modified buff back to original ip
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		MPI_Send(buff[i], nbuff[i], MPI_INT, G.smpi[lev].ip[i], 0, MPI_COMM_WORLD);
	}
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
                MPI_Probe(G.smpi[lev].ip[i], 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &nbuff[i]);
                buff[i] = malloc(nbuff[i] * sizeof(int));
                MPI_Recv(buff[i],nbuff[i],MPI_INT,G.smpi[lev].ip[i],0,MPI_COMM_WORLD,&stats);
	}

	//reassign sharing nodes from buff by removing buff == -1
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		k = 0;
		for(j=0; j<nbuff[i]; j++)
		{
			if(buff[i][j] != -1)
			{
				G.smpi[lev].node[i][k] = buff[i][j];
				k++;
			}
		}
		G.smpi[lev].nno[i] = k;
	}

	/*
	for(i=0; i<G.smpi[lev].nshareip; i++)
		for(j=0; j<G.smpi[lev].nno[i]; j++)
			printf("%d | %d %d %d\n",G.ip, G.smpi[lev].ip[i], j, G.smpi[lev].node[i][j]);
			*/
	
	int isize;

        G.smpi[lev].nshare_node = 0;
        isize = 1;
        G.smpi[lev].share_node = malloc(isize * sizeof(int));

        for(i=0; i<G.smpi[lev].nshareip; i++)
        {
                for(j=0; j<G.smpi[lev].nno[i]; j++)
                {
                        G.smpi[lev].share_node[ G.smpi[lev].nshare_node++ ] = G.smpi[lev].node[i][j];
                        if(G.smpi[lev].nshare_node >= isize)
                        {
                                isize *= 2;
                                G.smpi[lev].share_node = realloc(G.smpi[lev].share_node, isize * sizeof(int));
                        }
                }
        }
        G.smpi[lev].share_node = realloc(G.smpi[lev].share_node, G.smpi[lev].nshare_node * sizeof(int));
        sort(G.smpi[lev].share_node, G.smpi[lev].nshare_node);
        G.smpi[lev].nshare_node = remove_dup(G.smpi[lev].share_node, G.smpi[lev].nshare_node);


	for(i=0; i<G.smpi[lev].nshareip; i++)
		free(buff[i]);
	free(buff);
	free(nbuff);
	free(size);
}

void elements_for_sharing_eqs(lev)
	int lev;
{
	int e, i, inode, eq, j, a, d, k;

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].nel[i] = realloc(G.smpi[lev].nel[i], G.smpi[lev].n_veq[i] * sizeof(int));
		for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			G.smpi[lev].nel[i][j] = 0;

		for(a=0; a<NEN; a++)
		{
			G.smpi[lev].el[i][a] = realloc(G.smpi[lev].el[i][a], G.smpi[lev].n_veq[i] * sizeof(int));
		}
	}


	for(e=0; e<G.nel[lev]; e++)
	{
		for(a=0; a<VPTS; a++)
		{
			inode = node_index(G.ien[lev][a][e], lev);
			for(d=0; d<NSD; d++)
			{
				eq = G.id[lev][d][inode];
				for(i=0; i<G.smpi[lev].nshareip; i++)
				{
					j = value_in_sorted_array(eq, G.smpi[lev].veq[i], G.smpi[lev].n_veq[i]);
					if(j >= 0)
					{
						k = G.smpi[lev].nel[i][j];
						G.smpi[lev].el[i][k][j] = e;
						G.smpi[lev].nel[i][j]++;
					}
				}
			}
		}
	}

	/*
	for(i=0; i<G.smpi[lev].nshareip; i++)
	for(j=0; j<G.smpi[lev].n_veq[i]; j++)
	for(d=0; d<G.smpi[lev].nel[i][j]; d++)
		printf("%d | %d %d %d %d %d\n",G.ip,i,j,d,G.smpi[lev].el[i][d][j],G.smpi[lev].veq[i][j]);

	terminate();
	*/

	return;
}


void find_sharing_eqs(lev)
	int lev;
{
	int i, j, node, inode, d, Teq, veq;


	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].n_Teq[i] = 0;
		G.smpi[lev].Teq[i] = realloc(G.smpi[lev].Teq[i], G.smpi[lev].nno[i] * sizeof(int));

		G.smpi[lev].n_veq[i] = 0;
		G.smpi[lev].veq[i] = realloc(G.smpi[lev].veq[i], G.smpi[lev].nno[i] * 2 * sizeof(int));

		for(j=0; j<G.smpi[lev].nno[i]; j++)
		{
			node = G.smpi[lev].node[i][j];
			inode = node_index(node, lev);
			Teq = G.heat[lev].id[inode];
			if(Teq != -1)
			{
				G.smpi[lev].Teq[i][ G.smpi[lev].n_Teq[i] ] = Teq;
				G.smpi[lev].n_Teq[i]++;
			}
			for(d=0; d<NSD; d++)
			{
				veq = G.id[lev][d][inode];
				if(veq != -1)
				{
					G.smpi[lev].veq[i][ G.smpi[lev].n_veq[i] ] = veq;
					G.smpi[lev].n_veq[i]++;
				}
			}
		}
	}

	/*
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		sort(G.smpi[lev].veq[i], G.smpi[lev].n_veq[i]);
		sort(G.smpi[lev].Teq[i], G.smpi[lev].n_Teq[i]);
		for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			printf("%d %d %d\n",G.ip, G.smpi[lev].ip[i],G.smpi[lev].veq[i][j]);
	}

	terminate();
	*/

	return;
}


/*
void find_sharing_nodes_eqs_old(lev)
        int lev;
{
	int e, a, inode, node, i, j, eq, p, l, d;
	int ip[NEN], n;
	int *Tsize, *vsize;
	char *checked;

	checked = malloc(G.nno[lev] * sizeof(char));
	for(i=0; i<G.nno[lev]; i++)
		checked[i] = 'n';

	Tsize = malloc(G.smpi[lev].nshareip * sizeof(int));
	vsize = malloc(G.smpi[lev].nshareip * sizeof(int));

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		Tsize[i] = 1;
		vsize[i] = 1;

		G.smpi[lev].n_Teq[i] = 0;
		G.smpi[lev].Teq[i] = realloc(G.smpi[lev].Teq[i], Tsize[i] * sizeof(int));

		G.smpi[lev].n_veq[i] = 0;
		G.smpi[lev].veq[i] = realloc(G.smpi[lev].veq[i], vsize[i] * sizeof(int));
	}

	for(e=0; e<G.nel[lev]; e++)
	{
		for(a=0; a<VPTS; a++)
		{
			node = G.ien[lev][a][e];
			inode = node_index(node, lev);
			if(checked[inode] == 'n')
			{
				checked[inode] = 'y';
				n = what_ip_node_in(node, lev, ip);
				if(n == 0)
					continue;

				for(j=0; j<n; j++)
				{
					//for heat equations
					p = value_in_array(ip[j], G.smpi[lev].ip, G.smpi[lev].nshareip);
					eq = G.heat[lev].id[inode];
					if(eq != -1)
					{
						l = G.smpi[lev].n_Teq[p];
						G.smpi[lev].Teq[p][l] = eq;
						G.smpi[lev].n_Teq[p]++;

						if(G.smpi[lev].n_Teq[p] >= Tsize[p])
						{
							Tsize[p] *= 2;
							G.smpi[lev].Teq[p] = realloc(G.smpi[lev].Teq[p], Tsize[p] * sizeof(int));
						}
					}

					//for velocity equations
			                for(d=0; d<NSD; d++)
					{
						eq = G.id[lev][d][inode];

						if(eq != -1)
						{
							l = G.smpi[lev].n_veq[p];
							G.smpi[lev].veq[p][l] = eq;
							G.smpi[lev].n_veq[p]++;

							if(G.smpi[lev].n_veq[p] >= vsize[p])
							{
								vsize[p] *= 2;
								G.smpi[lev].veq[p] = realloc(G.smpi[lev].veq[p], vsize[p] * sizeof(int));
							}
						}
					}
				}
			}
		}
	}

	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		G.smpi[lev].Teq[i] = realloc(G.smpi[lev].Teq[i], G.smpi[lev].n_Teq[i] * sizeof(int));
		G.smpi[lev].veq[i] = realloc(G.smpi[lev].veq[i], G.smpi[lev].n_veq[i] * sizeof(int));
	}


	if(G.ip == 2)
	for(i=0; i<G.smpi[lev].nshareip; i++)
	{
		sort(G.smpi[lev].veq[i], G.smpi[lev].n_veq[i]);
		for(j=0; j<G.smpi[lev].n_veq[i]; j++)
			printf("%d %d %d %d\n",G.smpi[lev].ip[i],i,j,G.smpi[lev].veq[i][j]);
	}

	free(Tsize);
	free(vsize);
	free(checked);
	terminate();
        return;
};
*/

int what_ip_node_in(n,lev, ipval)
        int n,lev,*ipval;
{
        int i,nip;
        int ix[VPTS],iz[VPTS];
        int a,b;
        int e,ip;
        int ipx,ipz;
        int max_elx,max_elz,max_nel;


        max_elx = G.max_elx[G.max_level - 1];
        max_elz = G.max_elz[G.max_level - 1];
        max_nel = G.max_nel[G.max_level - 1];

        a = n / G.gnoz;
        b = n - a * G.gnoz;

        ix[0] = max(a - 1, 0);
        ix[1] = max(a - 1, 0);
        ix[2] = min(a, G.gnox - 2);
        ix[3] = min(a, G.gnox - 2);

        iz[0] = max(b - 1,0);
        iz[1] = min(b, G.gnoz - 2);
        iz[2] = min(b, G.gnoz - 2);
        iz[3] = max(b - 1,0);

        nip = 0;
        for(i=0; i<VPTS; i++)
        {
                ipx = ix[i] / max_elx;
                ipz = iz[i] / max_elz;

                if(iz[i] == G.gnoz - 1)
                        ipz--;
                if(ix[i] == G.gnox - 1)
                        ipx--;

                ip = ipz + ipx * G.npz;

                ix[i] %= max_elx;
                iz[i] %= max_elz;

                e = morton(ix[i],iz[i]) + max_nel * ip;

                for(ip=0; ip<G.np; ip++)
                {
                        if(ip == G.ip)
                                continue;
                        if(e>=G.global.min_Mcode[lev][ip] && e <=G.global.max_Mcode[lev][ip])
                        {
                                        ipval[nip] = ip;
                                        nip++;
                        }
                }
        }
        sort(ipval,nip);
        nip = remove_dup(ipval,nip);
        return nip;
}

void find_ip()
{
	return;
}



void build_interp_conn(lev)
	int lev;
{
        int llev,ulev,ip,level;
        int e,i,j,k,p,proc[VPTS],n,nedge,ncenter;
        int d,eq_l,eq_u,n1,n2,n3,n4,i1,i2,i3,i4;
        int *size,**remote_eq;
        int **remote_node,*node_amount,*node_size,**remote_node_n, *amount;
        char *check;

        llev = lev;
        ulev = lev + 1;

        check = malloc(G.nno[ulev] * sizeof(char));
        for(i=0; i<G.nno[ulev]; i++)
                check[i] = 'n';

	//for local eqs
        G.interp[ulev].local_eq_n = realloc(G.interp[ulev].local_eq_n, G.neq[ulev] * sizeof(int));
	for(d=0; d<VPTS; d++)
	        G.interp[ulev].local_eq[d] = realloc(G.interp[ulev].local_eq[d], G.neq[ulev] * sizeof(int));
        for(i=0; i<G.neq[ulev]; i++)
                G.interp[ulev].local_eq_n[i] = 0;

	G.interp[ulev].nip_send = 0;
	G.interp[ulev].ip_send = realloc(G.interp[ulev].ip_send, G.np * sizeof(int));//the IP I send data to

        amount = malloc(G.np * sizeof(int));
        remote_eq = malloc(G.np * sizeof(int *));
        size = malloc(G.np * sizeof(int));
        for(ip=0; ip<G.np; ip++)
        {
                amount[ip] = 0;
                size[ip] = 2;
                remote_eq[ip] = malloc(size[ip] * sizeof(int));
        }

        remote_node = malloc(G.np * sizeof(int *));
        remote_node_n = malloc(G.np * sizeof(int *));
        node_amount = malloc(G.np * sizeof(int));
        node_size = malloc(G.np * sizeof(int));
        for(ip=0; ip<G.np; ip++)
        {
                node_amount[ip] = 0;
                node_size[ip] = 1;
                remote_node[ip] = malloc(node_size[ip] * sizeof(int));
                remote_node_n[ip] = malloc(node_size[ip] * sizeof(int));
        }

	for(e=0; e<G.nel[llev]; e++)//iter on lower level elements because 1 higher level node comes from 1 or more lower level nodes
	{
		level = Mcode_lev(G.Mcode[llev][e]);
		if(llev >= level)//this says that the element is a leaf node
		{
			//first, check if nodes exist at center of edges, 2 lower nodes -> 1 higher nodes
			for(p=0; p<VPTS; p++)
			{
				n1 = G.ien[llev][p][e];
				n2 = G.ien[llev][(p+1) % VPTS ][e];
				nedge = (n1 + n2)/2;
				j = node_index(nedge, ulev);
                                if(j >= 0)//the node in the center of the edge exists in the upper level. This is local interpolation
                                {
                                        if(check[j] == 'y')//only need to check one, although a node is shared among elements 
                                                continue;
                                        check[j] = 'y';
                                        for(d=0; d<NSD; d++)//prepare eqs
                                        {
                                                eq_u = G.id[ulev][d][j];
                                                if(eq_u >= 0)
                                                {
                                                        i1 = node_index(n1, llev);
                                                        i2 = node_index(n2, llev);
                                                        G.interp[ulev].local_eq[0][eq_u] = G.id[llev][d][i1];
                                                        G.interp[ulev].local_eq[1][eq_u] = G.id[llev][d][i2];
                                                        G.interp[ulev].local_eq_n[eq_u] = 2;
                                                }
                                        }
                                }
				//no matter the edge node exist locally or not, need to share it with other processors
                		n = what_ip_node_in(nedge, ulev, proc);
				//temporary save these nodes exist remotely in relevant processor, using dynamic memeory 
				for(k=0; k<n; k++)
				{
					ip = proc[k];

					//save nodes
					remote_node[ip][node_amount[ip]] = nedge;
					remote_node_n[ip][node_amount[ip]] = 2;
					node_amount[ip]++;
					if(node_amount[ip] >= node_size[ip])
					{
						node_size[ip] *= 2;
						remote_node[ip] = realloc(remote_node[ip], node_size[ip] * sizeof(int));
						remote_node_n[ip] = realloc(remote_node_n[ip], node_size[ip] * sizeof(int));
					}
					//save eqs
					for(d=0; d<NSD; d++)
					{
						i1 = node_index(n1, llev);
						i2 = node_index(n2, llev);

						remote_eq[ip][amount[ip]  ] = G.id[llev][d][i1];
						remote_eq[ip][amount[ip]+1] = G.id[llev][d][i2];

						amount[ip] += 2;

						if(amount[ip] >= size[ip])
						{
							size[ip] *= 2;
							remote_eq[ip] = realloc(remote_eq[ip], size[ip] * sizeof(int));
						}
					}

					//record the IP where I send these remote data to 
					if(value_in_array(ip, G.interp[ulev].ip_send, G.interp[ulev].nip_send) == -1)//if not recorded already
						G.interp[ulev].ip_send[ G.interp[ulev].nip_send++ ] = ip;
				}
			}

			//second, check center node: 4 -> 1
                        ncenter = 0;
                        for(p=0; p<VPTS; p++)
                                ncenter += G.ien[llev][p][e];
                        ncenter /= VPTS;
                        j = node_index(ncenter,ulev);
                        if(j >= 0)
                        {
                                if(check[j] == 'y')
                                        continue;
                                check[j] = 'y';
                                for(d=0; d<NSD; d++)
                                {
                                        eq_u = G.id[ulev][d][j];
                                        if(eq_u >= 0)
                                        {
                                                for(p=0; p<VPTS; p++)
                                                {
                                                        i = node_index(G.ien[llev][p][e], llev);
                                                        G.interp[ulev].local_eq[p][eq_u] = G.id[llev][d][i];//use [p] at index, not 1&2 for edge nodes
                                                }
                                                G.interp[ulev].local_eq_n[eq_u] = VPTS;
                                        }
                                }
                        }
			//again, need to share with other ip, not matter locally exits or not
			n1 = G.ien[llev][0][e];
                        n2 = G.ien[llev][1][e];
                        n3 = G.ien[llev][2][e];
                        n4 = G.ien[llev][3][e];
                	n = what_ip_node_in(ncenter, ulev, proc);
			for(k=0; k<n; k++)
                        {
                                ip = proc[k];

                                remote_node[ip][node_amount[ip]] = ncenter;
                                remote_node_n[ip][node_amount[ip]] = VPTS;
                                node_amount[ip]++;
                                if(node_amount[ip] >= node_size[ip])
                                {
                                        node_size[ip] *= 2;
                                        remote_node[ip] = realloc(remote_node[ip], node_size[ip] * sizeof(int));
                                        remote_node_n[ip] = realloc(remote_node_n[ip], node_size[ip] * sizeof(int));
                                }

                                for(d=0; d<NSD; d++)
                                {
                                        i1 = node_index(n1, llev);
                                        i2 = node_index(n2, llev);
                                        i3 = node_index(n3, llev);
                                        i4 = node_index(n4, llev);

                                        remote_eq[ip][amount[ip]  ] = G.id[llev][d][i1];
                                        remote_eq[ip][amount[ip]+1] = G.id[llev][d][i2];
                                        remote_eq[ip][amount[ip]+2] = G.id[llev][d][i3];
                                        remote_eq[ip][amount[ip]+3] = G.id[llev][d][i4];

                                        amount[ip] += VPTS;

                                        if(amount[ip] >= size[ip])
                                        {
                                                size[ip] *= 2;
                                                remote_eq[ip] = realloc(remote_eq[ip], size[ip] * sizeof(int));
                                        }
                                }

				if(value_in_array(ip, G.interp[ulev].ip_send, G.interp[ulev].nip_send) == -1)//if not recorded already
					G.interp[ulev].ip_send[ G.interp[ulev].nip_send++ ] = ip;
                        }
		}
	}

	//third, if upper mesh and lower mesh have the same nodes, this is the easiest
	//this is not from element, but from nodes directly
	for(i=0; i<G.nno[llev]; i++)
	{
		j = node_index(G.node[llev][i], ulev);
                if(j >= 0)
                {
                        if(check[j] == 'y')
                                continue;
                        check[j] = 'y';
                        for(d=0; d<NSD; d++)
                        {
                                eq_u = G.id[ulev][d][j];
                                if(eq_u >= 0)
                                {
                                        eq_l = G.id[llev][d][i];
					k = G.interp[ulev].local_eq_n[eq_u];
                                        G.interp[ulev].local_eq[k][eq_u] = eq_l;
					G.interp[ulev].local_eq_n[eq_u]++;
                                }
                        }
                }


		//again, need to share with other ip, not matter this node exists locally or not
                n = what_ip_node_in(G.node[llev][i], ulev, proc);
                for(k=0; k<n; k++)
                {
                        ip = proc[k];

                        remote_node[ip][node_amount[ip]] = G.node[llev][i];
                        remote_node_n[ip][node_amount[ip]] = 1;
                        node_amount[ip]++;
                        if(node_amount[ip] >= node_size[ip])
                        {
                                node_size[ip] *= 2;
                                remote_node[ip] = realloc(remote_node[ip], node_size[ip] * sizeof(int));
                                remote_node_n[ip] = realloc(remote_node_n[ip], node_size[ip] * sizeof(int));
                        }

                        for(d=0; d<NSD; d++)
                        {
                                remote_eq[ip][amount[ip]  ] = G.id[llev][d][i];

                                amount[ip] += 1;

                                if(amount[ip] >= size[ip])
                                {
                                        size[ip] *= 2;
                                        remote_eq[ip] = realloc(remote_eq[ip], size[ip] * sizeof(int));
                                }
                        }
			if(value_in_array(ip, G.interp[ulev].ip_send, G.interp[ulev].nip_send) == -1)//if not recorded already
				G.interp[ulev].ip_send[ G.interp[ulev].nip_send++ ] = ip;
                }
	}

	//reassign memeory
        for(ip=0; ip<G.np; ip++)
        {
		if(amount[ip] > 0)
		{
			remote_eq[ip] = realloc(remote_eq[ip], amount[ip] * sizeof(int));
			remote_node[ip] = realloc(remote_node[ip], node_amount[ip] * sizeof(int));
		}
        }

//	G.interp[ulev].ip_send = realloc(G.interp[ulev].ip_send, G.interp[ulev].nip_send);//uncomment this line will crash the code. I do not know why
	sort(G.interp[ulev].ip_send, G.interp[ulev].nip_send);//may not necessary, but I like to keep things organized. Because nip_send is small, this sorting does not take time

	//setup interp framework for sharing across levels
	G.interp[ulev].node_amount = realloc(G.interp[ulev].node_amount, G.interp[ulev].nip_send * sizeof(int));
        G.interp[ulev].eq_amount = realloc(G.interp[ulev].eq_amount, G.interp[ulev].nip_send * sizeof(int));

        for(i=0; i<G.interp[ulev].nip_send; i++)
        {
                ip = G.interp[ulev].ip_send[i];
                G.interp[ulev].node_amount[i] = node_amount[ip];
                G.interp[ulev].remote_node[i] = realloc(G.interp[ulev].remote_node[i], node_amount[ip] * sizeof(int));
                G.interp[ulev].remote_node_n[i] = realloc(G.interp[ulev].remote_node_n[i], node_amount[ip] * sizeof(int));
                for(j=0; j<node_amount[ip]; j++)
                {
                        G.interp[ulev].remote_node[i][j] = remote_node[ip][j];
                        G.interp[ulev].remote_node_n[i][j] = remote_node_n[ip][j];
                }

                G.interp[ulev].eq_amount[i] = amount[ip];
                G.interp[ulev].remote_eq[i] = realloc(G.interp[ulev].remote_eq[i], amount[ip] * sizeof(int));
                for(j=0; j<amount[ip]; j++)
                        G.interp[ulev].remote_eq[i][j] = remote_eq[ip][j];
        }
	//Above, i have found what processors should I send data to	
	//now, i need to find what processors I will recv data from	
	//the approch is to first gather all sending ip information to all processors, 
	//and then check which ip is sending data to me (which is the same ip i recv data from). 

        int *counts, *disp, *global_ip_send, global_nip_send;

	counts = malloc(G.np * sizeof(int));
	disp = malloc(G.np * sizeof(int));

	MPI_Allreduce(&G.interp[ulev].nip_send,&global_nip_send,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	global_ip_send = malloc(global_nip_send * sizeof(int));

	MPI_Allgather(&G.interp[ulev].nip_send, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
	disp[0] = 0;
	for(i=1; i<G.np; i++)
		disp[i] = disp[i-1] + counts[i-1];

	MPI_Allgatherv(G.interp[ulev].ip_send, G.interp[ulev].nip_send, MPI_INT, global_ip_send, counts, disp, MPI_INT, MPI_COMM_WORLD);

	G.interp[ulev].nip_recv = 0;
	G.interp[ulev].ip_recv = realloc(G.interp[ulev].ip_recv, G.np * sizeof(int));//the IP I send data to

	k = 0;
	for(i=0; i<G.np; i++)
	{
		for(j=0; j<counts[i]; j++)
		{
			if(global_ip_send[k] == G.ip)
			{
				G.interp[ulev].ip_recv[ G.interp[ulev].nip_recv++ ] = i;
			}
			k++;
		}
	}
	sort(G.interp[ulev].ip_recv, G.interp[ulev].nip_recv);


	free(global_ip_send);
	free(counts);
	free(disp);

        free(amount);
        free(size);
        for(ip=0; ip<G.np; ip++)
        {
                free(remote_eq[ip]);
                free(remote_node[ip]);
                free(remote_node_n[ip]);
        }
        free(remote_eq);
        free(remote_node);
        free(remote_node_n);

	free(node_amount);
	free(node_size);

	return;
}

void build_project_conn(lev)
	int lev;
{
	int i, j, k, llev, ulev, d, eq_u, eq_l, n, proc[VPTS], ip;
        int **remote_node,*node_amount,*node_size,**remote_node_n, *amount, *size, **remote_eq;

	llev = lev;
	ulev = lev + 1;
	//similar but much simpler approch than building interp conn,
	//only 1 node needs to considered,
	//the difference is that we need to prepare data for sending under the upper level,
	//because data is project from upper level to lower level

	//for local eqs
        G.project[llev].local_eq_n = realloc(G.project[llev].local_eq_n, G.neq[llev] * sizeof(int));
	G.project[llev].local_eq[0] = realloc(G.project[llev].local_eq[0], G.neq[llev] * sizeof(int));//0 means 1->1, only one node is projected
        for(i=0; i<G.neq[llev]; i++)
                G.project[llev].local_eq_n[i] = 0;


        amount = malloc(G.np * sizeof(int));
        remote_eq = malloc(G.np * sizeof(int *));
        size = malloc(G.np * sizeof(int));
        for(ip=0; ip<G.np; ip++)
        {
                amount[ip] = 0;
                size[ip] = 1;
                remote_eq[ip] = malloc(size[ip] * sizeof(int));
        }

        remote_node = malloc(G.np * sizeof(int *));
        remote_node_n = malloc(G.np * sizeof(int *));
        node_amount = malloc(G.np * sizeof(int));
        node_size = malloc(G.np * sizeof(int));
        for(ip=0; ip<G.np; ip++)
        {
                node_amount[ip] = 0;
                node_size[ip] = 1;
                remote_node[ip] = malloc(node_size[ip] * sizeof(int));
                remote_node_n[ip] = malloc(node_size[ip] * sizeof(int));
        }


	for(i=0; i<G.nno[ulev]; i++)
	{
		//first, local nodes exist in both upper and lower levels
		j = node_index(G.node[ulev][i], llev);
		if(j >= 0)
		{
                        for(d=0; d<NSD; d++)
                        {
                                eq_l = G.id[llev][d][j];
                                if(eq_l >= 0)
                                {
                                        eq_u = G.id[ulev][d][i];
                                        G.project[llev].local_eq[0][eq_l] = eq_u;
					G.project[llev].local_eq_n[eq_l] = 1;
                                }
                        }
		}

		//if a node in lower level also exist in upper level for other processors, 
		//we need to tell the other processors that they need to send information to this current node,
		//in order to get eq values at lower level from other upper level ips
		
		//first, check what other ips this node exists in the lower level
                n = what_ip_node_in(G.node[ulev][i], llev, proc);

                for(k=0; k<n; k++)
                {
                        ip = proc[k];

			//record the node and related eq information, which will be sent to these processors later
			//The idea is that after these other processors prepare related eq values, they need to send them back to assign values to the current local eqs
                        remote_node[ip][node_amount[ip]] = G.node[ulev][i];
                        remote_node_n[ip][node_amount[ip]] = 1;
                        node_amount[ip]++;
                        if(node_amount[ip] >= node_size[ip])
                        {
                                node_size[ip] *= 2;
                                remote_node[ip] = realloc(remote_node[ip], node_size[ip] * sizeof(int));
                                remote_node_n[ip] = realloc(remote_node_n[ip], node_size[ip] * sizeof(int));
                        }

                        for(d=0; d<NSD; d++)
                        {
                                remote_eq[ip][amount[ip]] = G.id[ulev][d][i];
				amount[ip] += 1;
                                if(amount[ip] >= size[ip])
                                {
                                        size[ip] *= 2;
                                        remote_eq[ip] = realloc(remote_eq[ip], size[ip] * sizeof(int));
                                }
                        }
                }
	}

        //reassign memeory
        for(ip=0; ip<G.np; ip++)
        {
                remote_eq[ip] = realloc(remote_eq[ip], amount[ip] * sizeof(int));
                remote_node[ip] = realloc(remote_node[ip], node_amount[ip] * sizeof(int));
        }

        //setup interp framework for sharing across levels

	//number of ips data should be sent to from upper to lower level in projection is the same as the number of ips reced from lower levels, and the opposite is also true
	G.project[llev].nip_send = G.interp[ulev].nip_recv;
	G.project[llev].ip_send = realloc(G.project[llev].ip_send, G.project[llev].nip_send * sizeof(int));
	for(i=0; i<G.project[llev].nip_send; i++)
		G.project[llev].ip_send[i] = G.interp[ulev].ip_recv[i];

	//the opposite is also true
	G.project[llev].nip_recv = G.interp[ulev].nip_send;
	G.project[llev].ip_recv = realloc(G.project[llev].ip_recv, G.project[llev].nip_recv * sizeof(int));
	for(i=0; i<G.project[llev].nip_recv; i++)
		G.project[llev].ip_recv[i] = G.interp[ulev].ip_send[i];

	//for data
        G.project[llev].node_amount = realloc(G.project[llev].node_amount, G.project[llev].nip_send * sizeof(int));
        G.project[llev].eq_amount = realloc(G.project[llev].eq_amount, G.project[llev].nip_send * sizeof(int));

        for(i=0; i<G.project[llev].nip_send; i++)
        {
                ip = G.project[llev].ip_send[i];
                G.project[llev].node_amount[i] = node_amount[ip];

                G.project[llev].remote_node[i] = realloc(G.project[llev].remote_node[i], node_amount[ip] * sizeof(int));
                G.project[llev].remote_node_n[i] = realloc(G.project[llev].remote_node_n[i], node_amount[ip] * sizeof(int));
                for(j=0; j<node_amount[ip]; j++)
                {
                        G.project[llev].remote_node[i][j] = remote_node[ip][j];
                        G.project[llev].remote_node_n[i][j] = remote_node_n[ip][j];
                }

                G.project[llev].eq_amount[i] = amount[ip];
                G.project[llev].remote_eq[i] = realloc(G.project[llev].remote_eq[i], amount[ip] * sizeof(int));
                for(j=0; j<amount[ip]; j++)
                        G.project[llev].remote_eq[i][j] = remote_eq[ip][j];
        }

        free(amount);
        free(size);
        for(ip=0; ip<G.np; ip++)
        {
                free(remote_eq[ip]);
                free(remote_node[ip]);
                free(remote_node_n[ip]);
        }
        free(remote_eq);
        free(remote_node);
        free(remote_node_n);

	free(node_amount);
	free(node_size);

	return;
}


void setup_project_comm_frame(lev)
	int lev;
{
	int i,l,u,d,k,j,node,eq,nlocal,nremote;
	int number;
	int ip,loc,size;
	char *checked;
	int *buff_i = NULL;
	int llev, ulev, ml, mu;
	MPI_Status stats;

	llev = lev;
	ulev = lev + 1;
	ml = G.neq[llev];
	mu = G.neq[ulev];

	checked = malloc(ml * sizeof(char));
	for(l=0; l<ml; l++) 
		checked[l] = 'n';


	G.project[llev].local_l = realloc(G.project[llev].local_l, ml * sizeof(int));
	G.project[llev].local_u = realloc(G.project[llev].local_u, ml * sizeof(int));
	nlocal = 0;

	for(l=0; l<ml; l++)
	{
		if(G.project[llev].local_eq_n[l] == 1)
		{
			u = G.project[llev].local_eq[0][l];

			G.project[llev].local_l[nlocal] = l;
			G.project[llev].local_u[nlocal] = u;
			nlocal++;
//			res_l[l] = res_u[u];
			checked[l] = 'y';
		}
	}
	G.project[llev].local_l = realloc(G.project[llev].local_l, nlocal * sizeof(int));
	G.project[llev].local_u = realloc(G.project[llev].local_u, nlocal * sizeof(int));
	G.project[llev].nlocal = nlocal;


	for(i=0; i<G.project[llev].nip_send; i++)
	{
		nremote = 0;
		size = 1;

		G.project[llev].send_l[i] = realloc(G.project[llev].send_l[i], size * sizeof(int));
		G.project[llev].send_u[i] = realloc(G.project[llev].send_u[i], size * sizeof(int));

		ip = G.project[llev].ip_send[i];
		number = G.project[llev].node_amount[i];

		loc = 0; 
		for(j=0; j<G.project[llev].node_amount[i]; j++)
		{
			for(d=0; d<NSD; d++)
                        {
                                eq = j * NSD + d;
                                u = G.project[llev].remote_eq[i][loc++];
				G.project[llev].send_l[i][nremote] = eq;
				G.project[llev].send_u[i][nremote] = u;
				nremote++;
				if(nremote >= size)
				{
					size *= 2;
					G.project[llev].send_l[i] = realloc(G.project[llev].send_l[i], size * sizeof(int));
					G.project[llev].send_u[i] = realloc(G.project[llev].send_u[i], size * sizeof(int));
				}

//				if(u >= 0)
//	                                buff_d[eq] = res_u[u];
//				else 
//					buff_d[eq] = 0.0;
                        }
		}
		MPI_Send(G.project[llev].remote_node[i],number,MPI_INT,ip,0,MPI_COMM_WORLD);

		G.project[llev].send_l[i] = realloc(G.project[llev].send_l[i], nremote * sizeof(int));
		G.project[llev].send_u[i] = realloc(G.project[llev].send_u[i], nremote * sizeof(int));
		G.project[llev].nsend[i] = nremote;
	}

	for(i=0; i<G.project[llev].nip_recv; i++)
        {
		nremote = 0;
		size = 1;
		G.project[llev].recv_l[i] = realloc(G.project[llev].recv_l[i], size * sizeof(int));
		G.project[llev].recv_u[i] = realloc(G.project[llev].recv_u[i], size * sizeof(int));

                ip = G.project[llev].ip_recv[i];

                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &number);
                buff_i = malloc(number * sizeof(int));
                MPI_Recv(buff_i,number,MPI_INT,ip,0,MPI_COMM_WORLD,&stats);

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
						G.project[llev].recv_l[i][nremote] = l;
						G.project[llev].recv_u[i][nremote] = d + k * NSD;
						nremote++;
						if(nremote >= size)
						{
							size *= 2;
							G.project[llev].recv_l[i] = realloc(G.project[llev].recv_l[i], size * sizeof(int));
							G.project[llev].recv_u[i] = realloc(G.project[llev].recv_u[i], size * sizeof(int));
						}

		//				res_l[l] = buff_d[d + k * NSD];
						checked[l] = 'y';
					}
                                }
                        }
                }

		G.project[llev].recv_l[i] = realloc(G.project[llev].recv_l[i], nremote * sizeof(int));
		G.project[llev].recv_u[i] = realloc(G.project[llev].recv_u[i], nremote * sizeof(int));
		G.project[llev].nrecv[i] = nremote;

		free(buff_i);
        }

	free(checked);

	return;
}

void setup_interp_comm_frame(lev)
	int lev;
{
        int i,j,k,d,eq,node,a;
        int u,l,mu,ml,llev,ulev,nlocal,nremote,size;
        char *checked;
        int ip;
        int loc;
        int number;
        int *buff_i = NULL;
        MPI_Status stats;


        llev = lev;
        ulev = lev + 1;
        ml = G.neq[llev];
        mu = G.neq[ulev];

        checked = malloc(mu * sizeof(char));
        for(u=0; u<mu; u++)
                checked[u] = 'n';


        G.interp[ulev].local_u = realloc(G.interp[ulev].local_u, mu * sizeof(int));
        G.interp[ulev].local_nlv = realloc(G.interp[ulev].local_nlv, mu * sizeof(int));
	for(d=0; d<VPTS; d++)
	{
	        G.interp[ulev].local_lv[d] = realloc(G.interp[ulev].local_lv[d], mu * sizeof(int));
	}
        nlocal = 0;

        for(u=0; u<mu; u++)
        {
                if(G.interp[ulev].local_eq_n[u] > 0)
                {
			G.interp[ulev].local_u[nlocal] = u;
			d = 0;
                        for(i=0; i<G.interp[ulev].local_eq_n[u]; i++)
                        {
                                l = G.interp[ulev].local_eq[i][u];
                                if(l != -1)//boundary nodes
                                {
					G.interp[ulev].local_lv[d++][nlocal] = l;
                                }
                        }
			G.interp[ulev].local_nlv[nlocal] = d;

			nlocal++;
                        checked[u] = 'y';
                }
        }

        G.interp[ulev].local_nlv = realloc(G.interp[ulev].local_nlv, nlocal * sizeof(int));
	for(d=0; d<VPTS; d++)
        	G.interp[ulev].local_lv[d] = realloc(G.interp[ulev].local_lv[d], nlocal * sizeof(int));

        G.interp[ulev].local_u = realloc(G.interp[ulev].local_u, nlocal * sizeof(int));
        G.interp[ulev].nlocal = nlocal;


        for(i=0; i<G.interp[ulev].nip_send; i++)
        {
                ip = G.interp[ulev].ip_send[i];
                number = G.interp[ulev].node_amount[i];
                MPI_Send(G.interp[ulev].remote_node[i],number,MPI_INT,ip,0,MPI_COMM_WORLD);

		// the following section is commented out, because it is complex and does not increase the speed of calculation neither
		/*
		nremote = 0;
		size = G.interp[ulev].node_amount[i] * NSD;
		G.interp[ulev].send_u[i] = realloc(G.interp[ulev].send_u[i], size * sizeof(int));
		G.interp[ulev].send_nlv[i] = realloc(G.interp[ulev].send_nlv[i], size * sizeof(int));
		for(a=0; a<VPTS; a++)
			G.interp[ulev].send_lv[i][a] = realloc(G.interp[ulev].send_lv[i][a], size * sizeof(int));


                loc = 0;
                for(j=0; j<G.interp[ulev].node_amount[i]; j++)
                {
                        for(d=0; d<NSD; d++)
                        {
                                eq = j * NSD + d;
//				G.interp[ulev].send_u[i][nremote] = eq;
				a = 0;
                                for(k=0; k<G.interp[ulev].remote_node_n[i][j]; k++)
                                {
                                        l = G.interp[ulev].remote_eq[i][loc++];
                                        if(l != -1)
                                        {
//						G.interp[ulev].send_lv[i][a++][nremote] = l;
                                        }
                                }

				G.interp[ulev].send_nlv[i][nremote] = a;
				nremote++;
                        }
                }

			*/

		/*
		size = nremote;
		G.interp[ulev].send_u[i] = realloc(G.interp[ulev].send_u[i], size * sizeof(int));
		G.interp[ulev].send_nlv[i] = realloc(G.interp[ulev].send_nlv[i], size * sizeof(int));
		for(a=0; a<VPTS; a++)
			G.interp[ulev].send_lv[i][a] = realloc(G.interp[ulev].send_lv[i][a], size * sizeof(int));
		G.interp[ulev].nsend[i] = nremote;
		*/

        }

        for(i=0; i<G.interp[ulev].nip_recv; i++)
        {
                ip = G.interp[ulev].ip_recv[i];
                MPI_Probe(ip, 0, MPI_COMM_WORLD, &stats);
                MPI_Get_count(&stats, MPI_INT, &number);
                buff_i = malloc(number * sizeof(int));
                MPI_Recv(buff_i,number,MPI_INT,ip,0,MPI_COMM_WORLD,&stats);

		nremote = 0;
		size = number * NSD;
		G.interp[ulev].recv_u[i] = realloc(G.interp[ulev].recv_u[i], size * sizeof(int));
		G.interp[ulev].recv_l[i] = realloc(G.interp[ulev].recv_l[i], size * sizeof(int));

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
						G.interp[ulev].recv_u[i][nremote] = u;
						G.interp[ulev].recv_l[i][nremote] = d + k * NSD;
						nremote++;
                                                checked[u] = 'y';
                                        }
                                }
                        }
                }

		G.interp[ulev].recv_u[i] = realloc(G.interp[ulev].recv_u[i], nremote * sizeof(int));
		G.interp[ulev].recv_l[i] = realloc(G.interp[ulev].recv_l[i], nremote * sizeof(int));

                free(buff_i);
        }
	free(checked);
	return;
}

void find_nshare_for_eqs(lev)
	int lev;
{
	int i, j, m;
	int neq = G.neq_m[lev];
	int nsh, **seq, *nseq;

	G.smpi[lev].veq_nsh = realloc(G.smpi[lev].veq_nsh, neq * sizeof(int));

        for(i=0; i<neq; i++)
		G.smpi[lev].veq_nsh[i] = 1;

        nsh = G.smpi[lev].nshareip;
        nseq = G.smpi[lev].n_veq;
        seq = G.smpi[lev].veq;

        for(i=0; i<neq; i++)
        {
                for(j=0; j<nsh; j++)
                {
                        if(value_in_sorted_array(i, seq[j], nseq[j]) >= 0)
                        {
                                G.smpi[lev].veq_nsh[i] += 1;
                        }
                }
        }

	return;
}

void find_nshare_for_Teqs(lev)
	int lev;
{
	int i, j, m;
	int neq = G.heat[lev].neq_m;
	int nsh, **seq, *nseq;

	G.smpi[lev].Teq_nsh = realloc(G.smpi[lev].Teq_nsh, neq * sizeof(int));

        for(i=0; i<neq; i++)
		G.smpi[lev].Teq_nsh[i] = 1;

        nsh = G.smpi[lev].nshareip;
        nseq = G.smpi[lev].n_Teq;
        seq = G.smpi[lev].Teq;

        for(i=0; i<neq; i++)
        {
                for(j=0; j<nsh; j++)
                {
                        if(value_in_sorted_array(i, seq[j], nseq[j]) >= 0)
                        {
                                G.smpi[lev].Teq_nsh[i] += 1;
                        }
                }
        }

	return;
}
