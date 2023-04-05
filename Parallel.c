#include "global_variables.h"

void distribute_Mcode_among_processors()
{
	int lev;
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		distribute_Mcode_evenly_blocked_comm(lev);
		range_of_Mcode(lev);
	}

	return;
}

void distribute_Mcode_evenly(lev)
	int lev;
{
	int nel, *nelv, ip, e;
	int *enelv;
	int ave_nel,nel_remain;
	int *nsend, *nrecv, *size;
	int **send, **recv;
	int *tmp_Mcode, ntmp;
	MPI_Status stats;
	MPI_Request request[1000];
	int nir;

	//nel: total number of elements
	//nelv: accumulated number of element before distributed
	//enelv: accumulated number of element after distributed


	nelv = malloc(G.np * sizeof(int));
	MPI_Allgather(&G.nel[lev], 1, MPI_INT, nelv, 1, MPI_INT, MPI_COMM_WORLD); // send nel for each element to all elements

	for(ip=1; ip<G.np; ip++)
		nelv[ip] += nelv[ip - 1];//add previous nel up

	nel = nelv[G.np - 1];

	//update nel for each ip, evenly distribute
	ave_nel = (nel - nel % G.np) / G.np;
	nel_remain = nel - ave_nel * G.np;

	enelv = malloc(G.np * sizeof(int));
	for(ip = 0; ip<G.np; ip++)
	{
		enelv[ip] = ave_nel;
		if(ip < nel_remain)
			enelv[ip] += 1;
	}

	ntmp = 0;
	tmp_Mcode = malloc(enelv[G.ip] * sizeof(int));

	for(ip=1; ip<G.np; ip++)
		enelv[ip] += enelv[ip - 1];//add previous nel up

	//send and recv Mcodes
	//for convinience, I create memory for all processors. This does not waster much memory as nsend[ip] will be zero if there is no connection with ip.
	nsend = malloc(G.np * sizeof(int));
	size = malloc(G.np * sizeof(int));// for dynamic memory allocattion purpose only
	send = malloc(G.np * sizeof(int *));
	for(ip = 0; ip<G.np; ip++)
	{
		nsend[ip] = 0;
		size[ip] = 1;
		send[ip] = malloc(size[ip] * sizeof(int));
	}
	//prepare number of data and ip to send to
	for(e=0; e<G.nel[lev]; e++)//check every e within this ip
	{
		ip = find_element_ip(e + nelv[G.ip - 1], enelv);
		if(ip != G.ip) // does not belong to the current ip
		{
			send[ip][nsend[ip]] = G.Mcode[lev][e];
			nsend[ip]++;
			if(nsend[ip] >= size[ip])
			{
				size[ip] *= 2;
				send[ip] = realloc(send[ip], size[ip] * sizeof(int));
			}
		}
		else
		{
			tmp_Mcode[ntmp] = G.Mcode[lev][e];
			ntmp++;
		}
	}
	
	//send
	nir = 0;
	for(ip=0; ip<G.np; ip++)
	{
		if(nsend[ip] > 0)
		{
			MPI_Isend(send[ip], nsend[ip], MPI_INT, ip, 0, MPI_COMM_WORLD, &request[nir++]);
		//	printf("send %d | %d %d\n",G.ip, ip, nsend[ip]);
		}
	}

	//prepare number of data and ip to recv from
	nrecv = malloc(G.np * sizeof(int));
	recv = malloc(G.np * sizeof(int *));
	for(ip=0; ip<G.np; ip++)
	{
		nrecv[ip] = 0;
		recv[ip] = malloc(sizeof(int));
	}

	for(e=enelv[G.ip - 1]; e<enelv[G.ip]; e++)//check every e within this ip
	{
		ip = find_element_ip(e, nelv);
		if(ip != G.ip) // does not belong to the current ip
		{
			nrecv[ip]++;
		}
	}

	//recv
	for(ip=0; ip<G.np; ip++)
	{
		if(nrecv[ip] > 0)
		{
			recv[ip] = malloc(nrecv[ip] * sizeof(int));
			MPI_Irecv(recv[ip],nrecv[ip],MPI_INT,ip,0,MPI_COMM_WORLD,&request[nir++]);
//			printf("recv %d | %d %d\n",G.ip, ip, nrecv[ip]);

			//append recv Mcode
			for(e=0; e<nrecv[ip]; e++)
			{
				tmp_Mcode[ntmp] = recv[ip][e];
				ntmp++;
			}
		}
	}

	MPI_Waitall(nir, request, &stats);

	//update Mcode
	sort(tmp_Mcode, ntmp);
	G.nel[lev] = ntmp;
	G.Mcode[lev] = realloc(G.Mcode[lev], G.nel[lev] * sizeof(int));
	for(e=0; e<G.nel[lev]; e++)
		G.Mcode[lev][e] = tmp_Mcode[e];

	for(ip=0; ip<G.np; ip++)
	{
		free(send[ip]);
		free(recv[ip]);
	}


	free(send);
	free(recv);
	free(nsend);
	free(nrecv);

	free(tmp_Mcode);
	free(size);
	free(nelv);
	free(enelv);

	return;
}

void distribute_Mcode_evenly_blocked_comm(lev)
	int lev;
{
	int nel, *nelv, ip, e;
	int *enelv;
	int ave_nel,nel_remain;
	int *nsend, *nrecv, *size;
	int **send, **recv;
	int *tmp_Mcode, ntmp;
	MPI_Status stats;

	//nel: total number of elements
	//nelv: accumulated number of element before distributed
	//enelv: accumulated number of element after distributed


	nelv = malloc(G.np * sizeof(int));
	MPI_Allgather(&G.nel[lev], 1, MPI_INT, nelv, 1, MPI_INT, MPI_COMM_WORLD); // send nel for each element to all elements

	for(ip=1; ip<G.np; ip++)
		nelv[ip] += nelv[ip - 1];//add previous nel up

	nel = nelv[G.np - 1];

	//update nel for each ip, evenly distribute
	ave_nel = (nel - nel % G.np) / G.np;
	nel_remain = nel - ave_nel * G.np;

	enelv = malloc(G.np * sizeof(int));
	for(ip = 0; ip<G.np; ip++)
	{
		enelv[ip] = ave_nel;
		if(ip < nel_remain)
			enelv[ip] += 1;
	}

	ntmp = 0;
	tmp_Mcode = malloc(enelv[G.ip] * sizeof(int));

	for(ip=1; ip<G.np; ip++)
		enelv[ip] += enelv[ip - 1];//add previous nel up

	//send and recv Mcodes
	//for convinience, I create memory for all processors. This does not waster much memory as nsend[ip] will be zero if there is no connection with ip.
	nsend = malloc(G.np * sizeof(int));
	size = malloc(G.np * sizeof(int));// for dynamic memory allocattion purpose only
	send = malloc(G.np * sizeof(int *));
	for(ip = 0; ip<G.np; ip++)
	{
		nsend[ip] = 0;
		size[ip] = 1;
		send[ip] = malloc(size[ip] * sizeof(int));
	}
	//prepare number of data and ip to send to
	for(e=0; e<G.nel[lev]; e++)//check every e within this ip
	{
		ip = find_element_ip(e + nelv[G.ip - 1], enelv);
		if(ip != G.ip) // does not belong to the current ip
		{
			send[ip][nsend[ip]] = G.Mcode[lev][e];
			nsend[ip]++;
			if(nsend[ip] >= size[ip])
			{
				size[ip] *= 2;
				send[ip] = realloc(send[ip], size[ip] * sizeof(int));
			}
		}
		else
		{
			tmp_Mcode[ntmp] = G.Mcode[lev][e];
			ntmp++;
		}
	}
	
	//send
	for(ip=0; ip<G.np; ip++)
	{
		if(nsend[ip] > 0)
		{
			MPI_Send(send[ip], nsend[ip], MPI_INT, ip, 0, MPI_COMM_WORLD);
		//	printf("send %d | %d %d\n",G.ip, ip, nsend[ip]);
		}
	}

	//prepare number of data and ip to recv from
	nrecv = malloc(G.np * sizeof(int));
	recv = malloc(G.np * sizeof(int *));
	for(ip=0; ip<G.np; ip++)
	{
		nrecv[ip] = 0;
		recv[ip] = malloc(sizeof(int));
	}

	for(e=enelv[G.ip - 1]; e<enelv[G.ip]; e++)//check every e within this ip
	{
		ip = find_element_ip(e, nelv);
		if(ip != G.ip) // does not belong to the current ip
		{
			nrecv[ip]++;
		}
	}

	//recv
	for(ip=0; ip<G.np; ip++)
	{
		if(nrecv[ip] > 0)
		{
			recv[ip] = malloc(nrecv[ip] * sizeof(int));
			MPI_Recv(recv[ip],nrecv[ip],MPI_INT,ip,0,MPI_COMM_WORLD,&stats);
//			printf("recv %d | %d %d\n",G.ip, ip, nrecv[ip]);

			//append recv Mcode
			for(e=0; e<nrecv[ip]; e++)
			{
				tmp_Mcode[ntmp] = recv[ip][e];
				ntmp++;
			}
		}
	}

	//update Mcode
	sort(tmp_Mcode, ntmp);
	G.nel[lev] = ntmp;
	G.Mcode[lev] = realloc(G.Mcode[lev], G.nel[lev] * sizeof(int));
	for(e=0; e<G.nel[lev]; e++)
		G.Mcode[lev][e] = tmp_Mcode[e];

	for(ip=0; ip<G.np; ip++)
	{
		free(send[ip]);
		free(recv[ip]);
	}
	free(send);
	free(recv);
	free(nsend);
	free(nrecv);
	free(tmp_Mcode);
	free(size);
	free(nelv);
	free(enelv);

	return;
}


int find_element_ip(e, nel)
	int e, *nel;
{
	int ip;

	if(e < nel[0]) return 0;

	for(ip = 1; ip<G.np; ip++)
	{
		if(e >= nel[ip - 1] && e < nel[ip])
			return ip;
	}
}

void range_of_Mcode(lev)
	int lev;
{
	int Mcode,level,mul;

	Mcode = Mcode_loc(G.Mcode[lev][0]);
	MPI_Allgather(&Mcode, 1, MPI_INT, G.global.min_Mcode[lev], 1, MPI_INT, MPI_COMM_WORLD);

	Mcode = Mcode_loc(G.Mcode[lev][G.nel[lev]-1]);
	level = Mcode_lev(G.Mcode[lev][G.nel[lev]-1]);

	mul = (int)pow(2,level);
	mul = G.max_nel[G.max_level - 1]/mul/mul;
	Mcode += mul - 1;

	MPI_Allgather(&Mcode, 1, MPI_INT, G.global.max_Mcode[lev], 1, MPI_INT, MPI_COMM_WORLD);
	return;
}

void terminate(void)
{
        MPI_Finalize();
        fflush(stderr);
        fflush(stdout);
        exit(0);
        return;
}

void get_node_horiz_avg(data,lev,avg)
	double *data;
	double *avg;
	int lev;
{
	int e,ix[2],iz[2],node,i,k,a,index, index_arr[4];
	double xi,eta;
	double *global_avg;

	global_avg = malloc(G.gnoz * sizeof(double));


	for(i=0; i<G.gnoz; i++)
		avg[i] = 0.0;

	for(e=0; e<G.nel[lev]; e++)
	{
		for(a=0; a<VPTS; a++)
			index_arr[a] = node_index(G.ien[lev][a][e],lev);

		node = G.ien[lev][0][e];
                ix[0] = node / G.gnoz;
                iz[0] = node - ix[0] * G.gnoz;

		node = G.ien[lev][2][e];
                ix[1] = node / G.gnoz;
                iz[1] = node - ix[1] * G.gnoz;

		for(i=ix[0]; i<ix[1]; i++)
		{
			xi = (2 * i - ix[0] - ix[1])/(ix[1] - ix[0]);
			for(k=iz[0]; k<iz[1]; k++)
			{
				eta = (2 * k - iz[0] - iz[1])/(iz[1] - iz[0]);
				for(a=0; a<VPTS; a++)
				{
					index = index_arr[a];
					avg[k] += Na(a,xi,eta) * data[index];
				}
			}
		}

		if(iz[1] == G.gnoz -1)
		{
			k = iz[1];
			for(i=ix[0]; i<ix[1]; i++)
			{
				xi = (2 * i - ix[0] - ix[1])/(ix[1] - ix[0]);
				eta = 1.0;
				for(a=0; a<VPTS; a++)
				{
					index = index_arr[a];
					avg[k] += Na(a,xi,eta) * data[index];
				}
			}
		}

		if(ix[1] == G.gnox -1)
		{
			xi = 1.0;
			for(k=iz[0]; k<iz[1]; k++)
			{
				eta = (2 * k - iz[0] - iz[1])/(iz[1] - iz[0]);
				for(a=0; a<VPTS; a++)
				{
					index = index_arr[a];
					avg[k] += Na(a,xi,eta) * data[index];
				}
			}
		}

		if(node == G.gnox * G.gnoz - 1)
			avg[G.gnoz - 1] += data[G.nno[lev]-1];
	}

	MPI_Allreduce(avg,global_avg,G.gnoz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	for(i=0; i<G.gnoz; i++)
		avg[i] = global_avg[i] / G.gnox;

	free(global_avg);

	return;
}

void remove_node_horiz_avg(data,lev,avg)
	double *data,*avg;
	int lev;
{
	int i,node,ix,iz;

	for(i=0; i<G.nno[lev]; i++)
	{
		node = G.node[lev][i];
                ix = node / G.gnoz;
                iz = node - ix * G.gnoz;
		data[i] -= avg[iz];
	}

	return;
}


