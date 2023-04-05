#include "global_variables.h"

void construct_ien()
{
        int lev,e,a,b,step,loc,level,ipx,ipz,ip,code,node,gnoz;

        for(lev=G.min_level; lev<G.max_level; lev++)
		for(a=0; a<NEN; a++)
                	G.ien[lev][a] = realloc(G.ien[lev][a], G.nel[lev]*sizeof(int));

        gnoz = G.gnoz;

        for(lev=G.min_level; lev<G.max_level; lev++)
        {
                for(e=0; e<G.nel[lev]; e++)
                {
			level = Mcode_lev(G.Mcode[lev][e]);
			loc = Mcode_loc(G.Mcode[lev][e]);


			ip = loc/G.max_nel[G.max_level-1];
			ipx = ip/G.npz;
			ipz =  ip - ipx*G.npz;

			code = loc - ip * G.max_nel[G.max_level-1];
			demorton_no_lev(code,&a,&b);

			a += G.max_elx[G.max_level - 1] * ipx;
			b += G.max_elz[G.max_level - 1] * ipz;
			node = b + a * gnoz;

                	step = (int) pow(2, G.max_level - 1 - level);

                        G.ien[lev][0][e] = node;
                        G.ien[lev][1][e] = G.ien[lev][0][e] + step;
                        G.ien[lev][2][e] = G.ien[lev][1][e] + step * gnoz;
                        G.ien[lev][3][e] = G.ien[lev][0][e] + step * gnoz;
                }
        }
        return;
}

void node_index_from_ien()
{
        int lev,e,i,j;

        for(lev=G.min_level; lev<G.max_level; lev++)
        {
                G.nno[lev] = VPTS * G.nel[lev];
                G.node[lev] = realloc(G.node[lev], G.nno[lev] * sizeof(int));

                j = 0;
                for(e=0; e<G.nel[lev]; e++)
                for(i=0; i<VPTS; i++)
                        G.node[lev][j++] = G.ien[lev][i][e];

                sort(G.node[lev],G.nno[lev]);

                G.nno[lev] = remove_dup(G.node[lev],G.nno[lev]);
                G.node[lev] = realloc(G.node[lev], G.nno[lev] * sizeof(int));
        }

	if(G.control.use_hash)
		node_hashtable();

        return;
}
void construct_id()
{
	int lev;
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		construct_v_id(lev);
		construct_T_id(lev);
	}
	return;
}

void construct_v_id(lev)
	int lev;
{
	int i,A,P;

	for(i=0; i<NSD; i++)
		G.id[lev][i] = realloc(G.id[lev][i], G.nno[lev] * sizeof(int));


	G.veq_node[lev] = realloc(G.veq_node[lev], G.nno[lev] * 2 * sizeof(int));
	G.veq_dof[lev] = realloc(G.veq_dof[lev], G.nno[lev] * 2 * sizeof(int));

	P = 0;
	for(A=0; A<G.nno[lev]; A++)
	for(i=0; i<NSD; i++)
	{
		if(G.bc.v_bc[lev][i][A] == 'V')
		{
			G.id[lev][i][A] = -1;
		}
		else
		{
			G.id[lev][i][A] = P;
			G.veq_node[lev][P] = G.node[lev][A];
			G.veq_dof[lev][P] = i;
			P++;
		}
	}

	G.neq[lev] = P;
	G.veq_node[lev] = realloc(G.veq_node[lev], G.neq[lev] * sizeof(int));
	G.veq_dof[lev] = realloc(G.veq_dof[lev], G.neq[lev] * sizeof(int));

	G.neq_m[lev] = P;

	return;
}

void construct_T_id(lev)
	int lev;
{
	int P, A;

	G.heat[lev].id = realloc(G.heat[lev].id, G.nno[lev] * sizeof(int));
	G.heat[lev].eq_node = realloc(G.heat[lev].eq_node, G.nno[lev] * sizeof(int));

	P=0;
	for(A=0; A<G.nno[lev]; A++)
	{
		if(G.bc.T_bc[lev][A]=='T')
		{
			G.heat[lev].id[A] = -1;
		}
		else
		{
			G.heat[lev].id[A] = P;
			G.heat[lev].eq_node[P] = G.node[lev][A];
			P++;
		}
	}
	G.heat[lev].neq = P;
	G.heat[lev].eq_node = realloc(G.heat[lev].eq_node, G.heat[lev].neq * sizeof(int));
	G.heat[lev].neq_m = P;

	return;
}

void construct_lm()
{
	int lev;
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		construct_v_lm(lev);
		construct_T_lm(lev);
	}
	return;
}

void construct_v_lm(lev)
	int lev;
{
	int i,a,node,e;

	for(i=0; i<NSD; i++)
	{
		for(a=0; a<NEN; a++)
			G.lm[lev][i][a] = realloc(G.lm[lev][i][a], G.nel[lev] * sizeof(int));
	}

        for(i=0; i<NSD; i++)
	for(a=0; a<NEN; a++)
        for(e=0; e<G.nel[lev]; e++)
	{
		node = node_index(G.ien[lev][a][e],lev);
		G.lm[lev][i][a][e]=G.id[lev][i][node];
	}

	return;
}


void construct_T_lm(lev)
	int lev;
{
	int a,node,e;

	for(a=0; a<NEN; a++)
		G.heat[lev].lm[a] = realloc(G.heat[lev].lm[a], G.nel[lev] * sizeof(int));

        for(a=0; a<NEN; a++)
        for(e=0; e<G.nel[lev]; e++)
	{
		node = node_index(G.ien[lev][a][e], lev);
                G.heat[lev].lm[a][e]=G.heat[lev].id[node];
	}

	return;
}

void element_area()
{
	int e, lev, node1, node2;
	double dx, dz;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.eco[lev].area = realloc(G.eco[lev].area, G.nel[lev] * sizeof(double));
		for(e=0; e<G.nel[lev]; e++)
		{
			node1 = node_index(G.ien[lev][0][e], lev);
			node2 = node_index(G.ien[lev][2][e], lev);

			dx = G.X[lev][0][node2] - G.X[lev][0][node1];
			dz = G.X[lev][1][node2] - G.X[lev][1][node1];

			G.eco[lev].area[e] = dx * dz;
		}
	}
	return;
}

void what_elements_node_in()
{
	int lev, i, a, e;
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.p2n[lev].nel = realloc(G.p2n[lev].nel, G.nno[lev] * sizeof(int));
		for(a=0; a<NEN; a++)
		{
			G.p2n[lev].el[a] = realloc(G.p2n[lev].el[a], G.nno[lev] * sizeof(int));
		}

		for(i=0; i<G.nno[lev]; i++)
			G.p2n[lev].nel[i] = 0;

		for(e=0; e<G.nel[lev]; e++)
		{
			for(a=0; a<VPTS; a++)
			{
				i = node_index(G.ien[lev][a][e], lev);
				G.p2n[lev].el[G.p2n[lev].nel[i]++ ][i] = e;
			}
		}
	}
	return;
}
