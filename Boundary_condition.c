#include "global_variables.h"

void boundary_condition()
{
	temperature_bc();
	velocity_bc();
	return;
}

void temperature_bc()
{
	int lev,j;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.bc.T_bc[lev] = realloc(G.bc.T_bc[lev], G.nno[lev] * sizeof(char));
		G.bc.T_bc_val[lev] = realloc(G.bc.T_bc_val[lev], G.nno[lev] * sizeof(double));
	}

	for(lev=G.min_level; lev<G.max_level; lev++)
	for(j=0; j<G.nno[lev]; j++)
	{
		G.bc.T_bc[lev][j] = 'N';
		G.bc.T_bc_val[lev][j] = 0.0;
	}

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(j=0; j<G.nno[lev]; j++)
		{
			if((G.node[lev][j] + 1) % G.gnoz == 0)//top
			{
				G.bc.T_bc[lev][j] = G.bc.top_T_bc;
				G.bc.T_bc_val[lev][j] = G.bc.top_T_bc_val;
			}
			else if(G.node[lev][j] % G.gnoz == 0)//bottom
			{
				G.bc.T_bc[lev][j] = G.bc.bot_T_bc;
				G.bc.T_bc_val[lev][j] = G.bc.bot_T_bc_val;
			}
			else if(G.node[lev][j] < G.gnoz)//left
			{
				if(G.bc.T_bc[lev][j] != 'T')
					G.bc.T_bc[lev][j] = G.bc.lef_T_bc;
				G.bc.T_bc_val[lev][j] = G.bc.lef_T_bc_val;
			}
			else if(G.node[lev][j] > (G.gnox -1) * G.gnoz - 1)//right
			{
				if(G.bc.T_bc[lev][j] != 'T')
					G.bc.T_bc[lev][j] = G.bc.rig_T_bc;
				G.bc.T_bc_val[lev][j] = G.bc.rig_T_bc_val;
			}
		}
	}
	return;
}

void velocity_bc()
{
	int i,lev,j;

        for(lev=G.min_level; lev<G.max_level; lev++)
        {
		for(i=0; i<NSD; i++)
		{
			G.bc.v_bc[lev][i] = realloc(G.bc.v_bc[lev][i], G.nno[lev] * sizeof(char));
			G.bc.v_bc_val[lev][i] = realloc(G.bc.v_bc_val[lev][i], G.nno[lev] * sizeof(double));
		}
	}

	for(lev=G.min_level; lev<G.max_level; lev++)
	for(i=0; i<NSD; i++)
	for(j=0; j<G.nno[lev]; j++)
	{
		G.bc.v_bc[lev][i][j] = 'N';
		G.bc.v_bc_val[lev][i][j] = 0.0;
	}

	for(lev=G.min_level; lev<G.max_level; lev++)
	for(i=0; i<NSD; i++)
	{
		for(j=0; j<G.nno[lev]; j++)
		{
			if((G.node[lev][j] + 1) % G.gnoz == 0)//top
			{
				G.bc.v_bc[lev][i][j] = G.bc.top_v_bc[i];
				G.bc.v_bc_val[lev][i][j] = G.bc.top_v_bc_val[i];
			}
			else if(G.node[lev][j] % G.gnoz == 0)//bottom
			{
				G.bc.v_bc[lev][i][j] = G.bc.bot_v_bc[i];
				G.bc.v_bc_val[lev][i][j] = G.bc.bot_v_bc_val[i];
			}

			if(G.node[lev][j] < G.gnoz)//left
			{
				if(G.bc.v_bc[lev][i][j] != 'V')
					G.bc.v_bc[lev][i][j] = G.bc.lef_v_bc[i];
				G.bc.v_bc_val[lev][i][j] = G.bc.lef_v_bc_val[i];
			}
			else if(G.node[lev][j] > (G.gnox -1) * G.gnoz - 1)//right
			{
				if(G.bc.v_bc[lev][i][j] != 'V')
					G.bc.v_bc[lev][i][j] = G.bc.rig_v_bc[i];
				G.bc.v_bc_val[lev][i][j] = G.bc.rig_v_bc_val[i];
			}
		}
	}
	return;
}
