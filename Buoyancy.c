#include "global_variables.h"

void buoyancy_force()
{
        int i,lev;
	double *avg;

	avg = malloc(G.gnoz * sizeof(double));

        for(lev=G.min_level; lev<G.max_level; lev++)
	{
                G.buoyancy[lev] = realloc(G.buoyancy[lev], G.nno[lev] * sizeof(double));
		for(i=0;i<G.nno[lev];i++)
			G.buoyancy[lev][i] = G.Ra * G.T[lev][i];

		get_node_horiz_avg(G.buoyancy[lev], lev, avg);
		remove_node_horiz_avg(G.buoyancy[lev],lev,avg);
	}

	free(avg);
	return;
}
