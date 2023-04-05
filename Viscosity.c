#include "global_variables.h"

void viscosity()
{
	int i,lev;
	double z;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.vis[lev] = realloc(G.vis[lev], G.nno[lev] * sizeof(double));
		for(i=0; i< G.nno[lev]; i++)
		{
			G.vis[lev][i]=exp(G.viscosity.A * (G.viscosity.refT - G.T[lev][i]));
			z = G.X[lev][1][i];
			if(z < 0.77)
				G.vis[lev][i] *= G.viscosity.lm_jump;
		}

		if(G.control.grain_damage == 1)
			visc_from_grain_size(lev);
	}

	return;
}

