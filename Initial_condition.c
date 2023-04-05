#include "global_variables.h"

void initial_temperature()
{
	int lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
		G.T[lev] = realloc(G.T[lev], G.nno[lev] * sizeof(double));

	G.model_time = 0.0;

	switch (G.control.T_ini)
	{
		case 0:
			conductive_T();
			break;
		case 1:
			periodic_T();
			break;
		case 2:
			random_T_perturbation();
			break;
		case 3:
			periodic_T_with_random_perturb();
			break;
		case -1:
			read_T_from_file();
			break;
		default:
			fprintf(stderr,"Initial temperature is not set\n");
			terminate();
	}

	apply_T_bc_values();

	return;
}

void random_T_perturbation()
{
        double random;
        int i, lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(i=0; i<G.nno[lev]; i++)
		{
			random=(1.0*rand())/(1.0*RAND_MAX);
			G.T[lev][i] = G.control.T_mantle + G.control.T_mantle_perturb * random;
		}
	}
        return;
}

void conductive_T()
{
        double x,z;
        int i,lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(i=0; i<G.nno[lev]; i++)
		{
			x = G.X[lev][0][i];
			z = G.X[lev][1][i];
                        G.T[lev][i] = G.global.height - z 
                             + G.control.T_mantle_perturb
                             *sin(M_PI * (G.global.height - z)/G.global.height)
                             *cos(G.control.T_wave_number*M_PI*x/G.global.length);
		}
	}

        return;
}


void periodic_T()
{
        double x,z;
        int i,lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(i=0; i<G.nno[lev]; i++)
		{
			x = G.X[lev][0][i];
			z = G.X[lev][1][i];
                        G.T[lev][i] = G.control.T_mantle
                             + G.control.T_mantle_perturb
                             *sin(M_PI * (G.global.height - z)/G.global.height)
                             *cos(G.control.T_wave_number*M_PI*x/G.global.length);
		}
	}

        return;
}

void periodic_T_with_random_perturb()
{
        double x,z;
        double random;
        int i, lev;
	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		for(i=0; i<G.nno[lev]; i++)
		{
			random=(1.0*rand())/(1.0*RAND_MAX);
			x = G.X[lev][0][i];
			z = G.X[lev][1][i];
                        G.T[lev][i] = G.control.T_mantle
                             + G.control.T_mantle_perturb
                             *sin(M_PI * (G.global.height - z)/G.global.height)
                             *cos(G.control.T_wave_number*M_PI*x/G.global.length);

			G.T[lev][i] += random * 0.1 * G.control.T_mantle_perturb;
		}
	}

        return;
}

void apply_T_bc_values()
{
	int i,lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
        for(i=0; i<G.nno[lev]; i++)
        {
                if(G.bc.T_bc[lev][i] == 'T')
                {
                        G.T[lev][i] = G.bc.T_bc_val[lev][i];
                }
        }

	return;
}

void interpolate_T_from_file()
{
	int i, nline, lev;
	char filename[STRLEN];
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(filename,"%s/%07d/temperature.%d.%07d", G.old_T_file, G.old_T_timestep, G.ip, G.old_T_timestep);
	fp = fopen(filename, "r");
	if(fp == NULL)
	{
		fprintf(stderr,"Cannot find file: %s\n", filename);
		terminate();
	}

	nline = 0;
	while(fgets(tempstring, STRLEN, fp) != NULL)
		nline++;

	lev = G.old_T_level;

	printf("%d %d\n",nline, G.nno[lev]);

	if(nline != G.nno[lev])
	{
		fprintf(stderr,"The grid of the old T file does not match this model\n change old_T_level may help\n");
		terminate();
	}

	terminate();
	return;
}

void read_T_from_file()
{
	int i, lev;
	char filename[STRLEN];
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(filename,"%s/%07d/temperature.%d.%07d", G.old_T_file, G.old_T_timestep, G.ip, G.old_T_timestep);
	fp = fopen(filename, "r");
	if(fp == NULL)
	{
		fprintf(stderr,"Cannot find file: %s\n", filename);
		terminate();
	}

	lev = G.max_level-1;
	for(i=0; i<G.nno[lev]; i++)
	{
		fgets(tempstring, STRLEN, fp);
		sscanf(tempstring, "%lf", &G.T[lev][i]);
	}

	fclose(fp);
	return;
}
