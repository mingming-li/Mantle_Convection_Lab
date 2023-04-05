#include "global_variables.h"


void output()
{
	char tempstring[STRLEN];

//	if(!folder_exist(G.data_dir))
	if(access(G.data_dir, F_OK) == -1)
		mkdir(G.data_dir,0700);

	sprintf(tempstring,"%s/%07d",G.data_dir,G.timestep);
//	if(!folder_exist(tempstring))
	if(access(tempstring, F_OK) == -1)
		mkdir(tempstring,0700);

	MPI_Barrier(MPI_COMM_WORLD);

	output_coord();
	output_temperature();
	output_stress();
	output_viscosity();
	output_velocity();
	if(G.control.grain_damage==1)
		output_grain_size();

	return;
}


void output_coord()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/coord.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}

	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%e %e\n",G.X[lev][0][i],G.X[lev][1][i]);
	fclose(fp);

	return;
}

void output_temperature()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/temperature.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}

	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%e\n",G.T[lev][i]);
	fclose(fp);
	return;
}

void output_stress()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/stress.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}

	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%e\n",G.stress[lev][i]/1e6);
	fclose(fp);
	return;
}


void output_grain_size()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/grain_size.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}

	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%e\n",log10(G.gs[lev][i]));
	fclose(fp);
	return;
}

void output_viscosity()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/viscosity.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}

	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%e\n",G.vis[lev][i]);
	fclose(fp);
	return;
}

void output_velocity()
{
	int i;
	int lev = G.max_level - 1;
	char tempstring[STRLEN];
	FILE *fp;

	sprintf(tempstring,"%s/%07d/velocity.%d.%07d",G.data_dir,G.timestep,G.ip,G.timestep);
	fp = fopen(tempstring,"w");
	if(fp == NULL)
	{
		fprintf(stderr,"cannot open file %s\n",tempstring);
		terminate();
	}
	for(i=0; i<G.nno[lev]; i++)
		fprintf(fp,"%+e %+e\n",G.V[lev][0][i],G.V[lev][1][i]);
	fclose(fp);
	return;
}

int folder_exist(const char *path)
{
    struct stat stats;

    stat(path, &stats);

    if (S_ISDIR(stats.st_mode))
        return 1;

    return 0;
}
