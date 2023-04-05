#include "global_variables.h"

void construct_shape_function()
{
	static int been_here = 0;
	if(been_here == 0)
	{
		shape_function();
		dshape_function_dxi_deta();
		been_here = 1;
	}
	dx_dxi_dz_deta();
	jacobi_in_integration();
	return;
}

void shape_function()
{
        int a,d;
        double xi,eta;
        double BETA;

        BETA=1.0/sqrt(3.0);

        G.points[0][0]=-1.0;
        G.points[0][1]=-1.0;

        G.points[1][0]=-1.0;
        G.points[1][1]=1.0;

        G.points[2][0]=1.0;
        G.points[2][1]=1.0;

        G.points[3][0]=1.0;
        G.points[3][1]=-1.0;

        for(a=0;a<NEN;a++)
        {
                for(d=0;d<NSD;d++)
                {
                        G.gpoints[a][d]=G.points[a][d]*BETA;
                }
        }

        for(d=0;d<NEN;d++)
        {
                xi = G.points[d][0];
                eta= G.points[d][1];

                G.N[0].nodal[d]=0.25*(1-xi)*(1-eta);
                G.N[1].nodal[d]=0.25*(1-xi)*(1+eta);
                G.N[2].nodal[d]=0.25*(1+xi)*(1+eta);
                G.N[3].nodal[d]=0.25*(1+xi)*(1-eta);

                xi = G.gpoints[d][0];
                eta= G.gpoints[d][1];

                G.N[0].gauss[d]=0.25*(1-xi)*(1-eta);
                G.N[1].gauss[d]=0.25*(1-xi)*(1+eta);
                G.N[2].gauss[d]=0.25*(1+xi)*(1+eta);
                G.N[3].gauss[d]=0.25*(1+xi)*(1-eta);
        }
        G.N[0].center=0.25;
        G.N[1].center=0.25;
        G.N[2].center=0.25;
        G.N[3].center=0.25;


	return;
}

void dshape_function_dxi_deta()
{
        int d;
        double xi,eta;

        for(d=0;d<NEN;d++)
        {
                xi = G.points[d][0];
                eta= G.points[d][1];

                G.dN_dxi[0].nodal[d]=0.25*(-1.0)*(1.0-eta);
                G.dN_dxi[1].nodal[d]=0.25*(-1.0)*(1.0+eta);
                G.dN_dxi[2].nodal[d]=0.25*( 1.0)*(1.0+eta);
                G.dN_dxi[3].nodal[d]=0.25*( 1.0)*(1.0-eta);

                G.dN_deta[0].nodal[d]=0.25*(1.0-xi)*(-1.0);
                G.dN_deta[1].nodal[d]=0.25*(1.0-xi)*(+1.0);
                G.dN_deta[2].nodal[d]=0.25*(1.0+xi)*(+1.0);
                G.dN_deta[3].nodal[d]=0.25*(1.0+xi)*(-1.0);

                xi = G.gpoints[d][0];
                eta= G.gpoints[d][1];

                G.dN_dxi[0].gauss[d]=0.25*(-1.0)*(1-eta);
                G.dN_dxi[1].gauss[d]=0.25*(-1.0)*(1+eta);
                G.dN_dxi[2].gauss[d]=0.25*( 1.0)*(1+eta);
                G.dN_dxi[3].gauss[d]=0.25*( 1.0)*(1-eta);

                G.dN_deta[0].gauss[d]=0.25*(1.0-xi)*(-1.0);
                G.dN_deta[1].gauss[d]=0.25*(1.0-xi)*( 1.0);
                G.dN_deta[2].gauss[d]=0.25*(1.0+xi)*( 1.0);
                G.dN_deta[3].gauss[d]=0.25*(1.0+xi)*(-1.0);


        }
        xi = 0.0;
        eta= 0.0;

        G.dN_dxi[0].center=0.25*(-1.0)*(1-eta);
        G.dN_dxi[1].center=0.25*(-1.0)*(1+eta);
        G.dN_dxi[2].center=0.25*( 1.0)*(1+eta);
        G.dN_dxi[3].center=0.25*( 1.0)*(1-eta);

        G.dN_deta[0].center=0.25*(1.0-xi)*(-1.0);
        G.dN_deta[1].center=0.25*(1.0-xi)*( 1.0);
        G.dN_deta[2].center=0.25*(1.0+xi)*( 1.0);
        G.dN_deta[3].center=0.25*(1.0+xi)*(-1.0);

	return;
}

void dx_dxi_dz_deta()
{
        int e,lev,level;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.dx_dxi[lev] = realloc(G.dx_dxi[lev], G.nel[lev] * sizeof(double));
		G.dz_deta[lev] = realloc(G.dz_deta[lev], G.nel[lev] * sizeof(double));
	}

	for(lev=G.min_level; lev<G.max_level; lev++)
        for(e=0;e<G.nel[lev];e++)
        {
		level = Mcode_lev(G.Mcode[lev][e]);
                G.dx_dxi[lev][e]= 0.5*G.dx[level];
                G.dz_deta[lev][e]=0.5*G.dz[level];
        }

        G.dx_deta=0.0;
        G.dz_dxi=0.0;

	return;
}

void jacobi_in_integration()
{
	double temp[NSD][NSD];
	int lev,e;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.jacobi[lev] = realloc(G.jacobi[lev], G.nel[lev] * sizeof(double));
		for(e=0;e<G.nel[lev];e++)
		{
			temp[0][0]=G.dx_dxi[lev][e];
			temp[0][1]=G.dx_deta;
			temp[1][0]=G.dz_dxi;
			temp[1][1]=G.dz_deta[lev][e];
			G.jacobi[lev][e]=determinant(temp);
		}

	}

	return;
}

double determinant(a)
	double a[NSD][NSD];
{
	double value;
	value=a[0][0]*a[1][1]-a[0][1]*a[1][0];
	return value;
}

double Na(a,xi,eta)
	int a;
	double xi,eta;
{
	switch (a)
	{
		case 0:
			return 0.25*(1.0 - xi) * (1.0 - eta);
		case 1:
			return 0.25*(1.0 - xi) * (1.0 + eta);
		case 2:
			return 0.25*(1.0 + xi) * (1.0 + eta);
		case 3:
			return 0.25*(1.0 + xi) * (1.0 - eta);
		default:
			fprintf(stderr,"Na error: value of a should be 0-3\n");
			terminate();
	}

	return 0.0;
}
