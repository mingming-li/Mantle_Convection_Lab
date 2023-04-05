#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>
#include <stdbool.h>
#include "function_def.h"
#include <sys/stat.h>
#include <assert.h>
#include <unistd.h>


#define VPTS 4
#define CH_BIT 8
#define NSD 2
#define NEN 4
#define NEE 8
#define NED 2
#define NEG 4
#define NDF 2
#define LEVEL_BITS 4
#define LEAF_BITS 1
#define MAX_SHARE_IP 20
#define MAX_LEVEL 20
#define STRLEN 1000
#define MAX_T_NCOL 18
#define MAX_V_NCOL 36
#define MIN(a,b) ((a)<(b))?(a):(b)
#define MAX(a,b) ((a)>(b))?(a):(b)

double start_time,end_time;

struct TREE_NODE
{
	int loc;
        struct TREE_NODE *child[VPTS];
};

struct LIST_NODE
{
	int data;
	struct LIST_NODE *next;
};

struct HEAT
{
	int *Kncol;
	int **Krow_col;
	double **Krow_val;

	int *Cncol;
	int **Crow_col;
	double **Crow_val;

        double *F;
        double *d;
        double *ddot;
        double dt;
        double theta;
        double *tau[NSD];
	int *id;
	int *eq_node;
	int **lm;
        int neq, neq_m, neq_c, neq_lr;
};

struct P2N
{
	int *nel;
	int *el[NEN];
};

struct ECO
{
	double *area;
};



struct NODE_SYSTEM
{
        double nodal[NEN];
        double gauss[NEN];
        double center;
};

struct HASH_ITEM
{
	int data;
	int key;
};

struct SEND
{
	int nint;
	int *isend;
	int ndouble;
	double *dsend;
};


struct LIST
{
	int key;
	struct LIST *p;
};

struct GLOBAL_MESH
{
	int *nel;
	int **ip_nel;
	int **old_ip_nel;
	int **min_Mcode;
	int **max_Mcode;
	double height;
	double length;
};

struct INTERPOLATE
{
	int *node_amount;
	int *eq_amount;
	int *local_eq[VPTS];
	int *local_eq_n;
	int **remote_node;
	int **remote_node_n;
	int **remote_eq;

	int nip_send, *ip_send;
	int nip_recv, *ip_recv;


	int *local_l;
	int *local_u;
	int nlocal;

	int *send_l[MAX_SHARE_IP];
	int *send_u[MAX_SHARE_IP];
	int nsend[MAX_SHARE_IP];

	int *recv_l[MAX_SHARE_IP];
	int *recv_u[MAX_SHARE_IP];
	int nrecv[MAX_SHARE_IP];



	//the following are for interp, 
	//because of multiple to 1 node interp, 
	//we need another dimension

	int *local_nlv;
	int *local_lv[VPTS];
	int *send_nlv[MAX_SHARE_IP];
	int *send_lv[MAX_SHARE_IP][VPTS];

	int remote_project;
	int remote_interp;
};

struct MPI_SHARE
{
	int nshareip;
        int ip[MAX_SHARE_IP];

	int n_veq[MAX_SHARE_IP];
	int n_veq_m[MAX_SHARE_IP];
	int *veq[MAX_SHARE_IP];

	int *veq_nsh;

	int n_Teq[MAX_SHARE_IP];
	int n_Teq_m[MAX_SHARE_IP];
	int *Teq[MAX_SHARE_IP];


	int *Teq_nsh;

	int nno[MAX_SHARE_IP];
	int *node[MAX_SHARE_IP];
	int *nel[MAX_SHARE_IP];
	int **el[MAX_SHARE_IP];

	int *share_node;
	int nshare_node;

	int send_number[MAX_SHARE_IP];
	int *send_row[MAX_SHARE_IP];
	int *send_matr_row[MAX_SHARE_IP];
	int *send_matr_col[MAX_SHARE_IP];

	int recv_number[MAX_SHARE_IP];
	int *recv_row[MAX_SHARE_IP];
	int *recv_matr_ip[MAX_SHARE_IP];
	int *recv_matr_col[MAX_SHARE_IP];


	int *Tshare_node;
	int Tnshare_node;

	int Tsend_number[MAX_SHARE_IP];
	int *Tsend_row[MAX_SHARE_IP];
	int *Tsend_matr_row[MAX_SHARE_IP];
	int *Tsend_matr_col[MAX_SHARE_IP];

	int Trecv_number[MAX_SHARE_IP];
	int *Trecv_row[MAX_SHARE_IP];
	int *Trecv_matr_ip[MAX_SHARE_IP];
	int *Trecv_matr_col[MAX_SHARE_IP];
};

struct MAP_1
{
	int n;
	int *eq_l;
	int *eq_u;
};

struct HANGLE_NODES
{
	int ip;
	int node;
	int n[NSD];
};

struct HANGLE
{
	int nh;
	int *node;
	int *ni;
	int *cnode[NSD];
	int *ci[NSD];
};

struct HANGLE_EQ
{
	int node;
	int heq;
	int ceq[NSD];
	int cip[NSD];
};

struct BC
{
        char top_T_bc;
        char bot_T_bc;
        char lef_T_bc;
        char rig_T_bc;

        double top_T_bc_val;
        double bot_T_bc_val;
        double lef_T_bc_val;
        double rig_T_bc_val;

	char **T_bc;
	double **T_bc_val;

        char top_v_bc[NSD];
        char bot_v_bc[NSD];
        char lef_v_bc[NSD];
        char rig_v_bc[NSD];

        double top_v_bc_val[NSD];
        double bot_v_bc_val[NSD];
        double lef_v_bc_val[NSD];
        double rig_v_bc_val[NSD];

	char ***v_bc;
	double ***v_bc_val;

};

struct VISCOSITY
{
	double A;
	double refT;
	double lm_jump;
	double highest;
	double lowest;
};

struct CONTROL
{
	int T_ini;
	double T_mantle;
	double T_mantle_perturb;
	int T_wave_number;
	int use_hash;
	int use_dense_matrix;
	int mg_coarse_iteration;
	int mg_finest_iteration;
	int mg_cycle;

	int mpi_KF_method;

        double Q;
	double theta;
	int AMR;
	int cg_precondition;
	int Tcg_precondition;
	double accuracy;
	int auto_accuracy;
	double citcom_accu;

	int dev_mode;
	int check_v_only;
	int max_citcom_cycle;
	int max_cg_cycle;

	int grain_damage;
	int show_convg;
	int zero_P;
};

struct All_global_variables
{
	char data_dir[256];
	int max_level;
	int min_level;
	int *max_elx,*max_elz,*max_nel;
	int *max_nox,*max_noz,*max_nno;
	int gnox, gnoz;

	int *nel,*nno,*nno_prime;

	int ip,npx,npz,np,ipx,ipz;
	int timestep,max_timestep,solution_cycle,save_spacing;
	double model_time;

	int ***ien;

	int **Mcode;
	int **node;
	int **Mcode_ip;
	char **ntype;

	struct TREE_NODE *root;
	struct GLOBAL_MESH global;
	struct MPI_SHARE *smpi;

	double length,height;
	double *dx,*dz;
	double ***X;
	double **T;
	double **buoyancy;
	double **vis;
	double Ra;
	double ***V;
	double *P;
	double **stress;
	double **gs;

	char dimension[3];
	char geometry[6];
	char solver[3];

	struct BC bc;

	int ***id;
	int ****lm;
	int **veq_node;
	int **veq_dof;


	int *neq;
	int *neq_f;
	int *neq_m;
	int *neq_c;
	int *neq_lr;

	double points[NEN][NSD];
	double gpoints[NEN][NSD];

        struct NODE_SYSTEM N[NEN];
        struct NODE_SYSTEM dN_dxi[NEN];
        struct NODE_SYSTEM dN_deta[NEN];

	double **dx_dxi, dx_deta, dz_dxi, **dz_deta;
	double **jacobi;

	struct CONTROL control;
	struct VISCOSITY viscosity;

	FILE *fpout;

	struct HASH_ITEM ***hash_node;

	int **Kncol;
	int ***Krow_col;
	double ***Krow_val;
	double **Kdiag;

	int **Gncol;
	int ***Grow_col;
	double ***Grow_val;

	double **F;
	double **d;
	double **p;
	double **err;


	struct HANGLE *hangle;
	struct HANGLE *lg_hangle;
	struct INTERPOLATE *interp;
	struct INTERPOLATE *project;

	struct HEAT *heat;

	struct LIST_NODE *head;
	struct LIST_NODE *addition;

	double CPU_time;

	char old_T_file[STRLEN];
	int old_T_level, old_T_timestep;

	int lg_mul;

	struct P2N *p2n;
	struct ECO *eco;
}G;
