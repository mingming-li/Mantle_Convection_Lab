#include "global_variables.h"

void append();
void printlist();
int list_to_array();
void deletelist();

void deletelist(struct LIST_NODE** head_ref)
{
	/* deref head_ref to get the real head */
	struct LIST_NODE* current = *head_ref;
	struct LIST_NODE* next;

	while (current != NULL)
	{
	next = current->next;
	free(current);
	current = next;
	}

	/* deref head_ref to affect the real head back
	in the caller. */
	*head_ref = NULL;
}

int list_to_array(list, arr)
	struct LIST_NODE *list;
	int *arr;
{
	int i;
	i=0;
	while (list != NULL)
	{
		arr[i++] = list->data;
		list = list->next;
	}
	return i;
}


void append(struct LIST_NODE **head_ref, int new_data)
{
        struct LIST_NODE *new_node = malloc(sizeof(struct LIST_NODE *));

        struct LIST_NODE *last = *head_ref;

        new_node->data = new_data;
        new_node->next = NULL;

	//from GeeksForGeeks
        /* If the Linked List is empty, then make the new node as head */
        if (*head_ref == NULL)
        {
        	*head_ref = new_node;
        	return;
        }

        /* Else traverse till the last node */
        while (last->next != NULL)
                last = last->next;

        /* Change the next of last node */
        last->next = new_node;
        return;
}

void build_quad_tree()
{
	int lev;
	static int been_here=0;

	if(been_here == 0)//initial
	{
		for(lev=G.min_level; lev<G.max_level; lev++)
			G.nel[lev] = 0;
		initial_refine(G.max_level - 1);//initial setup, can be very flexible, user defined
		been_here = 1;
	}
	else
	{
		adaptive_mesh_refine(G.max_level - 1);
	}

	//project Mcode to lower levels for multigrid
	for(lev=G.max_level-1; lev>G.min_level; lev--)
	{
		G.nel[lev-1] = 0;
		G.Mcode[lev-1] = realloc(G.Mcode[lev-1], G.nel[lev] * sizeof(int));
		project_Mcode(lev,lev - 1);
		G.Mcode[lev-1] = realloc(G.Mcode[lev-1], G.nel[lev-1] * sizeof(int));
	}
	return;
}

void adaptive_mesh_refine(lev)
	int lev;
{
	return;
}

void initial_refine(lev)
	int lev;
{
	void refine_bottom();
	void refine_top();
	G.Mcode[lev] = realloc(G.Mcode[lev], G.max_nel[lev] * 2 * sizeof(int));//* 2 here is to secure enough memory

	G.Mcode[lev][0] = G.ip * G.max_nel[G.max_level-1]; // using global index for variable ip

	if(G.control.dev_mode == 0)
	{
		split_mul_lev(G.Mcode[lev][0], lev, G.max_level - 1, 0, 0);//even mesh
	}
	else if(G.control.dev_mode == 1)
	{
		if(G.ip % 2 == 0)
			split_mul_lev(G.Mcode[lev][0], lev, 1, 0, 0);
		else if(G.ip == 1)
			split_mul_lev(G.Mcode[lev][0], lev, 2, 0, 0);
		else
		{
			split_mul_lev(G.Mcode[lev][0], lev, 2, 0, 0);
			split_mul_lev(212, lev, 1, 2, 0);
			split_mul_lev(220, lev, 1, 2, 0);
			split_mul_lev(244, lev, 1, 2, 0);
			split_mul_lev(252, lev, 1, 2, 0);
		}
	}
	else if(G.control.dev_mode == 2)
	{
		split_mul_lev(G.Mcode[lev][0], lev, 5, 0, 0);//even mesh
		if(G.ip % 4 == 0)
			refine_bottom(lev);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else if(G.control.dev_mode == 3)
	{
		split_mul_lev(G.Mcode[lev][0], lev, 4, 0, 0);//even mesh
		if(G.ip % 8 == 0)
			refine_bottom(lev);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else if(G.control.dev_mode == 4)
	{
		if(G.ip % 2 == 0)
			split_mul_lev(G.Mcode[lev][0], lev, 2, 0, 0);//even mesh
		else
			split_mul_lev(G.Mcode[lev][0], lev, 1, 0, 0);//even mesh
	}
	else if(G.control.dev_mode == 5)
	{
		if(G.ip % 2 == 0)
			split_mul_lev(G.Mcode[lev][0], lev, 6, 0, 0);//even mesh
		else
			split_mul_lev(G.Mcode[lev][0], lev, 5, 0, 0);//even mesh
	}
	else if(G.control.dev_mode == 6)
	{
		int i;
		i = G.ip % G.npz;
		i = G.max_level - 1 - i;
		i = MAX(i, 4);
		split_mul_lev(G.Mcode[lev][0], lev, i, 0, 0);//even mesh
	}
	else if(G.control.dev_mode == 7)
	{
		if(G.ip % G.npz == 0)
			split_mul_lev(G.Mcode[lev][0], lev, G.max_level - 1, 0, 0);//even mesh
		else
			split_mul_lev(G.Mcode[lev][0], lev, G.max_level - 2, 0, 0);//even mesh
	}
	else if(G.control.dev_mode == 8)
	{
		split_mul_lev(G.Mcode[lev][0], lev, 2, 0, 0);
		split_mul_lev(64, lev, 1, 2, 0);
		split_mul_lev(80, lev, 1, 2, 0);
		split_mul_lev(96, lev, 1, 2, 0);
		split_mul_lev(112, lev, 1, 2, 0);
		split_mul_lev(192, lev, 1, 2, 0);
		split_mul_lev(208, lev, 1, 2, 0);
		split_mul_lev(224, lev, 1, 2, 0);
		split_mul_lev(240, lev, 1, 2, 0);
		split_mul_lev(212, lev, 1, 3, 0);
		split_mul_lev(220, lev, 1, 3, 0);
		split_mul_lev(244, lev, 1, 3, 0);
		split_mul_lev(252, lev, 1, 3, 0);
	}
	else if(G.control.dev_mode == 9)
	{
		split_mul_lev(G.Mcode[lev][0], lev, 5, 0, 0);//even mesh
		if(G.ip % 4 == 0)
			refine_bottom(lev);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else if(G.control.dev_mode == 10)
	{
		split_mul_lev(G.Mcode[lev][0], lev, 4, 0, 0);//even mesh
                if(G.ip % 8 == 0)
                        refine_bottom(lev);
		else if(G.ip % 8 == 7)
			refine_top(lev);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else if(G.control.dev_mode == 11)
	{
		int i, a, b, mul,width;
		split_mul_lev(G.Mcode[lev][0], lev, G.max_level-2, 0, 0);//even mesh

		mul = (int)pow(2, G.max_level-1);
		width = mul/8;

		for(i=0; i<mul*mul; i+=4)
		{
			demorton_no_lev(i, &a, &b);
			if((b<width || b >=mul-width) || (a<width || a >=mul-width))
				split_mul_lev(i, lev, 1, G.max_level-2, 0);//even mesh
		}

	}
	else if(G.control.dev_mode == 12)
	{
		int i, a, b, mul,width;
		split_mul_lev(G.Mcode[lev][0], lev, G.max_level-2, 0, 0);//even mesh

		mul = (int)pow(2, G.max_level-1);
		width = mul/8;

		for(i=0; i<mul*mul; i+=4)
		{
			demorton_no_lev(i, &a, &b);
			if(b >=mul-width)
				split_mul_lev(i, lev, 1, G.max_level-2, 0);//even mesh
		}

	}
	else if(G.control.dev_mode == 13)
	{
		int i, a, b, mul,width;
		split_mul_lev(G.Mcode[lev][0], lev, G.max_level-6, 0, 0);//even mesh
		mul = (int)pow(2, G.max_level-1);
		width = mul / 64;

		for(i=0; i<mul*mul; i+=1024)
		{
			demorton_no_lev(i, &a, &b);
			if(b <= width)
				split_mul_lev(i, lev, 5, G.max_level-6, 0);
			else if(b > 1*width && b<= 4*width)
				split_mul_lev(i, lev, 4, G.max_level-6, 0);
			else if(b > 4*width && b<= 8*width)
				split_mul_lev(i, lev, 3, G.max_level-6, 0);
			else if(b > 8*width && b<= 12*width)
				split_mul_lev(i, lev, 2, G.max_level-6, 0);
			else if(b > 12*width)
				split_mul_lev(i, lev, 1, G.max_level-6, 0);
		}
		if(G.ip == 0)
			fprintf(stderr,"done list\n");
	}
	else if(G.control.dev_mode == 14)
	{
		int i, a, b, mul,width;
		split_mul_lev(G.Mcode[lev][0], lev, G.max_level-2, 0, 0);//even mesh

		mul = (int)pow(2, G.max_level-1);
		width = mul/2+1;

		for(i=0; i<mul*mul; i+=4)
		{
			demorton_no_lev(i, &a, &b);
			if(b <= width)
				split_mul_lev(i, lev, 1, G.max_level-2, 0);//even mesh
		}

	}


	G.nel[lev] = list_to_array(G.head, G.Mcode[lev]);//from the next to head to escape the first 0
	deletelist(&G.head);//I do not need the list any more

	G.Mcode[lev] = realloc(G.Mcode[lev], G.nel[lev] * sizeof(int));
	if(G.ip == 0)
		printf("sort list\n");
	sort(G.Mcode[lev], G.nel[lev]);
	if(G.ip == 0)
		printf("remove list dup\n");
	G.nel[lev] = remove_dup(G.Mcode[lev], G.nel[lev]);

	/*
	for(int i=0;i<G.nel[lev]; i++)
		printf("%d %d\n",i,G.Mcode[lev][i]);
	terminate();
	*/

	append_level_to_Mcode(lev);

	return;
}



void split_mul_lev(node, lev, num, level, log)
	int node, lev, num, level, log;
{
	int i, child[VPTS];

	if(num == 0)
		return;

	split_Mcode(node, child, level);

	if(log == 0)
	{
		for(i=0; i<VPTS; i++)
			append(&G.head, child[i]);
	}
	else
	{
		for(i=0; i<VPTS; i++)
			append(&G.addition, child[i]);
	}

	for(i=0; i<VPTS; i++)
		split_mul_lev(child[i], lev, num - 1, level+1, log);//split child if num - 1 > 0
}

void split_Mcode(parent, child, level)
	int parent, level, child[VPTS];
{
	int mul,i;

	level += 1;
	for(i=0; i<VPTS; i++)
	{
		mul = G.max_nel[G.max_level - 1] >> 2 * level;
		child[i] = parent + i * mul;
	}
}

void append_level_to_Mcode(lev)
	int lev;
{
	int e, level;

	for(e=0; e<G.nel[lev]-1; e++)
	{
		level = log4(G.max_nel[lev] / (G.Mcode[lev][e+1] - G.Mcode[lev][e])); // level can be calculated from the diff of neiboring Mcodes


		G.Mcode[lev][e] = append_level(G.Mcode[lev][e], level);
	}

	e = G.nel[lev] - 1;

	level = log4(G.max_nel[lev] / (G.Mcode[lev][e] - Mcode_loc(G.Mcode[lev][e-1])));//last element, so use [e] - [e-1]
	G.Mcode[lev][e] = append_level(G.Mcode[lev][e], level);
	return;
}


void project_Mcode(ulev,llev)
	int ulev, llev;
{
	int e, level, loc;

	e = 0;
	while(e < G.nel[ulev])
	{
		level = Mcode_lev(G.Mcode[ulev][e]);
		loc = Mcode_loc(G.Mcode[ulev][e]);

		//The level below determines whether to skip 1 or 4
		if(level < ulev)
		{
			e++;
		}
		else
		{
			e += 4;
			level -= 1;
		}
		
		//no matter what, always copy the first Mcode. 
		G.Mcode[llev][G.nel[llev]] = append_level(loc, level);
		G.nel[llev]++;
	}

	return;
}

int Mcode_loc(z)
	int z;
{
	return z >> LEVEL_BITS;
}

int Mcode_lev(z)
	int z;
{
	return (z & ((1U << LEVEL_BITS)-1));
}

int append_level(code, lev)
        int lev,code;
{
        code = (code << LEVEL_BITS) | lev;
        return code;
}

int get_Mcode_level(z)
	int z;
{
	printf("not working here\n");
	return (z >> LEAF_BITS) & ((1U<<LEVEL_BITS) -1);
}

int get_Mcode_location(z)
	int z;
{
	printf("not working here\n");
	return z >> (LEAF_BITS + LEVEL_BITS);
}

int get_Mcode_leaf(z)
	int z;
{
	printf("not working here\n");
	return (z & ((1U << LEAF_BITS)-1));
}



void demorton_no_lev(z,x,y)
        int z,*x,*y;
{
        int i;

        *x=*y=0;
        for(i=0; i<sizeof(int)*CH_BIT; i+=2)
        {
                *x |= (z & 1U << (i+1)) >> (i/2+1);
                *y |= (z & 1U << i) >> i/2;
        }
        return;
}

int morton(x,z)
        int x,z;
{
        int i, loc=0;

        for(i=0; i<sizeof(int) * CH_BIT; i++)
                loc |= (x & 1U << i) << (i+1) | (z & 1U << i) << i;

        return loc;
}


void refine_bottom(lev)
	int lev;
{
	int start_Mcode, a, b;
	struct LIST_NODE *current = G.head;

	start_Mcode = G.ip * G.max_nel[G.max_level-1];

	if(G.control.dev_mode == 2)
	while (current != NULL)
	{
		demorton_no_lev(current->data-start_Mcode, &a, &b);
		if(b <= 112)
			split_mul_lev(current->data, lev, 4, 5, 1);
		else if(b > 112 && b <= 240)
			split_mul_lev(current->data, lev, 3, 5, 1);
		else if(b > 240 && b <= 368)
			split_mul_lev(current->data, lev, 2, 5, 1);
		else if(b > 368)
			split_mul_lev(current->data, lev, 1, 5, 1);
		current = current->next;
	}
	else if(G.control.dev_mode == 3)
	while (current != NULL)
	{
		demorton_no_lev(current->data-start_Mcode, &a, &b);
		if(b <= 48)
			split_mul_lev(current->data, lev, 4, 4, 1);
		else if(b > 48 && b <= 112)
			split_mul_lev(current->data, lev, 3, 4, 1);
		else if(b > 112 && b <= 240)
			split_mul_lev(current->data, lev, 2, 4, 1);
		else if(b > 240)
			split_mul_lev(current->data, lev, 1, 4, 1);
		current = current->next;
	}
	else if(G.control.dev_mode == 9)
	while (current != NULL)
	{
		demorton_no_lev(current->data-start_Mcode, &a, &b);
		if(b <= 48)
			split_mul_lev(current->data, lev, 4, 5, 1);
		else if(b > 48 && b <= 240)
			split_mul_lev(current->data, lev, 3, 5, 1);
		else if(b > 240 && b <= 368)
			split_mul_lev(current->data, lev, 2, 5, 1);
		else if(b > 368)
			split_mul_lev(current->data, lev, 1, 5, 1);
		current = current->next;
	}
	else if(G.control.dev_mode == 10)
	while (current != NULL)
	{
		demorton_no_lev(current->data-start_Mcode, &a, &b);
		if(b <= 48)
			split_mul_lev(current->data, lev, 4, 4, 1);
		else if(b > 48 && b <= 112)
			split_mul_lev(current->data, lev, 3, 4, 1);
		else if(b > 112 && b <= 240)
			split_mul_lev(current->data, lev, 2, 4, 1);
		else if(b > 240)
			split_mul_lev(current->data, lev, 1, 4, 1);
		current = current->next;
	}

	deletelist(&G.head);//I do not need the list any more
	G.head = G.addition;

	return;
}

void refine_top(lev)
	int lev;
{
	int start_Mcode, a, b;
	struct LIST_NODE *current = G.head;

	start_Mcode = G.ip * G.max_nel[G.max_level-1];

	while (current != NULL)
	{
		demorton_no_lev(current->data-start_Mcode, &a, &b);
		if(b < 48)
			split_mul_lev(current->data, lev, 1, 4, 1);
		else if(b >= 48 && b < 112)
			split_mul_lev(current->data, lev, 2, 4, 1);
		else if(b >= 112 && b < 240)
			split_mul_lev(current->data, lev, 3, 4, 1);
		else if(b >= 240)
			split_mul_lev(current->data, lev, 4, 4, 1);
		current = current->next;
	}



	deletelist(&G.head);//I do not need the list any more
	G.head = G.addition;

	return;
}
