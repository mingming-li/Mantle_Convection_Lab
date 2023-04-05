#include "global_variables.h"

void insertionSort();
void quickSort();
int partition();

int nextPrime();
bool isPrime();

void sort(arr,n)
	int *arr,n;
{
//	insertionSort(arr,n);
	quickSort(arr,0,n-1);

	return;
}

void insertionSort(arr, n)
        int *arr,n;
{
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

void swap_d(double* a, double* b)
{
	double t = *a;
	*a = *b;
	*b = t;
}

int partition (int arr[], int low, int high)
{
	int pivot = arr[high];
	int i = (low - 1),j;

	for (j = low; j <= high - 1; j++)
	{
		if (arr[j] < pivot)
		{
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
	if (low < high)
	{
		int pi = partition(arr, low, high);
		quickSort(arr, low, pi - 1);
		quickSort(arr, pi + 1, high);
	}
}

int remove_dup(arr, n)
	int *arr,n;
{
	int i,j,*temp;

	if (n==0 || n==1)
		return n;

	temp = malloc(n*sizeof(int));

	j = 0;
	for(i=0; i<n-1; i++)

	if (arr[i] != arr[i+1])
		temp[j++] = arr[i];

	temp[j++] = arr[n-1];

	for(i=0; i<j; i++)
		arr[i] = temp[i];

	free(temp);
	 
	return j;
}

int cmp(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int node_index(key,lev)
        int key,lev;
{
	int i;

	if(key == -1)
		return -1;

	if(G.control.use_hash)
	{
		struct HASH_ITEM *ptr = search(G.hash_node[lev], G.nno_prime[lev], key);
		i = ptr->data;
	}
	else
	{
		int *ptr = bsearch(&key, G.node[lev], G.nno[lev], sizeof(int), cmp);
		if(ptr != NULL)
			i = ((int)ptr - (int)G.node[lev])/sizeof(int);
		else
			i = -1;
	}

        return i;
}

void node_hashtable()
{
	int i,lev;

	for(lev=G.min_level; lev<G.max_level; lev++)
	{
		G.nno_prime[lev] = nextPrime(G.nno[lev]);
		G.hash_node[lev] = realloc(G.hash_node[lev], G.nno_prime[lev] * sizeof(struct HASH_ITEM *));
		for(i=0; i<G.nno_prime[lev]; i++)
			G.hash_node[lev][i] = NULL;

		for(i=0; i<G.nno[lev]; i++)
		{
			insert_hash_item(G.hash_node[lev], G.node[lev][i], i, G.nno_prime[lev]);
		}
	}

	/* test
	for(lev=G.min_level; lev<G.max_level; lev++)
	for(i=0; i<G.nno[lev]; i++)
	{
		printf("%d %d %d\n",G.node[lev][i], i, node_index(G.node[lev][i], lev));
	}
	terminate();
	*/

	return;
}

void insert_hash_item(arr, key, data, size)
	int key,data,size;
	struct HASH_ITEM **arr;
{
	struct HASH_ITEM *item = malloc(sizeof(struct HASH_ITEM *));
	int index;

	item->data = data;
	item->key = key;
	index = hash_code(key, size);
	while(arr[index] != NULL)
	{
		index++;
		index %= size;
	}
	arr[index] = item;
	return;
}

int hash_code(key,size)
	int key,size;
{
	return key % size;
}

struct HASH_ITEM *search(arr,size,key)
	struct HASH_ITEM **arr;
	int size,key;
{
	int index = hash_code(key,size);

	while(arr[index] != NULL)
	{
		if(arr[index]->key == key)
			return arr[index];

		index++;
		index %= size;
	}
	return NULL;
}



int nextPrime(int N)
{

    // Base case
    if (N <= 1)
        return 2;

    int prime = N;
    bool found = false;

    while (!found) {
        prime++;

        if (isPrime(prime))
            found = true;
    }

    return prime;
}

bool isPrime(int n)
{
	int i;
    // Corner cases
    if (n <= 1)  return false;
    if (n <= 3)  return true;

    // This is checked so that we can skip
    // middle five numbers in below loop
    if (n%2 == 0 || n%3 == 0) return false;

    for (i=5; i*i<=n; i=i+6)
        if (n%i == 0 || n%(i+2) == 0)
           return false;

    return true;
}

int max(a,b)
        int a,b;
{
        return a>b?a:b;
};

int min(a,b)
        int a,b;
{
        return a<b?a:b;
};


int log4(a)
	int a;
{
	int b=0;
	while (a >>= 2) b++;
	return b;
}

int value_in_array(a, arr, n)
	int a, *arr, n;
{
	int i;
	for(i=0; i<n; i++)
		if(a == arr[i])
			return i;
	return -1;
}

int value_in_sorted_array(a, arr, n)
	int a, *arr, n;
{
	int i;

	int *ptr = bsearch(&a, arr, n, sizeof(int), cmp);
	if(ptr != NULL)
		i = ((int)ptr - (int)arr)/sizeof(int);
	else
		i = -1;

	return i;
}

double global_v_prod(x, y, lev)
	double *x, *y;
	int lev;
{
	int i, j, k, neq_m, neq, *m;
	int *nsh, **seq, *nseq;
	double global_prod, prod, temp;

	neq = G.neq[lev];
	m = G.smpi[lev].veq_nsh;
	prod = 0.0;
	for(i=0; i<neq; i++)
		prod += x[i] * y[i] / m[i];

        MPI_Allreduce(&prod, &global_prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_prod;
}

double global_T_prod(x, y, lev)
	double *x, *y;
	int lev;
{
	int i, j, k, neq_m, neq, *m;
	int *nsh, **seq, *nseq;
	double global_prod, prod, temp;

	neq = G.heat[lev].neq;
	m = G.smpi[lev].Teq_nsh;

	prod = 0.0;
	for(i=0; i<neq; i++)
		prod += x[i] * y[i] / m[i];

        MPI_Allreduce(&prod, &global_prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_prod;
}

double global_p_prod(x, y, lev)
	double *x, *y;
	int lev;
{
	int i;
	double prod, global_prod;

	prod = 0.0;

        for(i=0; i<G.nel[lev]; i++)
        	prod += x[i] * y[i];

        MPI_Allreduce(&prod, &global_prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_prod;
}
