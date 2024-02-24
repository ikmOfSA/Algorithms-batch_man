#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gen_co_mat(int n, double co_mat[][n], double b_vec[], double alpha, double p, double K, double c, int m);
int solve_linear_eqs(int n, double co_mat[][n], double b_vec[], double J_vec[]);
int gauss_elim(int n, double co_mat[][n], double b_vec[]);
int back_elim(int n, double co_mat[][n], double b_vec[], double J_vec[]);
int part_pivot(int n, double co_mat[][n], double b_vec[], int j);
int row_swap(int n, double co_mat[][n], double b_vec[], int i, int j);
int print_mat(int n, double co_mat[][n], double b_vec[]);
int improv(int n, double J_vec[], double alpha, double p, double K, double c, int m);
int forward(int n, double J_vec[], double J_vec_new[], double alpha, double p, double K, double c);

int main()
{	
	int n = 10;		//number of states, starts from 0;
	double alpha = 0.95;
	double p = 0.5;
	double K = 20;
	double c = 1;

	double *J_vec = calloc(n, sizeof(double));	//current cost vector;
	double *J_vec_new = calloc(n, sizeof(double));	//next cost vector, i.e., the improvement step
	
	for (int j = 0; j < 1000; j++)
		forward(n, J_vec, J_vec_new, alpha, p, K, c);	//improvement step

	int m = forward(n, J_vec, J_vec_new, alpha, p, K, c);	//run again just for getting m, which is the policy threshold;

	for (int i = 0; i < n; i++)		//print cost and policy
		printf("%5.2f\t", J_vec[i]);
	printf("\nm = %d\n", m);

	return 0;    
}

/*the improvement step, i.e., get next cost vector;*/
int forward(int n, double J_vec[], double J_vec_new[], double alpha, double p, double K, double c)
{
	int m = n - 1;		//policy threshold, at least be n - 1 according to the problem;
	double proc, unproc;

	for (int i = 0; i < n - 1; i++) {		// i is the state;

		proc = K + alpha * (1 - p) * J_vec[0] + alpha * p * J_vec[1];		// cost if process the batch;
		unproc = c * i + alpha * (1 - p) * J_vec[i] + alpha * p * J_vec[i + 1];	// if not;

		if (proc < unproc) {
			J_vec_new[i] = proc;
			m = (m > i) ? i : m;	
		}
		else 
			J_vec_new[i] = unproc;	
	}
	
	J_vec_new[n-1] = K + alpha * (1-p) * J_vec[0] + alpha * p * J_vec[1];	//must process when reach n-1;

	for (int j = 0; j < n; j++)
		J_vec[j] = J_vec_new[j];
	
	return m;
}

int improv(int n, double J_vec[], double alpha, double p, double K, double c, int m)
{
	int i;
	int new_m = 0;	
	double proc, unproc;

	for(i = 1; i < n - 1; i++) {
		proc = K + alpha * (1 - p) * J_vec[0] + alpha * p * J_vec[1];
		unproc = c * i + alpha * (1 - p) * J_vec[i] + alpha * p * J_vec[i + 1];
		if (proc < unproc) {
			new_m = i;
			break;	
		}
	}	

	if (new_m == 0)
		new_m = n - 1;
	
	return new_m;
}

int gen_co_mat(int n, double co_mat[][n], double b_vec[], double alpha, double p, double K, double c, int m)
{
	int i, j;

	for(i = 0; i < n; i++) {

		for (j = 0; j < n; j++) {
			co_mat[i][j] = 0;	
		}

		b_vec[i] = 0;
	}
	
	for (i = 0; i < n; i++) {
		if (i < m) {
			co_mat[i][i] = 1 - alpha * (1 - p);
			co_mat[i][i+1] = -1 * alpha * p;
			b_vec[i] = c * i;
		} else {
			co_mat[i][0] = -1 * alpha * (1 - p);
			co_mat[i][1] = -1 * alpha * p;
			co_mat[i][i] += 1;			// in case i == 0 or 1;
			b_vec[i] = K;
		}
	}

	return 0;
}

int solve_linear_eqs(int n, double co_mat[][n], double b_vec[], double J_vec[])
{
	gauss_elim(n, co_mat, b_vec);
	print_mat(n, co_mat, b_vec);
	back_elim(n, co_mat, b_vec, J_vec);
}

int gauss_elim(int n, double co_mat[][n], double b_vec[])
{
	int i, j;
	double mult;

	for (j = 0; j < n; j++) {
			
		part_pivot(n, co_mat, b_vec, j);
		print_mat(n, co_mat, b_vec);

		for (i = j + 1; i < n; i++) {

			mult = co_mat[i][j] / co_mat[j][j];

			for (int k = 0; k < n; k++)
				co_mat[i][k] = co_mat[i][k] - mult * co_mat[j][k];

			b_vec[i] = b_vec[i] - mult * b_vec[j];
		}
		
		print_mat(n, co_mat, b_vec);
	}

	return 0;
}

int back_elim(int n, double co_mat[][n], double b_vec[], double J_vec[])
{	
	int i, j;
	double tmp;

	J_vec[n-1] = b_vec[n-1] / co_mat[n-1][n-1];

	for (i = n - 2; i >= 0; i--) {
		tmp = b_vec[i];

		for (j = i + 1; j < n; j++)
			tmp -= co_mat[i][j] * J_vec[j];

		J_vec[i] = tmp / co_mat[i][i];	
	}
	
	return 0;
}

int part_pivot(int n, double co_mat[][n], double b_vec[], int j)
{
	int max = j;
	for (int i = j + 1; i < n; i++)
		if (fabs(co_mat[i][j]) > fabs(co_mat[max][j]))
			max = i;
	row_swap(n, co_mat, b_vec, max, j);
	return 0;
}

int row_swap(int n, double co_mat[][n], double b_vec[], int i, int j)
{
	double *tmp_vec = malloc(sizeof(double) * n);
	double tmp_b;

	for (int k = 0; k < n; k++) {
		tmp_vec[k] = co_mat[i][k];	
		co_mat[i][k] = co_mat[j][k];
	}

	for (int h = 0; h < n; h++)
		co_mat[j][h] = tmp_vec[h];	

	tmp_b = b_vec[i];
	b_vec[i] = b_vec[j];
	b_vec[j] = tmp_b;	

	free(tmp_vec);
	return 0;
}	

int print_mat(int n, double co_mat[][n], double b_vec[])
{
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%5.2f\t", co_mat[i][j]);	
		}	
		printf("b = \t%5.2f\n", b_vec[i]);
	}
	printf("\n");

	return 0;
}
