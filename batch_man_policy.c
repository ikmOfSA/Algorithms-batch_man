#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int forward(int n, double J_vec[], double J_vec_new[], double alpha, double p, double K, double c);
int gen_co_mat(int n, double co_mat[][n], double b_vec[], double alpha, double p, double K, double c, int m);
int solve_linear_eqs(int n, double co_mat[][n], double b_vec[], double J_vec[]);
int gauss_elim(int n, double co_mat[][n], double b_vec[]);
int back_elim(int n, double co_mat[][n], double b_vec[], double J_vec[]);
int part_pivot(int n, double co_mat[][n], double b_vec[], int j);
int row_swap(int n, double co_mat[][n], double b_vec[], int i, int j);
int print_mat(int n, double co_mat[][n], double b_vec[]);
int improv(int n, double J_vec[], double alpha, double p, double K, double c, int m);
int run(int n, double alpha, double p, double K, double c, FILE * fp, int ite);
int pre_run(int n, int ite);

int main()
{
	int n[3] = {5, 10, 30};
	int ite[5] = {10, 30, 50, 70, 90};
	for (int i = 0; i < 3; i++)
		for ( int j = 0; j < 5; j++)
			pre_run(n[i], ite[j]);
	return 0;
}

int pre_run(int n, int ite)
{

	double alpha, p, K, c;

	char name[50];
	sprintf(name, "outcome_n_%d_ite_%d.cvs", n, ite);

	FILE *fp = fopen(name, "ab");

	fprintf(fp, "n,alpha,p,K,c,pol_m,val_m,");
	for (int x = 0; x < n; x++)
		fprintf(fp, "pol_%d,", x);

	for (int y = 0; y < n; y++)
		fprintf(fp, "val_%d,", y);

	fprintf(fp, "pol_avg, val_avg, avg_diff");
	fprintf(fp, "\n");

	for (alpha = 0.95; alpha > 0.5; alpha -= 0.1)
		for (p = 0.75; p > 0; p -= 0.25)
			for (K = 10; K < 500; K += 10)
				for (c = 1; c < 10; c++)
					run(n, alpha, p, K, c, fp, ite);
	fclose(fp);
	return 0;
}

int run(int n, double alpha, double p, double K, double c, FILE *fp, int ite)
{	
	int m = (n/2 < 1) ? 1 : n/2;			//initial policy threshold; when reaching m, process; states from 1 to n-1;

	double (* co_mat)[n] = calloc(n * n, sizeof(double)); //coefficient matrix of linear equations;
	double *b_vec = calloc(n, sizeof(double));	//output vector of linear equations;
	double *J_vec = calloc(n, sizeof(double));	//cost vector, index being state;

	int new_m;		//improved policy
	for (int k = 0; k < n - 1; k++)		//in case there is an infinite loop of the possible value of m, all with the same cost.
	{
		gen_co_mat(n, co_mat, b_vec, alpha, p, K, c, m);	//generate coeffcient matrix according to current policy
		solve_linear_eqs(n, co_mat, b_vec, J_vec);		//solve linear equations to get cost vector;
		new_m = improv(n, J_vec, alpha, p, K, c, m);		//get new policy, notice all the policy are thresholds.
		if (m == new_m)			//when improved policy is the same, optimal.
			break;
		m = new_m;
	}

	double *J_vec_val = calloc(n, sizeof(double));	//current cost vector;
	double *J_vec_val_new = calloc(n, sizeof(double));	//next cost vector, i.e., the improvement step
	
	for (int j = 0; j < ite; j++)
		forward(n, J_vec_val, J_vec_val_new, alpha, p, K, c);	//improvement step

	int m_val = forward(n, J_vec_val, J_vec_val_new, alpha, p, K, c);	//run again just for getting m, which is the policy threshold;
	
	fprintf(fp, "%d,%f,%f,%f,%f,%d,%d,", n, alpha, p, K, c, m, m_val);

	double pol_sum, pol_avg, val_sum, val_avg, avg_diff;

	pol_sum = 0;
	for (int x = 0; x < n; x++) {
		fprintf(fp, "%f,", J_vec[x]);
		pol_sum += J_vec[x];
	}
	pol_avg = pol_sum / n;
	
	val_sum = 0;
	for (int y = 0; y < n; y++) {
		fprintf(fp, "%f,", J_vec_val[y]);
		val_sum += J_vec_val[y];
	}
	val_avg = val_sum / n;	

	avg_diff = fabs((pol_avg - val_avg) / pol_avg);
	fprintf(fp, "%f,%f,%f", pol_avg, val_avg, avg_diff);
	
	fprintf(fp, "\n");

	return 0;    
}

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

/*get improved new policy*/
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

/*generate coefficient matrix according to current policy*/
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

/*code that solve the linear equation system, when there are single result;*/
int solve_linear_eqs(int n, double co_mat[][n], double b_vec[], double J_vec[])
{
	gauss_elim(n, co_mat, b_vec);
	//print_mat(n, co_mat, b_vec);
	back_elim(n, co_mat, b_vec, J_vec);
}

/*Gauss matrix elimination*/
int gauss_elim(int n, double co_mat[][n], double b_vec[])
{
	int i, j;
	double mult;

	for (j = 0; j < n; j++) {
			
		part_pivot(n, co_mat, b_vec, j);
		//print_mat(n, co_mat, b_vec);

		for (i = j + 1; i < n; i++) {

			mult = co_mat[i][j] / co_mat[j][j];

			for (int k = 0; k < n; k++)
				co_mat[i][k] = co_mat[i][k] - mult * co_mat[j][k];

			b_vec[i] = b_vec[i] - mult * b_vec[j];
		}
		
		//print_mat(n, co_mat, b_vec);
	}

	return 0;
}

/*obtain linear equation system result from echelon matrix*/
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

/*partial_pivot trick that enhence the gauss elimination*/
int part_pivot(int n, double co_mat[][n], double b_vec[], int j)
{
	int max = j;
	for (int i = j + 1; i < n; i++)
		if (fabs(co_mat[i][j]) > fabs(co_mat[max][j]))
			max = i;
	row_swap(n, co_mat, b_vec, max, j);
	return 0;
}

/*swap two row of the matrix*/
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

/*print matrix in format*/
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
