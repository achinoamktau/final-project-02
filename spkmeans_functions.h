#ifndef SPKMEANS_FUNCTIONS
#define SPKMEANS_FUNCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* structs: */
typedef struct ROTATION_MATRIX{
    double c;
    double s;
    int i;
    int j;
} rotation_mat;

typedef struct Jacobi_output{
    double *eigenValues;
    double **V;
} Jacobi_output;

double** alloc_nXm_matrix(int n, int m);
double** alloc_nXn_matrix(int n);
void free_contiguous_mat(double **mat);
double** create_identity_matrix(int n);
double** wam(double **dots, int d, int n);
double** ddg(double **W, int n);
double **gl(double **D, double **W, int n);
int *find_off_diag_max_abs_val(double **A, int n);
double **calc_A_tag(double **A_tag, double **A, int n, rotation_mat *P);
double calc_of_f_square(double **A, int n);
void multiply_rotation_matrix(double **V, rotation_mat *P,int n);
Jacobi_output* jacobi(double **A, int n);
double** create2DfromJacobi(Jacobi_output *jacobi_output, int n);
rotation_mat *calc_rotation_mat(rotation_mat *P,double **A, int n);
int compare_eigenStruct(const void *a, const void *b);
double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer);
void free2DMalloc(double **points,int n);

/* Kmeans functions: */
double calc_distance(const double *dot, const double *centroid, int d);
int find_nearest_centroid(double *dot, double **centroids, int k, int d);
void update_centroids(double **dots,const int *dots_location, double **new_centroids, int n, int d, int k);
int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num);
void print_2d_array(double **array, int row_num, int col_num);
void print_2d_array_transpose(double **array, int row_num, int col_num);
double **deepCopy2DArray(double **A, int row_num, int col_num);
void print_list(double *, int);
double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter);

#endif