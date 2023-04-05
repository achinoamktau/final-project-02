#include "spkmeans_functions.h"


void printSpectralClustering(double **dots, int d, int n, int k ){
    double **T,**final_centroids, **initialCentroids;
    T = getT(dots,d,n,&k);
    int max_iter = 300;

    initialCentroids = deepCopy2DArray(T,k,k); /* defining the first k 'dots' (from T) as the initial  centroids. Like in Ex1. */

    final_centroids = kmeans(T, initialCentroids, k, k, n, max_iter);
    print_2d_array(final_centroids,k,k);

    /* frees: */
    free_contiguous_mat(T);
    free_contiguous_mat(initialCentroids);
    /* no need to free final_centroids because it points to the same address as initialCentroids */

}

int main(int argc, char** argv){   
    FILE* fp;
    int d, n;
    n = 0;
    d = 0;
    double** laplacian;
    double** centroids;
    double** adj_matrix;
    double** degree_matrix;
    Jacobi_output* jacobi_res;
    char sepereator;
    char* goal;
    if(argc < 3){
        printf("not enough arguments, sory");
        return 1;
    }
    fp = fopen(argv[2], "r");
    if(fp == NULL){
        printf("file not found");
        return 1;
    }
    goal = argv[1];  
    sepereator = fgetc(fp);
    while(sepereator != '\n' && sepereator != EOF){
        if(sepereator == ','){
            d++;
        }
        sepereator = fgetc(fp);
    }
    d++;
    rewind(fp);
    for (sepereator = getc(fp); sepereator != EOF; sepereator = getc(fp)){
        if (sepereator == '\n'){ // Increment count if this character is newline
            n = n + 1;
        }
    }
    centroids = allocateMatrix(n, d);
    rewind(fp);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < d; j++){
            fscanf(fp, "%lf,", &centroids[i][j]);
        
        }
    }
    fclose(fp);

    if(strcmp("jacobi", goal) == 0){
        jacobi_res = jacobi(centroids, n);
        print_list(jacobi_res->eigenValues, n);
        print_2d_array(jacobi_res->V, n, n);
        free_contiguous_mat(jacobi_res->V);
        free(jacobi_res->eigenValues);
        free(jacobi_res);
    }
    adj_matrix = wam(centroids, n, d);
    degree_matrix = ddg(adj_matrix, n);
    if(strcmp("wam", goal) == 0){
        print_matrix(adj_matrix, n, n);
    }
    else if(strcmp("ddg", goal) == 0){
        print_matrix(degree_matrix, n, n);

    }
    else if(strcmp("gl", goal) == 0){
        laplacian = gl(adj_matrix, degree_matrix, n);
        print_matrix(laplacian, n, n);
        free_contiguous_mat(laplacian);
    }
    else if(strcmp("spk",goal)==0){
        printSpectralClustering(centroids, n, d, 0);    //k=0
    }
    else if(strcmp("jacobi", goal) == 1){
        printf("invalid goal");
    }
    if(strcmp("gl", goal) != 0){
        free_contiguous_mat(adj_matrix);
        free_contiguous_mat(degree_matrix);
    }
    free_contiguous_mat(centroids);     //maybe we dont want it contigous
    return 1;  
}


