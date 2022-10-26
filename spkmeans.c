#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "spkmeans.h"

double **T;
int K;
int d;
int N;
void print_matrix(double **matrix);
void print_diag(double **A);
int assert1(double **DataPts);
int assert2(double *DataPts);
int assert3(myStruct *array);
void k_means(double** T);
int update(double** clusters, double** copy_clusters);
void new_clus(double** group_clus, double** clusters);
int nearest_clus(double** clusters, double* point);
void calc_new_clus(double** clusters, double** arr, double** group_clus);
void copy_clus(double** arr1, double** arr2);
void spk_case(double** DataPts);
void find_k(double** V, double** A);
void find_U(double** V, double** U);
void find_T(double** U, double** T);
void updateVA(double** V, double** A, myStruct* array);
int find_max_diff(double** A);
void wam_case(double** DataPts);
void calc_weigh_adj_mat(double** DataPts, double **weigh_adj_mat);
void ddg_case(double** DataPts);
void calc_diag_deg_mat(double** weigh_adj_mat, double** diag_deg_mat);
void lnorm_case(double** DataPts);
void sqrt_diag_deg_mat(double** diag_deg_mat);
void mul_W_left(double** weigh_adj_mat, double** diag_deg_mat, double** lnorm);
void copy_mat(double** weigh_adj_mat, double** lnorm);
void mul_W_right(double** weigh_adj_mat, double** diag_deg_mat, double** lnorm);
void calc_norm_laplacian(double** lnorm);
void jacobi_case(double** A);
void jacobi_algo(double** A, double** AA, double** P, double** V);
double calc_offA(double** A);
int check_diag_mat(double** A, double** V);
void find_s_c_ind_val(double** A, double** P, double* index, int r);
double sign(double phi);
void calc_P(double** P, double* index);
void mult_PAP(double** A, double** AA, double* index);
void calc_V(double** P, double** V, int i);
int convergence(double** A, double* offA, int i, int j);
void transposeA(double** V, double** AA);




int main(int argc, char **argv) {
    double **DataPts;
    char line[500];
    char* num;
    FILE *fptr;
    int i;

    if (argc != 3){
        printf("An Error Has Occured");
        return 0;
    }
    fptr = fopen(argv[2],"r");

    DataPts = (double**) malloc(1000*sizeof(double*));
    assert1(DataPts);
    assert(DataPts!=NULL);

    fgets(line, sizeof(line), fptr);
    d = 0;
    DataPts[0] = (double*)malloc(1000*sizeof(double));
    assert2(DataPts[0]);
    assert(DataPts[0]!=NULL);
    num = strtok(line, ",");
    while (num != NULL){
        DataPts[0][d] = atof(num);
        d++;
        num = strtok(NULL, ",");
    }

    DataPts[0] = (double*)realloc(DataPts[0], d*sizeof(double));
    assert2(DataPts[0]);
    assert(DataPts[0]!=NULL);

    N = 1;
    while (fgets(line, sizeof(line), fptr) != NULL){
        num = strtok(line, ",");
        DataPts[N] = (double*) malloc(d*sizeof(double));
        assert2(DataPts[N]);
        assert(DataPts[N]!=NULL);
        for (i=0 ; i<d && num != NULL; i++){
            DataPts[N][i] = atof(num);
            num = strtok(NULL, ",");
        }
        N++;
    }

    DataPts = (double**) realloc(DataPts, N*sizeof(double*));
    assert1(DataPts);
    assert(DataPts!=NULL);
    for (i = 0; i < N; i++){
        assert2(DataPts[i]);
        assert(DataPts[i]!=NULL);
    }

    if(!strcmp(argv[1],"spk")){
        spk_case(DataPts);
        if (K>=N || K ==0 || K==1){
            printf("An Error Has Occured");
            return 0;
        }
        k_means(T);
    } else if(!strcmp(argv[1],"wam")){
        wam_case(DataPts);
    } else if(!strcmp(argv[1],"ddg")){
        ddg_case(DataPts);
    } else if(!strcmp(argv[1],"lnorm")){
        lnorm_case(DataPts);
    } else if(!strcmp(argv[1],"jacobi")){
        jacobi_case(DataPts);
    } else {
        printf("Invalid Input!\n");
        exit(0);
    }

    fclose(fptr);
    return 0;

}

/* prints an NxN matrix */
void print_matrix(double **matrix){
    int i,j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            if (matrix[i][j] < 0.0 && matrix[i][j] > -0.00005){
                printf("%.4f",-0.0000);
            } else {


                printf("%.4f",matrix[i][j]);
            }

            if(j < N-1){
                printf("%s",",");


            }
        }
        if(i < N-1){
            printf("\n");

        }
    }
    printf("\n");
}

/* prints the diameter of an NxN matrix*/
void print_diag(double **A){
    int i;
    for (i = 0 ; i < N ; i++){
        if (A[i][i] < 0.0 && A[i][i] > -0.00005){
            printf("%.4f",0.0000);


        } else {
            printf("%.4f",A[i][i]);


        }
        if(i < N-1){
            printf(",");

        } else {


            printf("\n");
        }
    }
    
}

/* this function checks if the allocation of a double array (using malloc) was successful */
int assert1(double **DataPts){
    if (DataPts != NULL)
        return 0;
    else {


        printf("An Error Has Occured");
        return 0;


    }
}
/* this function checks if the allocation of a double array at index i (using malloc) was successful */
int assert2(double *DataPts){
    if (DataPts != NULL)
        return 0;
    else {


        printf("An Error Has Occured");
        return 0;
    }
}
/* this function checks if the allocation of an array of structs (using malloc) was successful */
int assert3(myStruct *array){
    if (array != NULL)
        return 0;
    else {
        printf("An Error Has Occured");
        return 0;



    }
}

/*
 * the first case: Goal = spk
 *   - calculate the wam,ddg and lnorm matrices
 *   - find the eigen vectors and the eigen values 
 *   - if needed, find k
 *   - calculate the U,T matrices to implement the kmeans algorithm on T
 */
void spk_case(double** DataPts){
    int i;
    double **weigh_adj_mat;
    double **diag_deg_mat;
    double **lnorm;
    double **AA, **P, **V;
    double **U;
    weigh_adj_mat = (double**) malloc(N*sizeof(double*));
    assert1(weigh_adj_mat);
    assert(weigh_adj_mat!=NULL);
    diag_deg_mat = (double**) calloc(N, sizeof(double*));
    assert1(diag_deg_mat);
    assert(diag_deg_mat!=NULL);
    lnorm = (double**) malloc(N*sizeof(double*));
    assert1(lnorm);
    assert(lnorm!=NULL);
    AA = (double**) malloc(N*sizeof(double*));
    assert1(AA);
    assert(AA!=NULL);
    P = (double**) calloc(N, sizeof(double*));
    assert1(P);
    assert(P!=NULL);
    V = (double**) malloc(N*sizeof(double*));
    assert1(V);
    assert(V!=NULL);
    for (i = 0; i < N; i++){
        weigh_adj_mat[i] = (double*) malloc(N*sizeof(double));
        assert2(weigh_adj_mat[i]);
        assert(weigh_adj_mat[i]!=NULL);
        diag_deg_mat[i] = (double*) calloc(N, sizeof(double));
        assert2(diag_deg_mat[i]);
        assert(diag_deg_mat[i]!=NULL);
        lnorm[i] = (double*) malloc(N*sizeof(double));
        assert2(lnorm[i]);
        assert(lnorm[i]!=NULL);
        AA[i] = (double*) malloc(N*sizeof(double));
        assert2(AA[i]);
        assert(AA[i]!=NULL);
        P[i] = (double*) calloc(N, sizeof(double));
        assert2(P[i]);
        assert(P[i]!=NULL);
        V[i] = (double*) malloc(N*sizeof(double));
        assert2(V[i]);
        assert(V[i]!=NULL);
    }
    calc_weigh_adj_mat(DataPts, weigh_adj_mat);
    calc_diag_deg_mat(weigh_adj_mat, diag_deg_mat);
    sqrt_diag_deg_mat(diag_deg_mat);
    mul_W_left(weigh_adj_mat,diag_deg_mat, lnorm);
    copy_mat(weigh_adj_mat, lnorm);
    mul_W_right(weigh_adj_mat,diag_deg_mat, lnorm);
    calc_norm_laplacian(lnorm);
    jacobi_algo(lnorm, AA, P, V);
    find_k(V, lnorm);
    U = (double**) malloc(N*sizeof(double*));
    assert1(U);
    assert(U!=NULL);
    T = (double**) malloc(N*sizeof(double*));
    assert1(T);
    assert(T!=NULL);
    for (i = 0; i < N; i++){
        U[i] = (double*) malloc(K*sizeof(double));
        assert2(U[i]);
        assert(U[i]!=NULL);
        T[i] = (double*) malloc(K*sizeof(double));
        assert2(T[i]);
        assert(T[i]!=NULL);
    }
    find_U(V, U);
    find_T(U, T);
    d = K;
    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
        free(weigh_adj_mat[i]);
        free(diag_deg_mat[i]);
        free(lnorm[i]);
        free(AA[i]);
        free(P[i]);
        free(V[i]);
        free(U[i]);
    }
    free(DataPts);
    free(weigh_adj_mat);
    free(diag_deg_mat);
    free(lnorm);
    free(AA);
    free(P);
    free(V);
    free(U);
}

/* function that compares between struct variable types */
int Scomparator(const void *v1, const void *v2)
{
    const myStruct *p1 = (myStruct *)v1;
    const myStruct *p2 = (myStruct *)v2;
    if (p1->eigen_val < p2->eigen_val)
        return -1;
    else if (p1->eigen_val > p2->eigen_val)
        return +1;
    else if (p1->index_eigen_val < p2->index_eigen_val)
        return -1;
    else if (p1->index_eigen_val > p2->index_eigen_val)
        return +1;
    else
        return 0;
}

/* this function calculates K in the case of goal == spk and K==0.
 *  -use qsort algorithm to sort the eigen values to achieve a stable sort using the previous function (Scomparator) as a comparator
 *  -implement the Eigengap Heuristic algorithm
 *  -calculate the maximum difference between the first n\2 eigen values , this would be the value we're looking for (K)
 */
void find_k(double** V, double** A){
    int i,j;
    myStruct *array;
    array = (myStruct*)calloc(N, sizeof(myStruct));
    assert3(array);
    assert(array!=NULL);
    for (i = 0; i < N; i++) {
        array->index_eigen_val = i+1;
        array->eigen_val = A[i][i];
        array->eigen_vec = (double*) malloc(N*sizeof(double));
        assert2(array->eigen_vec);
        assert(array->eigen_vec!=NULL);
        for (j = 0; j < N; j++){
            array->eigen_vec[j] = V[j][i];
        }
        array++;
    }
    array = array -N;
    qsort(array, N, sizeof(myStruct), Scomparator);
    updateVA(V, A, array);
    if (K == 0)
        K = find_max_diff(A);
    for (i = 0; i < N; i++){
        free(array[i].eigen_vec);
    }
    free(array);
}

/* this function updates the martices V,A upon sorting the eigen values*/
void updateVA(double** V, double** A, myStruct* array){
    int i,j;
    for (i = 0; i < N; i++) {
        A[i][i] = array[i].eigen_val;
        for (j = 0; j < N; j++){
            V[j][i] = array[i].eigen_vec[j];
        }
    }
}

/* find the index that gives the maximum difference in the eigen values (after sorting them) */
int find_max_diff(double** A) {
    int i, index;
    double max;
    max = fabs(A[0][0] -A[1][1]);
    index = 1;

    for (i = 1; i <= N/2 -1; i++){
        if (N>2){
            if (max < fabs(A[i][i] -A[i+1][i+1])){
                index = i +1;
                max = fabs(A[i][i] -A[i+1][i+1]);
            }
        }
    }
    return index;
}

/* calculates U from V, U is an NxK matrix */
void find_U(double** V, double** U){
    int i,j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            U[i][j] = V[i][j];
        }
    }
}

/* calculate the matrix T from U by renormalizing each of U's rows */
void find_T(double** U, double** T){
    int i,j;
    double num;

    for (i = 0; i < N; i++) {
        num = 0;
        for (j = 0; j < K; j++) {
            num = num + pow(U[i][j],2);
        }
        if (num ==0){ 
              for (j = 0; j < K; j++) {
            T[i][j] = 0;
        }}
        else {

        num = pow(num, 0.5);
        for (j = 0; j < K; j++) {
            T[i][j] = U[i][j] / num;
        }}
    }
}

/*
 * the second case: goal = wam
 *     - calculate the wam matrix
 */
void wam_case(double** DataPts){
    int i;
    double **weigh_adj_mat;
    weigh_adj_mat = (double**) malloc(N*sizeof(double*));
    assert1(weigh_adj_mat);
    assert(weigh_adj_mat!=NULL);

    for (i = 0; i < N; i++){
        weigh_adj_mat[i] = (double*) malloc(N*sizeof(double));
        assert2(weigh_adj_mat[i]);
        assert(weigh_adj_mat[i]!=NULL);
    }

    calc_weigh_adj_mat(DataPts, weigh_adj_mat);
    print_matrix(weigh_adj_mat);
    

    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
        free(weigh_adj_mat[i]);
    }
    free(DataPts);
    free(weigh_adj_mat);

}

/* calculates the wam matrix using the given formula through the data points */
void calc_weigh_adj_mat(double** DataPts, double** weigh_adj_mat){
    int i,j,t;
    double candidate,num;

    for (i = 0; i < N; i++){
        for (j = i+1; j < N; j++) {
            candidate = 0;
            for (t = 0; t < d; t++) {
                num = DataPts[i][t] - DataPts[j][t];
                candidate = candidate + pow(num, 2);
            }

            candidate = pow(candidate, 0.5);
            candidate = candidate / 2.0;
            candidate = exp(-1 * candidate);
            weigh_adj_mat[i][j] = candidate;
            weigh_adj_mat[j][i] = candidate;
        }
    }

    for (i = 0; i < N; i++){
        weigh_adj_mat[i][i] = 0;
    }
}


/*
 * the third case: goal = ddg
 *     -calculate the ddg matrix
 */
void ddg_case(double** DataPts){
    int i;
    double **weigh_adj_mat;
    double **diag_deg_mat;
    weigh_adj_mat = (double**) malloc(N*sizeof(double*));
    assert1(weigh_adj_mat);
    assert(weigh_adj_mat!=NULL);
    diag_deg_mat = (double**) calloc(N, sizeof(double*));
    assert1(diag_deg_mat);
    assert(diag_deg_mat!=NULL);
    for (i = 0; i < N; i++){
        weigh_adj_mat[i] = (double*) malloc(N*sizeof(double));
        assert2(weigh_adj_mat[i]);
        assert(weigh_adj_mat[i]!=NULL);
        diag_deg_mat[i] = (double*) calloc(N, sizeof(double));
        assert2(diag_deg_mat[i]);
        assert(diag_deg_mat[i]!=NULL);
    }
    calc_weigh_adj_mat(DataPts, weigh_adj_mat);
    calc_diag_deg_mat(weigh_adj_mat, diag_deg_mat);
    print_matrix(diag_deg_mat);
    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
        free(weigh_adj_mat[i]);
        free(diag_deg_mat[i]);
    }
    free(DataPts);
    free(weigh_adj_mat);
    free(diag_deg_mat);
}

/* calculates the ddg matrix using the wam matrix */
void calc_diag_deg_mat(double** weigh_adj_mat, double** diag_deg_mat){
    int i,j;
    double num;
    for (i = 0; i < N; i++) {
        num = 0;
        for (j = 0; j < N; j++) {
            num = num + weigh_adj_mat[i][j];
        }
        diag_deg_mat[i][i] = num;
    }
}


/*
 * the fourth case: goal = lnorm
 *     -calculate the lnorm matrix
 */
void lnorm_case(double** DataPts){
    int i;
    double **weigh_adj_mat;
    double **diag_deg_mat;
    double **lnorm;
    weigh_adj_mat = (double**) malloc(N*sizeof(double*));
    assert1(weigh_adj_mat);
    assert(weigh_adj_mat!=NULL);
    diag_deg_mat = (double**) calloc(N, sizeof(double*));
    assert1(diag_deg_mat);
    assert(diag_deg_mat!=NULL);
    lnorm = (double**) malloc(N*sizeof(double*));
    assert1(lnorm);
    assert(lnorm!=NULL);
    for (i = 0; i < N; i++){
        weigh_adj_mat[i] = (double*) malloc(N*sizeof(double));
        assert2(weigh_adj_mat[i]);
        assert(weigh_adj_mat[i]!=NULL);
        diag_deg_mat[i] = (double*) calloc(N, sizeof(double));
        assert2(diag_deg_mat[i]);
        assert(diag_deg_mat[i]!=NULL);
        lnorm[i] = (double*) malloc(N*sizeof(double));
        assert2(lnorm[i]);
        assert(lnorm[i]!=NULL);
    }
    calc_weigh_adj_mat(DataPts, weigh_adj_mat);
    calc_diag_deg_mat(weigh_adj_mat, diag_deg_mat);
    sqrt_diag_deg_mat(diag_deg_mat);
    mul_W_left(weigh_adj_mat,diag_deg_mat, lnorm);
    copy_mat(weigh_adj_mat, lnorm);
    mul_W_right(weigh_adj_mat,diag_deg_mat, lnorm);
    calc_norm_laplacian(lnorm);
    print_matrix(lnorm);
    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
        free(weigh_adj_mat[i]);
        free(diag_deg_mat[i]);
        free(lnorm[i]);
    }
    free(DataPts);
    free(weigh_adj_mat);
    free(diag_deg_mat);
    free(lnorm);
}

/* calculates D^(-0.5) from D matrix */
void sqrt_diag_deg_mat(double** diag_deg_mat){
    int i;
    for (i = 0; i < N; i++){
        diag_deg_mat[i][i] = pow(diag_deg_mat[i][i], -0.5);
    }
}

/* calculates D^0.5W */
void mul_W_left(double** weigh_adj_mat, double** diag_deg_mat, double** lnorm){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            lnorm[i][j] = diag_deg_mat[i][i]*weigh_adj_mat[i][j];
        }
    }
}

/* this function copies NxN matrices */
void copy_mat(double** weigh_adj_mat, double** lnorm){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            weigh_adj_mat[i][j] = lnorm[i][j];
        }
    }
}

/* calculates D^(-0.5)WD^(-0.5) */
void mul_W_right(double** weigh_adj_mat, double** diag_deg_mat, double** lnorm){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            lnorm[i][j] = diag_deg_mat[j][j]*weigh_adj_mat[i][j];
        }
    }
}

/* calculates I - D^(-0.5)WD^(-0.5) */
void calc_norm_laplacian(double** lnorm){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j){
                lnorm[i][i] = 1 - lnorm[i][i];
            } else {
                lnorm[i][j]  = -1*lnorm[i][j];
            }
        }
    }
}

/*
 * the fifth case: goal = jacobi
*      - calculates the eigen values and the eigen vectors
 *     
 * note: the input of this function is a symmetrical matrix
 */
void jacobi_case(double** A){
    int i;
    double **AA, **P, **V;
    AA = (double**) malloc(N*sizeof(double*));
    assert1(AA);
    assert(AA!=NULL);
    P = (double**) calloc(N, sizeof(double*));
    assert1(P);
    assert(P!=NULL);
    V = (double**) malloc(N*sizeof(double*));
    assert1(V);
    assert(V!=NULL);
    for (i = 0; i < N; i++){
        AA[i] = (double*) malloc(N*sizeof(double));
        assert2(AA[i]);
        assert(AA[i]!=NULL);
        P[i] = (double*) calloc(N, sizeof(double));
        assert2(P[i]);
        assert(P[i]!=NULL);
        V[i] = (double*) malloc(N*sizeof(double));
        assert2(V[i]);
        assert(V[i]!=NULL);
    }
    jacobi_algo(A, AA, P, V);
    print_diag(A);
    print_matrix(V);
    for(i = 0 ; i < N  ; i++){
        free(A[i]);
        free(AA[i]);
        free(P[i]);
        free(V[i]);
    }
    free(A);
    free(AA);
    free(P);
    free(V);

}

/*
 * this function implements the jacobi algorithm:
 *    - find the max value out of the values "above" the diagonal in A
 *    - calculate the P matrix using s,c
 *    - calculate AA using A and P
 *    - calculate matrix V
 */
void jacobi_algo(double **A, double **AA, double **P, double **V){
    int i, flag;
    double *index, offA;
    i = 0;
    flag = 1;
    index = (double*) malloc(4*sizeof(double));
    assert2(index);
    assert(index!=NULL);
    offA = calc_offA(A);
    index[0] = offA;
    copy_mat(AA, A);
    flag = check_diag_mat(A, V);
    while(i < 100 && flag ){
        find_s_c_ind_val(A, P, index, i);
        calc_P(P, index);
        calc_V(P, V, i);
        mult_PAP(A, AA, index);
        flag = convergence(A, &offA, (int)index[0], (int)index[1]);
        copy_mat(A, AA);
        i++;

    }
    free(index);
}

/* calculates the offA using the given formula */
double calc_offA(double** A){
    int i,j;
    double num, diagonal;
    num = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            num += pow(A[i][j], 2);
        }
    }
   

    diagonal = 0;
    for (i = 0; i < N; i++) {
        diagonal += pow(A[i][i], 2);
    }
    return num-diagonal;
}

/* checks if the input matrix is symmetrical (this is used for the jacobi algorithm) */
int check_diag_mat(double** A, double** V){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = i+1; j < N; j++) {
            if (A[i][j] != 0.0){
                return 1;
            }
        }
    }

    for (i = 0; i < N; i++) {
        V[i][i] = 1;
    }
    return 0;
}

/* finds the index of the max value out of the values above the diagonal in A and the s,c values */
void find_s_c_ind_val(double** A, double** P, double* index, int r){
    int i, j, index1, index2;
    double phi, t, c, s;
    double max1;
    index1 = 0;
    index2 = 0;
    max1 = 0.0;
    if (r != 0){
        P[(int)index[0]][(int)index[1]] = 0;
        P[(int)index[1]][(int)index[0]] = 0;
    }
    for (i = 0; i < N; i++) {
        for (j = i+1; j < N; j++) {
            if (i == 0 && j == 1){
                index1 = i;
                index2 = j;
                max1 = pow(pow(A[i][j],2),0.5);

            } else if (pow(pow(A[i][j],2),0.5) > max1){
                index1 = i;
                index2 = j;
                max1 = pow(pow(A[i][j],2),0.5);
            }
        }
    }
    index[0] = (double)index1;
    index[1] = (double)index2;
    phi = (A[index2][index2]-A[index1][index1]) / (2.0*A[index1][index2]);
    t = sign(phi);
    c = 1.0 / (pow(pow(t,2) +1, 0.5));
    s = c*t;
    index[2] = c;
    index[3] = s;
}

/* calculates the value of phi */
double sign(double phi){
    if (phi >= 0.0){
        return 1 / (pow(pow(phi,2),0.5) + pow(pow(phi,2) +1, 0.5));
    } else {
        return -1 / (pow(pow(phi,2),0.5) + pow(pow(phi,2) +1, 0.5));
    }
}

/* calculates the P matrix */
void calc_P(double** P, double* index){
    int i;
    for (i = 0; i < N; i++){
        P[i][i] = 1;
    }
    P[(int)index[0]][(int)index[0]] = index[2];
    P[(int)index[1]][(int)index[1]] = index[2];
    P[(int)index[0]][(int)index[1]] = index[3];
    P[(int)index[1]][(int)index[0]] = -1*index[3];
}

/* calculates the matrix PAP using the given formula */
void mult_PAP(double** A, double** AA, double* index){
    int i,j,r;
    double c,s;
    i = index[0];
    j = index[1];
    c = index[2];
    s = index[3];
    for (r = 0; r < N; r++) {
        if(r != i && r != j){
            AA[r][i] = c*A[r][i] -s*A[r][j];
            AA[i][r] = c*A[r][i] -s*A[r][j];

            AA[r][j] = c*A[r][j] +s*A[r][i];
            AA[j][r] = c*A[r][j] +s*A[r][i];

        } else if(r == i){
            AA[r][r] = pow(c,2)*A[i][i] +pow(s,2)*A[j][j] -2*c*s*A[i][j];

        } else {
            AA[r][r] = pow(s,2)*A[i][i] +pow(c,2)*A[j][j] +2*c*s*A[i][j];
        }
    }
    AA[i][j] = 0;
    AA[j][i] = 0;
}

/* calculates V */
void calc_V(double** P, double** V, int t){
    if (t == 0){
        copy_mat(V, P);
    } else {
        int i,j,r;
        double num;
        double** tmp;
        tmp = (double**) malloc(N*sizeof(double*));
        assert1(tmp);
        assert(tmp!=NULL);

        for (i = 0; i < N; i++){
            tmp[i] = (double*) malloc(N*sizeof(double));
            assert2(tmp[i]);
            assert(tmp[i]!=NULL);
        }

        for (i = 0; i < N; i++){
            for (j = 0; j < N; j++) {
                num = 0;
                for (r = 0; r < N; r++) {
                    num = num + V[i][r] * P[r][j];
                }
                tmp[i][j] = num;
            }
        }

        copy_mat(V, tmp);
        for(i = 0 ; i < N  ; i++){
            free(tmp[i]);
        }
        free(tmp);
    }
}

/* checks if convergence occurs */
int convergence(double** A, double* offA, int i, int j){
    double offAA;
    offAA = *offA -2*pow(A[i][j],2);
    if(*offA - offAA <= pow(10,-5)){
        return 0;
    }
    *offA = offAA;
    return 1;
}

/* calculates the transpose of the V matrix */
void transposeA(double** V, double** AA){
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            AA[i][j] = V[j][i];
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            V[i][j] = AA[i][j];
        }
    }
}

/*
 * this function implements the kmeans algorithm from HW1
 */
void k_means(double** DataPts){
    int i;
    int j;
    int t;

    double **clusters;
    double **copy_clusters;
    double **group_clus;

    clusters = (double**) malloc(K*sizeof(double*));
    assert1(clusters);
    assert(clusters!=NULL);
    for (i = 0; i < K; i++){
        clusters[i] = (double*) malloc(d*sizeof(double));
        assert2(clusters[i]);
        assert(clusters[i]!=NULL);
    }

    for (i=0; i<K; i++){
        for (j = 0; j<d; j++){
            clusters[i][j] = DataPts[i][j];
        }
    }
    copy_clusters = (double**) malloc(K*sizeof(double*));
    assert1(copy_clusters);
    assert(copy_clusters!=NULL);
    for (i = 0; i< K; i++){
        copy_clusters[i] = (double*) malloc(d*sizeof(double));
        assert2(copy_clusters[i]);
        assert(copy_clusters[i]!=NULL);
    }
    group_clus = (double**) malloc(K*sizeof(double*));
    assert1(group_clus);
    assert(group_clus!=NULL);
    for (i = 0; i< K; i++){
        group_clus[i] = (double*) malloc((d+1)*sizeof(double));
        assert2(group_clus[i]);
        assert(group_clus[i]!=NULL);
    }
    i = 0;
    t = 1;
    while ((i < 300)&&(t == 1)){
        copy_clus(clusters, copy_clusters);
        calc_new_clus(clusters, DataPts, group_clus);
        new_clus(group_clus, clusters);
        t = update(clusters, copy_clusters);
        i++;
    }
    for(i = 0; i < K; i++){
        for(j = 0; j < d; j++){
            if (clusters[i][j] < 0.0 && clusters[i][j] > -0.00005){
                printf("%.4f",0.0000);
            } else {
                printf("%.4f",clusters[i][j]);
            }
            if(j < K-1){
                printf("%s",",");
            }
        }
        if(i < K-1){
            printf("\n");
        }
    }
    for(i = 0 ; i < N  ; i++){
        free(DataPts[i]);
    }
    free(DataPts);

    for(i = 0 ; i < K ; i++){
        free(group_clus[i]);
        free(copy_clusters[i]);
        free(clusters[i]);
    }
    free(group_clus);
    free(copy_clusters);
    free(clusters);
}


/* checks to which cluster the given point is nearest*/
int nearest_clus(double** clusters, double* point){
    double min;
    int index;
    int i;
    int j;
    double candidate;
    min=0;
    index=0;
    for (i=0; i<K; i++){
        candidate = 0;
        for (j=0; j<d; j++){
            candidate = candidate + (clusters[i][j] - point[j])*(clusters[i][j] - point[j]);
        }
        if ((candidate < min) || (i == 0)){
            min = candidate;
            index = i;
        }
    }
    return index;
}

/* checks if the two clusters are any different, if they differ in atleast one value the function returns 1, otherwise it returns 0 */
int update(double** clusters, double** copy_clusters){
    int i,j;
    for (i = 0; i < K; i++){
        for (j = 0; j < d; j++){
            if (clusters[i][j] != copy_clusters[i][j]){
                return 1;
            }
        }
    }
    return 0;
}

/* copies arr1 to arr2 */
void copy_clus(double** arr1, double** arr2){
    int i,j;
    for (i = 0; i < K; i++){
        for (j = 0; j < d; j++){
            arr2[i][j] = arr1[i][j];
        }
    }
}

/* calculates the new clusters */
void new_clus(double** group_clus, double** clusters){
    int i,j;
    for (i = 0; i < K; i++){
        for (j = 0; j < d; j++){
            clusters[i][j] = (group_clus[i][j]) / (group_clus[i][d]);
        }
    }
}



/* calculates the new clusters after we assign a cluster to each point */
void calc_new_clus(double** clusters, double** arr, double** group_clus){
    int i;
    int j;
    int index;
    for (i = 0; i < K; i++){
        for (j = 0; j < d+1; j++){
            group_clus[i][j]=0;
        }
    }
    for (i = 0; i < N; i++){
        index = nearest_clus(clusters, arr[i]);
        for (j = 0; j < d; j++){
            group_clus[index][j] = group_clus[index][j] + arr[i][j];
        }
        group_clus[index][d] += 1.0;
    }
}




