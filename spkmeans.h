#ifndef SPKMEANS_C_SPKMEANS_H
#define SPKMEANS_C_SPKMEANS_H

typedef struct myStruct{
    int index_eigen_val;
    double eigen_val;
    double *eigen_vec;
} myStruct;

void print_matrix(double **matrix);
void print_diag(double **A);
int assert1(double **DataPts);
int assert2(double *DataPts);
int assert3(myStruct *array);

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

int update(double** clusters, double** copy_clusters);
void new_clus(double** group_clus, double** clusters);
int nearest_clus(double** clusters, double* point);
void calc_new_clus(double** clusters, double** arr, double** group_clus);
void copy_clus(double** arr1, double** arr2);

#endif