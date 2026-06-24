#ifndef _MAIN
#define _MAIN

// for s-IGA
double E;             // ヤング率
double nu;            // ポアソン比
int Total_mesh;

int MAX_ORDER = 0;    // 基底関数の次数の最大値 + 1
int MAX_CP = 0;
int MAX_KNOT = 0;
int MAX_NO_CP_ON_ELEMENT;
int MAX_KIEL_SIZE;

int D_MATRIX_SIZE;
int N_STRAIN;
int N_STRESS;
int MAX_K_WHOLE_SIZE;
int K_Whole_Size;

int Total_connectivity;
int Total_connectivity_glo;
int Total_edge;
int Total_edge_glo;

double lame_lambda;
double lame_mu;

#endif