#ifndef _MAIN
#define _MAIN

// gauss array
int GP_1D;                                    // 1方向のガウス点数
int GP_2D;                                    // 2次元のガウス点数
int GP_3D;                                    // 3次元のガウス点数
int GP_ON_ELEMENT;                            // 要素内のガウス点数

// for s-IGA
double E;             // ヤング率
double nu;            // ポアソン比
int Total_mesh;

int MAX_ORDER = 0;         // 基底関数の次数の最大値 + 1 (各要素がサポートする制御点数)
int MAX_CP = 0;            // 
int MAX_KNOT = 0;          // 
int MAX_NO_CP_ON_ELEMENT;  // 
int MAX_KIEL_SIZE;         // 

int D_MATRIX_SIZE;         //
int N_STRAIN;              //
int N_STRESS;              //
int MAX_K_WHOLE_SIZE;      //
int K_Whole_Size;          //

int Total_connectivity;      //
int Total_connectivity_glo;  //
int Total_edge;              //
int Total_edge_glo;          //

// file pointer
FILE *fp;

#endif