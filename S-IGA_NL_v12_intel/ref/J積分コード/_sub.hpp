#ifndef _SUB
#define _SUB

#include <cstdio>

// gauss array
extern int GP_1D;                                    // 1方向のガウス点数
extern int GP_2D;                                    // 2次元のガウス点数
extern int GP_3D;                                    // 3次元のガウス点数
extern int GP_ON_ELEMENT;                            // 要素内のガウス点数

// for s-IGA
extern double E;             // ヤング率
extern double nu;            // ポアソン比
extern int Total_mesh;

extern int MAX_ORDER;        // 基底関数の次数の最大値 + 1
extern int MAX_CP;
extern int MAX_KNOT;
extern int MAX_NO_CP_ON_ELEMENT;
extern int MAX_KIEL_SIZE;

extern int D_MATRIX_SIZE;
extern int N_STRAIN;
extern int N_STRESS;
extern int MAX_K_WHOLE_SIZE;
extern int K_Whole_Size;

extern int Total_connectivity;
extern int Total_connectivity_glo;
extern int Total_edge;
extern int Total_edge_glo;

// file pointer
extern FILE *fp;

#endif