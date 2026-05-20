#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


// header
#include "_header.hpp"
#include "_sub.hpp"

using namespace std;

// memory allocation
void Allocation(const int num, information *info)
{
	if (num == 0)
	{
		info->Total_Knot_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));             // メッシュまでのノット数(メッシュ内のノット数は含まない)
		info->Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));			     // 各メッシュ上のパッチ数 （メッシュ内のパッチ数は含まない）   
		info->Total_Patch_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		     // メッシュまでのパッチ数（メッシュ内のパッチ数は含まない）
		info->Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	     // 各メッシュ上のコントロールポイント数
		info->Total_Control_Point_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));    // メッシュまでのコントロールポイント数(メッシュ内のコントロールポイント数は含まない)
		info->Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		     // 各メッシュ上の要素数
		info->Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		     // メッシュまでの要素数(メッシュ内の要素数は含まない)
		info->real_Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	     // 各メッシュ上の要素数(多重ノットを考慮した要素数) real_Total_Element_on_mesh[MAX_N_MESH]
		info->real_Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));     // 多重ノットを考慮したメッシュまでの要素数(メッシュ内の要素数は含まない) real_Total_Element_to_mesh[MAX_N_MESH + 1]
		info->Total_Load_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));             // 各メッシュまでの集中荷重の数
		info->Total_Constraint_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));       // 各メッシュまでの拘束条件の数
		info->Total_DistributeForce_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));  // 各メッシュまでの分布荷重の数
		if (info->Total_Knot_to_mesh == NULL || info->Total_Patch_on_mesh == NULL || info->Total_Patch_to_mesh == NULL || info->Total_Control_Point_on_mesh == NULL || info->Total_Control_Point_to_mesh == NULL || info->Total_Element_on_mesh == NULL || info->Total_Element_to_mesh == NULL || info->real_Total_Element_on_mesh == NULL || info->real_Total_Element_to_mesh == NULL || info->Total_Load_to_mesh == NULL || info->Total_Constraint_to_mesh == NULL || info->Total_DistributeForce_to_mesh == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		info->gauss_w = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                         // 
		info->gauss_point = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * MAX_DIMENSION);     // 
		info->gauss_w_ex = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                      // 
		info->gauss_point_ex = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * MAX_DIMENSION);  // 
		info->gauss_w_1D = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                      // 
		info->gauss_point_1D = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                  // 
		info->gauss_w_1D_ex = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                   // 
		info->gauss_point_1D_ex = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);               // 
		info->Gxi_1D = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                          // 
		info->w_1D = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                            // 
		info->Gxi = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * MAX_DIMENSION);             // 
		info->w = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND);                               // 
		if (info->gauss_w == NULL || info->gauss_point == NULL || info->gauss_w_ex == NULL || info->gauss_point_ex == NULL || info->gauss_w_1D == NULL || info->gauss_point_1D == NULL || info->gauss_w_1D_ex == NULL || info->gauss_point_1D_ex == NULL || info->Gxi_1D == NULL || info->w_1D == NULL || info->Gxi == NULL || info->w == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 1)
	{
		info->Order = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));								 // 各パッチの各方向の次数を格納する配列
		info->No_knot = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));							     // 各パッチとその各方向におけるノットベクトルの数を格納する配列
		info->Total_Control_Point_to_patch = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] + 1), sizeof(int));					     // 各パッチまでの総コントロールポイント数を格納する配列
		info->Total_Knot_to_patch_dim = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION + 1), sizeof(int));		     // 各パッチと方向におけるノットの累積数を格納する配列
		info->Position_Knots = (double *)malloc(sizeof(double) * info->Total_Knot_to_mesh[Total_mesh]);									     // 各ノットの位置を格納する配列
		info->No_Control_point = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));					 // 各パッチの各方向のコントロールポイント数を格納する配列  
		info->No_Control_point_in_patch = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));							     // 各パッチのコントロールポイントの総数を格納する配列       
		info->Patch_Control_point = (int *)malloc(sizeof(int) * MAX_CP);																     // 各パッチ内のコントロールポイントのコネクティビティを保持する配列   
		info->No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));						     // 各要素のコントロールポイント数を格納する配列        
		info->Node_Coordinate = (double *)malloc(sizeof(double) * (info->Total_Control_Point_to_mesh[Total_mesh] * (info->DIMENSION + 1)));  // 各ノードの座標を格納する配列(+1は重みの分)        
		info->Control_Coord_x = (double *)malloc(sizeof(double) * MAX_CP);																     // 各制御点のx座標を格納する配列        
		info->Control_Coord_y = (double *)malloc(sizeof(double) * MAX_CP);																     // 各制御点のy座標を格納する配列          
		info->Control_Coord_z = (double *)malloc(sizeof(double) * MAX_CP);																     // 各制御点のz座標を格納する配列             
		info->Control_Weight = (double *)malloc(sizeof(double) * MAX_CP);																     // 各制御点の重みを格納する配列          
		info->Constraint_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Constraint_to_mesh[Total_mesh] * 2));						     // 制約が適用されるノードとその方向を格納する配列
		info->Value_of_Constraint = (double *)calloc(info->Total_Constraint_to_mesh[Total_mesh], sizeof(double));						     // 制約の値を格納する配列
		// 荷重に関連する情報を格納するための配列群
		info->Load_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Load_to_mesh[Total_mesh] * 2));									  // 集中荷重をかける制御点のインデックスとその方向を規定する配列
		info->Value_of_Load = (double *)calloc(info->Total_Load_to_mesh[Total_mesh], sizeof(double));									  // 集中荷重の荷重量
		info->iPatch_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // 分布荷重をかけるパッチ番号を規定する配列                
		info->iCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // 分布荷重をかける面を規定する配列                                
		info->jCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  //                          
		info->type_load_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);							  // 分布荷重の方向を規定する配列
		info->val_Coord_array = (double *)calloc(info->Total_DistributeForce_to_mesh[Total_mesh], sizeof(double));						  // 分布荷重面に対して垂直方向の座標軸を規定する配列
		info->Range_Coord_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 2 * 2), sizeof(double));			  // 分布荷重を与え始めと終わりを規定する配列
		info->Coeff_Dist_Load_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 3 * 2), sizeof(double));		  // 二次関数を用いた分布荷重を規定する配列
		if (info->Order == NULL || info->No_knot == NULL || info->Total_Control_Point_to_patch == NULL || info->Total_Knot_to_patch_dim == NULL || info->Position_Knots == NULL || info->No_Control_point == NULL || info->No_Control_point_in_patch == NULL || info->Patch_Control_point == NULL || info->No_Control_point_ON_ELEMENT == NULL || info->Node_Coordinate == NULL || info->Control_Coord_x == NULL || info->Control_Coord_y == NULL || info->Control_Weight == NULL || info->Constraint_Node_Dir == NULL || info->Value_of_Constraint == NULL || info->Load_Node_Dir == NULL || info->Value_of_Load == NULL || info->iPatch_array == NULL || info->iCoord_array == NULL || info->jCoord_array == NULL || info->type_load_array == NULL || info->val_Coord_array == NULL || info->Range_Coord_array == NULL || info->Coeff_Dist_Load_array == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 2)
	{
		info->INC = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION));		     // 各パッチにおける各制御点のインデックス振り
		info->Controlpoint_of_Element = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * MAX_NO_CP_ON_ELEMENT));							     // 各要素における各制御点番号を格納
		info->Element_patch = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																     // 要素がどのパッチ内にあるかを示す配列
		info->Element_mesh = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																     // 要素がどのメッシュ内にあるかを示す配列
		info->line_No_real_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // ゼロエレメントではない要素列の数
		info->line_No_Total_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // ゼロエレメントを含むすべての要素列の数
		info->difference = (double *)calloc(info->Total_Knot_to_mesh[Total_mesh], sizeof(double));																     // 隣り合うノットベクトルの差
		info->Total_element_all_ID = (int *)calloc(info->Total_Element_to_mesh[Total_mesh], sizeof(int));														     // ゼロエレメントではない要素 = 1, ゼロエレメント = 0
		info->ENC = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION));														 // ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
		info->real_element_line = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION));  // ゼロエレメントではない要素列
		info->real_element = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																     // ゼロエレメントではない要素の番号
		info->real_El_No_on_mesh = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh]));				     // メッシュ上の実要素番号
		info->Equivalent_Nodal_Force = (double *)calloc(MAX_CP * info->DIMENSION, sizeof(double));															         // 等価節点荷重
		if (info->INC == NULL || info->Controlpoint_of_Element == NULL || info->Element_patch == NULL || info->Element_mesh == NULL || info->line_No_real_element == NULL || info->line_No_Total_element == NULL || info->difference == NULL || info->Total_element_all_ID == NULL || info->ENC == NULL || info->real_element_line == NULL || info->real_element == NULL || info->real_El_No_on_mesh == NULL || info->Equivalent_Nodal_Force == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 3)
	{
		info->a_matrix = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * info->DIMENSION * info->DIMENSION);                                    //      
		info->dSF = (double *)malloc(sizeof(double) * MAX_POW_NG * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);                                           // 基底関数の導関数        
		info->Gauss_Coordinate = (double *)calloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG * info->DIMENSION, sizeof(double));       // ガウス点の物理座標   
		info->Loc_parameter_on_Glo = (double *)malloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG * info->DIMENSION * sizeof(double));  // ローカルメッシュ上の要素のガウス点に対応するグローバルパッチ上のパラメトリック 
		info->opp_patch_num = (int *)malloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG * sizeof(int));                                 // ローカルメッシュ上の要素のガウス点に対応するグローバルパッチ番号を格納する配列 
		if (info->a_matrix == NULL || info->dSF == NULL || info->Gauss_Coordinate == NULL || info->Loc_parameter_on_Glo == NULL || info->opp_patch_num == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}

		info->Jac = (double *)malloc(sizeof(double) * info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG);                                       // ヤコビアン 
		info->B_Matrix = (double *)malloc(sizeof(double) * info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG * D_MATRIX_SIZE * MAX_KIEL_SIZE);  // 変位ひずみマトリクス 
		if (info->Jac == NULL || info->B_Matrix == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		if (MODE_EX == 0)
		{
			info->Jac_ex = (double *)malloc(sizeof(double) * info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG_EXTEND);                                       // 拡張されたヤコビアンを格納する配列
			info->B_Matrix_ex = (double *)malloc(sizeof(double) * info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG_EXTEND * D_MATRIX_SIZE * MAX_KIEL_SIZE);  // 拡張された変位ひずみマトリクスを格納する配列
			info->dSF_ex = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);                                             // 拡張された基底関数の導関数を格納する配列
			info->Gauss_Coordinate_ex = (double *)calloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG_EXTEND * info->DIMENSION, sizeof(double));         // 拡張されたガウス点の物理座標を格納する配列
			info->Loc_parameter_on_Glo_ex = (double *)malloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG_EXTEND * info->DIMENSION * sizeof(double));    // 
			info->opp_patch_num_ex = (int *)malloc(info->real_Total_Element_to_mesh[Total_mesh] * MAX_POW_NG_EXTEND * sizeof(int));                                   // 
			if (info->Jac_ex == NULL || info->B_Matrix_ex == NULL || info->dSF_ex == NULL || info->Gauss_Coordinate_ex == NULL || info->Loc_parameter_on_Glo_ex == NULL || info->opp_patch_num_ex == NULL)
			{
				printf("Cannot allocate memory\n");
				exit(1);
			}
		info->eoi.resize(info->Total_Element_to_mesh[Total_mesh]);             // 各要素と重なっている要素情報
		info->octree_subcell.resize(info->Total_Element_to_mesh[Total_mesh]);  // 
		}
	}

	else if (num == 4)
	{
		info->D = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * D_MATRIX_SIZE);                                          // 弾性マトリクス
		info->Index_Dof = (int *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));													     // Index_Dof[MAX_K_WHOLE_SIZE];
		info->K_Whole_Ptr = (int *)calloc(MAX_K_WHOLE_SIZE + 1, sizeof(int));											     // 全体剛性マトリクスの行ポインタ
		if (info->D == NULL || info->Index_Dof == NULL || info->K_Whole_Ptr == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 5)
	{
		info->Displacement = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE);     // Displacement[MAX_K_WHOLE_SIZE]
		if (info->Displacement == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
}


// set global variable
void Global_var(const int num, information *info)
{
	if (num == 0)
	{
		MAX_ORDER += 1;                                              // 次数の最大値 + 1  (要素内の一方向の制御点の数)
		MAX_NO_CP_ON_ELEMENT = pow_int(MAX_ORDER, info->DIMENSION);  // 要素内の制御点の数
		MAX_KIEL_SIZE = MAX_NO_CP_ON_ELEMENT * info->DIMENSION;      // 要素剛性マトリクスサイズ
	}
	else if (num == 1)
	{
		if (info->DIMENSION == 2)
		{
			D_MATRIX_SIZE = 3;

			// 0: xx, 1: yy, 2: xy, 3: zz
			N_STRAIN = 4;
			N_STRESS = 4;
		}
		else if (info->DIMENSION == 3)
		{
			D_MATRIX_SIZE = 6;

			// 0: xx, 1: yy, 2: zz, 3: xy, 4: yz, 5: xz
			N_STRAIN = 6;
			N_STRESS = 6;
		}
	}
	else if (num == 2)
	{
		MAX_K_WHOLE_SIZE = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;  // 全体剛性マトリクスのサイズ
	}
}


// Read file 1st time
void Get_Input_1(int tm, char **argv, information *info)
{
	char s[256];
	int temp_i;

	int i, j;

	fp = fopen(argv[tm + 1], "r");

	// DIMENSION
	fscanf(fp, "%d", &temp_i);
	info->DIMENSION = temp_i;
	fgets(s, 256, fp);
	printf("DIMENSION: %d\n", info->DIMENSION);
	if (info->DIMENSION != 2 && info->DIMENSION != 3)
	{
		printf("Error, wrong DIMENSION in input data\n");
		exit(1);
	}

	// 材料定数
	fscanf(fp, "%lf %lf", &E, &nu);
	fgets(s, 256, fp);
	printf("E: %le, nu: %le\n", E, nu);

	// パッチ数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;
	int *CP = (int *)malloc(sizeof(int) * temp_i * info->DIMENSION);
	printf("No_Patch: %d\n", temp_i);
	info->Total_Patch_on_mesh[tm] = temp_i;
	info->Total_Patch_to_mesh[tm + 1] = info->Total_Patch_to_mesh[tm] + temp_i;
	printf("Total_Patch_to_mesh[%d] = %d\n", tm, info->Total_Patch_to_mesh[tm]);

	// コントロールポイント数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	printf("Total_Control_Point: %d\n", temp_i);
	info->Total_Control_Point_on_mesh[tm] = temp_i;
	info->Total_Control_Point_to_mesh[tm + 1] = info->Total_Control_Point_to_mesh[tm] + temp_i;
	printf("Total_Control_Point_to_mesh[%d] = %d\n", tm + 1, info->Total_Control_Point_to_mesh[tm + 1]);

	// 各方向の次数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	info->Total_Knot_to_mesh[tm + 1] = info->Total_Knot_to_mesh[tm];
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Total_Knot_to_mesh[tm + 1] += temp_i;
		}
	}

	// 各パッチ各方向のコントロールポイント数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			CP[i * info->DIMENSION + j] = temp_i;
		}

		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		MAX_CP += temp_MAX_CP;
	}
	fgets(s, 256, fp);

	// パッチコネクティビティ(スキップ)
	for (i = 0; i < info->Total_Patch_on_mesh[tm]; i++)
	{
		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		for (j = 0; j < temp_MAX_CP; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// 境界条件
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	info->Total_Constraint_to_mesh[tm + 1] = info->Total_Constraint_to_mesh[tm] + Total_Constraint;
	info->Total_Load_to_mesh[tm + 1] = info->Total_Load_to_mesh[tm] + Total_Load;
	info->Total_DistributeForce_to_mesh[tm + 1] = info->Total_DistributeForce_to_mesh[tm] + Total_DistributeForce;

	printf("Total_Constraint = %d\n", Total_Constraint);
	printf("Total_Load = %d\n", Total_Load);
	printf("Total_DistributedForce = %d\n", Total_DistributeForce);

	printf("\n");

	fclose(fp);
	free(CP);
}


// Read file 2nd time
void Get_Input_2(int tm, char **argv, information *info)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j, k;

	fp = fopen(argv[tm + 1], "r");

	// DIMENSION(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);

	// 材料定数(スキップ)
	fscanf(fp, "%lf %lf", &temp_d, &temp_d);
	fgets(s, 256, fp);

	// パッチ数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;

	// コントロールポイント数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int Total_Control_Point = temp_i;

	// 各方向の次数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			if (MAX_ORDER < temp_i)
			{
				MAX_ORDER = temp_i;
			}
			// printf("Order[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			if (MAX_KNOT < temp_i)
			{
				MAX_KNOT = temp_i;
			}
			info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j + 1] = info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + temp_i;
			// printf("No_knot[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// 各パッチ各方向のコントロールポイント数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			// printf("No_Control_point[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < No_Patch; i++)
	{
		info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] = 1;
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] *= info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j];
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm] + 1] = info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]];
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			if (info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] != info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n", tm, i, j);
			}
		}
	}

	// for debug
	// for (i = 0; i < No_Patch; i++)
	// {
	// 	printf("No_Control_point_in_patch[%d] = %d\n", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]);
	// }

	printf("\n");

	// パッチコネクティビティ
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] = temp_i;
			if (tm > 0)
			{
				info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] += info->Total_Control_Point_to_mesh[tm];
			}
		}
	}
	fgets(s, 256, fp);

	// 境界条件(スキップ)
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	fgets(s, 256, fp);

	// ノットベクトルの読み込み
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]; k++)
			{
				fscanf(fp, "%lf", &temp_d);
				info->Position_Knots[info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + k] = temp_d;
			}
		}
	}
	fgets(s, 256, fp);

	int Total_Element = 0;
	for (i = 0; i < No_Patch; i++)
	{
		if (info->DIMENSION == 2)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1);
		}
		else if (info->DIMENSION == 3)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] + 1);
		}
	}
	printf("Total_Element = %d\n", Total_Element);
	info->Total_Element_on_mesh[tm] = Total_Element;
	info->Total_Element_to_mesh[tm + 1] = info->Total_Element_to_mesh[tm] + Total_Element;
	printf("Total_Element_on_mesh[%d] = %d\n", tm, info->Total_Element_on_mesh[tm]);

	// for debug
	// for (i = 0; i < No_Patch; i++)
	// {
	// 	printf("No_Control_point_ON_ELEMENT[%d] = %d\n", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]]);
	// }

	// 節点座標
	for (i = 0; i < Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &temp_i);
		for (j = 0; j < info->DIMENSION + 1; j++)
		{
			fscanf(fp, "%lf", &temp_d);
			info->Node_Coordinate[(temp_i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j] = temp_d;
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		int cp_to_patch = info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]];
		for (j = 0; j < info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]; j++)
		{
			int connectivity = info->Patch_Control_point[cp_to_patch + j] * (info->DIMENSION + 1);
			if (info->DIMENSION == 2)
			{
				info->Control_Coord_x[cp_to_patch + j] = info->Node_Coordinate[connectivity];
				info->Control_Coord_y[cp_to_patch + j] = info->Node_Coordinate[connectivity + 1];
				info->Control_Weight[cp_to_patch + j] = info->Node_Coordinate[connectivity + 2];
			}
			else if (info->DIMENSION == 3)
			{
				info->Control_Coord_x[cp_to_patch + j] = info->Node_Coordinate[connectivity];
				info->Control_Coord_y[cp_to_patch + j] = info->Node_Coordinate[connectivity + 1];
				info->Control_Coord_z[cp_to_patch + j] = info->Node_Coordinate[connectivity + 2];
				info->Control_Weight[cp_to_patch + j] = info->Node_Coordinate[connectivity + 3];
			}
		}
	}
	fgets(s, 256, fp);

	// 拘束
	for (i = 0; i < Total_Constraint; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0],
			   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1],
			   &info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]]);
		info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] = info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];
		// cout << "Constraint_Node_Dir[" << (i + info->Total_Constraint_to_mesh[tm]) * 2 + 0 << "] = " << info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] << " Constraint_Node_Dir[" << (i + info->Total_Constraint_to_mesh[tm]) * 2 + 1 << "] = " << info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1] << " Value_of_Constraint[" << i + info->Total_Constraint_to_mesh[tm] << "] = " << info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]] << endl;
	}
	fgets(s, 256, fp);

	// 荷重
	for (i = 0; i < Total_Load; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
			   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
			   &info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
		info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] = info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];

		printf("Load_Node_Dir[%d]= %d Load_Node_Dir[%d]= %d Value_of_Load[%d] = %le\n",
			   (i + info->Total_Load_to_mesh[tm]) * 2 + 0, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
			   (i + info->Total_Load_to_mesh[tm]) * 2 + 1, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
			   i + info->Total_Load_to_mesh[tm], info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	if (info->DIMENSION == 2)
	{
		int type_load, iPatch, iCoord;
		double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];

		for (i = 0; i < Total_DistributeForce; i++)
		{
			fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
			printf("Distibuted load number: %d\n", i);
			printf("type_load: %d iPatch: %d iCoord: %d val_Coord: %le Range_Coord: %le %le\nCoef_Dist_Load: %le %le %le\n",
					type_load, iPatch, iCoord, val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);

			// for S-IGA
			info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
			info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
			info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
			info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
		}
	}
	else if (info->DIMENSION == 3)
	{
		int type_load, iPatch, iCoord, jCoord;
		double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];

		for (i = 0; i < Total_DistributeForce; i++)
		{
			fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &jCoord, &val_Coord, &iRange_Coord[0], &iRange_Coord[1], &jRange_Coord[0], &jRange_Coord[1], &iCoeff_Dist_Load[0], &iCoeff_Dist_Load[1], &iCoeff_Dist_Load[2], &jCoeff_Dist_Load[0], &jCoeff_Dist_Load[1], &jCoeff_Dist_Load[2]);
			printf("Distibuted load number: %d\n", i);
			printf("type_load: %d iPatch: %d iCoord: %d jCoord: %d val_Coord: %le\nRange_Coord: %le %le %le %le\nCoef_Dist_Load: %le %le %le %le %le %le\n",
					type_load, iPatch, iCoord, jCoord, val_Coord,
					iRange_Coord[0], iRange_Coord[1], jRange_Coord[0], jRange_Coord[1],
					iCoeff_Dist_Load[0], iCoeff_Dist_Load[1], iCoeff_Dist_Load[2], jCoeff_Dist_Load[0], jCoeff_Dist_Load[1], jCoeff_Dist_Load[2]);

			// for S-IGA
			info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
			info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
			info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
			info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = jCoord;
			info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0] = iRange_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1] = iRange_Coord[1];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2] = jRange_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3] = jRange_Coord[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0] = iCoeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1] = iCoeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2] = iCoeff_Dist_Load[2];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3] = jCoeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4] = jCoeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5] = jCoeff_Dist_Load[2];
		}
	}
	fclose(fp);
}


// INC 等の作成
void Make_INC(information *info)
{
	int tm;

	// info->INC の計算(節点番号をξ, ηの番号で表す為の配列)
	for (tm = 0; tm < Total_mesh; tm++)
	{
		int b, B, e, h, i, j, k, l, n, o, p, q, x, y, z, ii, jj, kk, kkk, iiloc, jjloc, kkloc, r = 0;
		int type_load, iPatch, iCoord, jCoord;
		double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];
		int No_Patch = info->Total_Patch_on_mesh[tm];
		int Total_Patch_to_Now = info->Total_Patch_to_mesh[tm];
		int Total_Element = info->Total_Element_on_mesh[tm];
		int Total_Element_to_Now = info->Total_Element_to_mesh[tm];
		int Total_DistributeForce = info->Total_DistributeForce_to_mesh[tm + 1] - info->Total_DistributeForce_to_mesh[tm];

		if (info->DIMENSION == 2)
		{
			e = 0;  // 要素数をカウントするための変数
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;  // パッチ内の制御点の数をカウントするための変数
				int patch = l + Total_Patch_to_Now;  // 現在のパッチのグローバルインデックス
				for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)  // 現在のパッチのη方向の制御点数だけ処理を行う
				{
					for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)  // 現在のパッチのξ方向の制御点数だけ処理を行う
					{
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;   // INC : 各パッチの御点の各方向のインデックス振り 

						if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1])  // 現在の制御点が要素の範囲内にあるか
						{
							for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)   // jjloc : η方向の制御点のローカルインデックス分だけ処理を行う
							{
								for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)  // iiloc : ξ方向の制御点のローカルインデックス分だけ処理を行う
								{
									B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
									b = jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
									info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
								}
							}
							info->Element_patch[e + Total_Element_to_Now] = patch;
							info->Element_mesh[e + Total_Element_to_Now] = tm;
							e++;
						}
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (kk = 0; kk < info->No_Control_point[patch * info->DIMENSION + 2]; kk++)
				{
					for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
					{
						for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
						{
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 2] = kk;

							if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1] && kk >= info->Order[patch * info->DIMENSION + 2])
							{
								for (kkloc = 0; kkloc <= info->Order[patch * info->DIMENSION + 2]; kkloc++)
								{
									for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
									{
										for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
										{
											B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - kkloc * info->No_Control_point[patch * info->DIMENSION + 0] * info->No_Control_point[patch * info->DIMENSION + 1] - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
											b = kkloc * (info->Order[patch * info->DIMENSION + 0] + 1) * (info->Order[patch * info->DIMENSION + 1] + 1) + jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
											info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
										}
									}
								}
								info->Element_patch[e + Total_Element_to_Now] = patch;
								info->Element_mesh[e + Total_Element_to_Now] = tm;
								e++;
							}
							i++;
						}
					}
				}
			}
		}

		// for S-IGA line_No_real_elementの初期化
		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_real_element[patch * info->DIMENSION + j] = 0;
			}
		}

		//  line_No_real_elementの計算
		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_Total_element[patch * info->DIMENSION + j] = info->No_knot[patch * info->DIMENSION + j] - 2 * info->Order[patch * info->DIMENSION + j] - 1;

				for (kkk = info->Order[patch * info->DIMENSION + j]; kkk < info->No_knot[patch * info->DIMENSION + j] - info->Order[patch * info->DIMENSION + j] - 1; kkk++)
				{
					info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]]
						= info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk];
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]] != 0)
					{
						info->line_No_real_element[patch * info->DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		for (h = 0; h < Total_Element; h++)
		{
			info->Total_element_all_ID[h] = 0;
		}

		if (info->DIMENSION == 2)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
				{
					for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
					{
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (z = 0; z < info->line_No_Total_element[patch * info->DIMENSION + 2]; z++)
				{
					for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
					{
						for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
						{
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 2] = z;
							i++;
						}
					}
				}
			}
		}

		// 必要な要素の行と列の番号を求める
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (l = 0; l < No_Patch; l++)
			{
				e = 0;
				int patch = l + Total_Patch_to_Now;
				for (k = 0; k < info->line_No_Total_element[patch * info->DIMENSION + j]; k++)
				{
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + k] != 0)
					{
						info->real_element_line[patch * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + e * info->DIMENSION + j] = k;
						e++;
					}
				}
			}
		}

		// 必要な要素列上の要素のIDを1にする
		if (info->DIMENSION == 2)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; p++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 0])
					{
						for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; q++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 1])
							{
								info->Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = ele;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}
		else if (info->DIMENSION == 3)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (o = 0; o < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; o++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + o * info->DIMENSION + 0])
					{
						for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; p++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 1])
							{
								for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 2]; q++)
								{
									if (info->ENC[ele * info->DIMENSION + 2] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 2])
									{
										info->Total_element_all_ID[n]++;
									}
								}
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = ele;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1] * info->line_No_real_element[patch * info->DIMENSION + 2];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}

		// 分布荷重
		for (i = 0; i < Total_DistributeForce; i++)
		{
			if (info->DIMENSION == 2)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2];

				Setting_Dist_Load_2D(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
			}
			else if (info->DIMENSION == 3)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				jCoord = info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1];
				jRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2];
				jRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2];
				jCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3];
				jCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4];
				jCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5];

				Setting_Dist_Load_3D(tm, iPatch, iCoord, jCoord, val_Coord, iRange_Coord, jRange_Coord, type_load, iCoeff_Dist_Load, jCoeff_Dist_Load, info);
			}
		}
	}
}


// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info)
{
	int iii, jjj;
	int N_Seg_Load_Element_iDir = 0, jCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load;
	int iX, iY;
	int ic, ig;
	double val_jCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	static int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	static int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(0, info);

	// icoord : 分布荷重をかける面の情報
	// iCoord=0: Load on Eta=Constant
	// iCoord=1: Load on Xi=Constant
	if (iCoord == 0)
	{
		jCoord = 1;
	}
	else if (iCoord == 1)
	{
		jCoord = 0;
	}

	// 分布荷重の開始点と終点があるノットスパンを調べる
	for (iii = info->Order[iPatch * info->DIMENSION + iCoord]; iii < info->No_knot[iPatch * info->DIMENSION + iCoord] - info->Order[iPatch * info->DIMENSION + iCoord] - 1; iii++)   // ノットスパンを持ちうる範囲だけ処理
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii] - epsi <= Range_Coord[0])
			iPos[0] = iii;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii + 1] - epsi <= Range_Coord[1])
			iPos[1] = iii + 1;
	}

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(1);
	}

	// 分布荷重面に垂直な方向の荷重点があるノットスパンの検索
	for (jjj = info->Order[iPatch * info->DIMENSION + jCoord]; jjj < info->No_knot[iPatch * info->DIMENSION + jCoord] - info->Order[iPatch * info->DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj] - epsi <= val_Coord
			&& info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] + epsi > val_Coord)
		{
			jPos[0] = jjj;
			jPos[1] = jjj + 1;
			val_jCoord_Local = -1.0 + 2.0 * (val_Coord - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj])
							 / (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj]);
		}
	}

	//
	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(1);
	}

	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		N_Seg_Load_Element_iDir++;
	}

	iii = 0;
	if (iCoord == 1)  
	{
		iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0];
		for (iY = iPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < iPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			No_Element_for_Integration[iii] = SearchForElement_2D(mesh_n, iPatch, iX, iY, info);
			iii++;
		}
	}

	if (iCoord == 0)
	{
		iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1];
		for (iX = iPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < iPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			No_Element_for_Integration[iii] = SearchForElement_2D(mesh_n, iPatch, iX, iY, info);
			iii++;
		}
	}
	No_Element_For_Dist_Load = iii;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load; iii++)
	{
		if (info->Total_element_all_ID[No_Element_for_Integration[iii]] == 1)
		{
			iX = info->ENC[No_Element_for_Integration[iii] * info->DIMENSION + 0];
			iY = info->ENC[No_Element_for_Integration[iii] * info->DIMENSION + 1];

			for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
				iControlpoint[ic] = info->Controlpoint_of_Element[No_Element_for_Integration[iii] * MAX_NO_CP_ON_ELEMENT + ic];

			for (ig = 0; ig < GP_1D; ig++)
			{
				double Local_Coord[2], sfc, dxyzdge[3] = {0.0}, detJ, XiEtaCoordParen, valDistLoad;
				int icc;
				Local_Coord[jCoord] = val_jCoord_Local;
				Local_Coord[iCoord] = info->Gxi_1D[ig];

				ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii], info);
				XiEtaCoordParen = Position_Data_param[iCoord];
				valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen + Coeff_Dist_Load[2] * XiEtaCoordParen * XiEtaCoordParen;

				for (icc = 0; icc < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); icc++)
				{
					double temp_i = dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], info);

					dxyzdge[0] += temp_i * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
					dxyzdge[1] += temp_i * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
				}

				detJ = sqrt(dxyzdge[0] * dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);

				if (type_load < 2)
				{
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * info->w_1D[ig];
					}
				}
				else if (type_load == 2) // 法線方向
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1] / detJ;
					LoadDir[1] = -dxyzdge[0] / detJ;
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * info->w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * info->w_1D[ig];
					}
				}
				else if (type_load == 3)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[0] / detJ;
					LoadDir[1] = dxyzdge[1] / detJ;
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * info->w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * info->w_1D[ig];
					}
				}
			}
		}
	}
}


void Setting_Dist_Load_3D(int mesh_n, int iPatch, int iCoord, int jCoord, double val_Coord, double *iRange_Coord, double *jRange_Coord, int type_load, double *iCoeff_Dist_Load, double *jCoeff_Dist_Load, information *info)
{
	int iii, jjj, kkk;
	int N_Seg_Load_Element_iDir = 0, N_Seg_Load_Element_jDir = 0, kCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000}, kPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load_iDir, No_Element_For_Dist_Load_jDir;
	int iX, iY, iZ;
	int ic, ig_i, ig_j;
	double val_kCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	static int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh] * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	static int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(0, info);

	if (iCoord == 0 && jCoord == 1)
	{
		kCoord = 2;
	}
	else if (iCoord == 1 && jCoord == 2)
	{
		kCoord = 0;
	}
	else if (iCoord == 2 && jCoord == 0)
	{
		kCoord = 1;
	}
	else
	{
		printf("Error, incorrect input data at distributed load.\n");
		exit(1);
	}

	for (iii = info->Order[iPatch * info->DIMENSION + iCoord]; iii < info->No_knot[iPatch * info->DIMENSION + iCoord] - info->Order[iPatch * info->DIMENSION + iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii] - epsi <= iRange_Coord[0])
			iPos[0] = iii;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii + 1] - epsi <= iRange_Coord[1])
			iPos[1] = iii + 1;
	}
	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(1);
	}

	for (jjj = info->Order[iPatch * info->DIMENSION + jCoord]; jjj < info->No_knot[iPatch * info->DIMENSION + jCoord] - info->Order[iPatch * info->DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj] - epsi <= jRange_Coord[0])
			jPos[0] = jjj;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] - epsi <= jRange_Coord[1])
			jPos[1] = jjj + 1;
	}
	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(1);
	}

	for (kkk = info->Order[iPatch * info->DIMENSION + kCoord]; kkk < info->No_knot[iPatch * info->DIMENSION + kCoord] - info->Order[iPatch * info->DIMENSION + kCoord] - 1; kkk++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk] - epsi <= val_Coord
			&& info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk + 1] + epsi > val_Coord)
		{
			kPos[0] = kkk;
			kPos[1] = kkk + 1;
			val_kCoord_Local = - 1.0 + 2.0 * (val_Coord - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk])
							 / (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk]);
		}
	}
	if (kPos[0] < 0 || kPos[1] < 0)
	{
		printf("Error (Stop) kPos[0] = %d kPos[1] = %d\n", kPos[0], kPos[1]);
		exit(1);
	}

	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		N_Seg_Load_Element_iDir++;
	}

	for (jjj = jPos[0]; jjj < jPos[1]; jjj++)
	{
		N_Seg_Load_Element_jDir++;
	}

	iii = 0, jjj = 0;

	if (iCoord == 0 && jCoord == 1)
	{
		iZ = kPos[0] - info->Order[iPatch * info->DIMENSION + 2];
		for (iX = iPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < iPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			for (iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < jPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < jPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			jjj++;
		}
	}
	else if (iCoord == 1 && jCoord == 2)
	{
		iX = kPos[0] - info->Order[iPatch * info->DIMENSION + 0];
		for (iY = iPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < iPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			for (iZ = jPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < jPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iZ = jPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < jPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
		{
			jjj++;
		}
	}
	else if (iCoord == 2 && jCoord == 0)
	{
		iY = kPos[0] - info->Order[iPatch * info->DIMENSION + 1];
		for (iZ = iPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < iPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
		{
			for (iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < jPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < jPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			jjj++;
		}
	}

	No_Element_For_Dist_Load_iDir = iii;
	No_Element_For_Dist_Load_jDir = jjj;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load_iDir; iii++)
	{
		for (jjj = 0; jjj < No_Element_For_Dist_Load_jDir; jjj++)
		{
			if (info->Total_element_all_ID[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj]] == 1)
			{
				iX = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 0];
				iY = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 1];
				iZ = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 2];

				for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
					iControlpoint[ic] = info->Controlpoint_of_Element[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * MAX_NO_CP_ON_ELEMENT + ic];

				for (ig_i = 0; ig_i < GP_1D; ig_i++)
				{
					for (ig_j = 0; ig_j < GP_1D; ig_j++)
					{
						double Local_Coord[3], sfc, dxyzdgez_i[3] = {0.0}, dxyzdgez_j[3] = {0.0}, detJ, CoordParen_iDir, CoordParen_jDir, valDistLoad;
						int icc;
						Local_Coord[kCoord] = val_kCoord_Local;
						Local_Coord[iCoord] = info->Gxi_1D[ig_i];
						Local_Coord[jCoord] = info->Gxi_1D[ig_j];

						ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
						CoordParen_iDir = Position_Data_param[iCoord];

						ShapeFunc_from_paren(Position_Data_param, Local_Coord, jCoord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
						CoordParen_jDir = Position_Data_param[jCoord];

						valDistLoad = (iCoeff_Dist_Load[0] + iCoeff_Dist_Load[1] * CoordParen_iDir + iCoeff_Dist_Load[2] * CoordParen_iDir * CoordParen_iDir)
									* (jCoeff_Dist_Load[0] + jCoeff_Dist_Load[1] * CoordParen_jDir + jCoeff_Dist_Load[2] * CoordParen_jDir * CoordParen_jDir);

						for (icc = 0; icc < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); icc++)
						{
							double temp_i = dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
							double temp_j = dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);

							dxyzdgez_i[0] += temp_i * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
							dxyzdgez_i[1] += temp_i * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
							dxyzdgez_i[2] += temp_i * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 2];

							dxyzdgez_j[0] += temp_j * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
							dxyzdgez_j[1] += temp_j * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
							dxyzdgez_j[2] += temp_j * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 2];
						}

						detJ = (dxyzdgez_i[0] * dxyzdgez_j[1] - dxyzdgez_i[1] * dxyzdgez_j[0]) + (dxyzdgez_i[1] * dxyzdgez_j[2] - dxyzdgez_i[2] * dxyzdgez_j[1]) + (dxyzdgez_i[2] * dxyzdgez_j[0] - dxyzdgez_i[0] * dxyzdgez_j[2]);

						if (type_load < 3)
						{
							for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
							{
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * info->w_1D[ig_i] * info->w_1D[ig_j];
							}
						}
						else if (type_load == 3) // 法線方向
						{
							double LoadDir[3];
							LoadDir[0] = (dxyzdgez_i[1] * dxyzdgez_j[2] + dxyzdgez_i[2] * dxyzdgez_j[1]) / detJ;
							LoadDir[1] = (dxyzdgez_i[2] * dxyzdgez_j[0] + dxyzdgez_i[0] * dxyzdgez_j[2]) / detJ;
							LoadDir[2] = (dxyzdgez_i[0] * dxyzdgez_j[1] + dxyzdgez_i[1] * dxyzdgez_j[0]) / detJ;
							for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
							{
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] +=   LoadDir[0] * valDistLoad * sfc * detJ * info->w_1D[ig_i] * info->w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] +=   LoadDir[1] * valDistLoad * sfc * detJ * info->w_1D[ig_i] * info->w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 2] += - LoadDir[2] * valDistLoad * sfc * detJ * info->w_1D[ig_i] * info->w_1D[ig_j];
							}
						}
					}
				}
			}
		}
	}
}


int SearchForElement_2D(int mesh_n, int iPatch, int iX, int iY, information *info)
{
	int iii;

	for (iii = 0; iii < info->Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (info->Element_patch[iii + info->Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 0] && iY == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 1])
				goto loopend_2D;
		}
	}
	loopend_2D:

	return (iii);
}


int SearchForElement_3D(int mesh_n, int iPatch, int iX, int iY, int iZ, information *info)
{
	int iii;

	for (iii = 0; iii < info->Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (info->Element_patch[iii + info->Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 0] && iY == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 1] && iZ == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 2])
				goto loopend_3D;
		}
	}
	loopend_3D:

	return (iii);
}


// for IGA
void Preprocessing_IGA(information *info)
{
	int e, global_mesh_num = 0, m = 0;

	Make_gauss_array(m, info);
	#pragma omp parallel for
	for (e = 0; e < info->real_Total_Element_on_mesh[global_mesh_num]; e++)
	{
		Preprocessing(m, e, info);
	}
	#pragma omp barrier
}


// for S_IGA
void Preprocessing_S_IGA(information *info)
{
	// global patch
	Make_gauss_array(0, info);
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
	{
		Preprocessing(0, i, info);
	}
	#pragma omp barrier

	// local patch
	Make_gauss_array(0, info);
	#pragma omp parallel for
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		Preprocessing(0, i, info);

		for (int n = 0; n < GP_ON_ELEMENT; n++)
		{
			double output_para[MAX_DIMENSION] = {0.0};
			double data_result_shape[MAX_DIMENSION] = {0.0};
			for (int j = 0; j < info->DIMENSION; j++)
				data_result_shape[j] = info->Gauss_Coordinate[i * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + j];

			int itr_n = 0, p = 0;
			for (int j = 0; j < info->Total_Patch_on_mesh[0]; j++)
			{
				itr_n = calc_patch_parameter_coord(data_result_shape, j, output_para, info);
				if (itr_n != ERROR)
				{
					p = j;
					info->opp_patch_num[i * GP_ON_ELEMENT + n] = p;
					break;
				}
			}

			for (int j = 0; j < info->DIMENSION; j++)
				info->Loc_parameter_on_Glo[i * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + j] = output_para[j];

			// data_result_shapeがグローバルメッシュ上にないとき
			if (itr_n == ERROR)
			{
				printf("ERROR local patch is not inside of global patch\n");
				exit(1);
			}
		}
	}
	#pragma omp barrier

	// local patch
	if (MODE_EX == 0)
	{
		Make_gauss_array(1, info);
		#pragma omp parallel for
		for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			if (info->eoi[i].size() < 2)
				continue;

			Preprocessing(1, i, info);

			for (int n = 0; n < GP_ON_ELEMENT; n++)
			{
				double output_para[MAX_DIMENSION] = {0.0};
				double data_result_shape[MAX_DIMENSION] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					data_result_shape[j] = info->Gauss_Coordinate_ex[i * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + j];

				int itr_n = 0, p = 0;
				for (int j = 0; j < info->Total_Patch_on_mesh[0]; j++)
				{
					itr_n = calc_patch_parameter_coord(data_result_shape, j, output_para, info);
					if (itr_n != ERROR)
					{
						p = j;
						info->opp_patch_num_ex[i * GP_ON_ELEMENT + n] = p;
						break;
					}
				}

				for (int j = 0; j < info->DIMENSION; j++)
					info->Loc_parameter_on_Glo_ex[i * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + j] = output_para[j];

				// data_result_shapeがグローバルメッシュ上にないとき
				if (itr_n == ERROR)
				{
					printf("ERROR local patch is not inside of global patch\n");
					exit(1);
				}
			}
		}
		#pragma omp barrier
	}
}



// Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double *para_coord, information *info)
{
	int line[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int line_counter = 0;
		for (int j = 0; j < info->No_Control_point[patch_n * info->DIMENSION + i] - info->Order[patch_n * info->DIMENSION + i]; j++)
		{
			double knot_front = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j];
			double knot_back = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j + 1];

			if (j == 0 && fabs(knot_front - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_front;
				line[i] = line_counter;
				break;
			}

			if (fabs(knot_back - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_back;
				line[i] = line_counter;
				break;
			}

			if (knot_front < para_coord[i] && para_coord[i] < knot_back)
			{
				line[i] = line_counter;
				break;
			}

			line_counter++;
		}
	}

	int ele = 0;
	if (info->DIMENSION == 2)
		ele = line[0] + line[1] * info->line_No_Total_element[patch_n * info->DIMENSION + 0];
	else if (info->DIMENSION == 3)
		ele = line[0] + line[1] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] + line[2] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] * info->line_No_Total_element[patch_n * info->DIMENSION + 1];

	for (int i = 0; i < patch_n; i++)
	{
		int ele_in_patch = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			ele_in_patch *= (info->No_Control_point[i * info->DIMENSION + j] - info->Order[i * info->DIMENSION + j]);
		ele += ele_in_patch;
	}

	return ele;
}


void tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		xi_tilde[i] = -1.0 + 2.0 * (xi[i] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i]]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i]]);
	}
}


// 昇順ソート
void sort(int total, int *element_n_point)
{
	int i, j;
	int temp;

	for (i = 0; i < total; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (element_n_point[i] > element_n_point[j])
			{
				temp = element_n_point[i];
				element_n_point[i] = element_n_point[j];
				element_n_point[j] = temp;
			}
		}
	}
}


// Preprocessing
void Preprocessing(int m, int e, information *info)
{
	static double *a_matrix_parallel = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * info->DIMENSION * info->DIMENSION * omp_get_max_threads());
	static double *dsf_parallel = (double *)malloc(sizeof(double) * MAX_POW_NG_EXTEND * MAX_NO_CP_ON_ELEMENT * info->DIMENSION * omp_get_max_threads());
	double *a_matrix = &a_matrix_parallel[MAX_POW_NG_EXTEND * info->DIMENSION * info->DIMENSION * omp_get_thread_num()];
	double *dsf = &dsf_parallel[MAX_POW_NG_EXTEND * MAX_NO_CP_ON_ELEMENT * info->DIMENSION * omp_get_thread_num()];

	// ガウス点の物理座標を計算
	Make_Gauss_Coordinate(m, e, info);

	// ガウス点でのヤコビアン, Bマトリックスを計算
	Make_dSF(m, e, dsf, info);
	Make_Jac(m, e, a_matrix, dsf, info);
	Make_B_Matrix(m, e, a_matrix, dsf, info);
}


void Make_Gauss_Coordinate(int m, int e, information *info)
{
	int i, j, k;

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{	
		for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
		{
			double R = Shape_func(j, &info->Gxi[i * info->DIMENSION], e, info);

			for (k = 0; k < info->DIMENSION; k++)
			{
				if (m == 0)
				{
					info->Gauss_Coordinate[e * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + k] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				}
				else if (m == 1)
				{
					info->Gauss_Coordinate_ex[e * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + k] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				}
			}
		}
	}
}


void Make_dSF(int m, int e, double *dsf, information *info)
{
	int i, j, k;

	if (m == 0)
	{
		for (i = 0; i < GP_ON_ELEMENT; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + j * info->DIMENSION + k] = dShape_func(j, k, &info->Gxi[i * info->DIMENSION], e, info);
				}
			}
		}
	}
	else if (m == 1)
	{
		for (i = 0; i < GP_ON_ELEMENT; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + j * info->DIMENSION + k] = dShape_func(j, k, &info->Gxi[i * info->DIMENSION], e, info);
				}
			}
		}
	}
}


void Make_Jac(int m, int e, double *a_matrix, double *dsf, information *info)
{
	int i, j, k, l;
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->DIMENSION; k++)
			{
				if (info->DIMENSION == 2)
				{
					a_2x2[j][k] = 0.0;
					for (l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
					{
						if (m == 0)
						{
							a_2x2[j][k] += dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
						else if (m == 1)
						{
							a_2x2[j][k] += dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
					}
				}
				else if (info->DIMENSION == 3)
				{
					a_3x3[j][k] = 0.0;
					for (l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
					{
						if (m == 0)
						{
							a_3x3[j][k] += dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
						else if (m == 1)
						{
							a_3x3[j][k] += dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
					}
				}
			}
		}

		if (info->DIMENSION == 2)
		{
			J = InverseMatrix_2x2(a_2x2);

			for (j = 0; j < info->DIMENSION; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					a_matrix[i * info->DIMENSION * info->DIMENSION + j * info->DIMENSION + k] = a_2x2[j][k];
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			J = InverseMatrix_3x3(a_3x3);

			for (j = 0; j < info->DIMENSION; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					a_matrix[i * info->DIMENSION * info->DIMENSION + j * info->DIMENSION + k] = a_3x3[j][k];
				}
			}
		}

		if (m == 0)
		{
			info->Jac[e * GP_ON_ELEMENT + i] = J;
		}
		else if (m == 1)
		{
			info->Jac_ex[e * GP_ON_ELEMENT + i] = J;
		}

		if (J <= 0)
		{
			double jac;
			if (J == ERROR)
				jac = 0.0;
			else
				jac = J;

			printf("Error, J = %le\nJ must be positive value.", jac);
			exit(1);
		}
	}
}


void Make_B_Matrix(int m, int e, double *a_matrix, double *dsf, information *info)
{
	int i, j, k, l;
	static double *b_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	double *b = &b_parallel[info->DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
			{
				b[j * MAX_NO_CP_ON_ELEMENT + k] = 0.0;
				for (l = 0; l < info->DIMENSION; l++)
				{
					if (m == 0)
					{
						b[j * MAX_NO_CP_ON_ELEMENT + k] += a_matrix[i * info->DIMENSION * info->DIMENSION + l * info->DIMENSION + j] * dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + k * info->DIMENSION + l];
					}
					else if (m == 1)
					{
						b[j * MAX_NO_CP_ON_ELEMENT + k] += a_matrix[i * info->DIMENSION * info->DIMENSION + l * info->DIMENSION + j] * dsf[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + k * info->DIMENSION + l];
					}
				}
			}
		}

		if (info->DIMENSION == 2)
		{
			if (m == 0)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
					
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)]     = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
			else if (m == 1)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)]     = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			if (m == 0)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 2)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 1)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 2)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j)]     = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 2)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
			else if (m == 1)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 2)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 1)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 2)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j)]     = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[(long)e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 2)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
		}
	}
}


void Make_gauss_array(int select_GP, information *info)
{
	static int counter_GP = 0;

	static int GP_1D_NG = NG;
	static int GP_2D_NG = NG * NG;
	static int GP_3D_NG = NG * NG * NG;

	static int GP_1D_NG_EXTEND = NG_EXTEND;
	static int GP_2D_NG_EXTEND = NG_EXTEND * NG_EXTEND;
	static int GP_3D_NG_EXTEND = NG_EXTEND * NG_EXTEND * NG_EXTEND;

	if (counter_GP == 0)
	{
		Gauss_point(info, NG, info->gauss_w, info->gauss_point, info->gauss_w_1D, info->gauss_point_1D);
		Gauss_point(info, NG_EXTEND, info->gauss_w_ex, info->gauss_point_ex, info->gauss_w_1D_ex, info->gauss_point_1D_ex);
		counter_GP++;
	}

	if (select_GP == 0)
	{
		if (info->w[0] == info->gauss_w[0])
			return;

		GP_1D = GP_1D_NG;
		GP_2D = GP_2D_NG;
		GP_3D = GP_3D_NG;

		if (info->DIMENSION == 2)
			GP_ON_ELEMENT = GP_2D;
		else if (info->DIMENSION == 3)
			GP_ON_ELEMENT = GP_3D;

		info->w = info->gauss_w;
		info->Gxi = info->gauss_point;
		info->w_1D = info->gauss_w_1D;
		info->Gxi_1D = info->gauss_point_1D;
	}
	else if(select_GP == 1)
	{
		if (info->w[0] == info->gauss_w_ex[0])
			return;

		GP_1D = GP_1D_NG_EXTEND;
		GP_2D = GP_2D_NG_EXTEND;
		GP_3D = GP_3D_NG_EXTEND;

		if (info->DIMENSION == 2)
			GP_ON_ELEMENT = GP_2D;
		else if (info->DIMENSION == 3)
			GP_ON_ELEMENT = GP_3D;

		info->w = info->gauss_w_ex;
		info->Gxi = info->gauss_point_ex;
		info->w_1D = info->gauss_w_1D_ex;
		info->Gxi_1D = info->gauss_point_1D_ex;
	}
}


void Gauss_point(information *info, int temp_GP_1D, double *w_ptr, double *point_ptr, double *w_1D_ptr, double *point_1D_ptr)
{
	if (info->DIMENSION == 2)
	{
		if (temp_GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 5)
		{
			double A = pow((10.0 / 7.0), 0.5);
			double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
			double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
			double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
			double B = pow(70.0, 0.5);
			double w1 = 128.0 / 225.0;
			double w2 = (322.0 + 13.0 * B) / 900.0;
			double w3 = (322.0 - 13.0 * B) / 900.0;
			double w_vec[5] = {w3, w2, w1, w2, w3};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 6)
		{
			double G_vec[6];
			double w_vec[6];

			G_vec[0] = -0.9324695142031520278123;
			G_vec[1] = -0.661209386466264513661;
			G_vec[2] = -0.2386191860831969086305;
			G_vec[3] = 0.238619186083196908631;
			G_vec[4] = 0.661209386466264513661;
			G_vec[5] = 0.9324695142031520278123;

			w_vec[0] = 0.1713244923791703450403;
			w_vec[1] = 0.3607615730481386075698;
			w_vec[2] = 0.4679139345726910473899;
			w_vec[3] = 0.46791393457269104739;
			w_vec[4] = 0.3607615730481386075698;
			w_vec[5] = 0.1713244923791703450403;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 8)
		{
			double G_vec[8];
			double w_vec[8];

			G_vec[0] = -0.9602898564975362316836;
			G_vec[1] = -0.7966664774136267395916;
			G_vec[2] = -0.5255324099163289858177;
			G_vec[3] = -0.1834346424956498049395;
			G_vec[4] = 0.1834346424956498049395;
			G_vec[5] = 0.5255324099163289858177;
			G_vec[6] = 0.7966664774136267395916;
			G_vec[7] = 0.9602898564975362316836;

			w_vec[0] = 0.1012285362903762591525;
			w_vec[1] = 0.2223810344533744705444;
			w_vec[2] = 0.313706645877887287338;
			w_vec[3] = 0.3626837833783619829652;
			w_vec[4] = 0.3626837833783619829652;
			w_vec[5] = 0.313706645877887287338;
			w_vec[6] = 0.222381034453374470544;
			w_vec[7] = 0.1012285362903762591525;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 10)
		{
			double G_vec[10];
			double w_vec[10];

			G_vec[0] = -0.9739065285171717;
			G_vec[1] = -0.8650633666889845;
			G_vec[2] = -0.6794095682990244;
			G_vec[3] = -0.4333953941292472;
			G_vec[4] = -0.1488743389816312;
			G_vec[5] = 0.1488743389816312;
			G_vec[6] = 0.4333953941292472;
			G_vec[7] = 0.6794095682990244;
			G_vec[8] = 0.8650633666889845;
			G_vec[9] = 0.9739065285171717;

			w_vec[0] = 0.0666713443086881;
			w_vec[1] = 0.1494513491505804;
			w_vec[2] = 0.2190863625159820;
			w_vec[3] = 0.2692667193099965;
			w_vec[4] = 0.2955242247147530;
			w_vec[5] = 0.2955242247147530;
			w_vec[6] = 0.2692667193099965;
			w_vec[7] = 0.2190863625159820;
			w_vec[8] = 0.1494513491505804;
			w_vec[9] = 0.0666713443086881;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 20)
		{
			double G_vec[20];
			double w_vec[20];

			G_vec[0]  = -0.99312859918509492478612239;
			G_vec[1]  = -0.96397192727791379126766613;
			G_vec[2]  = -0.91223442825132590586775244;
			G_vec[3]  = -0.83911697182221882339452906;
			G_vec[4]  = -0.74633190646015079261430507;
			G_vec[5]  = -0.6360536807265150254528367;
			G_vec[6]  = -0.51086700195082709800436405;
			G_vec[7]  = -0.3737060887154195606725482;
			G_vec[8]  = -0.2277858511416450780804962;
			G_vec[9]  = -0.07652652113349733375464041;
			G_vec[10] = 0.07652652113349733375464041;
			G_vec[11] = 0.2277858511416450780804962;
			G_vec[12] = 0.37370608871541956067254818;
			G_vec[13] = 0.5108670019508270980043641;
			G_vec[14] = 0.6360536807265150254528367;
			G_vec[15] = 0.74633190646015079261430507;
			G_vec[16] = 0.83911697182221882339452906;
			G_vec[17] = 0.91223442825132590586775244;
			G_vec[18] = 0.9639719272779137912676661;
			G_vec[19] = 0.9931285991850949247861224;

			w_vec[0]  = 0.01761400713915211831186196;
			w_vec[1]  = 0.04060142980038694133103995;
			w_vec[2]  = 0.0626720483341090635695065;
			w_vec[3]  = 0.083276741576704748724758143;
			w_vec[4]  = 0.10193011981724043503675014;
			w_vec[5]  = 0.1181945319615184173123774;
			w_vec[6]  = 0.1316886384491766268984945;
			w_vec[7]  = 0.14209610931838205132929833;
			w_vec[8]  = 0.14917298647260374678782874;
			w_vec[9]  = 0.15275338713072585069808433;
			w_vec[10] = 0.1527533871307258506980843;
			w_vec[11] = 0.1491729864726037467878287;
			w_vec[12] = 0.14209610931838205132929833;
			w_vec[13] = 0.1316886384491766268984945;
			w_vec[14] = 0.11819453196151841731237738;
			w_vec[15] = 0.1019301198172404350367501;
			w_vec[16] = 0.08327674157670474872475814;
			w_vec[17] = 0.0626720483341090635695065;
			w_vec[18] = 0.04060142980038694133103995;
			w_vec[19] = 0.01761400713915211831186196;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
	}
	else if (info->DIMENSION == 3)
	{
		int temp_GP_2D = temp_GP_1D * temp_GP_1D;
		if (temp_GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 5)
		{
			double A = pow((10.0 / 7.0), 0.5);
			double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
			double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
			double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
			double B = pow(70.0, 0.5);
			double w1 = 128.0 / 225.0;
			double w2 = (322.0 + 13.0 * B) / 900.0;
			double w3 = (322.0 - 13.0 * B) / 900.0;
			double w_vec[5] = {w3, w2, w1, w2, w3};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 6)
		{
			double G_vec[6];
			double w_vec[6];

			G_vec[0] = -0.9324695142031520278123;
			G_vec[1] = -0.661209386466264513661;
			G_vec[2] = -0.2386191860831969086305;
			G_vec[3] = 0.238619186083196908631;
			G_vec[4] = 0.661209386466264513661;
			G_vec[5] = 0.9324695142031520278123;

			w_vec[0] = 0.1713244923791703450403;
			w_vec[1] = 0.3607615730481386075698;
			w_vec[2] = 0.4679139345726910473899;
			w_vec[3] = 0.46791393457269104739;
			w_vec[4] = 0.3607615730481386075698;
			w_vec[5] = 0.1713244923791703450403;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 8)
		{
			double G_vec[8];
			double w_vec[8];

			G_vec[0] = -0.9602898564975362316836;
			G_vec[1] = -0.7966664774136267395916;
			G_vec[2] = -0.5255324099163289858177;
			G_vec[3] = -0.1834346424956498049395;
			G_vec[4] = 0.1834346424956498049395;
			G_vec[5] = 0.5255324099163289858177;
			G_vec[6] = 0.7966664774136267395916;
			G_vec[7] = 0.9602898564975362316836;

			w_vec[0] = 0.1012285362903762591525;
			w_vec[1] = 0.2223810344533744705444;
			w_vec[2] = 0.313706645877887287338;
			w_vec[3] = 0.3626837833783619829652;
			w_vec[4] = 0.3626837833783619829652;
			w_vec[5] = 0.313706645877887287338;
			w_vec[6] = 0.222381034453374470544;
			w_vec[7] = 0.1012285362903762591525;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 10)
		{
			double G_vec[10];
			double w_vec[10];

			G_vec[0] = -0.9739065285171717;
			G_vec[1] = -0.8650633666889845;
			G_vec[2] = -0.6794095682990244;
			G_vec[3] = -0.4333953941292472;
			G_vec[4] = -0.1488743389816312;
			G_vec[5] = 0.1488743389816312;
			G_vec[6] = 0.4333953941292472;
			G_vec[7] = 0.6794095682990244;
			G_vec[8] = 0.8650633666889845;
			G_vec[9] = 0.9739065285171717;

			w_vec[0] = 0.0666713443086881;
			w_vec[1] = 0.1494513491505804;
			w_vec[2] = 0.2190863625159820;
			w_vec[3] = 0.2692667193099965;
			w_vec[4] = 0.2955242247147530;
			w_vec[5] = 0.2955242247147530;
			w_vec[6] = 0.2692667193099965;
			w_vec[7] = 0.2190863625159820;
			w_vec[8] = 0.1494513491505804;
			w_vec[9] = 0.0666713443086881;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 20)
		{
			double G_vec[20];
			double w_vec[20];

			G_vec[0]  = -0.99312859918509492478612239;
			G_vec[1]  = -0.96397192727791379126766613;
			G_vec[2]  = -0.91223442825132590586775244;
			G_vec[3]  = -0.83911697182221882339452906;
			G_vec[4]  = -0.74633190646015079261430507;
			G_vec[5]  = -0.6360536807265150254528367;
			G_vec[6]  = -0.51086700195082709800436405;
			G_vec[7]  = -0.3737060887154195606725482;
			G_vec[8]  = -0.2277858511416450780804962;
			G_vec[9]  = -0.07652652113349733375464041;
			G_vec[10] = 0.07652652113349733375464041;
			G_vec[11] = 0.2277858511416450780804962;
			G_vec[12] = 0.37370608871541956067254818;
			G_vec[13] = 0.5108670019508270980043641;
			G_vec[14] = 0.6360536807265150254528367;
			G_vec[15] = 0.74633190646015079261430507;
			G_vec[16] = 0.83911697182221882339452906;
			G_vec[17] = 0.91223442825132590586775244;
			G_vec[18] = 0.9639719272779137912676661;
			G_vec[19] = 0.9931285991850949247861224;

			w_vec[0]  = 0.01761400713915211831186196;
			w_vec[1]  = 0.04060142980038694133103995;
			w_vec[2]  = 0.0626720483341090635695065;
			w_vec[3]  = 0.083276741576704748724758143;
			w_vec[4]  = 0.10193011981724043503675014;
			w_vec[5]  = 0.1181945319615184173123774;
			w_vec[6]  = 0.1316886384491766268984945;
			w_vec[7]  = 0.14209610931838205132929833;
			w_vec[8]  = 0.14917298647260374678782874;
			w_vec[9]  = 0.15275338713072585069808433;
			w_vec[10] = 0.1527533871307258506980843;
			w_vec[11] = 0.1491729864726037467878287;
			w_vec[12] = 0.14209610931838205132929833;
			w_vec[13] = 0.1316886384491766268984945;
			w_vec[14] = 0.11819453196151841731237738;
			w_vec[15] = 0.1019301198172404350367501;
			w_vec[16] = 0.08327674157670474872475814;
			w_vec[17] = 0.0626720483341090635695065;
			w_vec[18] = 0.04060142980038694133103995;
			w_vec[19] = 0.01761400713915211831186196;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
	}
}


// K matrix
void Make_D_Matrix(information *info)
{
	int i, j;

	if (info->DIMENSION == 2)
	{
		if (DM == 0) // 平面応力状態
		{
			double Eone = E / (1.0 - nu * nu);
			double D1[3][3] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) * Eone / 2.0}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (DM == 1) // 平面ひずみ状態
		{
			double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			double D1[3][3] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2.0 / (1.0 - nu) * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
	}
	else if (info->DIMENSION == 3)
	{
		double E_ii = (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double E_ij = nu / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double G = E / (2.0 * (1.0 + nu));

		double D1[6][6] = {{E_ii, E_ij, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ii, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ij, E_ii, 0.0, 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  G  , 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  0.0, G  , 0.0},
						   { 0.0, 0.0,  0.0,  0.0, 0.0, G  }};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
	}
}


// BGマトリックスを求める
double Make_B_Matrix_anypoint(int El_No, double *B, double *Local_coord, information *info)
{
	double a_2x2[2][2], a_3x3[3][3];
	static double *b_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	double *b = &b_parallel[info->DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];

	double J = 0.0;
	int i, j, k;

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_2x2[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k][i] * dShape_func(j, k, Local_coord, El_No, info);
			}
		}

		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			B[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;

			B[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			B[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];

			B[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_3x3[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k][i] * dShape_func(j, k, Local_coord, El_No, info);
			}
		}

		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			B[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];

			B[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			B[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];

			B[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			B[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}

		return J;
	}

	return ERROR;
}


double Make_Jac_anypoint(int El_No, double *Local_coord, information *info)
{
	int i, j, k;
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_2x2[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_2x2(a_2x2);
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_3x3[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_3x3(a_3x3);
	}

	return J;
}

bool GaussianElimination2(double *sol, const double *diff, double *J, int dimension)
{
	const double EPSILON = std::numeric_limits<double>::epsilon();

	// initialize solution
	for (int i = 0; i < dimension; ++i)
		sol[i] = diff[i];

	// Forward elimination
	for (int i = 0; i < dimension; ++i)
	{
		// Find the pivot
		int pivot = i;
		for (int j = i + 1; j < dimension; ++j)
		{
			if (std::abs(J[j * dimension + i]) > std::abs(J[pivot * dimension + i]))
			{
				pivot = j;
			}
		}

		// Swap rows if necessary
		if (pivot != i)
		{
			for (int k = 0; k < dimension; ++k)
			{
				std::swap(J[i * dimension + k], J[pivot * dimension + k]);
			}
			std::swap(sol[i], sol[pivot]);
		}

		// Check for singular matrix
		if (std::abs(J[i * dimension + i]) < EPSILON)
		{
			return false; // Matrix is singular
		}

		// Eliminate column below pivot
		for (int j = i + 1; j < dimension; ++j)
		{
			double factor = J[j * dimension + i] / J[i * dimension + i];
			for (int k = i; k < dimension; ++k)
			{
				J[j * dimension + k] -= factor * J[i * dimension + k];
			}
			sol[j] -= factor * sol[i];
		}
	}

	// Back substitution
	for (int i = dimension - 1; i >= 0; --i)
	{
		sol[i] /= J[i * dimension + i];
		for (int j = 0; j < i; ++j)
		{
			sol[j] -= J[j * dimension + i] * sol[i];
		}
	}

	return true; // Success
}

void Search_ele_point_2D(const int &xi, const int &eta, const int &e_x_max, const int &e_y_max, int *p_x, int *p_y, int *e_x, int *e_y, int *point)
{
	int i;

	*e_x = xi / 2;
	*e_y = eta / 2;
	int p_num_x = xi % 2;
	int p_num_y = eta % 2;

	// 例外
	if (*e_x == e_x_max)
	{
		(*e_x) -= 1;
		p_num_x = 2;
	}
	if (*e_y == e_y_max)
	{
		(*e_y) -= 1;
		p_num_y = 2;
	}

	// ポイント番号を探索
	for (i = 0; i < 9; i++)
	{
		if (p_x[i] == p_num_x && p_y[i] == p_num_y)
		{
			*point = i;
			return;
		}
	}
}


void Search_ele_point_3D(const int &xi, const int &eta, const int &zeta, const int &e_x_max, const int &e_y_max, const int &e_z_max, int *p_x, int *p_y, int *p_z, int *e_x, int *e_y, int *e_z, int *point)
{
	int i;

	*e_x = xi / 2;
	*e_y = eta / 2;
	*e_z = zeta / 2;
	int p_num_x = xi % 2;
	int p_num_y = eta % 2;
	int p_num_z = zeta % 2;

	// 例外
	if (*e_x == e_x_max)
	{
		(*e_x) -= 1;
		p_num_x = 2;
	}
	if (*e_y == e_y_max)
	{
		(*e_y) -= 1;
		p_num_y = 2;
	}
	if (*e_z == e_z_max)
	{
		(*e_z) -= 1;
		p_num_z = 2;
	}

	// ポイント番号を探索
	for (i = 0; i < 27; i++)
	{
		if (p_x[i] == p_num_x && p_y[i] == p_num_y && p_z[i] == p_num_z)
		{
			*point = i;
			return;
		}
	}
}


// 変位場の読みとり
void Substitute_Displacement(information *info, char **argv)
{
	FILE *disp_fp;
	// ファイルの存在チェック
    if ((disp_fp = fopen(argv[Total_mesh + 1], "r")) == NULL) {
        fprintf(stderr, "Eror: No such displacement field file.\n");
        exit(1);
    }
	int temp_i = 0;
	double temp_d = 0;

	// 変位数の読み込み
	fscanf(disp_fp, "%d", &temp_i);
	if (temp_i != info->Total_Control_Point_to_mesh[Total_mesh])
	{
		printf("wrong input's displacements, check POST_ONLY\n");
		printf("number or displacement in files: %d\n", temp_i);
		printf("Total_Control_Point_to_mesh[Total_mesh]: %d\n", info->Total_Control_Point_to_mesh[Total_mesh]);
		exit(1);
	}

	// 次元数の読み込み
	fscanf(disp_fp, "%d", &temp_i);
	if (temp_i != info->DIMENSION)
	{
		printf("wrong input's Dimension, check POST_ONLY\n");
		exit(1);
	}

	// 変位場の読みとり
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		fscanf(disp_fp, "%lf", &temp_d);
		info->Displacement[i] = temp_d;
		// printf("info->Displacement[%d]: %lf\n", i, info->Displacement[i]);
	}
}


// 要素の重なり情報の読み込み
void Substitute_EOI(information *info, char **argv)
{
    FILE *fp_ele_ovelap;
    if ((fp_ele_ovelap = fopen(argv[Total_mesh + 2], "r")) == NULL) 
	{
        printf("Cannot open __EOI file\n");
        exit(1);
    }

    // ヘッダー情報の読み込み
    int total_elements, global_elements;
    if (fscanf(fp_ele_ovelap, "%d\n%d\n", &total_elements, &global_elements) != 2) {
        printf("Error reading header information\n");
        exit(1);
    }

    // 値の整合性チェック
    if (total_elements != info->Total_Element_to_mesh[Total_mesh] ||
        global_elements != info->Total_Element_to_mesh[1]) {
        printf("Error: Inconsistent mesh information in file\n");
        printf("File: total=%d, global=%d\n", total_elements, global_elements);
        printf("Mesh: total=%d, global=%d\n", 
               info->Total_Element_to_mesh[Total_mesh], 
               info->Total_Element_to_mesh[1]);
        exit(1);
    }

    // EOIの初期化
    info->eoi.clear();
    info->eoi.resize(total_elements);

    // EOI情報の読み込み
    for (int e = 0; e < total_elements; e++) {
        int element_num;
        size_t overlap_count;
        if (fscanf(fp_ele_ovelap, "%d %zu", &element_num, &overlap_count) != 2) {
            printf("Error reading element %d information\n", e);
            exit(1);
        }

        // 要素番号の整合性チェック
        if (element_num != e) {
            printf("Error: Inconsistent element numbering at element %d\n", e);
            exit(1);
        }

        // 重なり要素の読み込み
        for (size_t i = 0; i < overlap_count; i++) {
            int overlap_element;
            if (fscanf(fp_ele_ovelap, "%d", &overlap_element) != 1) {
                printf("Error reading overlapping element for element %d\n", e);
                exit(1);
            }
            info->eoi[e].push_back(overlap_element);
        }
    }
    fclose(fp_ele_ovelap);

    // データ検証
    for (int e = 0; e < total_elements; e++) {
        for (size_t i = 0; i < info->eoi[e].size(); i++) {
            if (info->eoi[e][i] >= total_elements) {
                printf("Invalid overlapping element number %d for element %d\n", info->eoi[e][i], e);
				exit(1);
            }
        }
    }
}


void Output_EOI(information *info)
{
    FILE *fp_ele_ovelap;

    // ヘッダー情報
	fp_ele_ovelap = fopen("__EOI.dat", "w");
    fprintf(fp_ele_ovelap, "%d\n", info->Total_Element_to_mesh[Total_mesh]);  // 全要素数
    fprintf(fp_ele_ovelap, "%d\n", info->Total_Element_to_mesh[1]);          // グローバル要素数

    // EOI情報の書き込み
    for (int e = 0; e < info->Total_Element_to_mesh[Total_mesh]; e++) {
        fprintf(fp_ele_ovelap, "%d %zu", e, info->eoi[e].size());
        for (size_t i = 0; i < info->eoi[e].size(); i++) {
            fprintf(fp_ele_ovelap, " %d", info->eoi[e][i]);
        }
        fprintf(fp_ele_ovelap, "\n");
    }

    fclose(fp_ele_ovelap);
}
