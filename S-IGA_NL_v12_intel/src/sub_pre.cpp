// header
#include "_header.hpp"
#include "_sub.hpp"

using namespace std;

// option file with index verification and error handling
void Get_Option(const char *s, vector<string> &opt_files)
{
	ifstream file(s);
	if (!file.is_open())
	{
		cerr << "Error: Cannot open file " << s << endl;
		return;
	}

	map<int, string> indexed_files;
	string line;
	int current_index = -1;
	int max_index = -1;

	while (getline(file, line))
	{
		if (line.empty())
			continue;

		// コメント行から番号を抽出
		if (line[0] == '#')
		{
			size_t bracket_start = line.find('[');
			size_t bracket_end = line.find(']');
			if (bracket_start != string::npos && bracket_end != string::npos)
			{
				string index_str = line.substr(bracket_start + 1, bracket_end - bracket_start - 1);
				
				// 空白を除去
				index_str.erase(remove_if(index_str.begin(), index_str.end(), ::isspace), index_str.end());
				
				// 数字のみかチェック
				if (!index_str.empty() && all_of(index_str.begin(), index_str.end(), ::isdigit))
				{
					try
					{
						current_index = stoi(index_str);
						max_index = max(max_index, current_index);
					}
					catch (const exception &e)
					{
						cerr << "Warning: Invalid index in line: " << line << endl;
						current_index = -1;
					}
				}
				else
				{
					// 数字でない場合はスキップ
					current_index = -1;
				}
			}
			continue;
		}

		// ファイル名行の処理
		size_t commentPos = line.find('#');
		if (commentPos != string::npos)
			line = line.substr(0, commentPos);

		istringstream iss(line);
		string filename;
		if (iss >> filename && current_index >= 0)
		{
			indexed_files[current_index] = filename;
		}
	}

	// vectorに番号順に格納(欠番は空文字列)
	opt_files.clear();
	opt_files.resize(max_index + 1, ""); // 空文字列で初期化

	for (const auto &pair : indexed_files)
		opt_files[pair.first] = pair.second;

	// output loaded files
	cout << "Option files loaded from " << s << " :" << endl;
	for (size_t i = 0; i < opt_files.size(); i++)
		cout << "  [" << i << "] " << (opt_files[i].empty() ? "(not specified)" : opt_files[i]) << endl;
}


// memory allocation
void Allocation(const int num, information *info)
{
	if (num == 0)
	{
		// general
		info->c.ANALYSIS_MODE_0               = any_cast<int>   (info->c.data_vec[ 0]);
		info->c.ANALYSIS_MODE_1               = any_cast<int>   (info->c.data_vec[ 1]);
		info->c.BASE_PATCH_SIGA               = any_cast<int>   (info->c.data_vec[ 2]);
		info->c.SOLVER                        = any_cast<int>   (info->c.data_vec[ 3]);
		info->c.M_MODE                        = any_cast<int>   (info->c.data_vec[ 4]);
		info->c.CALC_ON_GP                    = any_cast<int>   (info->c.data_vec[ 5]);
		info->c.CALC_ON_ELE_VERTEX            = any_cast<int>   (info->c.data_vec[ 6]);
		info->c.OUTPUT_SVG                    = any_cast<int>   (info->c.data_vec[ 7]);
		info->c.OUTPUT_DEFORMED               = any_cast<int>   (info->c.data_vec[ 8]);
		info->c.PARAVIEW_GLO_MODE             = any_cast<int>   (info->c.data_vec[ 9]);
		info->c.PARAVIEW_CRACK_REPRESENTATION = any_cast<int>   (info->c.data_vec[10]);
		info->c.BIN_MODE                      = any_cast<int>   (info->c.data_vec[11]);
		info->c.MODE_EX                       = any_cast<int>   (info->c.data_vec[12]);
		info->c.MAX_OCTREE_N                  = any_cast<int>   (info->c.data_vec[13]);
		info->c.DM                            = any_cast<int>   (info->c.data_vec[14]);
		info->c.NG                            = any_cast<int>   (info->c.data_vec[15]);
		info->c.NG_EXTEND                     = any_cast<int>   (info->c.data_vec[16]);
		info->c.BLENDING                      = any_cast<int>   (info->c.data_vec[17]);

		// B-bar method
		info->c.B_BAR                         = any_cast<int>   (info->c.data_vec[18]);

		// numerical parameters
		info->c.EPS                           = any_cast<double>(info->c.data_vec[19]);

		// nonlinear
		info->c.TOTAL_SECONDS                 = any_cast<double>(info->c.data_vec[20]);
		info->c.STEP_N                        = any_cast<int>   (info->c.data_vec[21]);
		info->c.CUTBACK                       = any_cast<int>   (info->c.data_vec[22]);
		info->c.NEWTON_RAPHSON_MAX_ITR        = any_cast<int>   (info->c.data_vec[23]);
		info->c.NEWTON_RAPHSON_EPS            = any_cast<double>(info->c.data_vec[24]);
		info->c.UPDATE_CONTROL_POINT_MAX_ITR  = any_cast<int>   (info->c.data_vec[25]);

		// fictitous
		info->c.FICT_COEF                     = any_cast<double>(info->c.data_vec[26]);

		info->Total_Knot_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));			  // 各メッシュ上のパッチ数
		info->Total_Patch_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		  // メッシュまでのパッチ数（メッシュ内のパッチ数は含まない）
		info->Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	  // 各メッシュ上のコントロールポイント数
		info->Total_Control_Point_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int)); // メッシュまでのコントロールポイント数(メッシュ内のコントロールポイント数は含まない)
		info->Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		  // 各メッシュ上の要素数
		info->Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		  // メッシュまでの要素数(メッシュ内の要素数は含まない)
		info->real_Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	  // real_Total_Element_on_mesh[MAX_N_MESH]
		info->real_Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));  // real_Total_Element_to_mesh[MAX_N_MESH + 1]
		info->Total_Load_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_Constraint_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_DistributeForce_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		if (info->Total_Knot_to_mesh == NULL || info->Total_Patch_on_mesh == NULL || info->Total_Patch_to_mesh == NULL || info->Total_Control_Point_on_mesh == NULL || info->Total_Control_Point_to_mesh == NULL || info->Total_Element_on_mesh == NULL || info->Total_Element_to_mesh == NULL || info->real_Total_Element_on_mesh == NULL || info->real_Total_Element_to_mesh == NULL || info->Total_Load_to_mesh == NULL || info->Total_Constraint_to_mesh == NULL || info->Total_DistributeForce_to_mesh == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}

		int max_pow_ng = 0;
		if (info->c.MODE_EX == 1)
			max_pow_ng = pow_int(info->c.NG, MAX_DIMENSION);
		else if (info->c.MODE_EX == 0)
			max_pow_ng = pow_int(info->c.NG_EXTEND, MAX_DIMENSION);
		else
		{
			printf("info->c.MODE_EX is not defined\n");
			exit(1);
		}

		info->gauss_w = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point = (double *)malloc(sizeof(double) * max_pow_ng * MAX_DIMENSION);
		info->gauss_w_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_ex = (double *)malloc(sizeof(double) * max_pow_ng * MAX_DIMENSION);
		info->gauss_w_1D = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_1D = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_w_1D_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_1D_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		if (info->gauss_w == NULL || info->gauss_point == NULL || info->gauss_w_ex == NULL || info->gauss_point_ex == NULL || info->gauss_w_1D == NULL || info->gauss_point_1D == NULL || info->gauss_w_1D_ex == NULL || info->gauss_point_1D_ex == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 1)
	{
		info->Order = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));								  // Order[MAX_N_PATCH][info->DIMENSION]
		info->No_knot = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));							  // No_knot[MAX_N_PATCH][info->DIMENSION]
		info->Total_Control_Point_to_patch = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] + 1), sizeof(int));					  // Total_Control_Point_to_patch[MAX_N_PATCH]
		info->Total_Knot_to_patch_dim = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION + 1), sizeof(int));		  // Total_Knot_to_patch_dim[MAX_N_PATCH][info->DIMENSION]
		info->Position_Knots = (double *)malloc(sizeof(double) * info->Total_Knot_to_mesh[Total_mesh]);									  // Position_Knots[MAX_N_PATCH][info->DIMENSION][MAX_N_KNOT];
		info->No_Control_point = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));					  // Order[MAX_N_PATCH][info->DIMENSION]
		info->No_Control_point_in_patch = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));							  // No_Control_point_in_patch[MAX_N_PATCH]
		info->Patch_Control_point = (int *)malloc(sizeof(int) * MAX_CP);																  // Patch_Control_point[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]
		info->No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));						  // No_Control_point_ON_ELEMENT[MAX_N_PATCH]
		info->Node_Coordinate = (double *)malloc(sizeof(double) * (info->Total_Control_Point_to_mesh[Total_mesh] * (info->DIMENSION + 1))); // Node_Coordinate[MAX_N_NODE][info->DIMENSION + 1];
		info->Control_Coord_x = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Coord_y = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Coord_z = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Weight = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Weight[MAX_N_NODE];
		info->Constraint_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Constraint_to_mesh[Total_mesh] * 2));						  // Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
		info->Value_of_Constraint = (double *)calloc(info->Total_Constraint_to_mesh[Total_mesh], sizeof(double));						  // Value_of_Constraint[MAX_N_CONSTRAINT];
		info->Load_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Load_to_mesh[Total_mesh] * 2));									  // Load_Node_Dir[MAX_N_LOAD][2];
		info->Value_of_Load = (double *)calloc(info->Total_Load_to_mesh[Total_mesh], sizeof(double));									  // Value_of_Load[MAX_N_LOAD];
		info->iPatch_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // iPatch_array[MAX_N_DISTRIBUTE_FORCE]
		info->iCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // iCoord_array[MAX_N_DISTRIBUTE_FORCE]
		info->jCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // jCoord_array[MAX_N_DISTRIBUTE_FORCE]
		info->type_load_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);							  // type_load_array[MAX_N_DISTRIBUTE_FORCE]
		info->val_Coord_array = (double *)calloc(info->Total_DistributeForce_to_mesh[Total_mesh], sizeof(double));						  // val_Coord_array[MAX_N_DISTRIBUTE_FORCE]
		info->Range_Coord_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 2 * 2), sizeof(double));			  // Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2 * 2] (info->DIMENSION == 3 のとき 2 -> 4)
		info->Coeff_Dist_Load_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 3 * 2), sizeof(double));		  // Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3 * 2] (info->DIMENSION == 3 のとき 3 -> 6)
		if (info->Order == NULL || info->No_knot == NULL || info->Total_Control_Point_to_patch == NULL || info->Total_Knot_to_patch_dim == NULL || info->Position_Knots == NULL || info->No_Control_point == NULL || info->No_Control_point_in_patch == NULL || info->Patch_Control_point == NULL || info->No_Control_point_ON_ELEMENT == NULL || info->Node_Coordinate == NULL || info->Control_Coord_x == NULL || info->Control_Coord_y == NULL || info->Control_Weight == NULL || info->Constraint_Node_Dir == NULL || info->Value_of_Constraint == NULL || info->Load_Node_Dir == NULL || info->Value_of_Load == NULL || info->iPatch_array == NULL || info->iCoord_array == NULL || info->jCoord_array == NULL || info->type_load_array == NULL || info->val_Coord_array == NULL || info->Range_Coord_array == NULL || info->Coeff_Dist_Load_array == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 2)
	{
		info->INC = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION));		 // INC[MAX_N_PATCH][MAX_N_NODE][info->DIMENSION]
		info->Controlpoint_of_Element = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * MAX_NO_CP_ON_ELEMENT));							 // Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CP_ON_ELEMENT]
		info->Element_patch = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // Element_patch[MAX_N_ELEMENT]
		info->Element_mesh = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // Element_mesh[MAX_N_ELEMENT] 要素がどのメッシュ内にあるかを示す配列
		info->line_No_real_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // line_No_real_element[MAX_N_PATCH][info->DIMENSION] ゼロエレメントではない要素列の数
		info->line_No_Total_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // line_No_Total_element[MAX_N_PATCH][info->DIMENSION] ゼロエレメントを含むすべての要素列の数
		info->difference = (double *)calloc(info->Total_Knot_to_mesh[Total_mesh], sizeof(double));																 // difference[MAX_N_PATCH][info->DIMENSION][MAX_N_KNOT] 隣り合うノットベクトルの差
		info->Total_element_all_ID = (int *)calloc(info->Total_Element_to_mesh[Total_mesh], sizeof(int));														 // Total_element_all_ID[MAX_N_ELEMENT] ゼロエレメントではない要素 = 1, ゼロエレメント = 0
		info->ENC = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION));														 // ENC[MAX_N_ELEMENT][info->DIMENSION] ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
		info->real_element_line = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION)); // real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][info->DIMENSION] ゼロエレメントではない要素列
		info->real_element = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // real_element[MAX_N_ELEMENT] ゼロエレメントではない要素の番号
		info->real_El_No_on_mesh = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh]));				 // real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT]
		info->Equivalent_Nodal_Force = (double *)calloc(MAX_CP * info->DIMENSION, sizeof(double));																 // Equivalent_Nodal_Force[MAX_N_NODE][info->DIMENSION] Equivalent nodal forces arising from the distributed load
		if (info->INC == NULL || info->Controlpoint_of_Element == NULL || info->Element_patch == NULL || info->Element_mesh == NULL || info->line_No_real_element == NULL || info->line_No_Total_element == NULL || info->difference == NULL || info->Total_element_all_ID == NULL || info->ENC == NULL || info->real_element_line == NULL || info->real_element == NULL || info->real_El_No_on_mesh == NULL || info->Equivalent_Nodal_Force == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 3)
	{
		info->gp.resize(info->Total_Element_to_mesh[Total_mesh]);
		info->eoi.resize(info->Total_Element_to_mesh[Total_mesh]);
		info->octree_subcell.resize(info->Total_Element_to_mesh[Total_mesh]);
	}
	else if (num == 4)
	{
		info->D = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * D_MATRIX_SIZE);
		info->Index_Dof = (int *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));													   // Index_Dof[MAX_K_WHOLE_SIZE];
		info->K_Whole_Ptr = (int *)calloc(MAX_K_WHOLE_SIZE + 1, sizeof(int));											   // K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1]
		if (info->D == NULL || info->Index_Dof == NULL || info->K_Whole_Ptr == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 5)
	{
		info->K_Whole_Col = (int *)malloc(sizeof(int) * info->K_Whole_Ptr[K_Whole_Size]);	  // K_Whole_Col[MAX_NON_ZERO]
		info->K_Whole_Val = (double *)calloc(info->K_Whole_Ptr[K_Whole_Size], sizeof(double)); // K_Whole_Val[MAX_NON_ZERO]
		if (info->K_Whole_Col == NULL || info->K_Whole_Val == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
		{
			info->full_K = (double *)calloc((long)K_Whole_Size * K_Whole_Size, sizeof(double));
			if (info->full_K == NULL)
			{
				printf("Cannot allocate memory\n");
				exit(1);
			}
		}
	}
	else if (num == 6)
	{
		info->sol_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double)); // sol_vec[MAX_K_WHOLE_SIZE]
		info->rhs_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double)); // rhs_vec[MAX_K_WHOLE_SIZE]
		if (info->sol_vec == NULL || info->rhs_vec == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 7)
	{
		int gp_n = pow_int(info->c.NG, info->DIMENSION);
		info->Displacement = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE);
		info->Displacement_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * gp_n * info->DIMENSION, sizeof(double));
		info->Strain_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * gp_n * N_STRAIN, sizeof(double));
		info->Stress_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * gp_n * N_STRESS, sizeof(double));
		info->Displacement_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(2, info->DIMENSION) * info->DIMENSION, sizeof(double));
		info->Strain_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(2, info->DIMENSION) * N_STRAIN, sizeof(double));
		info->Stress_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(2, info->DIMENSION) * N_STRESS, sizeof(double));
		info->ReactionForce = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double));
		if (info->Displacement == NULL || info->Strain_at_GP == NULL || info->Stress_at_GP == NULL || info->ReactionForce == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 8)
	{
		int total_array = 0, check_n = 0, FE_n = 0, opp_n = 0, point_on_element = pow_int(3, info->DIMENSION);
		int total_edge = 0, edge_on_ele = 0;
		if (info->DIMENSION == 2)
		{
			for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
			{
				int x = info->line_No_Total_element[i * info->DIMENSION] * 2 + 1, y = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;
				total_array += 2 * (x + y);
			}
			check_n = 16;
			FE_n = 32;
			opp_n = 4;
			edge_on_ele = 2;
		}
		else if (info->DIMENSION == 3)
		{
			for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
			{
				int x = info->line_No_Total_element[i * info->DIMENSION] * 2 + 1, y = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1, z = info->line_No_Total_element[i * info->DIMENSION + 2] * 2 + 1;
				total_array += 2 * (x * y + y * z + x * z);
			}
			check_n = 24;
			FE_n = 36;
			opp_n = 6;
			edge_on_ele = 4;
		}
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				total_edge += info->line_No_Total_element[i * info->DIMENSION + j] * edge_on_ele * 3;
			}
		}
		info->Connectivity = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_ele = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_point = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_coord = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * info->DIMENSION, sizeof(double));
		info->Patch_check = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * check_n);
		info->Patch_array = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * total_array);
		info->Face_Edge_info = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * FE_n);
		info->Opponent_patch_num = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * opp_n);
		info->Edge_coord = (double *)calloc(total_edge * info->DIMENSION, sizeof(double));
		info->disp_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * info->DIMENSION, sizeof(double));
		info->strain_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRAIN, sizeof(double));
		info->stress_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRESS, sizeof(double));
		info->vm_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		info->hs_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh] * point_on_element; i++)
		{
			info->Connectivity[i] = -1;
		}
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh] * FE_n; i++)
		{
			info->Face_Edge_info[i] = -1;
		}
		if (info->Connectivity == NULL || info->Connectivity_ele == NULL || info->Connectivity_point == NULL || info->Connectivity_coord == NULL || info->Patch_check == NULL || info->Patch_array == NULL || info->Face_Edge_info == NULL || info->Opponent_patch_num == NULL || info->disp_at_connectivity == NULL || info->strain_at_connectivity == NULL || info->stress_at_connectivity == NULL || info->vm_at_connectivity == NULL || info->hs_at_connectivity == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		if (info->c.OUTPUT_DEFORMED == 0)
		{
			info->deformed_Connectivity_coord = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * info->DIMENSION, sizeof(double));
			info->deformed_Edge_coord = (double *)calloc(total_edge * info->DIMENSION, sizeof(double));
			if (info->deformed_Connectivity_coord == NULL || info->deformed_Edge_coord == NULL)
			{
				printf("Cannot allocate memory\n");
				exit(1);
			}
		}

		info->Connectivity_all_ele.resize(info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->ecn.resize(info->Total_Element_to_mesh[Total_mesh]);
		info->elastic_strain_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRAIN, sizeof(double));
		info->equivalent_plastic_strain_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		info->yield_stress_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		info->back_stress_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRESS, sizeof(double));
		if (info->elastic_strain_at_connectivity == NULL || info->equivalent_plastic_strain_at_connectivity == NULL || info->yield_stress_at_connectivity == NULL || info->back_stress_at_connectivity == NULL)
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
		MAX_ORDER += 1; // 次数の最大値 + 1
		MAX_NO_CP_ON_ELEMENT = pow_int(MAX_ORDER, info->DIMENSION);
		MAX_KIEL_SIZE = MAX_NO_CP_ON_ELEMENT * info->DIMENSION;
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
		MAX_K_WHOLE_SIZE = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;
	}
}


// Read file 1st time
void Get_Input_1(int tm, const char *filename, information *info)
{
	char s[256];
	int temp_i;

	int i, j;

	FILE *fp = fopen(filename, "r");

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
	vector<int> CP(temp_i * info->DIMENSION);
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
	printf("Total_Control_Point_to_mesh[%d] = %d\n", tm, info->Total_Control_Point_to_mesh[tm]);

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
}


// Read file 2nd time
void Get_Input_2(int tm, const char *filename, information *info)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j, k;

	FILE *fp = fopen(filename, "r");

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
			printf("Order[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
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
			printf("No_knot[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
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
			printf("No_Control_point[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
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

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_in_patch[%d] = %d\t", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]);
	}
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

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_ON_ELEMENT[%d] = %d\n", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]]);
	}

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

	// init
	info->Controlpoint_of_Element_in_patch.resize(info->Total_Element_to_mesh[Total_mesh] * MAX_NO_CP_ON_ELEMENT);

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
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
					{
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;

						if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1])
						{
							for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
								{
									B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
									b = jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
									info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;

									info->Controlpoint_of_Element_in_patch[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = i - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc;
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

											info->Controlpoint_of_Element_in_patch[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = i - kkloc * info->No_Control_point[patch * info->DIMENSION + 0] * info->No_Control_point[patch * info->DIMENSION + 1] - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc;
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

				// Setting_Dist_Load_2D(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
				setDistLoad(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
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

				// Setting_Dist_Load_3D(tm, iPatch, iCoord, jCoord, val_Coord, iRange_Coord, jRange_Coord, type_load, iCoeff_Dist_Load, jCoeff_Dist_Load, info);
				setDistLoad(tm, iPatch, iCoord, jCoord, val_Coord, iRange_Coord, jRange_Coord, type_load, iCoeff_Dist_Load, jCoeff_Dist_Load, info);
			}
		}
	}
}


#if 0
// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info)
{
	int GP_1D = info->c.NG;
	int iii, jjj;
	int jCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load;
	int iX, iY;
	int ic, ig;
	double val_jCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	static int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	static int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(info);

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

	for (iii = info->Order[iPatch * info->DIMENSION + iCoord]; iii < info->No_knot[iPatch * info->DIMENSION + iCoord] - info->Order[iPatch * info->DIMENSION + iCoord] - 1; iii++)
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

	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(1);
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
				Local_Coord[iCoord] = info->gauss_point_1D[ig];

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
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * info->gauss_w_1D[ig];
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
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig];
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
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig];
					}
				}
			}
		}
	}
}


void Setting_Dist_Load_3D(int mesh_n, int iPatch, int iCoord, int jCoord, double val_Coord, double *iRange_Coord, double *jRange_Coord, int type_load, double *iCoeff_Dist_Load, double *jCoeff_Dist_Load, information *info)
{
	int GP_1D = info->c.NG;
	int iii, jjj, kkk;
	int kCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000}, kPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load_iDir, No_Element_For_Dist_Load_jDir;
	int iX, iY, iZ;
	int ic, ig_i, ig_j;
	double val_kCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	static int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh] * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	static int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(info);

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
						Local_Coord[iCoord] = info->gauss_point_1D[ig_i];
						Local_Coord[jCoord] = info->gauss_point_1D[ig_j];

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
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * info->gauss_w_1D[ig_i] * info->gauss_w_1D[ig_j];
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
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] +=   LoadDir[0] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig_i] * info->gauss_w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] +=   LoadDir[1] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig_i] * info->gauss_w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 2] += - LoadDir[2] * valDistLoad * sfc * detJ * info->gauss_w_1D[ig_i] * info->gauss_w_1D[ig_j];
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
#endif


void setDistLoad(int current_mesh, int patch, int coord, double target_knot, double *range, int type_load, double *dist_load_coeff, information *info)
{
	Make_gauss_array(info);

	// check coordinates and set load surface
	int target_coord = -1;
	if (coord == 0)
		target_coord = 1;
	else if (coord == 1)
		target_coord = 0;
	else
	{
		printf("Error, incorrect input data at distributed load.\n");
		exit(1);
	}

	// set target ENC
	int current_patch = patch + info->Total_Patch_to_mesh[current_mesh];
	int id = current_patch * info->DIMENSION + target_coord;
	int start_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->Order[id]];
	int end_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->No_knot[id] - info->Order[id] - 1];

	bool isStart = false;
	int target_ENC = -1;
	if (fabs(target_knot - start_knot) < MERGE_ERROR)
	{
		target_ENC = 0;
		isStart = true;
	}
	else if (fabs(target_knot - end_knot) < MERGE_ERROR)
	{
		target_ENC = info->No_Control_point[id] - info->Order[id] - 1;
		isStart = false;
	}
	else
	{
		printf("Error, incorrect input data at distributed load.\ntarget_knot is invalid.\n");
		exit(1);
	}
	
	// set coord ENC
	vector<int> ENC(2);
	bool errorFlag[2] = {false, false};
	int id_coord = current_patch * info->DIMENSION + coord;
	int knot_id_coord = info->Total_Knot_to_patch_dim[id_coord] + info->Order[id_coord];
	for (int i = 0; i < info->No_Control_point[id_coord] - info->Order[id_coord]; i++)
	{
		if (fabs(range[0] - info->Position_Knots[knot_id_coord + i]) < MERGE_ERROR)
		{
			ENC[0] = i;
			errorFlag[0] |= true;
		}
		if (fabs(range[1] - info->Position_Knots[knot_id_coord + i + 1]) < MERGE_ERROR)
		{
			ENC[1] = i;
			errorFlag[1] |= true;
		}
	}
	if (!errorFlag[0] && !errorFlag[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// serach element
	vector<int> ele_list;
	for (int i = 0; i < info->Total_Element_on_mesh[current_mesh]; i++)
	{
		int e = i + info->Total_Element_to_mesh[current_mesh];
		if (info->Element_patch[e] == current_patch)
		{
			bool isTargetEle = true;
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int current_ENC = info->ENC[e * info->DIMENSION + j];

				// target coord
				if (j == target_coord)
				{
					if (current_ENC == target_ENC)
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range
				else if (j == coord)
				{
					if (ENC[0] <= current_ENC && current_ENC <= ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}
			}

			if (isTargetEle)
				ele_list.emplace_back(e);
		}
	}

	// parametric coordinates
	int gp_1d = info->c.NG;
	vector<double> gp_para(gp_1d * info->DIMENSION);
	vector<double> gp_w(gp_1d);
	for (int i = 0; i < gp_1d; i++)
	{
		gp_para[i * info->DIMENSION + coord] = info->gauss_point_1D[i];
		gp_para[i * info->DIMENSION + target_coord] = isStart ? -1.0 : 1.0;
		gp_w[i] = info->gauss_w_1D[i];
	}

	// integration for distributed load
	for (size_t i = 0; i < ele_list.size(); i++)
	{
		int e = ele_list[i];
		for (int j = 0; j < gp_1d; j++)
		{
			double *para = gp_para.data() + j * info->DIMENSION;
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
			shape_and_dshape(R.data(), dR.data(), para, e, true, info);

			// jacobian matrix
			vector<double> jac(info->DIMENSION, 0.0);
			for (int k = 0; k < info->DIMENSION; k++)
				for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
					for (int l = 0; l < info->DIMENSION; l++)
						if (l == coord)
							jac[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];

			// jacobian ||normal_J||
			vector<double> normal_J(info->DIMENSION);
			normal_J[0] = jac[1];
			normal_J[1] = -jac[0];
			double J = 0.0;
			for (int k = 0; k < info->DIMENSION; k++)
			{
				J += normal_J[k] * normal_J[k];
			}
			J = sqrt(J);

			// integrand
			double *c = dist_load_coeff;
			vector<double> current_knot(info->DIMENSION);
			ShapeFunc_from_paren(current_knot.data(), para, coord, e, info);
			double x = current_knot[coord];

			double val = (c[0] + c[1] * x + c[2] * x * x);

			// type 0, 1
			if (type_load < 3)
			{
				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + type_load;
					info->Equivalent_Nodal_Force[index] += coeff * R[k];
				}
			}

			// type 3 (normal direction)
			else if (type_load == 2)
			{
				vector<double> normal(info->DIMENSION);
				normal[0] = normal_J[0] / J;
				normal[1] = normal_J[1] / J;

				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION;
					info->Equivalent_Nodal_Force[index + 0] += coeff * normal[0] * R[k];
					info->Equivalent_Nodal_Force[index + 1] += coeff * normal[1] * R[k];
				}
			}

			else
			{
				printf("Error, incorrect input data at distributed load.\nChange the type of the distributed load.\n");
				exit(1);
			}
		}
	}
}


void setDistLoad(int current_mesh, int patch, int coord_i, int coord_j, double target_knot, double *range_i, double *range_j, int type_load, double *dist_load_coeff_i, double *dist_load_coeff_j, information *info)
{
	Make_gauss_array(info);

	// check coordinates and set load surface
	int target_coord = -1;
	if (coord_i == 0 && coord_j == 1)
		target_coord = 2;
	else if (coord_i == 1 && coord_j == 2)
		target_coord = 0;
	else if (coord_i == 2 && coord_j == 0)
		target_coord = 1;
	else
	{
		printf("Error, incorrect input data at distributed load.\n");
		exit(1);
	}

	// set target ENC
	int current_patch = patch + info->Total_Patch_to_mesh[current_mesh];
	int id = current_patch * info->DIMENSION + target_coord;
	int start_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->Order[id]];
	int end_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->No_knot[id] - info->Order[id] - 1];

	bool isStart = false;
	int target_ENC = -1;
	if (fabs(target_knot - start_knot) < MERGE_ERROR)
	{
		target_ENC = 0;
		isStart = true;
	}
	else if (fabs(target_knot - end_knot) < MERGE_ERROR)
	{
		target_ENC = info->No_Control_point[id] - info->Order[id] - 1;
		isStart = false;
	}
	else
	{
		printf("Error, incorrect input data at distributed load.\ntarget_knot is invalid.\n");
		exit(1);
	}

	// set coord_i ENC
	vector<int> i_ENC(2);
	bool errorFlag_i[2] = {false, false};
	int id_i = current_patch * info->DIMENSION + coord_i;
	int knot_id_i = info->Total_Knot_to_patch_dim[id_i] + info->Order[id_i];
	for (int i = 0; i < info->No_Control_point[id_i] - info->Order[id_i]; i++)
	{
		if (fabs(range_i[0] - info->Position_Knots[knot_id_i + i]) < MERGE_ERROR)
		{
			i_ENC[0] = i;
			errorFlag_i[0] |= true;
		}
		if (fabs(range_i[1] - info->Position_Knots[knot_id_i + i + 1]) < MERGE_ERROR)
		{
			i_ENC[1] = i;
			errorFlag_i[1] |= true;
		}
	}
	if (!errorFlag_i[0] && !errorFlag_i[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// set coord_j ENC
	vector<int> j_ENC(2);
	bool errorFlag_j[2] = {false, false};
	int id_j = current_patch * info->DIMENSION + coord_j;
	int knot_id_j = info->Total_Knot_to_patch_dim[id_j] + info->Order[id_j];
	for (int i = 0; i < info->No_Control_point[id_j] - info->Order[id_j]; i++)
	{
		if (fabs(range_j[0] - info->Position_Knots[knot_id_j + i]) < MERGE_ERROR)
		{
			j_ENC[0] = i;
			errorFlag_j[0] |= true;
		}
		if (fabs(range_j[1] - info->Position_Knots[knot_id_j + i + 1]) < MERGE_ERROR)
		{
			j_ENC[1] = i;
			errorFlag_j[1] |= true;
		}
	}
	if (!errorFlag_j[0] && !errorFlag_j[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// serach element
	vector<int> ele_list;
	for (int i = 0; i < info->Total_Element_on_mesh[current_mesh]; i++)
	{
		int e = i + info->Total_Element_to_mesh[current_mesh];
		if (info->Element_patch[e] == current_patch)
		{
			bool isTargetEle = true;
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int current_ENC = info->ENC[e * info->DIMENSION + j];

				// target coord
				if (j == target_coord)
				{
					if (current_ENC == target_ENC)
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range i
				else if (j == coord_i)
				{
					if (i_ENC[0] <= current_ENC && current_ENC <= i_ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range j
				else if (j == coord_j)
				{
					if (j_ENC[0] <= current_ENC && current_ENC <= j_ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}
			}

			if (isTargetEle)
				ele_list.emplace_back(e);
		}
	}

	// parametric coordinates
	int gp_1d = info->c.NG;
	int gp_2d = info->c.NG * info->c.NG;
	vector<double> gp_para(gp_2d * info->DIMENSION);
	vector<double> gp_w(gp_2d);
	for (int i = 0; i < gp_1d; i++)
		for (int j = 0; j < gp_1d; j++)
		{
			gp_para[(i * gp_1d + j) * info->DIMENSION + coord_i] = info->gauss_point_1D[i];
			gp_para[(i * gp_1d + j) * info->DIMENSION + coord_j] = info->gauss_point_1D[j];
			gp_para[(i * gp_1d + j) * info->DIMENSION + target_coord] = isStart ? -1.0 : 1.0;
			gp_w[i * gp_1d + j] = info->gauss_w_1D[i] * info->gauss_w_1D[j];
		}

	// integration for distributed load
	for (size_t i = 0; i < ele_list.size(); i++)
	{
		int e = ele_list[i];
		for (int j = 0; j < gp_2d; j++)
		{
			double *para = gp_para.data() + j * info->DIMENSION;
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
			shape_and_dshape(R.data(), dR.data(), para, e, true, info);

			// jacobian matrix
			vector<double> jac_i(info->DIMENSION, 0.0);
			vector<double> jac_j(info->DIMENSION, 0.0);
			for (int k = 0; k < info->DIMENSION; k++)
				for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						if (l == coord_i)
							jac_i[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
						else if (l == coord_j)
							jac_j[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
					}

			// jacobian ||cross(jac_i, jac_j)||
			double J = 0.0;
			vector<double> normal_J(info->DIMENSION);
			for (int k = 0; k < info->DIMENSION; k++)
			{
				normal_J[k] = jac_i[(k + 1) % info->DIMENSION] * jac_j[(k + 2) % info->DIMENSION] - jac_i[(k + 2) % info->DIMENSION] * jac_j[(k + 1) % info->DIMENSION];
				J += normal_J[k] * normal_J[k];
			}
			J = sqrt(J);

			// integrand
			double *ci = dist_load_coeff_i;
			double *cj = dist_load_coeff_j;
			vector<double> current_knot(info->DIMENSION);
			ShapeFunc_from_paren(current_knot.data(), para, coord_i, e, info);
			ShapeFunc_from_paren(current_knot.data(), para, coord_j, e, info);
			double x = current_knot[coord_i];
			double y = current_knot[coord_j];

			double val = (ci[0] + ci[1] * x + ci[2] * x * x) * (cj[0] + cj[1] * y + cj[2] * y * y);

			// type 0, 1, 2
			if (type_load < 3)
			{
				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + type_load;
					info->Equivalent_Nodal_Force[index] += coeff * R[k];
				}
			}

			// type 3 (normal direction)
			else if (type_load == 3)
			{
				vector<double> normal(info->DIMENSION);
				normal[0] = normal_J[0] / J;
				normal[1] = normal_J[1] / J;
				normal[2] = normal_J[2] / J;

				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION;
					info->Equivalent_Nodal_Force[index + 0] += coeff * normal[0] * R[k];
					info->Equivalent_Nodal_Force[index + 1] += coeff * normal[1] * R[k];
					info->Equivalent_Nodal_Force[index + 2] += coeff * normal[2] * R[k];
				}
			}

			else
			{
				printf("Error, incorrect input data at distributed load.\nChange the type of the distributed load.\n");
				exit(1);
			}
		}
	}

	// exit(0);
}


// search element
int ele_check(int patch_n, double *para_coord, information *info)
{
	int line[MAX_DIMENSION] = {0};
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int line_counter = 0;
		for (int j = 0; j < info->No_Control_point[patch_n * info->DIMENSION + i] - info->Order[patch_n * info->DIMENSION + i]; j++)
		{
			double knot_front = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j];
			double knot_back = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j + 1];

			if (fabs(knot_front - para_coord[i]) < MERGE_ERROR)
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


// gauss points
void Make_gauss_array(information *info)
{
	if (!info->isSetGaussArray)
	{
		Gauss_point(info, info->c.NG, info->gauss_w, info->gauss_point, info->gauss_w_1D, info->gauss_point_1D);
		if (info->c.MODE_EX == 0)
			Gauss_point(info, info->c.NG_EXTEND, info->gauss_w_ex, info->gauss_point_ex, info->gauss_w_1D_ex, info->gauss_point_1D_ex);
		info->isSetGaussArray = true;
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
		else
		{
			printf("Error, incorrect input data at gauss points.\nChange the number of the gauss points.\n");
			exit(1);
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
		else
		{
			printf("Error, incorrect input data at gauss points.\nChange the number of the gauss points.\n");
			exit(1);
		}
	}
}


void Make_Gauss_points(bool isSinglePatch, information *info)
{
	int n[2] = {0};
	n[0] = pow_int(info->c.NG, info->DIMENSION), n[1] = pow_int(info->c.NG_EXTEND, info->DIMENSION);

	Make_gauss_array(info);
	if (isSinglePatch)
	{
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			info->gp[i].setVar(false, false, p, i, n[0], info);
		}
	}
	else
	{
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];

			// small deformation problem
			if (info->c.MODE_EX == 0 && info->c.ANALYSIS_MODE_1 == 0)
			{
				if (info->eoi[i].size() > 1 && info->c.BASE_PATCH_SIGA == 0)
				{
					if (i < info->Total_Element_to_mesh[1])
						info->gp[i].setVar(false, true, p, i, n[0], info);
					else
						info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				}
				else if (info->eoi[i].size() > 1 && info->c.BASE_PATCH_SIGA == 1)
				{
					if (i < info->Total_Element_to_mesh[1])
						info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
					else
						info->gp[i].setVar(false, true, p, i, n[0], info);
				}
				else
					info->gp[i].setVar(false, true, p, i, n[0], info);
			}

			// large deformation problem
			else if (info->c.MODE_EX == 0 && info->c.ANALYSIS_MODE_1 == 1)
			{
				// ローカル重複部
				#if 0
				if (i >= info->Total_Element_to_mesh[1] && info->eoi[i].size() > 1)
					info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				else
					info->gp[i].setVar(false, true, p, i, n[0], info);
				#endif

				// ローカル重複部とグローバル重複部
				#if 0
				if (info->eoi[i].size() > 1)
					info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				else
					info->gp[i].setVar(false, true, p, i, n[0], info);
				#endif

				// ローカル重複部(2要素以上)とグローバル重複部(1要素以上)
				#if 1
				if (i < info->Total_Element_to_mesh[1] && info->eoi[i].size() > 0)
					info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				else if (i >= info->Total_Element_to_mesh[1] && info->eoi[i].size() > 1)
					info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				else
					info->gp[i].setVar(false, true, p, i, n[0], info);
				#endif

				// ローカル全てとグローバル重複部
				#if 0
				if (i >= info->Total_Element_to_mesh[1] || info->eoi[i].size() > 1)
					info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				else
					info->gp[i].setVar(false, true, p, i, n[0], info);
				#endif

				// 全て
				#if 0
				info->gp[i].setVar(true, true, p, i, n[0], n[1], info);
				#endif
			}

			// less gauss points
			else if (info->c.MODE_EX == 1)
				info->gp[i].setVar(false, true, p, i, n[0], info);

			// error message
			else
			{
				printf("info->c.MODE_EX is not correct.\n");
				exit(1);
			}
		}
	}
}


void gp_switch(bool flag, information *info)
{
	// if flag == true, switch to the extended gauss points. if flag == false, switch to the original gauss points.
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		info->gp[i].switch_gp(flag);
	}
}


// K matrix
void Make_D_Matrix(information *info)
{
	int i, j;

	if (info->DIMENSION == 2)
	{
		if (info->c.DM == 0) // 平面応力状態
		{
			double Eone = E / (1.0 - nu * nu);
			double D1[3][3] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) * Eone / 2.0}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (info->c.DM == 1) // 平面ひずみ状態
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


// 拘束されている行数を省いた行列の番号の制作
void Make_Index_Dof(information *info)
{
	int i, k = 0;

	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		info->Index_Dof[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		if (info->Index_Dof[i] != ERROR)
		{
			info->Index_Dof[i] = k;
			k++;
		}
	}

	K_Whole_Size = k;
	printf("\nK_Whole_Size = %d\n\n", k);
}


#if 1
void Make_K_Whole_Ptr_Col(information *info, int mode_select)
{
	// mode_select == 0: Ptr を作成
	// mode_select == 1: Col を作成
	
	vector<set<int>> &ptr = info->ptr;

	if (mode_select == 0)
	{
		ptr.resize(K_Whole_Size);
		vector<std::mutex> mtx(K_Whole_Size);

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					if (row >= 0)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[p]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
							{
								int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + m];
								if (col >= row)
								{
									std::lock_guard<std::mutex> lock(mtx[row]);
									ptr[row].insert(col);
								}
							}
				}

			// coupling element
			if (i >= info->Total_Element_to_mesh[1])
				for (size_t j = 0; j < info->eoi[i].size(); j++)
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							if (row >= 0)
								for (int m = 0; m < info->No_Control_point_ON_ELEMENT[p]; m++)
									for (int n = 0; n < info->DIMENSION; n++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + n];
										if (col >= row)
										{
											std::lock_guard<std::mutex> lock(mtx[row]);
											ptr[row].insert(col);
										}
									}
						}
		}

		// substitute ptr to K_Whole_Ptr
		info->K_Whole_Ptr[0] = 0;
		for (int i = 0; i < K_Whole_Size; i++)
			info->K_Whole_Ptr[i + 1] = info->K_Whole_Ptr[i] + ptr[i].size();
	}
	
	else if (mode_select == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < K_Whole_Size; i++)
		{
			int count = 0;
			for (auto itr = ptr[i].begin(); itr != ptr[i].end(); ++itr)
			{
				info->K_Whole_Col[info->K_Whole_Ptr[i] + count] = *itr;
				count++;
			}
		}
	}
}
#else
void Make_K_Whole_Ptr_Col(information *info, int mode_select)
{
	// mode_select == 0: Ptr を作成
	// mode_select == 1: Col を作成

	constexpr int K_DIVISION_LENGE = 10;
	// Node_To_Node = (int *)malloc(sizeof(int) * K_DIVISION_LENGE * info->Total_Control_Point_to_mesh[Total_mesh]); // Node_To_Node[K_DIVISION_LENGE][10000]
	// Total_Control_Point_To_Node = (int *)malloc(sizeof(int) * K_DIVISION_LENGE);
	static vector<int> Node_To_Node(K_DIVISION_LENGE * info->Total_Control_Point_to_mesh[Total_mesh]);
	static vector<int> Total_Control_Point_To_Node(K_DIVISION_LENGE);

	if (mode_select == 0)
	{
		// 初期化
		for (int i = 0; i < K_Whole_Size + 1; i++)
			info->K_Whole_Ptr[i] = 0;

		// 大きく分割するためのループ
		for (int N = 0; N < info->Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
		{
			// 各節点に接する節点を取得
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				Total_Control_Point_To_Node[i] = 0;
			}
			for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (int ii = 0; ii < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; ii++)
				{
					int NE = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + ii] - N;
					if (0 <= NE && NE < K_DIVISION_LENGE)
					{
						// ローカル要素
						for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
						{
							// 数字がない時
							if (Total_Control_Point_To_Node[NE] == 0)
							{
								// 節点番号を取得
								Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								Total_Control_Point_To_Node[NE]++;
							}
							// 同じものがあったら, k > 0 以降の取得, kのカウント
							int k;
							for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
							{
								if (Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j])
									break;
							}
							// 未設定のNode_To_Node取得
							if (k == Total_Control_Point_To_Node[NE])
							{
								Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								Total_Control_Point_To_Node[NE]++;
							}
						}
						// 別メッシュとの重なりを考慮
						if (info->eoi[i].size() > 0)
						{
							for (size_t jj = 0; jj < info->eoi[i].size(); jj++)
							{
								// ローカル要素
								for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][jj]]]; j++)
								{
									// 数字がない時
									if (Total_Control_Point_To_Node[NE] == 0)
									{
										// 節点番号を取得
										Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j];
										Total_Control_Point_To_Node[NE]++;
									}

									// 同じものがあったら, k > 0 以降の取得, kのカウント
									int k;
									for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
									{
										if (Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j])
											break;
									}
									// 未設定のNode_To_Node取得
									if (k == Total_Control_Point_To_Node[NE])
									{
										Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j];
										Total_Control_Point_To_Node[NE]++;
									}
								}
							}
						}
					}
				}
			}

			// 順番に並び替える
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
				{
					for (int j = 0; j < Total_Control_Point_To_Node[i]; j++)
					{
						int Min = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
						for (int k = j; k < Total_Control_Point_To_Node[i]; k++)
						{
							if (Min > Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k])
							{
								Min = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k];
								No = k;
							}
						}
						for (int k = No; k > j; k--)
						{
							Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k] = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k - 1];
						}
						Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
					}
				}
			}

			// 節点からcol ptrを求める
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				for (int ii = 0; ii < info->DIMENSION; ii++)
				{
					if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
					{
						int i_index = info->Index_Dof[(N + i) * info->DIMENSION + ii];
						if (i_index >= 0)
						{
							info->K_Whole_Ptr[i_index + 1] = info->K_Whole_Ptr[i_index];
							for (int j = 0; j < Total_Control_Point_To_Node[i]; j++)
							{
								for (int jj = 0; jj < info->DIMENSION; jj++)
								{
									int j_index = info->Index_Dof[Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] * info->DIMENSION + jj];
									if (j_index >= 0 && j_index >= i_index)
										info->K_Whole_Ptr[i_index + 1]++;
								}
							}
						}
					}
				}
			}
		}
	}
	else if (mode_select == 1)
	{
		// 大きく分割するためのループ
		for (int N = 0; N < info->Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
		{
			// 各節点に接する節点を取得
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				Total_Control_Point_To_Node[i] = 0;
			}
			for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (int ii = 0; ii < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; ii++)
				{
					int NE = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + ii] - N;
					if (0 <= NE && NE < K_DIVISION_LENGE)
					{
						// ローカル要素
						for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
						{
							// 数字がない時
							if (Total_Control_Point_To_Node[NE] == 0)
							{
								// 節点番号を取得
								Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								Total_Control_Point_To_Node[NE]++;
							}
							// 同じものがあったら, k > 0 以降の取得, kのカウント
							int k;
							for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
							{
								if (Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j])
									break;
							}
							// 未設定のNode_To_Node取得
							if (k == Total_Control_Point_To_Node[NE])
							{
								Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								Total_Control_Point_To_Node[NE]++;
							}
						}
						// 別メッシュとの重なりを考慮
						if (info->eoi[i].size() > 0)
						{
							for (size_t jj = 0; jj < info->eoi[i].size(); jj++)
							{
								// ローカル要素
								for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][jj]]]; j++)
								{
									// 数字がない時
									if (Total_Control_Point_To_Node[NE] == 0)
									{
										// 節点番号を取得
										Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j];
										Total_Control_Point_To_Node[NE]++;
									}

									// 同じものがあったら, k > 0 以降の取得, kのカウント
									int k;
									for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
									{
										if (Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j])
										{
											break;
										}
									}
									// 未設定のNode_To_Node取得
									if (k == Total_Control_Point_To_Node[NE])
									{
										Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[info->eoi[i][jj] * MAX_NO_CP_ON_ELEMENT + j];
										Total_Control_Point_To_Node[NE]++;
									}
								}
							}
						}
					}
				}
			}

			// 順番に並び替える
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
				{
					for (int j = 0; j < Total_Control_Point_To_Node[i]; j++)
					{
						int Min = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
						for (int k = j; k < Total_Control_Point_To_Node[i]; k++)
						{
							if (Min > Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k])
							{
								Min = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k];
								No = k;
							}
						}
						for (int k = No; k > j; k--)
						{
							Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k] = Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k - 1];
						}
						Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
					}
				}
			}

			// 節点からcol ptrを求める
			for (int i = 0; i < K_DIVISION_LENGE; i++)
			{
				for (int ii = 0; ii < info->DIMENSION; ii++)
				{
					if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
					{
						int i_index = info->Index_Dof[(N + i) * info->DIMENSION + ii];
						int k = 0;
						if (i_index >= 0)
						{
							for (int j = 0; j < Total_Control_Point_To_Node[i]; j++)
							{
								for (int jj = 0; jj < info->DIMENSION; jj++)
								{
									int j_index = info->Index_Dof[Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] * info->DIMENSION + jj];
									if (j_index >= 0 && j_index >= i_index)
									{
										info->K_Whole_Col[info->K_Whole_Ptr[i_index] + k] = j_index;
										k++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
#endif


// K matrix の値を求める
void Make_K_Whole_Val(information *info)
{
	gp_switch(false, info);
	#pragma omp parallel for schedule(dynamic)
	for (int re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		thread_local double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		int i = info->real_element[re];

		// 各要素のK_ELを求める
		Make_K_EL(i, K_EL, info);

		// Valを求める
		for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
			for (int j2 = 0; j2 < info->DIMENSION; j2++)
			{
				long row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
				if (row >= 0)
					for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
						for (int k2 = 0; k2 < info->DIMENSION; k2++)
						{
							long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
							if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
							{
								if (col >= row)
								{
									int id = RowCol_to_icount(row, col, info);
									double val = K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->K_Whole_Val[id] += val;
								}
							}
							else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
								if (col >= 0)
								{
									long id = (long)row * K_Whole_Size + col;
									double val = K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->full_K[id] += val;
								}
						}
			}
	}

	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		thread_local double *coupled_K_EL_1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		int i = info->real_element[re];

		// ローカルメッシュ上の要素について, 重なっている要素が存在するとき
		if (info->c.BASE_PATCH_SIGA == 0)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
			{
				for (size_t j = 0; j < info->eoi[i].size(); j++)
				{
					// 各要素のcoupled_K_ELを求める
					Make_coupled_K_EL_locbase(i, info->eoi[i][j], coupled_K_EL_1, info);

					// Valを求める
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
										{
											if (col >= row)
											{
												int id = RowCol_to_icount(row, col, info);
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->K_Whole_Val[id] += val;
											}
										}
										else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
										{
											if (col >= 0)
											{
												long id_1 = (long)row * K_Whole_Size + col;
												long id_2 = (long)col * K_Whole_Size + row;
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->full_K[id_1] += val;
												#pragma omp atomic
												info->full_K[id_2] += val;
											}
										}
									}
						}
				}
			}
		}
		else if (info->c.BASE_PATCH_SIGA == 1)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[i] == 0 && info->eoi[i].size() > 0))
			{
				for (size_t j = 0; j < info->eoi[i].size(); j++)
				{
					// 各要素のcoupled_K_ELを求める
					Make_coupled_K_EL_globase(info->eoi[i][j], i, coupled_K_EL_1, info);

					// Valを求める
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
										{
											if (col >= row)
											{
												int id = RowCol_to_icount(row, col, info);
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->K_Whole_Val[id] += val;
												if (id == -1)
												{
													cout << "error here\n" << endl;
													exit(1);
												}
											}
										}
										else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
										{
											if (col >= 0)
											{
												long id_1 = (long)row * K_Whole_Size + col;
												long id_2 = (long)col * K_Whole_Size + row;
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->full_K[id_1] += val;
												#pragma omp atomic
												info->full_K[id_2] += val;
											}
										}
									}
						}
				}
			}
		}
	}
}


void Make_K_Whole_Val_octree(information *info)
{
	gp_switch(false, info);
	#pragma omp parallel for schedule(dynamic)
	for (int re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		thread_local double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		int i = info->real_element[re];

		// 各要素のK_ELを求める
		Make_K_EL(i, K_EL, info);

		// Valを求める
		for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
			for (int j2 = 0; j2 < info->DIMENSION; j2++)
			{
				long row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
				if (row >= 0)
					for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
						for (int k2 = 0; k2 < info->DIMENSION; k2++)
						{
							long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
							if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
							{
								if (col >= row)
								{
									int id = RowCol_to_icount(row, col, info);
									double val = K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->K_Whole_Val[id] += val;
								}
							}
							else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
								if (col >= 0)
								{
									long id = (long)row * K_Whole_Size + col;
									double val = K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->full_K[id] += val;
								}
						}
			}
	}

	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		thread_local double *coupled_K_EL_1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		int i = info->real_element[re];

		// ローカルメッシュ上の要素について, 重なっている要素が存在するとき
		if (info->c.BASE_PATCH_SIGA == 0)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
			{
				for (size_t j = 0; j < info->octree_subcell[i].size(); j++)
				{
					for (size_t k = 0; k < info->octree_subcell[i][j].overlapping_ele.size(); k++)
					{
						// 各要素のcoupled_K_ELを求める
						Make_coupled_K_EL_octree_locbase(info->octree_subcell[i][j], info->octree_subcell[i][j].overlapping_ele[k], coupled_K_EL_1, info);

						// Valを求める
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->octree_subcell[i][j].overlapping_ele[k]]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[info->octree_subcell[i][j].overlapping_ele[k] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
											{
												if (col >= row)
												{
													int id = RowCol_to_icount(row, col, info);
													double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

													#pragma omp atomic
													info->K_Whole_Val[id] += val;
												}
											}
											else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
											{
												if (col >= 0)
												{
													long id_1 = (long)row * K_Whole_Size + col;
													long id_2 = (long)col * K_Whole_Size + row;
													double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

													#pragma omp atomic
													info->full_K[id_1] += val;
													#pragma omp atomic
													info->full_K[id_2] += val;
												}
											}
										}
							}
					}
				}
			}
		}
		else if (info->c.BASE_PATCH_SIGA == 1)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[i] == 0 && info->eoi[i].size() > 0))
			{
				for (size_t j = 0; j < info->eoi[i].size(); j++)
				{
					// 各要素のcoupled_K_ELを求める
					Make_coupled_K_EL_globase(info->eoi[i][j], i, coupled_K_EL_1, info);

					// Valを求める
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (info->c.SOLVER == 0 || info->c.SOLVER == 1 || info->c.SOLVER == 4)
										{
											if (col >= row)
											{
												int id = RowCol_to_icount(row, col, info);
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->K_Whole_Val[id] += val;
												if (id == -1)
												{
													cout << "error here\n" << endl;
													exit(1);
												}
											}
										}
										else if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
										{
											if (col >= 0)
											{
												long id_1 = (long)row * K_Whole_Size + col;
												long id_2 = (long)col * K_Whole_Size + row;
												double val = coupled_K_EL_1[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

												#pragma omp atomic
												info->full_K[id_1] += val;
												#pragma omp atomic
												info->full_K[id_2] += val;
											}
										}
									}
						}
				}
			}
		}
	}
}


// 要素合成マトリックス
void Make_K_EL(int El_No, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]] * info->DIMENSION;

	thread_local double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	thread_local double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No].n(); i++)
	{
		double J = Make_B_Matrix_anypoint(El_No, B, &info->gp[El_No].para()[i * info->DIMENSION], info);

		BDBJ(KIEL_SIZE, B, J, K1, info);
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


// 結合要素剛性マトリックス, ローカル要素ベース積分
void Make_coupled_K_EL_locbase(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	thread_local double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	thread_local double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	thread_local double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{
		// check target element
		if (!(info->gp[El_No_loc].isOverlay()[i] && info->gp[El_No_loc].opp_ele()[i] == El_No_glo))
			continue;

		double J = Make_B_Matrix_anypoint(El_No_loc, B, &info->gp[El_No_loc].para()[i * info->DIMENSION], info);
		Make_B_Matrix_anypoint(El_No_glo, BG, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], info);

		// BGTDBLJの計算
		coupled_BDBJ(KIEL_SIZE_glo, KIEL_SIZE_loc, B, BG, J, K1, info);
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No_loc].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


// 結合要素剛性マトリックス, ローカル要素ベース積分
void Make_coupled_K_EL_octree_locbase(subcell &subcell_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int BDBJ_flag = 0;
	int El_No_loc = subcell_loc.e;
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	thread_local double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	thread_local double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	thread_local double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int k = 0; k < info->gp[El_No_loc].n(); k++)
	{
		// J, B, BG の作成
		double para[MAX_DIMENSION];
		double coord[MAX_DIMENSION];
		double out_para[MAX_DIMENSION];
		double G_Gxi[MAX_DIMENSION];

		subcell_loc.get_physical_coord(para, &info->gp[El_No_loc].para()[k * info->DIMENSION], coord, info);
		double J = Make_octree_B_Matrix_anypoint(subcell_loc.subcell_jac_1D, El_No_loc, B, para, info);

		int itr_n = 0, g_p = 0;
		for (int n = 0; n < info->Total_Patch_on_mesh[0]; n++)
		{
			itr_n = calc_patch_parameter_coord(coord, n, out_para, info);
			if (itr_n != ERROR)
			{
				g_p = n;
				break;
			}
		}

		// 要素内外判定
		if (info->DIMENSION == 2)
		{
			if (out_para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 0] + info->Order[g_p * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]] &&
				out_para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 0] + info->Order[g_p * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] &&
				out_para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 1] + info->Order[g_p * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]] &&
				out_para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 1] + info->Order[g_p * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1])
			{
				BDBJ_flag = 1;

				// 親要素座標の算出
				tilde_coord(G_Gxi, out_para, g_p, El_No_glo, info);
			}
			else
			{
				BDBJ_flag = 0;
			}
		}
		else if (info->DIMENSION == 3)
		{
			if (out_para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 0] + info->Order[g_p * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]] &&
				out_para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 0] + info->Order[g_p * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] &&
				out_para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 1] + info->Order[g_p * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]] &&
				out_para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 1] + info->Order[g_p * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1] &&
				out_para[2] >= info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 2] + info->Order[g_p * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2]] &&
				out_para[2] <  info->Position_Knots[info->Total_Knot_to_patch_dim[g_p * info->DIMENSION + 2] + info->Order[g_p * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2] + 1])
			{
				BDBJ_flag = 1;

				// 親要素座標の算出
				tilde_coord(G_Gxi, out_para, g_p, El_No_glo, info);
			}
			else
			{
				BDBJ_flag = 0;
			}
		}

		// 要素内であるとき, 結合要素剛性マトリックス計算
		if (BDBJ_flag)
		{
			// 重なるグローバル要素のBマトリックス
			Make_B_Matrix_anypoint(El_No_glo, BG, G_Gxi, info);

			// BGTDBLJの計算
			coupled_BDBJ(KIEL_SIZE_glo, KIEL_SIZE_loc, B, BG, J, K1, info);
			for (int l = 0; l < KIEL_SIZE_glo; l++)
				for (int m = 0; m < KIEL_SIZE_loc; m++)
					coupled_K_EL[l * MAX_KIEL_SIZE + m] += info->gp[El_No_loc].w()[k] * K1[l * MAX_KIEL_SIZE + m];
		}
	}
}


// 結合要素剛性マトリックス, グローバル要素ベース積分
void Make_coupled_K_EL_globase(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	cout << "Not yet implemented in func 'Make_coupled_K_EL_globase'." << endl;
	exit(1);

	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	static double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	static double *BL = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	static double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_glo].n(); i++)
	{
		if (!(info->gp[El_No_glo].isOverlay()[i] && (info->gp[El_No_glo].opp_ele()[i] == El_No_loc)))
			continue;

		double J = Make_B_Matrix_anypoint(El_No_glo, B, &info->gp[El_No_glo].para()[i * info->DIMENSION], info);
		Make_B_Matrix_anypoint(El_No_loc, BL, &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION], info);

		// BGTDBLJの計算
		coupled_BDBJ(KIEL_SIZE_glo, KIEL_SIZE_loc, BL, B, J, K1, info);
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No_glo].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


// BGマトリックスを求める
double Make_B_Matrix_anypoint(int El_No, double *B, double *Local_coord, information *info)
{
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];
	thread_local double *b = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);

	int i, j, k;

	thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	thread_local vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	shape_and_dshape(R.data(), dR.data(), Local_coord, El_No, true, info);

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_2x2[i][j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k][i] * dR[j * info->DIMENSION + k];
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
					a_3x3[i][j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k][i] * dR[j * info->DIMENSION + k];
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


double Make_octree_B_Matrix_anypoint(double one_over_sub_1D, int El_No, double *B, double *Local_coord, information *info)
{
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];
	thread_local double *b = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);

	int i, j, k;

	thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	thread_local vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	shape_and_dshape(R.data(), dR.data(), Local_coord, El_No, true, info);

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_2x2[i][j] += dR[k * info->DIMENSION + j] * one_over_sub_1D * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k][i] * dR[j * info->DIMENSION + k] * one_over_sub_1D;
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
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_3x3[i][j] += dR[k * info->DIMENSION + j] * one_over_sub_1D * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k][i] * dR[j * info->DIMENSION + k] * one_over_sub_1D;
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
	}

	return J;
}


void BDBJ(int KIEL_SIZE, double *B, double J, double *K_EL, information *info)
{
	int i, j, k;
	thread_local double *BD = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (i = 0; i < KIEL_SIZE * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	// [B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
		for (j = 0; j < D_MATRIX_SIZE; j++)
			for (k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * info->D[k * D_MATRIX_SIZE + j];

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
	}
}


void coupled_BDBJ(int KIEL_SIZE_glo, int KIEL_SIZE_loc, double *B, double *BG, double J, double *K_EL, information *info)
{
	int i, j, k;
	thread_local double *BD = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (i = 0; i < KIEL_SIZE_glo * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	//[B]GT[D][B]LJの計算
	for (i = 0; i < KIEL_SIZE_glo; i++)
		for (j = 0; j < D_MATRIX_SIZE; j++)
			for (k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * info->D[k * D_MATRIX_SIZE + j];

	for (i = 0; i < KIEL_SIZE_glo; i++)
		for (j = 0; j < KIEL_SIZE_loc; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
}


double Make_Jac_anypoint(int El_No, double *Local_coord, information *info)
{
	int i, j, k;
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];

	thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	thread_local vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	shape_and_dshape(R.data(), dR.data(), Local_coord, El_No, true, info);

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
					a_2x2[i][j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
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
					a_3x3[i][j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}
		}

		J = InverseMatrix_3x3(a_3x3);
	}

	return J;
}


// F vector
void Make_F_Vec(information *info)
{
	int i, index;

	for (i = 0; i < info->Total_Load_to_mesh[Total_mesh]; i++)
	{
		index = info->Index_Dof[info->Load_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Load_Node_Dir[i * 2 + 1]];
		if (index >= 0)
			info->rhs_vec[index] += info->Value_of_Load[i];
	}
}


// 強制変位
void Make_F_Vec_disp_const(information *info)
{
	gp_switch(false, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->real_Total_Element_to_mesh[Total_mesh]; i++)
	{
		thread_local double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		int e = info->real_element[i];

		int disp_const_counter = 0;
		for (int j = 0; j < info->DIMENSION; j++)
		{
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
			{
				int id = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + j];
				if (id < 0)
					disp_const_counter++;
			}
		}

		if (disp_const_counter > 0)
		{
			Make_K_EL(e, K_EL, info);
			for (int j = 0; j < info->DIMENSION; j++)
			{
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int id_row = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + j;
					int row = info->Index_Dof[id_row];
					if (row >= 0)
					{
						int row_local = k * info->DIMENSION + j;
						for (int l = 0; l < info->DIMENSION; l++)
						{
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
							{
								int id_col = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l;
								int col = info->Index_Dof[id_col];
								if (col < 0)
								{
									int col_local = m * info->DIMENSION + l;
									for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
									{
										if (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] == info->Constraint_Node_Dir[n * 2 + 0] && l == info->Constraint_Node_Dir[n * 2 + 1])
										{
											double val = -K_EL[row_local * MAX_KIEL_SIZE + col_local] * info->Value_of_Constraint[n];

											#pragma omp atomic
											info->rhs_vec[row] += val;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->real_Total_Element_to_mesh[Total_mesh]; i++)
	{
		thread_local double *coupled_K_EL_1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		int e = info->real_element[i];

		// for S-IGA
		if (info->c.BASE_PATCH_SIGA == 0)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[e] > 0 && info->eoi[e].size() > 0))
				for (size_t j = 0; j < info->eoi[e].size(); j++)
				{
					int loc_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
						{
							int loc_id = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (loc_id < 0)
								loc_disp_const_counter++;
						}
					}

					int glo_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; l++)
						{
							int glo_id = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (glo_id < 0)
								glo_disp_const_counter++;
						}
					}

					if (loc_disp_const_counter > 0)
					{
						Make_coupled_K_EL_locbase(e, info->eoi[e][j], coupled_K_EL_1, info);
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														double val = -coupled_K_EL_1[row_local * MAX_KIEL_SIZE + col_local] * info->Value_of_Constraint[n];

														#pragma omp atomic
														info->rhs_vec[row] += val;
													}
												}
											}
										}
								}
							}
					}

					if (glo_disp_const_counter > 0)
					{
						Make_coupled_K_EL_locbase(e, info->eoi[e][j], coupled_K_EL_1, info);
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														double val = -coupled_K_EL_1[col_local * MAX_KIEL_SIZE + row_local] * info->Value_of_Constraint[n];

														#pragma omp atomic
														info->rhs_vec[row] += val;
													}
												}
											}
										}
								}
							}
					}
				}
		}
		else if (info->c.BASE_PATCH_SIGA == 1)
		{
			if (Total_mesh >= 2 && (info->Element_mesh[e] == 0 && info->eoi[e].size() > 0))
				for (size_t j = 0; j < info->eoi[e].size(); j++)
				{
					int loc_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; l++)
						{
							int loc_id = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (loc_id < 0)
								loc_disp_const_counter++;
						}
					}

					int glo_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
						{
							int glo_id = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (glo_id < 0)
								glo_disp_const_counter++;
						}
					}

					if (loc_disp_const_counter > 0)
					{
						Make_coupled_K_EL_globase(info->eoi[e][j], e, coupled_K_EL_1, info);
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														double val = -coupled_K_EL_1[row_local * MAX_KIEL_SIZE + col_local] * info->Value_of_Constraint[n];

														#pragma omp atomic
														info->rhs_vec[row] += val;
													}
												}
											}
										}
								}
							}
					}

					if (glo_disp_const_counter > 0)
					{
						Make_coupled_K_EL_globase(info->eoi[e][j], e, coupled_K_EL_1, info);
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[e][j]]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[e][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														double val = -coupled_K_EL_1[col_local * MAX_KIEL_SIZE + row_local] * info->Value_of_Constraint[n];

														#pragma omp atomic
														info->rhs_vec[row] += val;
													}
												}
											}
										}
								}
							}
					}
				}
		}
	}
}


// 分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Force_to_F_Vec(information *info)
{
	int i, j, index;
	for (j = 0; j < info->DIMENSION; j++)
	{
		for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		{
			index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
				info->rhs_vec[index] += info->Equivalent_Nodal_Force[i * info->DIMENSION + j];
		}
	}
}