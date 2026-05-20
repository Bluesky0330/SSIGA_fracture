// header
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <any>
#include "_class.hpp"
#include "_func.hpp"
#include "_sub.hpp"

void constant::get_data(const char *s, information *info)
{
	ifstream file(s);
	if (!file)
	{
		cerr << "Cannot open file." << endl;
		exit(1);
	}

	string line;
	while (getline(file, line))
	{
		if (line.empty() || line[0] == '#')
			continue;

		size_t commentPos = line.find('#');
		if (commentPos != string::npos)
			line = line.substr(0, commentPos);

		istringstream iss(line);
		string variableName;
		string valueStr;

		if (iss >> variableName >> valueStr)
		{
			cout << variableName << ": " << valueStr << endl;
			// 値に小数点が含まれているかチェック
			if (valueStr.find('.') != string::npos)
			{
				// 浮動小数点数として処理
				double value = stod(valueStr);
				info->c.data_vec.emplace_back(value);
			}
			else
			{
				// 整数として処理
				int value = stoi(valueStr);
				info->c.data_vec.emplace_back(value);
			}
		}
	}
}


void target_domain::setPara(int octree_level, int *octree_position, information *info)
{
	int division = pow_int(2, octree_level);

	int p = info->Element_patch[e];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int offset = info->ENC[e * info->DIMENSION + i];
		double s = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + i] + info->Order[p * info->DIMENSION + i] + offset];
		double e = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + i] + info->Order[p * info->DIMENSION + i] + offset + 1];
		double id = (double)octree_position[i];
		para_start[i] = (s + (e - s) * id / (double)division);
		para_end[i] = (s + (e - s) * (id + 1.0) / (double)division);
	}
}


octree_node::octree_node(int *input_octree_position, int input_level, information *info) : is_leaf(true), level(input_level)
{
	for (int i = 0; i < info->DIMENSION; i++)
		octree_position[i] = input_octree_position[i];
}


void octree_node::recieve_para(target_domain &input_td, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		para_start[i] = input_td.para_start[i];
		para_end[i] = input_td.para_end[i];
	}
}


subcell::subcell(int input_e, octree_node &input_octree_node, information *info) : e(input_e), patch(info->Element_patch[input_e]), level(input_octree_node.level)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		para_start[i] = input_octree_node.para_start[i];
		para_end[i] = input_octree_node.para_end[i];
	}
	for (size_t i = 0; i < input_octree_node.overlapping_ele.size(); i++)
		overlapping_ele.emplace_back(input_octree_node.overlapping_ele[i]);

	subcell_jac_1D = 1.0 / static_cast<double>(pow_int(2, level));
}


void subcell::get_tilde_para(double *tilde_para, double *child_tilde_para, information *info)
{
	double para[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		para[i] = para_start[i] + 0.5 * (child_tilde_para[i] + 1.0) * (para_end[i] - para_start[i]);
	tilde_coord(tilde_para, para, patch, e, info);
}


void subcell::get_physical_coord(double *tilde_para, double *child_tilde_para, double *coord, information *info)
{
	get_tilde_para(tilde_para, child_tilde_para, info);

	for (int i = 0; i < info->DIMENSION; i++)
		coord[i] = 0.0;

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	shape_and_dshape(R.data(), tilde_para, e, info);
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
			coord[j] += R[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	}
}


void gauss_points::set(int input_patch, int input_ele, int input_n, bool opp_flag, information *info)
{
	patch = input_patch;
	ele = input_ele;
	n = input_n;

	int gp_n = pow_int(info->c.NG, info->DIMENSION);

	if (input_n == gp_n)
	{
		para = info->gauss_point;
		w = info->gauss_w;
	}
	else
	{
		para = info->gauss_point_ex;
		w = info->gauss_w_ex;
	}

	coord.resize(n * info->DIMENSION);

	makeCoord(info);

	strain_energy_density.resize(n, 0.0);
	delta_strain_energy_density.resize(n, 0.0);

	// overlap variables
	isOverlay.resize(n, false);
	// if (opp_flag)
	{
		opp_patch.resize(n);
		opp_ele.resize(n);
		opp_para.resize(n * info->DIMENSION, 0.0);
		opp_para_tilde.resize(n * info->DIMENSION, 0.0);

		makeOpp(info);
	}

	// nonlinear variables
	deformation_gradient.resize(n * info->DIMENSION * info->DIMENSION, 0.0);
	current_deformation_gradient.resize(n * info->DIMENSION * info->DIMENSION, 0.0);
	elastic_strain.resize(n * D_MATRIX_SIZE, 0.0);
	elastic_strain_trial.resize(n * D_MATRIX_SIZE);
	current_elastic_strain.resize(n * D_MATRIX_SIZE);
	back_stress.resize(n * D_MATRIX_SIZE, 0.0);
	current_back_stress.resize(n * D_MATRIX_SIZE);
	stress.resize(n * D_MATRIX_SIZE, 0.0);
	current_stress.resize(n * D_MATRIX_SIZE);
	equivalent_stress.resize(n);
	yield_stress.resize(n);
	current_yield_stress.resize(n);
	equivalent_plastic_strain.resize(n, 0.0);
	equivalent_plastic_strain_increment.resize(n, 0.0);

	previous_elastic_strain.resize(n * D_MATRIX_SIZE, 0.0);
	previous_stress.resize(n * D_MATRIX_SIZE, 0.0);
}


void gauss_points::reset(information *info)
{
	// recalculate
	makeCoord(info);

	// overlap variables
	makeOpp(info);
}


void gauss_points::makeCoord(information *info)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
			coord[i * info->DIMENSION + j] = 0.0;

		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), &para[i * info->DIMENSION], ele, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
				coord[i * info->DIMENSION + k] += R[j] * info->Node_Coordinate[info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
	}
}


void gauss_points::makeOpp(information *info)
{
	int start = 0, end = 0;
	if (patch >= info->Total_Patch_to_mesh[1])
	{
		start = 0, end = info->Total_Patch_to_mesh[1];
	}
	else
	{
		start = info->Total_Patch_to_mesh[1], end = info->Total_Patch_to_mesh[Total_mesh];
	}

	for (int i = 0; i < n; i++)
	{
		isOverlay[i] = false;
		double temp_para[MAX_DIMENSION] = {0.0};

		int itr_n = 0;
		for (int j = start; j < end; j++)
		{
			itr_n = calc_patch_parameter_coord(&coord[i * info->DIMENSION], j, temp_para, info);
			if (itr_n != ERROR)
			{
				isOverlay[i] = true;
				opp_patch[i] = j;
				for (int k = 0; k < info->DIMENSION; k++)
					opp_para[i * info->DIMENSION + k] = temp_para[k];
				break;
			}
		}

		if (isOverlay[i])
		{
			opp_ele[i] = ele_check(opp_patch[i], &opp_para[i * info->DIMENSION], info);
			tilde_coord(&opp_para_tilde[i * info->DIMENSION], &opp_para[i * info->DIMENSION], opp_patch[i], opp_ele[i], info);
		}
	}
}


void gauss_points::setBlendingFunction(information *info)
{
	if (info->aux != nullptr && ele < info->Total_Element_to_mesh[1])
		for (int i = 0; i < n; i++)
		{
			double bf = std::min(blendingFunction(ele, patch, &para[i * info->DIMENSION], info), blendingFunction(ele, patch, &para[i * info->aux->DIMENSION], info->aux.get()));
			blending_function.emplace_back(bf);
			if (bf < MERGE_ERROR)
				isSkip.emplace_back(true);
			else
				isSkip.emplace_back(false);
		}
	else if (info->aux == nullptr && ele < info->Total_Element_to_mesh[1])
	{
		constexpr double bf = ERROR;
		for (int i = 0; i < n; i++)
		{
			blending_function.emplace_back(bf);
			isSkip.emplace_back(true);
		}
	}
	else
		for (int i = 0; i < n; i++)
		{
			double bf = blendingFunction(ele, patch, &para[i * info->DIMENSION], info);
			blending_function.emplace_back(bf);
			if (bf < MERGE_ERROR)
				isSkip.emplace_back(true);
			else
				isSkip.emplace_back(false);
		}
}


J_point_info::J_point_info(information *info, int patch, int dir1, int dir2, int dir3, int dir4) : p(patch), crack_dir(dir1), r_dir(dir2), tangent_dir{dir3, dir4}
{
	// init
	for (int i = 0; i < info->DIMENSION; i++)
		normal.emplace_back(0.0);

	// set ele_in_normal_patch
	int current_e = 0;
	for (int i = 0; i < p; i++)
	{
		int temp_e = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			temp_e *= info->No_Control_point[i * info->DIMENSION + j] - info->Order[i * info->DIMENSION + j];
		current_e += temp_e;
	}
	int ele_n = 1;
	for (int i = 0; i < info->DIMENSION; i++)
		ele_n *= info->No_Control_point[p * info->DIMENSION + i] - info->Order[p * info->DIMENSION + i];

	for (int i = 0; i < ele_n; i++)
		ele_in_normal_patch.emplace_back(i + current_e);
}


void J_point_info::init(J_info &J, information *info)
{
	int size = J.e_list.size() * J.sub_n;

	coord.clear();
	normal.clear();
	J_val.clear();
	J_integral_area.clear();
	delta_arc_length.clear();
	delta_arc_length_at_half_para.clear();
	coord.resize(size, vector<double>(info->DIMENSION, 0.0));
	normal.resize(size, vector<double>(info->DIMENSION, 0.0));
	J_val.resize(size, 0.0);
	J_integral_area.resize(size, 0.0);
	delta_arc_length.resize(size, 0.0);
	delta_arc_length_at_half_para.resize(size, 0.0);

	// interaction integral method
	IIM_K_val.clear();
	IIM_K_val.resize(size, vector<double>(3, 0.0)); // 3 MODES
	Q.clear();
	Q.resize(size, vector<vector<double>>(info->DIMENSION, vector<double>(info->DIMENSION, 0.0)));
}


void sub_ele_J::init(information *info)
{
	// allocate internal state variables
	elastic_strain.resize(gp_num, vector<double>(D_MATRIX_SIZE, 0.0));
	stress.resize(gp_num, vector<double>(D_MATRIX_SIZE, 0.0));
	previous_elastic_strain.resize(gp_num, vector<double>(D_MATRIX_SIZE, 0.0));
	previous_stress.resize(gp_num, vector<double>(D_MATRIX_SIZE, 0.0));
	yield_stress.resize(gp_num, get_hardening_stress(0.0, info));
	equivalent_plastic_strain.resize(gp_num, 0.0);
	previous_W.resize(gp_num, 0.0);

	// make para_tilde
	double delta = 1.0 / (static_cast<double>(sub_n * face_n));
	para_tilde.resize(gp_num * info->DIMENSION);
	for (int i = 0; i < gp_num; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
		{
			if (j == crack_dir)
			{
				double normalized_local_para = (gp[i * info->DIMENSION + j] + 1.0) / 2.0;
				double normalized_para = (static_cast<double>(slice * face_n + face) + normalized_local_para) * delta;
				para_tilde[i * info->DIMENSION + j] = normalized_para * 2.0 - 1.0;
			}
			else
				para_tilde[i * info->DIMENSION + j] = gp[i * info->DIMENSION + j];
		}
	}
	
	coord.resize(gp_num * info->DIMENSION, 0.0);
	for (int i = 0; i < gp_num; i++)
	{
		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		double *para_tilde_ptr = para_tilde.data() + i * info->DIMENSION;
		shape_and_dshape(R.data(), para_tilde_ptr, e, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
				coord[i * info->DIMENSION + k] += R[j] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
	}
	
	// for overlay
	overlay.resize(gp_num, false);
	if (Total_mesh > 1)
	{
		overlay_ele.resize(gp_num, -1);

		opp_para_tilde.resize(gp_num * info->DIMENSION, 0.0);
		for (int i = 0; i < gp_num; i++)
		{
			double *coord_ptr = coord.data() + i * info->DIMENSION;
			double *opp_para_tilde_ptr = opp_para_tilde.data() + i * info->DIMENSION;
			int itr_n = 0, g_p = 0;
			vector<double> out_natural_coord(info->DIMENSION, 0.0);
			for (int j = 0; j < info->Total_Patch_on_mesh[0]; j++)
			{
				itr_n = calc_patch_parameter_coord(coord_ptr, j, out_natural_coord.data(), info);
				if (itr_n != ERROR)
				{
					g_p = j;
					overlay[i] = true;
					overlay_ele[i] = ele_check(g_p, out_natural_coord.data(), info);
					tilde_coord(opp_para_tilde_ptr, out_natural_coord.data(), g_p, overlay_ele[i], info);
					break;
				}
			}
		}
	}
}