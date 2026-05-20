// header
#include "_header.hpp"
#include "_sub.hpp"

using namespace std;

// search overlapping element
void searchOverlappingEle(information *info)
{
	vector<vector<int>> ecn(info->Total_Element_to_mesh[Total_mesh]);  // element connectivity
	vector<vector<int>> &eoi = info->eoi; // element overlapping information
	vector<vector<int>> enoi(info->Total_Element_to_mesh[Total_mesh]); // element non-overlapping information

	int zero_vec[MAX_DIMENSION] = {0};

	// func ptr
	bool (*search)(target_domain &, target_domain &, information *) = (info->DIMENSION == 2) ? search2D : search3D;

	// make element connectivity
	elementConnectivity(info, ecn);

	// set loc as base
	#pragma omp parallel for schedule(dynamic)
	for (int e = info->Total_Element_to_mesh[1]; e < info->Total_Element_to_mesh[Total_mesh]; e++)
	{
		target_domain td_base(e);
		vector<target_domain> td, td_swap; // target domain
		td_base.setPara(0, zero_vec, info);

		// set start element
		if (int start_e = setStartEle(e, info); start_e != ERROR)
			td.emplace_back(start_e);
		else
		{
			printf("Error: start element is not found, local element is not inside of global patch \n");
			exit(1);
		}

		// search overlapping element
		bool break_flag = false, inner_flag = false;
		for (int i = 0; ; i++)
		{
			for (size_t j = 0; j < td.size(); j++)
			{
				td[j].setPara(0, zero_vec, info);

				auto overlapping = search(td_base, td[j], info);
				if (overlapping)
					eoi[e].emplace_back(td[j].e);
				else
					enoi[e].emplace_back(td[j].e);
			}

			// sort eoi and enoi
			sortUniqueErase(eoi[e]);
			sortUniqueErase(enoi[e]);

			// set next element
			td_swap.clear();
			for (size_t j = 0; j < td.size(); j++)
			{
				auto it = std::find(eoi[e].begin(), eoi[e].end(), td[j].e);
				if (it != eoi[e].end())
					for (size_t k = 0; k < ecn[td[j].e].size(); k++)
					{
						int next_e = ecn[td[j].e][k];
						auto it_eoi = std::find(eoi[e].begin(), eoi[e].end(), next_e);
						auto it_enoi = std::find(enoi[e].begin(), enoi[e].end(), next_e);
						if (it_eoi == eoi[e].end() && it_enoi == enoi[e].end())
							td_swap.emplace_back(next_e);
					}
			}

			// first time, if there is no intersection, check if L element is inside G element or outside
			if (i == 0 && eoi[e].size() == 0 && !checkInner(td_base, td[0], info))
			{
				eoi[e].emplace_back(td[0].e);
				break_flag = true;
				inner_flag = true;
			}
			else if (eoi[e].size() == 0)
			{
				// add and sort enoi
				for (size_t j = 0; j < td.size(); j++)
					enoi[e].emplace_back(td[j].e);
				sortUniqueErase(enoi[e]);

				for (size_t j = 0; j < td.size(); j++)
				{
					for (size_t k = 0; k < ecn[td[j].e].size(); k++)
					{
						int next_e = ecn[td[j].e][k];
						auto it_enoi = std::find(enoi[e].begin(), enoi[e].end(), next_e);
						if (it_enoi == enoi[e].end())
							td_swap.emplace_back(next_e);
					}
				}
			}

			if (td_swap.size() == 0)
				break_flag = true;

			// sort and swap td_swap
			sortUniqueErase(td_swap);
			std::swap(td, td_swap);

			if (break_flag)
				break;
		}

		// add inner element
		if (!inner_flag)
			addInnerEle(eoi[e], enoi[e], td_base, info);
	}

	// ここでeoi結果を出力
	// for (int e = info->Total_Element_to_mesh[1]; e < info->Total_Element_to_mesh[Total_mesh]; e++)
	// {
	//     printf("Local element %d overlaps with global elements: ", e);
	//     for (size_t i = 0; i < eoi[e].size(); i++)
	//     {
	//         printf("%d ", eoi[e][i]);
	//     }
	//     printf("(total: %zu elements)\n", eoi[e].size());
	// }

	// make global eoi
	for (int e = info->Total_Element_to_mesh[1]; e < info->Total_Element_to_mesh[Total_mesh]; e++)
		for (size_t i = 0; i < eoi[e].size(); i++)
			eoi[eoi[e][i]].emplace_back(e);
	#pragma omp parallel for
	for (int e = 0; e < info->Total_Element_to_mesh[1]; e++)
		sortUniqueErase(eoi[e]);

	// ここでeoi結果を出力
	// for (int e = 0; e < info->Total_Element_to_mesh[1]; e++)
	// {
	//     if (eoi[e].size() > 0)
	//     {
	//         printf("Global element %d overlaps with local elements: ", e);
	//         for (size_t i = 0; i < eoi[e].size(); i++)
	//         {
	//             printf("%d ", eoi[e][i]);
	//         }
	//         printf("(total: %zu elements)\n", eoi[e].size());
	//     }
	// }

	return;
}


// search start element
int setStartEle(const int e, information *info)
{
	// start coordinate
	double para[MAX_DIMENSION] = {0.0};
	double out_coord[MAX_DIMENSION] = {0.0};
	double out_para[MAX_DIMENSION] = {0.0};
	int global_patch_num = info->Global_local_patch;
	int local_patch_num = info->Element_patch[e];

	// transform to patch parameter coordinate from element parameter coordinate
	double para_patch[MAX_DIMENSION] = {0.0};
	trans_ele_patch_coord(para_patch, para, local_patch_num, e, info);

	// search local element and patch for geometry
	int geo_patch_num = local_patch_num - info->Total_Patch_on_mesh[0];
	int geo_e = geo_ele_check(geo_patch_num, para_patch, info);

	// transform to element parameter coordinate from patch parameter coordinate
	double para_ele[MAX_DIMENSION] = {0.0};
	geo_tilde_coord(para_ele, para_patch, geo_patch_num, geo_e, info);

	// physical coordinate
    geo_parameter_coord(geo_e, para_ele, out_coord, info);

	// search global point
	int itr_n = 0, g_p = 0;
	itr_n = calc_global_patch_parameter_coord(out_coord, global_patch_num, out_para, info);
	if (itr_n != ERROR)
	{
		g_p = info->Global_local_patch;
	}

	if (itr_n == ERROR){
		return ERROR;
	}

	return ele_check(g_p, out_para, info);
}


bool checkInner(target_domain &td_base, target_domain &td, information *info)
{
	// check inner
	double vertex_para[MAX_DIMENSION] = {-1.0, -1.0, -1.0};
	double out_coord[MAX_DIMENSION] = {0.0};
	double out_para[MAX_DIMENSION] = {0.0};
	trans_ele_patch_coord(out_coord, vertex_para, info->Element_patch[td.e], td.e, info);

	// search local point
	bool inner_flag = false;
	int itr_n = 0;
	for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		itr_n = calc_global_patch_parameter_coord(out_coord, i, out_para, info);
		if (itr_n != ERROR)
		{
			if (info->Element_patch[td_base.e] != i)
				return false;
			inner_flag = true;
			break;
		}
	}

	// if vertex locates in local patch
	if (inner_flag)
		if (!checkPosition(td_base, out_para, info))
			inner_flag = false;

	return inner_flag;
}


bool checkInnerCenter(target_domain &td_base, int e_center, information *info)
{
	// check inner
	double center_para[MAX_DIMENSION] = {0.0, 0.0, 0.0};
	double out_coord[MAX_DIMENSION] = {0.0};
	double out_para[MAX_DIMENSION] = {0.0};
	trans_ele_patch_coord(out_coord, center_para, info->Element_patch[e_center], e_center, info);

	// search local point
	bool inner_flag = false;
	int itr_n = 0;
	for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		itr_n = calc_local_patch_parameter_coord(out_coord, i, out_para, info);
		if (itr_n != ERROR)
		{
			if (info->Element_patch[td_base.e] != i)
				return false;
			inner_flag = true;
			break;
		}
	}

	// if a point locates in local patch
	if (inner_flag)
		if (!checkPosition(td_base, out_para, info))
			inner_flag = false;

	return inner_flag;
}


bool checkInnerCenterEachPatch(target_domain &td_base, target_domain &td, information *info)
{
	// check inner
	double center_para[MAX_DIMENSION] = {0.0};
	for (int i = 0; i < info->DIMENSION; i++)
		center_para[i] = (td.para_start[i] + td.para_end[i]) / 2.0;

	double out_para[MAX_DIMENSION] = {0.0};

	int patch = 0;
	double start_end[2];
	if (info->Element_mesh[td.e] == 0)
	{
		patch = info->Element_patch[td_base.e];
		start_end[0] = info->Total_Patch_to_mesh[1];
		start_end[1] = info->Total_Patch_to_mesh[Total_mesh];
	}
	else
	{
		patch = info->Element_patch[td.e];
		start_end[0] = 0;
		start_end[1] = info->Total_Patch_to_mesh[1];
	}

	// search point
	bool inner_flag = false;
	int itr_n = 0;
	for (int i = start_end[0]; i < start_end[1]; i++)
	{
		itr_n = calc_local_patch_parameter_coord(center_para, i, out_para, info);
		if (itr_n != ERROR)
		{
			if (patch != i)
				return false;
			inner_flag = true;
			break;
		}
	}

	// if a point locates in patch
	if (inner_flag)
		if (!checkPosition(td_base, out_para, info))
			inner_flag = false;
	
	return inner_flag;
}


bool checkPosition(target_domain &td_base, double *check_para, information *info)
{
	// check position
	for (int i = 0; i < info->DIMENSION; i++)
	{
		if (fabs(check_para[i] - td_base.para_start[i]) < MERGE_ERROR || fabs(check_para[i] - td_base.para_end[i]) < MERGE_ERROR || (td_base.para_start[i] < check_para[i] && check_para[i] < td_base.para_end[i]))
			continue;
		else
			return false;
	}
	return true;
}


void addInnerEle(vector<int> &eoi_e, vector<int> &enoi_e, target_domain &td_base, information *info)
{
	vector<max_min> max_min_val;
	for (int i = 0; i < 2; i++)
	{
		vector<int> &temp = (i == 0) ? eoi_e : enoi_e;
		for (size_t j = 0; j < temp.size(); j++)
		{
			int patch_num = info->Element_patch[temp[j]];
			auto it = find_if(max_min_val.begin(), max_min_val.end(), [patch_num](max_min &obj) { return obj.patch_num == patch_num; });
			if (it == max_min_val.end())
			{
				max_min_val.emplace_back(patch_num);
				for (int k = 0; k < info->DIMENSION; k++)
				{
					max_min_val.back().min[k] = info->ENC[temp[j] * info->DIMENSION + k];
					max_min_val.back().max[k] = info->ENC[temp[j] * info->DIMENSION + k];
				}
			}
			else
			{
				for (int k = 0; k < info->DIMENSION; k++)
				{
					if (max_min_val[it - max_min_val.begin()].min[k] > info->ENC[temp[j] * info->DIMENSION + k])
						max_min_val[it - max_min_val.begin()].min[k] = info->ENC[temp[j] * info->DIMENSION + k];
					if (max_min_val[it - max_min_val.begin()].max[k] < info->ENC[temp[j] * info->DIMENSION + k])
						max_min_val[it - max_min_val.begin()].max[k] = info->ENC[temp[j] * info->DIMENSION + k];
				}
			}
		}
	}

	for (size_t i = 0; i < max_min_val.size(); i++)
	{
		int e_to_patch = 0;
		int p = max_min_val[i].patch_num;
		for (int j = 0; j < p; j++)
		{
			int temp_e_to_patch = 1;
			for (int k = 0; k < info->DIMENSION; k++)
				temp_e_to_patch *= (info->No_Control_point[j * info->DIMENSION + k] - info->Order[j * info->DIMENSION + k]);
			e_to_patch += temp_e_to_patch;
		}

		if (info->DIMENSION == 2)
		{
			for (int j = max_min_val[i].min[0]; j <= max_min_val[i].max[0]; j++)
				for (int k = max_min_val[i].min[1]; k <= max_min_val[i].max[1]; k++)
				{
					int e_center = e_to_patch + j + k * (info->No_Control_point[p * info->DIMENSION] - info->Order[p * info->DIMENSION]);
					if (checkInnerCenter(td_base, e_center, info))
						eoi_e.emplace_back(e_center);
				}
		}
		else if (info->DIMENSION == 3)
		{
			for (int j = max_min_val[i].min[0]; j <= max_min_val[i].max[0]; j++)
				for (int k = max_min_val[i].min[1]; k <= max_min_val[i].max[1]; k++)
					for (int l = max_min_val[i].min[2]; l <= max_min_val[i].max[2]; l++)
					{
						int e_center = e_to_patch + j + k * (info->No_Control_point[p * info->DIMENSION] - info->Order[p * info->DIMENSION]) + l * (info->No_Control_point[p * info->DIMENSION] - info->Order[p * info->DIMENSION]) * (info->No_Control_point[p * info->DIMENSION + 1] - info->Order[p * info->DIMENSION + 1]);
						if (checkInnerCenter(td_base, e_center, info))
							eoi_e.emplace_back(e_center);
					}
		}
	}
	sortUniqueErase(eoi_e);
	return;
}


template <typename T>
void sortUniqueErase(vector<T> &vec)
{
	std::sort(vec.begin(), vec.end());
	auto last = std::unique(vec.begin(), vec.end());
	vec.erase(last, vec.end());
}


void physical_coord(const int e, double *para, double *out_coord, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
		out_coord[i] = 0.0;

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	shape_and_dshape(R.data(), para, e, info);
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			out_coord[j] += R[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	return;
}


// ローカル形状表現の写像される、グローバルパラメータ座標
void geo_parameter_coord(const int e, double *para, double *out_coord, information *info)
{   
	for (int i = 0; i < info->DIMENSION; i++)
		out_coord[i] = 0.0;

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	geo_shape_and_dshape(R.data(), para, e, info);
	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e]]; i++)
        for (int j = 0; j < info->DIMENSION; j++)
            out_coord[j] += R[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
    return;
}


void physical_coord_R(const int e, double *R, double *out_coord, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
		out_coord[i] = 0.0;

	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			out_coord[j] += R[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	return;
}


// ローカル形状表現の写像される、グローバルパラメータ座標
void geo_parameter_coord_R(const int e_local, double *R_local, double *out_coord, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
		out_coord[i] = 0.0;

	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e_local]]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			out_coord[j] += R_local[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e_local * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	return;
}


void makeChildPosition(int *child_position, int i, int *parent_position, information *info)
{
	constexpr int div_1D = 2;
	if (info->DIMENSION == 2)
	{
		child_position[0] = parent_position[0] * div_1D + i % div_1D;
		child_position[1] = parent_position[1] * div_1D + i / div_1D;
		return;
	}
	else if (info->DIMENSION == 3)
	{
		constexpr int div_2D = div_1D * div_1D;
		child_position[0] = parent_position[0] * div_1D + i % div_1D;
		child_position[1] = parent_position[1] * div_1D + (i / div_1D) % div_1D;
		child_position[2] = parent_position[2] * div_1D + i / div_2D;
		return;
	}
	return;
}


// search overlapping element
bool search2D(target_domain &td_base, target_domain &td, information *info)
{
	constexpr int max_itr = 8;
	constexpr int edge_n = 4;

	vector<double> R_base(MAX_NO_CP_ON_ELEMENT);
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR_base(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);

	// fix axis loop
	for (int i = 0; i < edge_n; i++)
	{
		for (int j = 0; j < edge_n; j++)
		{
			int target_axis_base = 0;
			int target_axis = 0;
			int patch_base = info->Element_patch[td_base.e] - info->Total_Patch_on_mesh[0]; // ローカル形状表現のパッチ番号
			double para_base[MAX_DIMENSION] = {0.0};
			double para[MAX_DIMENSION] = {0.0};
			double tilde_para_base[MAX_DIMENSION] = {0.0};
			double dR_simple[MAX_DIMENSION] = {0.0};

			setEdge2D(i, &target_axis_base, para_base, td_base); // このエッジはローカル変位表現のものを使用
			setEdge2D(j, &target_axis, para, td);

			// Newton-Raphson loop
			bool singular_flag = true;
			int singular_axis = -1;
			for (int k = 0; k < max_itr; k++)
			{
				double coord_base[MAX_DIMENSION] = {0.0};
				double coord[MAX_DIMENSION] = {0.0};

				vector<double> diff(MAX_DIMENSION);
				vector<double> sol(MAX_DIMENSION);
				vector<double> J(MAX_DIMENSION * MAX_DIMENSION, 0.0);

				if (check_nonfinite(para_base, info->DIMENSION) || check_nonfinite(para, info->DIMENSION))
					break;
				int e_base = geo_ele_check(patch_base, para_base, info); // ローカル形状表現の要素番号を取得
				geo_tilde_coord(tilde_para_base, para_base, patch_base, e_base, info);
				geo_shape_and_dshape(R_base.data(), dR_base.data(), tilde_para_base, e_base, false, info);
				geo_parameter_coord_R(e_base, R_base.data(), coord_base, info);

				// グローバルのパラメータ座標系における微分を1に設定(dξ/dξ = 1)
				std::fill(dR_simple, dR_simple + info->DIMENSION, 0.0);
				dR_simple[target_axis] = 1.0;

				// 非ベースのcoord配列には、グローバルパラメータ座標を直接代入
				std::copy(para, para + info->DIMENSION, coord);

				double r = 0.0;
				for (int l = 0; l < info->DIMENSION; l++)
				{
					diff[l] = coord_base[l] - coord[l];
					r += diff[l] * diff[l];
				}
				r = sqrt(r);

				// check position
				if (!isfinite(r))
					break;
				if (r < MERGE_ERROR)
					return true;

                // make J
                for (int l = 0; l < info->DIMENSION; l++)
                    for (int n = 0; n < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e_base]]; n++)
                        J[l * info->DIMENSION + 0] += dR_base[n * info->DIMENSION + target_axis_base] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e_base * MAX_NO_CP_ON_ELEMENT + n] * (info->DIMENSION + 1) + l];
                for (int l = 0; l < info->DIMENSION; l++)
                    J[l * info->DIMENSION + 1] -= dR_simple[l];

				// solve [J]{sol} = {diff}
				singular_flag = GaussianElimination2(sol.data(), diff.data(), J.data(), info->DIMENSION);

				if (!singular_flag && k == 0)
				{
					for (int l = 0; l < info->DIMENSION; l++)
						if (fabs(coord_base[l] - coord[l]) < MERGE_ERROR)
							singular_axis = (l == 0) ? 1 : 0;
					if (singular_axis == -1)
					{
						singular_flag = true;
						continue;
					}
					break;
				}
				else if (!singular_flag && k != 0)
				{
					singular_flag = true;
					continue;
				}

				// update
				para_base[target_axis_base] -= sol[0];
				para[target_axis] -= sol[1];

				if (para_base[target_axis_base] > td_base.para_end[target_axis_base])
					para_base[target_axis_base] = td_base.para_end[target_axis_base];
				if (para_base[target_axis_base] < td_base.para_start[target_axis_base])
					para_base[target_axis_base] = td_base.para_start[target_axis_base];

				if (para[target_axis] > td.para_end[target_axis])
					para[target_axis] = td.para_end[target_axis];
				if (para[target_axis] < td.para_start[target_axis])
					para[target_axis] = td.para_start[target_axis];
			}

			// reduce variable because there are infinite solutions
			if (!singular_flag)
			{
				setEdge2D(i, &target_axis_base, para_base, td_base);
				setEdge2D(j, &target_axis, para, td);

				double coord_base[MAX_DIMENSION] = {0.0};

				if (check_nonfinite(para_base, info->DIMENSION))
					break;
				int e_base = geo_ele_check(patch_base, para_base, info);
				geo_tilde_coord(tilde_para_base, para_base, patch_base, e_base, info);
				geo_parameter_coord(e_base, tilde_para_base, coord_base, info);

				for (int k = 0; k < max_itr; k++)
				{
					double coord[MAX_DIMENSION] = {0.0};

					double diff;
					double sol;
					double J = 0.0;

					if (check_nonfinite(para, info->DIMENSION))
						break;
					std::fill(dR_simple, dR_simple + info->DIMENSION, 0.0);
					dR_simple[target_axis] = 1.0;
					std::copy(para, para + info->DIMENSION, coord);

					double r = 0.0;
					diff = coord_base[singular_axis] - coord[singular_axis];
					r = fabs(diff);

					// check position
					if (!isfinite(r))
						break;
					if (r < MERGE_ERROR)
						return true;

					// make J
					J -= dR_simple[target_axis];

					// solve [J]{sol} = {diff}
					sol = diff / J;

					// update
					para[target_axis] -= sol;

					if (para[target_axis] > td.para_end[target_axis])
						para[target_axis] = td.para_end[target_axis];
					if (para[target_axis] < td.para_start[target_axis])
						para[target_axis] = td.para_start[target_axis];
				}

				setEdge2D(i, &target_axis_base, para_base, td_base);
				setEdge2D(j, &target_axis, para, td);

				double coord[MAX_DIMENSION] = {0.0};
				if (check_nonfinite(para, info->DIMENSION))
					break;
				std::copy(para, para + info->DIMENSION, coord);

				for (int k = 0; k < max_itr; k++)
				{
					double coord_base[MAX_DIMENSION] = {0.0};

					double diff;
					double sol;
					double J = 0.0;

					if (check_nonfinite(para_base, info->DIMENSION))
						break;
					int e_base = geo_ele_check(patch_base, para_base, info);
					geo_tilde_coord(tilde_para_base, para_base, patch_base, e_base, info);
					geo_shape_and_dshape(R_base.data(), dR_base.data(), tilde_para_base, e_base, false, info);
					geo_parameter_coord_R(e_base, R_base.data(), coord_base, info);

					double r = 0.0;
					diff = coord_base[singular_axis] - coord[singular_axis];
					r = fabs(diff);

					// check position
					if (!isfinite(r))
						break;
					if (r < MERGE_ERROR)
						return true;

					// make J
					for (int n = 0; n < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e_base]]; n++)
						J += dR_base[n * info->DIMENSION + target_axis_base] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e_base * MAX_NO_CP_ON_ELEMENT + n] * (info->DIMENSION + 1) + singular_axis];

					// solve [J]{sol} = {diff}
					sol = diff / J;

					// update
					para_base[target_axis_base] -= sol;

					if (para_base[target_axis_base] > td_base.para_end[target_axis_base])
						para_base[target_axis_base] = td_base.para_end[target_axis_base];
					if (para_base[target_axis_base] < td_base.para_start[target_axis_base])
						para_base[target_axis_base] = td_base.para_start[target_axis_base];
				}
			}
		}
	}
	return false;
}


void setEdge2D(int edge_num, int *target_axis, double *para, target_domain &temp)
{
	// 左辺
	if (edge_num == 0)
	{
		*target_axis = 1;
		para[0] = temp.para_start[0];
		para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
	}
	// 右辺
	else if (edge_num == 1)
	{
		*target_axis = 1;
		para[0] = temp.para_end[0];
		para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
	}
	// 下辺
	else if (edge_num == 2)
	{
		*target_axis = 0;
		para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
		para[1] = temp.para_start[1];
	}
	// 上辺
	else if (edge_num == 3)
	{
		*target_axis = 0;
		para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
		para[1] = temp.para_end[1];
	}
	return;
}


// search overlapping element
bool search3D(target_domain &td_base, target_domain &td, information *info)
{
	constexpr int max_itr = 5;
	constexpr int edge_n = 12;
	constexpr int face_n = 6;

	vector<double> R_base(MAX_NO_CP_ON_ELEMENT);
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR_base(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);

	// fix axis loop (using local edge and global face)
	for (int i = 0; i < edge_n; i++)
	{
		for (int j = 0; j < face_n; j++)
		{
			int target_axis_base = 0;
			int target_axis[2] = {0, 0};
			int patch_base = info->Element_patch[td_base.e] - info->Total_Patch_on_mesh[0]; // ローカル形状表現のパッチ番号
			double para_base[MAX_DIMENSION] = {0.0};
			double para[MAX_DIMENSION] = {0.0};
			double tilde_para_base[MAX_DIMENSION] = {0.0};
			double dR_simple[MAX_DIMENSION * MAX_DIMENSION] = {0.0};

			setEdge3D(0, i, &target_axis_base, para_base, td_base);
			setEdge3D(1, j, target_axis, para, td);

			// Newton-Raphson loop
			bool singular_flag = true;
			for (int k = 0; k < max_itr; k++)
			{
				double coord_base[MAX_DIMENSION] = {0.0};
				double coord[MAX_DIMENSION] = {0.0};

				vector<double> diff(MAX_DIMENSION);
				vector<double> sol(MAX_DIMENSION);
				vector<double> J(MAX_DIMENSION * MAX_DIMENSION, 0.0);

				if (check_nonfinite(para_base, info->DIMENSION) || check_nonfinite(para, info->DIMENSION))
					break;
				int e_base = geo_ele_check(patch_base, para_base, info);
				geo_tilde_coord(tilde_para_base, para_base, patch_base, e_base, info);
				geo_shape_and_dshape(R_base.data(), dR_base.data(), tilde_para_base, e_base, false, info);
				geo_parameter_coord_R(e_base, R_base.data(), coord_base, info);

				std::fill(dR_simple, dR_simple + info->DIMENSION * info->DIMENSION, 0.0);
				for (int m = 0; m < info->DIMENSION - 1; m++)
				    dR_simple[target_axis[m] * info->DIMENSION + 1 + m] = 1.0;
				std::copy(para, para + info->DIMENSION, coord);

				double r = 0.0;
				for (int l = 0; l < info->DIMENSION; l++)
				{
					diff[l] = coord_base[l] - coord[l];
					r += diff[l] * diff[l];
				}
				r = sqrt(r);

				// check position
				if (!isfinite(r))
					break;
				if (r < MERGE_ERROR)
					return true;

                // make J
				for (int l = 0; l < info->DIMENSION; l++)
					for (int n = 0; n < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e_base]]; n++)
						J[l * info->DIMENSION + 0] += dR_base[n * info->DIMENSION + target_axis_base] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e_base * MAX_NO_CP_ON_ELEMENT + n] * (info->DIMENSION + 1) + l];
				for (int l = 0; l < info->DIMENSION; l++)
					for (int m = 0; m < info->DIMENSION; m++)
							J[m * info->DIMENSION + l] -= dR_simple[m * info->DIMENSION + l];

				// solve [J]{sol} = {diff}
				singular_flag = GaussianElimination2(sol.data(), diff.data(), J.data(), info->DIMENSION);

				if (!singular_flag)
					continue;

				// update
				para_base[target_axis_base] -= sol[0];
				para[target_axis[0]] -= sol[1];
				para[target_axis[1]] -= sol[2];

				if (para_base[target_axis_base] > td_base.para_end[target_axis_base])
					para_base[target_axis_base] = td_base.para_end[target_axis_base];
				if (para_base[target_axis_base] < td_base.para_start[target_axis_base])
					para_base[target_axis_base] = td_base.para_start[target_axis_base];

				for (int l = 0; l < 2; l++)
				{
					if (para[target_axis[l]] > td.para_end[target_axis[l]])
						para[target_axis[l]] = td.para_end[target_axis[l]];
					if (para[target_axis[l]] < td.para_start[target_axis[l]])
						para[target_axis[l]] = td.para_start[target_axis[l]];
				}
			}
		}
	}

	// fix axis loop (using global edge and local face)
	for (int i = 0; i < edge_n; i++)
	{
		for (int j = 0; j < face_n; j++)
		{
			int target_axis_base = 0;
			int target_axis[2] = {0, 0};
			int patch = info->Element_patch[td_base.e] - info->Total_Patch_on_mesh[0]; // ローカル形状表現のパッチ番号
			double para_base[MAX_DIMENSION] = {0.0};
			double para[MAX_DIMENSION] = {0.0};
			double tilde_para[MAX_DIMENSION] = {0.0};
			double dR_base_simple[MAX_DIMENSION] = {0.0};

			setEdge3D(0, i, &target_axis_base, para_base, td);
			setEdge3D(1, j, target_axis, para, td_base);

			// Newton-Raphson loop
			bool singular_flag = true;
			for (int k = 0; k < max_itr; k++)
			{
				double coord_base[MAX_DIMENSION] = {0.0};
				double coord[MAX_DIMENSION] = {0.0};

				vector<double> diff(MAX_DIMENSION);
				vector<double> sol(MAX_DIMENSION);
				vector<double> J(MAX_DIMENSION * MAX_DIMENSION, 0.0);

				if (check_nonfinite(para_base, info->DIMENSION) || check_nonfinite(para, info->DIMENSION))
					break;
				int e = geo_ele_check(patch, para, info);
				geo_tilde_coord(tilde_para, para, patch, e, info);
				geo_shape_and_dshape(R.data(), dR.data(), tilde_para, e, false, info);
				geo_parameter_coord_R(e, R.data(), coord, info);

				std::fill(dR_base_simple, dR_base_simple + info->DIMENSION, 0.0);
				dR_base_simple[target_axis_base] = 1.0;
				std::copy(para_base, para_base + info->DIMENSION, coord_base);

				double r = 0.0;
				for (int l = 0; l < info->DIMENSION; l++)
				{
					diff[l] = coord_base[l] - coord[l];
					r += diff[l] * diff[l];
				}
				r = sqrt(r);

				// check position
				if (!isfinite(r))
					break;
				if (r < MERGE_ERROR)
					return true;

                // make J
				for (int l = 0; l < info->DIMENSION; l++)
					J[l * info->DIMENSION + 0] += dR_base_simple[l];
				for (int l = 0; l < 2; l++)
					for (int m = 0; m < info->DIMENSION; m++)
						for (int n = 0; n < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[e]]; n++)
							J[m * info->DIMENSION + l + 1] -= dR[n * info->DIMENSION + target_axis[l]] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + n] * (info->DIMENSION + 1) + m];

				// solve [J]{sol} = {diff}
				singular_flag = GaussianElimination2(sol.data(), diff.data(), J.data(), info->DIMENSION);

				if (!singular_flag)
					continue;

				// update
				para_base[target_axis_base] -= sol[0];
				para[target_axis[0]] -= sol[1];
				para[target_axis[1]] -= sol[2];

				if (para_base[target_axis_base] > td.para_end[target_axis_base])
					para_base[target_axis_base] = td.para_end[target_axis_base];
				if (para_base[target_axis_base] < td.para_start[target_axis_base])
					para_base[target_axis_base] = td.para_start[target_axis_base];

				for (int l = 0; l < 2; l++)
				{
					if (para[target_axis[l]] > td_base.para_end[target_axis[l]])
						para[target_axis[l]] = td_base.para_end[target_axis[l]];
					if (para[target_axis[l]] < td_base.para_start[target_axis[l]])
						para[target_axis[l]] = td_base.para_start[target_axis[l]];
				}
			}
		}
	}

	return false;
}


void setEdge3D(int edge_or_face, int edge_or_face_num, int *target_axis, double *para, target_domain &temp)
{
	// edge
	if (edge_or_face == 0)
	{
		*target_axis = edge_or_face_num / 4;
		if (edge_or_face_num == 0)
		{
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_start[1];
			para[2] = temp.para_start[2];
		}
		else if (edge_or_face_num == 1)
		{
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_start[1];
			para[2] = temp.para_end[2];
		}
		else if (edge_or_face_num == 2)
		{
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_end[1];
			para[2] = temp.para_start[2];
		}
		else if (edge_or_face_num == 3)
		{
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_end[1];
			para[2] = temp.para_end[2];
		}
		else if (edge_or_face_num == 4)
		{
			para[0] = temp.para_start[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_start[2];
		}
		else if (edge_or_face_num == 5)
		{
			para[0] = temp.para_start[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_end[2];
		}
		else if (edge_or_face_num == 6)
		{
			para[0] = temp.para_end[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_start[2];
		}
		else if (edge_or_face_num == 7)
		{
			para[0] = temp.para_end[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_end[2];
		}
		else if (edge_or_face_num == 8)
		{
			para[0] = temp.para_start[0];
			para[1] = temp.para_start[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 9)
		{
			para[0] = temp.para_start[0];
			para[1] = temp.para_end[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 10)
		{
			para[0] = temp.para_end[0];
			para[1] = temp.para_start[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 11)
		{
			para[0] = temp.para_end[0];
			para[1] = temp.para_end[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		return;
	}

	// face
	else if (edge_or_face == 1)
	{
		if (edge_or_face_num == 0)
		{
			target_axis[0] = 1, target_axis[1] = 2;
			para[0] = temp.para_start[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 1)
		{
			target_axis[0] = 1, target_axis[1] = 2;
			para[0] = temp.para_end[0];
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 2)
		{
			target_axis[0] = 0, target_axis[1] = 2;
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_start[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 3)
		{
			target_axis[0] = 0, target_axis[1] = 2;
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = temp.para_end[1];
			para[2] = (temp.para_start[2] + temp.para_end[2]) * 0.5;
		}
		else if (edge_or_face_num == 4)
		{
			target_axis[0] = 0, target_axis[1] = 1;
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_start[2];
		}
		else if (edge_or_face_num == 5)
		{
			target_axis[0] = 0, target_axis[1] = 1;
			para[0] = (temp.para_start[0] + temp.para_end[0]) * 0.5;
			para[1] = (temp.para_start[1] + temp.para_end[1]) * 0.5;
			para[2] = temp.para_end[2];
		}
		return;
	}
}


// element connectivitiy
void elementConnectivity(information *info, vector<vector<int>> &ecn)
{
	int total_array = 0, check_n = 0, FE_n = 0, opp_n = 0, point_on_element = pow_int(3, info->DIMENSION);
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
	}

	vector<int> Connectivity(info->Total_Element_to_mesh[Total_mesh] * point_on_element);
	vector<int> Connectivity_ele(info->Total_Element_to_mesh[Total_mesh] * point_on_element);
	vector<int> Connectivity_point(info->Total_Element_to_mesh[Total_mesh] * point_on_element);
	vector<int> Patch_check(info->Total_Patch_to_mesh[Total_mesh] * check_n);
	vector<int> Patch_array(info->Total_Patch_to_mesh[Total_mesh] * total_array);
	vector<int> Face_Edge_info(info->Total_Patch_to_mesh[Total_mesh] * FE_n, -1);
	vector<int> Opponent_patch_num(info->Total_Patch_to_mesh[Total_mesh] * opp_n);

	vector<vector<int>> Connectivity_all_ele(info->Total_Element_to_mesh[Total_mesh] * point_on_element);

	// Make Patch_check
	int Patch_check_counter = 0;
	if (info->DIMENSION == 2)
    {
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];

			// 辺0 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺0 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺2 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺2 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺4 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺4 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺6 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺6 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点0
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点1
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];
			int temp = info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] * (info->No_Control_point[i * info->DIMENSION + 2] - 1);
			int temp_Patch_check_counter;
			int temp_Patch_check[8];

			// 面0 点0
			temp_Patch_check[0] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 面0 点1
			temp_Patch_check[1] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面0 点2
			temp_Patch_check[2] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面0 点3
			temp_Patch_check[3] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面1 点0
			temp_Patch_check[4] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp];
			Patch_check_counter++;
			// 面1 点1
			temp_Patch_check[5] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面1 点2
			temp_Patch_check[6] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面1 点3
			temp_Patch_check[7] = Patch_check_counter;
			Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面2 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点1
			temp_Patch_check_counter = temp_Patch_check[3];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点2
			temp_Patch_check_counter = temp_Patch_check[7];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点1
			temp_Patch_check_counter = temp_Patch_check[1];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点2
			temp_Patch_check_counter = temp_Patch_check[5];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点0
			temp_Patch_check_counter = temp_Patch_check[1];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点3
			temp_Patch_check_counter = temp_Patch_check[5];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点0
			temp_Patch_check_counter = temp_Patch_check[3];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点3
			temp_Patch_check_counter = temp_Patch_check[7];
			Patch_check[Patch_check_counter] = Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
		}
	}

	// Check
	int patch_start = 0;
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
			{
				patch_start = i;
			}

			for (int j = patch_start; j < i; j++)
			{
				int own = i * 16;
				int opp = j * 16;
				for (int k = 0; k < 4; k++)
				{
					int kk = 2 * k;
					for (int l = 0; l < 8; l++)
					{
						// 辺が一致している場合 Face_Edge_info を True
						if (Patch_check[own + kk * 2] == Patch_check[opp + l * 2] && Patch_check[own + kk * 2 + 1] == Patch_check[opp + l * 2 + 1])
						{
							Face_Edge_info[i * 32 + k * 8 + l] = 1;
							Opponent_patch_num[i * 4 + k] = j;
							goto end_loop_a;
						}
					}
				}
				end_loop_a:;
			}
		}
	}
	else if (info->DIMENSION == 3)
	{
		int check_a[4], check_b[4];
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
			{
				patch_start = i;
			}

			for (int j = patch_start; j < i; j++)
			{
				int own = i * 24;
				int opp = j * 24;
				for (int k = 0; k < 6; k++)
				{
					int kk = k * 4;
					for (int m = 0; m < 4; m++)
					{
						check_a[m] = Patch_check[own + kk + m];
					}
					sort(4, check_a);

					for (int l = 0; l < 6; l++)
					{
						int ll = l * 4;
						for (int m = 0; m < 4; m++)
						{
							check_b[m] = Patch_check[opp + ll + m];
						}
						sort(4, check_b);

						// 面が一致している場合 Face_Edge_info を Mode 番号 (m) に
						if (check_a[0] == check_b[0] && check_a[1] == check_b[1] && check_a[2] == check_b[2] && check_a[3] == check_b[3])
						{
							for (int m = 0; m < 4; m++)
							{
								if (Patch_check[own + kk] == Patch_check[opp + ll + m])
								{
									Face_Edge_info[i * 36 + k * 6 + l] = m;
									Opponent_patch_num[i * 6 + k] = j;
									goto end_loop_b;
								}
							}
						}
					}
				}
				end_loop_b:;
			}
		}
	}

	// Make connectivity
	int counter = 0;
	int point_counter = 0;
	int ele_counter = 0;
	if (info->DIMENSION == 2)
	{
		int p_x[9] = {0, 2, 2, 0, 1, 2, 1, 0, 1};
		int p_y[9] = {0, 0, 2, 2, 0, 1, 2, 1, 1};
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (int j = 0; j < i; j++)
			{
				int a, b;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				own += 2 * (a + b);
			}

			int Edge[4] = {0};
    		int Edge_counter = i * 32;
			int temp_n = 0;

			int ii = info->line_No_Total_element[i * info->DIMENSION + 0] * 2 + 1;
			int jj = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;

			// 重なっている辺の Patch_array を作成
			for (int j = 0; j < 4; j++)
			{
				// j == 0 ではなにもしない
				if (j == 1)
				{
					own += ii;
				}
				else if (j == 2)
				{
					own += jj;
				}
				else if (j == 3)
				{
					own += ii;
				}

				for (int k = 0; k < 8; k++)
				{
					if (Face_Edge_info[Edge_counter + j * 8 + k] == 1)
					{
						int opp = 0;
						for (int l = 0; l < Opponent_patch_num[i * 4 + j]; l++)
						{
							int a, b;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							opp += 2 * (a + b);
						}

						Edge[j] = 1;

						int p = k / 2;
						int q = k % 2;

						int temp_ii = info->line_No_Total_element[Opponent_patch_num[i * 4 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[Opponent_patch_num[i * 4 + j] * info->DIMENSION + 1] * 2 + 1;

						if (p == 0)
						{
							temp_n = temp_ii;
						}
						else if (p == 1)
						{
							temp_n = temp_jj;
							opp += temp_ii;
						}
						else if (p == 2)
						{
							temp_n = temp_ii;
							opp += temp_ii + temp_jj;
						}
						else if (p == 3)
						{
							temp_n = temp_jj;
							opp += 2 * temp_ii + temp_jj;
						}

						if (q == 0)
						{
							for (int l = 0; l < temp_n; l++)
							{
								Patch_array[own + l] = Patch_array[opp + l];
							}
							break;
						}
						else if (q == 1)
						{
							for (int l = 0; l < temp_n; l++)
							{
								Patch_array[own + l] = Patch_array[opp + (temp_n - 1) - l];
							}
							break;
						}
					}
				}
			}

			// コネクティビティを作成
			own = 0;
			for (int j = 0; j < i; j++)
			{
				int a, b;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				own += 2 * (a + b);
			}

			int e_x_max = info->line_No_Total_element[i * info->DIMENSION + 0];
			int e_y_max = info->line_No_Total_element[i * info->DIMENSION + 1];
			for (int y = 0; y < e_y_max; y++)
			{
				for (int x = 0; x < e_x_max; x++)
				{
					for (int point = 0; point < 9; point++)
					{
						// point 0
						if (point == 0)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + xi];
							}
							else if (x > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 1];
							}
							else if (y > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 3];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 1
						else if (point == 1)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + xi];
							}
							else if (y > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 2];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 2
						else if (point == 2)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + jj + xi];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 3
						else if (point == 3)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + jj + xi];
							}
							else if (x > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 2];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 4
						else if (point == 4)
						{
							if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + xi];
							}
							else if (y > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 6];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 5
						else if (point == 5)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + eta];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 6
						else if (point == 6)
						{
							if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + ii + jj + xi];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 7
						else if (point == 7)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Patch_array[own + 2 * ii + jj + eta];
							}
							else if (x > 0)
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 5];
							}
							else
							{
								Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 8
						else if (point == 8)
						{
							Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
							Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
							Connectivity_point[counter] = point;
							counter++;
						}

						Connectivity_all_ele[Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point]].emplace_back(ele_counter + y * e_x_max + x);
					}
				}
			}

			// Patch_array の作ってない分を作成
			int point = 0;
			for (int eta = 0; eta < jj; eta++)
			{
				for (int xi = 0; xi < ii; xi++)
				{
					int e_x, e_y;

					if (eta == 0 && Edge[0] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						Patch_array[own + xi] = Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (eta == jj - 1 && Edge[2] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						Patch_array[own + ii + jj + xi] = Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == 0 && Edge[3] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						Patch_array[own + 2 * ii + jj + eta] = Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == ii - 1 && Edge[1] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						Patch_array[own + ii + eta] = Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
				}
			}
			point_counter += e_x_max * e_y_max * 9;
			ele_counter += e_x_max * e_y_max;

			if (i == info->Total_Patch_to_mesh[1] - 1)
			{
				Total_connectivity_glo = counter;
			}
		}
		Total_connectivity = counter;
	}
	else if (info->DIMENSION == 3)
	{
		int p_x[27] = {0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 1, 1, 1, 1, 1};
		int p_y[27] = {0, 0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0, 0, 2, 2, 1, 1, 0, 2, 1, 1, 1};
		int p_z[27] = {0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 1};
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (int j = 0; j < i; j++)
			{
				int a, b, c;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				c = info->line_No_Total_element[j * info->DIMENSION + 2] * 2 + 1;
				own += 2 * (a * b + b * c + a * c);
			}

			int Face[6] = {0};
			int Face_counter = i * 36;
			int own_a = 0, own_b = 0;

			int ii = info->line_No_Total_element[i * info->DIMENSION + 0] * 2 + 1;
			int jj = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;
			int kk = info->line_No_Total_element[i * info->DIMENSION + 2] * 2 + 1;

			// 重なっている面の Patch_array 配列を作成
			for (int j = 0; j < 6; j++)
			{
				if (j == 0)
				{
					own_a = ii;
					own_b = jj;
				}
				else if (j == 1)
				{
					own += ii * jj;
					own_a = ii;
					own_b = jj;
				}
				else if (j == 2)
				{
					own += ii * jj;
					own_a = jj;
					own_b = kk;
				}
				else if (j == 3)
				{
					own += jj * kk;
					own_a = ii;
					own_b = kk;
				}
				else if (j == 4)
				{
					own += ii * kk;
					own_a = jj;
					own_b = kk;
				}
				else if (j == 5)
				{
					own += jj * kk;
					own_a = ii;
					own_b = kk;
				}

				for (int k = 0; k < 6; k++)
				{
					if (Face_Edge_info[Face_counter + j * 6 + k] >= 0)
					{
						int opp = 0;
						for (int l = 0; l < Opponent_patch_num[i * 6 + j]; l++)
						{
							int a, b, c;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							c = info->line_No_Total_element[l * info->DIMENSION + 2] * 2 + 1;
							opp += 2 * (a * b + b * c + a * c);
						}

						Face[j] = 1;

						int temp_ii = info->line_No_Total_element[Opponent_patch_num[i * 6 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[Opponent_patch_num[i * 6 + j] * info->DIMENSION + 1] * 2 + 1;
						int temp_kk = info->line_No_Total_element[Opponent_patch_num[i * 6 + j] * info->DIMENSION + 2] * 2 + 1;

						for (int l = 1; l <= k; l++)
						{
							if (l == 1)
							{
								opp += temp_ii * temp_jj;
							}
							else if (l == 2)
							{
								opp += temp_ii * temp_jj;
							}
							else if (l == 3)
							{
								opp += temp_jj * temp_kk;
							}
							else if (l == 4)
							{
								opp += temp_ii * temp_kk;
							}
							else if (l == 5)
							{
								opp += temp_jj * temp_kk;
							}
						}
						
						if (((j == 0 || j == 2 || j == 5) && (k == 1 || k == 3 || k == 4)) || ((k == 0 || k == 2 || k == 5) && (j == 1 || j == 3 || j == 4)))
						{
							if (Face_Edge_info[Face_counter + j * 6 + k] == 0)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + l * own_a + m];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k]  == 1)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + m * own_b + ((own_b - 1) - l)];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k]  == 2)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + ((own_b - 1) - l) * own_a + ((own_a - 1) - m)];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k]  == 3)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + ((own_a - 1) - m) * own_b + l];
								break;
							}
						}
						else
						{
							if (Face_Edge_info[Face_counter + j * 6 + k] == 0)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + m * own_b + l];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k] == 1)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + ((own_a - 1) - m)] = Patch_array[opp + l * own_a + m];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k] == 2)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + l * own_a + m] = Patch_array[opp + ((own_a - 1) - m) * own_b + ((own_b - 1) - l)];
								break;
							}
							else if (Face_Edge_info[Face_counter + j * 6 + k] == 3)
							{
								for (int l = 0; l < own_b; l++)
									for (int m = 0; m < own_a; m++)
										Patch_array[own + ((own_b - 1) - l) * own_a + m] = Patch_array[opp + l * own_a + m];
								break;
							}
						}
					}
				}
			}

			// コネクティビティを作成
			own = 0;
			for (int j = 0; j < i; j++)
			{
				int a, b, c;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				c = info->line_No_Total_element[j * info->DIMENSION + 2] * 2 + 1;
				own += 2 * (a * b + b * c + a * c);
			}

			int temp1, temp2, temp3, temp4, temp5;
			temp1 = ii * jj;
			temp2 = temp1 + ii * jj;
			temp3 = temp2 + jj * kk;
			temp4 = temp3 + ii * kk;
			temp5 = temp4 + jj * kk;

			int e_x_max = info->line_No_Total_element[i * info->DIMENSION + 0];
			int e_y_max = info->line_No_Total_element[i * info->DIMENSION + 1];
			int e_z_max = info->line_No_Total_element[i * info->DIMENSION + 2];
			for (int z = 0; z < e_z_max; z++)
			{
				for (int y = 0; y < e_y_max; y++)
				{
					for (int x = 0; x < e_x_max; x++)
					{
						for (int point = 0; point < 27; point++)
						{
							int temp = point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + point;

							// point 0
							if (point == 0)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 1];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 3];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 4];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 1
							else if (point == 1)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 2];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 5];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 2
							else if (point == 2)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 3
							else if (point == 3)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 2];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 4
							else if (point == 4)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 5];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 5
							else if (point == 5)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 6
							else if (point == 6)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 7
							else if (point == 7)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 6];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 8
							else if (point == 8)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 10];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 12];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 9
							else if (point == 9)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 13];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 10
							else if (point == 10)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 11
							else if (point == 11)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 9];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 15];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 12
							else if (point == 12)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 13
							else if (point == 13)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 14
							else if (point == 14)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 15
							else if (point == 15)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 13];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 16
							else if (point == 16)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 17];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 19];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 17
							else if (point == 17)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 18];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 18
							else if (point == 18)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 19
							else if (point == 19)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 18];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 20
							else if (point == 20)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (x > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 21];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 21
							else if (point == 21)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp4 + zeta * jj + eta];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 22
							else if (point == 22)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 23];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 23
							else if (point == 23)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									Connectivity[temp] = Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 24
							else if (point == 24)
							{
								if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									Connectivity[temp] = Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 25];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 25
							else if (point == 25)
							{
								if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									Connectivity[temp] = Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									Connectivity[temp] = counter;
									Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 26
							else if (point == 26)
							{
								Connectivity[temp] = counter;
								Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
								Connectivity_point[counter] = point;
								counter++;
							}

							Connectivity_all_ele[Connectivity[temp]].emplace_back(ele_counter + z * e_y_max * e_x_max + y * e_x_max + x);
						}
					}
				}
			}

			// Patch_array の作ってない分を作成
			int point = 0;
			for (int zeta = 0; zeta < kk; zeta++)
			{
				for (int eta = 0; eta < jj; eta++)
				{
					for (int xi = 0; xi < ii; xi++)
					{
						int e_x, e_y, e_z;

						if (zeta == 0 && Face[0] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + eta * ii + xi] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (zeta == kk - 1 && Face[1] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + temp1 + eta * ii + xi] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == 0 && Face[3] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + temp3 + zeta * ii + xi] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == jj - 1 && Face[5] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + temp5 + zeta * ii + xi] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == 0 && Face[2] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + temp2 + zeta * jj + eta] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == ii - 1 && Face[4] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							Patch_array[own + temp4 + zeta * jj + eta] = Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
					}
				}
			}
			point_counter += e_x_max * e_y_max * e_z_max * 27;
			ele_counter += e_x_max * e_y_max * e_z_max;

			if (i == info->Total_Patch_to_mesh[1] - 1)
			{
				Total_connectivity_glo = counter;
			}
		}
		Total_connectivity = counter;
	}

	// Make ecn
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < point_on_element; j++)
		{
			for (size_t k = 0; k < Connectivity_all_ele[Connectivity[i * point_on_element + j]].size(); k++)
			{
				ecn[i].emplace_back(Connectivity_all_ele[Connectivity[i * point_on_element + j]][k]);
			}
		}
		// sort
		sortUniqueErase(ecn[i]);
	}
}