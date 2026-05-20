// header
#include "_header.hpp"
#include "_sub.hpp"
#include "_MI3D.hpp"

using namespace std;

void Fracture_analysis(information *info)
{
	// allocation J point info
	Allocation(10, info);

	// set Displacement to disp
	info->disp = info->Displacement;

	// set 0 step displacement for J calculation
	Init_increment_field(info);

	// calculate J integral
	calc_J(info);
}


void Init_increment_field(information *info)
{
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		info->disp_increment[i] = 0.0;
		info->disp_overlay_increment[i] = 0.0;
	}

	gp_switch(true, info);

	#pragma omp parallel for collapse(1)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < info->gp[i].n(); j++)
			info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
	}
}

vector<debug_J_vtk> debug_J;
vector<int> e_index;

bool get_info_J(vector<J_point_info> &jp_list, vector<J_info> &J_list, information *info, vector<double> &circle_r)
{
	if (info->c.INTEGRAL_DOMAIN_TYPE == 0) // 特異パッチ
	{
		int temp_p = info->Total_Patch_to_mesh[1];
	
		jp_list.emplace_back(info, temp_p + 0, 1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 4, 1, 0, 1, 2);
	
		J_list.emplace_back(temp_p + 0, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 1, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 2, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 3, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 4, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 5, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 6, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 7, 1, 0, 1, 10 /*sub_n*/);
	}
	else if (info->c.INTEGRAL_DOMAIN_TYPE == 1) // SQUARE_PATCH
	{
		jp_list.emplace_back(info, 0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);

		J_list.emplace_back(0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
		J_list.emplace_back(1 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
		J_list.emplace_back(2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
		J_list.emplace_back(3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	}
	else if (info->c.INTEGRAL_DOMAIN_TYPE == 2) // 卒研でつかったローカル
	{
		int temp_p = info->Total_Patch_to_mesh[1];
	
		jp_list.emplace_back(info, temp_p + 0, 1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 1, 1, 0, 1, 2);
	
		J_list.emplace_back(temp_p + 0, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 2, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 4, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 6, 1, 0, 0, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 1, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 3, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 5, 1, 0, 1, 10 /*sub_n*/);
		J_list.emplace_back(temp_p + 7, 1, 0, 1, 10 /*sub_n*/);
	}
	else if (info->c.INTEGRAL_DOMAIN_TYPE == 3) // Surface Center Crack at the hole in a Plate
	{
		int temp_p = info->Total_Patch_to_mesh[1];
	
		jp_list.emplace_back(info, temp_p + 0,  1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 8,  1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 16, 1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 24, 1, 0, 1, 2);
	
		J_list.emplace_back(temp_p + 0, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 1, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 2, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 3, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 4, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 5, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 6, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 7, 1, 0, 0,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 8, 1, 0, 1,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 9, 1, 0, 1,  3/*sub_n*/);
		J_list.emplace_back(temp_p + 10, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 11, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 12, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 13, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 14, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 15, 1, 0, 1, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 16, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 17, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 18, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 19, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 20, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 21, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 22, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 23, 1, 0, 2, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 24, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 25, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 26, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 27, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 28, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 29, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 30, 1, 0, 3, 3 /*sub_n*/);
		J_list.emplace_back(temp_p + 31, 1, 0, 3, 3 /*sub_n*/);
	}
	else if (info->c.INTEGRAL_DOMAIN_TYPE == 4) //矩形パッチ
	{
		// #define SQUARE_PATCH
		int temp_p = info->Total_Patch_to_mesh[1];
	
		jp_list.emplace_back(info, temp_p + 0,  1, 0, 1, 2);
		jp_list.emplace_back(info, temp_p + 4,  1, 0, 1, 2);
	
		J_list.emplace_back(temp_p + 0, 1, 0, 0, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 1, 1, 0, 0, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 2, 1, 0, 0, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 3, 1, 0, 0, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 4, 1, 0, 1, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 5, 1, 0, 1, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 6, 1, 0, 1, 3/*sub_n*/);
		J_list.emplace_back(temp_p + 7, 1, 0, 1, 3/*sub_n*/);
	}
	else
	{
		printf("Error: Invalid INTEGRAL_DOMAIN_TYPE\n");
		exit(1);
	}
	return true;
}


// calc J integral
void calc_J(information *info)
{
	cout << "start calc_J" << endl;

	#if 0
	return;
	#endif

	constexpr int sub_gp_n_1d = 4;
	constexpr int sub_gp_n = sub_gp_n_1d * sub_gp_n_1d * sub_gp_n_1d;

	constexpr int face_n = 2;

	static vector<J_point_info> jp_list;
	static vector<J_info> J_list;

	static vector<sub_ele_J> sub_ele;

	static vector<double> sub_gp_para(sub_gp_n * info->DIMENSION, 0.0);
	static vector<double> sub_gp_weight(sub_gp_n, 0.0);
	static vector<double> sub_gp_para_1d_dummy(sub_gp_n * info->DIMENSION, 0.0);
	static vector<double> sub_gp_weight_1d_dummy(sub_gp_n, 0.0);

	static vector<double> circle_r(2, 1.0);

	static bool isInitialized = false;
	if (!isInitialized)
	{
		Gauss_point(info, sub_gp_n_1d, sub_gp_weight.data(), sub_gp_para.data(), sub_gp_weight_1d_dummy.data(), sub_gp_para_1d_dummy.data());

		isInitialized = get_info_J(jp_list, J_list, info, circle_r);

		// set e_list
		for (size_t i = 0; i < J_list.size(); i++)
		{
			int current_p = J_list[i].p;
			int current_e = 0;
			for (int j = 0; j < current_p; j++)
			{
				int temp_e = 1;
				for (int k = 0; k < info->DIMENSION; k++)
					temp_e *= info->No_Control_point[j * info->DIMENSION + k] - info->Order[j * info->DIMENSION + k];
				current_e += temp_e;
			}
			int ele_n = 1;
			for (int j = 0; j < info->DIMENSION; j++)
				ele_n *= info->No_Control_point[current_p * info->DIMENSION + j] - info->Order[current_p * info->DIMENSION + j];

			// resize e_list
			int dir = J_list[i].crack_dir;
			int size = info->No_Control_point[current_p * info->DIMENSION + dir] - info->Order[current_p * info->DIMENSION + dir];
			J_list[i].e_list.resize(size);
	
			for (int j = 0; j < ele_n; j++)
			{
				int e = j + current_e;
				int enc = info->ENC[e * info->DIMENSION + J_list[i].crack_dir];
				J_list[i].e_list[enc].emplace_back(e);
			}
		}

		// init jp_list
		for (size_t i = 0; i < jp_list.size(); i++)
		{
			size_t index = 0;
			for (size_t j = 0; j < J_list.size(); j++)
			{
				if (J_list[j].p == jp_list[i].p)
				{
					index = j;
					break;
				}
			}
			J_info &temp = J_list[index];
			jp_list[i].init(temp, info);
		}

		// debug output
		#if 0
		for (size_t i = 0; i < J_list.size(); i++)
		{
			printf("J_list[%d]: %d %d %d\n", i, J_list[i].p, J_list[i].crack_dir, J_list[i].r_dir);
			for (int j = 0; j < J_list[i].e_list.size(); j++)
			{
				printf("\t");
				for (int k = 0; k < J_list[i].e_list[j].size(); k++)
					printf("%d ", J_list[i].e_list[j][k]);
				printf("\n");
			}
		}
		#endif

		// debug vtk output
		#if 0
		// init debug_J_vtk
		e_index.resize(info->Total_Element_to_mesh[Total_mesh], -1);
		size_t total_num_J_vtk = 0;
		int pow_ng = pow_int(info->c.NUM_GAUSS_POINTS, info->DIMENSION);
		const int face_n = 2;
		int e_counter = 0;
		for (size_t i = 0; i < J_list.size(); i++)
		{
			for (size_t j = 0; j < J_list[i].e_list.size(); j++)
			{
				total_num_J_vtk += J_list[i].e_list[j].size() * J_list[i].sub_n * face_n * pow_ng;
				
				for (size_t k = 0; k < J_list[i].e_list[j].size(); k++)
				{
					int e = J_list[i].e_list[j][k];
					if (e_index[e] == -1)
					{
						e_index[e] = e_counter;
						e_counter++;
					}
				}
			}
		}
		debug_J.reserve(total_num_J_vtk);
		for (size_t i = 0; i < total_num_J_vtk; i++)
		{
			debug_J.emplace_back(info->DIMENSION);
		}
		#endif
		
		// calc vertual crack extension area
		for (size_t i = 0; i < jp_list.size(); i++)
		{
			J_point_info &jp = jp_list[i];
			int index = 0;
			for (size_t j = 0; j < J_list.size(); j++)
			{
				if (J_list[j].p == jp_list[i].p)
				{
					index = j;
					break;
				}
			}

			// ローカルかグローバルかの判定
			bool islocal = false;
			if (i < (size_t)info->Total_Element_on_mesh[1])
				islocal = true;

			for (size_t j = 0; j < J_list[index].e_list.size(); j++)
			{
				// search curve element
				int e_on_curve = -1;
				int crack_dir = jp.crack_dir;
				int target_e = J_list[index].e_list[j][0]; // e_list has same ENC
				for (size_t k = 0; k < jp.ele_in_normal_patch.size(); k++)
				{
					int e = jp.ele_in_normal_patch[k];
					if (info->ENC[e * info->DIMENSION + crack_dir] == info->ENC[target_e * info->DIMENSION + crack_dir])
					{
						e_on_curve = e;
						break;
					}
				}
				if (e_on_curve == -1)
				{
					printf("error: e_on_curve not found\n");
					exit(1);
				}

				// calc vertual crack extension area
				make_virtual_crack_extension_area(jp, J_list[i], j, e_on_curve, islocal, info);
			}
		}

		// alloc sub_ele
		cout << "start allocation" << endl;
		#pragma omp parallel
		{
			// スレッドローカルの一時バッファ
			vector<sub_ele_J> local_sub_ele;
			local_sub_ele.reserve(1000); // 必要に応じて適当な見積もり

			#pragma omp for nowait schedule(dynamic)
			for (size_t i = 0; i < J_list.size(); i++)
				for (size_t j = 0; j < J_list[i].e_list.size(); j++)
					for (size_t k = 0; k < J_list[i].e_list[j].size(); k++)
						for (int l = 0; l < J_list[i].sub_n; l++)
							for (int m = 0; m < face_n; m++)
							{
								int jp_list_index = J_list[i].normal_info_num;
								int e = J_list[i].e_list[j][k];
								local_sub_ele.emplace_back(sub_gp_n_1d, sub_gp_n, e, info->Element_patch[e], J_list[i].sub_n, face_n, j, l, m, i, jp_list_index, J_list[i].crack_dir, sub_gp_para.data(), sub_gp_weight.data());
							}

			#pragma omp critical (sub_ele_merge_lock)
			{
				// 各スレッドの結果を結合
				sub_ele.insert(sub_ele.end(), local_sub_ele.begin(), local_sub_ele.end());
			}
		}
		
		// init sub_ele
		cout << "start init" << endl;
		#pragma omp parallel for schedule(dynamic)
		for (size_t i = 0; i < sub_ele.size(); i++)
			sub_ele[i].init(info);
	}

	// calc J
	J_sub_element(J_list, jp_list, sub_ele, info);

	// debug output SIF in MODE1
	std::ofstream output_sif_file("output_sif.txt");
	if (!output_sif_file.is_open())
	{
		std::cerr << "Error: Cannot open output_sif.txt" << std::endl;
	}
	else
	{
		output_sif_file << std::scientific << std::setprecision(15);
		output_sif_file << "patches, points, x, y, theta_deg, K, J\n";
	}
	for (size_t i = 0; i < jp_list.size(); i++)
	{
		for (size_t j = 0; j < jp_list[i].J_val.size(); j++)
		{	
			double J = jp_list[i].J_val[j];
			if (info->c.CRACK_TYPE == 1) // symmetric
				J *= 2.0;
			double E_prime = E / (1.0 - nu * nu);
			double K = sqrt(E_prime * J);
			double theta_deg = atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI;
			printf("i: %zu j: %zu (x,y,θ,K,J): %.15e %.15e %.15e %.15e %.15e\n", i, j, jp_list[i].coord[j][0], jp_list[i].coord[j][1], theta_deg, K, J);
			if (output_sif_file.is_open())
				output_sif_file << i << " " << j << " " << jp_list[i].coord[j][0] << " " << jp_list[i].coord[j][1] << " " << theta_deg << " " << K << " " << J << "\n";

		}
	}
	printf("\n");
	if (output_sif_file.is_open())
		output_sif_file.close();

	// debug vtk output by using debug_J
	#if 0
	writeVTKFile(debug_J, "debug_J.vtk", info->DIMENSION);
	#endif
}


void writeVTKFile(const std::vector<debug_J_vtk>& debug_J, const std::string& filename, int dim)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // VTKファイルヘッダー
    file << "# vtk DataFile Version 3.0\n";
    file << "Debug J VTK output\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n\n";
    
    // ポイント数
    size_t numPoints = debug_J.size();
    file << "POINTS " << numPoints << " double\n";
    
    // 座標データの出力
    for (const auto& point : debug_J) {
        for (int i = 0; i < dim; i++) {
            file << std::fixed << std::setprecision(6) << point.coord[i];
            if (i < dim - 1) file << " ";
        }
        // 2Dの場合は3次元目に0を追加
        if (dim == 2) {
            file << " 0.0";
        }
        file << "\n";
    }

    // セルデータ（各点を独立したセルとして扱う）
    file << "\nCELLS " << numPoints << " " << numPoints * 2 << "\n";
    for (size_t i = 0; i < numPoints; i++) {
        file << "1 " << i << "\n";
    }

    // セルタイプ（VTK_VERTEX = 1）
    file << "\nCELL_TYPES " << numPoints << "\n";
    for (size_t i = 0; i < numPoints; i++) {
        file << "1\n";
    }

    // ポイントデータ
    file << "\nPOINT_DATA " << numPoints << "\n";

    // スカラーデータ q
    file << "SCALARS q double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& point : debug_J) {
        file << std::fixed << std::setprecision(6) << point.q << "\n";
    }

    // ベクトルデータ q_grad
    file << "\nVECTORS q_grad double\n";
    for (const auto& point : debug_J) {
        for (int i = 0; i < dim; i++) {
            file << std::fixed << std::setprecision(6) << point.q_grad[i];
            if (i < dim - 1) file << " ";
        }
        // 2Dの場合は3次元目に0を追加
        if (dim == 2) {
            file << " 0.0";
        }
        file << "\n";
    }

    // ベクトルデータ normal
    file << "\nVECTORS normal double\n";
    for (const auto& point : debug_J) {
        for (int i = 0; i < dim; i++) {
            file << std::fixed << std::setprecision(6) << point.normal[i];
            if (i < dim - 1) file << " ";
        }
        // 2Dの場合は3次元目に0を追加
        if (dim == 2) {
            file << " 0.0";
        }
        file << "\n";
    }

    file.close();
    std::cout << "VTK file written to: " << filename << std::endl;
}


void make_virtual_crack_extension_area(J_point_info &jp, J_info &J, size_t line, int e_on_curve, bool islocal, information *info)
{
	vector<double> &J_integral_area = jp.J_integral_area;
	vector<vector<double>> &normal = jp.normal;
	vector<vector<double>> &coord = jp.coord;

	double delta = 1.0 / static_cast<double>(J.sub_n);

	for (int i = 0; i < J.sub_n; i++)
	{
		// calc para_disp
		vector<double> para_disp_c(info->DIMENSION, 0.0);
		vector<double> para_disp_s(info->DIMENSION, 0.0);
		vector<double> para_disp_e(info->DIMENSION, 0.0);
		vector<double> para_disp_n(info->DIMENSION, 0.0);
		for (int m = 0; m < info->DIMENSION; m++)
		{
			// para
			if (m == J.crack_dir)
			{
				double normalized_para_c = (static_cast<double>(i) + 0.5) * delta;
				double normalized_para_s = (static_cast<double>(i) + 0.0) * delta;
				double normalized_para_e = (static_cast<double>(i) + 1.0) * delta;
				double normalized_para_n = (static_cast<double>(i) + 0.0) * delta;
				para_disp_c[m] = normalized_para_c * 2.0 - 1.0;
				para_disp_s[m] = normalized_para_s * 2.0 - 1.0;
				para_disp_e[m] = normalized_para_e * 2.0 - 1.0;
				para_disp_n[m] = normalized_para_n * 2.0 - 1.0;
			}
			else
			{
				para_disp_c[m] = -1.0;
				para_disp_s[m] = -1.0;
				para_disp_e[m] = -1.0;
				para_disp_n[m] = -1.0;
			}
		}
		para_disp_n[J.r_dir] = 1.0;
		
		vector<double> coord_c(info->DIMENSION, 0.0);
		vector<double> coord_s(info->DIMENSION, 0.0);
		vector<double> coord_e(info->DIMENSION, 0.0);
		vector<double> coord_n(info->DIMENSION, 0.0);
		if (!islocal)
		{
			vector<double> R_c(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_s(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_e(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_n(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R_c.data(), para_disp_c.data(), e_on_curve, info);
			shape_and_dshape(R_s.data(), para_disp_s.data(), e_on_curve, info);
			shape_and_dshape(R_e.data(), para_disp_e.data(), e_on_curve, info);
			shape_and_dshape(R_n.data(), para_disp_n.data(), e_on_curve, info);
			for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e_on_curve]]; m++)
			{
				for (int n = 0; n < info->DIMENSION; n++)
				{
					int id = info->Controlpoint_of_Element[e_on_curve * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + n;
					coord_c[n] += R_c[m] * info->Node_Coordinate[id];
					coord_s[n] += R_s[m] * info->Node_Coordinate[id];
					coord_e[n] += R_e[m] * info->Node_Coordinate[id];
					coord_n[n] += R_n[m] * info->Node_Coordinate[id];
				}
			}
		}
		else
		{
			// calc global para
			vector<double> para_glo_c(info->DIMENSION, 0.0);
			vector<double> para_glo_s(info->DIMENSION, 0.0);
			vector<double> para_glo_e(info->DIMENSION, 0.0);
			vector<double> para_glo_n(info->DIMENSION, 0.0);
			int ele_glo_c = trans_local_para_to_global_para(e_on_curve, para_disp_c.data(), para_glo_c.data(), info);
			int ele_glo_s = trans_local_para_to_global_para(e_on_curve, para_disp_s.data(), para_glo_s.data(), info);
			int ele_glo_e = trans_local_para_to_global_para(e_on_curve, para_disp_e.data(), para_glo_e.data(), info);
			int ele_glo_n = trans_local_para_to_global_para(e_on_curve, para_disp_n.data(), para_glo_n.data(), info);

			// calc local coord (global para)
			vector<double> R_c(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_s(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_e(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_n(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R_c.data(), para_glo_c.data(), ele_glo_c, info);
			shape_and_dshape(R_s.data(), para_glo_s.data(), ele_glo_s, info);
			shape_and_dshape(R_e.data(), para_glo_e.data(), ele_glo_e, info);
			shape_and_dshape(R_n.data(), para_glo_n.data(), ele_glo_n, info);
			auto accumulate_coord = [&](int ele, const vector<double>& R, vector<double>& local_coord)
			{
			    int patch = info->Element_patch[ele];
			    int cp_n = info->No_Control_point_ON_ELEMENT[patch];
			
			    for (int m = 0; m < cp_n; ++m)
			    {
			        int cp = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + m];
			        int base = cp * (info->DIMENSION + 1);
			        for (int n = 0; n < info->DIMENSION; ++n)
			            local_coord[n] += R[m] * info->Node_Coordinate[base + n];
			    }
			};
			accumulate_coord(ele_glo_c, R_c, coord_c);
			accumulate_coord(ele_glo_s, R_s, coord_s);
			accumulate_coord(ele_glo_e, R_e, coord_e);
			accumulate_coord(ele_glo_n, R_n, coord_n);
		}

		for (int m = 0; m < info->DIMENSION; m++)
			coord[line * J.sub_n + i][m] = coord_c[m];

		vector<double> vec_a(info->DIMENSION, 0.0);
		vector<double> vec_b(info->DIMENSION, 0.0);
		vector<double> vec_c(info->DIMENSION, 0.0);
		vector<double> vec_e(info->DIMENSION, 0.0);
		for (int m = 0; m < info->DIMENSION; m++)
		{
			vec_a[m] = coord_c[m] - coord_s[m];
			vec_b[m] = coord_e[m] - coord_c[m];
			vec_c[m] = coord_e[m] - coord_s[m];
			vec_e[m] = coord_n[m] - coord_c[m];
		}

		// norm of vec_c
		// 今回のq関数では、仮想亀裂進展面積ΔAは∫_{h_0} q dx_3 = h_0 / 2.0
		double norm = 0.0;
		for (int m = 0; m < info->DIMENSION; m++)
			norm += vec_c[m] * vec_c[m];
		J_integral_area[line * J.sub_n + i] = sqrt(norm) * 0.5;

		// cross (a, b)
		vector<double> cross_vec_ab(info->DIMENSION, 0.0);
		bool isValid = calc_normalized_cross(cross_vec_ab.data(), vec_a.data(), vec_b.data());

		// cross (e, c)
		vector<double> cross_vec_ec(info->DIMENSION, 0.0);
		calc_normalized_cross(cross_vec_ec.data(), vec_e.data(), vec_c.data());

		// if cross (a, b) are parallel, use cross (c, ec)
		if (!isValid)
		{
			// cross (c, ec)
			vector<double> &current_normal = normal[line * J.sub_n + i];
			calc_normalized_cross(current_normal.data(), vec_c.data(), cross_vec_ec.data());
		}
		else
		{
			// cross (c, ab)
			vector<double> &current_normal = normal[line * J.sub_n + i];
			calc_normalized_cross(current_normal.data(), vec_c.data(), cross_vec_ab.data());

			// check direction
			double dot = 0.0;
			for (int m = 0; m < info->DIMENSION; m++)
				dot += current_normal[m] * vec_a[m];
			if (dot < 0.0)
				for (int m = 0; m < info->DIMENSION; m++)
					current_normal[m] *= -1.0;
		}

		// cout << "line: " << line << " i: " << i << " normal: ";
		// for (int m = 0; m < info->DIMENSION; m++)			cout << normal[line * J.sub_n + i][m] << " ";
		// cout << endl;
		// fflush(stdout);
		// exit(0);

		#if defined(INTERACTION_INTEGRAL_METHOD)
		// opening direction (n_open): crack plane normal
		// use cross(vec_e, vec_c) which is normal to the plane spanned by propagation (vec_e) and front tangent (vec_c)
		vector<double> n_open(info->DIMENSION, 0.0);
		for (int m = 0; m < info->DIMENSION; m++) n_open[m] = cross_vec_ec[m];
		vector_normalize(n_open.data(), info->DIMENSION);

		// propagation direction (b_prop): use computed normal[]
		vector<double> b_prop(info->DIMENSION, 0.0);
		for (int m = 0; m < info->DIMENSION; m++) b_prop[m] = normal[line * J.sub_n + i][m];
		vector_normalize(b_prop.data(), info->DIMENSION);

		// tangent (t): build from (prop, open) to keep a right-handed basis
		// columns are (x', y', z') = (prop, open, tangent) so z' = x' x y'
		vector<double> t(info->DIMENSION, 0.0);
		calc_normalized_cross(t.data(), b_prop.data(), n_open.data());

		for (int m = 0; m < info->DIMENSION; m++)
		{
			jp.Q[line * J.sub_n + i][m][0] = b_prop[m]; // column 0: propagation (x')
			jp.Q[line * J.sub_n + i][m][1] = n_open[m]; // column 1: opening (y')
			jp.Q[line * J.sub_n + i][m][2] = t[m];      // column 2: front tangent (z')
		}
		#endif

		#if defined(YONTENMAGE)
		constexpr int sub_gp_n_1d = 4;
		constexpr double sub_gp_para_1d[4] = { -0.861136311594052575224, -0.3399810435848562648027, 0.3399810435848562648027, 0.861136311594052575224 };
		constexpr double sub_gp_weight_1d[4] = { 0.3478548451374538573731, 0.6521451548625461426269, 0.6521451548625461426269, 0.3478548451374538573731 };

		for (int i = 0; i < J.sub_n; i++)
		{
			// delta_arc_length
			double current_delta_arc_length = 0.0;
			double current_delta_arc_length_at_half_para = 0.0;
			for (int j = 0; j < sub_gp_n_1d; j++)
			{
				// make para
				vector<double> current_para(info->DIMENSION, 0.0);
				vector<double> current_para_hp(info->DIMENSION, 0.0);
				for (int k = 0; k < info->DIMENSION; k++)
				{
					if (k == J.crack_dir)
					{
						double normalized_para = (static_cast<double>(i) + (sub_gp_para_1d[j] + 1.0) * 0.5) * delta;
						double normalized_para_jp = (static_cast<double>(i) + (sub_gp_para_1d[j] + 1.0) * 0.25) * delta;
						current_para[k] = normalized_para * 2.0 - 1.0;
						current_para_hp[k] = normalized_para_jp * 2.0 - 1.0;
					}
					else
					{
						current_para[k] = -1.0;
						current_para_hp[k] = -1.0;
					}
				}

				// calc Jacobian
				double J = Make_Jac_anypoint(e_on_curve, current_para.data(), info);
				double J_hp = Make_Jac_anypoint(e_on_curve, current_para_hp.data(), info);
				current_delta_arc_length += sub_gp_weight_1d[j] * J;
				current_delta_arc_length_at_half_para += sub_gp_weight_1d[j] * J_hp * 0.5; // det|∂ξ/∂ξ'|=0.5

			}
			jp.delta_arc_length[line * J.sub_n + i] = current_delta_arc_length;
			jp.delta_arc_length_at_half_para[line * J.sub_n + i] = current_delta_arc_length_at_half_para;
		}
		#endif

		#if 0
		#pragma omp critical
		{
			// debug output
			printf("J_point[%d][%d] Q matrix:\n", line, i);
			for (int m = 0; m < info->DIMENSION; m++)
			{
				printf("\t");
				for (int n = 0; n < info->DIMENSION; n++)
				{
					printf("%.6e ", jp.Q[line * J.sub_n + i][m][n]);
				}
				printf("\n");
			}
			exit(0);
		}
		#endif
	}

	for (int i = 0; i < J.sub_n; i++)
	{
		// debug output
		#if 0
		printf("J_integral_area[%d]: %.6e\n", i, J_integral_area[line * J.sub_n + i]);
		#endif
	}
}


void J_sub_element(vector<J_info> &J_list, vector<J_point_info> &jp_list, vector<sub_ele_J> &sub_ele, information *info)
{
	cout << "start J_sub_element" << endl;

	const int face_n = 2;

	vector<double> D_e(D_MATRIX_SIZE * D_MATRIX_SIZE, 0.0);
    for (int i = 0; i < D_MATRIX_SIZE; i++){
        for (int j = 0; j < D_MATRIX_SIZE; j++)
            D_e[i * D_MATRIX_SIZE + j] = info->D[i * D_MATRIX_SIZE + j];
	}

	// init J_val
	for (size_t i = 0; i < jp_list.size(); i++){
		for (size_t j = 0; j < jp_list[i].J_val.size(); j++)
			jp_list[i].J_val[j] = 0.0;
	}
	
	// init IIM_K_val
	#if defined(INTERACTION_INTEGRAL_METHOD)
	for (size_t i = 0; i < jp_list.size(); i++)
		for (size_t j = 0; j < jp_list[i].IIM_K_val.size(); j++)
			for (size_t k = 0; k < jp_list[i].IIM_K_val[j].size(); k++)
				jp_list[i].IIM_K_val[j][k] = 0.0;
	#endif

	#pragma omp parallel
	{
		// cache J_val
		vector<vector<double>> local_J_val(jp_list.size());
		for (size_t i = 0; i < jp_list.size(); i++)
			local_J_val[i].resize(jp_list[i].J_val.size(), 0.0);

		// cache IIM_K_val
		#if defined(INTERACTION_INTEGRAL_METHOD)
		vector<vector<vector<double>>> local_IIM_K_val(jp_list.size());
		for (size_t i = 0; i < jp_list.size(); i++)
			local_IIM_K_val[i].resize(jp_list[i].IIM_K_val.size(), vector<double>(3, 0.0));
		#endif

		#pragma omp for schedule(dynamic) nowait
		for (size_t i = 0; i < sub_ele.size(); i++)
		{
			sub_ele_J &current_sub_ele = sub_ele[i];
			int e = current_sub_ele.e;

			J_info &current_J_list = J_list[current_sub_ele.J_list_index];
			const vector<double> &J_integral_area = jp_list[current_sub_ele.jp_list_index].J_integral_area;
			const vector<vector<double>> &normal = jp_list[current_sub_ele.jp_list_index].normal;
			int index_J = current_sub_ele.line * current_J_list.sub_n + current_sub_ele.slice;
			double delta = 1.0 / (static_cast<double>(current_J_list.sub_n * face_n));

			for (int j = 0; j < current_sub_ele.gp_num; j++)
			{
				// calc para
				vector<double> disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> strain_trial(D_MATRIX_SIZE, 0.0);
				vector<double> bl(D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
				vector<double> b_disp_grad(info->DIMENSION * info->DIMENSION * MAX_KIEL_SIZE, 0.0);
				vector<double> u(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
				double *para_ptr = current_sub_ele.para_tilde.data() + j * info->DIMENSION;

				Make_B_Linear(e, para_ptr, bl.data(), info);
				{
					vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
					if (e < info->Total_Element_to_mesh[1])
						Make_B_component(e, para_ptr, b.data(), info);
					else
						Make_B_component_for_SSIGA(e, para_ptr, b.data(), info);

					for (int cp = 0; cp < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; cp++)
						for (int comp = 0; comp < info->DIMENSION; comp++)
							for (int dir = 0; dir < info->DIMENSION; dir++)
							{
								int offset = comp * info->DIMENSION + dir;
								b_disp_grad[offset * MAX_KIEL_SIZE + cp * info->DIMENSION + comp] = b[dir * MAX_NO_CP_ON_ELEMENT + cp];
							}
				}

				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
						u[k * info->DIMENSION + l] = info->disp[id];
					}
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							strain_trial[k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
						}
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int offset = l * info->DIMENSION + k;
						for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]] * info->DIMENSION; m++)
							disp_grad[k * info->DIMENSION + l] += u[m] * b_disp_grad[offset * MAX_KIEL_SIZE + m];
					}
				// overlay
				if (Total_mesh > 1 && current_sub_ele.overlay[j])
				{
					// calc opp_para
					double *opp_para_ptr = current_sub_ele.opp_para_tilde.data() + j * info->DIMENSION;
					int overlay_ele = current_sub_ele.overlay_ele[j];
					vector<double> opp_bl(D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
					vector<double> opp_b_disp_grad(info->DIMENSION * info->DIMENSION * MAX_KIEL_SIZE, 0.0);
					vector<double> opp_u(info->DIMENSION * MAX_KIEL_SIZE, 0.0);

					Make_B_Linear(overlay_ele, opp_para_ptr, opp_bl.data(), info);
					{
						vector<double> opp_b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
						if (overlay_ele < info->Total_Element_to_mesh[1])
							Make_B_component(overlay_ele, opp_para_ptr, opp_b.data(), info);
						else
							Make_B_component_for_SSIGA(overlay_ele, opp_para_ptr, opp_b.data(), info);

						for (int cp = 0; cp < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; cp++)
							for (int comp = 0; comp < info->DIMENSION; comp++)
								for (int dir = 0; dir < info->DIMENSION; dir++)
								{
									int offset = dir * info->DIMENSION + comp;
									opp_b_disp_grad[offset * MAX_KIEL_SIZE + cp * info->DIMENSION + comp] = opp_b[dir * MAX_NO_CP_ON_ELEMENT + cp];
								}
					}

					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int id = info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
							opp_u[k * info->DIMENSION + l] = info->disp[id];
						}
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
							{
								strain_trial[k] += opp_bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * opp_u[l * info->DIMENSION + m];
							}
					for (int k = 0; k < info->DIMENSION; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int offset = k * info->DIMENSION + l;
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]] * info->DIMENSION; m++)
								disp_grad[k * info->DIMENSION + l] += opp_u[m] * opp_b_disp_grad[offset * MAX_KIEL_SIZE + m];
						}
				}

				// Calculate trial stresses: {sigma}^trial = [D] * {epsilon}^trial (small deformation)
				vector<double> stress_trial(D_MATRIX_SIZE, 0.0);
				for (int k = 0; k < D_MATRIX_SIZE; k++)
				for (int l = 0; l < D_MATRIX_SIZE; l++)
				stress_trial[k] += D_e[k * D_MATRIX_SIZE + l] * strain_trial[l];
				
				// 使わないかも
				// double equivalent_stress_trial = calc_equivalent_stress(stress_trial.data());

				// update current trial
				current_sub_ele.elastic_strain[j] = strain_trial;
				current_sub_ele.stress[j] = stress_trial;

				#if 0
				// elastoplastic correction (disabled for linear elastic analysis)
				if (equivalent_stress_trial > current_sub_ele.yield_stress[j])
				{
					// Calculate equivalent plastic strain increment
					double equivalent_plastic_strain_increment = calc_equivalent_plastic_strain_increment(equivalent_stress_trial, current_sub_ele.equivalent_plastic_strain[j], current_sub_ele.yield_stress[j], info);

					// Calculate hardening stress increment
					double hardening_stress_increment = get_hardening_stress(current_sub_ele.equivalent_plastic_strain[j] + equivalent_plastic_strain_increment, info) - get_hardening_stress(current_sub_ele.equivalent_plastic_strain[j], info);

					// Calculate trial relative deviatoric stresses
					double current_relative_hydrostatic_stress = (1.0 / 3.0) * (stress_trial[0] + stress_trial[1] + stress_trial[2]);
					for (int k = 0; k < info->DIMENSION; k++)
						stress_trial[k] -= current_relative_hydrostatic_stress;

					// Calculate final elastic strains
					double coef_1 = equivalent_plastic_strain_increment * 1.5 / equivalent_stress_trial;
					for (int k = 0; k < info->DIMENSION; k++)
						current_sub_ele.elastic_strain[j][k] -= coef_1 * stress_trial[k];
					for (int k = info->DIMENSION; k < D_MATRIX_SIZE; k++)
						current_sub_ele.elastic_strain[j][k] -= 2.0 * coef_1 * stress_trial[k];

					// Calculate final yield stress
					current_sub_ele.yield_stress[j] = current_sub_ele.yield_stress[j] + (1.0 - info->kinematic_hardening_fraction) * hardening_stress_increment;

					// Calculate final stresses
					double coef_2 = current_sub_ele.yield_stress[j] / equivalent_stress_trial;
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						current_sub_ele.stress[j][k] = coef_2 * stress_trial[k];
					for (int k = 0; k < info->DIMENSION; k++)
						current_sub_ele.stress[j][k] += current_relative_hydrostatic_stress;

					// update equivalent plastic strain
					current_sub_ele.equivalent_plastic_strain[j] += equivalent_plastic_strain_increment;
				}
				#endif

				// W (single-step energy density)
				double W = 0.0;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					W += 0.5 * stress_trial[k] * strain_trial[k];

				// (∂q/∂x_i) = ∂q/∂ξ^ * ∂ξ^/∂ξ * ∂ξ/∂x_i
				vector<double> q_grad(info->DIMENSION, 0.0);
				vector<double> q_xitildei(info->DIMENSION, 0.0);
				vector<double> a_inv(info->DIMENSION * info->DIMENSION, 0.0);

				#ifndef SQUARE_PATCH
				calc_q_gradient(e, para_ptr, current_sub_ele.slice, current_sub_ele.face, current_J_list, q_xitildei, info);
				#else
				calc_q_gradient_square_patch(e, para_ptr, current_sub_ele.slice, current_sub_ele.face, current_J_list, q_xitildei, info);
				#endif

				double jac;
				if (e < info->Total_Element_to_mesh[1])
				{
					jac = Make_a_q(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, delta, true);
					Make_a_q(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, 0.0, false);
				}
				else
				{
					jac = Make_a_q_for_SSIGA(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, delta, true);
					Make_a_q_for_SSIGA(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, 0.0, false);
				}

				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						q_grad[k] += q_xitildei[l] * a_inv[l * info->DIMENSION + k];
	
				double coef = current_sub_ele.w[j] * jac / J_integral_area[index_J];

				vector<double> W_I(info->DIMENSION * info->DIMENSION, 0.0);
				for (int k = 0; k < info->DIMENSION; k++)
					W_I[k * info->DIMENSION + k] = W;
	
				vector<double> stress_tensor(info->DIMENSION * info->DIMENSION, 0.0);
				stress_tensor[0] = current_sub_ele.stress[j][0];
				stress_tensor[1] = current_sub_ele.stress[j][3];
				stress_tensor[2] = current_sub_ele.stress[j][5];
				stress_tensor[3] = current_sub_ele.stress[j][3];
				stress_tensor[4] = current_sub_ele.stress[j][1];
				stress_tensor[5] = current_sub_ele.stress[j][4];
				stress_tensor[6] = current_sub_ele.stress[j][5];
				stress_tensor[7] = current_sub_ele.stress[j][4];
				stress_tensor[8] = current_sub_ele.stress[j][2];
	
				vector<double> stress_disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							// 変更前（正しい？）
							stress_disp_grad[k * info->DIMENSION + l] += stress_tensor[k * info->DIMENSION + m] * disp_grad[l * info->DIMENSION + m];

							// 変更後
							// stress_disp_grad[l * info->DIMENSION + k] += stress_tensor[k * info->DIMENSION + m] * disp_grad[m * info->DIMENSION + l];
						}

				vector<double> J_int(info->DIMENSION, 0.0);
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						J_int[k] += (stress_disp_grad[l * info->DIMENSION + k] - W_I[k * info->DIMENSION + l]) * q_grad[l] * coef;
					}
	
				for (int k = 0; k < info->DIMENSION; k++)
				{
					double val = J_int[k] * normal[index_J][k];
					local_J_val[current_sub_ele.jp_list_index][index_J] += val;
				}

				// interaction integral method
				#if defined(INTERACTION_INTEGRAL_METHOD)
				vector<double> aux_W_I(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> aux_disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> aux_strain_tensor(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> aux_stress_tensor(info->DIMENSION * info->DIMENSION, 0.0);

				vector<double> strain_tensor(info->DIMENSION * info->DIMENSION, 0.0);
				strain_tensor[0] = current_sub_ele.elastic_strain[j][0];
				strain_tensor[1] = 0.5 * current_sub_ele.elastic_strain[j][3];
				strain_tensor[2] = 0.5 * current_sub_ele.elastic_strain[j][5];
				strain_tensor[3] = 0.5 * current_sub_ele.elastic_strain[j][3];
				strain_tensor[4] = current_sub_ele.elastic_strain[j][1];
				strain_tensor[5] = 0.5 * current_sub_ele.elastic_strain[j][4];
				strain_tensor[6] = 0.5 * current_sub_ele.elastic_strain[j][5];
				strain_tensor[7] = 0.5 * current_sub_ele.elastic_strain[j][4];
				strain_tensor[8] = current_sub_ele.elastic_strain[j][2];

				const double *coord = jp_list[current_sub_ele.jp_list_index].coord[index_J].data();
				const double *current_sub_ele_coord = &current_sub_ele.coord[j * info->DIMENSION];
				const vector<vector<vector<double>>> &Q = jp_list[current_sub_ele.jp_list_index].Q;

				// calc three K modes
				double r = 0.0, theta = 0.0;

				calc_polar_coord_for_auxiliary_fields(r, theta, Q[index_J], coord, current_sub_ele_coord);
				for (int k = 0; k < 3; k++)
				{
					// calc aux stress tensor and aux disp grad for mode k
					calc_auxiliary_fields(aux_disp_grad, aux_strain_tensor, aux_stress_tensor, &r, &theta, Q[index_J], k, D_e);

					// aux W_I
					double aux_W = 0.0;
					for (int l = 0; l < info->DIMENSION * info->DIMENSION; l++)
						aux_W += 0.5 * (aux_stress_tensor[l] * strain_tensor[l] + stress_tensor[l] * aux_strain_tensor[l]);
					for (int l = 0; l < info->DIMENSION; l++)
						aux_W_I[l * info->DIMENSION + l] = aux_W;

					// second term: aux_stress_tensor * disp_grad
					vector<double> aux_stress_disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
					for (int kdir = 0; kdir < info->DIMENSION; kdir++)
						for (int ldir = 0; ldir < info->DIMENSION; ldir++)
							for (int mcomp = 0; mcomp < info->DIMENSION; mcomp++)
							{
								// aux_stress_disp_grad[kdir * info->DIMENSION + ldir] += aux_stress_tensor[kdir * info->DIMENSION + mcomp] * disp_grad[mcomp * info->DIMENSION + ldir];
								aux_stress_disp_grad[kdir * info->DIMENSION + ldir] += aux_stress_tensor[kdir * info->DIMENSION + mcomp] * disp_grad[ldir * info->DIMENSION + mcomp];
							}

					// third term: stress_tensor * aux_disp_grad
					vector<double> stress_aux_disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
					for (int kdir = 0; kdir < info->DIMENSION; kdir++)
						for (int ldir = 0; ldir < info->DIMENSION; ldir++)
							for (int mcomp = 0; mcomp < info->DIMENSION; mcomp++)
							{
								// stress_aux_disp_grad[kdir * info->DIMENSION + ldir] += stress_tensor[kdir * info->DIMENSION + mcomp] * aux_disp_grad[mcomp * info->DIMENSION + ldir];
								stress_aux_disp_grad[kdir * info->DIMENSION + ldir] += stress_tensor[kdir * info->DIMENSION + mcomp] * aux_disp_grad[ldir * info->DIMENSION + mcomp];
							}

					vector<double> IIM_int(info->DIMENSION, 0.0);
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							IIM_int[l] += (stress_aux_disp_grad[m * info->DIMENSION + l] + aux_stress_disp_grad[m * info->DIMENSION + l] - aux_W_I[l * info->DIMENSION + m]) * q_grad[m] * coef;
						}

					// add to IIM_K_val
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double val = IIM_int[l] * normal[index_J][l] * normal_projection_sign;
						local_IIM_K_val[current_sub_ele.jp_list_index][index_J][k] += val;
					}
				}
				#endif
			}
		}
		#pragma omp critical (J_val_merge_lock)
		{
			for (size_t i = 0; i < jp_list.size(); i++)
				for (size_t j = 0; j < jp_list[i].J_val.size(); j++)
					jp_list[i].J_val[j] += local_J_val[i][j];
			
			#if defined(INTERACTION_INTEGRAL_METHOD)
			for (size_t i = 0; i < jp_list.size(); i++)
				for (size_t j = 0; j < jp_list[i].IIM_K_val.size(); j++)
					for (size_t k = 0; k < jp_list[i].IIM_K_val[j].size(); k++)
						jp_list[i].IIM_K_val[j][k] += local_IIM_K_val[i][j][k];
			#endif
		}
	}
}


void normal_vector(J_info &J, double *para, int e, J_point_info &jp, vector<vector<double>> &normal, information *info)
{
	// calc by adopting cross product for tangent vectors
	vector<vector<double>> tangent_vector(2, vector<double>(info->DIMENSION, 0.0));

	double para_on_curve[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		if (i != J.crack_dir)
		{
			// para_on_curve[i] = -1.0;
			para_on_curve[i] = -0.99; // -1.0だと亀裂先端の要素では縮退されていて法線ベクトルを求められない
		}
		else
		{
			para_on_curve[i] = para[i];
			// para_on_curve[i] = 0.0;
		}
	}

	vector<double> a(info->DIMENSION * info->DIMENSION, 0.0);
	// Make_a(e, para, a.data(), info); // 要修正かも！！亀裂先端の要素では縮退されていてこの方法では法線ベクトルを求められていない
	Make_a(e /*e_on_curve*/, para_on_curve, a.data(), info); // 要修正かも！！亀裂先端の要素では縮退されていてこの方法では法線ベクトルを求められていない
	// Make_a(e+1 /*e_on_curve*/, para_on_curve, a.data(), info); // 要修正かも！！亀裂先端の要素では縮退されていてこの方法では法線ベクトルを求められていない
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			tangent_vector[i][j] = a[jp.tangent_dir[i] * info->DIMENSION + j];

	vector<double> temp_normal(info->DIMENSION, 0.0);
	vector_normalize(tangent_vector[0].data(), info->DIMENSION);
	vector_normalize(tangent_vector[1].data(), info->DIMENSION);
	calc_normalized_cross(temp_normal.data(), tangent_vector[0].data(), tangent_vector[1].data());

	for (int i = 0; i < info->DIMENSION; i++)
	{
		// normal[i][0] = temp_normal[i];
		// normal[i][1] = tangent_vector[0][i];
		// normal[i][2] = tangent_vector[1][i];
		normal[i][0] = temp_normal[i];
		normal[i][1] = tangent_vector[1][i];
		normal[i][2] = - tangent_vector[0][i];
	}

	// debug output
	#if 1
	printf("normal vector for element %d:\n", e);
	printf("para: ");
	for (int i = 0; i < info->DIMENSION; i++)
		printf("%e ", para[i]);
	printf("\n");
	for (int i = 0; i < info->DIMENSION; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
			printf("%e ", normal[i][j]);
		printf("\n");
	}
	exit(0);
	#endif
}


bool calc_normalized_cross(double *cross, double *a, double *b)
{
	cross[0] = a[1] * b[2] - a[2] * b[1];
	cross[1] = a[2] * b[0] - a[0] * b[2];
	cross[2] = a[0] * b[1] - a[1] * b[0];

	double norm = 0.0;
	for (int i = 0; i < 3; i++)
		norm += cross[i] * cross[i];
	norm = sqrt(norm);

	if (norm > MERGE_ERROR)
	{
		for (int i = 0; i < 3; i++)
			cross[i] /= norm;
		return true;
	}

	return false; // if a and b are parallel, return false
}


void vector_normalize(double *vec, int dim)
{
	double norm = 0.0;
	for (int i = 0; i < dim; i++)
		norm += vec[i] * vec[i];
	norm = sqrt(norm);
	if (norm > MERGE_ERROR) // avoid division by zero
	{
		for (int i = 0; i < dim; i++)
			vec[i] /= norm;
	}
	else
	{
		for (int i = 0; i < dim; i++)
			vec[i] = 0.0; // if norm is zero, set vector to zero
	}
}

void calc_q_gradient(int e, double *tilde_coord, int line, int face, J_info &J_list, vector<double> &q_grad, information *info)
{
	int r_dir = J_list.r_dir;
	int crack_dir = J_list.crack_dir;
	int patch = J_list.p;
	double *knot_ptr_r_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + r_dir] + info->Order[patch * info->DIMENSION + r_dir]];
	double *knot_ptr_crack_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + crack_dir] + info->Order[patch * info->DIMENSION + crack_dir]];

	vector<double> natural_coord(2, 0.0);

	// r_dir natural coord
	int current_knot_r_dir = info->ENC[e * info->DIMENSION + r_dir];
	natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir] + 0.5 * (tilde_coord[r_dir] + 1.0) * (knot_ptr_r_dir[current_knot_r_dir + 1] - knot_ptr_r_dir[current_knot_r_dir]);
	if (fabs(natural_coord[0] - knot_ptr_r_dir[current_knot_r_dir]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir];
	else if (fabs(natural_coord[0] - knot_ptr_r_dir[current_knot_r_dir + 1]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir + 1];

	// crack_dir natural coord
	int current_knot_crack_dir = info->ENC[e * info->DIMENSION + crack_dir];
	natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir] + 0.5 * (tilde_coord[crack_dir] + 1.0) * (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]);
	if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir];
	else if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir + 1]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir + 1];

	// cout << "natural_coord: " << natural_coord[0] << " " << natural_coord[1] << endl;

	vector<double> r_coord(2, 0.0);
	vector<double> crack_coord(2, 0.0);

	// r_coord
	int ele_n_r_dir = info->No_Control_point[patch * info->DIMENSION + r_dir] - info->Order[patch * info->DIMENSION + r_dir];
	r_coord[0] = knot_ptr_r_dir[0];
	r_coord[1] = knot_ptr_r_dir[ele_n_r_dir];

	// crack_coord
	// cout << "tilde_coord[crack_dir]: " << tilde_coord[crack_dir] << endl;
	double delta = (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]) / static_cast<double>(J_list.sub_n);
	vector<double> temp(2, 0.0);
	temp[0] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * line;
	temp[1] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * (line + 1);
	if (face == 1)
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[1];
	}
	else
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[0];
	}

	// cout << "r_coord: " << r_coord[0] << " " << r_coord[1] << endl;
	// cout << "crack_coord: " << crack_coord[0] << " " << crack_coord[1] << endl;

	// init
	for (int i = 0; i < info->DIMENSION; i++)
		q_grad[i] = 0.0;

	// r direction
	q_grad[r_dir] = (-1.0 / (r_coord[1] - r_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]));
	// q_grad[crack_dir] = (-1.0 / (r_coord[1] - r_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]));
	// q_grad[r_dir] = coef[0] * (-1.0 / (r_coord[1] - r_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]));
	// q_grad[crack_dir] = coef[0] * (-1.0 / (r_coord[1] - r_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]));
	// q_grad[r_dir] = (-1.0 / (r_coord[1] - r_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]));

	// crack direction
	q_grad[crack_dir] = (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));
	// q_grad[r_dir] = (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));
	// q_grad[crack_dir] = coef[1] * (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));;
	// q_grad[r_dir] = coef[1] * (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));;
	// q_grad[crack_dir] = (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));;

	// cout << "r_dir: " << r_dir << " crack_dir: " << crack_dir << " q_grad " << q_grad[0] << " " << q_grad[1] << " " << q_grad[2] << endl;
	// exit(0);

	if ((natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0]) < 0 || (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]) < 0)
	{
		printf("Error: negative gradient in calc_q_gradient\n");
		exit(0);
	}
}

void calc_q_gradient_square_patch(int e, double *tilde_coord, int line, int face, J_info &J_list, vector<double> &q_grad, information *info)
{
	int quadrant = J_list.quadrant;

	int r1_dir = J_list.r_dir;
	int crack_dir = J_list.crack_dir;
	int patch = J_list.p;
	double *knot_ptr_r1_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + r1_dir] + info->Order[patch * info->DIMENSION + r1_dir]];
	double *knot_ptr_crack_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + crack_dir] + info->Order[patch * info->DIMENSION + crack_dir]];

	// r2 direction
	int r2_dir = -1;
	double *knot_ptr_r2_dir = nullptr;
	for (int i = 0; i < info->DIMENSION; i++)
		if (i != r1_dir && i != crack_dir)
		{
			r2_dir = i;
			knot_ptr_r2_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + r2_dir] + info->Order[patch * info->DIMENSION + r2_dir]];
			break;
		}

	vector<double> natural_coord(3, 0.0);

	// r1_dir natural coord
	int current_knot_r1_dir = info->ENC[e * info->DIMENSION + r1_dir];
	natural_coord[0] = knot_ptr_r1_dir[current_knot_r1_dir] + 0.5 * (tilde_coord[r1_dir] + 1.0) * (knot_ptr_r1_dir[current_knot_r1_dir + 1] - knot_ptr_r1_dir[current_knot_r1_dir]);
	if (fabs(natural_coord[0] - knot_ptr_r1_dir[current_knot_r1_dir]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r1_dir[current_knot_r1_dir];
	else if (fabs(natural_coord[0] - knot_ptr_r1_dir[current_knot_r1_dir + 1]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r1_dir[current_knot_r1_dir + 1];

	// crack_dir natural coord
	int current_knot_crack_dir = info->ENC[e * info->DIMENSION + crack_dir];
	natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir] + 0.5 * (tilde_coord[crack_dir] + 1.0) * (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]);
	if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir];
	else if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir + 1]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir + 1];

	// r2_dir natural coord
	int current_knot_r2_dir = info->ENC[e * info->DIMENSION + r2_dir];
	natural_coord[2] = knot_ptr_r2_dir[current_knot_r2_dir] + 0.5 * (tilde_coord[r2_dir] + 1.0) * (knot_ptr_r2_dir[current_knot_r2_dir + 1] - knot_ptr_r2_dir[current_knot_r2_dir]);
	if (fabs(natural_coord[2] - knot_ptr_r2_dir[current_knot_r2_dir]) < MERGE_ERROR)
		natural_coord[2] = knot_ptr_r2_dir[current_knot_r2_dir];
	else if (fabs(natural_coord[2] - knot_ptr_r2_dir[current_knot_r2_dir + 1]) < MERGE_ERROR)
		natural_coord[2] = knot_ptr_r2_dir[current_knot_r2_dir + 1];

	vector<double> coef_vec(info->DIMENSION, 1.0);
	if (quadrant == 1)
	{
		natural_coord[0] = natural_coord[0]; // r1_dir
		natural_coord[1] = natural_coord[1]; // crack_dir
		natural_coord[2] = natural_coord[2]; // r2_dir
	}
	else if (quadrant == 2)
	{
		coef_vec[r1_dir] = -1.0;
		natural_coord[0] = knot_ptr_r1_dir[info->No_Control_point[patch * info->DIMENSION + r1_dir] - info->Order[patch * info->DIMENSION + r1_dir]] - natural_coord[0]; // r1_dir
		natural_coord[1] = natural_coord[1]; // crack_dir
		natural_coord[2] = natural_coord[2]; // r2_dir
	}
	else if (quadrant == 3)
	{
		coef_vec[r1_dir] = -1.0;
		coef_vec[r2_dir] = -1.0;
		natural_coord[0] = knot_ptr_r1_dir[info->No_Control_point[patch * info->DIMENSION + r1_dir] - info->Order[patch * info->DIMENSION + r1_dir]] - natural_coord[0]; // r1_dir
		natural_coord[1] = natural_coord[1]; // crack_dir
		natural_coord[2] = knot_ptr_r2_dir[info->No_Control_point[patch * info->DIMENSION + r2_dir] - info->Order[patch * info->DIMENSION + r2_dir]] - natural_coord[2]; // r2_dir
	}
	else if (quadrant == 4)
	{
		coef_vec[r2_dir] = -1.0;
		natural_coord[0] = natural_coord[0]; // r1_dir
		natural_coord[1] = natural_coord[1]; // crack_dir
		natural_coord[2] = knot_ptr_r2_dir[info->No_Control_point[patch * info->DIMENSION + r2_dir] - info->Order[patch * info->DIMENSION + r2_dir]] - natural_coord[2]; // r2_dir
	}

	// cout << "natural_coord: " << natural_coord[0] << " " << natural_coord[1] << endl;

	vector<double> r1_coord(2, 0.0);
	vector<double> crack_coord(2, 0.0);
	vector<double> r2_coord(2, 0.0);

	// r1_coord
	int ele_n_r1_dir = info->No_Control_point[patch * info->DIMENSION + r1_dir] - info->Order[patch * info->DIMENSION + r1_dir];
	r1_coord[0] = knot_ptr_r1_dir[0];
	r1_coord[1] = knot_ptr_r1_dir[ele_n_r1_dir];

	// crack_coord
	// cout << "tilde_coord[crack_dir]: " << tilde_coord[crack_dir] << endl;
	double delta = (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]) / static_cast<double>(J_list.sub_n);
	vector<double> temp(2, 0.0);
	temp[0] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * line;
	temp[1] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * (line + 1);
	if (face == 1)
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[1];
	}
	else
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[0];
	}

	// r2_coord
	int ele_n_r2_dir = info->No_Control_point[patch * info->DIMENSION + r2_dir] - info->Order[patch * info->DIMENSION + r2_dir];
	r2_coord[0] = knot_ptr_r2_dir[0];
	r2_coord[1] = knot_ptr_r2_dir[ele_n_r2_dir];

	// cout << "r1_coord: " << r1_coord[0] << " " << r1_coord[1] << endl;
	// cout << "r2_coord: " << r2_coord[0] << " " << r2_coord[1] << endl;
	// cout << "crack_coord: " << crack_coord[0] << " " << crack_coord[1] << endl;

	// init
	for (int i = 0; i < info->DIMENSION; i++)
		q_grad[i] = 0.0;

	// r direction
	q_grad[r1_dir] = (-1.0 / (r1_coord[1] - r1_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[2] - r2_coord[0]) / (r2_coord[1] - r2_coord[0]));

	// crack direction
	q_grad[crack_dir] = (1.0 - (natural_coord[0] - r1_coord[0]) / (r1_coord[1] - r1_coord[0])) * (-1.0 / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[2] - r2_coord[0]) / (r2_coord[1] - r2_coord[0]));

	// r2 direction
	q_grad[r2_dir] = (1.0 - (natural_coord[0] - r1_coord[0]) / (r1_coord[1] - r1_coord[0])) * (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0])) * (-1.0 / (r2_coord[1] - r2_coord[0]));

	for (int i = 0; i < info->DIMENSION; i++)
		q_grad[i] *= coef_vec[i];
}


double calc_virtual_crack_extension_area(J_info &J, int e, const int face_n, int sub_num, information *info)
{
	const int gp_1d = info->c.NUM_GAUSS_POINTS;
	double delta = 1.0 / (static_cast<double>(J.sub_n) * static_cast<double>(face_n));

	double virtual_crack_extension_area = 0.0;

	// calc length of element
	// double length = 0.0;
	for (int i = 0; i < face_n; i++)
	{
		for (int j = 0; j < gp_1d; j++)
		{
			double para_on_curve[MAX_DIMENSION];
			for (int k = 0; k < info->DIMENSION; k++)
			{
				if (k != J.crack_dir)
					para_on_curve[k] = -1.0;
				else
				{
					double normalized_local_para = (info->gauss_point_1D[j] + 1.0) / 2.0;
					double normalized_para = (static_cast<double>(sub_num * face_n + i) + normalized_local_para) * delta;
					para_on_curve[k] = normalized_para * 2.0 - 1.0;
				}
			}

			vector<double> jac_mat(info->DIMENSION * info->DIMENSION, 0.0);
			Make_a_tilde_prime(e, para_on_curve, jac_mat.data(), info, J.crack_dir, delta);

			vector<double> x_xi(info->DIMENSION, 0.0);
			for (int k = 0; k < info->DIMENSION; k++)
			{
				// x_xi[k] = jac_mat[J.crack_dir * info->DIMENSION + k];
				x_xi[k] = jac_mat[k * info->DIMENSION + J.crack_dir];
			}

			double val = 0.0;
			for (int k = 0; k < info->DIMENSION; k++)
				val += x_xi[k] * x_xi[k];
			val = sqrt(val);
			// length += val * info->gauss_w_1D[j];

			virtual_crack_extension_area += val * info->gauss_w_1D[j] * calc_q_1D(info->gauss_point_1D[j], i);
			// cout << calc_q_1D(info->gauss_point_1D[j], i) << " " << info->gauss_point_1D[j] << endl;
		}
		// cout << endl;
	}

	// calc subdivision area VCE [ΔA = Δ * q / 2 (q = 1)]
	#if 0
	double VCE = length * 0.5;
	#else
	double VCE = virtual_crack_extension_area;
	#endif

	#if 0
	cout << "length: " << length << " VCE: " << VCE << endl;
	cout << "virtual_crack_extension_area: " << virtual_crack_extension_area << endl;
	exit(0);
	#endif

	#if 0
	static vector<int> e_list;
	static vector<int> sub_list;
	// if not in both e_list and sub_list
	static std::vector<std::pair<int, int>> e_sub_list;
	auto it = std::find(e_sub_list.begin(), e_sub_list.end(), std::make_pair(e, sub_num));
	if (it == e_sub_list.end())
	{
		e_sub_list.emplace_back(e, sub_num);
		cout << "\033[33m";
		cout << "e: " << e << " sub_num: " << sub_num << " length: " << length << endl;
		cout << "\033[0m";
	}
	#endif

	return VCE;
}


double calc_q_1D(double normalized_para_1D, int face)
{
	double current_para = (normalized_para_1D + 1.0) / 2.0; // convert to [0, 1] range

	if (face == 1)
	{
		return (1.0 - current_para);
	}
	else
	{
		return current_para;
	}

	return 1.0;
}


double calc_q(int e, vector<double> &tilde_coord, int line, int face, J_info &J_list, information *info)
{
	int r_dir = J_list.r_dir;
	int crack_dir = J_list.crack_dir;
	int patch = J_list.p;
	double *knot_ptr_r_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + r_dir] + info->Order[patch * info->DIMENSION + r_dir]];
	double *knot_ptr_crack_dir = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + crack_dir] + info->Order[patch * info->DIMENSION + crack_dir]];

	vector<double> natural_coord(2, 0.0);

	// r_dir natural coord
	int current_knot_r_dir = info->ENC[e * info->DIMENSION + r_dir];
	natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir] + 0.5 * (tilde_coord[r_dir] + 1.0) * (knot_ptr_r_dir[current_knot_r_dir + 1] - knot_ptr_r_dir[current_knot_r_dir]);
	if (fabs(natural_coord[0] - knot_ptr_r_dir[current_knot_r_dir]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir];
	else if (fabs(natural_coord[0] - knot_ptr_r_dir[current_knot_r_dir + 1]) < MERGE_ERROR)
		natural_coord[0] = knot_ptr_r_dir[current_knot_r_dir + 1];

	// crack_dir natural coord
	int current_knot_crack_dir = info->ENC[e * info->DIMENSION + crack_dir];
	natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir] + 0.5 * (tilde_coord[crack_dir] + 1.0) * (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]);
	if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir];
	else if (fabs(natural_coord[1] - knot_ptr_crack_dir[current_knot_crack_dir + 1]) < MERGE_ERROR)
		natural_coord[1] = knot_ptr_crack_dir[current_knot_crack_dir + 1];

	vector<double> r_coord(2, 0.0);
	vector<double> crack_coord(2, 0.0);

	// r_coord
	int ele_n_r_dir = info->No_Control_point[patch * info->DIMENSION + r_dir] - info->Order[patch * info->DIMENSION + r_dir];
	r_coord[0] = knot_ptr_r_dir[0];
	r_coord[1] = knot_ptr_r_dir[ele_n_r_dir];

	// crack_coord
	double delta = (knot_ptr_crack_dir[current_knot_crack_dir + 1] - knot_ptr_crack_dir[current_knot_crack_dir]) / static_cast<double>(J_list.sub_n);
	vector<double> temp(2, 0.0);
	temp[0] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * line;
	temp[1] = knot_ptr_crack_dir[current_knot_crack_dir] + delta * (line + 1);
	if (face == 1)
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[1];
	}
	else
	{
		crack_coord[0] = 0.5 * (temp[0] + temp[1]);
		crack_coord[1] = temp[0];
	}

	// r direction
	double q = (1.0 - (natural_coord[1] - crack_coord[0]) / (crack_coord[1] - crack_coord[0])) * (1.0 - (natural_coord[0] - r_coord[0]) / (r_coord[1] - r_coord[0]));

	return q;
}


// Stress
double calc_equivalent_stress(double *stress)
{
	double a = (stress[0] - stress[1]) * (stress[0] - stress[1]);
	double b = (stress[1] - stress[2]) * (stress[1] - stress[2]);
	double c = (stress[2] - stress[0]) * (stress[2] - stress[0]);
	double d = stress[3] * stress[3];
	double e = stress[4] * stress[4];
	double f = stress[5] * stress[5];

	return sqrt(0.5 * (a + b + c + 6.0 * (d + e + f)));
}


double Make_a_q(int ele, double *para, double *a_inv, information *info, int crack_dir, double delta, bool coef_flag)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, coef_flag, info);
	vector<double> coef(info->DIMENSION, 1.0);
	if (coef_flag)
		coef[crack_dir] = delta;

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif
					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}	
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a_inv[i * info->DIMENSION + j] = a_2x2[i * info->DIMENSION + j];

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif
					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a_inv[i * info->DIMENSION + j] = a_3x3[i * info->DIMENSION + j];

		return J;
	}

	return J;
}


double Make_a_q_for_SSIGA(int ele_disp, double *para_disp, double *a_inv, information *info, int crack_dir, double delta, bool coef_flag)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];
	double b_2x2[4], b_3x3[9];
	double c_2x2[4], c_3x3[9];

	// パラメータ座標再定義　(para_dispはローカル変位要素のパラメータ座標)
	vector<double> para_geo(info->DIMENSION); // ローカル形状要素パラメータ座標
	vector<double> para_glo(info->DIMENSION); // グローバル要素パラメータ座標

	// 要素番号の取得と座標系の変換
	int ele_geo = trans_local_para_to_local_geo_para(ele_disp, para_disp, para_geo.data(), info);
	int ele_glo = trans_local_para_to_global_para(ele_disp, para_disp, para_glo.data(), info);

	vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
	vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR_geo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	vector<double> dR_glo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	geo_shape_and_dshape(R_geo.data(), dR_geo.data(), para_geo.data(), ele_geo, false, info);
	shape_and_dshape(R_glo.data(), dR_glo.data(), para_glo.data(), ele_glo, false, info);

	// 変位要素空間の微分によって生じる定数
	if (coef_flag)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			double coef = dShapeFunc_from_paren(i, ele_disp, info);
			for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; j++)
				dR_geo[j * info->DIMENSION + i] *= coef;
		}
	}

	// サブ要素空間の微分によって生じる定数
	vector<double> coef(info->DIMENSION, 1.0);
	if (coef_flag)
		coef[crack_dir] = delta;
	for (int i = 0; i < info->DIMENSION; i++)
	{
		for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; j++)
			dR_geo[j * info->DIMENSION + i] *= coef[i];
	}

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; k++)
				{
					int id0 = info->Geo_Controlpoint_of_Element[ele_geo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_2x2[i * 2 + j] += dR_geo[k * info->DIMENSION + j] * info->Geo_Node_Coordinate[id0];
				}
			}
		}
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				b_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_glo]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele_glo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					b_2x2[i * 2 + j] += dR_glo[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
				}	
			}
		}
        for (int i = 0; i < 2; i++)
        {
			for (int j = 0; j < 2; j++)
            {
				c_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < 2; k++)
                {
                    c_2x2[i * 2 + j] += b_2x2[i * 2 + k] * a_2x2[k * 2 + j];
                }
            }
        }

		J = InverseMatrix_2x2(c_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a_inv[i * info->DIMENSION + j] = c_2x2[i * info->DIMENSION + j];

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; k++)
				{
					int id0 = info->Geo_Controlpoint_of_Element[ele_geo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_3x3[i * 3 + j] += dR_geo[k * info->DIMENSION + j] * info->Geo_Node_Coordinate[id0];
				}
			}
		}
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				b_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_glo]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele_glo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					b_3x3[i * 3 + j] += dR_glo[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
				}
			}
		}
        for (int i = 0; i < 3; i++)
        {
			for (int j = 0; j < 3; j++)
            {
				c_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < 3; k++)
                {
                    c_3x3[i * 3 + j] += b_3x3[i * 3 + k] * a_3x3[k * 3 + j];
                }
            }
        }

		J = InverseMatrix_3x3(c_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a_inv[i * info->DIMENSION + j] = c_3x3[i * info->DIMENSION + j];

		return J;
	}

	return J;
}


void Make_a(int ele, double *para, double *a, information *info)
{
	double a_2x2[4], a_3x3[9];

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, false, info);
	// shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif

					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a[i * info->DIMENSION + j] = a_2x2[j * info->DIMENSION + i];
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif

					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a[i * info->DIMENSION + j] = a_3x3[j * info->DIMENSION + i];
	}
}


void Make_a_tilde_prime(int ele, double *para, double *a, information *info, int crack_dir, double delta)
{
	double a_2x2[4], a_3x3[9];

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	// shape_and_dshape(R.data(), dR.data(), para, ele, false, info);
	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

	vector<double> coef(info->DIMENSION, 1.0);
	coef[crack_dir] = delta;

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif

					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a[i * info->DIMENSION + j] = a_2x2[i * info->DIMENSION + j];
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					
					a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					
					#if 0
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
					if (Total_mesh < 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					else if (Total_mesh >= 2)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					#endif

					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a[i * info->DIMENSION + j] = a_3x3[i * info->DIMENSION + j];
	}
}


