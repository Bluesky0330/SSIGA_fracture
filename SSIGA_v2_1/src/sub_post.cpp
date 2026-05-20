// header
#include "_header.hpp"
#include "_sub.hpp"
// HDF5 C++ API
#include <H5Cpp.h>
using namespace H5;

using namespace std;

// Postprocessing
void Make_Displacement(information *info)
{
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
		info->Displacement[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = info->Value_of_Constraint[i];

	#pragma omp parallel for
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
		{
			int index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
				info->Displacement[i * info->DIMENSION + j] = info->sol_vec[index];
		}
	}
}


void Output_Displacement(information *info)
{
	FILE *disp_fp;
	disp_fp = fopen("__Displacement.dat", "w");
	fprintf(disp_fp, "%d %d\n", info->Total_Control_Point_to_mesh[Total_mesh], info->DIMENSION);
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
		{
			if (j == 0)
				fprintf(disp_fp, "% .20e", info->Displacement[i * info->DIMENSION + j]);
			else
				fprintf(disp_fp, " % .20e", info->Displacement[i * info->DIMENSION + j]);
		}
		fprintf(disp_fp, "\n");
	}
}


void Calc_on_Element_Vertex(information *info)
{	
	FILE *fp_disp, *fp_strain, *fp_stress, *fp_pc;
	fp_disp   = fopen("_Displacement_overlay_at_ele_vertex.dat", "w");
	fp_strain = fopen("_Strain_overlay_at_ele_vertex.dat", "w");
	fp_stress = fopen("_Stress_overlay_at_ele_vertex.dat", "w");
	fp_pc     = fopen("_PhysicalCoordinate_at_ele_vertex.dat", "w");
	if (info->DIMENSION == 2)
	{
		fprintf(fp_disp,   "要素番号\t頂点番号\tdisp_x\tdisp_y\n");
		fprintf(fp_strain, "要素番号\t頂点番号\txx\tyy\txy\tzz\n");
		fprintf(fp_stress, "要素番号\t頂点番号\txx\tyy\txy\tzz\n");
		fprintf(fp_pc,     "要素番号\t頂点番号\tx\ty\n");
	}
	else if (info->DIMENSION == 3)
	{
		fprintf(fp_disp,   "要素番号\t頂点番号\tdisp_x\tdisp_y\tdisp_z\n");
		fprintf(fp_strain, "要素番号\t頂点番号\txx\tyy\tzz\txy\tyz\txz\n");
		fprintf(fp_stress, "要素番号\t頂点番号\txx\tyy\tzz\txy\tyz\txz\n");
		fprintf(fp_pc,     "要素番号\t頂点番号\tx\ty\tz\n");
	}

	int vertex_n = pow_int(2, info->DIMENSION);
	vector<double> point_array(vertex_n * info->DIMENSION);

	int counter = 0;
	if (info->DIMENSION == 2)
	{
		point_array[counter] = -1.0;	point_array[counter + 1] = -1.0;	counter += 2;
		point_array[counter] =  1.0;	point_array[counter + 1] = -1.0;	counter += 2;
		point_array[counter] = -1.0;	point_array[counter + 1] =  1.0;	counter += 2;
		point_array[counter] =  1.0;	point_array[counter + 1] =  1.0;
	}
	else if (info->DIMENSION == 3)
	{
		point_array[counter] = -1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] = -1.0;	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] = -1.0;	counter += 3;
		point_array[counter] = -1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] = -1.0;	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] = -1.0;	counter += 3;

		point_array[counter] = -1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] = 	1.0;	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] = 	1.0; 	counter += 3;
		point_array[counter] = -1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] = 	1.0; 	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] = 	1.0;
	}

	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < vertex_n; j++)
		{
			double temp_point[MAX_DIMENSION] = {0.0};
			double temp_point_glo[MAX_DIMENSION] = {0.0};
			double temp_para_glo[MAX_DIMENSION] = {0.0};
			double temp_point_loc[MAX_DIMENSION] = {0.0};
			double temp_para_loc[MAX_DIMENSION] = {0.0};

			int element = i;
			int point = j;

			// for SSIGA local
			double geo_coord_tilde[MAX_DIMENSION] = {0.0};	//ローカルジオメトリ表現の要素パラメータ座標
			int geo_element = 0;							//ローカルジオメトリ表現の要素番号

			// make temp_point
			for (int j = 0; j < info->DIMENSION; j++)
				temp_point[j] = point_array[point * info->DIMENSION + j];
			
			if (i < info->Total_Element_on_mesh[0]) // IGA or SSIGA global
			{
				vector<double> R(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0);
				shape_and_dshape(R.data(), temp_point, element, info);
				Make_B_Linear(element, temp_point, bl.data(), info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->PhysicalCoordinate_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] += R[k] * info->Node_Coordinate[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
						info->Displacement_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] += R[k] * d;
						u[k * info->DIMENSION + l] = d;
					}
				// strain
				for (int k= 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
			}
			if (i >= info->Total_Element_on_mesh[0] && Total_mesh > 1)
			{
				int patch_num = info->Element_patch[element];							//要素番号からパッチ番号を取得
				int geo_patch_num = patch_num - info->Total_Patch_on_mesh[0];			//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい
				double coord[MAX_DIMENSION] = {0.0};
				double out_coord[MAX_DIMENSION] = {0.0};									//ローカルのパッチパラメータ座標
				double coord_tilde[MAX_DIMENSION] = {0.0};								//ローカルの要素パラメータ座標
				double glo_coord_tilde[MAX_DIMENSION] = {0.0};							//グローバルの要素パラメータ座標
				double Parameter_coord_at_ele_vertex[MAX_DIMENSION] = {0.0};			//要素のパラメータ座標	
			
				std::copy(temp_point, temp_point + info->DIMENSION, coord_tilde);				//要素のパラメータ座標をガウス点座標に変換
				trans_ele_patch_coord(coord, coord_tilde, patch_num, element, info);			//要素のパラメータ座標をパッチのパラメータ座標に変換
				geo_element = geo_ele_check(geo_patch_num, coord, info);						//パラメータ座標から、要素を探索
				geo_tilde_coord(geo_coord_tilde, coord, geo_patch_num, geo_element, info);		//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換
			
				vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
				vector<double> R_disp(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
				geo_shape_and_dshape(R_geo.data(), geo_coord_tilde, geo_element, info);
				Bspline_shape_and_dshape(R_disp.data(), temp_point, element, info);
				Make_B_Linear(element, temp_point, bl.data(), info);
			
				// physical coordinate
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[geo_element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						Parameter_coord_at_ele_vertex[l] += R_geo[k] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_element * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
				int ele_glo = ele_check(info->Global_local_patch, Parameter_coord_at_ele_vertex, info);
				tilde_coord(glo_coord_tilde, Parameter_coord_at_ele_vertex, info->Global_local_patch, ele_glo, info);
				physical_coord(ele_glo, glo_coord_tilde, out_coord, info);
				for (int l = 0; l < info->DIMENSION; l++)
					info->PhysicalCoordinate_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] = out_coord[l];
			
				// displacement
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] += R_disp[k] * d;
						u[k * info->DIMENSION + l] = d;
					}
				
				// strain
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
			}
			// overlay global
			int status_glo_overlay = 0;
			if (i < info->Total_Element_on_mesh[0] && Total_mesh > 1 && info->Element_patch[element] == info->Global_local_patch)
			{
				double temp_point_patch[MAX_DIMENSION] = {0.0};
				trans_ele_patch_coord(temp_point_patch, temp_point, info->Element_patch[element], element, info);

				// make temp_point_glo
				int itr_n = 0, loc_patch = 0;
				for (int k = info->Total_Patch_to_mesh[1]; k < info->Total_Patch_to_mesh[Total_mesh]; k++)
				{
					itr_n = calc_local_patch_parameter_coord(temp_point_patch, k, temp_para_loc, info);
					loc_patch = k;
					if (itr_n != ERROR)
					{
						status_glo_overlay = 1;
						break;
					}
				}

				int element_loc = 0;
				if (status_glo_overlay)
				{
					element_loc = ele_check(loc_patch, temp_para_loc, info);
					tilde_coord(temp_point_loc, temp_para_loc, loc_patch, element_loc, info);
				}

				// overlay displacement
				if (status_glo_overlay)
				{
					vector<double> R_loc(MAX_NO_CP_ON_ELEMENT);
					vector<double> bl_loc(D_MATRIX_SIZE * MAX_KIEL_SIZE);
					vector<double> u_loc(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
					Bspline_shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
					Make_B_Linear(element_loc, temp_point_loc, bl_loc.data(), info);
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							double d_loc = info->Displacement[info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							info->Displacement_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] += R_loc[k] * d_loc;
							u_loc[k * info->DIMENSION + l] = d_loc;
						}

					// strain
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
								info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k] += bl_loc[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u_loc[l * info->DIMENSION + m];
				}
			}

			// overlay local
			if (i >= info->Total_Element_on_mesh[0] && Total_mesh > 1)
			{
				geo_parameter_coord(geo_element, geo_coord_tilde, temp_para_glo, info);
				int element_glo = ele_check(info->Global_local_patch, temp_para_glo, info);
				tilde_coord(temp_point_glo, temp_para_glo, info->Global_local_patch, element_glo, info);

				// overlay displacement
				vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl_glo(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u_glo(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
				shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
				Make_B_Linear(element_glo, temp_point_glo, bl_glo.data(), info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_ele_vertex[element * vertex_n * info->DIMENSION + point * info->DIMENSION + l] += R_glo[k] * d_glo;
						u_glo[k * info->DIMENSION + l] = d_glo;
					}
				
				// strain
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k] += bl_glo[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u_glo[l * info->DIMENSION + m];
			}

			// stress
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				for (int l = 0; l < D_MATRIX_SIZE; l++)
					info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k] += info->D[k * D_MATRIX_SIZE + l] * info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRESS + l];

			// 2D specific adjustments
			if (info->DIMENSION == 2)
			{
				// plane stress condition
				if (info->c.ANALYSIS_MODE == 0)
				{
					info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 3] = -1.0 * nu / E * (info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 0] + info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 1]);
					info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 3] = 0.0;
				}

				// plane strain condition
				else if (info->c.ANALYSIS_MODE == 1)
				{
					info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 3] = 0.0;
					info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 3] = nu * (info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 0] + info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + 1]);
				}
			}

			// output file
			// displacement
			fprintf(fp_disp, "%d\t%d", element, point);
			for (int k = 0; k < info->DIMENSION; k++)
			fprintf(fp_disp, "\t%.15e", info->Displacement_at_ele_vertex[(element * vertex_n + point) * info->DIMENSION + k]);
			fprintf(fp_disp, "\n");
			// printf("wrote displacement\n");
			
			// strain
			fprintf(fp_strain, "%d\t%d", element, point);
			for (int k = 0; k < N_STRAIN; k++)
			fprintf(fp_strain, "\t%.15e", info->Strain_at_ele_vertex[(element * vertex_n + point) * N_STRAIN + k]);
			fprintf(fp_strain, "\n");
			// printf("wrote strain\n");
			
			// stress
			fprintf	(fp_stress, "%d\t%d", element, point);
			for (int k = 0; k < N_STRESS; k++)
			fprintf(fp_stress, "\t%.15e", info->Stress_at_ele_vertex[(element * vertex_n + point) * N_STRESS + k]);
			fprintf(fp_stress, "\n");
			// printf("wrote stress\n");

			// physical coordinate
			fprintf(fp_pc, "%d\t%d", element, point);
			for (int k = 0; k < info->DIMENSION; k++)
			fprintf(fp_pc, "\t%.15e", info->PhysicalCoordinate_at_ele_vertex[(element * vertex_n + point) * info->DIMENSION + k]);
			fprintf(fp_pc, "\n");
			// printf("wrote physical coordinate\n");
		}
	}

	if (fp_disp) 
	    fclose(fp_disp);

	if (fp_strain)
	    fclose(fp_strain);

	if (fp_stress)
	    fclose(fp_stress);

	if (fp_pc)
	    fclose(fp_pc);
}


void Calc_on_Gauss_Point(information *info)
{	
	const int n_gp = pow_int(info->c.NUM_GAUSS_POINTS, info->DIMENSION);
	const int n_gp_ex = pow_int(info->c.NUM_GAUSS_POINTS_EXTENDED, info->DIMENSION);
	const int gp_stride = (n_gp_ex > n_gp) ? n_gp_ex : n_gp;

	FILE *fp_disp, *fp_strain, *fp_stress, *fp_rf, *fp_pc;
	fp_disp   = fopen("_Displacement_overlay_at_GP.dat", "w");
	fp_strain = fopen("_Strain_overlay_at_GP.dat", "w");
	fp_stress = fopen("_Stress_overlay_at_GP.dat", "w");
	fp_rf     = fopen("_ReactionForce.dat", "w");
	fp_pc     = fopen("_PhysicalCoordinate_at_GP.dat", "w");
	if (info->DIMENSION == 2)
	{
		fprintf(fp_disp,   "要素番号\tガウス点番号\tdisp_x\tdisp_y\n");
		fprintf(fp_strain, "要素番号\tガウス点番号\txx\tyy\txy\tzz\n");
		fprintf(fp_stress, "要素番号\tガウス点番号\txx\tyy\txy\tzz\n");
		fprintf(fp_rf,     "コントロールポイント番号\trf_x\trf_y\n");
		fprintf(fp_pc,     "要素番号\tガウス点番号\tx\ty\n");
	}
	else if (info->DIMENSION == 3)
	{
		fprintf(fp_disp,   "要素番号\tガウス点番号\tdisp_x\tdisp_y\tdisp_z\n");
		fprintf(fp_strain, "要素番号\tガウス点番号\txx\tyy\tzz\txy\tyz\txz\n");
		fprintf(fp_stress, "要素番号\tガウス点番号\txx\tyy\tzz\txy\tyz\txz\n");
		fprintf(fp_rf,     "コントロールポイント番号\trf_x\trf_y\trf_z\n");
		fprintf(fp_pc,     "要素番号\tガウス点番号\tx\ty\tz\n");
	}

	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < info->gp[i].n(); j++)
		{
			double temp_point[MAX_DIMENSION] = {0.0};
			double temp_point_glo[MAX_DIMENSION] = {0.0};
			double temp_para_glo[MAX_DIMENSION] = {0.0};
			double temp_point_loc[MAX_DIMENSION] = {0.0};
			double temp_para_loc[MAX_DIMENSION] = {0.0};

			int element = i;
			int point = j;
			int gp_base = element * gp_stride + point;

			// for SSIGA local
			double geo_coord_tilde[MAX_DIMENSION] = {0.0};	//ローカルジオメトリ表現の要素パラメータ座標
			int geo_element = 0;							//ローカルジオメトリ表現の要素番号

			// make temp_point
			for (int k = 0; k < info->DIMENSION; k++)
				temp_point[k] = info->gp[element].para()[point * info->DIMENSION + k];
			
			if (i < info->Total_Element_on_mesh[0]) // IGA or SSIGA global
			{
				vector<double> R(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0);
				shape_and_dshape(R.data(), temp_point, element, info);
				Make_B_Linear(element, temp_point, bl.data(), info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->PhysicalCoordinate_at_GP[gp_base * info->DIMENSION + l] += R[k] * info->Node_Coordinate[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
						info->Displacement_at_GP[gp_base * info->DIMENSION + l] += R[k] * d;
						u[k * info->DIMENSION + l] = d;
					}
				// strain
				for (int k= 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_GP[gp_base * N_STRESS + k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
			}
			if (i >= info->Total_Element_on_mesh[0] && Total_mesh > 1)
			{
				int patch_num = info->Element_patch[element];							//要素番号からパッチ番号を取得
				int geo_patch_num = patch_num - info->Total_Patch_on_mesh[0];			//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい
				double coord[MAX_DIMENSION] = {0.0};
				double out_coord[MAX_DIMENSION] = {0.0};									//ローカルのパッチパラメータ座標
				double coord_tilde[MAX_DIMENSION] = {0.0};								//ローカルの要素パラメータ座標
				double glo_coord_tilde[MAX_DIMENSION] = {0.0};							//グローバルの要素パラメータ座標
				double Parameter_coord_at_GP[MAX_DIMENSION] = {0.0};			//要素のパラメータ座標	
			
				std::copy(temp_point, temp_point + info->DIMENSION, coord_tilde);				//要素のパラメータ座標をガウス点座標に変換
				trans_ele_patch_coord(coord, coord_tilde, patch_num, element, info);			//要素のパラメータ座標をパッチのパラメータ座標に変換
				geo_element = geo_ele_check(geo_patch_num, coord, info);						//パラメータ座標から、要素を探索
				geo_tilde_coord(geo_coord_tilde, coord, geo_patch_num, geo_element, info);		//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換
			
				vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
				vector<double> R_disp(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
				geo_shape_and_dshape(R_geo.data(), geo_coord_tilde, geo_element, info);
				Bspline_shape_and_dshape(R_disp.data(), temp_point, element, info);
				Make_B_Linear(element, temp_point, bl.data(), info);
			
				// physical coordinate
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[geo_element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						Parameter_coord_at_GP[l] += R_geo[k] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_element * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
				int ele_glo = ele_check(info->Global_local_patch, Parameter_coord_at_GP, info);
				tilde_coord(glo_coord_tilde, Parameter_coord_at_GP, info->Global_local_patch, ele_glo, info);
				physical_coord(ele_glo, glo_coord_tilde, out_coord, info);
				for (int l = 0; l < info->DIMENSION; l++)
					info->PhysicalCoordinate_at_GP[gp_base * info->DIMENSION + l] = out_coord[l];
			
				// displacement
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_GP[gp_base * info->DIMENSION + l] += R_disp[k] * d;
						u[k * info->DIMENSION + l] = d;
					}
				
				// strain
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_GP[gp_base * N_STRESS + k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
			}
			// overlay global
			int status_glo_overlay = 0;
			if (i < info->Total_Element_on_mesh[0] && Total_mesh > 1 && info->Element_patch[element] == info->Global_local_patch)
			{
				double temp_point_patch[MAX_DIMENSION] = {0.0};
				trans_ele_patch_coord(temp_point_patch, temp_point, info->Element_patch[element], element, info);

				// make temp_point_glo
				int itr_n = 0, loc_patch = 0;
				for (int k = info->Total_Patch_to_mesh[1]; k < info->Total_Patch_to_mesh[Total_mesh]; k++)
				{
					itr_n = calc_local_patch_parameter_coord(temp_point_patch, k, temp_para_loc, info);
					loc_patch = k;
					if (itr_n != ERROR)
					{
						status_glo_overlay = 1;
						break;
					}
				}

				int element_loc = 0;
				if (status_glo_overlay)
				{
					element_loc = ele_check(loc_patch, temp_para_loc, info);
					tilde_coord(temp_point_loc, temp_para_loc, loc_patch, element_loc, info);
				}

				// overlay displacement
				if (status_glo_overlay)
				{
					vector<double> R_loc(MAX_NO_CP_ON_ELEMENT);
					vector<double> bl_loc(D_MATRIX_SIZE * MAX_KIEL_SIZE);
					vector<double> u_loc(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
					Bspline_shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
					Make_B_Linear(element_loc, temp_point_loc, bl_loc.data(), info);
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							double d_loc = info->Displacement[info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							info->Displacement_at_GP[gp_base * info->DIMENSION + l] += R_loc[k] * d_loc;
							u_loc[k * info->DIMENSION + l] = d_loc;
						}

					// strain
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
								info->Strain_at_GP[gp_base * N_STRESS + k] += bl_loc[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u_loc[l * info->DIMENSION + m];
				}
			}

			// overlay local
			if (i >= info->Total_Element_on_mesh[0] && Total_mesh > 1)
			{
				geo_parameter_coord(geo_element, geo_coord_tilde, temp_para_glo, info);
				int element_glo = ele_check(info->Global_local_patch, temp_para_glo, info);
				tilde_coord(temp_point_glo, temp_para_glo, info->Global_local_patch, element_glo, info);

				// overlay displacement
				vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl_glo(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u_glo(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
				shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
				Make_B_Linear(element_glo, temp_point_glo, bl_glo.data(), info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_GP[gp_base * info->DIMENSION + l] += R_glo[k] * d_glo;
						u_glo[k * info->DIMENSION + l] = d_glo;
					}
				
				// strain
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->Strain_at_GP[gp_base * N_STRESS + k] += bl_glo[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u_glo[l * info->DIMENSION + m];
			}

			// stress
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				for (int l = 0; l < D_MATRIX_SIZE; l++)
					info->Stress_at_GP[gp_base * N_STRESS + k] += info->D[k * D_MATRIX_SIZE + l] * info->Strain_at_GP[gp_base * N_STRESS + l];

			// 2D specific adjustments
			if (info->DIMENSION == 2)
			{
				// plane stress condition
				if (info->c.ANALYSIS_MODE == 0)
				{
					info->Stress_at_GP[gp_base * N_STRESS + 3] = -1.0 * nu / E * (info->Stress_at_GP[gp_base * N_STRESS + 0] + info->Stress_at_GP[gp_base * N_STRESS + 1]);
					info->Stress_at_GP[gp_base * N_STRESS + 3] = 0.0;
				}

				// plane strain condition
				else if (info->c.ANALYSIS_MODE == 1)
				{
					info->Stress_at_GP[gp_base * N_STRESS + 3] = 0.0;
					info->Stress_at_GP[gp_base * N_STRESS + 3] = nu * (info->Stress_at_GP[gp_base * N_STRESS + 0] + info->Stress_at_GP[gp_base * N_STRESS + 1]);
				}
			}

			// output file
			// displacement
			fprintf(fp_disp, "%d\t%d", element, point);
			for (int k = 0; k < info->DIMENSION; k++)
			fprintf(fp_disp, "\t%.15e", info->Displacement_at_GP[gp_base * info->DIMENSION + k]);
			fprintf(fp_disp, "\n");
			// printf("wrote displacement\n");
			
			// strain
			fprintf(fp_strain, "%d\t%d", element, point);
			for (int k = 0; k < N_STRAIN; k++)
			fprintf(fp_strain, "\t%.15e", info->Strain_at_GP[gp_base * N_STRAIN + k]);
			fprintf(fp_strain, "\n");
			// printf("wrote strain\n");
			
			// stress
			fprintf	(fp_stress, "%d\t%d", element, point);
			for (int k = 0; k < N_STRESS; k++)
				fprintf(fp_stress, "\t%.15e", info->Stress_at_GP[gp_base * N_STRESS + k]);
			fprintf(fp_stress, "\n");
			// printf("wrote stress\n");

			// physical coordinate
			fprintf(fp_pc, "%d\t%d", element, point);
			for (int k = 0; k < info->DIMENSION; k++)
				fprintf(fp_pc, "\t%.15e", info->PhysicalCoordinate_at_GP[gp_base * info->DIMENSION + k]);
			fprintf(fp_pc, "\n");
			// printf("wrote physical coordinate\n");
		}
	}

	if (fp_disp) 
	    fclose(fp_disp);

	if (fp_strain)
	    fclose(fp_strain);

	if (fp_stress)
	    fclose(fp_stress);

	if (fp_pc)
	    fclose(fp_pc);
}


// グローバル解析モデルの自然座標を出力する関数
void output_global_parameters(information *info)
{
	Init_viewer_info(info);
	Make_info_for_viewer_by_shape_func_at_global_parameter_space(info);
	Make_boundary_line(info);

	static int point_on_element = 0;
	static int point_on_edge = 3;
	static char geometry_type[256];
	static char element_type[256];

	char str_glo[256] = "global_patch_at_global_parameter_space.xmf2";
	static FILE *fp_glo = fopen(str_glo, "w");
	char str_loc[256] = "local_patch_at_global_parameter_space.xmf2";
	char str_loc_bl[256] = "local_patch_boundary_line_at_global_parameter_space.xmf2";
	char str_loc_cp[256] = "local_patch_control_point_at_global_parameter_space.xmf2";
	static FILE *fp_loc = fopen(str_loc, "w");
	static FILE *fp_loc_bl = fopen(str_loc_bl, "w");
	static FILE *fp_loc_cp = fopen(str_loc_cp, "w");

	
	// set point_on_element, geometry_type, element_type
	static int count = 0;
	if (count++ == 0)
	{
		if (info->DIMENSION == 2)
		{
			point_on_element = 9;
			strcpy(geometry_type, "XY");
        	strcpy(element_type, "Quadrilateral_9");
		}
		else if (info->DIMENSION == 3)
		{
			point_on_element = 27;
			strcpy(geometry_type, "XYZ");
			strcpy(element_type, "HEXAHEDRON_27");
		}
	}

	// set step
	static double time = 0.0;
	static int step = 0;

	// first step
	fprintf(fp_glo, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp_glo, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp_glo, "  <Domain>\n");
	fprintf(fp_glo, "    <Grid Name=\"global_patch_at_global_parameter_space.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

	fprintf(fp_loc, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp_loc, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp_loc, "  <Domain>\n");
	fprintf(fp_loc, "    <Grid Name=\"local_patch_at_global_parameter_space.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
	
	fprintf(fp_loc_bl, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp_loc_bl, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp_loc_bl, "  <Domain>\n");
	fprintf(fp_loc_bl, "    <Grid Name=\"local_patch_boundary_line_at_global_parameter_space.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
	
	fprintf(fp_loc_cp, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp_loc_cp, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp_loc_cp, "  <Domain>\n");
	fprintf(fp_loc_cp, "    <Grid Name=\"local_patch_control_point_at_global_parameter_space.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

	// middle step
	#pragma region middle step
	{
		// HDF5 dataset creation helper: apply chunking (~512KB) and gzip compression level 6
		auto make_plist = [&](const hsize_t* dims, int rank, const DataType& dtype) {
			DSetCreatPropList plist;
			size_t type_size = dtype.getSize();
			const size_t target_bytes = 512u * 1024u; // ~512KB per chunk
			if (rank == 1)
			{
				hsize_t chunk[1];
				size_t elems = type_size ? (target_bytes / type_size) : target_bytes;
				if (elems < 1) elems = 1;
				chunk[0] = dims[0] < (hsize_t)elems ? dims[0] : (hsize_t)elems;
				plist.setChunk(1, chunk);
			}
			else if (rank == 2)
			{
				hsize_t chunk[2];
				chunk[1] = dims[1]; // take full minor dimension
				size_t denom = type_size * (chunk[1] ? (size_t)chunk[1] : (size_t)1);
				size_t rows = denom ? (target_bytes / denom) : target_bytes;
				if (rows < 1) rows = 1;
				chunk[0] = dims[0] < (hsize_t)rows ? dims[0] : (hsize_t)rows;
				plist.setChunk(2, chunk);
			}
			// enable gzip compression
			// plist.setDeflate(0); // set to compression level 0 (no compression)
			plist.setDeflate(1); // set to compression level 1 (fastest)
			// plist.setDeflate(6);
			// plist.setDeflate(9); // set to maximum compression level 9
			return plist;
		};

		#pragma region middle step global patch
		{
			fprintf(fp_glo, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
			fprintf(fp_glo, "        <Time Value=\"%.20e\"/>\n", time);
			fprintf(fp_glo, "        <Grid Name=\"ien\">\n");
			fprintf(fp_glo, "          <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n", element_type, info->Total_Element_to_mesh[1]);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Topology>\n");
			fprintf(fp_glo, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Geometry>\n");
			fprintf(fp_glo, "          <Attribute Name=\"displacement\" AttributeType=\"Vector\" Dimensions=\"%d %d\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;displacement&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"stress\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/xyz\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_TRUNC);
				hsize_t dims_xyz[2] = { (hsize_t)Total_connectivity_glo, (hsize_t)info->DIMENSION };
				DataSpace dataspace_xyz(2, dims_xyz);
				DSetCreatPropList plist_xyz = make_plist(dims_xyz, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_xyz, plist_xyz);
				vector<double> buf_xyz((size_t)Total_connectivity_glo * info->DIMENSION);
				for (int i = 0; i < Total_connectivity_glo; i++)
					for (int j = 0; j < info->DIMENSION; j++)
						buf_xyz[i * info->DIMENSION + j] = info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j];
				dset_xyz.write(buf_xyz.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Element_to_mesh[1], point_on_element);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/ien\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_ien[2] = { (hsize_t)info->Total_Element_to_mesh[1], (hsize_t)point_on_element };
				DataSpace dataspace_ien(2, dims_ien);
				DSetCreatPropList plist_ien = make_plist(dims_ien, 2, PredType::NATIVE_INT);
				DataSet dset_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_ien, plist_ien);
				vector<int> buf_ien((size_t)info->Total_Element_to_mesh[1] * point_on_element);
				for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
					for (int j = 0; j < point_on_element; j++)
						buf_ien[i * point_on_element + j] = info->Connectivity[i * point_on_element + j];
				dset_ien.write(buf_ien.data(), PredType::NATIVE_INT);
			}
			fprintf(fp_glo, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/displacement\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_disp[2] = { (hsize_t)Total_connectivity_glo, (hsize_t)info->DIMENSION };
				DataSpace dataspace_disp(2, dims_disp);
				DSetCreatPropList plist_disp = make_plist(dims_disp, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_disp = file.createDataSet("/displacement", PredType::NATIVE_DOUBLE, dataspace_disp, plist_disp);
				vector<double> buf_disp((size_t)Total_connectivity_glo * info->DIMENSION);
				for (int i = 0; i < Total_connectivity_glo; i++)
					for (int j = 0; j < info->DIMENSION; j++)
						buf_disp[i * info->DIMENSION + j] = info->disp_at_connectivity[i * info->DIMENSION + j];
				dset_disp.write(buf_disp.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"stress\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/stress\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_stress[2] = { (hsize_t)Total_connectivity_glo, 6 };
				DataSpace dataspace_stress(2, dims_stress);
				DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
				vector<double> buf_stress((size_t)Total_connectivity_glo * 6);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 2];
						buf_stress[i * 6 + 2] = 0.0;
						buf_stress[i * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 6 + 4] = 0.0;
						buf_stress[i * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 3];
						buf_stress[i * 6 + 2] = info->stress_at_connectivity[i * N_STRESS + 4];
						buf_stress[i * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 6 + 4] = info->stress_at_connectivity[i * N_STRESS + 5];
						buf_stress[i * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 2];
					}
				}
				dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/elastic_strain\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_es[2] = { (hsize_t)Total_connectivity_glo, 6 };
				DataSpace dataspace_es(2, dims_es);
				DSetCreatPropList plist_es = make_plist(dims_es, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_es = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_es, plist_es);
				vector<double> buf_es((size_t)Total_connectivity_glo * 6);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 2];
						buf_es[i * 6 + 2] = 0.0;
						buf_es[i * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 6 + 4] = 0.0;
						buf_es[i * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 3];
						buf_es[i * 6 + 2] = info->strain_at_connectivity[i * N_STRAIN + 4];
						buf_es[i * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 6 + 4] = info->strain_at_connectivity[i * N_STRAIN + 5];
						buf_es[i * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 2];
					}
				}
				dset_es.write(buf_es.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/eqs\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_eqs[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_eqs(1, dims_eqs);
				DSetCreatPropList plist_eqs = make_plist(dims_eqs, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_eqs = file.createDataSet("/eqs", PredType::NATIVE_DOUBLE, dataspace_eqs, plist_eqs);
				vector<double> buf_eqs((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_eqs[i] = info->vm_at_connectivity[i];
				dset_eqs.write(buf_eqs.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/hs\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_hs[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_hs(1, dims_hs);
				DSetCreatPropList plist_hs = make_plist(dims_hs, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_hs = file.createDataSet("/hs", PredType::NATIVE_DOUBLE, dataspace_hs, plist_hs);
				vector<double> buf_hs((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_hs[i] = info->hs_at_connectivity[i];
				dset_hs.write(buf_hs.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "        </Grid>\n");
			fprintf(fp_glo, "      </Grid>\n");
		}
		#pragma endregion

		// local patch
		if (Total_mesh > 1)
		{
			#pragma region middle step local patch
			{
				fprintf(fp_loc, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc, "          <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n", element_type, info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Topology>\n");
				fprintf(fp_loc, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Geometry>\n");
				fprintf(fp_loc, "          <Attribute Name=\"displacement\" AttributeType=\"Vector\" Dimensions=\"%d %d\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;displacement&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"stress\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_xyz[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), (hsize_t)info->DIMENSION };
					DataSpace dataspace_xyz(2, dims_xyz);
					DSetCreatPropList plist_xyz = make_plist(dims_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_xyz, plist_xyz);
					vector<double> buf_xyz((size_t)(Total_connectivity - Total_connectivity_glo) * info->DIMENSION);
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_xyz[(i - Total_connectivity_glo) * info->DIMENSION + j] = info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j];
					dset_xyz.write(buf_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1], point_on_element);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_ien[2] = { (hsize_t)(info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]), (hsize_t)point_on_element };
					DataSpace dataspace_ien(2, dims_ien);
					DSetCreatPropList plist_ien = make_plist(dims_ien, 2, PredType::NATIVE_INT);
					DataSet dset_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_ien, plist_ien);
					vector<int> buf_ien((size_t)(info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]) * point_on_element);
					for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
						for (int j = 0; j < point_on_element; j++)
							buf_ien[(i - info->Total_Element_to_mesh[1]) * point_on_element + j] = info->Connectivity[i * point_on_element + j] - Total_connectivity_glo;
					dset_ien.write(buf_ien.data(), PredType::NATIVE_INT);
				}
				fprintf(fp_loc, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/displacement\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_disp[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), (hsize_t)info->DIMENSION };
					DataSpace dataspace_disp(2, dims_disp);
					DSetCreatPropList plist_disp = make_plist(dims_disp, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_disp = file.createDataSet("/displacement", PredType::NATIVE_DOUBLE, dataspace_disp, plist_disp);
					vector<double> buf_disp((size_t)(Total_connectivity - Total_connectivity_glo) * info->DIMENSION);
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_disp[(i - Total_connectivity_glo) * info->DIMENSION + j] = info->disp_at_connectivity[i * info->DIMENSION + j];
					dset_disp.write(buf_disp.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"stress\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_stress[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 6 };
					DataSpace dataspace_stress(2, dims_stress);
					DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
					vector<double> buf_stress((size_t)(Total_connectivity - Total_connectivity_glo) * 6);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 2];
							buf_stress[(i - Total_connectivity_glo) * 6 + 2] = 0.0;
							buf_stress[(i - Total_connectivity_glo) * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 6 + 4] = 0.0;
							buf_stress[(i - Total_connectivity_glo) * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 3];
							buf_stress[(i - Total_connectivity_glo) * 6 + 2] = info->stress_at_connectivity[i * N_STRESS + 4];
							buf_stress[(i - Total_connectivity_glo) * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 6 + 4] = info->stress_at_connectivity[i * N_STRESS + 5];
							buf_stress[(i - Total_connectivity_glo) * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 2];
						}
					}
					dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/elastic_strain\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_elastic_strain[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 6 };
					DataSpace dataspace_elastic_strain(2, dims_elastic_strain);
					DSetCreatPropList plist_elastic = make_plist(dims_elastic_strain, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_elastic_strain = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_elastic_strain, plist_elastic);
					vector<double> buf_elastic_strain((size_t)(Total_connectivity - Total_connectivity_glo) * 6);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 2];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 2] = 0.0;
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 4] = 0.0;
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 3];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 2] = info->strain_at_connectivity[i * N_STRAIN + 4];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 4] = info->strain_at_connectivity[i * N_STRAIN + 5];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 2];
						}
					}
					dset_elastic_strain.write(buf_elastic_strain.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/equivalent_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_eqs[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_eqs(1, dims_eqs);
					DSetCreatPropList plist_eqs = make_plist(dims_eqs, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_eqs = file.createDataSet("/equivalent_stress", PredType::NATIVE_DOUBLE, dataspace_eqs, plist_eqs);
					vector<double> buf_eqs((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_eqs[i - Total_connectivity_glo] = info->vm_at_connectivity[i];
					dset_eqs.write(buf_eqs.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/hydrostatic_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_hs[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_hs(1, dims_hs);
					DSetCreatPropList plist_hs = make_plist(dims_hs, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_hs = file.createDataSet("/hydrostatic_stress", PredType::NATIVE_DOUBLE, dataspace_hs, plist_hs);
					vector<double> buf_hs((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_hs[i - Total_connectivity_glo] = info->hs_at_connectivity[i];
					dset_hs.write(buf_hs.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "        </Grid>\n");
				fprintf(fp_loc, "      </Grid>\n");
			}
			#pragma endregion

			#pragma region middle step local control point
			{
				fprintf(fp_loc_cp, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc_cp, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc_cp, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc_cp, "          <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\">\n", info->Geo_Total_Control_Point_on_mesh);
				fprintf(fp_loc_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_cp, "          </Topology>\n");
				fprintf(fp_loc_cp, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_cp, "          </Geometry>\n");
				fprintf(fp_loc_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", info->Geo_Total_Control_Point_on_mesh, info->DIMENSION);
				fprintf(fp_loc_cp, "            %s/local_patch_cp_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_cp_step") + to_string(step) + ".h5";
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_cp_xyz[2] = { (hsize_t)(info->Geo_Total_Control_Point_on_mesh), (hsize_t)(info->DIMENSION) };
					DataSpace dataspace_cp_xyz(2, dims_cp_xyz);
					DSetCreatPropList plist_cp_xyz = make_plist(dims_cp_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_cp_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_cp_xyz, plist_cp_xyz);
					vector<double> buf_cp_xyz((size_t)(info->Geo_Total_Control_Point_on_mesh) * info->DIMENSION);
					if (Total_mesh < 2)
						for (int i = 0; i < info->Geo_Total_Control_Point_on_mesh; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[i * info->DIMENSION + j] = info->Geo_Node_Coordinate[i * (info->DIMENSION + 1) + j]; // deformation of control points is not considered
					else
						for (int i = 0; i < info->Geo_Total_Control_Point_on_mesh; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[i * info->DIMENSION + j] = info->Geo_Node_Coordinate[i * (info->DIMENSION + 1) + j];
					dset_cp_xyz.write(buf_cp_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Geo_Total_Control_Point_on_mesh);
				fprintf(fp_loc_cp, "            %s/local_patch_cp_gps_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_cp_gps_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_cp_ien[2] = { (hsize_t)(info->Geo_Total_Control_Point_on_mesh), 1 };
					DataSpace dataspace_cp_ien(2, dims_cp_ien);
					DSetCreatPropList plist_cp_ien = make_plist(dims_cp_ien, 2, PredType::NATIVE_INT);
					DataSet dset_cp_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_cp_ien, plist_cp_ien);
					vector<int> buf_cp_ien((size_t)(info->Geo_Total_Control_Point_on_mesh));
					for (int i = 0; i < info->Geo_Total_Control_Point_on_mesh; i++)
						buf_cp_ien[i] = i;
					dset_cp_ien.write(buf_cp_ien.data(), PredType::NATIVE_INT);
				}

				fprintf(fp_loc_cp, "        </Grid>\n");
				fprintf(fp_loc_cp, "      </Grid>\n");
			}

			#pragma region middle step local patch boundary line
			{
				fprintf(fp_loc_bl, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc_bl, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc_bl, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc_bl, "          <Topology TopologyType=\"EDGE_3\" NumberOfElements=\"%d\">\n", Total_edge - Total_edge_glo);
				fprintf(fp_loc_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_bl, "          </Topology>\n");
				fprintf(fp_loc_bl, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_bl, "          </Geometry>\n");
				fprintf(fp_loc_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", (Total_edge - Total_edge_glo) * point_on_edge, info->DIMENSION);
				fprintf(fp_loc_bl, "            %s/local_boundary_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_boundary_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_boundary_xyz[2] = { (hsize_t)((Total_edge - Total_edge_glo) * point_on_edge), (hsize_t)(info->DIMENSION) };
					DataSpace dataspace_boundary_xyz(2, dims_boundary_xyz);
					DSetCreatPropList plist_bxyz = make_plist(dims_boundary_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_boundary_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_boundary_xyz, plist_bxyz);
					vector<double> buf_boundary_xyz((size_t)(Total_edge - Total_edge_glo) * point_on_edge * info->DIMENSION);
					for (int i = Total_edge_glo; i < Total_edge; i++)
						for (int j = 0; j < point_on_edge; j++)
							for (int k = 0; k < info->DIMENSION; k++)
								buf_boundary_xyz[(i - Total_edge_glo) * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k];
					dset_boundary_xyz.write(buf_boundary_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", Total_edge - Total_edge_glo, point_on_edge);
				fprintf(fp_loc_bl, "            %s/local_boundary_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_boundary_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_boundary_ien[2] = { (hsize_t)(Total_edge - Total_edge_glo), (hsize_t)(point_on_edge) };
					DataSpace dataspace_boundary_ien(2, dims_boundary_ien);
					DSetCreatPropList plist_bien = make_plist(dims_boundary_ien, 2, PredType::NATIVE_INT);
					DataSet dset_boundary_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_boundary_ien, plist_bien);
					vector<int> buf_boundary_ien((size_t)(Total_edge - Total_edge_glo) * point_on_edge);
					for (int i = 0; i < Total_edge - Total_edge_glo; i++)
						for (int j = 0; j < point_on_edge; j++)
							buf_boundary_ien[i * point_on_edge + j] = i * point_on_edge + j;
					dset_boundary_ien.write(buf_boundary_ien.data(), PredType::NATIVE_INT);
				}

				fprintf(fp_loc_bl, "        </Grid>\n");
				fprintf(fp_loc_bl, "      </Grid>\n");
			}
			#pragma endregion
		}
	}
	#pragma endregion

	fprintf(fp_glo, "    </Grid>\n");
	fprintf(fp_glo, "  </Domain>\n");
	fprintf(fp_glo, "</Xdmf>\n");

	fprintf(fp_loc, "    </Grid>\n");
	fprintf(fp_loc, "  </Domain>\n");
	fprintf(fp_loc, "</Xdmf>\n");
	
	fprintf(fp_loc_bl, "    </Grid>\n");
	fprintf(fp_loc_bl, "  </Domain>\n");
	fprintf(fp_loc_bl, "</Xdmf>\n");

	fprintf(fp_loc_cp, "    </Grid>\n");
	fprintf(fp_loc_cp, "  </Domain>\n");
	fprintf(fp_loc_cp, "</Xdmf>\n");

	step++;
}


// paraview
void output_for_paraview_timestep(information *info, bool isGeometryOnly, double factor)
{
	Init_viewer_info(info);
	Make_info_for_viewer_by_shape_func(info);
	Make_boundary_line(info);

	static int point_on_element = 0;
	static int point_on_edge = 3;
	static char geometry_type[256];
	static char element_type[256];

	char str_glo[256] = "global_patch.xmf2";
	char str_glo_cp[256] = "global_patch_control_point.xmf2";
	char str_glo_bl[256] = "global_patch_boundary_line.xmf2";
	static FILE *fp_glo = fopen(str_glo, "w");
	static FILE *fp_glo_cp = fopen(str_glo_cp, "w");
	static FILE *fp_glo_bl = fopen(str_glo_bl, "w");
	char str_loc[256] = "local_patch.xmf2";
	char str_loc_bl[256] = "local_patch_boundary_line.xmf2";
	static FILE *fp_loc = fopen(str_loc, "w");
	static FILE *fp_loc_bl = fopen(str_loc_bl, "w");

	// ローカル制御点は、グローバルパラメータ空間に写像するために使用され、物理座標で意味を持たない
	// char str_loc_cp[256] = "local_patch_control_point.xmf2";
	// static FILE *fp_loc_cp = fopen(str_loc_cp, "w");
	
	// set point_on_element, geometry_type, element_type
	static int count = 0;
	if (count++ == 0)
	{
		if (info->DIMENSION == 2)
		{
			point_on_element = 9;
			strcpy(geometry_type, "XY");
        	strcpy(element_type, "Quadrilateral_9");
		}
		else if (info->DIMENSION == 3)
		{
			point_on_element = 27;
			strcpy(geometry_type, "XYZ");
			strcpy(element_type, "HEXAHEDRON_27");
		}
	}

	// for debug
	// for (int i = 0; i < Total_connectivity; i++)
	// {
	// 	cout << "displacement at " << i << ": " 
	// 		 << info->disp_at_connectivity[i * info->DIMENSION + 0] << ", "
	// 		 << info->disp_at_connectivity[i * info->DIMENSION + 1] << endl;
	// }

	// set time
	static const double total_time = 10.0;
	double time = total_time * factor;
	if (fabs(time) < MERGE_ERROR)
		time = 0.0;
	else if (fabs(time - total_time) < MERGE_ERROR)
		time = total_time;

	// set step
	static int step = 0;

	// first step
	if (time == 0.0)
	{
		fprintf(fp_glo, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_glo, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_glo, "  <Domain>\n");
		fprintf(fp_glo, "    <Grid Name=\"global_patch.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

		fprintf(fp_glo_cp, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_glo_cp, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_glo_cp, "  <Domain>\n");
		fprintf(fp_glo_cp, "    <Grid Name=\"global_patch_control_point.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

		fprintf(fp_glo_bl, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_glo_bl, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_glo_bl, "  <Domain>\n");
		fprintf(fp_glo_bl, "    <Grid Name=\"global_patch_boundary_line.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

		fprintf(fp_loc, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_loc, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_loc, "  <Domain>\n");
		fprintf(fp_loc, "    <Grid Name=\"local_patch.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
		
		fprintf(fp_loc_bl, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_loc_bl, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_loc_bl, "  <Domain>\n");
		fprintf(fp_loc_bl, "    <Grid Name=\"local_patch_boundary_line.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
		
		// ローカル制御点はグローバルパラメータ空間内に配置されるため、物理座標へ可視化する必要がない。
		// fprintf(fp_loc_cp, "<?xml version=\"1.0\" ?>\n");
		// fprintf(fp_loc_cp, "<Xdmf Version=\"2.0\">\n");
		// fprintf(fp_loc_cp, "  <Domain>\n");
		// fprintf(fp_loc_cp, "    <Grid Name=\"local_patch_control_point.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
	}

	// middle step
	#pragma region middle step
	{
		// HDF5 dataset creation helper: apply chunking (~512KB) and gzip compression level 6
		auto make_plist = [&](const hsize_t* dims, int rank, const DataType& dtype) {
			DSetCreatPropList plist;
			size_t type_size = dtype.getSize();
			const size_t target_bytes = 512u * 1024u; // ~512KB per chunk
			if (rank == 1)
			{
				hsize_t chunk[1];
				size_t elems = type_size ? (target_bytes / type_size) : target_bytes;
				if (elems < 1) elems = 1;
				chunk[0] = dims[0] < (hsize_t)elems ? dims[0] : (hsize_t)elems;
				plist.setChunk(1, chunk);
			}
			else if (rank == 2)
			{
				hsize_t chunk[2];
				chunk[1] = dims[1]; // take full minor dimension
				size_t denom = type_size * (chunk[1] ? (size_t)chunk[1] : (size_t)1);
				size_t rows = denom ? (target_bytes / denom) : target_bytes;
				if (rows < 1) rows = 1;
				chunk[0] = dims[0] < (hsize_t)rows ? dims[0] : (hsize_t)rows;
				plist.setChunk(2, chunk);
			}
			// enable gzip compression
			// plist.setDeflate(0); // set to compression level 0 (no compression)
			plist.setDeflate(1); // set to compression level 1 (fastest)
			// plist.setDeflate(6);
			// plist.setDeflate(9); // set to maximum compression level 9
			return plist;
		};

		#pragma region middle step global patch
		{
			fprintf(fp_glo, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
			fprintf(fp_glo, "        <Time Value=\"%.20e\"/>\n", time);
			fprintf(fp_glo, "        <Grid Name=\"ien\">\n");
			fprintf(fp_glo, "          <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n", element_type, info->Total_Element_to_mesh[1]);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Topology>\n");
			fprintf(fp_glo, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Geometry>\n");
			fprintf(fp_glo, "          <Attribute Name=\"displacement\" AttributeType=\"Vector\" Dimensions=\"%d %d\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;displacement&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"stress\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/xyz\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_TRUNC);
				hsize_t dims_xyz[2] = { (hsize_t)Total_connectivity_glo, (hsize_t)info->DIMENSION };
				DataSpace dataspace_xyz(2, dims_xyz);
				DSetCreatPropList plist_xyz = make_plist(dims_xyz, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_xyz, plist_xyz);
				vector<double> buf_xyz((size_t)Total_connectivity_glo * info->DIMENSION);
				for (int i = 0; i < Total_connectivity_glo; i++)
					for (int j = 0; j < info->DIMENSION; j++)
						buf_xyz[i * info->DIMENSION + j] = info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j];
				dset_xyz.write(buf_xyz.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Element_to_mesh[1], point_on_element);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/ien\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_ien[2] = { (hsize_t)info->Total_Element_to_mesh[1], (hsize_t)point_on_element };
				DataSpace dataspace_ien(2, dims_ien);
				DSetCreatPropList plist_ien = make_plist(dims_ien, 2, PredType::NATIVE_INT);
				DataSet dset_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_ien, plist_ien);
				vector<int> buf_ien((size_t)info->Total_Element_to_mesh[1] * point_on_element);
				for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
					for (int j = 0; j < point_on_element; j++)
						buf_ien[i * point_on_element + j] = info->Connectivity[i * point_on_element + j];
				dset_ien.write(buf_ien.data(), PredType::NATIVE_INT);
			}
			fprintf(fp_glo, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/displacement\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_disp[2] = { (hsize_t)Total_connectivity_glo, (hsize_t)info->DIMENSION };
				DataSpace dataspace_disp(2, dims_disp);
				DSetCreatPropList plist_disp = make_plist(dims_disp, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_disp = file.createDataSet("/displacement", PredType::NATIVE_DOUBLE, dataspace_disp, plist_disp);
				vector<double> buf_disp((size_t)Total_connectivity_glo * info->DIMENSION);
				for (int i = 0; i < Total_connectivity_glo; i++)
					for (int j = 0; j < info->DIMENSION; j++)
						buf_disp[i * info->DIMENSION + j] = info->disp_at_connectivity[i * info->DIMENSION + j];
				dset_disp.write(buf_disp.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"stress\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/stress\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_stress[2] = { (hsize_t)Total_connectivity_glo, 6 };
				DataSpace dataspace_stress(2, dims_stress);
				DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
				vector<double> buf_stress((size_t)Total_connectivity_glo * 6);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 2];
						buf_stress[i * 6 + 2] = 0.0;
						buf_stress[i * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 6 + 4] = 0.0;
						buf_stress[i * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 3];
						buf_stress[i * 6 + 2] = info->stress_at_connectivity[i * N_STRESS + 4];
						buf_stress[i * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 6 + 4] = info->stress_at_connectivity[i * N_STRESS + 5];
						buf_stress[i * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 2];
					}
				}
				dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/elastic_strain\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_es[2] = { (hsize_t)Total_connectivity_glo, 6 };
				DataSpace dataspace_es(2, dims_es);
				DSetCreatPropList plist_es = make_plist(dims_es, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_es = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_es, plist_es);
				vector<double> buf_es((size_t)Total_connectivity_glo * 6);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 2];
						buf_es[i * 6 + 2] = 0.0;
						buf_es[i * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 6 + 4] = 0.0;
						buf_es[i * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 3];
						buf_es[i * 6 + 2] = info->strain_at_connectivity[i * N_STRAIN + 4];
						buf_es[i * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 6 + 4] = info->strain_at_connectivity[i * N_STRAIN + 5];
						buf_es[i * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 2];
					}
				}
				dset_es.write(buf_es.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/eqs\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_eqs[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_eqs(1, dims_eqs);
				DSetCreatPropList plist_eqs = make_plist(dims_eqs, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_eqs = file.createDataSet("/eqs", PredType::NATIVE_DOUBLE, dataspace_eqs, plist_eqs);
				vector<double> buf_eqs((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_eqs[i] = info->vm_at_connectivity[i];
				dset_eqs.write(buf_eqs.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/hs\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_hs[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_hs(1, dims_hs);
				DSetCreatPropList plist_hs = make_plist(dims_hs, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_hs = file.createDataSet("/hs", PredType::NATIVE_DOUBLE, dataspace_hs, plist_hs);
				vector<double> buf_hs((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_hs[i] = info->hs_at_connectivity[i];
				dset_hs.write(buf_hs.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "        </Grid>\n");
			fprintf(fp_glo, "      </Grid>\n");
		}
		#pragma endregion

		#pragma region middle step global control point
		{
			fprintf(fp_glo_cp, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
			fprintf(fp_glo_cp, "        <Time Value=\"%.20e\"/>\n", time);
			fprintf(fp_glo_cp, "        <Grid Name=\"ien\">\n");
			fprintf(fp_glo_cp, "          <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\">\n", info->Total_Control_Point_to_mesh[1]);
			fprintf(fp_glo_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo_cp, "          </Topology>\n");
			fprintf(fp_glo_cp, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
			fprintf(fp_glo_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo_cp, "          </Geometry>\n");
			fprintf(fp_glo_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", info->Total_Control_Point_to_mesh[1], info->DIMENSION);
			fprintf(fp_glo_cp, "            %s/global_cp_step%d.h5:/cp_xyz\n", "./bin", step);
			fprintf(fp_glo_cp, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_cp_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_TRUNC);
				hsize_t dims_cp_xyz[2] = { (hsize_t)info->Total_Control_Point_to_mesh[1], (hsize_t)info->DIMENSION };
				DataSpace dataspace_cp_xyz(2, dims_cp_xyz);
				DSetCreatPropList plist_cp_xyz = make_plist(dims_cp_xyz, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_cp_xyz = file.createDataSet("/cp_xyz", PredType::NATIVE_DOUBLE, dataspace_cp_xyz, plist_cp_xyz);
				vector<double> buf_cp_xyz((size_t)info->Total_Control_Point_to_mesh[1] * info->DIMENSION);
				if (Total_mesh < 2)
					for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_cp_xyz[i * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->Displacement[i * info->DIMENSION + j]);
				else
					for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_cp_xyz[i * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j];
				dset_cp_xyz.write(buf_cp_xyz.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Control_Point_to_mesh[1]);
			fprintf(fp_glo_cp, "            %s/global_cp_step%d.h5:/cp_ien\n", "./bin", step);
			fprintf(fp_glo_cp, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_cp_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_cp_ien[2] = { (hsize_t)info->Total_Control_Point_to_mesh[1], 1 };
				DataSpace dataspace_cp_ien(2, dims_cp_ien);
				DSetCreatPropList plist_cp_ien = make_plist(dims_cp_ien, 2, PredType::NATIVE_INT);
				DataSet dset_cp_ien = file.createDataSet("/cp_ien", PredType::NATIVE_INT, dataspace_cp_ien, plist_cp_ien);
				vector<int> buf_cp_ien((size_t)info->Total_Control_Point_to_mesh[1]);
				for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
					buf_cp_ien[i] = i;
				dset_cp_ien.write(buf_cp_ien.data(), PredType::NATIVE_INT);
			}
			fprintf(fp_glo_cp, "        </Grid>\n");
			fprintf(fp_glo_cp, "      </Grid>\n");
		}
		#pragma endregion

		#pragma region middle step global patch boundary line
		{
			fprintf(fp_glo_bl, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
			fprintf(fp_glo_bl, "        <Time Value=\"%.20e\"/>\n", time);
			fprintf(fp_glo_bl, "        <Grid Name=\"ien\">\n");
			fprintf(fp_glo_bl, "          <Topology TopologyType=\"EDGE_3\" NumberOfElements=\"%d\">\n", Total_edge_glo);
			fprintf(fp_glo_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo_bl, "          </Topology>\n");
			fprintf(fp_glo_bl, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
			fprintf(fp_glo_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo_bl, "          </Geometry>\n");
			fprintf(fp_glo_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_edge_glo * point_on_edge, info->DIMENSION);
			fprintf(fp_glo_bl, "            %s/global_boundary_step%d.h5:/boundary_xyz\n", "./bin", step);
			fprintf(fp_glo_bl, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_boundary_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_TRUNC);
				hsize_t dims_b_xyz[2] = { (hsize_t)(Total_edge_glo * point_on_edge), (hsize_t)info->DIMENSION };
				DataSpace dataspace_b_xyz(2, dims_b_xyz);
				DSetCreatPropList plist_b_xyz = make_plist(dims_b_xyz, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_b_xyz = file.createDataSet("/boundary_xyz", PredType::NATIVE_DOUBLE, dataspace_b_xyz, plist_b_xyz);
				vector<double> buf_b_xyz((size_t)(Total_edge_glo * point_on_edge) * info->DIMENSION);
				for (int i = 0; i < Total_edge_glo; i++)
					for (int j = 0; j < point_on_edge; j++)
						for (int k = 0; k < info->DIMENSION; k++)
							buf_b_xyz[(i * point_on_edge + j) * info->DIMENSION + k] = info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k];
				dset_b_xyz.write(buf_b_xyz.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", Total_edge_glo, point_on_edge);
			fprintf(fp_glo_bl, "            %s/global_boundary_step%d.h5:/boundary_ien\n", "./bin", step);
			fprintf(fp_glo_bl, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_boundary_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_b_ien[2] = { (hsize_t)Total_edge_glo, (hsize_t)point_on_edge };
				DataSpace dataspace_b_ien(2, dims_b_ien);
				DSetCreatPropList plist_b_ien = make_plist(dims_b_ien, 2, PredType::NATIVE_INT);
				DataSet dset_b_ien = file.createDataSet("/boundary_ien", PredType::NATIVE_INT, dataspace_b_ien, plist_b_ien);
				vector<int> buf_b_ien((size_t)Total_edge_glo * point_on_edge);
				for (int i = 0; i < Total_edge_glo; i++)
					for (int j = 0; j < point_on_edge; j++)
						buf_b_ien[i * point_on_edge + j] = i * 3 + j;
				dset_b_ien.write(buf_b_ien.data(), PredType::NATIVE_INT);
			}
			fprintf(fp_glo_bl, "        </Grid>\n");
			fprintf(fp_glo_bl, "      </Grid>\n");
		}
		#pragma endregion

		// local patch
		if (Total_mesh > 1)
		{
			#pragma region middle step local patch
			{
				fprintf(fp_loc, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc, "          <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n", element_type, info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Topology>\n");
				fprintf(fp_loc, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Geometry>\n");
				fprintf(fp_loc, "          <Attribute Name=\"displacement\" AttributeType=\"Vector\" Dimensions=\"%d %d\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;displacement&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"stress\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Tensor6\" Dimensions=\"%d 6\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_xyz[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), (hsize_t)info->DIMENSION };
					DataSpace dataspace_xyz(2, dims_xyz);
					DSetCreatPropList plist_xyz = make_plist(dims_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_xyz, plist_xyz);
					vector<double> buf_xyz((size_t)(Total_connectivity - Total_connectivity_glo) * info->DIMENSION);
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_xyz[(i - Total_connectivity_glo) * info->DIMENSION + j] = info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j];
					dset_xyz.write(buf_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1], point_on_element);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_ien[2] = { (hsize_t)(info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]), (hsize_t)point_on_element };
					DataSpace dataspace_ien(2, dims_ien);
					DSetCreatPropList plist_ien = make_plist(dims_ien, 2, PredType::NATIVE_INT);
					DataSet dset_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_ien, plist_ien);
					vector<int> buf_ien((size_t)(info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]) * point_on_element);
					for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
						for (int j = 0; j < point_on_element; j++)
							buf_ien[(i - info->Total_Element_to_mesh[1]) * point_on_element + j] = info->Connectivity[i * point_on_element + j] - Total_connectivity_glo;
					dset_ien.write(buf_ien.data(), PredType::NATIVE_INT);
				}
				fprintf(fp_loc, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/displacement\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_disp[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), (hsize_t)info->DIMENSION };
					DataSpace dataspace_disp(2, dims_disp);
					DSetCreatPropList plist_disp = make_plist(dims_disp, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_disp = file.createDataSet("/displacement", PredType::NATIVE_DOUBLE, dataspace_disp, plist_disp);
					vector<double> buf_disp((size_t)(Total_connectivity - Total_connectivity_glo) * info->DIMENSION);
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_disp[(i - Total_connectivity_glo) * info->DIMENSION + j] = info->disp_at_connectivity[i * info->DIMENSION + j];
					dset_disp.write(buf_disp.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"stress\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_stress[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 6 };
					DataSpace dataspace_stress(2, dims_stress);
					DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
					vector<double> buf_stress((size_t)(Total_connectivity - Total_connectivity_glo) * 6);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 2];
							buf_stress[(i - Total_connectivity_glo) * 6 + 2] = 0.0;
							buf_stress[(i - Total_connectivity_glo) * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 6 + 4] = 0.0;
							buf_stress[(i - Total_connectivity_glo) * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 6 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 6 + 1] = info->stress_at_connectivity[i * N_STRESS + 3];
							buf_stress[(i - Total_connectivity_glo) * 6 + 2] = info->stress_at_connectivity[i * N_STRESS + 4];
							buf_stress[(i - Total_connectivity_glo) * 6 + 3] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 6 + 4] = info->stress_at_connectivity[i * N_STRESS + 5];
							buf_stress[(i - Total_connectivity_glo) * 6 + 5] = info->stress_at_connectivity[i * N_STRESS + 2];
						}
					}
					dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 6\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/elastic_strain\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_elastic_strain[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 6 };
					DataSpace dataspace_elastic_strain(2, dims_elastic_strain);
					DSetCreatPropList plist_elastic = make_plist(dims_elastic_strain, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_elastic_strain = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_elastic_strain, plist_elastic);
					vector<double> buf_elastic_strain((size_t)(Total_connectivity - Total_connectivity_glo) * 6);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 2];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 2] = 0.0;
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 4] = 0.0;
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 0] = info->strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 1] = info->strain_at_connectivity[i * N_STRAIN + 3];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 2] = info->strain_at_connectivity[i * N_STRAIN + 4];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 3] = info->strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 4] = info->strain_at_connectivity[i * N_STRAIN + 5];
							buf_elastic_strain[(i - Total_connectivity_glo) * 6 + 5] = info->strain_at_connectivity[i * N_STRAIN + 2];
						}
					}
					dset_elastic_strain.write(buf_elastic_strain.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/equivalent_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_eqs[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_eqs(1, dims_eqs);
					DSetCreatPropList plist_eqs = make_plist(dims_eqs, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_eqs = file.createDataSet("/equivalent_stress", PredType::NATIVE_DOUBLE, dataspace_eqs, plist_eqs);
					vector<double> buf_eqs((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_eqs[i - Total_connectivity_glo] = info->vm_at_connectivity[i];
					dset_eqs.write(buf_eqs.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/hydrostatic_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_hs[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_hs(1, dims_hs);
					DSetCreatPropList plist_hs = make_plist(dims_hs, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_hs = file.createDataSet("/hydrostatic_stress", PredType::NATIVE_DOUBLE, dataspace_hs, plist_hs);
					vector<double> buf_hs((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_hs[i - Total_connectivity_glo] = info->hs_at_connectivity[i];
					dset_hs.write(buf_hs.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "        </Grid>\n");
				fprintf(fp_loc, "      </Grid>\n");
			}
			#pragma endregion

			// ローカル制御点はグローバルパラメータ空間内に配置されるため、物理座標へ可視化する必要がない。
			#if 0
			#pragma region middle step local control point
			{
				fprintf(fp_loc_cp, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc_cp, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc_cp, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc_cp, "          <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]);
				fprintf(fp_loc_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_cp, "          </Topology>\n");
				fprintf(fp_loc_cp, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc_cp, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_cp, "          </Geometry>\n");
				fprintf(fp_loc_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1], info->DIMENSION);
				fprintf(fp_loc_cp, "            %s/local_patch_cp_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_cp_step") + to_string(step) + ".h5";
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_cp_xyz[2] = { (hsize_t)(*info->Total_Control_Point_on_mesh), (hsize_t)(info->DIMENSION) };
					DataSpace dataspace_cp_xyz(2, dims_cp_xyz);
					DSetCreatPropList plist_cp_xyz = make_plist(dims_cp_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_cp_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_cp_xyz, plist_cp_xyz);
					vector<double> buf_cp_xyz((size_t)(*info->Total_Control_Point_on_mesh) * info->DIMENSION);
					if (Total_mesh < 2)
						for (int i = 0; i < *info->Total_Control_Point_on_mesh; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[i * info->DIMENSION + j] = info->Geo_Node_Coordinate[i * (info->DIMENSION + 1) + j]; // deformation of control points is not considered
					else
						for (int i = 0; i < *info->Total_Control_Point_on_mesh; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[i * info->DIMENSION + j] = info->Geo_Node_Coordinate[i * (info->DIMENSION + 1) + j];
					dset_cp_xyz.write(buf_cp_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", *info->Total_Control_Point_on_mesh);
				fprintf(fp_loc_cp, "            %s/local_patch_cp_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_cp_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_cp_ien[2] = { (hsize_t)(*info->Total_Control_Point_on_mesh), 1 };
					DataSpace dataspace_cp_ien(2, dims_cp_ien);
					DSetCreatPropList plist_cp_ien = make_plist(dims_cp_ien, 2, PredType::NATIVE_INT);
					DataSet dset_cp_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_cp_ien, plist_cp_ien);
					vector<int> buf_cp_ien((size_t)(*info->Total_Control_Point_on_mesh));
					for (int i = 0; i < *info->Total_Control_Point_on_mesh; i++)
						buf_cp_ien[i] = i;
					dset_cp_ien.write(buf_cp_ien.data(), PredType::NATIVE_INT);
				}

				fprintf(fp_loc_cp, "        </Grid>\n");
				fprintf(fp_loc_cp, "      </Grid>\n");
			}
			#pragma endregion
			#endif

			#pragma region middle step local patch boundary line
			{
				fprintf(fp_loc_bl, "      <Grid Name=\"step%d\" GridType=\"Collection\">\n", step);
				fprintf(fp_loc_bl, "        <Time Value=\"%.20e\"/>\n", time);
				fprintf(fp_loc_bl, "        <Grid Name=\"ien\">\n");
				fprintf(fp_loc_bl, "          <Topology TopologyType=\"EDGE_3\" NumberOfElements=\"%d\">\n", Total_edge - Total_edge_glo);
				fprintf(fp_loc_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_bl, "          </Topology>\n");
				fprintf(fp_loc_bl, "          <Geometry GeometryType=\"%s\">\n", geometry_type);
				fprintf(fp_loc_bl, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc_bl, "          </Geometry>\n");
				fprintf(fp_loc_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", (Total_edge - Total_edge_glo) * point_on_edge, info->DIMENSION);
				fprintf(fp_loc_bl, "            %s/local_boundary_step%d.h5:/xyz\n", "./bin", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_boundary_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_TRUNC);
					hsize_t dims_boundary_xyz[2] = { (hsize_t)((Total_edge - Total_edge_glo) * point_on_edge), (hsize_t)(info->DIMENSION) };
					DataSpace dataspace_boundary_xyz(2, dims_boundary_xyz);
					DSetCreatPropList plist_bxyz = make_plist(dims_boundary_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_boundary_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_boundary_xyz, plist_bxyz);
					vector<double> buf_boundary_xyz((size_t)(Total_edge - Total_edge_glo) * point_on_edge * info->DIMENSION);
					for (int i = Total_edge_glo; i < Total_edge; i++)
						for (int j = 0; j < point_on_edge; j++)
							for (int k = 0; k < info->DIMENSION; k++)
								buf_boundary_xyz[(i - Total_edge_glo) * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k];
					dset_boundary_xyz.write(buf_boundary_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", Total_edge - Total_edge_glo, point_on_edge);
				fprintf(fp_loc_bl, "            %s/local_boundary_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_boundary_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_boundary_ien[2] = { (hsize_t)(Total_edge - Total_edge_glo), (hsize_t)(point_on_edge) };
					DataSpace dataspace_boundary_ien(2, dims_boundary_ien);
					DSetCreatPropList plist_bien = make_plist(dims_boundary_ien, 2, PredType::NATIVE_INT);
					DataSet dset_boundary_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_boundary_ien, plist_bien);
					vector<int> buf_boundary_ien((size_t)(Total_edge - Total_edge_glo) * point_on_edge);
					for (int i = 0; i < Total_edge - Total_edge_glo; i++)
						for (int j = 0; j < point_on_edge; j++)
							buf_boundary_ien[i * point_on_edge + j] = i * point_on_edge + j;
					dset_boundary_ien.write(buf_boundary_ien.data(), PredType::NATIVE_INT);
				}

				fprintf(fp_loc_bl, "        </Grid>\n");
				fprintf(fp_loc_bl, "      </Grid>\n");
			}
			#pragma endregion
		}
	}
	#pragma endregion

	// final step
	if (isGeometryOnly || time == total_time)
	{
		fprintf(fp_glo, "    </Grid>\n");
		fprintf(fp_glo, "  </Domain>\n");
		fprintf(fp_glo, "</Xdmf>\n");

		fprintf(fp_glo_cp, "    </Grid>\n");
		fprintf(fp_glo_cp, "  </Domain>\n");
		fprintf(fp_glo_cp, "</Xdmf>\n");

		fprintf(fp_glo_bl, "    </Grid>\n");
		fprintf(fp_glo_bl, "  </Domain>\n");
		fprintf(fp_glo_bl, "</Xdmf>\n");

		fprintf(fp_loc, "    </Grid>\n");
		fprintf(fp_loc, "  </Domain>\n");
		fprintf(fp_loc, "</Xdmf>\n");

		
		fprintf(fp_loc_bl, "    </Grid>\n");
		fprintf(fp_loc_bl, "  </Domain>\n");
		fprintf(fp_loc_bl, "</Xdmf>\n");

		// ローカル制御点はグローバルパラメータ空間内に配置されるため、物理座標へ可視化する必要がない。
		// fprintf(fp_loc_cp, "    </Grid>\n");
		// fprintf(fp_loc_cp, "  </Domain>\n");
		// fprintf(fp_loc_cp, "</Xdmf>\n");
	}

	step++;
}


void Make_connectivity(information *info)
{
	int i, j, k, l, m;

	// Make Patch_check
	int Patch_check_counter = 0;
	if (info->DIMENSION == 2)
    {
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];

			// 辺0 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺0 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺2 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺2 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺4 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺4 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺6 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺6 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];
			int temp = info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] * (info->No_Control_point[i * info->DIMENSION + 2] - 1);
			int temp_Patch_check_counter;
			int temp_Patch_check[8];

			// 面0 点0
			temp_Patch_check[0] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 面0 点1
			temp_Patch_check[1] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面0 点2
			temp_Patch_check[2] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面0 点3
			temp_Patch_check[3] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面1 点0
			temp_Patch_check[4] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp];
			Patch_check_counter++;
			// 面1 点1
			temp_Patch_check[5] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面1 点2
			temp_Patch_check[6] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面1 点3
			temp_Patch_check[7] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面2 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点1
			temp_Patch_check_counter = temp_Patch_check[3];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点2
			temp_Patch_check_counter = temp_Patch_check[7];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点1
			temp_Patch_check_counter = temp_Patch_check[1];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点2
			temp_Patch_check_counter = temp_Patch_check[5];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点0
			temp_Patch_check_counter = temp_Patch_check[1];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点3
			temp_Patch_check_counter = temp_Patch_check[5];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点0
			temp_Patch_check_counter = temp_Patch_check[3];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点3
			temp_Patch_check_counter = temp_Patch_check[7];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
		}
	}

	// Check
	int patch_start = 0;
	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
				patch_start = i;

			for (j = patch_start; j < i; j++)
			{
				// crack check
				int continue_flag = 0;
				int mesh_num = 0;
				for (k = 0; k <= Total_mesh; k++)
					if (info->Total_Patch_to_mesh[k] > i)
					{
						mesh_num = k - 1;
						break;
					}

				for (k = info->crack_pair_n_to_mesh[mesh_num]; k < info->crack_pair_n_to_mesh[mesh_num + 1]; k++)
					if (((info->crack_pair[k * 2] == i - patch_start) && (info->crack_pair[k * 2 + 1] == j - patch_start)) || ((info->crack_pair[k * 2] == j - patch_start) && (info->crack_pair[k * 2 + 1] == i - patch_start)))
						continue_flag = 1;

				if (continue_flag)
					continue;

				int own = i * 16;
				int opp = j * 16;
				for (k = 0; k < 4; k++)
				{
					int kk = 2 * k;
					for (l = 0; l < 8; l++)
					{
						// 辺が一致している場合 Face_Edge_info を True
						if (info->Patch_check[own + kk * 2] == info->Patch_check[opp + l * 2] && info->Patch_check[own + kk * 2 + 1] == info->Patch_check[opp + l * 2 + 1])
						{
							info->Face_Edge_info[i * 32 + k * 8 + l] = 1;
							info->Opponent_patch_num[i * 4 + k] = j;
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
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
				patch_start = i;

			for (j = patch_start; j < i; j++)
			{
				// crack check
				int continue_flag = 0;
				int mesh_num = 0;
				for (k = 0; k <= Total_mesh; k++)
					if (info->Total_Patch_to_mesh[k] > i)
					{
						mesh_num = k - 1;
						break;
					}

				for (k = info->crack_pair_n_to_mesh[mesh_num]; k < info->crack_pair_n_to_mesh[mesh_num + 1]; k++)
					if (((info->crack_pair[k * 2] == i - patch_start) && (info->crack_pair[k * 2 + 1] == j - patch_start)) || ((info->crack_pair[k * 2] == j - patch_start) && (info->crack_pair[k * 2 + 1] == i - patch_start)))
						continue_flag = 1;

				if (continue_flag)
					continue;

				int own = i * 24;
				int opp = j * 24;
				for (k = 0; k < 6; k++)
				{
					int kk = k * 4;
					for (m = 0; m < 4; m++)
					{
						check_a[m] = info->Patch_check[own + kk + m];
					}
					sort(4, check_a);

					for (l = 0; l < 6; l++)
					{
						int ll = l * 4;
						for (m = 0; m < 4; m++)
						{
							check_b[m] = info->Patch_check[opp + ll + m];
						}
						sort(4, check_b);

						// 面が一致している場合 Face_Edge_info を Mode 番号 (m) に
						if (check_a[0] == check_b[0] && check_a[1] == check_b[1] && check_a[2] == check_b[2] && check_a[3] == check_b[3])
						{
							for (m = 0; m < 4; m++)
							{
								if (info->Patch_check[own + kk] == info->Patch_check[opp + ll + m])
								{
									info->Face_Edge_info[i * 36 + k * 6 + l] = m;
									info->Opponent_patch_num[i * 6 + k] = j;
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
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (j = 0; j < i; j++)
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
			for (j = 0; j < 4; j++)
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

				for (k = 0; k < 8; k++)
				{
					if (info->Face_Edge_info[Edge_counter + j * 8 + k] == 1)
					{
						int opp = 0;
						for (l = 0; l < info->Opponent_patch_num[i * 4 + j]; l++)
						{
							int a, b;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							opp += 2 * (a + b);
						}

						Edge[j] = 1;

						int p = k / 2;
						int q = k % 2;

						int temp_ii = info->line_No_Total_element[info->Opponent_patch_num[i * 4 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[info->Opponent_patch_num[i * 4 + j] * info->DIMENSION + 1] * 2 + 1;

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
							for (l = 0; l < temp_n; l++)
							{
								info->Patch_array[own + l] = info->Patch_array[opp + l];
							}
							break;
						}
						else if (q == 1)
						{
							for (l = 0; l < temp_n; l++)
							{
								info->Patch_array[own + l] = info->Patch_array[opp + (temp_n - 1) - l];
							}
							break;
						}
					}
				}
			}

			// コネクティビティを作成
			int x, y;

			own = 0;
			for (j = 0; j < i; j++)
			{
				int a, b;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				own += 2 * (a + b);
			}

			int e_x_max = info->line_No_Total_element[i * info->DIMENSION + 0];
			int e_y_max = info->line_No_Total_element[i * info->DIMENSION + 1];
			for (y = 0; y < e_y_max; y++)
			{
				for (x = 0; x < e_x_max; x++)
				{
					for (int point = 0; point < 9; point++)
					{
						// point 0
						if (point == 0)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 1];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 3];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 1
						else if (point == 1)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 2];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 2
						else if (point == 2)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 3
						else if (point == 3)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 2];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 4
						else if (point == 4)
						{
							if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 6];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 5
						else if (point == 5)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 6
						else if (point == 6)
						{
							if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 7
						else if (point == 7)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 5];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 8
						else if (point == 8)
						{
							info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
							info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
							info->Connectivity_point[counter] = point;
							counter++;
						}

						info->Connectivity_all_ele[info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point]].emplace_back(ele_counter + y * e_x_max + x);
					}
				}
			}

			// Patch_array の作ってない分を作成
			int point = 0;
			int xi, eta;
			for (eta = 0; eta < jj; eta++)
			{
				for (xi = 0; xi < ii; xi++)
				{
					int e_x, e_y;

					if (eta == 0 && Edge[0] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + xi] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (eta == jj - 1 && Edge[2] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + ii + jj + xi] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == 0 && Edge[3] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + 2 * ii + jj + eta] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == ii - 1 && Edge[1] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + ii + eta] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
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
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (j = 0; j < i; j++)
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
			for (j = 0; j < 6; j++)
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

				for (k = 0; k < 6; k++)
				{
					if (info->Face_Edge_info[Face_counter + j * 6 + k] >= 0)
					{
						int opp = 0;
						for (l = 0; l < info->Opponent_patch_num[i * 6 + j]; l++)
						{
							int a, b, c;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							c = info->line_No_Total_element[l * info->DIMENSION + 2] * 2 + 1;
							opp += 2 * (a * b + b * c + a * c);
						}

						Face[j] = 1;

						int temp_ii = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 1] * 2 + 1;
						int temp_kk = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 2] * 2 + 1;

						for (l = 1; l <= k; l++)
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
							if (info->Face_Edge_info[Face_counter + j * 6 + k] == 0)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + l * own_a + m];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 1)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + m * own_b + ((own_b - 1) - l)];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 2)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + ((own_b - 1) - l) * own_a + ((own_a - 1) - m)];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 3)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + ((own_a - 1) - m) * own_b + l];
								break;
							}
						}
						else
						{
							if (info->Face_Edge_info[Face_counter + j * 6 + k] == 0)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + m * own_b + l];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k] == 1)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + ((own_a - 1) - m)] = info->Patch_array[opp + l * own_a + m];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k] == 2)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + ((own_a - 1) - m) * own_b + ((own_b - 1) - l)];
								break;
							}
							else if (info->Face_Edge_info[Face_counter + j * 6 + k] == 3)
							{
								for (l = 0; l < own_b; l++)
									for (m = 0; m < own_a; m++)
										info->Patch_array[own + ((own_b - 1) - l) * own_a + m] = info->Patch_array[opp + l * own_a + m];
								break;
							}
						}
					}
				}
			}

			// コネクティビティを作成
			int x, y, z;

			own = 0;
			for (j = 0; j < i; j++)
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
			for (z = 0; z < e_z_max; z++)
			{
				for (y = 0; y < e_y_max; y++)
				{
					for (x = 0; x < e_x_max; x++)
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 1];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 3];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 4];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 2];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 5];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 2];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 5];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 10];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 12];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 13];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 9];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 15];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 13];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 17];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 19];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 18];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 18];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 21];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 23];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 25];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
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
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 26
							else if (point == 26)
							{
								info->Connectivity[temp] = counter;
								info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}

							info->Connectivity_all_ele[info->Connectivity[temp]].emplace_back(ele_counter + z * e_y_max * e_x_max + y * e_x_max + x);
						}
					}
				}
			}

			// Patch_array の作ってない分を作成
			int point = 0;
			int xi, eta, zeta;
			for (zeta = 0; zeta < kk; zeta++)
			{
				for (eta = 0; eta < jj; eta++)
				{
					for (xi = 0; xi < ii; xi++)
					{
						int e_x, e_y, e_z;

						if (zeta == 0 && Face[0] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + eta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (zeta == kk - 1 && Face[1] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp1 + eta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == 0 && Face[3] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp3 + zeta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == jj - 1 && Face[5] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp5 + zeta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == 0 && Face[2] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp2 + zeta * jj + eta] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == ii - 1 && Face[4] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp4 + zeta * jj + eta] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
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

	// sort Connectivity_all_ele
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Total_connectivity; i++)
		sortUniqueErase(info->Connectivity_all_ele[i]);

	// sort ecn
	int point_on_element = pow_int(3, info->DIMENSION);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<int> temp_ecn;
		for (int j = 0; j < point_on_element; j++)
		{
			int cnt = info->Connectivity[i * point_on_element + j];
			for (size_t k = 0; k < info->Connectivity_all_ele[cnt].size(); k++)
			{
				temp_ecn.emplace_back(info->Connectivity_all_ele[cnt][k]);
			}
		}

		// sort
		sortUniqueErase(temp_ecn);

		info->ecn[i] = std::move(temp_ecn);
	}
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


void Init_viewer_info(information *info)
{
	#pragma omp parallel for
	for (int i = 0; i < Total_connectivity; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
		{
			info->Connectivity_coord[i * info->DIMENSION + j] = 0.0;
			info->disp_at_connectivity[i * info->DIMENSION + j] = 0.0;
		}
		for (int j = 0; j < D_MATRIX_SIZE; j++)
		{
			info->strain_at_connectivity[i * D_MATRIX_SIZE + j] = 0.0;
			info->stress_at_connectivity[i * D_MATRIX_SIZE + j] = 0.0;
		}
		info->vm_at_connectivity[i] = 0.0;
		info->hs_at_connectivity[i] = 0.0;
	}
}


void Make_average_points(vector<average_points> &ap, information *info)
{
	int point_on_element = pow_int(3, info->DIMENSION);
	vector<double> point_array(point_on_element * info->DIMENSION);

	// Make point in element
	int counter = 0;
	if (info->DIMENSION == 2)
	{
		// 0, 2, 8, 6, 1, 5, 7, 3, 4

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 4
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 6
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;
	}
	else if (info->DIMENSION == 3)
	{
		// 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 13

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 4
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 6
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 9
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 10
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 11
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 12
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 13
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 14
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 15
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 16
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 17
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 18
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 19
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 20
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 21
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 22
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 23
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 24
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 25
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 26
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;
	}

	#pragma omp parallel for
	for (int i = 0; i < Total_connectivity; i++)
	{
		int ce = info->Connectivity_ele[i];
		int point = info->Connectivity_point[i];

		int NG_mod = info->c.NUM_GAUSS_POINTS % 2;
		int point_in_ele = 0;
		// 頂点
		if (0 <= point && point <= 7)
		{
			point_in_ele = 1;
		}
		// 辺中間点
		else if (8 <= point && point <= 13)
		{
			if (NG_mod == 0)
			{
				point_in_ele = 2;
			}
			else if (NG_mod == 1)
			{
				point_in_ele = 1;
			}
		}
		// 面中間点
		else if (14 <= point && point <= 25)
		{
			if (NG_mod == 0)
			{
				point_in_ele = 4;
			}
			else if (NG_mod == 1)
			{
				point_in_ele = 1;
			}
		}
		// 中間点
		else if (point == 26)
		{
			if (NG_mod == 0)
			{
				point_in_ele = 8;
			}
			else if (NG_mod == 1)
			{
				point_in_ele = 1;
			}
		}

		int average_point_n = point_in_ele * info->Connectivity_all_ele[i].size();

		// make Connectivity_coord
		double temp_point[MAX_DIMENSION] = {0.0};
		for (int j = 0; j < info->DIMENSION; j++)
		{
			info->Connectivity_coord[i * info->DIMENSION + j] = 0.0;
			temp_point[j] = point_array[point * info->DIMENSION + j];
		}

		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), temp_point, ce, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ce]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
				info->Connectivity_coord[i * info->DIMENSION + k] += R[j] * info->Node_Coordinate[info->Controlpoint_of_Element[ce * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];

		vector<distance_norm> dn;
		for (size_t j = 0; j < info->Connectivity_all_ele[i].size(); j++)
		{
			int ele = info->Connectivity_all_ele[i][j];
			for (int k = 0; k < info->gp[ele].n(); k++)
			{
				distance_norm temp_dn;
				temp_dn.e = ele;
				temp_dn.gp_num = k;
				for (int l = 0; l < info->DIMENSION; l++)
				{
					temp_dn.r += pow(info->gp[ele].coord()[k * info->DIMENSION + l] - info->Connectivity_coord[i * info->DIMENSION + l], 2);
				}
				temp_dn.r = sqrt(temp_dn.r);
				dn.emplace_back(temp_dn);
			}
		}

		// sort
		sort(dn.begin(), dn.end());

		// make average points
		ap[i].n = average_point_n;
		double total_w = 0.0;
		for (int j = 0; j < ap[i].n; j++)
		{
			// total_w += 1.0 / dn[j].r;
			total_w += 1.0 / pow(dn[j].r, 2);
		}
		for (int j = 0; j < ap[i].n; j++)
		{
			ap[i].target_e.emplace_back(dn[j].e);
			ap[i].target_gp_num.emplace_back(dn[j].gp_num);
			// ap[i].w.emplace_back(1.0 / dn[j].r / total_w);
			ap[i].w.emplace_back(1.0 / pow(dn[j].r, 2) / total_w);
		}
	}
}


void Make_info_for_viewer_by_shape_func(information *info)
{
	int point_on_element = pow_int(3, info->DIMENSION);
	vector<double> point_array(point_on_element * info->DIMENSION);

	// Make point in element
	int counter = 0;
	if (info->DIMENSION == 2)
	{
		// 0, 2, 8, 6, 1, 5, 7, 3, 4

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 4
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 6
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;
	}
	else if (info->DIMENSION == 3)
	{
		// 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 13

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 4
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 6
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 9
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 10
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 11
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 12
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 13
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 14
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 15
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 16
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 17
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 18
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 19
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 20
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 21
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 22
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 23
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 24
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 25
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 26
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Total_connectivity; i++)
	{
		double temp_point[MAX_DIMENSION] = {0.0};
		double temp_point_glo[MAX_DIMENSION] = {0.0};
		double temp_para_glo[MAX_DIMENSION] = {0.0};
		double temp_point_loc[MAX_DIMENSION] = {0.0};
		double temp_para_loc[MAX_DIMENSION] = {0.0};

		int element = info->Connectivity_ele[i];
		int point = info->Connectivity_point[i];

		// for SSIGA local
		double geo_coord_tilde[MAX_DIMENSION] = {0.0};	//ローカルジオメトリ表現の要素パラメータ座標
		int geo_element = 0;							//ローカルジオメトリ表現の要素番号

		// make temp_point
		for (int j = 0; j < info->DIMENSION; j++)
			temp_point[j] = point_array[point * info->DIMENSION + j];

		// make physical coordinate, disp_at_connectivity
		if (i < Total_connectivity_glo)
		{
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
			vector<double> u(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
			shape_and_dshape(R.data(), temp_point, element, info);
			Make_B_Linear(element, temp_point, bl.data(), info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					info->Connectivity_coord[i * info->DIMENSION + k] += R[j] * info->Node_Coordinate[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
					info->disp_at_connectivity[i * info->DIMENSION + k] += R[j] * d;
					u[j * info->DIMENSION + k] = d;
				}

			// strain
			for (int j = 0; j < D_MATRIX_SIZE; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						info->strain_at_connectivity[i * N_STRESS + j] += bl[j * MAX_KIEL_SIZE + k * info->DIMENSION + l] * u[k * info->DIMENSION + l];
		}
		if (i >= Total_connectivity_glo && Total_mesh > 1)
		{
			int patch_num = info->Element_patch[element];							//要素番号からパッチ番号を取得
			int geo_patch_num = patch_num - info->Total_Patch_on_mesh[0];			//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい
			double coord[MAX_DIMENSION] = {0.0};									//ローカルのパッチパラメータ座標
			double coord_tilde[MAX_DIMENSION] = {0.0};								//ローカルの要素パラメータ座標
			double glo_coord_tilde[MAX_DIMENSION] = {0.0};							//グローバルの要素パラメータ座標
			double Connectivity_coord_para[MAX_DIMENSION] = {0.0};					//要素のパラメータ座標	

			std::copy(temp_point, temp_point + info->DIMENSION, coord_tilde);				//要素のパラメータ座標をガウス点座標に変換
			trans_ele_patch_coord(coord, coord_tilde, patch_num, element, info);			//要素のパラメータ座標をパッチのパラメータ座標に変換
			geo_element = geo_ele_check(geo_patch_num, coord, info);						//パラメータ座標から、要素を探索
			geo_tilde_coord(geo_coord_tilde, coord, geo_patch_num, geo_element, info);		//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換

			vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_disp(MAX_NO_CP_ON_ELEMENT);
			vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
			vector<double> u(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
			geo_shape_and_dshape(R_geo.data(), geo_coord_tilde, geo_element, info);
			Bspline_shape_and_dshape(R_disp.data(), temp_point, element, info);
			Make_B_Linear(element, temp_point, bl.data(), info);

			// physical coordinate
			for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[geo_element]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
					Connectivity_coord_para[k] += R_geo[j] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_element * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
			int ele_glo = ele_check(info->Global_local_patch, Connectivity_coord_para, info);
			tilde_coord(glo_coord_tilde, Connectivity_coord_para, info->Global_local_patch, ele_glo, info);
			physical_coord(ele_glo, glo_coord_tilde, &info->Connectivity_coord[i * info->DIMENSION], info);

			// displacement
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					info->disp_at_connectivity[i * info->DIMENSION + k] += R_disp[j] * d;
					u[j * info->DIMENSION + k] = d;
				}

			// strain
			for (int j = 0; j < D_MATRIX_SIZE; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						info->strain_at_connectivity[i * N_STRESS + j] += bl[j * MAX_KIEL_SIZE + k * info->DIMENSION + l] * u[k * info->DIMENSION + l];
		}
		
		// overlay global
		int status_glo_overlay = 0;
		if (i < Total_connectivity_glo && Total_mesh > 1 && info->Element_patch[element] == info->Global_local_patch)
		{
			double temp_point_patch[MAX_DIMENSION] = {0.0};
			trans_ele_patch_coord(temp_point_patch, temp_point, info->Element_patch[element], element, info);

			// make temp_point_glo
			int itr_n = 0, loc_patch = 0;
			for (int j = info->Total_Patch_to_mesh[1]; j < info->Total_Patch_to_mesh[Total_mesh]; j++)
			{
				itr_n = calc_local_patch_parameter_coord(temp_point_patch, j, temp_para_loc, info);
				loc_patch = j;
				if (itr_n != ERROR)
				{
					status_glo_overlay = 1;
					break;
				}
			}

			int element_loc = 0;
			if (status_glo_overlay)
			{
				element_loc = ele_check(loc_patch, temp_para_loc, info);
				tilde_coord(temp_point_loc, temp_para_loc, loc_patch, element_loc, info);
			}

			// overlay displacement
			if (status_glo_overlay)
			{
				vector<double> R_loc(MAX_NO_CP_ON_ELEMENT);
				vector<double> bl_loc(D_MATRIX_SIZE * MAX_KIEL_SIZE);
				vector<double> u_loc(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
				Bspline_shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
				Make_B_Linear(element_loc, temp_point_loc, bl_loc.data(), info);
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; j++)
					for (int k = 0; k < info->DIMENSION; k++)
					{
						double d_loc = info->Displacement[info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
						info->disp_at_connectivity[i * info->DIMENSION + k] += R_loc[j] * d_loc;
						u_loc[j * info->DIMENSION + k] = d_loc;
					}

				// strain
				for (int j = 0; j < D_MATRIX_SIZE; j++)
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
							info->strain_at_connectivity[i * N_STRESS + j] += bl_loc[j * MAX_KIEL_SIZE + k * info->DIMENSION + l] * u_loc[k * info->DIMENSION + l];
			}
		}

		// overlay local
		if (i >= Total_connectivity_glo && Total_mesh > 1)
		{
			geo_parameter_coord(geo_element, geo_coord_tilde, temp_para_glo, info);
			int element_glo = ele_check(info->Global_local_patch, temp_para_glo, info);
			tilde_coord(temp_point_glo, temp_para_glo, info->Global_local_patch, element_glo, info);

			// overlay displacement
			vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
			vector<double> bl_glo(D_MATRIX_SIZE * MAX_KIEL_SIZE);
			vector<double> u_glo(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
			shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
			Make_B_Linear(element_glo, temp_point_glo, bl_glo.data(), info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					info->disp_at_connectivity[i * info->DIMENSION + k] += R_glo[j] * d_glo;
					u_glo[j * info->DIMENSION + k] = d_glo;
				}
			
			// strain
			for (int j = 0; j < D_MATRIX_SIZE; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						info->strain_at_connectivity[i * N_STRESS + j] += bl_glo[j * MAX_KIEL_SIZE + k * info->DIMENSION + l] * u_glo[k * info->DIMENSION + l];
		}

		// stress
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				info->stress_at_connectivity[i * N_STRESS + j] += info->D[j * D_MATRIX_SIZE + k] * info->strain_at_connectivity[i * N_STRESS + k];

		// 2D specific adjustments
		if (info->DIMENSION == 2)
		{
			// plane stress condition
			if (info->c.ANALYSIS_MODE == 0)
			{
				info->strain_at_connectivity[i * N_STRESS + 3] = -1.0 * nu / E * (info->stress_at_connectivity[i * N_STRESS + 0] + info->stress_at_connectivity[i * N_STRESS + 1]);
				info->stress_at_connectivity[i * N_STRESS + 3] = 0.0;
			}

			// plane strain condition
			else if (info->c.ANALYSIS_MODE == 1)
			{
				info->strain_at_connectivity[i * N_STRESS + 3] = 0.0;
				info->stress_at_connectivity[i * N_STRESS + 3] = nu * (info->stress_at_connectivity[i * N_STRESS + 0] + info->stress_at_connectivity[i * N_STRESS + 1]);
			}
		}

		auto calc_equivalent_stress = [&](double *stress) {
			double s_eq = sqrt(0.5 * (pow(stress[0] - stress[1], 2) + pow(stress[1] - stress[2], 2) + pow(stress[2] - stress[0], 2) + 6.0 * (pow(stress[3], 2) + pow(stress[4], 2) + pow(stress[5], 2))));
			return s_eq;
		};

		// stress 6 components
		vector<double> stress_6(6, 0.0);
		if (info->DIMENSION == 2)
		{
			double *ptr = &info->stress_at_connectivity[i * N_STRESS];
			stress_6[0] = ptr[0]; stress_6[1] = ptr[1]; stress_6[2] = ptr[3]; stress_6[3] = ptr[2];
		}
		else if (info->DIMENSION == 3)
			for (int j = 0; j < N_STRESS; j++)
				stress_6[j] = info->stress_at_connectivity[i * N_STRESS + j];

		// equivalent stress
		info->vm_at_connectivity[i] = calc_equivalent_stress(stress_6.data());

		// hydrostatic stress
		double hs = (stress_6[0] + stress_6[1] + stress_6[2]) / 3.0;
		info->hs_at_connectivity[i] = hs;
	}
}


// グローバルパラメータ空間による形状関数の計算
void Make_info_for_viewer_by_shape_func_at_global_parameter_space(information *info)
{
	int point_on_element = pow_int(3, info->DIMENSION);
	vector<double> point_array(point_on_element * info->DIMENSION);

	// Make point in element
	int counter = 0;
	if (info->DIMENSION == 2)
	{
		// 0, 2, 8, 6, 1, 5, 7, 3, 4

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 4
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 6
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;
	}
	else if (info->DIMENSION == 3)
	{
		// 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 13

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 4
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 6
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 9
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 10
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 11
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 12
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 13
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 14
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 15
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 16
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 17
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 18
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 19
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 20
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 21
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 22
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 23
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 24
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 25
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 26
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Total_connectivity; i++)
	{
		double temp_point[MAX_DIMENSION] = {0.0};

		int element = info->Connectivity_ele[i];
		int point = info->Connectivity_point[i];

		// for SSIGA local
		double geo_coord_tilde[MAX_DIMENSION] = {0.0};	//ローカルジオメトリ表現の要素パラメータ座標
		int geo_element = 0;							//ローカルジオメトリ表現の要素番号

		// make temp_point
		for (int j = 0; j < info->DIMENSION; j++)
			temp_point[j] = point_array[point * info->DIMENSION + j];

		// make Connectivity_coord
		if (i < Total_connectivity_glo)
		{
			vector<double> para(info->DIMENSION);
			trans_ele_patch_coord(para.data(), temp_point, info->Element_patch[element], element, info);
			for (int j = 0; j < info->DIMENSION; j++)
				info->Connectivity_coord[i * info->DIMENSION + j] += para[j];
		}
		if (i >= Total_connectivity_glo && Total_mesh > 1)
		{
			int patch_num = info->Element_patch[element];							//要素番号からパッチ番号を取得
			int geo_patch_num = patch_num - info->Total_Patch_on_mesh[0];			//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい
			double coord[MAX_DIMENSION] = {0.0};									//ローカルのパッチパラメータ座標
			double coord_tilde[MAX_DIMENSION] = {0.0};								//ローカルの要素パラメータ座標
			
			std::copy(temp_point, temp_point + info->DIMENSION, coord_tilde);				//要素のパラメータ座標をガウス点座標に変換
			trans_ele_patch_coord(coord, coord_tilde, patch_num, element, info);			//要素のパラメータ座標をパッチのパラメータ座標に変換
			geo_element = geo_ele_check(geo_patch_num, coord, info);						//パラメータ座標から、要素を探索
			geo_tilde_coord(geo_coord_tilde, coord, geo_patch_num, geo_element, info);		//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換

			vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
			vector<double> R_disp(MAX_NO_CP_ON_ELEMENT);
			vector<double> bl(D_MATRIX_SIZE * MAX_KIEL_SIZE);
			vector<double> u(MAX_KIEL_SIZE * info->DIMENSION, 0.0);
			geo_shape_and_dshape(R_geo.data(), geo_coord_tilde, geo_element, info);
			Bspline_shape_and_dshape(R_disp.data(), temp_point, element, info);
			Make_B_Linear(element, temp_point, bl.data(), info);

			// physical coordinate
			for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[geo_element]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
					info->Connectivity_coord[i * info->DIMENSION + k] += R_geo[j] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_element * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
		}
	}
}

void Make_boundary_line(information *info)
{
	int i, j, k;
	int point_on_edge = 3;

	int counter = 0;
	double *E_ptr, *C_ptr, *D_ptr;
	E_ptr = info->Edge_coord;
	C_ptr = info->Connectivity_coord;
	D_ptr = info->disp_at_connectivity;

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		if (info->DIMENSION == 2)
		{
			int point_on_element = 9;
			int x = info->ENC[i * info->DIMENSION + 0];
			int y = info->ENC[i * info->DIMENSION + 1];
			int x_max = info->line_No_Total_element[info->Element_patch[i] * info->DIMENSION + 0] - 1;
			int y_max = info->line_No_Total_element[info->Element_patch[i] * info->DIMENSION + 1] - 1;

			if (x == 0)
			{
				int POINT[3] = {0, 3, 7};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == x_max)
			{
				int POINT[3] = {1, 2, 5};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == 0)
			{
				int POINT[3] = {0, 1, 4};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == y_max)
			{
				int POINT[3] = {3, 2, 6};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
		}
		else if (info->DIMENSION == 3)
		{
			int point_on_element = 27;
			int x = info->ENC[i * info->DIMENSION + 0];
			int y = info->ENC[i * info->DIMENSION + 1];
			int z = info->ENC[i * info->DIMENSION + 2];
			int x_max = info->line_No_Total_element[info->Element_patch[i] * info->DIMENSION + 0] - 1;
			int y_max = info->line_No_Total_element[info->Element_patch[i] * info->DIMENSION + 1] - 1;
			int z_max = info->line_No_Total_element[info->Element_patch[i] * info->DIMENSION + 2] - 1;

			if (x == 0 && y == 0)
			{
				int POINT[3] = {0, 4, 16};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == x_max && y == 0)
			{
				int POINT[3] = {1, 5, 17};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == 0 && y == y_max)
			{
				int POINT[3] = {3, 7, 19};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == x_max && y == y_max)
			{
				int POINT[3] = {2, 6, 18};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == 0 && z == 0)
			{
				int POINT[3] = {0, 3, 11};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == x_max && z == 0)
			{
				int POINT[3] = {1, 2, 9};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == 0 && z == z_max)
			{
				int POINT[3] = {4, 7, 15};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (x == x_max && z == z_max)
			{
				int POINT[3] = {5, 6, 13};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == 0 && z == 0)
			{
				int POINT[3] = {0, 1, 8};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == y_max && z == 0)
			{
				int POINT[3] = {3, 2, 10};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == 0 && z == z_max)
			{
				int POINT[3] = {4, 5, 12};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
			if (y == y_max && z == z_max)
			{
				int POINT[3] = {7, 6, 14};
				for (j = 0; j < point_on_edge; j++)
					for (k = 0; k < info->DIMENSION; k++)
						E_ptr[counter * point_on_edge * info->DIMENSION + j * info->DIMENSION + k] = C_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k] + D_ptr[info->Connectivity[i * point_on_element + POINT[j]] * info->DIMENSION + k];
				counter++;
			}
		}

		if (i == info->Total_Element_to_mesh[1] - 1)
			Total_edge_glo = counter;
	}
	Total_edge = counter;
}