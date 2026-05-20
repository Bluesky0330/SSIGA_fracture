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
	int i, j;

	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
		info->Displacement[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = info->Value_of_Constraint[i];

	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			int index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
				info->Displacement[i * info->DIMENSION + j] = info->sol_vec[index];
		}
	}
}


void Substitute_Displacement(information *info, char **argv)
{
	FILE *disp_fp;
	disp_fp = fopen(argv[Total_mesh + 1], "r");

	int temp_i;
	double temp_d;

	fscanf(disp_fp, "%d", &temp_i);
	if (temp_i != info->Total_Control_Point_to_mesh[Total_mesh])
	{
		printf("wrong input, check POST_ONLY");
		exit(1);
	}

	fscanf(disp_fp, "%d", &temp_i);
	if (temp_i != info->DIMENSION)
	{
		printf("wrong input, check POST_ONLY");
		exit(1);
	}

	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		fscanf(disp_fp, "%lf", &temp_d);
		info->Displacement[i] = temp_d;
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


void calc_at_GP(information *info)
{
	gp_switch(false, info);
	int gp_n = pow_int(info->c.NG, info->DIMENSION);

	int i, j, k, l, m;

	double data_result_shape[MAX_DIMENSION];
	double temp_point_glo[MAX_DIMENSION];
	double temp_para_glo[MAX_DIMENSION];

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *temp_disp = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	double *temp_disp_glo = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	FILE *fp_disp, *fp_strain, *fp_stress, *fp_rf, *fp_pc;
	fp_disp   = fopen("_Displacement_overlay_at_GP.dat", "w");
	fp_strain = fopen("_Strain_overlay_at_GP.dat", "w");
	fp_stress = fopen("_Stress_overlay_at_GP.dat", "w");
	fp_rf     = fopen("_ReactoinForce.dat", "w");
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

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int e = info->real_element[i], element_glo = 0;

		for (j = 0; j < info->gp[e].n(); j++)
		{
			// B matrix
			Make_B_Matrix_anypoint(e, B, &info->gp[e].para()[j * info->DIMENSION], info);

			for (k = 0; k < info->DIMENSION; k++)
				data_result_shape[k] = 0.0;

			// make displacement, temp_disp
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &info->gp[e].para()[j * info->DIMENSION], e, info);
			for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				for (l = 0; l < info->DIMENSION; l++)
				{
					double d = info->Displacement[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
					data_result_shape[l] += R[k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
					info->Displacement_at_GP[e * gp_n * info->DIMENSION + j * info->DIMENSION + l] += R[k] * d;
					temp_disp[k * info->DIMENSION + l] = d;
				}

			// overlay
			if (i >= info->Total_Element_on_mesh[0])
			{
				// make temp_point_glo
				int itr_n = 0, glo_patch = 0;
				for (k = 0; k < info->Total_Patch_on_mesh[0]; k++)
				{
					itr_n = calc_patch_parameter_coord(data_result_shape, k, temp_para_glo, info);
					glo_patch = k;
					if (itr_n != ERROR)
						break;
				}

				element_glo = ele_check(glo_patch, temp_para_glo, info);
				tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

				// BG matrix
				Make_B_Matrix_anypoint(element_glo, BG, temp_point_glo, info);

				// overlay displacement, make temp_disp_glo
				vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
				shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; k++)
					for (l = 0; l < info->DIMENSION; l++)
					{
						double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_GP[e * gp_n * info->DIMENSION + j * info->DIMENSION + l] += R_glo[k] * d_glo;
						temp_disp_glo[k * info->DIMENSION + l] = d_glo;
					}
			}

			int ii = e * gp_n + j;

			// strain
			int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[e]] * info->DIMENSION;
			for (k = 0; k < D_MATRIX_SIZE; k++)
				for (l = 0; l < KIEL_SIZE; l++)
					info->Strain_at_GP[ii * N_STRAIN + k] += B[k * MAX_KIEL_SIZE + l] * temp_disp[l];

			// overlay strain
			if (i >= info->Total_Element_on_mesh[0])
			{
				int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]] * info->DIMENSION;
				for (k = 0; k < D_MATRIX_SIZE; k++)
					for (l = 0; l < KIEL_SIZE_glo; l++)
						info->Strain_at_GP[ii * N_STRAIN + k] += BG[k * MAX_KIEL_SIZE + l] * temp_disp_glo[l];
			}

			// stress
			for (k = 0; k < D_MATRIX_SIZE; k++)
				for (l = 0; l < D_MATRIX_SIZE; l++)
					info->Stress_at_GP[ii * N_STRESS + k] += info->D[k * D_MATRIX_SIZE + l] * info->Strain_at_GP[ii * N_STRAIN + l];

			// if DIMENSION == 2, make strain zz, stress zz
			if (info->DIMENSION == 2)
			{
				// 平面応力状態
				if (info->c.DM == 0)
					info->Strain_at_GP[ii * N_STRAIN + 3] = - 1.0 * nu / E * (info->Stress_at_GP[ii * N_STRESS + 0] + info->Stress_at_GP[ii * N_STRESS + 1]);
				// 平面ひずみ状態
				else if (info->c.DM == 1)
					info->Stress_at_GP[ii * N_STRESS + 3] = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu) * (info->Strain_at_GP[ii * N_STRAIN + 0] + info->Strain_at_GP[ii * N_STRAIN + 1]);
			}

			// Reaction Force
			double J = Make_Jac_anypoint(e, &info->gp[e].para()[j * info->DIMENSION], info);
			for (k = 0; k < D_MATRIX_SIZE; k++)
				for (l = 0; l < info->DIMENSION; l++)
					for (m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
						info->ReactionForce[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l] += B[k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->Stress_at_GP[ii * N_STRESS + k] * info->gp[e].w()[j] * J;

			// output file
			// displacement
			fprintf(fp_disp, "%d\t%d", e, j);
			for (k = 0; k < info->DIMENSION; k++)
				fprintf(fp_disp, "\t%.15e", info->Displacement_at_GP[ii * info->DIMENSION + k]);
			fprintf(fp_disp, "\n");

			// strain
			fprintf(fp_strain, "%d\t%d", e, j);
			for (k = 0; k < N_STRAIN; k++)
				fprintf(fp_strain, "\t%.15e", info->Strain_at_GP[ii * N_STRAIN + k]);
			fprintf(fp_strain, "\n");

			// stress
			fprintf(fp_stress, "%d\t%d", e, j);
			for (k = 0; k < N_STRESS; k++)
				fprintf(fp_stress, "\t%.15e", info->Stress_at_GP[ii * N_STRESS + k]);
			fprintf(fp_stress, "\n");

			// physical coordinate
			fprintf(fp_pc, "%d\t%d", e, j);
			for (k = 0; k < info->DIMENSION; k++)
				fprintf(fp_pc, "\t%.15e", data_result_shape[k]);
			fprintf(fp_pc, "\n");
		}
	}

	// Reaction Force
	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp_rf, "%d", i);
		for (j = 0; j < info->DIMENSION; j++)
			fprintf(fp_rf, "\t%.15e", info->ReactionForce[i * info->DIMENSION + j]);
		fprintf(fp_rf, "\n");
	}
}


void calc_at_ele_vertex(information *info)
{
	int vertex_n = pow_int(2, info->DIMENSION);
	double *point_array = (double *)malloc(sizeof(double) * vertex_n * info->DIMENSION);

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
		point_array[counter] = -1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] =  1.0;	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] = -1.0;	point_array[counter + 2] =  1.0;	counter += 3;
		point_array[counter] = -1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] =  1.0;	counter += 3;
		point_array[counter] =  1.0;	point_array[counter + 1] =  1.0;	point_array[counter + 2] =  1.0;
	}

	int i, j, k, l;

	double temp_point[MAX_DIMENSION];
	double data_result_shape[MAX_DIMENSION];
	double temp_point_glo[MAX_DIMENSION];
	double temp_para_glo[MAX_DIMENSION];

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *temp_disp = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	double *temp_disp_glo = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

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

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int e = info->real_element[i], element_glo = 0;

		for (j = 0; j < vertex_n; j++)
		{
			// make temp_point
			for (k = 0; k < info->DIMENSION; k++)
				temp_point[k] = point_array[j * info->DIMENSION + k];

			// B matrix
			Make_B_Matrix_anypoint(e, B, temp_point, info);

			for (k = 0; k < info->DIMENSION; k++)
				data_result_shape[k] = 0.0;

			// make displacement, temp_disp
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), temp_point, e, info);
			for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				for (l = 0; l < info->DIMENSION; l++)
				{
					double d = info->Displacement[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
					data_result_shape[l] += R[k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
					info->Displacement_at_ele_vertex[e * vertex_n * info->DIMENSION + j * info->DIMENSION + l] += R[k] * d;
					temp_disp[k * info->DIMENSION + l] = d;
				}

			// overlay
			if (i >= info->Total_Element_on_mesh[0])
			{
				// make temp_point_glo
				int itr_n = 0, glo_patch = 0;
				for (k = 0; k < info->Total_Patch_on_mesh[0]; k++)
				{
					itr_n = calc_patch_parameter_coord(data_result_shape, k, temp_para_glo, info);
					glo_patch = k;
					if (itr_n != ERROR)
						break;
				}

				element_glo = ele_check(glo_patch, temp_para_glo, info);
				tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

				// BG matrix
				Make_B_Matrix_anypoint(element_glo, BG, temp_point_glo, info);

				// overlay displacement, make temp_disp_glo
				vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
				shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; k++)
					for (l = 0; l < info->DIMENSION; l++)
					{
						double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						info->Displacement_at_ele_vertex[e * vertex_n * info->DIMENSION + j * info->DIMENSION + l] += R_glo[k] * d_glo;
						temp_disp_glo[k * info->DIMENSION + l] = d_glo;
					}
			}

			int ii = e * vertex_n + j;

			// strain
			int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[e]] * info->DIMENSION;
			for (k = 0; k < D_MATRIX_SIZE; k++)
				for (l = 0; l < KIEL_SIZE; l++)
					info->Strain_at_ele_vertex[ii * N_STRAIN + k] += B[k * MAX_KIEL_SIZE + l] * temp_disp[l];

			// overlay strain
			if (i >= info->Total_Element_on_mesh[0])
			{
				int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]] * info->DIMENSION;
				for (k = 0; k < D_MATRIX_SIZE; k++)
					for (l = 0; l < KIEL_SIZE_glo; l++)
						info->Strain_at_ele_vertex[ii * N_STRAIN + k] += BG[k * MAX_KIEL_SIZE + l] * temp_disp_glo[l];
			}

			// stress
			for (k = 0; k < D_MATRIX_SIZE; k++)
				for (l = 0; l < D_MATRIX_SIZE; l++)
					info->Stress_at_ele_vertex[ii * N_STRESS + k] += info->D[k * D_MATRIX_SIZE + l] * info->Strain_at_ele_vertex[ii * N_STRAIN + l];

			// if DIMENSION == 2, make strain zz, stress zz
			if (info->DIMENSION == 2)
			{
				// 平面応力状態
				if (info->c.DM == 0)
					info->Strain_at_ele_vertex[ii * N_STRAIN + 3] = - 1.0 * nu / E * (info->Stress_at_ele_vertex[ii * N_STRESS + 0] + info->Stress_at_ele_vertex[ii * N_STRESS + 1]);
				// 平面ひずみ状態
				else if (info->c.DM == 1)
					info->Stress_at_ele_vertex[ii * N_STRESS + 3] = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu) * (info->Strain_at_ele_vertex[ii * N_STRAIN + 0] + info->Strain_at_ele_vertex[ii * N_STRAIN + 1]);
			}

			// output file
			// displacement
			fprintf(fp_disp, "%d\t%d", e, j);
			for (k = 0; k < info->DIMENSION; k++)
				fprintf(fp_disp, "\t%.15e", info->Displacement_at_ele_vertex[ii * info->DIMENSION + k]);
			fprintf(fp_disp, "\n");

			// strain
			fprintf(fp_strain, "%d\t%d", e, j);
			for (k = 0; k < N_STRAIN; k++)
				fprintf(fp_strain, "\t%.15e", info->Strain_at_ele_vertex[ii * N_STRAIN + k]);
			fprintf(fp_strain, "\n");

			// stress
			fprintf(fp_stress, "%d\t%d", e, j);
			for (k = 0; k < N_STRESS; k++)
				fprintf(fp_stress, "\t%.15e", info->Stress_at_ele_vertex[ii * N_STRESS + k]);
			fprintf(fp_stress, "\n");

			// physical coordinate
			fprintf(fp_pc, "%d\t%d", e, j);
			for (k = 0; k < info->DIMENSION; k++)
				fprintf(fp_pc, "\t%.15e", data_result_shape[k]);
			fprintf(fp_pc, "\n");
		}
	}
}


// output SVG K matrix
void K_output_svg(information *info)
{
	// [K] = [[K^G, K^GL], [K^GL, K^L]]

	int i, j;
	int ndof = K_Whole_Size;

	char color_vec[2][10] = {"#f5f5f5", "#ee82ee"};
	// 0	whitesmoke
	// 1	violet
	// https://www.colordic.org/

	double space = 3.0, scale = 1000.0 / (((double)ndof) + 2.0 * space);

	double width = (((double)ndof) + 2.0 * space) * scale;
	double height = width;

	char str[256] = "K_matrix.svg";
	FILE *fp = fopen(str, "w");

	fprintf(fp, "<?xml version='1.0'?>\n");
	fprintf(fp, "<svg width='%.16e' height='%.16e' version='1.1' style='background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height);

	double xx = space * scale;
	double yy = space * scale;
	double ww = ndof * scale;
	double hh = ndof * scale;
	fprintf(fp, "<rect x='%.16e' y='%.16e' width='%.16e' height='%.16e' fill='%s' />\n", xx, yy, ww, hh, color_vec[0]);

	vector<int> K_bool(ndof); // 一行分保存する

	// 各行の成分を抽出
	for (i = 0; i < ndof; i++)
	{
		for (j = 0; j < ndof; j++)
			K_bool[j] = 0;

		for (j = 0; j < ndof; j++)
		{
			int temp_count = 0;
			if (i <= j)
				temp_count = RowCol_to_icount(i, j, info);
			else if (i > j)
				temp_count = RowCol_to_icount(j, i, info);

			if (temp_count != -1)
				K_bool[j] = 1;
		}

		for (j = 0; j < ndof; j++)
		{
			double x = (((double)j) + space) * scale;
			double y = (((double)i) + space) * scale;

			if (K_bool[j] == 1)
				fprintf(fp, "<rect x='%.16e' y='%.16e' width='%.16e' height='%.16e' fill='%s' />\n", x, y, scale, scale, color_vec[1]);
		}
	}

	fprintf(fp, "</svg>");
	fclose(fp);
}


// output octree subcell
void Output_Octree_subcell(information *info)
{
	int point_on_element = pow_int(3, info->DIMENSION);
	double *point_array = (double *)malloc(sizeof(double) * point_on_element * info->DIMENSION);

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

	int total_subcell = 0;
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		total_subcell += info->octree_subcell[i].size();

	vector<vector<double>> out_octree_coord(info->Total_Element_to_mesh[Total_mesh]);

	#pragma omp parallel for
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		for (size_t j = 0; j < info->octree_subcell[i].size(); j++)
			for (int k = 0; k < point_on_element; k++)
			{
				double temp_para[MAX_DIMENSION], temp_coord[MAX_DIMENSION];
				info->octree_subcell[i][j].get_physical_coord(temp_para, &point_array[k * info->DIMENSION], temp_coord, info);

				for (int l = 0; l < info->DIMENSION; l++)
					out_octree_coord[i].emplace_back(temp_coord[l]);
			}

	char str_octree[256] = "octree_subcell.xmf2";
	FILE *fp = fopen(str_octree, "w");

	fprintf(fp, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp, "  <Domain>\n");
	fprintf(fp, "    <Grid Name=\"ien\">\n");
	if (info->DIMENSION == 2)
		fprintf(fp, "      <Topology TopologyType=\"Quadrilateral_9\" NumberOfElements=\"%d\">\n", total_subcell);
	else if (info->DIMENSION == 3)
		fprintf(fp, "      <Topology TopologyType=\"HEXAHEDRON_27\" NumberOfElements=\"%d\">\n", total_subcell);
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n");
	fprintf(fp, "      </Topology>\n");
	if (info->DIMENSION == 2)
		fprintf(fp, "      <Geometry GeometryType=\"XY\">\n");
	else if (info->DIMENSION == 3)
		fprintf(fp, "      <Geometry GeometryType=\"XYZ\">\n");
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n");
	fprintf(fp, "      </Geometry>\n");
	fprintf(fp, "      <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", total_subcell * point_on_element, info->DIMENSION);
	fprintf(fp, "        bin/octree_subcell_xyz.bin\n");
	fprintf(fp, "      </DataItem>\n");
	fstream file_octree_subcell_xyz("./bin/octree_subcell_xyz.bin", ios::binary | ios::out);
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		for (size_t j = 0; j < info->octree_subcell[i].size(); j++)
			for (int k = 0; k < point_on_element; k++)
				for (int l = 0; l < info->DIMENSION; l++)
				{
					double cast_temp = out_octree_coord[i][(j * point_on_element + k) * info->DIMENSION + l];
					file_octree_subcell_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
	file_octree_subcell_xyz.close();
	fprintf(fp, "      <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", total_subcell, point_on_element);
	fprintf(fp, "        bin/octree_subcell_ien.bin\n");
	fprintf(fp, "      </DataItem>\n");
	fstream file_octree_subcell_ien("./bin/octree_subcell_ien.bin", ios::binary | ios::out);
	int counter_ien = 0;
	for (int i = 0; i < total_subcell; i++)
		for (int j = 0; j < point_on_element; j++)
		{
			int cast_temp = counter_ien++;
			file_octree_subcell_ien.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
		}
	file_octree_subcell_ien.close();
	fprintf(fp, "    </Grid>\n");
	fprintf(fp, "  </Domain>\n");
	fprintf(fp, "</Xdmf>\n");
}


// paraview
void output_for_paraview_timestep(information *info, double factor)
{
	Init_viewer_info(info);
	Make_info_for_viewer(info);
	Make_boundary_line(info);

	static int point_on_element = 0;
	static int point_on_edge = 3;
	static int digit = 0;
	static char geometry_type[256];
	static char element_type[256];

	char str_glo[256] = "global_patch.xmf2";
	char str_glo_cp[256] = "global_patch_control_point.xmf2";
	char str_glo_bl[256] = "global_patch_boundary_line.xmf2";
	static FILE *fp_glo = fopen(str_glo, "w");
	static FILE *fp_glo_cp = fopen(str_glo_cp, "w");
	static FILE *fp_glo_bl = fopen(str_glo_bl, "w");
	char str_loc[256] = "local_patch.xmf2";
	char str_loc_cp[256] = "local_patch_control_point.xmf2";
	char str_loc_bl[256] = "local_patch_boundary_line.xmf2";
	static FILE *fp_loc = fopen(str_loc, "w");
	static FILE *fp_loc_cp = fopen(str_loc_cp, "w");
	static FILE *fp_loc_bl = fopen(str_loc_bl, "w");

	// set digit, point_on_element, geometry_type, element_type
	static int count = 0;
	if (count++ == 0)
	{
		int num = info->Total_Control_Point_to_mesh[Total_mesh];
		while (num != 0)
		{
			num = num / 10;
			digit++;
		}

		if (info->DIMENSION == 2)
		{
			point_on_element = 9;
			geometry_type[0] = 'X'; geometry_type[1] = 'Y'; geometry_type[2] = '\0';
			element_type[0] = 'Q'; element_type[1] = 'u'; element_type[2] = 'a'; element_type[3] = 'd'; element_type[4] = 'r'; element_type[5] = 'i'; element_type[6] = 'l'; element_type[7] = 'a'; element_type[8] = 't'; element_type[9] = 'e'; element_type[10] = 'r'; element_type[11] = 'a';  element_type[12] = 'l'; element_type[13] = '_'; element_type[14] = '9'; element_type[15] = '\0';
		}
		else if (info->DIMENSION == 3)
		{
			point_on_element = 27;
			geometry_type[0] = 'X'; geometry_type[1] = 'Y'; geometry_type[2] = 'Z'; geometry_type[3] = '\0'; 
			element_type[0] = 'H'; element_type[1] = 'E'; element_type[2] = 'X'; element_type[3] = 'A'; element_type[4] = 'H'; element_type[5] = 'E'; element_type[6] = 'D'; element_type[7] = 'R'; element_type[8] = 'O'; element_type[9] = 'N'; element_type[10] = '_'; element_type[11] = '2';  element_type[12] = '7'; element_type[13] = '\0';
		}
	}

	// set time
	double time = info->c.TOTAL_SECONDS * factor;
	if (fabs(time) < MERGE_ERROR)
		time = 0.0;
	else if (fabs(time - info->c.TOTAL_SECONDS) < MERGE_ERROR)
		time = info->c.TOTAL_SECONDS;
	
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

		fprintf(fp_loc_cp, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_loc_cp, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_loc_cp, "  <Domain>\n");
		fprintf(fp_loc_cp, "    <Grid Name=\"local_patch_control_point.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");

		fprintf(fp_loc_bl, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp_loc_bl, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp_loc_bl, "  <Domain>\n");
		fprintf(fp_loc_bl, "    <Grid Name=\"local_patch_boundary_line.xmf2\" CollectionType=\"Temporal\" GridType=\"Collection\">\n");
	}

	// output text
	if (info->c.BIN_MODE == 0)
	{
		// middle step global patch
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
			fprintf(fp_glo, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo, info->DIMENSION);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				for (int j = 0; j < info->DIMENSION; j++)
					fprintf(fp_glo, "%.16e ", info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j]);
				fprintf(fp_glo, "\n");
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\">\n", info->Total_Element_to_mesh[1], point_on_element);
			for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
			{
				fprintf(fp_glo, "            ");
				for (int j = 0; j < point_on_element; j++)
					fprintf(fp_glo, "%*d ", -digit, info->Connectivity[i * point_on_element + j]);
				fprintf(fp_glo, "\n");
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo, info->DIMENSION);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				for (int j = 0; j < info->DIMENSION; j++)
					fprintf(fp_glo, "%.16e ", info->disp_at_connectivity[i * info->DIMENSION + j]);
				fprintf(fp_glo, "\n");
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			if (info->DIMENSION == 2)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 3]);
				}
			else if (info->DIMENSION == 3)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 2]);
				}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			if (info->DIMENSION == 2)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->elastic_strain_at_connectivity[i * N_STRAIN + 0], info->elastic_strain_at_connectivity[i * N_STRAIN + 1], info->elastic_strain_at_connectivity[i * N_STRAIN + 3]);
				}
			else if (info->DIMENSION == 3)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->elastic_strain_at_connectivity[i * N_STRAIN + 0], info->elastic_strain_at_connectivity[i * N_STRAIN + 1], info->elastic_strain_at_connectivity[i * N_STRAIN + 2]);
				}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				fprintf(fp_glo, "%.16e\n", info->vm_at_connectivity[i]);
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				fprintf(fp_glo, "%.16e\n", info->hs_at_connectivity[i]);
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				fprintf(fp_glo, "%.16e\n", info->equivalent_plastic_strain_at_connectivity[i]);
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				fprintf(fp_glo, "            ");
				fprintf(fp_glo, "%.16e\n", info->yield_stress_at_connectivity[i]);
			}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity_glo);
			if (info->DIMENSION == 2)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->back_stress_at_connectivity[i * N_STRESS + 0], info->back_stress_at_connectivity[i * N_STRESS + 1], info->back_stress_at_connectivity[i * N_STRESS + 3]);
				}
			else if (info->DIMENSION == 3)
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					fprintf(fp_glo, "            ");
					fprintf(fp_glo, "%.16e %.16e %.16e\n", info->back_stress_at_connectivity[i * N_STRESS + 0], info->back_stress_at_connectivity[i * N_STRESS + 1], info->back_stress_at_connectivity[i * N_STRESS + 2]);
				}
			fprintf(fp_glo, "          </DataItem>\n");
			fprintf(fp_glo, "        </Grid>\n");
			fprintf(fp_glo, "      </Grid>\n");
		}

		// middle step global control point
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
			fprintf(fp_glo_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", info->Total_Control_Point_to_mesh[1], info->DIMENSION);
			if (Total_mesh < 2)
				for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
				{
					fprintf(fp_glo_cp, "            ");
					for (int j = 0; j < info->DIMENSION; j++)
						fprintf(fp_glo_cp, "%.16e ", info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + info->disp[i * info->DIMENSION + j]);
					fprintf(fp_glo_cp, "\n");
				}
			else
				for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
				{
					fprintf(fp_glo_cp, "            ");
					for (int j = 0; j < info->DIMENSION; j++)
						fprintf(fp_glo_cp, "%.16e ", info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + info->disp_overlay[i * info->DIMENSION + j]);
					fprintf(fp_glo_cp, "\n");
				}
			fprintf(fp_glo_cp, "          </DataItem>\n");
			fprintf(fp_glo_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" NumberType=\"Int\" Format=\"XML\">\n", info->Total_Control_Point_to_mesh[1]);
			for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
				fprintf(fp_glo_cp, "            %d\n", i);
			fprintf(fp_glo_cp, "          </DataItem>\n");
			fprintf(fp_glo_cp, "        </Grid>\n");
			fprintf(fp_glo_cp, "      </Grid>\n");
		}

		// middle step global patch boundary line
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
			fprintf(fp_glo_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_edge_glo * point_on_edge, info->DIMENSION);
			for (int i = 0; i < Total_edge_glo; i++)
			{
				for (int j = 0; j < point_on_edge; j++)
				{
					fprintf(fp_glo_bl, "            ");
					for (int k = 0; k < info->DIMENSION; k++)
						fprintf(fp_glo_bl, "%.16e ", info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k]);
					fprintf(fp_glo_bl, "\n");
				}
			}
			fprintf(fp_glo_bl, "          </DataItem>\n");
			fprintf(fp_glo_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\">\n", Total_edge_glo, point_on_edge);
			for (int i = 0; i < Total_edge_glo; i++)
				fprintf(fp_glo_bl, "            %*d %*d %*d\n", -digit, i * 3, -digit, i * 3 + 1, -digit, i * 3 + 2);
			fprintf(fp_glo_bl, "          </DataItem>\n");
			fprintf(fp_glo_bl, "        </Grid>\n");
			fprintf(fp_glo_bl, "      </Grid>\n");
		}

		if (Total_mesh > 1)
		{
			// middle step local patch
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
				fprintf(fp_loc, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					for (int j = 0; j < info->DIMENSION; j++)
						fprintf(fp_loc, "%.16e ", info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j]);
					fprintf(fp_loc, "\n");
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\">\n", info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1], point_on_element);
				for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
				{
					fprintf(fp_loc, "            ");
					for (int j = 0; j < point_on_element; j++)
						fprintf(fp_loc, "%*d ", -digit, info->Connectivity[i * point_on_element + j] - Total_connectivity_glo);
					fprintf(fp_loc, "\n");
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					for (int j = 0; j < info->DIMENSION; j++)
						fprintf(fp_loc, "%.16e ", info->disp_at_connectivity[i * info->DIMENSION + j]);
					fprintf(fp_loc, "\n");
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				if (info->DIMENSION == 2)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 3]);
					}
				else if (info->DIMENSION == 3)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 2]);
					}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				if (info->DIMENSION == 2)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->elastic_strain_at_connectivity[i * N_STRAIN + 0], info->elastic_strain_at_connectivity[i * N_STRAIN + 1], info->elastic_strain_at_connectivity[i * N_STRAIN + 3]);
					}
				else if (info->DIMENSION == 3)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->elastic_strain_at_connectivity[i * N_STRAIN + 0], info->elastic_strain_at_connectivity[i * N_STRAIN + 1], info->elastic_strain_at_connectivity[i * N_STRAIN + 2]);
					}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					fprintf(fp_loc, "%.16e\n", info->vm_at_connectivity[i]);
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					fprintf(fp_loc, "%.16e\n", info->hs_at_connectivity[i]);
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					fprintf(fp_loc, "%.16e\n", info->equivalent_plastic_strain_at_connectivity[i]);
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					fprintf(fp_loc, "            ");
					fprintf(fp_loc, "%.16e\n", info->yield_stress_at_connectivity[i]);
				}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", Total_connectivity - Total_connectivity_glo);
				if (info->DIMENSION == 2)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->back_stress_at_connectivity[i * N_STRESS + 0], info->back_stress_at_connectivity[i * N_STRESS + 1], info->back_stress_at_connectivity[i * N_STRESS + 3]);
					}
				else if (info->DIMENSION == 3)
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						fprintf(fp_loc, "            ");
						fprintf(fp_loc, "%.16e %.16e %.16e\n", info->back_stress_at_connectivity[i * N_STRESS + 0], info->back_stress_at_connectivity[i * N_STRESS + 1], info->back_stress_at_connectivity[i * N_STRESS + 2]);
					}
				fprintf(fp_loc, "          </DataItem>\n");
				fprintf(fp_loc, "        </Grid>\n");
				fprintf(fp_loc, "      </Grid>\n");
			}

			// middle step local control point
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
				fprintf(fp_loc_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1], info->DIMENSION);
				if (Total_mesh < 2)
					for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
					{
						fprintf(fp_loc_cp, "            ");
						for (int j = 0; j < info->DIMENSION; j++)
							fprintf(fp_loc_cp, "%.16e ", info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + info->disp[i * info->DIMENSION + j]);
						fprintf(fp_loc_cp, "\n");
					}
				else
					for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
					{
						fprintf(fp_loc_cp, "            ");
						for (int j = 0; j < info->DIMENSION; j++)
							fprintf(fp_loc_cp, "%.16e ", info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + info->disp_overlay[i * info->DIMENSION + j]);
						fprintf(fp_loc_cp, "\n");
					}
				fprintf(fp_loc_cp, "          </DataItem>\n");
				fprintf(fp_loc_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" NumberType=\"Int\" Format=\"XML\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]);
				for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]; i++)
					fprintf(fp_loc_cp, "            %d\n", i);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				fprintf(fp_loc_cp, "        </Grid>\n");
				fprintf(fp_loc_cp, "      </Grid>\n");
			}

			// middle step local patch boundary line
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
				fprintf(fp_loc_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", (Total_edge - Total_edge_glo) * point_on_edge, info->DIMENSION);
				for (int i = Total_edge_glo; i < Total_edge; i++)
				{
					for (int j = 0; j < point_on_edge; j++)
					{
						fprintf(fp_loc_bl, "            ");
						for (int k = 0; k < info->DIMENSION; k++)
							fprintf(fp_loc_bl, "%.16e ", info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k]);
						fprintf(fp_loc_bl, "\n");
					}
				}
				fprintf(fp_loc_bl, "          </DataItem>\n");
				fprintf(fp_loc_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\">\n", Total_edge - Total_edge_glo, point_on_edge);
				for (int i = 0; i < Total_edge - Total_edge_glo; i++)
					fprintf(fp_loc_bl, "            %*d %*d %*d\n", -digit, i * 3, -digit, i * 3 + 1, -digit, i * 3 + 2);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				fprintf(fp_loc_bl, "        </Grid>\n");
				fprintf(fp_loc_bl, "      </Grid>\n");
			}
		}
	}

	// output binary
	else if (info->c.BIN_MODE == 1)
	{
		// middle step global patch
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
			fprintf(fp_glo, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            bin/global_patch_xyz%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_xyz = "./bin/global_patch_xyz" + to_string(step) + ".bin";
			fstream file_glo_xyz(s_file_glo_xyz, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
				for (int j = 0; j < info->DIMENSION; j++)
				{
					double cast_temp = (info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j]);
					file_glo_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			file_glo_xyz.close();
			fprintf(fp_glo, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Element_to_mesh[1], point_on_element);
			fprintf(fp_glo, "            bin/global_patch_ien%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_ien = "./bin/global_patch_ien" + to_string(step) + ".bin";
			fstream file_glo_ien(s_file_glo_ien, ios::binary | ios::out);
			for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
				for (int j = 0; j < point_on_element; j++)
				{
					int cast_temp = info->Connectivity[i * point_on_element + j];
					file_glo_ien.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			file_glo_ien.close();
			fprintf(fp_glo, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo, info->DIMENSION);
			fprintf(fp_glo, "            bin/global_patch_disp%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_disp = "./bin/global_patch_disp" + to_string(step) + ".bin";
			fstream file_glo_disp(s_file_glo_disp, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
				for (int j = 0; j < info->DIMENSION; j++)
				{
					double cast_temp = info->disp_at_connectivity[i * info->DIMENSION + j];
					file_glo_disp.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			file_glo_disp.close();
			fprintf(fp_glo, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_stress%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_stress = "./bin/global_patch_stress" + to_string(step) + ".bin";
			fstream file_glo_stress(s_file_glo_stress, ios::binary | ios::out);
			if (info->DIMENSION == 2)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->stress_at_connectivity[i * N_STRESS];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->stress_at_connectivity[i * N_STRESS + 1];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->stress_at_connectivity[i * N_STRESS + 3];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			else if (info->DIMENSION == 3)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->stress_at_connectivity[i * N_STRESS];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->stress_at_connectivity[i * N_STRESS + 1];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->stress_at_connectivity[i * N_STRESS + 2];
					file_glo_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			file_glo_stress.close();
			fprintf(fp_glo, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_elastic_strain%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_elastic_strain = "./bin/global_patch_elastic_strain" + to_string(step) + ".bin";
			fstream file_glo_elastic_strain(s_file_glo_elastic_strain, ios::binary | ios::out);
			if (info->DIMENSION == 2)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 3];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			else if (info->DIMENSION == 3)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 2];
					file_glo_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			file_glo_elastic_strain.close();
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_eqs%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_eqs = "./bin/global_patch_eqs" + to_string(step) + ".bin";
			fstream file_glo_eqs(s_file_glo_eqs, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				double cast_temp = info->vm_at_connectivity[i];
				file_glo_eqs.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
			}
			fprintf(fp_glo, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_hs%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_hs = "./bin/global_patch_hs" + to_string(step) + ".bin";
			fstream file_glo_hs(s_file_glo_hs, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				double cast_temp = info->hs_at_connectivity[i];
				file_glo_hs.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
			}
			file_glo_eqs.close();
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_eps%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_eps = "./bin/global_patch_eps" + to_string(step) + ".bin";
			fstream file_glo_eps(s_file_glo_eps, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				double cast_temp = info->equivalent_plastic_strain_at_connectivity[i];
				file_glo_eps.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
			}
			file_glo_eps.close();
			fprintf(fp_glo, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_ys%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_ys = "./bin/global_patch_ys" + to_string(step) + ".bin";
			fstream file_glo_ys(s_file_glo_ys, ios::binary | ios::out);
			for (int i = 0; i < Total_connectivity_glo; i++)
			{
				double cast_temp = info->yield_stress_at_connectivity[i];
				file_glo_ys.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
			}
			file_glo_ys.close();
			fprintf(fp_glo, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            bin/global_patch_bs%d.bin\n", step);
			fprintf(fp_glo, "          </DataItem>\n");
			string s_file_glo_back_stress = "./bin/global_patch_bs" + to_string(step) + ".bin";
			fstream file_glo_back_stress(s_file_glo_back_stress, ios::binary | ios::out);
			if (info->DIMENSION == 2)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 1];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 3];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			else if (info->DIMENSION == 3)
			{
				double cast_temp;
				for (int i = 0; i < Total_connectivity_glo; i++)
				{
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 1];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 2];
					file_glo_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			}
			file_glo_back_stress.close();
			fprintf(fp_glo, "        </Grid>\n");
			fprintf(fp_glo, "      </Grid>\n");
		}

		// middle step global control point
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
			fprintf(fp_glo_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Control_Point_to_mesh[1], info->DIMENSION);
			fprintf(fp_glo_cp, "            bin/global_patch_cp_xyz%d.bin\n", step);
			fprintf(fp_glo_cp, "          </DataItem>\n");
			string s_file_glo_cp_xyz = "./bin/global_patch_cp_xyz" + to_string(step) + ".bin";
			fstream file_glo_cp_xyz(s_file_glo_cp_xyz, ios::binary | ios::out);
			if (Total_mesh < 2)
				for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
					for (int j = 0; j < info->DIMENSION; j++)
					{
						double cast_temp = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp[i * info->DIMENSION + j]);
						file_glo_cp_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
			else
				for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
					for (int j = 0; j < info->DIMENSION; j++)
					{
						double cast_temp = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp_overlay[i * info->DIMENSION + j]);
						file_glo_cp_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
			file_glo_cp_xyz.close();
			fprintf(fp_glo_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Control_Point_to_mesh[1]);
			fprintf(fp_glo_cp, "            bin/global_patch_cp_ien%d.bin\n", step);
			fprintf(fp_glo_cp, "          </DataItem>\n");
			string s_file_glo_cp_ien = "./bin/global_patch_cp_ien" + to_string(step) + ".bin";
			fstream file_glo_cp_ien(s_file_glo_cp_ien, ios::binary | ios::out);
			for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
				file_glo_cp_ien.write((reinterpret_cast<char *>(&i)), sizeof(i));
			file_glo_cp_ien.close();
			fprintf(fp_glo_cp, "        </Grid>\n");
			fprintf(fp_glo_cp, "      </Grid>\n");
		}

		// middle step global patch boundary line
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
			fprintf(fp_glo_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_edge_glo * point_on_edge, info->DIMENSION);
			fprintf(fp_glo_bl, "            bin/global_patch_boundary_xyz%d.bin\n", step);
			fprintf(fp_glo_bl, "          </DataItem>\n");
			string s_file_glo_boundary_xyz = "./bin/global_patch_boundary_xyz" + to_string(step) + ".bin";
			fstream file_glo_boundary_xyz(s_file_glo_boundary_xyz, ios::binary | ios::out);
			for (int i = 0; i < Total_edge_glo; i++)
				for (int j = 0; j < point_on_edge; j++)
					for (int k = 0; k < info->DIMENSION; k++)
					{
						double cast_temp = info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k];
						file_glo_boundary_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
			file_glo_boundary_xyz.close();
			fprintf(fp_glo_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", Total_edge_glo, point_on_edge);
			fprintf(fp_glo_bl, "            bin/global_patch_boundary_ien%d.bin\n", step);
			fprintf(fp_glo_bl, "          </DataItem>\n");
			string s_file_glo_boundary_ien = "./bin/global_patch_boundary_ien" + to_string(step) + ".bin";
			fstream file_glo_boundary_ien(s_file_glo_boundary_ien, ios::binary | ios::out);
			for (int i = 0; i < Total_edge_glo; i++)
				for (int j = 0; j < 3; j++)
				{
					int cast_temp = i * 3 + j;
					file_glo_boundary_ien.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
			file_glo_boundary_ien.close();
			fprintf(fp_glo_bl, "        </Grid>\n");
			fprintf(fp_glo_bl, "      </Grid>\n");
		}

		if (Total_mesh > 1)
		{
			// middle step local patch
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
				fprintf(fp_loc, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            bin/local_patch_xyz%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_xyz = "./bin/local_patch_xyz" + to_string(step) + ".bin";
				fstream file_loc_xyz(s_file_loc_xyz, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					for (int j = 0; j < info->DIMENSION; j++)
					{
						double cast_temp = (info->Connectivity_coord[i * info->DIMENSION + j] + info->disp_at_connectivity[i * info->DIMENSION + j]);
						file_loc_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				file_loc_xyz.close();
				fprintf(fp_loc, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1], point_on_element);
				fprintf(fp_loc, "            bin/local_patch_ien%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_ien = "./bin/local_patch_ien" + to_string(step) + ".bin";
				fstream file_loc_ien(s_file_loc_ien, ios::binary | ios::out);
				for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
					for (int j = 0; j < point_on_element; j++)
					{
						int cast_temp = info->Connectivity[i * point_on_element + j] - Total_connectivity_glo;
						file_loc_ien.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				file_loc_ien.close();
				fprintf(fp_loc, "          <DataItem Name=\"displacement\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
				fprintf(fp_loc, "            bin/local_patch_disp%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_disp = "./bin/local_patch_disp" + to_string(step) + ".bin";
				fstream file_loc_disp(s_file_loc_disp, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					for (int j = 0; j < info->DIMENSION; j++)
					{
						double cast_temp = info->disp_at_connectivity[i * info->DIMENSION + j];
						file_loc_disp.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				file_loc_disp.close();
				fprintf(fp_loc, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_stress%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_stress = "./bin/local_patch_stress" + to_string(step) + ".bin";
				fstream file_loc_stress(s_file_loc_stress, ios::binary | ios::out);
				if (info->DIMENSION == 2)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->stress_at_connectivity[i * N_STRESS];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->stress_at_connectivity[i * N_STRESS + 1];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->stress_at_connectivity[i * N_STRESS + 3];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				else if (info->DIMENSION == 3)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->stress_at_connectivity[i * N_STRESS];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->stress_at_connectivity[i * N_STRESS + 1];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->stress_at_connectivity[i * N_STRESS + 2];
						file_loc_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				file_loc_stress.close();
				fprintf(fp_loc, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_elastic_strain%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_elastic_strain = "./bin/local_patch_elastic_strain" + to_string(step) + ".bin";
				fstream file_loc_elastic_strain(s_file_loc_elastic_strain, ios::binary | ios::out);
				if (info->DIMENSION == 2)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 3];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				else if (info->DIMENSION == 3)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->elastic_strain_at_connectivity[i * N_STRAIN + 2];
						file_loc_elastic_strain.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				file_loc_elastic_strain.close();
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_eqs%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_eqs = "./bin/local_patch_eqs" + to_string(step) + ".bin";
				fstream file_loc_eqs(s_file_loc_eqs, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					double cast_temp = info->vm_at_connectivity[i];
					file_loc_eqs.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
				file_loc_eqs.close();
				fprintf(fp_loc, "          <DataItem Name=\"hydrostatic_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_hs%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_hs = "./bin/local_patch_hs" + to_string(step) + ".bin";
				fstream file_loc_hs(s_file_loc_hs, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					double cast_temp = info->hs_at_connectivity[i];
					file_loc_hs.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
				file_loc_hs.close();
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_eps%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_eps = "./bin/local_patch_eps" + to_string(step) + ".bin";
				fstream file_loc_eps(s_file_loc_eps, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					double cast_temp = info->equivalent_plastic_strain_at_connectivity[i];
					file_loc_eps.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
				file_loc_eps.close();
				fprintf(fp_loc, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_ys%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_ys = "./bin/local_patch_ys" + to_string(step) + ".bin";
				fstream file_loc_ys(s_file_loc_ys, ios::binary | ios::out);
				for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
				{
					double cast_temp = info->yield_stress_at_connectivity[i];
					file_loc_ys.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
				}
				file_loc_ys.close();
				fprintf(fp_loc, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            bin/local_patch_bs%d.bin\n", step);
				fprintf(fp_loc, "          </DataItem>\n");
				string s_file_loc_back_stress = "./bin/local_patch_bs" + to_string(step) + ".bin";
				fstream file_loc_back_stress(s_file_loc_back_stress, ios::binary | ios::out);
				if (info->DIMENSION == 2)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 1];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 3];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				else if (info->DIMENSION == 3)
				{
					double cast_temp;
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
					{
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 1];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						cast_temp = info->back_stress_at_connectivity[i * N_STRESS + 2];
						file_loc_back_stress.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				}
				file_loc_back_stress.close();
				fprintf(fp_loc, "        </Grid>\n");
				fprintf(fp_loc, "      </Grid>\n");
			}

			// middle step local control point
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
				fprintf(fp_loc_cp, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1], info->DIMENSION);
				fprintf(fp_loc_cp, "            bin/local_patch_cp_xyz%d.bin\n", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				string s_file_loc_cp_xyz = "./bin/local_patch_cp_xyz" + to_string(step) + ".bin";
				fstream file_loc_cp_xyz(s_file_loc_cp_xyz, ios::binary | ios::out);
				if (Total_mesh < 2)
					for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
						for (int j = 0; j < info->DIMENSION; j++)
						{
							double cast_temp = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp[i * info->DIMENSION + j]);
							file_loc_cp_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						}
				else
					for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
						for (int j = 0; j < info->DIMENSION; j++)
						{
							double cast_temp = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp_overlay[i * info->DIMENSION + j]);
							file_loc_cp_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						}
				file_loc_cp_xyz.close();
				fprintf(fp_loc_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]);
				fprintf(fp_loc_cp, "            bin/local_patch_cp_ien%d.bin\n", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				string s_file_loc_cp_ien = "./bin/local_patch_cp_ien" + to_string(step) + ".bin";
				fstream file_loc_cp_ien(s_file_loc_cp_ien, ios::binary | ios::out);
				for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]; i++)
					file_loc_cp_ien.write((reinterpret_cast<char *>(&i)), sizeof(i));
				file_loc_cp_ien.close();
				fprintf(fp_loc_cp, "        </Grid>\n");
				fprintf(fp_loc_cp, "      </Grid>\n");
			}

			// middle step local patch boundary line
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
				fprintf(fp_loc_bl, "          <DataItem Name=\"xyz\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">\n", (Total_edge - Total_edge_glo) * point_on_edge, info->DIMENSION);
				fprintf(fp_loc_bl, "            bin/local_patch_boundary_xyz%d.bin\n", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				string s_file_loc_boundary_xyz = "./bin/local_patch_boundary_xyz" + to_string(step) + ".bin";
				fstream file_loc_boundary_xyz(s_file_loc_boundary_xyz, ios::binary | ios::out);
				for (int i = Total_edge_glo; i < Total_edge; i++)
					for (int j = 0; j < point_on_edge; j++)
						for (int k = 0; k < info->DIMENSION; k++)
						{
							double cast_temp = info->Edge_coord[i * point_on_edge * info->DIMENSION + j * info->DIMENSION + k];
							file_loc_boundary_xyz.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
						}
				file_loc_boundary_xyz.close();
				fprintf(fp_loc_bl, "          <DataItem Name=\"ien\" Dimensions=\"%d %d\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\">\n", Total_edge - Total_edge_glo, point_on_edge);
				fprintf(fp_loc_bl, "            bin/local_patch_boundary_ien%d.bin\n", step);
				fprintf(fp_loc_bl, "          </DataItem>\n");
				string s_file_loc_boundary_ien = "./bin/local_patch_boundary_ien" + to_string(step) + ".bin";
				fstream file_loc_boundary_ien(s_file_loc_boundary_ien, ios::binary | ios::out);
				for (int i = 0; i < Total_edge - Total_edge_glo; i++)
					for (int j = 0; j < 3; j++)
					{
						int cast_temp = i * 3 + j;
						file_loc_boundary_ien.write((reinterpret_cast<char *>(&cast_temp)), sizeof(cast_temp));
					}
				file_loc_boundary_ien.close();
				fprintf(fp_loc_bl, "        </Grid>\n");
				fprintf(fp_loc_bl, "      </Grid>\n");
			}
		}
	}

	// output HDF5
	else if (info->c.BIN_MODE == 2)
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

		// middle step global patch
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
			fprintf(fp_glo, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
			fprintf(fp_glo, "          </Attribute>\n");
			fprintf(fp_glo, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
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
			fprintf(fp_glo, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/stress\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_stress[2] = { (hsize_t)Total_connectivity_glo, 3 };
				DataSpace dataspace_stress(2, dims_stress);
				DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
				vector<double> buf_stress((size_t)Total_connectivity_glo * 3);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 3 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 3 + 1] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 3 + 2] = info->stress_at_connectivity[i * N_STRESS + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_stress[i * 3 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
						buf_stress[i * 3 + 1] = info->stress_at_connectivity[i * N_STRESS + 1];
						buf_stress[i * 3 + 2] = info->stress_at_connectivity[i * N_STRESS + 2];
					}
				}
				dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/elastic_strain\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_es[2] = { (hsize_t)Total_connectivity_glo, 3 };
				DataSpace dataspace_es(2, dims_es);
				DSetCreatPropList plist_es = make_plist(dims_es, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_es = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_es, plist_es);
				vector<double> buf_es((size_t)Total_connectivity_glo * 3);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 3 + 0] = info->elastic_strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 3 + 1] = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 3 + 2] = info->elastic_strain_at_connectivity[i * N_STRAIN + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_es[i * 3 + 0] = info->elastic_strain_at_connectivity[i * N_STRAIN + 0];
						buf_es[i * 3 + 1] = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
						buf_es[i * 3 + 2] = info->elastic_strain_at_connectivity[i * N_STRAIN + 2];
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
			fprintf(fp_glo, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/eps\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_eps[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_eps(1, dims_eps);
				DSetCreatPropList plist_eps = make_plist(dims_eps, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_eps = file.createDataSet("/eps", PredType::NATIVE_DOUBLE, dataspace_eps, plist_eps);
				vector<double> buf_eps((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_eps[i] = info->equivalent_plastic_strain_at_connectivity[i];
				dset_eps.write(buf_eps.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/ys\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_ys[1] = { (hsize_t)Total_connectivity_glo };
				DataSpace dataspace_ys(1, dims_ys);
				DSetCreatPropList plist_ys = make_plist(dims_ys, 1, PredType::NATIVE_DOUBLE);
				DataSet dset_ys = file.createDataSet("/ys", PredType::NATIVE_DOUBLE, dataspace_ys, plist_ys);
				vector<double> buf_ys((size_t)Total_connectivity_glo);
				for (int i = 0; i < Total_connectivity_glo; i++)
					buf_ys[i] = info->yield_stress_at_connectivity[i];
				dset_ys.write(buf_ys.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity_glo);
			fprintf(fp_glo, "            %s/global_patch_step%d.h5:/back_stress\n", "./bin", step);
			fprintf(fp_glo, "          </DataItem>\n");
			{
				string h5name = string("./bin/global_patch_step") + to_string(step) + string(".h5");
				H5File file(h5name, H5F_ACC_RDWR);
				hsize_t dims_bs[2] = { (hsize_t)Total_connectivity_glo, 3 };
				DataSpace dataspace_bs(2, dims_bs);
				DSetCreatPropList plist_bs = make_plist(dims_bs, 2, PredType::NATIVE_DOUBLE);
				DataSet dset_bs = file.createDataSet("/back_stress", PredType::NATIVE_DOUBLE, dataspace_bs, plist_bs);
				vector<double> buf_bs((size_t)Total_connectivity_glo * 3);
				if (info->DIMENSION == 2)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_bs[i * 3 + 0] = info->back_stress_at_connectivity[i * N_STRESS + 0];
						buf_bs[i * 3 + 1] = info->back_stress_at_connectivity[i * N_STRESS + 1];
						buf_bs[i * 3 + 2] = info->back_stress_at_connectivity[i * N_STRESS + 3];
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (int i = 0; i < Total_connectivity_glo; i++)
					{
						buf_bs[i * 3 + 0] = info->back_stress_at_connectivity[i * N_STRESS + 0];
						buf_bs[i * 3 + 1] = info->back_stress_at_connectivity[i * N_STRESS + 1];
						buf_bs[i * 3 + 2] = info->back_stress_at_connectivity[i * N_STRESS + 2];
					}
				}
				dset_bs.write(buf_bs.data(), PredType::NATIVE_DOUBLE);
			}
			fprintf(fp_glo, "        </Grid>\n");
			fprintf(fp_glo, "      </Grid>\n");
		}

		// middle step global control point
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
							buf_cp_xyz[i * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp[i * info->DIMENSION + j]);
				else
					for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
						for (int j = 0; j < info->DIMENSION; j++)
							buf_cp_xyz[i * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp_overlay[i * info->DIMENSION + j]);
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

		// middle step global patch boundary line
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

		// local patch
		if (Total_mesh > 1)
		{
			// middle step local patch
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
				fprintf(fp_loc, "          <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"elastic_strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;elastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"hydrostatic_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;hydrostatic_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"equivalent_plastic_strain\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;equivalent_plastic_strain&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"yield_stress\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;yield_stress&quot;]</DataItem>\n", step + 1);
				fprintf(fp_loc, "          </Attribute>\n");
				fprintf(fp_loc, "          <Attribute Name=\"back_stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/Grid[%d]//Grid/DataItem[@Name=&quot;back_stress&quot;]</DataItem>\n", step + 1);
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
				fprintf(fp_loc, "          <DataItem Name=\"stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_stress[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 3 };
					DataSpace dataspace_stress(2, dims_stress);
					DSetCreatPropList plist_stress = make_plist(dims_stress, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_stress = file.createDataSet("/stress", PredType::NATIVE_DOUBLE, dataspace_stress, plist_stress);
					vector<double> buf_stress((size_t)(Total_connectivity - Total_connectivity_glo) * 3);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 3 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 3 + 1] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 3 + 2] = info->stress_at_connectivity[i * N_STRESS + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_stress[(i - Total_connectivity_glo) * 3 + 0] = info->stress_at_connectivity[i * N_STRESS + 0];
							buf_stress[(i - Total_connectivity_glo) * 3 + 1] = info->stress_at_connectivity[i * N_STRESS + 1];
							buf_stress[(i - Total_connectivity_glo) * 3 + 2] = info->stress_at_connectivity[i * N_STRESS + 2];
						}
					}
					dset_stress.write(buf_stress.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"elastic_strain\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/elastic_strain\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_elastic_strain[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 3 };
					DataSpace dataspace_elastic_strain(2, dims_elastic_strain);
					DSetCreatPropList plist_elastic = make_plist(dims_elastic_strain, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_elastic_strain = file.createDataSet("/elastic_strain", PredType::NATIVE_DOUBLE, dataspace_elastic_strain, plist_elastic);
					vector<double> buf_elastic_strain((size_t)(Total_connectivity - Total_connectivity_glo) * 3);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 0] = info->elastic_strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 1] = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 2] = info->elastic_strain_at_connectivity[i * N_STRAIN + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 0] = info->elastic_strain_at_connectivity[i * N_STRAIN + 0];
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 1] = info->elastic_strain_at_connectivity[i * N_STRAIN + 1];
							buf_elastic_strain[(i - Total_connectivity_glo) * 3 + 2] = info->elastic_strain_at_connectivity[i * N_STRAIN + 2];
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
				fprintf(fp_loc, "          <DataItem Name=\"equivalent_plastic_strain\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/equivalent_plastic_strain\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_eps[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_eps(1, dims_eps);
					DSetCreatPropList plist_eps = make_plist(dims_eps, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_eps = file.createDataSet("/equivalent_plastic_strain", PredType::NATIVE_DOUBLE, dataspace_eps, plist_eps);
					vector<double> buf_eps((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_eps[i - Total_connectivity_glo] = info->equivalent_plastic_strain_at_connectivity[i];
					dset_eps.write(buf_eps.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"yield_stress\" Dimensions=\"%d\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "			%s/local_patch_step%d.h5:/yield_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_ys[1] = { (hsize_t)(Total_connectivity - Total_connectivity_glo) };
					DataSpace dataspace_ys(1, dims_ys);
					DSetCreatPropList plist_ys = make_plist(dims_ys, 1, PredType::NATIVE_DOUBLE);
					DataSet dset_ys = file.createDataSet("/yield_stress", PredType::NATIVE_DOUBLE, dataspace_ys, plist_ys);
					vector<double> buf_ys((size_t)(Total_connectivity - Total_connectivity_glo));
					for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						buf_ys[i - Total_connectivity_glo] = info->yield_stress_at_connectivity[i];
					dset_ys.write(buf_ys.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc, "          <DataItem Name=\"back_stress\" Dimensions=\"%d 3\" ItemType=\"Uniform\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", Total_connectivity - Total_connectivity_glo);
				fprintf(fp_loc, "            %s/local_patch_step%d.h5:/back_stress\n", "./bin", step);
				fprintf(fp_loc, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_patch_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_back_stress[2] = { (hsize_t)(Total_connectivity - Total_connectivity_glo), 3 };
					DataSpace dataspace_back_stress(2, dims_back_stress);
					DSetCreatPropList plist_bs = make_plist(dims_back_stress, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_back_stress = file.createDataSet("/back_stress", PredType::NATIVE_DOUBLE, dataspace_back_stress, plist_bs);
					vector<double> buf_back_stress((size_t)(Total_connectivity - Total_connectivity_glo) * 3);
					if (info->DIMENSION == 2)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 0] = info->back_stress_at_connectivity[i * N_STRESS + 0];
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 1] = info->back_stress_at_connectivity[i * N_STRESS + 1];
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 2] = info->back_stress_at_connectivity[i * N_STRESS + 3];
						}
					}
					else if (info->DIMENSION == 3)
					{
						for (int i = Total_connectivity_glo; i < Total_connectivity; i++)
						{
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 0] = info->back_stress_at_connectivity[i * N_STRESS + 0];
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 1] = info->back_stress_at_connectivity[i * N_STRESS + 1];
							buf_back_stress[(i - Total_connectivity_glo) * 3 + 2] = info->back_stress_at_connectivity[i * N_STRESS + 2];
						}
					}
					dset_back_stress.write(buf_back_stress.data(), PredType::NATIVE_DOUBLE);
				}
	
				fprintf(fp_loc, "        </Grid>\n");
				fprintf(fp_loc, "      </Grid>\n");
			}
	
			// middle step local control point
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
					hsize_t dims_cp_xyz[2] = { (hsize_t)(info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]), (hsize_t)(info->DIMENSION) };
					DataSpace dataspace_cp_xyz(2, dims_cp_xyz);
					DSetCreatPropList plist_cp_xyz = make_plist(dims_cp_xyz, 2, PredType::NATIVE_DOUBLE);
					DataSet dset_cp_xyz = file.createDataSet("/xyz", PredType::NATIVE_DOUBLE, dataspace_cp_xyz, plist_cp_xyz);
					vector<double> buf_cp_xyz((size_t)(info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION);
					if (Total_mesh < 2)
						for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[(i - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp[i * info->DIMENSION + j]);
					else
						for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
							for (int j = 0; j < info->DIMENSION; j++)
								buf_cp_xyz[(i - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j] = info->Node_Coordinate[i * (info->DIMENSION + 1) + j] + (info->disp_overlay[i * info->DIMENSION + j]);
					dset_cp_xyz.write(buf_cp_xyz.data(), PredType::NATIVE_DOUBLE);
				}
				fprintf(fp_loc_cp, "          <DataItem Name=\"ien\" Dimensions=\"%d 1\" ItemType=\"Uniform\" NumberType=\"Int\" Format=\"HDF\">\n", info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]);
				fprintf(fp_loc_cp, "            %s/local_patch_cp_step%d.h5:/ien\n", "./bin", step);
				fprintf(fp_loc_cp, "          </DataItem>\n");
				{
					string h5name = string("./bin/local_cp_step") + to_string(step) + string(".h5");
					H5File file(h5name, H5F_ACC_RDWR);
					hsize_t dims_cp_ien[2] = { (hsize_t)(info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]), 1 };
					DataSpace dataspace_cp_ien(2, dims_cp_ien);
					DSetCreatPropList plist_cp_ien = make_plist(dims_cp_ien, 2, PredType::NATIVE_INT);
					DataSet dset_cp_ien = file.createDataSet("/ien", PredType::NATIVE_INT, dataspace_cp_ien, plist_cp_ien);
					vector<int> buf_cp_ien((size_t)(info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]));
					for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1]; i++)
						buf_cp_ien[i] = i;
					dset_cp_ien.write(buf_cp_ien.data(), PredType::NATIVE_INT);
				}
	
				fprintf(fp_loc_cp, "        </Grid>\n");
				fprintf(fp_loc_cp, "      </Grid>\n");
			}
	
			// middle step local patch boundary line
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
		}
	}

	// final step
	if (time == info->c.TOTAL_SECONDS)
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

		fprintf(fp_loc_cp, "    </Grid>\n");
		fprintf(fp_loc_cp, "  </Domain>\n");
		fprintf(fp_loc_cp, "</Xdmf>\n");

		fprintf(fp_loc_bl, "    </Grid>\n");
		fprintf(fp_loc_bl, "  </Domain>\n");
		fprintf(fp_loc_bl, "</Xdmf>\n");
	}

	step++;

	printf("Finish 'output_for_paraview_timestep'\n\n");
}


#if 0
void output_quantity_timestep(information *info, int step)
{
	// disp
	// elastic_strain
	// stress
	// back_stress
	// von_mises_stress
	// equivalent_plastic_strain
	// yield_stress 
}
#endif


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
				if (info->c.PARAVIEW_CRACK_REPRESENTATION)
				{
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
				}

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
				if (info->c.PARAVIEW_CRACK_REPRESENTATION)
				{
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
				}

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
			info->disp_at_connectivity[i * info->DIMENSION + j] = 0.0;
		}
		for (int j = 0; j < D_MATRIX_SIZE; j++)
		{
			info->elastic_strain_at_connectivity[i * D_MATRIX_SIZE + j] = 0.0;
			info->stress_at_connectivity[i * D_MATRIX_SIZE + j] = 0.0;
			info->back_stress_at_connectivity[i * D_MATRIX_SIZE + j] = 0.0;
		}
		info->vm_at_connectivity[i] = 0.0;
		info->hs_at_connectivity[i] = 0.0;
		info->equivalent_plastic_strain_at_connectivity[i] = 0.0;
		info->yield_stress_at_connectivity[i] = 0.0;
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

		int NG_mod = info->c.NG % 2;
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


void Make_info_for_viewer(information *info)
{
	gp_switch(true, info);

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

	// make average points and weights, Connectivity_coord
	static bool isInitialised = false;
	static vector<average_points> ap(Total_connectivity);
	if (!isInitialised)
	{
		Make_average_points(ap, info);
		isInitialised = true;
	}

	#pragma omp parallel for
	for (int i = 0; i < Total_connectivity; i++)
	{
		double temp_point[MAX_DIMENSION] = {0.0};
		double temp_point_glo[MAX_DIMENSION] = {0.0};
		double temp_para_glo[MAX_DIMENSION] = {0.0};
		double temp_point_loc[MAX_DIMENSION] = {0.0};
		double temp_para_loc[MAX_DIMENSION] = {0.0};

		int element = info->Connectivity_ele[i];
		int point = info->Connectivity_point[i];

		// make temp_point
		for (int j = 0; j < info->DIMENSION; j++)
			temp_point[j] = point_array[point * info->DIMENSION + j];

		// make physical coordinate, disp_at_connectivity
		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), temp_point, element, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				double d = info->disp[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
				info->disp_at_connectivity[i * info->DIMENSION + k] += R[j] * d;
			}

		// overlay global
		int status_glo_overlay = 0;
		if (i < Total_connectivity_glo && Total_mesh > 1)
		{
			// make temp_point_glo
			int itr_n = 0, loc_patch = 0;
			for (int j = info->Total_Patch_to_mesh[1]; j < info->Total_Patch_to_mesh[Total_mesh]; j++)
			{
				itr_n = calc_patch_parameter_coord(&info->Connectivity_coord[i * info->DIMENSION], j, temp_para_loc, info);
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
				shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; j++)
					for (int k = 0; k < info->DIMENSION; k++)
					{
						double d_loc = info->disp[info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
						info->disp_at_connectivity[i * info->DIMENSION + k] += R_loc[j] * d_loc;
					}
			}
		}

		// overlay local
		if (i >= Total_connectivity_glo && Total_mesh > 1)
		{
			// make temp_point_glo
			int itr_n = 0, glo_patch = 0;
			for (int j = 0; j < info->Total_Patch_on_mesh[0]; j++)
			{
				itr_n = calc_patch_parameter_coord(&info->Connectivity_coord[i * info->DIMENSION], j, temp_para_glo, info);
				glo_patch = j;
				if (itr_n != ERROR)
					break;
			}

			int element_glo = ele_check(glo_patch, temp_para_glo, info);
			tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

			// overlay displacement
			vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					double d_glo = info->disp[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					info->disp_at_connectivity[i * info->DIMENSION + k] += R_glo[j] * d_glo;
				}
		}

		for (int j = 0; j < ap[i].n; j++)
		{
			int e = ap[i].target_e[j];
			int gp_num = ap[i].target_gp_num[j];

			for (int l = 0; l < D_MATRIX_SIZE; l++)
			{
				info->elastic_strain_at_connectivity[i * N_STRESS + l] += info->gp[e].elastic_strain()[gp_num * D_MATRIX_SIZE + l] * ap[i].w[j];
				info->stress_at_connectivity[i * N_STRESS + l] += info->gp[e].stress()[gp_num * D_MATRIX_SIZE + l] * ap[i].w[j];
				info->back_stress_at_connectivity[i * N_STRESS + l] += info->gp[e].back_stress()[gp_num * D_MATRIX_SIZE + l] * ap[i].w[j];
			}

			info->vm_at_connectivity[i] += info->gp[e].equivalent_stress()[gp_num] * ap[i].w[j];
			info->equivalent_plastic_strain_at_connectivity[i] += info->gp[e].equivalent_plastic_strain()[gp_num] * ap[i].w[j];
			info->yield_stress_at_connectivity[i] += info->gp[e].yield_stress()[gp_num] * ap[i].w[j];
		}

		// hydrostatic stress
		if (info->DIMENSION == 2)
		{
			cout << "Error: 2D is not implemented for hydrostatic stress." << endl;
			exit(1);
		}
		else if (info->DIMENSION == 3)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				info->hs_at_connectivity[i] += info->stress_at_connectivity[i * N_STRESS + j];
			info->hs_at_connectivity[i] /= 3.0;
		}
	}
}


void Make_boundary_line(information *info)
{
	int i, j, k;
	int point_on_edge = 3;

	for (int loop = 0; loop < (info->c.OUTPUT_DEFORMED == 0 ? 2 : 1); loop++)
	{
		int counter = 0;
		double *E_ptr, *C_ptr, *D_ptr;
		if (loop == 0)
		{
			E_ptr = info->Edge_coord;
			C_ptr = info->Connectivity_coord;
			D_ptr = info->disp_at_connectivity;
		}
		else if (loop == 1)
		{
			E_ptr = info->deformed_Edge_coord;
			C_ptr = info->deformed_Connectivity_coord;
		}

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
}


// nonlinear postprocessing
void NL_post_information(information *info)
{
	cout << "start NL_post_information" << endl;
	gp_switch(true, info);

	// define line
	// double point_start[MAX_DIMENSION] = {0.5, 3.5, 3.0};
	// double point_end[MAX_DIMENSION] = {0.5, 6.0, 3.0};
	// double point_start[MAX_DIMENSION] = {1.0, 1.0, 0.0};
	// double point_end[MAX_DIMENSION] = {1.0, 3.0, 0.0};

	// hole
	// double point_start[MAX_DIMENSION] = {0.0, 1.0, 0.0};
	// double point_end[MAX_DIMENSION] = {0.0, 3.0, 0.0};
	// double point_start[MAX_DIMENSION] = {0.0, 1.0, 0.0};
	// double point_end[MAX_DIMENSION] = {0.0, 5.0, 0.0};

	// crack
	// double point_start[MAX_DIMENSION] = {0.25, 3.5, 3.0};
	// double point_end[MAX_DIMENSION] = {0.25, 6.0, 3.0};
	double point_start[MAX_DIMENSION] = {0.0, 0.5, 0.0};
	double point_end[MAX_DIMENSION] = {0.0, 3.0, 0.0};

	// double point_start[MAX_DIMENSION] = {0.5, 1.0, 0.0};
	// double point_end[MAX_DIMENSION] = {0.5, 3.0, 0.0};
	double line_vector[MAX_DIMENSION] = {0.0};

	// discretization of line
	double line_discretization = 501.0;

	for (int i = 0; i < info->DIMENSION; i++)
	{
		line_vector[i] = point_end[i] - point_start[i];
	}

	double *point_array = (double *)malloc(sizeof(double) * info->DIMENSION * (int)line_discretization);
	for (int i = 0; i < (int)line_discretization; i++)
	{
		for (int j = 0; j < info->DIMENSION; j++)
			point_array[i * info->DIMENSION + j] = point_start[j] + line_vector[j] * i / (line_discretization - 1.0);
	}

	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_1.txt", 1);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_2.txt", 2);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_4.txt", 4);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_8.txt", 8);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_16.txt", 16);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_32.txt", 32);
		}
	}
	for (int i = 0; i < Total_mesh; i++)
	{
		for (int j = 0; j < (int)line_discretization; j++)
		{
			data_from_point(i, &point_array[j * info->DIMENSION], j, info, "line_data_64.txt", 64);
		}
	}

	cout << "end NL_post_information" << endl;
}


void data_from_point(int line_overlaying_mesh_num, double *point, int delta, information *info, const string &filename, const int average_points)
{
	// make connectivityで作成した要素の情報から，任意の点の情報を補間する
	int itr_n = 0, line_overlaying_patch = 0, status = 0;
	double para[MAX_DIMENSION];
	for (int i = info->Total_Patch_to_mesh[line_overlaying_mesh_num]; i < info->Total_Patch_to_mesh[line_overlaying_mesh_num + 1]; i++)
	{
		itr_n = calc_patch_parameter_coord(point, i, para, info);
		line_overlaying_patch = i;
		if (itr_n != ERROR)
		{
			status = 1;
			break;
		}
	}

	if (!status)
		return;
	
	int e = ele_check(line_overlaying_patch, para, info);

	vector<distance_norm> dn;
	#pragma omp parallel for
	for (size_t i = 0; i < info->ecn[e].size(); i++)
	{
		int target_e = info->ecn[e][i];

		for (int j = 0; j < info->gp[target_e].n(); j++)
		{
			distance_norm temp_dn;
			temp_dn.e = target_e;
			temp_dn.gp_num = j;
			double GP_coord[MAX_DIMENSION] = {0.0};
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &info->gp[target_e].para()[j * info->DIMENSION], target_e, info);
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[target_e]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
					GP_coord[l] += R[k] * info->Node_Coordinate[info->Controlpoint_of_Element[target_e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
			for (int k = 0; k < info->DIMENSION; k++)
				temp_dn.r += (GP_coord[k] - point[k]) * (GP_coord[k] - point[k]);
			temp_dn.r = sqrt(temp_dn.r);

			#pragma omp critical
			{
				dn.emplace_back(temp_dn);
			}
		}
	}

	// sort
	sort(dn.begin(), dn.end());

	// average points
	int ap_n = (average_points <= static_cast<int>(dn.size())) ? average_points : static_cast<int>(dn.size());
	double total_w = 0.0;

	vector<double> weight(ap_n);
	for (int i = 0; i < ap_n; i++)
	{
		// total_w += 1.0 / dn[i].r;
		total_w += 1.0 / pow(dn[i].r, 2.0);
	}
	for (int i = 0; i < ap_n; i++)
	{
		// weight[i] = 1.0 / dn[i].r / total_w;
		weight[i] = 1.0 / pow(dn[i].r, 2.0) / total_w;
	}

	double equivalent_stress = 0.0;
	double equivalent_plastic_strain = 0.0;
	double stress_11 = 0.0;
	double stress_22 = 0.0;
	double stress_33 = 0.0;
	double strain_11 = 0.0;
	double strain_22 = 0.0;
	double strain_33 = 0.0;

	for (int i = 0; i < ap_n; i++)
	{
		equivalent_stress += info->gp[dn[i].e].equivalent_stress()[dn[i].gp_num] * weight[i];
		equivalent_plastic_strain += info->gp[dn[i].e].equivalent_plastic_strain()[dn[i].gp_num] * weight[i];
		stress_11 += info->gp[dn[i].e].stress()[dn[i].gp_num * D_MATRIX_SIZE + 0] * weight[i];
		stress_22 += info->gp[dn[i].e].stress()[dn[i].gp_num * D_MATRIX_SIZE + 1] * weight[i];
		stress_33 += info->gp[dn[i].e].stress()[dn[i].gp_num * D_MATRIX_SIZE + 2] * weight[i];
		strain_11 += info->gp[dn[i].e].elastic_strain()[dn[i].gp_num * D_MATRIX_SIZE + 0] * weight[i];
		strain_22 += info->gp[dn[i].e].elastic_strain()[dn[i].gp_num * D_MATRIX_SIZE + 1] * weight[i];
		strain_33 += info->gp[dn[i].e].elastic_strain()[dn[i].gp_num * D_MATRIX_SIZE + 2] * weight[i];
	}

	// write data
	FILE *fp_line = fopen(filename.c_str(), "a");
	if (delta == 0)
		fprintf(fp_line, "mesh, Point_x, y, z, delta, Equivalent_stress, Equivalent_plastic_strain, stress_11, 22, 33, strain_11, 22, 33\n");
	fprintf(fp_line, "%d ", line_overlaying_mesh_num);
	for (int i = 0; i < info->DIMENSION; i++)
	{
		fprintf(fp_line, "%.20e ", point[i]);
	}
	fprintf(fp_line, "%d ", delta);
	fprintf(fp_line, "%.20e %.20e ", equivalent_stress, equivalent_plastic_strain);
	fprintf(fp_line, "%.20e %.20e %.20e ", stress_11, stress_22, stress_33);
	fprintf(fp_line, "%.20e %.20e %.20e\n", strain_11, strain_22, strain_33);
	fclose(fp_line);
}