// header
#include "_header.hpp"
#include "_sub.hpp"
#include "_MI3D.hpp"

using namespace std;

// When a cutback occurs, record the step width and treat it
// as the maximum allowed width for subsequent analyses.
static double g_locked_max_delta = 1.0e300;
static bool g_locked_max_set = false;

static inline void set_locked_max_delta(double d)
{
	if (!g_locked_max_set)
	{
		g_locked_max_delta = d;
		g_locked_max_set = true;
	}
	else if (d < g_locked_max_delta)
	{
		// keep the smallest observed cutback width as the maximum
		g_locked_max_delta = d;
	}
}

// Nonlinear
void NL_loop(information *info, char *option_file, char *global_file)
{
	if (info->DIMENSION == 2)
	{
		printf("change dimension, this program is not addopted to 'plane stress' or 'plane strain'\n");
		exit(1);
	}

	// moment boundary condition
	#if 0
		Make_gauss_array(1, info);
		double moment_sigma[9];
		double coord[MAX_DIMENSION] = {0.0};
		double *B_linear = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);

		// init
		for (int i = 0; i < K_Whole_Size; i++)
		{
			info->rhs_vec[i] = 0.0;
		}
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < GP_ON_ELEMENT; j++)
			{
				// make B
				double J = Make_Jac_anypoint(i, &info->Gxi[j * info->DIMENSION], info);
				Make_B_Matrix_anypoint(i, B_linear, &info->Gxi[j * info->DIMENSION], info);

				// init
				for (int k = 0; k < info->DIMENSION; k++)
					coord[k] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
				{
					double R = Shape_func(k, &info->Gxi[j * info->DIMENSION], i, info);
					for (int l = 0; l < info->DIMENSION; l++)
						coord[l] += R * info->Node_Coordinate[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + l];
				}

				// make moment sigma
				double y = coord[2];
				double L = 200.0;
				double young_E = E;
				double pi = 3.14159265358979311600;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					moment_sigma[k] = 0.0;
				moment_sigma[1] = (2.0 * pi * young_E) * y / L;

				// make rhs vec
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < D_MATRIX_SIZE; m++)
						{
							int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							if (id >= 0)
								info->rhs_vec[id] += B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * moment_sigma[m] * J * info->w[j];
						}
			}
		}

		for (int i = 0; i < K_Whole_Size; i++)
		{
			cout << i << " " << info->rhs_vec[i] << endl;;
		}
		Make_gauss_array(0, info);
	#endif

	// get material information
	get_material(info);

	// get load disp curve calculation performed surface
	get_load_surface(info);

	// allocation
	Allocation_init_NL(info);

	// get blending data
	if (info->c.BLENDING == 1 && Total_mesh >= 2)
	{
		get_blending_data(info);
		get_auxiliary(info, option_file, global_file);
		setBlendingArea(info);

		// debug blending function
		#if 0
		debugBlendingFunctionGlo(info, 0);
		debugBlendingFunctionLoc(info, 0);
		debugBlendingFunctionAux(info, 0);
		printf("finish 'debugBlendingFunction'\n\n");
		#endif

		// area
		#if 0
		{
			gp_switch(true, info);

			double volume = 0.0;
			double volume_glo = 0.0;
			double volume_loc = 0.0;
			double volume_aux = 0.0;
			// glo
			#pragma omp parallel for reduction(+:volume_glo) collapse(1)
			for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
			{
				for (int j = 0; j < info->gp[i].n(); j++)
				{
					if (info->gp[i].isSkip()[j])
						continue;

					double J = Make_Jac_anypoint(i, &info->gp[i].para()[j * info->DIMENSION], info);
					double b = info->gp[i].blending_function()[j];
					double v = J * info->gp[i].w()[j] * b;
					#pragma omp atomic
					volume_glo += v;
				}
			}
			cout << "volume_glo = " << volume_glo << endl;

			// loc
			#pragma omp parallel for reduction(+:volume_loc) collapse(1)
			for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (int j = 0; j < info->gp[i].n(); j++)
				{
					if (info->gp[i].isSkip()[j])
						continue;

					double J = Make_Jac_anypoint(i, &info->gp[i].para()[j * info->DIMENSION], info);
					double b = info->gp[i].blending_function()[j];
					double v = J * info->gp[i].w()[j] * b;
					#pragma omp atomic
					volume_loc += v;
				}
			}
			cout << "volume_loc = " << volume_loc << endl;

			// aux
			information *aux = info->aux.get();
			gp_switch(true, aux);
			#pragma omp parallel for reduction(+:volume_aux) collapse(1)
			for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (int j = 0; j < aux->gp[i].n(); j++)
				{
					if (aux->gp[i].isSkip()[j])
						continue;

					double J = Make_Jac_anypoint(i, &aux->gp[i].para()[j * aux->DIMENSION], aux);
					double b = aux->gp[i].blending_function()[j];
					double v = J * aux->gp[i].w()[j] * b;
					#pragma omp atomic
					volume_aux += v;
				}
			}
			cout << "volume_aux = " << volume_aux << endl;

			volume = volume_glo + volume_loc + volume_aux;

			printf("volume = %.20e\n", volume);

			double true_volume = 4.10730091830127541641e+00;
			printf("diff = %.20e\n", volume - true_volume);
			printf("error = %.20e\n", (volume - true_volume) / true_volume);
		}
		#endif
	}

	// get decreased order data
	if (info->c.B_BAR == 1)
	{
		get_decreased_order_data(info, option_file);
		cout << "Finish 'get_decreased_order_data'" << endl;
	}

	// for output
	Allocation(7, info);
	Allocation(8, info);
	Make_connectivity(info);

	printf("Finish 'Make_connectivity'\n\n");

	std::vector<std::shared_ptr<LoadStepNode>> step_roots;
	int N = info->c.STEP_N;

	for (int i = 0; i < N; i++)
	{
		double s = (double)i / N;
		double e = (double)(i + 1) / N;
		step_roots.push_back(std::make_shared<LoadStepNode>(s, e, i));
	}

	const double min_step_size = 1e-4;

	int timestep_counter = 0;
	output_for_paraview_timestep(info, 0.0);

	cout << "\033[36m";
	printf("load step %d / %d\n", ++timestep_counter - 1, info->c.STEP_N);
	cout << "\033[0m";

	// solve
	for (auto& root : step_roots)
	{
		solve_tree(root, info, min_step_size);
		std::cout << "\033[36m";
		printf("load step %d / %d\n", ++timestep_counter - 1, info->c.STEP_N);
		std::cout << "\033[0m";
	}

	#if 1
	NL_post_information(info);
	#endif

	#if 0
	double stress_exact[6] = {1.559481185031401e+02, 1.559481185031401e+02, 3.638789431739937e+02, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00};
	double elast_exact[6]  = {0.0};
	elast_exact[2] = 0.5 * log(9.0 / 4.0);

	double *BL = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double error_norm_a = 0.0;
	double error_norm_b = 0.0;
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < GP_ON_ELEMENT; j++)
		{
			double J = Make_B_Linear(i, &info->Gxi[j * info->DIMENSION], BL, info);
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				error_norm_a += (info->elastic_strain[(long)i * GP_ON_ELEMENT * D_MATRIX_SIZE + j * D_MATRIX_SIZE + k] - elast_exact[k]) * (info->stress[(long)i * GP_ON_ELEMENT * D_MATRIX_SIZE + j * D_MATRIX_SIZE + k] - stress_exact[k]) * J * info->w[j];
				error_norm_b += info->elastic_strain[(long)i * GP_ON_ELEMENT * D_MATRIX_SIZE + j * D_MATRIX_SIZE + k] * info->stress[(long)i * GP_ON_ELEMENT * D_MATRIX_SIZE + j * D_MATRIX_SIZE + k] * J * info->w[j];
			}
		}
	}
	
	error_norm_a = sqrt(error_norm_a);
	error_norm_b = sqrt(error_norm_b);

	cout << "\033[31m";
	printf("%.20e\n", error_norm_a);
	printf("%.20e\n", error_norm_b);
	printf("%.20e\n", error_norm_a / error_norm_b);
	cout << "\033[0m";
	#endif

	#if 0

	double *BL = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double a = 0.0;
	double b = 0.0;

	double J = 0.0;
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < GP_ON_ELEMENT; j++)
		{
			double Jac = Make_B_Linear(i, &info->Gxi[j * info->DIMENSION], BL, info);
			double det_F = calc_3x3_determinant(&info->gp[i].deformation_gradient()[j * info->DIMENSION * info->DIMENSION]);
			printf("%d %d %.20e %.20e %.20e %.20e\n", i, j, info->equivalent_stress[i * GP_ON_ELEMENT + j], info->yield_stress[i * GP_ON_ELEMENT + j], info->equivalent_plastic_strain[i * GP_ON_ELEMENT + j], det_F);

			a += info->equivalent_stress[i * GP_ON_ELEMENT + j] * Jac * info->w[j];
			b += info->equivalent_plastic_strain[i * GP_ON_ELEMENT + j] * Jac * info->w[j];
			J += Jac * info->w[j];
		}
	}

	a /= J;
	b /= J;

	cout << "\033[31m";
	printf("%.20e\n", a);
	printf("%.20e\n", b);
	// printf("%.20e\n", error_norm_a / error_norm_b);
	cout << "\033[0m";

	#endif

	// end analysis
	return;
}


void Allocation_init_NL(information *info)
{
	// memory allocation
	info->delta_F = (double *)calloc(K_Whole_Size, sizeof(double));
	info->F = (double *)calloc(K_Whole_Size, sizeof(double));
	info->residual_vec = (double *)calloc(K_Whole_Size, sizeof(double));
	info->internal_force = (double *)malloc(sizeof(double) * K_Whole_Size);
	info->external_force = (double *)malloc(sizeof(double) * K_Whole_Size);
	info->previous_external_force = (double *)calloc(K_Whole_Size, sizeof(double));
	info->rhs_vec_initial = (double *)calloc(sizeof(double) * K_Whole_Size, sizeof(double));
	info->forced_disp_T = (double *)calloc(K_Whole_Size, sizeof(double));
	info->disp_increment = (double *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(double));
	info->disp_overlay = (double *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(double));
	info->disp_overlay_increment = (double *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(double));
	info->disp = (double *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(double));
	info->disp_previous = (double *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(double));

	// B-bar
	info->b_average = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION * MAX_NO_CP_ON_ELEMENT, sizeof(double));

	if (info->delta_F == NULL || info->F == NULL || info->residual_vec == NULL || info->internal_force == NULL || info->external_force == NULL || info->rhs_vec_initial == NULL || info->previous_external_force == NULL || info->forced_disp_T == NULL || info->disp_increment == NULL || info->disp_overlay == NULL || info->disp_overlay_increment == NULL || info->disp == NULL || info->disp_previous == NULL || info->b_average == NULL)
	{
		printf("Cannot allocate memory\n");
		exit(1);
	}

	// init
	gp_switch(true, info);

	#pragma omp parallel for
	for(int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		// init deformation_gradient
		for (int j = 0; j < info->gp[i].n(); j++)
			for (int k = 0; k < info->DIMENSION; k++)
				info->gp[i].deformation_gradient()[j * info->DIMENSION * info->DIMENSION + k * info->DIMENSION + k] = 1.0;

		// init yield stress
		for (int j = 0; j < info->gp[i].n(); j++)
			info->gp[i].yield_stress()[j] = get_hardening_stress(0.0, info);

		// init equivalent plastic strain
		for (int j = 0; j < info->gp[i].n(); j++)
		{
			info->gp[i].equivalent_plastic_strain()[j] = 0.0;
			info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
		}
	}

	// init delta_F
	for (int i = 0; i < K_Whole_Size; i++)
		info->delta_F[i] = info->rhs_vec[i];

	// set Displacement to disp
	info->Displacement = info->disp;

	// B-bar init
	#if 0
	if (info->c.B_BAR == 1)
	{
		cout << "start B-bar init" << endl;
		#if 1
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			double *b_ave = &info->b_average[i * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];
			vector<double> b_temp(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
			double ele_volume = 0.0;
			for (int j = 0; j < info->gp[i].n(); j++)
			{
				double J = Make_B_component(i, &info->gp[i].para()[j * info->DIMENSION], b_temp.data(), info);
				double coef = J * info->gp[i].w()[j];
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; l++)
						b_ave[k * MAX_NO_CP_ON_ELEMENT + l] += b_temp[k * MAX_NO_CP_ON_ELEMENT + l] * coef;
				ele_volume += coef;
			}
			double ele_volume_inv = 1.0 / ele_volume;
			for (int j = 0; j < info->DIMENSION; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					b_ave[j * MAX_NO_CP_ON_ELEMENT + k] *= ele_volume_inv;
		}
		#elif 0
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			double *b_ave = &info->b_average[i * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];
			vector<double> b_temp(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
			int gp_counter = 0;
			for (int j = 0; j < info->gp[i].n(); j++)
			{
				Make_B_component(i, &info->gp[i].para()[j * info->DIMENSION], b_temp.data(), info);
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; l++)
						b_ave[k * MAX_NO_CP_ON_ELEMENT + l] += b_temp[k * MAX_NO_CP_ON_ELEMENT + l];
				gp_counter++;
			}
			double gp_counter_inv = 1.0 / static_cast<double>(gp_counter);
			for (int j = 0; j < info->DIMENSION; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					b_ave[j * MAX_NO_CP_ON_ELEMENT + k] *= gp_counter_inv;
		}
		#elif 0
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			double *b_ave = &info->b_average[i * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];
			vector<double> b_temp(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
			double ele_volume = 0.0;
	
			double para[3] = {0.0};
			double J = Make_B_component(i, para, b_temp.data(), info);
			for (int k = 0; k < info->DIMENSION; k++)
				for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; l++)
					b_ave[k * MAX_NO_CP_ON_ELEMENT + l] += b_temp[k * MAX_NO_CP_ON_ELEMENT + l];
		}
		#endif
		cout << "finish B-bar init" << endl;
	}
	#endif
}


void Update_FR(information *info, double load_factor)
{
	// update F
	for (int i = 0; i < K_Whole_Size; i++)
		info->F[i] = info->delta_F[i] * load_factor;

	// update residual_vec
	for (int i = 0; i < K_Whole_Size; i++)
		info->residual_vec[i] = info->delta_F[i] * load_factor;
}


void Update_increment_field(information *info)
{
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		info->disp_previous[i] = info->disp[i];
		info->disp[i] += info->disp_increment[i];
		info->disp_increment[i] = 0.0;
		info->disp_overlay[i] += info->disp_overlay_increment[i];
		info->disp_overlay_increment[i] = 0.0;
	}

	gp_switch(true, info);

	#pragma omp parallel for collapse(1)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < info->gp[i].n(); j++)
		{
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				info->gp[i].previous_elastic_strain()[j * D_MATRIX_SIZE + k] = info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k];
				info->gp[i].previous_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].stress()[j * D_MATRIX_SIZE + k];
			}

			for (int k = 0; k < (info->DIMENSION * info->DIMENSION); k++)
				info->gp[i].deformation_gradient()[j * info->DIMENSION * info->DIMENSION + k] = info->gp[i].current_deformation_gradient()[j * info->DIMENSION * info->DIMENSION + k];

			for (int k = 0; k < D_MATRIX_SIZE; k++)
				info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k] = info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k];

			for (int k = 0; k < D_MATRIX_SIZE; k++)
				info->gp[i].stress()[j * D_MATRIX_SIZE + k] = info->gp[i].current_stress()[j * D_MATRIX_SIZE + k];
			info->gp[i].equivalent_stress()[j] = calc_equivalent_stress(&info->gp[i].stress()[j * D_MATRIX_SIZE]);

			info->gp[i].equivalent_plastic_strain()[j] += info->gp[i].equivalent_plastic_strain_increment()[j];
			info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;

			info->gp[i].yield_stress()[j] = info->gp[i].current_yield_stress()[j];

			for (int k = 0; k < D_MATRIX_SIZE; k++)
				info->gp[i].back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k];
		}
	}

	// auxiliary patch
	if (Total_mesh >= 2 && info->c.BLENDING == 1)
	{
		information *aux = info->aux.get();
		for (int i = 0; i < aux->Total_Control_Point_to_mesh[1] * aux->DIMENSION; i++)
		{
			aux->disp_previous[i] = aux->disp[i];
			aux->disp[i] += aux->disp_increment[i];
			aux->disp_increment[i] = 0.0;
			aux->disp_overlay[i] += info->disp_overlay_increment[i];
			aux->disp_overlay_increment[i] = 0.0;
		}
		for (int i = aux->Total_Control_Point_to_mesh[1] * aux->DIMENSION; i < aux->Total_Control_Point_to_mesh[Total_mesh] * aux->DIMENSION; i++)
		{
			aux->disp_previous[i] = aux->disp[i];
			aux->disp[i] = 0.0;
			aux->disp_increment[i] = 0.0;
			aux->disp_overlay[i] += aux->disp_overlay_increment[i];
			aux->disp_overlay_increment[i] = 0.0;
		}

		gp_switch(true, aux);

		#pragma omp parallel for collapse(1)
		for (int i = 0; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < aux->gp[i].n(); j++)
			{
				for (int k = 0; k < D_MATRIX_SIZE; k++)
				{
					aux->gp[i].previous_elastic_strain()[j * D_MATRIX_SIZE + k] = aux->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k];
					aux->gp[i].previous_stress()[j * D_MATRIX_SIZE + k] = aux->gp[i].stress()[j * D_MATRIX_SIZE + k];
				}

				for (int k = 0; k < (aux->DIMENSION * aux->DIMENSION); k++)
					aux->gp[i].deformation_gradient()[j * aux->DIMENSION * aux->DIMENSION + k] = aux->gp[i].current_deformation_gradient()[j * aux->DIMENSION * aux->DIMENSION + k];

				for (int k = 0; k < D_MATRIX_SIZE; k++)
					aux->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k] = aux->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k];

				for (int k = 0; k < D_MATRIX_SIZE; k++)
					aux->gp[i].stress()[j * D_MATRIX_SIZE + k] = aux->gp[i].current_stress()[j * D_MATRIX_SIZE + k];
				aux->gp[i].equivalent_stress()[j] = calc_equivalent_stress(&aux->gp[i].stress()[j * D_MATRIX_SIZE]);

				aux->gp[i].equivalent_plastic_strain()[j] += aux->gp[i].equivalent_plastic_strain_increment()[j];
				aux->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;

				aux->gp[i].yield_stress()[j] = aux->gp[i].current_yield_stress()[j];

				for (int k = 0; k < D_MATRIX_SIZE; k++)
					aux->gp[i].back_stress()[j * D_MATRIX_SIZE + k] = aux->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k];
			}
		}
	}
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

	// auxiliary patch
	if (Total_mesh >= 2 && info->c.BLENDING == 1)
	{
		information *aux = info->aux.get();
		for (int i = 0; i < aux->Total_Control_Point_to_mesh[1] * aux->DIMENSION; i++)
		{
			aux->disp_increment[i] = 0.0;
			aux->disp_overlay_increment[i] = 0.0;
		}
		for (int i = aux->Total_Control_Point_to_mesh[1] * aux->DIMENSION; i < aux->Total_Control_Point_to_mesh[Total_mesh] * aux->DIMENSION; i++)
		{
			aux->disp[i] = 0.0;
			aux->disp_increment[i] = 0.0;
			aux->disp_overlay_increment[i] = 0.0;
		}

		gp_switch(true, aux);

		#pragma omp parallel for collapse(1)
		for (int i = 0; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < aux->gp[i].n(); j++)
				aux->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
		}
	}
}


void get_load_surface(information *info)
{
	if (info->opt_files[3].empty())
		return;

	int temp_i;

	std::ifstream ifs(info->opt_files[3]);
	if (!ifs)
	{
		std::cerr << "Cannot open load surface file: " << info->opt_files[3] << std::endl;
		exit(1);
	}

	// load disp curve calculation performed surface
	ifs >> temp_i;
	int surface_n = temp_i;

	// allocation
	info->load_disp_element = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);
	info->load_disp_element_direction = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);
	info->load_disp = (int *)calloc(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, sizeof(int));
	vector<int> patch_num(surface_n);
	vector<int> load_disp_cp(info->Total_Control_Point_to_mesh[Total_mesh]);

	ifs >> temp_i; int direction = temp_i;
	ifs >> temp_i; int position = temp_i;
	ifs >> temp_i; int physical_direction = temp_i;

	for (int i = 0; i < surface_n; i++)
	{
		ifs >> temp_i;
		patch_num[i] = temp_i;
	}

	ifs.close();

	// make element
	int cp_counter = 0;
	info->load_disp_element_n = 0;
	for (int i = 0; i < surface_n; i++)
	{
		int start_element = 0;
		for (int j = 0; j < patch_num[i]; j++)
		{
			int temp_ele_n = 1;
			for (int k = 0; k < info->DIMENSION; k++)
				temp_ele_n *= info->No_Control_point[j * info->DIMENSION + k] - info->Order[j * info->DIMENSION + k];
			start_element += temp_ele_n;
		}

		int end_element = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			end_element *= info->No_Control_point[patch_num[i] * info->DIMENSION + j] - info->Order[patch_num[i] * info->DIMENSION + j];
		end_element += start_element;

		for (int j = start_element; j < end_element; j++)
		{
			int ele_posi = 0;
			if (position == 0)
				ele_posi = 0;
			else if (position == 1)
				ele_posi = info->No_Control_point[patch_num[i] * info->DIMENSION + direction] - info->Order[patch_num[i] * info->DIMENSION + direction] - 1;
			else
			{
				printf("please set correct position in file load_surface.\n");
				exit(1);
			}
			
			if (info->ENC[j * info->DIMENSION + direction] == ele_posi)
			{
				int cp_posi = 0;
				if (position == 0)
					cp_posi = 0;
				else if (position == 1)
					cp_posi = info->No_Control_point[patch_num[i] * info->DIMENSION + direction] - 1;
	
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[j]]; k++)
				{
					int cp = info->Controlpoint_of_Element[j * MAX_NO_CP_ON_ELEMENT + k];
					if (info->INC[info->Element_patch[j] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + cp * info->DIMENSION + direction] == cp_posi)
					{
						load_disp_cp[cp_counter] = cp;
						cp_counter++;
					}
				}

				info->load_disp_element_direction[info->load_disp_element_n] = direction;
				info->load_disp_element[info->load_disp_element_n] = j;
				info->load_disp_element_n++;
			}
		}
	}

	for (int i = 0; i < cp_counter; i++)
	{
		info->load_disp[load_disp_cp[i] * info->DIMENSION + physical_direction] = 1;
	}
}


void Load_disp_curve(information *info, int step)
{
	if (info->opt_files[3].empty())
		return;

	gp_switch(true, info);

	static bool check_bool = 0;
	static double max_forced_disp = 0.0;
	if (!check_bool)
	{
		for (int i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
		{
			if (max_forced_disp < info->Value_of_Constraint[i])
				max_forced_disp = info->Value_of_Constraint[i];
		}
		check_bool = true;
	}

	double delta_disp = max_forced_disp * (double)(step + 1) / (double)info->c.STEP_N;
	// double delta_disp = max_forced_disp * factor;

	double f = 0.0;
	#if 0
	// static double *f_vec = (double *)malloc(sizeof(double) * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION);
	static double *K_linear = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
	static double *K_nonlinear = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	// make f
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int row_count = 0;
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->load_disp[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
				if (id == 1)
					row_count++;
			}

		if (row_count > 0)
		{
			// K_linear
			Calc_K_linear_EL(i, K_linear, info);

			// K_nonlinear
			if (info->c.ANALYSIS_MODE_1 == 1)
			{
				Calc_K_nonlinear_EL(i, K_nonlinear, info);
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
						K_linear[j * MAX_KIEL_SIZE + k] += K_nonlinear[j * MAX_KIEL_SIZE + k];
			}

			for (int j = 0; j < info->DIMENSION; j++)
			{
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
				{
					int ii_id = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + j;
					int ii_local = k * info->DIMENSION + j;
					if (info->load_disp[ii_id] == 1)
					{
						for (int l = 0; l < info->DIMENSION; l++)
						{
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
							{
								int jj_id = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + m;
								int jj_local = m * info->DIMENSION + l;
								f += K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->disp[jj_id] + info->disp_increment[jj_id]);
							}
						}
					}
				}
			}
		}
	}
	#endif

	double temp_f = 0.0;

	#pragma omp parallel for reduction(+:temp_f) collapse(1) schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		thread_local double *B_linear = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);

		for (int j = 0; j < info->gp[i].n(); j++)
		{
			double J = Make_B_Linear(i, &info->gp[i].para()[j * info->DIMENSION], info->gp[i].isOverlay()[j], info->gp[i].opp_ele()[j], &info->gp[i].opp_para_tilde()[j * info->DIMENSION], B_linear, info);

			// B^T stress
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
					for (int m = 0; m < D_MATRIX_SIZE; m++)
						if (info->load_disp[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l] == 1)
								temp_f += B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].stress()[j * D_MATRIX_SIZE + m] * info->gp[i].w()[j] * J;
		}
	}

	FILE *fp_ols = fopen("output_load_surface.txt", "a");
	fprintf(fp_ols, "%.20e %.20e %.20e\n", delta_disp, f, temp_f);
	fclose(fp_ols);
}


void Update_previous_external_force(information *info)
{
	double *swap_ptr;

	swap_ptr = info->external_force;
	info->external_force = info->previous_external_force;
	info->previous_external_force = swap_ptr;
}


void solve_tree(std::shared_ptr<LoadStepNode> node, information *info, double min_size)
{
	if (!node->is_leaf())
	{
		for (auto &child : node->children)
			solve_tree(child, info, min_size);
		return;
	}

	// Check if this node exceeds the locked maximum delta
	double node_size = node->size();
	if (g_locked_max_set && node_size > g_locked_max_delta)
	{
		// This node is too large; skip trial and subdivide immediately
		std::cout << "\033[33m";
		printf("Node [%.6e, %.6e] (size=%.6e) exceeds locked max delta (%.6e). Skipping trial and subdividing.\n",
		       node->start, node->end, node_size, g_locked_max_delta);
		std::cout << "\033[0m";

		Init_increment_field(info);

		if (node_size < min_size)
		{
			std::cerr << "ERROR: step too small\n";
			return;
		}

		node->subdivide();
		for (auto &child : node->children)
			solve_tree(child, info, min_size);
		return;
	}

	// --- 外力ベクトル更新（区間に応じた外力） ---
	Update_FR(info, node->end); // 外力ベクトルに load_factor を反映

	// --- Newton-Raphson ---
	bool is_diverged = Newton_Raphson(node->step_id, node->end, node->size(), info);

	if (!is_diverged)
	{
		node->converged = true;

		// --- 成功したら更新処理 ---
		Update_increment_field(info);
		Update_previous_external_force(info);
		output_for_paraview_timestep(info, node->end);

		// debug blending 等
		#if 0
		if (info->c.ANALYSIS_MODE_1 == 1)
		{
			debugBlendingFunctionGlo(info, node->step_id + 1);
			debugBlendingFunctionLoc(info, node->step_id + 1);
			debugBlendingFunctionAux(info, node->step_id + 1);
			printf("finish 'debugBlendingFunction'\n\n");
		}
		#endif

		// calc J
		calc_J(info);
	}
	else
	{
		if (info->c.CUTBACK == 0)
		{
			cout << "INFO: CUTBACK is disabled (set to 0), Newton-Raphson did not converge. Proceeding without cutback." << endl;
			node->converged = false;

			// --- 更新処理 ---
			Update_increment_field(info);
			Update_previous_external_force(info);
			output_for_paraview_timestep(info, node->end);

			// debug blending 等
			#if 0
			if (info->c.ANALYSIS_MODE_1 == 1)
			{
				debugBlendingFunctionGlo(info, node->step_id + 1);
				debugBlendingFunctionLoc(info, node->step_id + 1);
				debugBlendingFunctionAux(info, node->step_id + 1);
				printf("finish 'debugBlendingFunction'\n\n");
			}
			#endif

			// calc J
			calc_J(info);

			return;
		}
	
		// Record that a cutback happened: lock future maximum delta to the attempted width
		set_locked_max_delta(node->size());

		Init_increment_field(info);

		if (node->size() < min_size)
		{
			std::cerr << "ERROR: step too small\n";
			return;
		}
	
		node->subdivide();
		for (auto &child : node->children)
			solve_tree(child, info, min_size);
	}
}


bool Newton_Raphson(int step, double factor, double delta_factor, information *info)
{
	int stop = 0;
	vector<bool> ep_flag(2, false);
	vector<bool> ep_flag_aux(2, false);
	bool is_diverged = false;

	bool aux_flag = (Total_mesh >= 2 && info->c.BLENDING == 1) ? true : false;
	cout << "aux_flag = " << aux_flag << endl;

	for (int i = 0; i < info->c.NEWTON_RAPHSON_MAX_ITR; i++)
	{
		if (info->c.ANALYSIS_MODE_0 == 0)
		{
			// update internal force
			cout << "start update internal force" << endl;
			Update_internal_force_linear(ep_flag, info);
			if (aux_flag)
				Update_internal_force_linear(ep_flag_aux, info->aux.get());
			cout << "finish update internal force" << endl;

			if (i == 1)
				break;
		}
		else if (info->c.ANALYSIS_MODE_0 == 1)
		{
			#if 1
			// update internal force
			cout << "start update internal force" << endl;
			Update_internal_force(ep_flag, info);
			if (aux_flag)
				Update_internal_force(ep_flag_aux, info->aux.get());
			cout << "finish update internal force" << endl;
			#else
			// update internal force, aux stored values are interpolated from info's global
			cout << "start update internal force" << endl;
			Update_internal_force(ep_flag, info);
			if (aux_flag)
				Update_internal_force_interpolate(info, info->aux.get());
			cout << "finish update internal force" << endl;
			#endif
		}
		else
		{
			printf("please set correct ANALYSIS_MODE_0\n");
			exit(1);
		}

		// update external force
		Update_external_force(info);

		// check convergence
		is_diverged = Convergence(ep_flag, ep_flag_aux, aux_flag, info, i, &stop, step, factor, delta_factor);
		if (is_diverged)
		{
			printf("Newton-Raphson diverged at iteration %d\n", i);
			return true;
		}

		if (stop)
		{
			printf("Newton-Raphson converged at iteration %d\n", i);
			break;
		}
		else if (i == info->c.NEWTON_RAPHSON_MAX_ITR - 1)
		{
			printf("Newton-Raphson didn't converge\n");
			is_diverged = true;
			break;
		}

		// generate coefficient matrix
		Calc_K(info);
		if (aux_flag)
			Calc_K(info->aux.get());

		// update Neumann force vector
		cout << "start Update_force_vector" << endl;
		Update_force_vector(info, factor, i);
		if (aux_flag)
			Update_force_vector(info->aux.get(), factor, i);
		cout << "finish Update_force_vector" << endl;

		// solve by PCG
		if (info->c.SOLVER == 0)
			PCG_Solver(K_Whole_Size, info->c.EPS, info);
		// solve by GMRES
		else if (info->c.SOLVER == 1)
			GMRES(K_Whole_Size, info->c.EPS, info);
		else if (info->c.SOLVER == 2)
			GaussianElimination(info->sol_vec, info->rhs_vec, info->full_K, K_Whole_Size);
		// else if (info->c.SOLVER == 3)
		// 	LU();
		else if (info->c.SOLVER == 4)
			intel_PARDISO(info->sol_vec, info->rhs_vec, K_Whole_Size, info);

		// update displaccement increments
		Update_disp_increment(delta_factor, info);

		if (Total_mesh >= 2 && info->c.ANALYSIS_MODE_1 == 1)
		{
			#if 0
			// ordinary update using overdetermined system
			printf("start cp update\n");
			// Update_Control_Point_glo(info);
			// Update_Control_Point_loc(info);
			Update_Control_Point_glo_old(info);
			Update_Control_Point_loc_old(info);
			printf("finish cp update\n");

			if (aux_flag)
			{
				printf("start aux cp update\n");
				// Update_Control_Point_glo(info);
				// Update_Control_Point_loc(info);
				Update_Control_Point_substitute_glo_aux(info, info->aux.get());
				Update_Control_Point_aux_old(info->aux.get());
				printf("finish aux cp update\n");
			}
			#endif
			#if 0
			// L2 projection
			cout << "start glo cp update using L2 projection\n";
			Update_Control_Point_glo_using_L2_projection(info);
			// Update_Control_Point_glo_using_L2_projection_and_aux(info);
			cout << "start loc cp update using L2 projection\n";
			Update_Control_Point_loc_using_L2_projection(info);
			if (aux_flag)
			{
				Update_Control_Point_substitute_glo_aux(info, info->aux.get());
				cout << "start aux cp update using L2 projection\n";
				Update_Control_Point_aux_using_L2_projection(info->aux.get());
			}
			#endif
		}

		// if (info->c.ANALYSIS_MODE_1 == 1 && Total_mesh >= 2)
		// 	Update_Gauss_points(info);
	}
	return is_diverged;
}


bool Convergence(vector<bool> &ep_flag, vector<bool> &ep_flag_aux, bool &aux_flag, information *info, int itr, int *stop, int step, double factor, double delta_factor)
{
	bool is_diverged = false;
	double residual = 0.0;
	double siga_residual_norm[2] = {0.0};

	// Calculate norms
	double f_norm  = 0.0;
	double df_norm = 0.0;
	double r_norm  = 0.0;
	double siga_r_norm[2] = {0.0, 0.0};

	int index = 0;
	for (int i = info->Total_Control_Point_to_mesh[1] * info->DIMENSION; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		index = info->Index_Dof[i];
		if (index >= 0)
			break;
	}

	for (int i = 0; i < K_Whole_Size; i++)
	{
		f_norm += calc_square(info->external_force[i]);
		df_norm += calc_square(info->external_force[i] - info->previous_external_force[i]);

		double temp = calc_square(info->external_force[i] - info->internal_force[i]);
		r_norm += temp;
		if (i < index)
			siga_r_norm[0] += temp;
		else
			siga_r_norm[1] += temp;
	}
	f_norm = sqrt(f_norm);
	df_norm = sqrt(df_norm);
	r_norm = sqrt(r_norm);
	siga_r_norm[0] = sqrt(siga_r_norm[0]);
	siga_r_norm[1] = sqrt(siga_r_norm[1]);

	// Calculate reference norm
	double reference_norm = fmax(f_norm, df_norm);

	// only Dirichlet boundary
	if (reference_norm <= MERGE_ERROR)
	{
		reference_norm = 0.0;
		r_norm = 0.0;
		siga_r_norm[0] = 0.0;
		siga_r_norm[1] = 0.0;

		for (int i = 0; i < K_Whole_Size; i++)
		{
			reference_norm += calc_square(info->rhs_vec_initial[i]);

			double temp = calc_square(info->external_force[i] - info->internal_force[i]);
			r_norm += temp;
			if (i < index)
				siga_r_norm[0] += temp;
			else
				siga_r_norm[1] += temp;
		}
		reference_norm = sqrt(reference_norm);
		r_norm = sqrt(r_norm);
		siga_r_norm[0] = sqrt(siga_r_norm[0]);
		siga_r_norm[1] = sqrt(siga_r_norm[1]);
	}

	if (itr == 0)
	{
		residual = 1.0;
		for (int i = 0; i < 2; i++)
			siga_residual_norm[i] = 1.0;
	}
	else
	{
		residual = r_norm / reference_norm;
		for (int i = 0; i < 2; i++)
			siga_residual_norm[i] = siga_r_norm[i] / reference_norm;
	}

	cout << "\033[33m";
	printf("\nNewton Raphson (factor = %le, delta_factor = %le)\n", factor, delta_factor);
	if (g_locked_max_set)
	{
		cout << "\033[35m";
		printf("[Locked max delta: %.6e]\n", g_locked_max_delta);
		cout << "\033[33m";
	}
	printf("itr = %d, residual = %1.5le", itr, residual);
	if (Total_mesh > 1)
		printf(", residual glo = %1.5le, residual loc = %1.5le", siga_residual_norm[0], siga_residual_norm[1]);
	printf("\n");
	cout << "\033[0m";

	if (!is_value_finite_robust(residual) || residual > 1.0e3)
	{
		cout << "\033[31m";
		printf("residual is NaN or Inf, Newton-Raphson diverged\n");
		cout << "\033[0m";
		is_diverged = true;
	}

	double sol_norm = 0.0;
	double siga_sol_norm[2] = {0.0};
	for (int i = 0; i < K_Whole_Size; i++)
	{
		double temp = calc_square(info->sol_vec[i]);
		sol_norm += temp;
		if (i < index)
			siga_sol_norm[0] += temp;
		else
			siga_sol_norm[1] += temp;
	}
	sol_norm = sqrt(sol_norm);
	siga_sol_norm[0] = sqrt(siga_sol_norm[0]);
	siga_sol_norm[1] = sqrt(siga_sol_norm[1]);

	cout << "\033[33m";
	printf("itr = %d, sol norm = %1.5le", itr, sol_norm);
	if (Total_mesh > 1)
		printf(", sol norm glo = %1.5le, sol norm loc = %1.5le", siga_sol_norm[0], siga_sol_norm[1]);
	printf("\n\n");
	cout << "\033[0m";

	string s = "";
	if (ep_flag[0] == true)
		s += "e";
	if (ep_flag[1] == true)
		s += "p";
	if (aux_flag == true)
	{
		s += "/";
		if (ep_flag_aux[0] == true)
			s += "e";
		if (ep_flag_aux[1] == true)
			s += "p";
	}

	FILE *fp_rc = fopen("residual_convergence.txt", "a");
	static int count = 0;
	if (Total_mesh == 1)
	{
		if (count++ == 0)
			fprintf(fp_rc, "step itr factor residual_norm sol_norm \n");
		fprintf(fp_rc, "%d %d %le %.20e %.20e %s\n", step, itr, factor, residual, sol_norm, s.c_str());
	}
	else if (Total_mesh > 1)
	{
		if (count++ == 0)
			fprintf(fp_rc, "step itr factor residual_norm sol_norm residual_norm_glo residual_norm_loc sol_norm_glo sol_norm_loc %s\n", s.c_str());
		fprintf(fp_rc, "%d %d %le %.20e %.20e %.20e %.20e %.20e %.20e %s\n", step, itr, factor, residual, sol_norm, siga_residual_norm[0], siga_residual_norm[1], siga_sol_norm[0], siga_sol_norm[1], s.c_str());
	}
	fclose(fp_rc);

	if (residual <= info->c.NEWTON_RAPHSON_EPS)
		*stop = 1;

	static double previous_residual = 0.0;
	if (info->c.NEWTON_RAPHSON_MAX_ITR != 1 && itr == info->c.NEWTON_RAPHSON_MAX_ITR - 1)
	{
		double ratio = residual / previous_residual;

		if (ratio > 1.0)
		{
			cout << "\033[31m";
			printf("residual increased by %1.5le, Newton-Raphson diverged\n", ratio);
			cout << "\033[0m";
			is_diverged = true;
		}
	}
	previous_residual = residual;

	return is_diverged;
}


void Calc_K(information *info)
{
	cout << "start calc K" << endl;

	// auxilary patch
	bool isAuxPatch = false;
	if (Total_mesh > 1 && info->c.BLENDING == 1 && info->aux == nullptr)
		isAuxPatch = true;
	cout << "isAuxPatch = " << isAuxPatch << endl;

	// init
	if (isAuxPatch == false)
	{
		if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
		{
			#pragma omp parallel for
			for (int i = 0; i < K_Whole_Size * K_Whole_Size; i++)
				info->full_K[i] = 0.0;
		}
		else
		{
			#pragma omp parallel for
			for (int i = 0; i < info->K_Whole_Ptr[K_Whole_Size]; i++)
				info->K_Whole_Val[i] = 0.0;
		}
	}

	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<double> K_linear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		vector<double> K_nonlinear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		vector<double> K_linear_conditional(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		vector<double> K_nonlinear_conditional(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		if (Total_mesh == 1 || (Total_mesh > 1 && info->c.BLENDING == 0))
		{
			// K_linear
			Calc_K_linear_EL(i, K_linear.data(), info);

			// K_nonlinear
			if (info->c.ANALYSIS_MODE_1 == 1)
			{
				Calc_K_nonlinear_EL(i, K_nonlinear.data(), info);
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
					{
						int id = j * MAX_KIEL_SIZE + k;
						double val = K_nonlinear[id];

						K_linear[id] += val;
					}
			}
		}
		else if (Total_mesh > 1 && info->c.BLENDING == 1)
		{
			// [K_L]^GG + [K_NL]^GG by global, local and axiliary elements
			if (i < info->Total_Element_to_mesh[1])
			{
				// [K_L]^GG + [K_NL]^GG by global elements
				if (isAuxPatch == false)
				{
					int dummy_element = 0;

					// K_linear based on global elements
					Calc_K_linear_EL_conditional_for_glo(false, i, dummy_element, K_linear.data(), info);

					// K_nonlinear based on global elements
					if (info->c.ANALYSIS_MODE_1 == 1)
					{
						Calc_K_nonlinear_EL_conditional_for_glo(false, i, dummy_element, K_nonlinear.data(), info);
						for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
							for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
							{
								int id = j * MAX_KIEL_SIZE + k;
								double val = K_nonlinear[id];

								K_linear[id] += val;
							}
					}
				}
				else if (isAuxPatch == true)
				{
					// init
					for (int j = 0; j < MAX_KIEL_SIZE * MAX_KIEL_SIZE; j++)
						K_linear[j] = 0.0;
				}

				// [K_L]^GG + [K_NL]^GG by local or auxiliary elements
				for (size_t j = 0; j < info->eoi[i].size(); j++)
				{
					// K_linear based on local or auxiliary elements
					Calc_K_linear_EL_conditional_for_glo(true, i, info->eoi[i][j], K_linear_conditional.data(), info);

					// K_nonlinear based on local or auxiliary elements
					if (info->c.ANALYSIS_MODE_1 == 1)
					{
						Calc_K_nonlinear_EL_conditional_for_glo(true, i, info->eoi[i][j], K_nonlinear_conditional.data(), info);
						for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
							for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
							{
								int id = k * MAX_KIEL_SIZE + l;
								double val = K_nonlinear_conditional[id];

								K_linear_conditional[id] += val;
							}
					}

					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
						{
							int id = k * MAX_KIEL_SIZE + l;
							double val = K_linear_conditional[id];

							K_linear[id] += val;
						}
				}
			}

			// [K_L]^LL + [K_NL]^LL by local elements
			else if (i >= info->Total_Element_to_mesh[1] && isAuxPatch == false)
			{
				// K_linear
				Calc_K_linear_EL(i, K_linear.data(), info);

				// K_nonlinear
				if (info->c.ANALYSIS_MODE_1 == 1)
				{
					Calc_K_nonlinear_EL(i, K_nonlinear.data(), info);
					for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
						for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
						{
							int id = j * MAX_KIEL_SIZE + k;
							double val = K_nonlinear[id];

							K_linear[id] += val;
						}
				}
			}

			else if (i >= info->Total_Element_to_mesh[1] && isAuxPatch == true)
				continue;
		}

		// K_EL add to K_Whole_Val
		if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
		{
			for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
				for (int j2 = 0; j2 < info->DIMENSION; j2++)
				{
					long row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
					if (row >= 0)
						for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
							for (int k2 = 0; k2 < info->DIMENSION; k2++)
							{
								long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
								if (col >= 0)
								{
									int id = row * K_Whole_Size + col;
									double val = K_linear[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->full_K[id] += val;
								}
							}
				}
		}
		else
		{
			for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
				for (int j2 = 0; j2 < info->DIMENSION; j2++)
				{
					long row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
					if (row >= 0)
						for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
							for (int k2 = 0; k2 < info->DIMENSION; k2++)
							{
								long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
								if (col >= row)
								{
									int id = RowCol_to_icount(row, col, info);
									double val = K_linear[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

									#pragma omp atomic
									info->K_Whole_Val[id] += val;
								}
							}
				}
		}
	}

	if (Total_mesh < 2)
	{
		cout << "end calc K (IGA)" << endl;
		return;
	}

	if (isAuxPatch)
	{
		cout << "end calc K" << endl;
		return;
	}

	// for coupling
	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<double> K_linear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		vector<double> K_nonlinear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		// coupling
		if (Total_mesh > 1 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
		{
			for (size_t j = 0; j < info->eoi[i].size(); j++)
			{
				// K_linear
				Calc_coupled_K_linear_EL(i, info->eoi[i][j], K_linear_coupled.data(), info);

				// K_nonlinear
				if (info->c.ANALYSIS_MODE_1 == 1)
				{
					Calc_coupled_K_nonlinear_EL(i, info->eoi[i][j], K_nonlinear_coupled.data(), info);
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]] * info->DIMENSION; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
						{
							int id = k * MAX_KIEL_SIZE + l;
							double val = K_nonlinear_coupled[id];

							K_linear_coupled[id] += val;
						}
				}

				// coupled_K_EL add to K_Whole_Val
				if (info->c.SOLVER == 2 || info->c.SOLVER == 3)
				{
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							long row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (col >= 0)
										{
											long id_1 = row * K_Whole_Size + col;
											long id_2 = col * K_Whole_Size + row;
											double val = K_linear_coupled[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];

											#pragma omp atomic
											info->full_K[id_1] += val;
											#pragma omp atomic
											info->full_K[id_2] += val;
										}
									}
						}
				}
				else
				{
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							long row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (col >= row)
										{
											int id = RowCol_to_icount(row, col, info);
											double val = K_linear_coupled[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];
											
											#pragma omp atomic
											info->K_Whole_Val[id] += val;
										}
									}
						}
				}
			}
		}
	}

	cout << "end calc K" << endl;
}


void Calc_K_linear_EL(int El_No, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]] * info->DIMENSION;
	vector<double> D(D_MATRIX_SIZE * D_MATRIX_SIZE);
	vector<double> B_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
	
	// information *dec_order= info->dec_order.get(); 
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	// int p = info->Element_patch[El_No];
	// double D_dil = D_elastic_dil_coef();
	vector<double> K2(MAX_KIEL_SIZE * MAX_KIEL_SIZE, 0.0);

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No].n(); i++)
	{
		// Generate linear B matrix
		double J = Make_B_Linear(El_No, &info->gp[El_No].para()[i * info->DIMENSION], info->gp[El_No].isOverlay()[i], info->gp[El_No].opp_ele()[i], &info->gp[El_No].opp_para_tilde()[i * info->DIMENSION], B_linear.data(), info);

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[El_No].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[El_No].equivalent_plastic_strain()[i], info->gp[El_No].equivalent_plastic_strain_increment()[i], El_No, i, info);
		
		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[El_No].current_stress()[i * D_MATRIX_SIZE], &info->gp[El_No].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[El_No].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// calc K_EL
		BDBJ_NL(KIEL_SIZE, B_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No].w()[i] * K1[j * MAX_KIEL_SIZE + k];

		#if 0
		// additional terms for B_bar
		if (info->c.B_BAR == 1)
		{
			double coef = D_dil * info->gp[El_No].w()[i] * J;
			shape_and_dshape(R.data(), &info->gp[El_No].para()[i * info->DIMENSION], El_No, dec_order);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int row = j * info->DIMENSION + k;
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[p]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							int col = l * info->DIMENSION + m;
							for (int n = 0; n < dec_order->No_Control_point_ON_ELEMENT[p]; n++)
							{
								int current_cp = dec_order->Controlpoint_of_Element_in_patch[El_No * MAX_NO_CP_ON_ELEMENT + n];
								// K_EL[row * MAX_KIEL_SIZE + col] += R[n] * coef * dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col]; // ほら貝
								K_EL[row * MAX_KIEL_SIZE + col] += coef * dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col]; // リーゼント
								// K2[row * MAX_KIEL_SIZE + col] += dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col];
							}
						}
				}

			// for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			// 	for (int k = 0; k < info->DIMENSION; k++)
			// 	{
			// 		int row = j * info->DIMENSION + k;
			// 		for (int n = 0; n < dec_order->No_Control_point_ON_ELEMENT[p]; n++)
			// 		{
			// 			int current_cp = dec_order->Controlpoint_of_Element_in_patch[El_No * MAX_NO_CP_ON_ELEMENT + n];
			// 			// K_EL[row * MAX_KIEL_SIZE + col] += R[n] * coef * dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col];
			// 			K_EL[row * MAX_KIEL_SIZE + col] += dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col];
			// 		}
			// 	}
		}
		#endif
	}

	#if 0
	// additional terms for B_bar
	if (info->c.B_BAR == 1)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[p]; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int row = i * info->DIMENSION + j;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int col = k * info->DIMENSION + l;
						for (int m = 0; m < dec_order->No_Control_point_ON_ELEMENT[p]; m++)
						{
							int current_cp = dec_order->Controlpoint_of_Element_in_patch[El_No * MAX_NO_CP_ON_ELEMENT + m];
							K_EL[row * MAX_KIEL_SIZE + col] += D_dil * dec_order->L2_G[p][current_cp][row] * dec_order->L2_M_inv_G[p][current_cp][col];
						}
					}
			}
	}
	#endif
}


void Calc_K_linear_EL_conditional_for_glo(bool isLocal, int El_No_glo, int El_No_loc, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;

	vector<double> D(D_MATRIX_SIZE * D_MATRIX_SIZE);
	vector<double> B_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	int target_El_No = isLocal ? El_No_loc : El_No_glo;

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[target_El_No].n(); i++)
	{
		if (info->gp[target_El_No].isSkip()[i])
			continue;

		double J = 0.0;
		if (!isLocal)
		{
			double *para = &info->gp[El_No_glo].para()[i * info->DIMENSION];

			// Generate linear B matrix
			J = Make_B_Linear(El_No_glo, para, info->gp[El_No_glo].isOverlay()[i], info->gp[El_No_glo].opp_ele()[i], &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION], B_linear.data(), info);
		}
		else if (isLocal)
		{
			if (info->gp[El_No_loc].opp_ele()[i] != El_No_glo)
				continue;
			double *para_glo = &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION];
			double *para_loc = &info->gp[El_No_loc].para()[i * info->DIMENSION];

			// Generate linear B matrix
			J = Make_updated_jacobian(El_No_loc, para_loc, info->gp[El_No_loc].isOverlay()[i], El_No_glo, para_glo, info);
			Make_B_Linear(El_No_glo, para_glo, info->gp[El_No_loc].isOverlay()[i], El_No_loc, para_loc, B_linear.data(), info);
		}

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[target_El_No].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[target_El_No].equivalent_plastic_strain()[i], info->gp[target_El_No].equivalent_plastic_strain_increment()[i], target_El_No, i, info);

		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[target_El_No].current_stress()[i * D_MATRIX_SIZE], &info->gp[target_El_No].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[target_El_No].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// calc K_EL
		double coef = info->gp[target_El_No].blending_function()[i] * info->gp[target_El_No].w()[i];
		BDBJ_NL(KIEL_SIZE, B_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += coef * K1[j * MAX_KIEL_SIZE + k];
	}
}


#if 0
void Calc_K_linear_EL_conditional_for_loc(bool isGlobal, int El_No_loc, int El_No_glo, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	vector<double> D(D_MATRIX_SIZE * D_MATRIX_SIZE);
	vector<double> B_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	int target_El_No = isGlobal ? El_No_glo : El_No_loc;

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[target_El_No].n(); i++)
	{
		if (info->gp[target_El_No].isSkip()[i])
			continue;

		double J = 0.0;
		if (!isGlobal)
		{
			double *para = &info->gp[El_No_loc].para()[i * info->DIMENSION];

			// Generate linear B matrix
			J = Make_B_Linear(El_No_loc, para, B_linear.data(), info);
		}
		else if (isGlobal)
		{
			if (!info->gp[El_No_glo].isOverlay()[i])
				continue;
			if (info->gp[El_No_glo].opp_ele()[i] != El_No_loc)
				continue;
			double *para_loc = &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION];
			double *para_glo = &info->gp[El_No_glo].para()[i * info->DIMENSION];

			// Generate linear B matrix
			J = Make_updated_jacobian(El_No_glo, para_glo, info);
			Make_B_Linear(El_No_loc, para_loc, B_linear.data(), info);
		}

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[target_El_No].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[target_El_No].equivalent_plastic_strain()[i], info->gp[target_El_No].equivalent_plastic_strain_increment()[i], target_El_No, i, info);

		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[target_El_No].current_stress()[i * D_MATRIX_SIZE], &info->gp[target_El_No].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[target_El_No].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// calc K_EL
		double coef = info->gp[target_El_No].blending_function()[i] * info->gp[target_El_No].w()[i];
		BDBJ_NL(KIEL_SIZE, B_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += coef * K1[j * MAX_KIEL_SIZE + k];
	}
}
#endif


void Calc_coupled_K_linear_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	vector<double> D(D_MATRIX_SIZE * D_MATRIX_SIZE);
	vector<double> BL_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> BG_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{
		if (!(info->gp[El_No_loc].isOverlay()[i] && (info->gp[El_No_loc].opp_ele()[i] == El_No_glo)))
			continue;

		// Generate local linear B matrix
		double J = Make_B_Linear(El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], info->gp[El_No_loc].isOverlay()[i], El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], BL_linear.data(), info);

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[El_No_loc].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[El_No_loc].equivalent_plastic_strain()[i], info->gp[El_No_loc].equivalent_plastic_strain_increment()[i], El_No_loc, i, info);

		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[El_No_loc].current_stress()[i * D_MATRIX_SIZE], &info->gp[El_No_loc].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[El_No_loc].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// Generate global linear B matrix
		Make_B_Linear(El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], info->gp[El_No_loc].isOverlay()[i], El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], BG_linear.data(), info);

		coupled_BDBJ_NL(KIEL_SIZE_glo, KIEL_SIZE_loc, BL_linear.data(), BG_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No_loc].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


#if 0
void Calc_coupled_K_linear_EL_conditional(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	vector<double> D(D_MATRIX_SIZE * D_MATRIX_SIZE);
	vector<double> BL_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> BG_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	// based on global element
	for (int i = 0; i < info->gp[El_No_glo].n(); i++)
	{
		if (!(info->gp[El_No_glo].isOverlay()[i] && (info->gp[El_No_glo].opp_ele()[i] == El_No_loc)))
			continue;

		if (info->gp[El_No_glo].isSkip()[i])
			continue;

		// Generate global linear B matrix
		double J = Make_B_Linear(El_No_glo, &info->gp[El_No_glo].para()[i * info->DIMENSION], BG_linear.data(), info);
		double coef = info->gp[El_No_glo].blending_function()[i] * info->gp[El_No_glo].w()[i];

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[El_No_glo].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[El_No_glo].equivalent_plastic_strain()[i], info->gp[El_No_glo].equivalent_plastic_strain_increment()[i], El_No_glo, i, info);

		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[El_No_glo].current_stress()[i * D_MATRIX_SIZE], &info->gp[El_No_glo].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[El_No_glo].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// Generate local linear B matrix
		Make_B_Linear(El_No_loc, &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION], BL_linear.data(), info);

		coupled_BDBJ_NL(KIEL_SIZE_glo, KIEL_SIZE_loc, BL_linear.data(), BG_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += coef * K1[j * MAX_KIEL_SIZE + k];
	}

	// based on local element
	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{
		if (!(info->gp[El_No_loc].isOverlay()[i] && (info->gp[El_No_loc].opp_ele()[i] == El_No_glo)))
			continue;

		if (info->gp[El_No_loc].isSkip()[i])
			continue;

		// Generate local linear B matrix
		double J = Make_B_Linear(El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], BL_linear.data(), info);
		double coef = info->gp[El_No_loc].blending_function()[i] * info->gp[El_No_loc].w()[i];

		// Generate D matrix
		D_elastic_smart(D.data(), info);
		if (info->gp[El_No_loc].equivalent_plastic_strain_increment()[i] > 0.0)
			D_elastoplastic(D.data(), info->gp[El_No_loc].equivalent_plastic_strain()[i], info->gp[El_No_loc].equivalent_plastic_strain_increment()[i], El_No_loc, i, info);

		// Modify D matrix with finite strains
		if (info->c.ANALYSIS_MODE_1 == 1)
			Modify_D_matrix(D.data(), &info->gp[El_No_loc].current_stress()[i * D_MATRIX_SIZE], &info->gp[El_No_loc].elastic_strain_trial()[i * D_MATRIX_SIZE], &info->gp[El_No_loc].current_deformation_gradient()[i * info->DIMENSION * info->DIMENSION]);

		// Generate global linear B matrix
		Make_B_Linear(El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], BG_linear.data(), info);

		coupled_BDBJ_NL(KIEL_SIZE_glo, KIEL_SIZE_loc, BL_linear.data(), BG_linear.data(), D.data(), J, K1.data());
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += coef * K1[j * MAX_KIEL_SIZE + k];
	}
}
#endif


void Calc_K_nonlinear_EL(int El_No, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]] * info->DIMENSION;
	constexpr int S_size = 9;

	vector<double> Sigma(S_size * S_size);
	vector<double> B_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BTS(S_size * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No].n(); i++)
	{
		// Generate nonlinear B matrix
		double J = Make_B_Nonlinear(El_No, &info->gp[El_No].para()[i * info->DIMENSION], info->gp[El_No].isOverlay()[i], info->gp[El_No].opp_ele()[i], &info->gp[El_No].opp_para_tilde()[i * info->DIMENSION], B_nonlinear.data(), info);

		// Generate S matrix
		generate_S_matrix(Sigma.data(), &info->gp[El_No].current_stress()[i * D_MATRIX_SIZE]);

		// Calculate integral of [B_NL]^T * [S] * [B_NL] * dv
		for (int j = 0; j < S_size * MAX_KIEL_SIZE; j++)
			BTS[j] = 0.0;
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K1[j * MAX_KIEL_SIZE + k] = 0.0;

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < S_size; k++)
					BTS[j * S_size + k] += B_nonlinear[l * MAX_KIEL_SIZE + j] * Sigma[l * S_size + k];

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < KIEL_SIZE; k++)
					K1[j * MAX_KIEL_SIZE + k] += BTS[j * S_size + l] * B_nonlinear[l * MAX_KIEL_SIZE + k];

		double temp_val = J * info->gp[El_No].w()[i];
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += K1[j * MAX_KIEL_SIZE + k] * temp_val;
	}
}


void Calc_K_nonlinear_EL_conditional_for_glo(bool isLocal, int El_No_glo, int El_No_loc, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	constexpr int S_size = 9;

	vector<double> Sigma(S_size * S_size);
	vector<double> B_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BTS(S_size * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	int target_El_No = isLocal ? El_No_loc : El_No_glo;

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[target_El_No].n(); i++)
	{
		if (info->gp[target_El_No].isSkip()[i])
			continue;

		double J = 0.0;
		if (!isLocal)
		{
			double *para = &info->gp[El_No_glo].para()[i * info->DIMENSION];

			// Generate nonlinear B matrix
			J = Make_B_Nonlinear(El_No_glo, para, info->gp[El_No_glo].isOverlay()[i], info->gp[El_No_glo].opp_ele()[i], &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION], B_nonlinear.data(), info);
		}
		else if (isLocal)
		{
			if (info->gp[El_No_loc].opp_ele()[i] != El_No_glo)
				continue;
			double *para_glo = &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION];
			double *para_loc = &info->gp[El_No_loc].para()[i * info->DIMENSION];

			// Generate nonlinear B matrix
			J = Make_updated_jacobian(El_No_loc, para_loc, info->gp[El_No_loc].isOverlay()[i], El_No_glo, para_glo, info);
			Make_B_Nonlinear(El_No_glo, para_glo, info->gp[El_No_loc].isOverlay()[i], El_No_loc, para_loc, B_nonlinear.data(), info);
		}

		// Generate S matrix
		generate_S_matrix(Sigma.data(), &info->gp[target_El_No].current_stress()[i * D_MATRIX_SIZE]);

		// Calculate integral of [B_NL]^T * [S] * [B_NL] * dv
		for (int j = 0; j < S_size * MAX_KIEL_SIZE; j++)
			BTS[j] = 0.0;
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K1[j * MAX_KIEL_SIZE + k] = 0.0;

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < S_size; k++)
					BTS[j * S_size + k] += B_nonlinear[l * MAX_KIEL_SIZE + j] * Sigma[l * S_size + k];

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < KIEL_SIZE; k++)
					K1[j * MAX_KIEL_SIZE + k] += BTS[j * S_size + l] * B_nonlinear[l * MAX_KIEL_SIZE + k];

		double temp_val = J * info->gp[target_El_No].w()[i] * info->gp[target_El_No].blending_function()[i];
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += K1[j * MAX_KIEL_SIZE + k] * temp_val;
	}
}


#if 0
void Calc_K_nonlinear_EL_conditional_for_loc(bool isGlobal, int El_No_loc, int El_No_glo, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;
	constexpr int S_size = 9;

	vector<double> Sigma(S_size * S_size);
	vector<double> B_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BTS(S_size * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	int target_El_No = isGlobal ? El_No_glo : El_No_loc;

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[target_El_No].n(); i++)
	{
		if (info->gp[target_El_No].isSkip()[i])
			continue;

		double J = 0.0;
		if (!isGlobal)
		{
			double *para = &info->gp[El_No_loc].para()[i * info->DIMENSION];

			// Generate nonlinear B matrix
			J = Make_B_Nonlinear(El_No_loc, para, B_nonlinear.data(), info);
		}
		else if (isGlobal)
		{
			if (!info->gp[El_No_glo].isOverlay()[i])
				continue;
			if (info->gp[El_No_glo].opp_ele()[i] != El_No_loc)
				continue;
			double *para_loc = &info->gp[El_No_glo].opp_para_tilde()[i * info->DIMENSION];
			double *para_glo = &info->gp[El_No_glo].para()[i * info->DIMENSION];

			// Generate nonlinear B matrix
			J = Make_updated_jacobian(El_No_glo, para_glo, info);
			Make_B_Nonlinear(El_No_loc, para_loc, B_nonlinear.data(), info);
		}

		// Generate S matrix
		generate_S_matrix(Sigma.data(), &info->gp[target_El_No].current_stress()[i * D_MATRIX_SIZE]);

		// Calculate integral of [B_NL]^T * [S] * [B_NL] * dv
		for (int j = 0; j < S_size * MAX_KIEL_SIZE; j++)
			BTS[j] = 0.0;
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K1[j * MAX_KIEL_SIZE + k] = 0.0;

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < S_size; k++)
					BTS[j * S_size + k] += B_nonlinear[l * MAX_KIEL_SIZE + j] * Sigma[l * S_size + k];

		for (int j = 0; j < KIEL_SIZE; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < KIEL_SIZE; k++)
					K1[j * MAX_KIEL_SIZE + k] += BTS[j * S_size + l] * B_nonlinear[l * MAX_KIEL_SIZE + k];

		double temp_val = J * info->gp[target_El_No].w()[i] * info->gp[target_El_No].blending_function()[i];
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += K1[j * MAX_KIEL_SIZE + k] * temp_val;
	}
}
#endif


void Calc_coupled_K_nonlinear_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;
	constexpr int S_size = 9;

	vector<double> Sigma(S_size * S_size);
	vector<double> BL_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BG_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BTS(S_size * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{	
		if (!(info->gp[El_No_loc].isOverlay()[i] && (info->gp[El_No_loc].opp_ele()[i] == El_No_glo)))
			continue;
		
		// Generate local nonlinear B matrix
		double J = Make_B_Nonlinear(El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], info->gp[El_No_loc].isOverlay()[i], El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], BL_nonlinear.data(), info);

		// Generate S matrix
		generate_S_matrix(Sigma.data(), &info->gp[El_No_loc].current_stress()[i * D_MATRIX_SIZE]);

		// Generate global nonlinear B matrix
		Make_B_Nonlinear(El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], info->gp[El_No_loc].isOverlay()[i], El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], BG_nonlinear.data(), info);

		// Calculate integral of [B_NL]^T * [S] * [B_NL] * dv
		for (int j = 0; j < S_size * MAX_KIEL_SIZE; j++)
			BTS[j] = 0.0;
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				K1[j * MAX_KIEL_SIZE + k] = 0.0;

		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < S_size; k++)
					BTS[j * S_size + k] += BG_nonlinear[l * MAX_KIEL_SIZE + j] * Sigma[l * S_size + k];

		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < KIEL_SIZE_loc; k++)
					K1[j * MAX_KIEL_SIZE + k] += BTS[j * S_size + l] * BL_nonlinear[l * MAX_KIEL_SIZE + k];

		double temp_val = J * info->gp[El_No_loc].w()[i];
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += K1[j * MAX_KIEL_SIZE + k] * temp_val;
	}
}


#if 0
/* 未着手 */
void Calc_coupled_K_nonlinear_EL_conditional(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;
	constexpr int S_size = 9;

	vector<double> Sigma(S_size * S_size);
	vector<double> BL_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BG_nonlinear(S_size * MAX_KIEL_SIZE);
	vector<double> BTS(S_size * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{	
		if (!(info->gp[El_No_loc].isOverlay()[i] && (info->gp[El_No_loc].opp_ele()[i] == El_No_glo)))
			continue;
		
		// Generate local nonlinear B matrix
		double J = Make_B_Nonlinear(El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], BL_nonlinear.data(), info);

		// Generate S matrix
		generate_S_matrix(Sigma.data(), &info->gp[El_No_loc].current_stress()[i * D_MATRIX_SIZE]);

		// Generate global nonlinear B matrix
		Make_B_Nonlinear(El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], BG_nonlinear.data(), info);

		// Calculate integral of [B_NL]^T * [S] * [B_NL] * dv
		for (int j = 0; j < S_size * MAX_KIEL_SIZE; j++)
			BTS[j] = 0.0;
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				K1[j * MAX_KIEL_SIZE + k] = 0.0;

		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < S_size; k++)
					BTS[j * S_size + k] += BG_nonlinear[l * MAX_KIEL_SIZE + j] * Sigma[l * S_size + k];

		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int l = 0; l < S_size; l++)
				for (int k = 0; k < KIEL_SIZE_loc; k++)
					K1[j * MAX_KIEL_SIZE + k] += BTS[j * S_size + l] * BL_nonlinear[l * MAX_KIEL_SIZE + k];

		double temp_val = J * info->gp[El_No_loc].w()[i];
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += K1[j * MAX_KIEL_SIZE + k] * temp_val;
	}
}
#endif


void generate_S_matrix(double *Sigma, const double *current_stress)
{
	// init
	for (int i = 0; i < 9 * 9; i++)
		Sigma[i] = 0.0;

	// Substitute current stresses for S matrix
	for (int i = 0; i < 3; i++)
	{
		Sigma[i * 9 + i]     = current_stress[0];
		Sigma[i * 9 + i + 3] = current_stress[3];
		Sigma[i * 9 + i + 6] = current_stress[5];

		Sigma[(i + 3) * 9 + i]     = current_stress[3];
		Sigma[(i + 3) * 9 + i + 3] = current_stress[1];
		Sigma[(i + 3) * 9 + i + 6] = current_stress[4];

		Sigma[(i + 6) * 9 + i]     = current_stress[5];
		Sigma[(i + 6) * 9 + i + 3] = current_stress[4];
		Sigma[(i + 6) * 9 + i + 6] = current_stress[2];
	}
}


void BDBJ_NL(int KIEL_SIZE, double *B, double *D, double J, double *K_EL)
{
	vector<double> BD(MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (int i = 0; i < KIEL_SIZE * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	// [B]T[D][B]Jの計算
	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
}


void coupled_BDBJ_NL(int KIEL_SIZE_glo, int KIEL_SIZE_loc, double *B, double *BG, double *D, double J, double *K_EL)
{
	vector<double> BD(MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	//[B]GT[D][B]LJの計算
	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
}


void Update_disp_increment(double delta_factor, information *info)
{
	for (int i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		int id = info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1];
		info->disp_increment[id] = info->Value_of_Constraint[i] * delta_factor;
		// info->disp_increment[id] = info->Value_of_Constraint[i] / (double)(info->c.STEP_N);
	}

	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
		{
			int index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
			{
				// info->disp_increment[i * info->DIMENSION + j] = info->sol_vec[index];
				info->disp_increment[i * info->DIMENSION + j] += info->sol_vec[index];
			}
		}

	// auxiliary patch
	if (Total_mesh >= 2 && info->c.BLENDING == 1)
	{
		information *aux = info->aux.get();
		for (int i = 0; i < aux->Total_Constraint_to_mesh[1]; i++)
		{
			int id = aux->Constraint_Node_Dir[i * 2 + 0] * aux->DIMENSION + aux->Constraint_Node_Dir[i * 2 + 1];
			aux->disp_increment[id] = aux->Value_of_Constraint[i] * delta_factor;
			// aux->disp_increment[id] = aux->Value_of_Constraint[i] / (double)(aux->c.STEP_N);
		}

		// global patch
		for (int i = 0; i < aux->Total_Control_Point_to_mesh[1]; i++)
			for (int j = 0; j < aux->DIMENSION; j++)
			{
				int index = aux->Index_Dof[i * aux->DIMENSION + j];
				if (index >= 0)
				{
					// aux->disp_increment[i * aux->DIMENSION + j] = info->sol_vec[index];
					aux->disp_increment[i * aux->DIMENSION + j] += info->sol_vec[index];
				}
			}

		// auxiliary patch
		for (int i = aux->Total_Control_Point_to_mesh[1]; i < aux->Total_Control_Point_to_mesh[Total_mesh]; i++)
			for (int j = 0; j < aux->DIMENSION; j++)
				aux->disp_increment[i * aux->DIMENSION + j] = 0.0;
	}
}


// calculate B
double Make_B_component_small_strain(int ele, double *para, double *out_B_component, information *info)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];

	static thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k * 2 + i] * dR[j * info->DIMENSION + k];
			}
		}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k * 3 + i] * dR[j * info->DIMENSION + k];
			}
		}

		return J;
	}

	return J;
}


double Make_B_component(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, double *out_B_component, information *info)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];

	static thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k * 2 + i] * dR[j * info->DIMENSION + k];
			}
		}

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		if (Total_mesh < 2)
		{
			for (int i = 0; i < info->DIMENSION; i++)
			{
				for (int j = 0; j < info->DIMENSION; j++)
				{
					a_3x3[i * 3 + j] = 0.0;
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
					{
						int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
						int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
	
						if (info->c.ANALYSIS_MODE_1 == 0)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
						else if (info->c.ANALYSIS_MODE_1 == 1)
						{
							if (Total_mesh < 2)
								a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
							else if (Total_mesh >= 2)
								a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
						}
	
						// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
					}
				}
			}
	
			J = InverseMatrix_3x3(a_3x3);
	
			for (int i = 0; i < info->DIMENSION; i++)
			{
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
				{
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
					for (int k = 0; k < info->DIMENSION; k++)
						out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k * 3 + i] * dR[j * info->DIMENSION + k];
				}
			}
		}
		else
		{
			//***********************
			// Modified on 2024/06/15
			//***********************
			// ヤコビアンの計算方法を変更
			// 形状アップデートを行わずに，∂N/∂ξ^G=∑[∂N/∂ξ^G][J^G]^(-1), ∂N/∂ξ^L=∑[∂N/∂ξ^L][J^L]^(-1)を計算する
			// J^G=∂x/∂ξ^G=∂(X^G+u^G)/∂ξ^G+(∂u^L/∂ξ^L)(∂ξ^L/∂ξ^G), J^L=∂x/∂ξ^L=∂(X^L+u^L)/∂ξ^L+(∂u^G/∂ξ^G)(∂ξ^G/∂ξ^L)
			// ここで，(∂ξ^L/∂ξ^G)および(∂ξ^G/∂ξ^L)はそれぞれ，マッピングテンソルTとして計算する．
			// Tは初期配置を参照して，(∂ξ^L/∂ξ^G)=[J^L_0]^(-1)[J^G_0]，(∂ξ^G/∂ξ^L)=[J^G_0]^(-1)[J^L_0]で計算する．
			//***********************
	
			//    | global patch    | local patch
			// ---+-----------------+----------------
			// a: | ∂(X^G+u^G)/∂ξ^G | ∂(X^L+u^L)/∂ξ^L
			// b: | ∂(u^L)/∂ξ^L     | ∂(u^G)/∂ξ^G
			// T: | ∂ξ^L/∂ξ^G       | ∂ξ^G/∂ξ^L
	
			// Total_Jac = a + b * T
			// B_component = Total_Jac^(-1) * dR^{G or L}
	
			double *Total_Jac = a_3x3;
			double J_initial_coord[9] = {0.0};
	
			for (int i = 0; i < info->DIMENSION; i++)
				for (int j = 0; j < info->DIMENSION; j++)
				{
					a_3x3[i * 3 + j] = 0.0;
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
					{
						int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
						int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
	
						if (info->c.ANALYSIS_MODE_1 == 0)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
						else if (info->c.ANALYSIS_MODE_1 == 1)
						{
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
							J_initial_coord[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
						}
					}
				}
	
			// overlay
			if (Total_mesh >= 2 && info->c.ANALYSIS_MODE_1 == 1 && isOverlay)
			{
				static thread_local vector<double> R_opp(MAX_NO_CP_ON_ELEMENT);
				static thread_local vector<double> dR_opp(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	
				shape_and_dshape(R_opp.data(), dR_opp.data(), opp_para, opp_ele, true, info);
	
				// calc b_3x3 and opp_J_initial_coord
				double b_3x3[9] = {0.0};
				double opp_J_initial_coord[9] = {0.0};
				for (int i = 0; i < info->DIMENSION; i++)
					for (int j = 0; j < info->DIMENSION; j++)
					{
						b_3x3[i * 3 + j] = 0.0;
						opp_J_initial_coord[i * 3 + j] = 0.0;
						for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[opp_ele]]; k++)
						{
							int id0 = info->Controlpoint_of_Element[opp_ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
							int id1 = info->Controlpoint_of_Element[opp_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
							b_3x3[i * 3 + j] += dR_opp[k * info->DIMENSION + j] * (info->disp[id1] + info->disp_increment[id1]);
							opp_J_initial_coord[i * 3 + j] += dR_opp[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
						}
					}
	
				// calc T
				double T[9] = {0.0};
				InverseMatrix_3x3(opp_J_initial_coord);
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int k = 0; k < 3; k++)
						{
							T[i * 3 + k] += opp_J_initial_coord[i * 3 + j] * J_initial_coord[j * 3 + k];
							// T[i * 3 + k] += opp_J_initial_coord[j * 3 + i] * J_initial_coord[j * 3 + k];
						}
	
				double b_T[9] = {0.0};
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int k = 0; k < 3; k++)
							b_T[i * 3 + k] += b_3x3[i * 3 + j] * T[j * 3 + k];
	
				for (int i = 0; i < 9; i++)
					Total_Jac[i] += b_T[i];
			}
	
			J = InverseMatrix_3x3(Total_Jac);
	
			for (int i = 0; i < info->DIMENSION; i++)
			{
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
				{
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
					for (int k = 0; k < info->DIMENSION; k++)
						out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += Total_Jac[k * 3 + i] * dR[j * info->DIMENSION + k];
				}
			}
		}

		return J;
	}

	return J;
}


double Make_B_component(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, double *out_B_component, vector<double> &R, vector<double> &dR, information *info)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];

	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_2x2(a_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k * 2 + i] * dR[j * info->DIMENSION + k];
			}
		}

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		#if 0
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_3x3(a_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k * 3 + i] * dR[j * info->DIMENSION + k];
			}
		}
		#else

		//***********************
		// Modified on 2024/06/15
		//***********************
		// ヤコビアンの計算方法を変更
		// 形状アップデートを行わずに，∂N/∂ξ^G=∑[∂N/∂ξ^G][J^G]^(-1), ∂N/∂ξ^L=∑[∂N/∂ξ^L][J^L]^(-1)を計算する
		// J^G=∂x/∂ξ^G=∂(X^G+u^G)/∂ξ^G+(∂u^L/∂ξ^L)(∂ξ^L/∂ξ^G), J^L=∂x/∂ξ^L=∂(X^L+u^L)/∂ξ^L+(∂u^G/∂ξ^G)(∂ξ^G/∂ξ^L)
		// ここで，(∂ξ^L/∂ξ^G)および(∂ξ^G/∂ξ^L)はそれぞれ，マッピングテンソルTとして計算する．
		// Tは初期配置を参照して，(∂ξ^L/∂ξ^G)=[J^L_0]^(-1)[J^G_0]，(∂ξ^G/∂ξ^L)=[J^G_0]^(-1)[J^L_0]で計算する．
		//***********************

		//    | global patch    | local patch
		// ---+-----------------+----------------
		// a: | ∂(X^G+u^G)/∂ξ^G | ∂(X^L+u^L)/∂ξ^L
		// b: | ∂(u^L)/∂ξ^L     | ∂(u^G)/∂ξ^G
		// T: | ∂ξ^L/∂ξ^G       | ∂ξ^G/∂ξ^L

		// Total_Jac = a + b * T
		// B_component = Total_Jac^(-1) * dR^{G or L}

		double *Total_Jac = a_3x3;
		double J_initial_coord[9] = {0.0};

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						J_initial_coord[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					}
				}
			}

		// overlay
		if (Total_mesh >= 2 && info->c.ANALYSIS_MODE_1 == 1 && isOverlay)
		{
			static thread_local vector<double> R_opp(MAX_NO_CP_ON_ELEMENT);
			static thread_local vector<double> dR_opp(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

			shape_and_dshape(R_opp.data(), dR_opp.data(), opp_para, opp_ele, true, info);

			// calc b_3x3 and opp_J_initial_coord
			double b_3x3[9] = {0.0};
			double opp_J_initial_coord[9] = {0.0};
			for (int i = 0; i < info->DIMENSION; i++)
				for (int j = 0; j < info->DIMENSION; j++)
				{
					b_3x3[i * 3 + j] = 0.0;
					opp_J_initial_coord[i * 3 + j] = 0.0;
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[opp_ele]]; k++)
					{
						int id0 = info->Controlpoint_of_Element[opp_ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
						int id1 = info->Controlpoint_of_Element[opp_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;
						b_3x3[i * 3 + j] += dR_opp[k * info->DIMENSION + j] * (info->disp[id1] + info->disp_increment[id1]);
						opp_J_initial_coord[i * 3 + j] += dR_opp[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					}
				}

			// calc T
			double T[9] = {0.0};
			InverseMatrix_3x3(opp_J_initial_coord);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						T[i * 3 + k] += opp_J_initial_coord[j * 3 + i] * J_initial_coord[j * 3 + k];

			double b_T[9] = {0.0};
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						b_T[i * 3 + k] += b_3x3[i * 3 + j] * T[j * 3 + k];

			for (int i = 0; i < 9; i++)
				Total_Jac[i] += b_T[i];
		}

		J = InverseMatrix_3x3(Total_Jac);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += Total_Jac[k * 3 + i] * dR[j * info->DIMENSION + k];
			}
		}
		#endif

		return J;
	}

	return J;
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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

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
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * coef[j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				a[i * info->DIMENSION + j] = a_3x3[i * info->DIMENSION + j];
	}
}


double Make_updated_jacobian(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component(ele, para, isOverlay, opp_ele, opp_para, b.data(), info);
	return J;
}


double Make_B_L_NL(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, double *out_BL, double *out_BNL, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component(ele, para, isOverlay, opp_ele, opp_para, b.data(), info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];

			// NonLinear
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];

			// Nonlinear
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
		}

		// B-bar method
		if (info->c.B_BAR == 1)
		{
			#if 0
			double *b_ave = &info->b_average[ele * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];

			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}

			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}
			#else
			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}

			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}
			#endif
		}

		return J;
	}
	return J;
}


double Make_B_Linear(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, double *out_BL, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component(ele, para, isOverlay, opp_ele, opp_para, b.data(), info);
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		
		// B-bar method
		if (info->c.B_BAR == 1)
		{
			#if 0
			double *b_ave = &info->b_average[ele * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];

			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];

				// out_BL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				// out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[0];
				// out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[0];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i]     += diff[1];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] += diff[1];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i]     += diff[2];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] += diff[2];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}
			#else
			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}
			#endif
		}

		return J;
	}
	return J;
}


double Make_B_Nonlinear(int ele, double *para, bool isOverlay, int opp_ele, double *opp_para, double *out_BNL, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component(ele, para, isOverlay, opp_ele, opp_para, b.data(), info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// NonLinear
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Nonlinear
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
		}

		// B-bar method
		if (info->c.B_BAR == 1)
		{
			#if 0
			double *b_ave = &info->b_average[ele * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];

			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];

				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[0];
				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[0];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     += diff[1];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] += diff[1];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     += diff[2];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] += diff[2];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}
			#else
			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}
			#endif
		}

		return J;
	}
	return J;
}


double Make_B_Linear_small_strain(int ele, double *para, double *out_BL, information *info)
{
	// vector<double> debug_point = {-0.861136, -0.993057, -0.861136};
	// double eps = 1e-6;
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component_small_strain(ele, para, b.data(), info);

	// if (ele == 8 && abs(para[0] - debug_point[0]) < eps && abs(para[1] - debug_point[1]) < eps && abs(para[2] - debug_point[2]) < eps)
	// {
	// 	cout << "b" << endl;
	// 	for (int i = 0; i < info->DIMENSION * MAX_NO_CP_ON_ELEMENT; i++)
	// 		cout << b[i] << " ";
	// 	exit(0);
	// }

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		
		// B-bar method
		if (info->c.B_BAR == 1)
		{
			#if 0
			double *b_ave = &info->b_average[ele * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];

			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];

				// out_BL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				// out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[0];
				// out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[0];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i]     += diff[1];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				// out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] += diff[1];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i]     += diff[2];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] += diff[2];
				// out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}
			#else
			// defference linear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}
			#endif
		}

		// if (ele == 8 && abs(para[0] - debug_point[0]) < eps && abs(para[1] - debug_point[1]) < eps && abs(para[2] - debug_point[2]) < eps)
		// {
		// 	cout << "out_BL_old" << endl;
		// 	for (int i = 0; i < D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT; i++)
		// 		cout << out_BL[i] << " ";
		// 	exit(0);
		// }

		return J;
	}
	return J;
}


double Make_B_Nonlinear_small_strain(int ele, double *para, double *out_BNL, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component_small_strain(ele, para, b.data(), info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// NonLinear
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Nonlinear
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
		}

		// B-bar method
		if (info->c.B_BAR == 1)
		{
			#if 0
			double *b_ave = &info->b_average[ele * info->DIMENSION * MAX_NO_CP_ON_ELEMENT];

			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double diff[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					diff[j] = factor * (b_ave[j * MAX_NO_CP_ON_ELEMENT + i] - b[j * MAX_NO_CP_ON_ELEMENT + i]);

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];

				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     += diff[0];
				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] += diff[0];
				// out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] += diff[0];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     += diff[1];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] += diff[1];
				// out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] += diff[1];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     += diff[2];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] += diff[2];
				// out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] += diff[2];
			}
			#else
			// defference nonlinear
			for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
			{
				constexpr double factor = 1.0 / 3.0;
				double dilatational_term[3] = {0.0};
				for (int j = 0; j < info->DIMENSION; j++)
					dilatational_term[j] = factor * b[j * MAX_NO_CP_ON_ELEMENT + i];

				out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     -= dilatational_term[0];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] -= dilatational_term[1];
				out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] -= dilatational_term[2];
			}
			#endif
		}

		return J;
	}
	return J;
}


#if 0
void Make_dBduo_component(int ele, double *para, vector<vector<double>> &out_dBduo_component, information *info)
{
	int p = info->Element_patch[ele];

	double a_inv_2x2[4], a_inv_3x3[9];

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_inv_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_inv_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_inv_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_inv_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_inv_2x2[i * 2 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_2x2(a_inv_2x2);
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_inv_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id0 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					int id1 = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + i;

					if (info->c.ANALYSIS_MODE_1 == 0)
						a_inv_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id0];
					else if (info->c.ANALYSIS_MODE_1 == 1)
					{
						if (Total_mesh < 2)
							a_inv_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
						else if (Total_mesh >= 2)
							a_inv_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * (info->Node_Coordinate[id0] + info->disp_overlay[id1] + info->disp_overlay_increment[id1]);
					}

					// a_inv_3x3[i * 3 + j] += dShape_func(k, j, para, ele, info) * (info->Node_Coordinate[id0] + info->disp[id1] + info->disp_increment[id1]);
				}
			}
		}

		J = InverseMatrix_3x3(a_inv_3x3);
	}

	// [dB/duo] = [dN/dxi] * (- [A]^-1 * [dA/duo] * [A]^-1), [A] = [dxi/dx(uo)]
	vector<vector<double>> dA_duo(info->DIMENSION * info->DIMENSION, vector<double>(info->No_Control_point_ON_ELEMENT[p] * info->DIMENSION, 0.0));
	vector<vector<double>> temp_tensor_1(info->DIMENSION * info->DIMENSION, vector<double>(info->No_Control_point_ON_ELEMENT[p] * info->DIMENSION), 0.0);
	vector<vector<double>> temp_tensor_2(info->DIMENSION * info->DIMENSION, vector<double>(info->No_Control_point_ON_ELEMENT[p] * info->DIMENSION), 0.0);

	// [dA/duo]
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[p]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			dA_duo[j * info->DIMENSION + j][i * info->DIMENSION + j] = dR[i * info->DIMENSION + j];
	
	// (- [A]^-1 * [dA/duo] * [A]^-1)
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[p]; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				for (int k = 0; k < info->DIMENSION; k++)
					for (int m = 0; m < info->DIMENSION; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							temp_tensor_1[k * info->DIMENSION + l][i * info->DIMENSION + j] -= a_inv_2x2[k * 2 + m] * dA_duo[m * info->DIMENSION + l][i * info->DIMENSION + j];
				for (int k = 0; k < info->DIMENSION; k++)
					for (int m = 0; m < info->DIMENSION; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							temp_tensor_2[k * info->DIMENSION + l][i * info->DIMENSION + j] += temp_tensor_1[k * info->DIMENSION + m][i * info->DIMENSION + j] * a_inv_2x2[m * 2 + l];
			}
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[p]; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				for (int k = 0; k < info->DIMENSION; k++)
					for (int m = 0; m < info->DIMENSION; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							temp_tensor_1[k * info->DIMENSION + l][i * info->DIMENSION + j] -= a_inv_3x3[k * 3 + m] * dA_duo[m * info->DIMENSION + l][i * info->DIMENSION + j];
				for (int k = 0; k < info->DIMENSION; k++)
					for (int m = 0; m < info->DIMENSION; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							temp_tensor_2[k * info->DIMENSION + l][i * info->DIMENSION + j] += temp_tensor_1[k * info->DIMENSION + m][i * info->DIMENSION + j] * a_inv_3x3[m * 3 + l];
			}
	}

	// [dN/dxi] * (- [A]^-1 * [dA/duo] * [A]^-1)
	for (int i = 0; i < info->DIMENSION; i++)
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
					for (int m = 0; m < info->DIMENSION; m++)
						out_dBduo_component[i * MAX_NO_CP_ON_ELEMENT + j][k * info->DIMENSION + l] += dR[]
}


void Make_dBduo_L_NL(int ele, double *para, vector<vector<double>> &out_dBduo_L, vector<vector<double>> &out_dBduo_NL, information *info)
{
	vector<vector<double>> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, vector<double>(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0));
	double J = Make_dBduo_component(ele, para, b, info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];

			// NonLinear
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[2 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];

			// Nonlinear
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[4 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[6 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BNL[7 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BNL[8 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
		}
	}
}
#endif


// Internal Force
void Update_internal_force(vector<bool> &ep_flag, information *info)
{
	// auxilary patch
	bool isAuxPatch = false;
	if (Total_mesh > 1 &&  info->c.BLENDING == 1 && info->aux == nullptr)
		isAuxPatch = true;
	cout << "isAuxPatch: " << isAuxPatch << endl;

	// init
	if (isAuxPatch == false)
	{
		for (int i = 0; i < K_Whole_Size; i++)
			info->internal_force[i] = 0.0;
	}

	// information *dec_order= info->dec_order.get();

	// Generate elastic D matrix
	vector<double> D_e(D_MATRIX_SIZE * D_MATRIX_SIZE);
	D_elastic(D_e.data(), info);

	long elastic_counter = 0;
	long plastic_counter = 0;
	long unloading_counter = 0;

	gp_switch(true, info);

	#pragma omp parallel for schedule(dynamic) reduction(+:elastic_counter, plastic_counter, unloading_counter) collapse(1)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		if (isAuxPatch && i < info->Total_Element_to_mesh[1])
			continue;

		double temp = 0.0;

		thread_local double *dF = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *B_linear = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
		thread_local double *elastic_strain_tensor = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *elastic_strain_tensor_trial = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *B_e = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *B_e_trial = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *temp_matrix = (double *)malloc(sizeof(double) * info->DIMENSION * info->DIMENSION);
		thread_local double *relative_stress_trial = (double *)malloc(sizeof(double) * D_MATRIX_SIZE);
		thread_local double *B_linear_overlay = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);

		for (int j = 0; j < info->gp[i].n(); j++)
		{
			bool overlay = false;
			int overlay_ele = -1;

			// search element
			if (Total_mesh > 1)
			{
				overlay = info->gp[i].isOverlay()[j];
				if (overlay)
				{
					overlay_ele = info->gp[i].opp_ele()[j];
					Make_B_Linear(overlay_ele, &info->gp[i].opp_para_tilde()[j * info->DIMENSION], overlay, i, &info->gp[i].para()[j * info->DIMENSION], B_linear_overlay, info);
				}
			}

			// Generate linear B matrix
			double J = Make_B_Linear(i, &info->gp[i].para()[j * info->DIMENSION], info->gp[i].isOverlay()[j], info->gp[i].opp_ele()[j], &info->gp[i].opp_para_tilde()[j * info->DIMENSION], B_linear, info);

			// 有限変形
			if (info->c.ANALYSIS_MODE_1 == 1)
			{
				/* debug 1 */
				{
					// Calculate relative deformation gradients: [dF] = ([I] - [d(du)/d(x + u + du)])^-1
					identify_tensor(info->DIMENSION, dF);
					for (int k = 0; k < info->DIMENSION; k++)
						for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
							for (int l = 0; l < info->DIMENSION; l++)
								dF[k * info->DIMENSION + l] -= B_linear[l * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + k];
								// dF[k * info->DIMENSION + l] -= B_linear[l * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_overlay_increment[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + k];
					if (overlay)
					{
						for (int k = 0; k < info->DIMENSION; k++)
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; m++)
								for (int l = 0; l < info->DIMENSION; l++)
									dF[k * info->DIMENSION + l] -= B_linear_overlay[l * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + k];
									// dF[k * info->DIMENSION + l] -= B_linear_overlay[l * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_overlay_increment[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + k];
					}

					if (info->DIMENSION == 2)
						InverseMatrix_2x2(dF);
					else if (info->DIMENSION == 3)
						InverseMatrix_3x3(dF);
				}

				// Calculate deformation gradients: [F] = [dF] * [F]
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						info->gp[i].current_deformation_gradient()[j * info->DIMENSION * info->DIMENSION + k * info->DIMENSION + l] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < info->DIMENSION; m++)
							info->gp[i].current_deformation_gradient()[j * info->DIMENSION * info->DIMENSION + k * info->DIMENSION + m] += dF[k * info->DIMENSION + l] * info->gp[i].deformation_gradient()[j * info->DIMENSION * info->DIMENSION + l * info->DIMENSION + m];

				// Calculate elastic left Cauchy-Green deformations: [B]^e = exp(2 * {epsilon}^e)
				if (info->DIMENSION == 2)
				{
					// voigt to tensor
					elastic_strain_tensor[0] = 2.0 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 0];
					elastic_strain_tensor[1] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 2];
					elastic_strain_tensor[2] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 2];
					elastic_strain_tensor[3] = 2.0 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 1];
				}
				else if (info->DIMENSION == 3)
				{
					// voigt to tensor
					elastic_strain_tensor[0] = 2.0 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 0];
					elastic_strain_tensor[1] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 3];
					elastic_strain_tensor[2] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 5];
					elastic_strain_tensor[3] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 3];
					elastic_strain_tensor[4] = 2.0 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 1];
					elastic_strain_tensor[5] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 4];
					elastic_strain_tensor[6] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 5];
					elastic_strain_tensor[7] = 2.0 * 0.5 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 4];
					elastic_strain_tensor[8] = 2.0 * info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + 2];
				}
				tensor_exponential(info->DIMENSION, elastic_strain_tensor, B_e);

				// Calculate trial elastic left Cauchy-Green deformations: [B]^trial = [dF] * [B]^e * [dF]^T
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double a = 0.0;
						for (int m = 0; m < info->DIMENSION; m++)
							a += dF[k * info->DIMENSION + m] * B_e[m * info->DIMENSION + l];
						temp_matrix[k * info->DIMENSION + l] = a;
					}
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double a = 0.0;
						for (int m = 0; m < info->DIMENSION; m++)
							a += temp_matrix[k * info->DIMENSION + m] * dF[l * info->DIMENSION + m];
						B_e_trial[k * info->DIMENSION + l] = a;
					}

				// Calculate trial elastic strains: {epsilon}^trial = 1/2 * ln([B]^trial)
				tensor_logarithm(info->DIMENSION, B_e_trial, elastic_strain_tensor_trial);
				if (info->DIMENSION == 2)
				{
					// tensor to voigt
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 0] = 0.5 * elastic_strain_tensor_trial[0];
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 1] = 0.5 * elastic_strain_tensor_trial[3];
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 2] = 0.5 * (elastic_strain_tensor_trial[1] + elastic_strain_tensor_trial[2]);
				}
				else if (info->DIMENSION == 3)
				{
					// tensor to voigt
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 0] = 0.5 * elastic_strain_tensor_trial[0];
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 1] = 0.5 * elastic_strain_tensor_trial[4];
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 2] = 0.5 * elastic_strain_tensor_trial[8];
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 3] = 0.5 * (elastic_strain_tensor_trial[1] + elastic_strain_tensor_trial[3]);
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 4] = 0.5 * (elastic_strain_tensor_trial[5] + elastic_strain_tensor_trial[7]);
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + 5] = 0.5 * (elastic_strain_tensor_trial[2] + elastic_strain_tensor_trial[6]);
				}
			}
			// 微小変形
			else if (info->c.ANALYSIS_MODE_1 == 0)
			{
				// Calculate trial elastic strains: {epsilon}^trial = {epsilon} + [B] * {du}
				for (int k = 0; k < D_MATRIX_SIZE; k++)
				{
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] = info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k];
					for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] += B_linear[k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];
				}
				if (overlay)
				{
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; m++)
							for (int l = 0; l < info->DIMENSION; l++)
								info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] += B_linear_overlay[k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];
				}

				// if (info->c.B_BAR == 1)
				// {
				// 	int p = info->Element_patch[i];
				// 	vector<double> R(MAX_NO_CP_ON_ELEMENT);
				// 	vector<vector<double>> R_L2_M_inv_G(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->DIMENSION, 0.0));
				// 	shape_and_dshape(R.data(), &info->gp[i].para()[j * info->DIMENSION], p, dec_order);
				// 	for (int k = 0; k < dec_order->No_Control_point_ON_ELEMENT[p]; k++)
				// 	{
				// 		int current_cp = dec_order->Controlpoint_of_Element[p * MAX_NO_CP_ON_ELEMENT + k];
				// 		for (int m = 0; m < info->No_Control_point_ON_ELEMENT[p]; m++)
				// 			for (int l = 0; l < info->DIMENSION; l++)
				// 				R_L2_M_inv_G[m][l] += R[k] * dec_order->L2_M_inv_G[p][current_cp][m * info->DIMENSION + l];
				// 	}

				// 	// Calculate mean strain increment: Delta{epsilon}_mean = 1/3 * sum({R_L2_M_inv_G}^T * {du})
				// 	double mean_strain_increment = 0.0;
				// 	for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
				// 		for (int l = 0; l < info->DIMENSION; l++)
				// 			mean_strain_increment += R_L2_M_inv_G[m][l] * info->disp_increment[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];

				// 	if (overlay)
				// 	{
				// 		int p_overlay = info->Element_patch[overlay_ele];
				// 		vector<double> R_overlay(MAX_NO_CP_ON_ELEMENT);
				// 		vector<vector<double>> R_L2_M_inv_G_overlay(info->No_Control_point_ON_ELEMENT[p_overlay], vector<double>(info->DIMENSION, 0.0));
				// 		shape_and_dshape(R_overlay.data(), &info->gp[i].opp_para_tilde()[j * info->DIMENSION], p_overlay, dec_order);
				// 		for (int k = 0; k < dec_order->No_Control_point_ON_ELEMENT[p_overlay]; k++)
				// 		{
				// 			int current_cp = dec_order->Controlpoint_of_Element[p_overlay * MAX_NO_CP_ON_ELEMENT + k];
				// 			for (int m = 0; m < info->No_Control_point_ON_ELEMENT[p_overlay]; m++)
				// 				for (int l = 0; l < info->DIMENSION; l++)
				// 					R_L2_M_inv_G_overlay[m][l] += R_overlay[k] * dec_order->L2_M_inv_G[p_overlay][current_cp][m * info->DIMENSION + l];
				// 		}

				// 		// Add overlay part
				// 		for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
				// 			for (int l = 0; l < info->DIMENSION; l++)
				// 				mean_strain_increment += R_L2_M_inv_G_overlay[m][l] * info->disp_increment[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];
				// 	}

				// 	mean_strain_increment /= 3.0;

				// 	// Subtract mean strain increment from trial elastic strains
				// 	for (int k = 0; k < info->DIMENSION; k++)
				// 		info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] += mean_strain_increment;
				// }
			}

			for (int k = 0; k < D_MATRIX_SIZE; k++)
				info->gp[i].elastic_strain_trial()[j * D_MATRIX_SIZE + k] = info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k];

			// Calculate trial stresses: {sigma}^trial = [D] * {epsilon}^trial
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				temp = 0.0;
				for (int l = 0; l < D_MATRIX_SIZE; l++)
					temp += D_e[k * D_MATRIX_SIZE + l] * info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + l];
				info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] = temp;
			}

			// Calculate trial relative stresses: {sigma}^trial = {sigma}^trial - {beta}
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				relative_stress_trial[k] = info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] - info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];

			// Calculate trial relative equivalent stress: sigma^trial_e
			double relative_equivalent_stress_trial = calc_equivalent_stress(relative_stress_trial);

			#if 0
			// Calculate previous_relative_stress
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				previous_relative_stress[k] = gp[i].previous_relative_stress_trial()[j * D_MATRIX_SIZE + k] - info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k];
				// previous_relative_stress[k] = info->gp[i].stress()[j * D_MATRIX_SIZE + k] - info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			#endif
			#if 0
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				gp[i].previous_relative_stress_trial()[j * D_MATRIX_SIZE + k] = relative_stress_trial[k] + info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			#endif
			#if 0
			// Calculate trial relative equivalent stress increment
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				increment_relative_stress_trial[k] = relative_stress_trial[k] - previous_relative_stress[k];

			// Calculate trial relative deviatoric stresses
			double hydrostatic_relative_stress_trial = (1.0 / 3.0) * (relative_stress_trial[0] + relative_stress_trial[1] + relative_stress_trial[2]);
			for (int k = 0; k < info->DIMENSION; k++)
				deviatoric_relative_stress_trial[k] += relative_stress_trial[k];
			for (int k = 0; k < info->DIMENSION; k++)
				deviatoric_relative_stress_trial[k] -= hydrostatic_relative_stress_trial;

			// Calculate unloading_check
			double unloading_check = 0.0;
			for (int k = 0; k < info->DIMENSION; k++)
				unloading_check += increment_relative_stress_trial[k] * deviatoric_relative_stress_trial[k];
			for (int k = info->DIMENSION; k < D_MATRIX_SIZE; k++)
				unloading_check += 2.0 * increment_relative_stress_trial[k] * deviatoric_relative_stress_trial[k];

			#endif
			if (relative_equivalent_stress_trial <= info->gp[i].yield_stress()[j])
			{
				elastic_counter++;

				info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j];
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			}
			#if 0
			else if (unloading_check < 0)
			{
				#pragma omp atomic
				unloading_counter++;

				info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j];
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			}
			#endif
			#if 0
			else if (unloading_check < 0)
			{
				#pragma omp atomic
				unloading_counter++;

				// Calculate equivalent plastic strain increment
				info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;

				// Calculate hardening stress increment
				double hardening_stress_increment = 0.0;

				// Calculate trial relative deviatoric stresses
				double current_relative_hydrostatic_stress = (1.0 / 3.0) * (relative_stress_trial[0] + relative_stress_trial[1] + relative_stress_trial[2]);
				for (int k = 0; k < info->DIMENSION; k++)
					relative_stress_trial[k] -= current_relative_hydrostatic_stress;

				// Calculate final elastic strains
				temp = info->gp[i].equivalent_plastic_strain_increment()[j] * 1.5 / relative_equivalent_stress_trial;
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= temp * relative_stress_trial[k];
				for (int k = info->DIMENSION; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= 2.0 * temp * relative_stress_trial[k];

				// Calculate final yield stress
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j] + (1.0 - info->kinematic_hardening_fraction) * hardening_stress_increment;

				// Calculate final back stresses
				temp = info->kinematic_hardening_fraction * hardening_stress_increment / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];

				// Calculate final stresses
				temp = info->gp[i].current_yield_stress()[j] / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] += current_relative_hydrostatic_stress;
			}
			#endif
			#if 0
			else if (relative_equivalent_stress_trial_increment < 0.0)
			{
				#pragma omp atomic
				unloading_counter++;

				// Calculate equivalent plastic strain increment
				info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;

				// Calculate hardening stress increment
				double hardening_stress_increment = 0.0;

				// Calculate trial relative deviatoric stresses
				double current_relative_hydrostatic_stress = (1.0 / 3.0) * (relative_stress_trial[0] + relative_stress_trial[1] + relative_stress_trial[2]);
				for (int k = 0; k < info->DIMENSION; k++)
					relative_stress_trial[k] -= current_relative_hydrostatic_stress;

				// Calculate final elastic strains
				temp = info->gp[i].equivalent_plastic_strain_increment()[j] * 1.5 / relative_equivalent_stress_trial;
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= temp * relative_stress_trial[k];
				for (int k = info->DIMENSION; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= 2.0 * temp * relative_stress_trial[k];

				// Calculate final yield stress
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j] + (1.0 - info->kinematic_hardening_fraction) * hardening_stress_increment;

				// Calculate final back stresses
				temp = info->kinematic_hardening_fraction * hardening_stress_increment / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];

				// Calculate final stresses
				temp = info->gp[i].current_yield_stress()[j] / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] += current_relative_hydrostatic_stress;
			}
			#endif
			else
			{
				plastic_counter++;

				// Calculate equivalent plastic strain increment
				info->gp[i].equivalent_plastic_strain_increment()[j] = calc_equivalent_plastic_strain_increment(relative_equivalent_stress_trial, info->gp[i].equivalent_plastic_strain()[j], info->gp[i].yield_stress()[j], info);

				// Calculate hardening stress increment
				double hardening_stress_increment = get_hardening_stress(info->gp[i].equivalent_plastic_strain()[j] + info->gp[i].equivalent_plastic_strain_increment()[j], info) - get_hardening_stress(info->gp[i].equivalent_plastic_strain()[j], info);

				// Calculate trial relative deviatoric stresses
				double current_relative_hydrostatic_stress = (1.0 / 3.0) * (relative_stress_trial[0] + relative_stress_trial[1] + relative_stress_trial[2]);
				for (int k = 0; k < info->DIMENSION; k++)
					relative_stress_trial[k] -= current_relative_hydrostatic_stress;

				// Calculate final elastic strains
				temp = info->gp[i].equivalent_plastic_strain_increment()[j] * 1.5 / relative_equivalent_stress_trial;
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= temp * relative_stress_trial[k];
				for (int k = info->DIMENSION; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] -= 2.0 * temp * relative_stress_trial[k];

				// Calculate final yield stress
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j] + (1.0 - info->kinematic_hardening_fraction) * hardening_stress_increment;

				// Calculate final back stresses
				temp = info->kinematic_hardening_fraction * hardening_stress_increment / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];

				// Calculate final stresses
				temp = info->gp[i].current_yield_stress()[j] / relative_equivalent_stress_trial;
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] + temp * relative_stress_trial[k];
				for (int k = 0; k < info->DIMENSION; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] += current_relative_hydrostatic_stress;
			}

			// Calculate Cauchy stresses from Kirchhoff stresses
			if (info->c.ANALYSIS_MODE_1 == 1)
			{
				double inverse_volume_change = 1.0 / calc_3x3_determinant(&info->gp[i].current_deformation_gradient()[j * info->DIMENSION * info->DIMENSION]);

				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] *= inverse_volume_change;
			}

			// unloading check
			#if 0
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				current_relative_stress[k] = info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] - info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k];
				previous_relative_stress[k] = info->gp[i].stress()[j * D_MATRIX_SIZE + k] - info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			}
			double relative_equivalent_stress_trial_increment = calc_equivalent_stress(current_relative_stress) - calc_equivalent_stress(previous_relative_stress);
			if (relative_equivalent_stress_trial_increment < 0.0)
			{
				unloading_counter++;
				info->gp[i].equivalent_plastic_strain_increment()[j] = 0.0;
				info->gp[i].current_yield_stress()[j] = info->gp[i].yield_stress()[j];
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					info->gp[i].current_back_stress()[j * D_MATRIX_SIZE + k] = info->gp[i].back_stress()[j * D_MATRIX_SIZE + k];
			}
			#endif

			// Internal force
			if (Total_mesh == 1)
			{
				// Calculate integral of [B]^T * {sigma} * dV
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < D_MATRIX_SIZE; m++)
						{
							int index = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
							int id = info->Index_Dof[index];
							if (id >= 0)
							{
								double val = B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * J * info->gp[i].w()[j];

								#pragma omp atomic
								info->internal_force[id] += val;
							}
						}
			}
			else if (Total_mesh > 1)
			{
				// skip gauss point which bf is zero
				if (info->c.BLENDING == 1 && info->gp[i].isSkip()[j])
					continue;

				double coef = (info->c.BLENDING == 1) ? J * info->gp[i].w()[j] * info->gp[i].blending_function()[j] : J * info->gp[i].w()[j];

				// Calculate integral of [B]^T * {sigma} * dV
				if (isAuxPatch == false)
				{
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
							for (int m = 0; m < D_MATRIX_SIZE; m++)
							{
								int index = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
								int id = info->Index_Dof[index];
								if (id >= 0)
								{
									double val = B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * coef;

									#pragma omp atomic
									info->internal_force[id] += val;
								}
							}
				}

				// Calculate f^G in Omega^L by using local elements
				if (i >= info->Total_Element_to_mesh[1])
				{
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
							for (int m = 0; m < D_MATRIX_SIZE; m++)
							{
								int index = info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
								int id = info->Index_Dof[index];
								if (id >= 0)
								{
									double val = B_linear_overlay[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * coef;

									#pragma omp atomic
									info->internal_force[id] += val;
								}
							}
				}
			}

			#if 0
			// Calculate integral of [B]^T * {sigma} * dV
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
					for (int m = 0; m < D_MATRIX_SIZE; m++)
					{
						int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						if (id >= 0)
						{
							#pragma omp critical
							{
								info->internal_force[id] += B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * J * info->gp[i].w()[j];
							}
						}
					}
			#endif
		}
	}

	/* debug stress zz */
	// for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	// {
	// 	for (int j = 0; j < GP_ON_ELEMENT; j++)
	// 	{
	// 		if (j == 0)
	// 		{
	// 			cout << "\033[31m";
	// 			cout << i << "\t" << info->gp[i].current_stress()[j * D_MATRIX_SIZE + 2] << endl;
	// 			cout << "\033[0m";
	// 		}
	// 	}
	// }

	// init
	for (auto it = ep_flag.begin(); it < ep_flag.end(); ++it)
		*it = false;

	if (elastic_counter > 0)
	{
		cout << "\033[31m";
		cout << "elastic" << endl;
		cout << "\033[0m";
		ep_flag[0] = true;
	}
	if (unloading_counter > 0)
	{
		cout << "\033[31m";
		cout << "unloading" << endl;
		cout << "\033[0m";
	}
	if (plastic_counter > 0)
	{
		cout << "\033[31m";
		cout << "plastic" << endl;
		cout << "\033[0m";
		ep_flag[1] = true;
	}

	#if 0
		FILE *fp_stress = fopen("test_Stress_ZZ.dat", "a");
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < info->gp[i].n(); j++)
			{
				fprintf(fp_stress, "%d %d %.20e\n", i, j, info->gp[i].current_stress()[j * D_MATRIX_SIZE + 2]);
			}
		}
		fprintf(fp_stress, "\n");
		fclose(fp_stress);
	#endif

	#if 0
		FILE *fp_IF = fopen("test_InternalForce.dat", "a");
		for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int id = info->Index_Dof[i * info->DIMENSION + j];
				if (id >= 0)
				{
					// glo
					if (i < info->Total_Control_Point_to_mesh[1])
						fprintf(fp_IF, "glo %d %.20e\n", i * info->DIMENSION + j, info->internal_force[id]);
					// loc
					else
						fprintf(fp_IF, "loc %d %.20e\n", i * info->DIMENSION + j, info->internal_force[id]);
				}
			}
		}
		fprintf(fp_IF, "\n");
		fclose(fp_IF);
	#endif
}

#if 0
void Update_internal_force_interpolate(information *info, information *aux)
{
	class alignas(64) interpolate_points {
	public:
		int le;           // aux's local element
		int legp;         // aux's local element gauss point

		int n;
		vector<double> w; // weight
		vector<int> ge;   // info's global element
		vector<int> gegp; // info's global element gauss point
	};

	static vector<vector<interpolate_points>> ip(aux->Total_Element_to_mesh[Total_mesh] - aux->Total_Element_to_mesh[1]);

	gp_switch(true, info);
	gp_switch(true, aux);

	// first time
	static bool isInitialized = false;
	if (!isInitialized)
	{
		// number of interpolation points
		// const int n = 1;
		const int n = pow_int((info->c.NG / 2) + (info->c.NG % 2), info->DIMENSION);
		// const int n = pow_int(info->c.NG, info->DIMENSION);

		// init
		#pragma omp parallel for
		for (size_t i = 0; i < ip.size(); i++)
		{
			ip[i].resize(aux->gp[i].n());
			for (size_t j = 0; j < ip[i].size(); j++)
			{
				ip[i][j].le = i + aux->Total_Element_to_mesh[1];
				ip[i][j].legp = j;
			}
		}

		// interpolate from info's global gauss points to aux's local gauss points
		#pragma omp parallel for schedule(dynamic)
		for (size_t i = 0; i < ip.size(); i++)
		{
			for (size_t j = 0; j < ip[i].size(); j++)
			{
				int aux_e = ip[i][j].le;
				int info_e = aux->gp[aux_e].opp_ele()[j];
				vector<distance_norm> dn;
				for (size_t k = 0; k < info->ecn[info_e].size(); k++)
				{
					int ge = info->ecn[info_e][k];
					for (int l = 0; l < info->gp[ge].n(); l++)
					{
						distance_norm temp;
						temp.e = ge;
						temp.gp_num = l;
						temp.r = 0.0;
						for (int m = 0; m < info->DIMENSION; m++)
							temp.r += pow(info->gp[ge].coord()[l * info->DIMENSION + m] - aux->gp[aux_e].coord()[j * info->DIMENSION + m], 2.0);
						temp.r = sqrt(temp.r);
						dn.emplace_back(temp);
					}
				}
				// search info's local gauss points
				vector<int> done;
				for (size_t k = 0; k < info->eoi[info_e].size(); k++)
				{
					int current_info_le = info->eoi[info_e][k];
					for (size_t l = 0; l < info->ecn[current_info_le].size(); l++)
					{
						int info_le = info->ecn[current_info_le][l];
						if (find(done.begin(), done.end(), info_le) != done.end())
							continue;
						done.emplace_back(info_le);
						for (int m = 0; m < info->gp[info_le].n(); m++)
						{
							distance_norm temp;
							temp.e = info_le;
							temp.gp_num = m;
							temp.r = 0.0;
							for (int nn = 0; nn < info->DIMENSION; nn++)
								temp.r += pow(info->gp[info_le].coord()[m * info->DIMENSION + nn] - aux->gp[aux_e].coord()[j * info->DIMENSION + nn], 2.0);
							temp.r = sqrt(temp.r);
							dn.emplace_back(temp);
						}
					}
				}

				// sort
				sort(dn.begin(), dn.end());

				// set
				ip[i][j].n = n;
				ip[i][j].w.resize(n);
				ip[i][j].ge.resize(n);
				ip[i][j].gegp.resize(n);
				for (int k = 0; k < n; k++)
				{
					ip[i][j].ge[k] = dn[k].e;
					ip[i][j].gegp[k] = dn[k].gp_num;
				}

				// weight
				double sum = 0.0;
				for (int k = 0; k < n; k++)
				{
					// ip[i][j].w[k] = 1.0 / dn[k].r;
					ip[i][j].w[k] = 1.0 / pow(dn[k].r, 2.0);
					sum += ip[i][j].w[k];
				}

				// normalize
				for (int k = 0; k < n; k++)
					ip[i][j].w[k] /= sum;
			}
		}

		isInitialized = true;
	}

	// interpolate from info's global gauss points to aux's local gauss points
	#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < ip.size(); i++)
	{
		for (size_t j = 0; j < ip[i].size(); j++)
		{
			int le = ip[i][j].le;
			int legp = ip[i][j].legp;

			// init
			aux->gp[le].current_yield_stress()[legp] = 0.0;
			aux->gp[le].equivalent_plastic_strain()[legp] = 0.0;
			aux->gp[le].equivalent_plastic_strain_increment()[legp] = 0.0;
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				aux->gp[le].elastic_strain_trial()[legp * D_MATRIX_SIZE + k] = 0.0;
				aux->gp[le].current_elastic_strain()[legp * D_MATRIX_SIZE + k] = 0.0;
				aux->gp[le].current_back_stress()[legp * D_MATRIX_SIZE + k] = 0.0;
				aux->gp[le].back_stress()[legp * D_MATRIX_SIZE + k] = 0.0;
				aux->gp[le].current_stress()[legp * D_MATRIX_SIZE + k] = 0.0;
			}

			// stored values
			for (int k = 0; k < ip[i][j].n; k++)
			{
				int ge = ip[i][j].ge[k];
				int gegp = ip[i][j].gegp[k];
				double w = ip[i][j].w[k];

				aux->gp[le].current_yield_stress()[legp] += info->gp[ge].current_yield_stress()[gegp] * w;
				aux->gp[le].equivalent_plastic_strain()[legp] += info->gp[ge].equivalent_plastic_strain()[gegp] * w;
				aux->gp[le].equivalent_plastic_strain_increment()[legp] += info->gp[ge].equivalent_plastic_strain_increment()[gegp] * w;
				for (int m = 0; m < D_MATRIX_SIZE; m++)
				{
					aux->gp[le].elastic_strain_trial()[legp * D_MATRIX_SIZE + m] += info->gp[ge].elastic_strain_trial()[gegp * D_MATRIX_SIZE + m] * w;
					aux->gp[le].current_elastic_strain()[legp * D_MATRIX_SIZE + m] += info->gp[ge].current_elastic_strain()[gegp * D_MATRIX_SIZE + m] * w;
					aux->gp[le].current_back_stress()[legp * D_MATRIX_SIZE + m] += info->gp[ge].current_back_stress()[gegp * D_MATRIX_SIZE + m] * w;
					aux->gp[le].back_stress()[legp * D_MATRIX_SIZE + m] += info->gp[ge].back_stress()[gegp * D_MATRIX_SIZE + m] * w;
					aux->gp[le].current_stress()[legp * D_MATRIX_SIZE + m] += info->gp[ge].current_stress()[gegp * D_MATRIX_SIZE + m] * w;
				}
			}
		}
	}

	// internal force
	#pragma omp parallel for schedule(dynamic)
	for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (int j = 0; j < aux->gp[i].n(); j++)
		{
			vector<double> B_linear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
			vector<double> B_linear_overlay(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

			int overlay_ele = aux->gp[i].opp_ele()[j];
			Make_B_Linear(overlay_ele, &aux->gp[i].opp_para_tilde()[j * aux->DIMENSION], B_linear_overlay.data(), aux);
			double J = Make_B_Linear(i, &aux->gp[i].para()[j * aux->DIMENSION], B_linear.data(), aux);
			double coef = J * aux->gp[i].w()[j] * aux->gp[i].blending_function()[j];

			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
					for (int m = 0; m < D_MATRIX_SIZE; m++)
					{
						int id = info->Index_Dof[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
						if (id >= 0)
						{
							double val = B_linear_overlay[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * coef;

							#pragma omp atomic
							info->internal_force[id] += val;
						}
					}
		}
	}

	#if 1
	// vtk file output, only points
	FILE* fp_IF = fopen("interpolate_points.vtk", "w");
	if (!fp_IF) {
		perror("Failed to open file");
		return;
	}

	fprintf(fp_IF, "# vtk DataFile Version 2.0\n");
	fprintf(fp_IF, "interpolate_points\n");
	fprintf(fp_IF, "ASCII\n");
	fprintf(fp_IF, "DATASET POLYDATA\n");

	// 全ポイント数を計算
	int totalPoints = 0;
	for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++) {
		totalPoints += aux->gp[i].n();
	}
	fprintf(fp_IF, "POINTS %d float\n", totalPoints);

	// ポイント座標の出力
	for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++) {
		for (int j = 0; j < aux->gp[i].n(); j++) {
			fprintf(fp_IF, "%.20e %.20e %.20e\n", 
					aux->gp[i].coord()[j * aux->DIMENSION], 
					aux->gp[i].coord()[j * aux->DIMENSION + 1], 
					aux->gp[i].coord()[j * aux->DIMENSION + 2]);
		}
	}

	// CELLS
	fprintf(fp_IF, "VERTICES %d %d\n", totalPoints, totalPoints * 2);
	for (int i = 0; i < totalPoints; i++) {
		fprintf(fp_IF, "1 %d\n", i);
	}

	// POINT_DATA セクション
	fprintf(fp_IF, "POINT_DATA %d\n", totalPoints);

	// current_yield_stress データ
	fprintf(fp_IF, "SCALARS current_yield_stress float\n");
	fprintf(fp_IF, "LOOKUP_TABLE default\n");
	for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++) {
		for (int j = 0; j < aux->gp[i].n(); j++) {
			fprintf(fp_IF, "%.20e\n", aux->gp[i].current_yield_stress()[j]);
		}
	}

	// stress z データ
	fprintf(fp_IF, "SCALARS current_stress float\n");
	fprintf(fp_IF, "LOOKUP_TABLE default\n");
	for (int i = aux->Total_Element_to_mesh[1]; i < aux->Total_Element_to_mesh[Total_mesh]; i++) {
		for (int j = 0; j < aux->gp[i].n(); j++) {
			fprintf(fp_IF, "%.20e\n", aux->gp[i].current_stress()[j * D_MATRIX_SIZE + 2]);
		}
	}

	fclose(fp_IF);
	#endif
}
#endif

void Update_internal_force_linear(vector<bool> &ep_flag, information *info)
{
	if (info->c.ANALYSIS_MODE_1 == 1)
	{
		cerr << "Error: Update_internal_force_linear cannot be used in large deformation analysis." << endl;
		exit(1);
	}

	// auxilary patch
	bool isAuxPatch = false;
	if (Total_mesh > 1 &&  info->c.BLENDING == 1 && info->aux == nullptr)
		isAuxPatch = true;
	cout << "isAuxPatch: " << isAuxPatch << endl;

	// init
	if (isAuxPatch == false)
	{
		for (int i = 0; i < K_Whole_Size; i++)
			info->internal_force[i] = 0.0;
	}

	// Generate elastic D matrix
	vector<double> D_e(D_MATRIX_SIZE * D_MATRIX_SIZE);
	D_elastic_smart(D_e.data(), info);

	long elastic_counter = 0;
	long plastic_counter = 0;
	long unloading_counter = 0;

	gp_switch(true, info);

	#pragma omp parallel for schedule(dynamic) reduction(+:elastic_counter, plastic_counter, unloading_counter) collapse(1)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		if (isAuxPatch && i < info->Total_Element_to_mesh[1])
			continue;

		double temp = 0.0;

		thread_local double *B_linear = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
		thread_local double *B_linear_overlay = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);

		for (int j = 0; j < info->gp[i].n(); j++)
		{
			bool overlay = false;
			int overlay_ele = -1;

			// search element
			if (Total_mesh > 1)
			{
				overlay = info->gp[i].isOverlay()[j];
				if (overlay)
				{
					overlay_ele = info->gp[i].opp_ele()[j];
					Make_B_Linear(overlay_ele, &info->gp[i].opp_para_tilde()[j * info->DIMENSION], overlay, i, &info->gp[i].para()[j * info->DIMENSION], B_linear_overlay, info);
				}
			}

			// Generate linear B matrix
			double J = Make_B_Linear(i, &info->gp[i].para()[j * info->DIMENSION], info->gp[i].isOverlay()[j], info->gp[i].opp_ele()[j], &info->gp[i].opp_para_tilde()[j * info->DIMENSION], B_linear, info);

			// Calculate trial elastic strains: {epsilon}^trial = {epsilon} + [B] * {du}
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] = info->gp[i].elastic_strain()[j * D_MATRIX_SIZE + k];
				for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
					for (int l = 0; l < info->DIMENSION; l++)
						info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] += B_linear[k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];
			}
			if (overlay)
			{
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; m++)
						for (int l = 0; l < info->DIMENSION; l++)
							info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + k] += B_linear_overlay[k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->disp_increment[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l];
			}

			// Calculate trial stresses: {sigma}^trial = [D] * {epsilon}^trial
			for (int k = 0; k < D_MATRIX_SIZE; k++)
			{
				temp = 0.0;
				for (int l = 0; l < D_MATRIX_SIZE; l++)
					temp += D_e[k * D_MATRIX_SIZE + l] * info->gp[i].current_elastic_strain()[j * D_MATRIX_SIZE + l];
				info->gp[i].current_stress()[j * D_MATRIX_SIZE + k] = temp;
			}

			elastic_counter++;

			// Internal force
			if (Total_mesh == 1)
			{
				// Calculate integral of [B]^T * {sigma} * dV
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < D_MATRIX_SIZE; m++)
						{
							int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							if (id >= 0)
							{
								double val = B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * J * info->gp[i].w()[j];

								#pragma omp atomic
								info->internal_force[id] += val;
							}
						}
			}
			else if (Total_mesh > 1)
			{
				// skip gauss point which bf is zero
				if (info->c.BLENDING == 1 && info->gp[i].isSkip()[j])
					continue;

				double coef = (info->c.BLENDING == 1) ? J * info->gp[i].w()[j] * info->gp[i].blending_function()[j] : J * info->gp[i].w()[j];

				// Calculate integral of [B]^T * {sigma} * dV
				if (isAuxPatch == false)
				{
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
							for (int m = 0; m < D_MATRIX_SIZE; m++)
							{
								int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
								if (id >= 0)
								{
									double val = B_linear[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * coef;

									#pragma omp atomic
									info->internal_force[id] += val;
								}
							}
				}

				// Calculate f^G in Omega^L by using local elements
				if (i >= info->Total_Element_to_mesh[1])
				{
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
							for (int m = 0; m < D_MATRIX_SIZE; m++)
							{
								int id = info->Index_Dof[info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
								if (id >= 0)
								{
									double val = B_linear_overlay[m * MAX_KIEL_SIZE + k * info->DIMENSION + l] * info->gp[i].current_stress()[j * D_MATRIX_SIZE + m] * coef;

									#pragma omp atomic
									info->internal_force[id] += val;
								}
							}
				}
			}
		}
	}

	// init
	for (auto it = ep_flag.begin(); it < ep_flag.end(); ++it)
		*it = false;

	if (elastic_counter > 0)
	{
		cout << "\033[31m";
		cout << "elastic" << endl;
		cout << "\033[0m";
		ep_flag[0] = true;
	}
	if (unloading_counter > 0)
	{
		cout << "\033[31m";
		cout << "unloading" << endl;
		cout << "\033[0m";
	}
	if (plastic_counter > 0)
	{
		cout << "\033[31m";
		cout << "plastic" << endl;
		cout << "\033[0m";
		ep_flag[1] = true;
	}
}


// External Force
void Update_external_force(information *info)
{
	#pragma omp parallel for
	for (int i = 0; i < K_Whole_Size; i++)
		info->external_force[i] = info->F[i];
}


// Force vector
void Update_force_vector(information *info, double factor, const int itr)
{
	// auxilary patch
	bool isAuxPatch = false;
	if (Total_mesh > 1 && info->c.BLENDING == 1 && info->aux == nullptr)
		isAuxPatch = true;
	cout << "isAuxPatch: " << isAuxPatch << endl;

	// {f_N} = {f^ext} - {f^int}
	if (isAuxPatch == false)
	{
		#pragma omp parallel for
		for (int i = 0; i < K_Whole_Size; i++)
			info->rhs_vec[i] = info->external_force[i] - info->internal_force[i];
	}

	// {f_N} -= [K_ND] * {u_D}
	if (itr == 0)
	{
		gp_switch(true, info);
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			vector<double> K_linear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
			vector<double> K_nonlinear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
			vector<double> K_linear_conditional(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
			vector<double> K_nonlinear_conditional(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

			int disp_const_counter = 0;
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					if (id < 0)
					{
						#pragma omp atomic
						disp_const_counter++;
					}
				}

			if (disp_const_counter > 0)
			{
				if (Total_mesh == 1 || (Total_mesh > 1 && info->c.BLENDING == 0))
				{
					// K_linear
					Calc_K_linear_EL(i, K_linear.data(), info);

					// K_nonlinear
					if (info->c.ANALYSIS_MODE_1 == 1)
					{
						Calc_K_nonlinear_EL(i, K_nonlinear.data(), info);
						for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
							for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
							{
								int id = j * MAX_KIEL_SIZE + k;
								double val = K_nonlinear[id];

								K_linear[id] += val;
							}
					}
				}
				else if (Total_mesh > 1 && info->c.BLENDING == 1)
				{
					// [K_L]^GG + [K_NL]^GG by global, local and axiliary elements
					if (i < info->Total_Element_to_mesh[1])
					{
						// [K_L]^GG + [K_NL]^GG by global elements
						if (isAuxPatch == false)
						{
							int dummy_element = 0;

							// K_linear based on global elements
							Calc_K_linear_EL_conditional_for_glo(false, i, dummy_element, K_linear.data(), info);

							// K_nonlinear based on global elements
							if (info->c.ANALYSIS_MODE_1 == 1)
							{
								Calc_K_nonlinear_EL_conditional_for_glo(false, i, dummy_element, K_nonlinear.data(), info);
								for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
									for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
									{
										int id = j * MAX_KIEL_SIZE + k;
										double val = K_nonlinear[id];

										K_linear[id] += val;
									}
							}
						}
						else if (isAuxPatch == true)
						{
							// init
							for (int j = 0; j < MAX_KIEL_SIZE * MAX_KIEL_SIZE; j++)
								K_linear[j] = 0.0;
						}

						// [K_L]^GG + [K_NL]^GG by local or auxiliary elements
						for (size_t j = 0; j < info->eoi[i].size(); j++)
						{
							// K_linear based on local elements
							Calc_K_linear_EL_conditional_for_glo(true, i, info->eoi[i][j], K_linear_conditional.data(), info);

							// K_nonlinear based on local elements
							if (info->c.ANALYSIS_MODE_1 == 1)
							{
								Calc_K_nonlinear_EL_conditional_for_glo(true, i, info->eoi[i][j], K_nonlinear_conditional.data(), info);
								for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
									for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
									{
										int id = k * MAX_KIEL_SIZE + l;
										double val = K_nonlinear_conditional[id];

										K_linear_conditional[id] += val;
									}
							}

							for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
								for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
								{
									int id = k * MAX_KIEL_SIZE + l;
									double val = K_linear_conditional[id];

									K_linear[id] += val;
								}
						}
					}

					// [K_L]^LL + [K_NL]^LL by local elements
					else if (i >= info->Total_Element_to_mesh[1] && isAuxPatch == false)
					{
						// K_linear
						Calc_K_linear_EL(i, K_linear.data(), info);

						// K_nonlinear
						if (info->c.ANALYSIS_MODE_1 == 1)
						{
							Calc_K_nonlinear_EL(i, K_nonlinear.data(), info);
							for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; j++)
								for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; k++)
								{
									int id = j * MAX_KIEL_SIZE + k;
									double val = K_nonlinear[id];

									K_linear[id] += val;
								}
						}
					}

					else if (i >= info->Total_Element_to_mesh[1] && isAuxPatch == true)
						continue;
				}

				for (int j = 0; j < info->DIMENSION; j++)
				{
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
					{
						int ii = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + j;
						int b = info->Index_Dof[ii];
						if (b >= 0)
						{
							int ii_local = k * info->DIMENSION + j;
							for (int l = 0; l < info->DIMENSION; l++)
							{
								for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
								{
									int jj = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l;
									int bb = info->Index_Dof[jj];
									if (bb < 0)
									{
										int jj_local = m * info->DIMENSION + l;
										for (int kk_const = 0; kk_const < info->Total_Constraint_to_mesh[Total_mesh]; kk_const++)
										{
											if (info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] == info->Constraint_Node_Dir[kk_const * 2 + 0] && l == info->Constraint_Node_Dir[kk_const * 2 + 1])
											{
												int id = info->Constraint_Node_Dir[kk_const * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[kk_const * 2 + 1];

												#if 1
													double val = -K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] * factor - (info->disp[id] + info->disp_increment[id]));
													// double val = -K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));

													#pragma omp atomic
													info->rhs_vec[b] += val;
												#endif
												#if 0
													if (Total_mesh < 2)
													{
														#pragma omp atomic
														info->rhs_vec[b] -= K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));
													}
													else
													{
														#pragma omp atomic
														info->rhs_vec[b] -= K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp_overlay[id] + info->disp_overlay_increment[id]));
													}
												#endif
												// info->rhs_vec[b] -= K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));
												// info->rhs_vec[b] -= K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * (info->Value_of_Constraint[kk_const] * (double)(1.0) / (double)(info->c.STEP_N));
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

	if (isAuxPatch)
		return;

	if (itr == 0 && Total_mesh > 1)
	{
		// for coupling
		gp_switch(true, info);
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			vector<double> K_linear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
			vector<double> K_nonlinear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

			if (Total_mesh >= 2 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
			{
				for (size_t j = 0; j < info->eoi[i].size(); j++)
				{
					int loc_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; l++)
						{
							int loc_id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (loc_id < 0)
							{
								#pragma omp atomic
								loc_disp_const_counter++;
							}
						}
					}

					int glo_disp_const_counter = 0;
					for (int k = 0; k < info->DIMENSION; k++)
					{
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; l++)
						{
							int glo_id = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
							if (glo_id < 0)
							{
								#pragma omp atomic
								glo_disp_const_counter++;
							}
						}
					}

					if (loc_disp_const_counter > 0 || glo_disp_const_counter > 0)
					{
						// K_linear
						Calc_coupled_K_linear_EL(i, info->eoi[i][j], K_linear_coupled.data(), info);

						// K_nonlinear
						if (info->c.ANALYSIS_MODE_1 == 1)
						{
							Calc_coupled_K_nonlinear_EL(i, info->eoi[i][j], K_nonlinear_coupled.data(), info);
							for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]] * info->DIMENSION; k++)
								for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION; l++)
								{
									int id = k * MAX_KIEL_SIZE + l;
									double val = K_nonlinear_coupled[id];

									K_linear_coupled[id] += val;
								}
						}
					}

					if (loc_disp_const_counter > 0)
					{
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														int id = info->Constraint_Node_Dir[n * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[n * 2 + 1];
														#if 1
															double val = -K_linear_coupled[row_local * MAX_KIEL_SIZE + col_local] * (info->Value_of_Constraint[n] * factor - (info->disp[id] + info->disp_increment[id]));
															// double val = -K_linear_coupled[row_local * MAX_KIEL_SIZE + col_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));

															#pragma omp atomic
															info->rhs_vec[row] += val;
														#endif
														#if 0
															if (Total_mesh < 2)
															{
																#pragma omp atomic
																info->rhs_vec[row] -= K_linear_coupled[row_local * MAX_KIEL_SIZE + col_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));
															}
															else
															{
																#pragma omp atomic
																info->rhs_vec[row] -= K_linear_coupled[row_local * MAX_KIEL_SIZE + col_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp_overlay[id] + info->disp_overlay_increment[id]));
															}
														#endif
													}
												}
											}
										}
								}
							}
					}

					if (glo_disp_const_counter > 0)
					{
						for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
							for (int j2 = 0; j2 < info->DIMENSION; j2++)
							{
								int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
								if (row >= 0)
								{
									int row_local = j1 * info->DIMENSION + j2;
									for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k1++)
										for (int k2 = 0; k2 < info->DIMENSION; k2++)
										{
											int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
											if (col < 0)
											{
												int col_local = k1 * info->DIMENSION + k2;
												for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
												{
													if (info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
													{
														int id = info->Constraint_Node_Dir[n * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[n * 2 + 1];
														#if 1
															double val = -K_linear_coupled[col_local * MAX_KIEL_SIZE + row_local] * (info->Value_of_Constraint[n] * factor - (info->disp[id] + info->disp_increment[id]));
															// double val = -K_linear_coupled[col_local * MAX_KIEL_SIZE + row_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));

															#pragma omp atomic
															info->rhs_vec[row] += val;
														#endif
														#if 0
															if (Total_mesh < 2)
															{
																#pragma omp atomic
																info->rhs_vec[row] -= K_linear_coupled[col_local * MAX_KIEL_SIZE + row_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp[id] + info->disp_increment[id]));
															}
															else
															{
																#pragma omp atomic
																info->rhs_vec[row] -= K_linear_coupled[col_local * MAX_KIEL_SIZE + row_local] * (info->Value_of_Constraint[n] * (double)(step + 1.0) / (double)(info->c.STEP_N) - (info->disp_overlay[id] + info->disp_increment[id]));
															}
														#endif
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

	// if (itr == 0 && step == 0)
	if (itr == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < K_Whole_Size; i++)
		{
			info->rhs_vec_initial[i] = info->rhs_vec[i];
			// info->rhs_vec_initial[i] = info->external_force[i] - info->internal_force[i];
		}
	}
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


double calc_equivalent_plastic_strain_increment(const double relative_equivalent_stress_trial, const double equivalent_plastic_strain, const double yield_stress, information *info)
{
	constexpr int max_iteratoin = 1000;
	constexpr double tolerance = TOLERANCE;

	const double three_times_shear_moduls = 3.0 * 0.5 * E / (1.0 + nu);
	const double hardening_stress = get_hardening_stress(equivalent_plastic_strain, info);

	double residual, residual_gradient;
	double d_ep = 0.0;

	// Solve 3 * G * dep - (se - sy - (sh (ep + dep) - sh (ep))) = 0
	for (int i = 0; i < max_iteratoin; i++)
	{
		double current_hardening_stress = get_hardening_stress(equivalent_plastic_strain + d_ep, info);
		residual = three_times_shear_moduls * d_ep - (relative_equivalent_stress_trial - yield_stress - (current_hardening_stress - hardening_stress));

		if (fabs(residual) <= tolerance)
			return d_ep;

		double current_hardening_modulus = get_hardening_modulus(equivalent_plastic_strain + d_ep, info);
		residual_gradient = three_times_shear_moduls + current_hardening_modulus;
		d_ep -= residual / residual_gradient;
	}

	return nan("");
}


void get_material(information *info)
{
	int temp_i;
	double temp_d;

	FILE *fp = fopen(info->opt_files[1].c_str(), "r");

	// kinematic_hardening_fraction
	fscanf(fp, "%lf", &temp_d);
	info->kinematic_hardening_fraction = temp_d;

	fscanf(fp, "%d", &temp_i);
	info->ss_curve_num = temp_i;
	info->ss_curve_stress = (double *)malloc(sizeof(double) * info->ss_curve_num);
	info->ss_curve_plastic_strain = (double *)malloc(sizeof(double) * info->ss_curve_num);
	for (int i = 0; i < info->ss_curve_num; i++)
	{
		fscanf(fp, "%lf", &temp_d);
		info->ss_curve_plastic_strain[i] = temp_d;
		fscanf(fp, "%lf", &temp_d);
		info->ss_curve_stress[i] = temp_d;
	}

	if (info->c.PARAVIEW_CRACK_REPRESENTATION)
	{
		fp = fopen(info->opt_files[2].c_str(), "r");
		info->crack_pair_n_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));

		int total_crack_pair_n = 0;
		for(int i = 0; i < Total_mesh; i++)
		{
			fscanf(fp, "%d", &temp_i);
			info->crack_pair_n_to_mesh[i] = total_crack_pair_n;
			total_crack_pair_n += temp_i;
		}
		info->crack_pair_n_to_mesh[Total_mesh] = total_crack_pair_n;

		info->crack_pair = (int *)malloc(sizeof(int) * total_crack_pair_n * 2);
		for(int i = 0; i < total_crack_pair_n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->crack_pair[i * 2 + j] = temp_i;
			}
		}
	}
}


double get_hardening_stress(const double equivalent_plastic_strain, information *info)
{
	if (info->ss_curve_num == 1)
		return info->ss_curve_stress[0];
	
	if (equivalent_plastic_strain < info->ss_curve_plastic_strain[0])
	{
		double sy_i0 = info->ss_curve_stress[0];
		double sy_i1 = info->ss_curve_stress[1];

		double d_ep = info->ss_curve_plastic_strain[1] - info->ss_curve_plastic_strain[0];
		double d_ep_i0 = equivalent_plastic_strain - info->ss_curve_plastic_strain[0];
		double d_ep_i1 = info->ss_curve_plastic_strain[1] - equivalent_plastic_strain;

		return (d_ep_i1 * sy_i0 + d_ep_i0 * sy_i1) / d_ep;
	}

	if (equivalent_plastic_strain > info->ss_curve_plastic_strain[info->ss_curve_num - 1])
	{
		double sy_i0 = info->ss_curve_stress[info->ss_curve_num - 2];
		double sy_i1 = info->ss_curve_stress[info->ss_curve_num - 1];

		double d_ep = info->ss_curve_plastic_strain[info->ss_curve_num - 1] - info->ss_curve_plastic_strain[info->ss_curve_num - 2];
		double d_ep_i0 = equivalent_plastic_strain - info->ss_curve_plastic_strain[info->ss_curve_num - 2];
		double d_ep_i1 = info->ss_curve_plastic_strain[info->ss_curve_num - 1] - equivalent_plastic_strain;

		return (d_ep_i1 * sy_i0 + d_ep_i0 * sy_i1) / d_ep;
	}

	for (int i = 0; i < info->ss_curve_num; i++)
		if (equivalent_plastic_strain >= info->ss_curve_plastic_strain[i] && equivalent_plastic_strain <= info->ss_curve_plastic_strain[i + 1])
		{
			double sy_i0 = info->ss_curve_stress[i];
			double sy_i1 = info->ss_curve_stress[i + 1];

			double d_ep = info->ss_curve_plastic_strain[i + 1] - info->ss_curve_plastic_strain[i];
			double d_ep_i0 = equivalent_plastic_strain - info->ss_curve_plastic_strain[i];
			double d_ep_i1 = info->ss_curve_plastic_strain[i + 1] - equivalent_plastic_strain;

			return (d_ep_i1 * sy_i0 + d_ep_i0 * sy_i1) / d_ep;
		}
	
	return nan("");
}


double get_hardening_modulus(const double equivalent_plastic_strain, information *info)
{
	if (info->ss_curve_num == 0)
		return 0.0;

	if (info->ss_curve_num == 1)
		return 0.0;

	if (equivalent_plastic_strain < info->ss_curve_plastic_strain[0])
	{
		double dsy = info->ss_curve_stress[1] - info->ss_curve_stress[0];
		double dep = info->ss_curve_plastic_strain[1] - info->ss_curve_plastic_strain[0];

		return dsy / dep;
	}

	if (equivalent_plastic_strain > info->ss_curve_plastic_strain[info->ss_curve_num - 1])
	{
		double dsy = info->ss_curve_stress[info->ss_curve_num - 1] - info->ss_curve_stress[info->ss_curve_num - 2];
		double dep = info->ss_curve_plastic_strain[info->ss_curve_num - 1] - info->ss_curve_plastic_strain[info->ss_curve_num - 2];

		return dsy / dep;
	}

	for (int i = 0; i < info->ss_curve_num - 1; i++)
		if (equivalent_plastic_strain >= info->ss_curve_plastic_strain[i] && equivalent_plastic_strain <= info->ss_curve_plastic_strain[i + 1])
		{
			double dsy = info->ss_curve_stress[i + 1] - info->ss_curve_stress[i];
			double dep = info->ss_curve_plastic_strain[i + 1] - info->ss_curve_plastic_strain[i];

			return dsy / dep;
		}

	return nan("");
}


// D matrix
void D_elastic_smart(double *D_e, information *info)
{
	thread_local bool is_initialized = false;
	thread_local vector<double> D_e_const(D_MATRIX_SIZE * D_MATRIX_SIZE, 0.0);
	thread_local vector<double> D_dil_const(D_MATRIX_SIZE * D_MATRIX_SIZE, 0.0);

	if (!is_initialized)
	{
		D_elastic(D_e_const.data(), info);
		
		if (info->c.B_BAR)
			D_dil(D_dil_const.data(), info);

		is_initialized = true;
	}

	for (int i = 0; i < D_MATRIX_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			D_e[i * D_MATRIX_SIZE + j] = D_e_const[i * D_MATRIX_SIZE + j] - D_dil_const[i * D_MATRIX_SIZE + j];

	return;
}


double D_elastic_dil_coef()
{
	thread_local bool is_initialized = false;
	thread_local double D_dil_const;

	if (!is_initialized)
	{
		double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
		double mu = E / (2.0 * (1.0 + nu));
		D_dil_const = lambda + (2.0 / 3.0) * mu;
		is_initialized = true;
	}

	return D_dil_const;
}


void D_elastic(double *D_e, information *info)
{
	int array_2[6] = {1, 1, 2, 2, 1, 2};
	int array_3[12] = {1, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1};
	int *array_ptr = (info->DIMENSION == 2 ? array_2 : array_3);
	int x[4];

	double K = E / (3.0 * (1.0 - 2.0 * nu));
	double G = E / (2.0 * (1.0 + nu));

	for (int i = 0; i < D_MATRIX_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
		{
			x[0] = array_ptr[i * 2];
			x[1] = array_ptr[i * 2 + 1];
			x[2] = array_ptr[j * 2];
			x[3] = array_ptr[j * 2 + 1];

			// Hencky model
			double d_ij = (x[0] == x[1] ? 1.0 : 0.0), d_kl = (x[2] == x[3] ? 1.0 : 0.0);
			double d_ik = (x[0] == x[2] ? 1.0 : 0.0), d_jl = (x[1] == x[3] ? 1.0 : 0.0);
			double d_il = (x[0] == x[3] ? 1.0 : 0.0), d_jk = (x[1] == x[2] ? 1.0 : 0.0);
			D_e[i * D_MATRIX_SIZE + j] = G * (d_ik * d_jl + d_il * d_jk) + (K - (2.0 / 3.0) * G) * d_ij * d_kl;
		}
}


void D_dil(double *D_dil, information *info)
{
	double K = E / (3.0 * (1.0 - 2.0 * nu));

	for (int i = 0; i < D_MATRIX_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			D_dil[i * D_MATRIX_SIZE + j] = 0.0;

	for (int i = 0; i < info->DIMENSION; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			D_dil[i * D_MATRIX_SIZE + j] = K;
}


void D_elastoplastic(double *D_ep, const double equivalent_plastic_strain, const double equivalent_plastic_strain_increment, const int ele, const int gauss_point, information *info)
{
	thread_local double *trial_relative_stresses = (double *)malloc(sizeof(double) * D_MATRIX_SIZE);

	// static const double bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));
	thread_local const double shear_modulus = 0.5 * E / (1.0 + nu);
	const double hardening_modulus = get_hardening_modulus(equivalent_plastic_strain + equivalent_plastic_strain_increment, info);

	// Calculate trial volumetric strain
	double trial_volumetric_strain = 0.0;
	for (int i = 0; i < info->DIMENSION; i++)
		trial_volumetric_strain += info->gp[ele].elastic_strain_trial()[gauss_point * D_MATRIX_SIZE + i];

	// Calculate trial relative deviatoric stresses
	for (int i = 0; i < info->DIMENSION; i++)
		trial_relative_stresses[i] = 2.0 * shear_modulus * (info->gp[ele].elastic_strain_trial()[gauss_point * D_MATRIX_SIZE + i] - trial_volumetric_strain / 3.0) - info->gp[ele].back_stress()[gauss_point * D_MATRIX_SIZE + i];
	for (int i = info->DIMENSION; i < D_MATRIX_SIZE; i++)
		trial_relative_stresses[i] = 2.0 * shear_modulus * 0.5 * info->gp[ele].elastic_strain_trial()[gauss_point * D_MATRIX_SIZE + i] - info->gp[ele].back_stress()[gauss_point * D_MATRIX_SIZE + i];

	// Calculate trial relative equivalent stress
	double trial_relative_equivalent_stress = calc_equivalent_stress(trial_relative_stresses);

	// Calculate elastic-plastic D matrix
	double temp = -equivalent_plastic_strain_increment * 6.0 * shear_modulus * shear_modulus / trial_relative_equivalent_stress;
	for (int i = 0; i < info->DIMENSION; i++)
		D_ep[i * D_MATRIX_SIZE + i] += temp;
	for (int i = info->DIMENSION; i < D_MATRIX_SIZE; i++)
		D_ep[i * D_MATRIX_SIZE + i] += 0.5 * temp;
	for (int i = 0; i < info->DIMENSION; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			D_ep[i * D_MATRIX_SIZE + j] -= temp / 3.0;
	double coefficient = 9.0 * shear_modulus * shear_modulus * (equivalent_plastic_strain_increment / trial_relative_equivalent_stress - 1.0 / (3.0 * shear_modulus + hardening_modulus));
	for (int i = 0; i < D_MATRIX_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			D_ep[i * D_MATRIX_SIZE + j] += coefficient * (trial_relative_stresses[i] / trial_relative_equivalent_stress) * (trial_relative_stresses[j] / trial_relative_equivalent_stress);
}


void Modify_D_matrix(double *D, double *current_stress, double *elastic_strain_trial, double *current_deformation_gradient)
{
	thread_local int tensor_size = pow_int(3, 4);
	thread_local int matrix_size = pow_int(3, 2);

	thread_local double *consistent_d_tensor = (double *)malloc(sizeof(double) * tensor_size);
	thread_local double *d_tensor = (double *)malloc(sizeof(double) * tensor_size);
	thread_local double *l_tensor = (double *)malloc(sizeof(double) * tensor_size);
	thread_local double *b_tensor = (double *)malloc(sizeof(double) * tensor_size);
	thread_local double *current_stress_tensor = (double *)malloc(sizeof(double) * matrix_size);
	thread_local double *trial_elastic_left_cauchy_green_deformations = (double *)malloc(sizeof(double) * matrix_size);
	thread_local double *trial_elastic_strain_tensor = (double *)malloc(sizeof(double) * matrix_size);
	thread_local double *identity_tensor = (double *)malloc(sizeof(double) * matrix_size);

	identify_tensor(3, identity_tensor);

	// Convert D matrix to D tensor
	convertSymmetric4thOrderMatrixToTensor(d_tensor, D);

	// Calculate trial elastic left Cauchy-Green deformations: [B]^trial = exp(2 * {epsilon}^trial)
	trial_elastic_strain_tensor[0 * 3 + 0] = 2.0 * elastic_strain_trial[0];
	trial_elastic_strain_tensor[0 * 3 + 1] = 2.0 * 0.5 * elastic_strain_trial[3];
	trial_elastic_strain_tensor[0 * 3 + 2] = 2.0 * 0.5 * elastic_strain_trial[5];
	trial_elastic_strain_tensor[1 * 3 + 0] = 2.0 * 0.5 * elastic_strain_trial[3];
	trial_elastic_strain_tensor[1 * 3 + 1] = 2.0 * elastic_strain_trial[1];
	trial_elastic_strain_tensor[1 * 3 + 2] = 2.0 * 0.5 * elastic_strain_trial[4];
	trial_elastic_strain_tensor[2 * 3 + 0] = 2.0 * 0.5 * elastic_strain_trial[5];
	trial_elastic_strain_tensor[2 * 3 + 1] = 2.0 * 0.5 * elastic_strain_trial[4];
	trial_elastic_strain_tensor[2 * 3 + 2] = 2.0 * elastic_strain_trial[2];

	tensor_exponential(3, trial_elastic_strain_tensor, trial_elastic_left_cauchy_green_deformations);

	// Calculate [L] = d(ln([B]^trial))/d[B]^trial
	tensor_logarithm_derivative(3, trial_elastic_left_cauchy_green_deformations, l_tensor);

	// Calculate B_ijkl = delta_ik * B^trial_jl + delta_jk * B^trial_il
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
					b_tensor[i * 27 + j * 9 + k * 3 + l] = identity_tensor[i * 3 + k] * trial_elastic_left_cauchy_green_deformations[j * 3 + l] + identity_tensor[j * 3 + k] * trial_elastic_left_cauchy_green_deformations[i * 3 + l];

	// Convert current stresses to current stress tensor
	current_stress_tensor[0 * 3 + 0] = current_stress[0];
	current_stress_tensor[0 * 3 + 1] = current_stress[3];
	current_stress_tensor[0 * 3 + 2] = current_stress[5];
	current_stress_tensor[1 * 3 + 0] = current_stress[3];
	current_stress_tensor[1 * 3 + 1] = current_stress[1];
	current_stress_tensor[1 * 3 + 2] = current_stress[4];
	current_stress_tensor[2 * 3 + 0] = current_stress[5];
	current_stress_tensor[2 * 3 + 1] = current_stress[4];
	current_stress_tensor[2 * 3 + 2] = current_stress[2];

	// Calculate inverse of volume change
	double inverse_volume_change = 1.0 / calc_3x3_determinant(current_deformation_gradient);

	// Calculate 1 / (2 * J) * ([D] : [L] : [B])_ijkl - (sigma_il * delta_jk + sigma_jl * delta_ik)
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
				{
					double d_i_j_k_l = 0.0;
					for (int ii = 0; ii < 3; ii++)
						for (int jj = 0; jj < 3; jj++)
						{
							double d_ii_jj_kk_ll_times_b_kk_ll_k_l = 0.0;
							for (int kk = 0; kk < 3; kk++)
								for (int ll = 0; ll < 3; ll++)
									d_ii_jj_kk_ll_times_b_kk_ll_k_l += l_tensor[ii * 27 + jj * 9 + kk * 3 + ll] * b_tensor[kk * 27 + ll * 9 + k * 3 + l];
							d_i_j_k_l += d_tensor[i * 27 + j * 9 + ii * 3 + jj] * d_ii_jj_kk_ll_times_b_kk_ll_k_l;
						}
					d_i_j_k_l *= 0.5 * inverse_volume_change;

					d_i_j_k_l -= current_stress_tensor[i * 3 + l] * identity_tensor[j * 3 + k];
					d_i_j_k_l -= current_stress_tensor[j * 3 + l] * identity_tensor[i * 3 + k];

					consistent_d_tensor[i * 27 + j * 9 + k * 3 + l] = d_i_j_k_l;
				}

	// Convert D tensor to D matrix
	convertSymmetric4thOrderTensorToMatrix(D, consistent_d_tensor);
}


int voigtIndex(const int i, const int j)
{
	if (i == j)
	{
		return i;
	}
	else if ((i == 1 && j == 2) || (i == 2 && j == 1))
	{
		return 3;
	}
	else if ((i == 0 && j == 2) || (i == 2 && j == 0))
	{
		return 4;
	}
	else if ((i == 0 && j == 1) || (i == 1 && j == 0))
	{
		return 5;
	}
	else
	{
		return -1;
	}
	
	return ERROR;
}


void convertSymmetric4thOrderMatrixToTensor(double *tensor, double *matrix)
{
	tensor[0 * 27 + 0 * 9 + 0 * 3 + 0] = matrix[0 * 6 + 0];
	tensor[0 * 27 + 0 * 9 + 0 * 3 + 1] = 0.5 * (matrix[0 * 6 + 3] + matrix[3 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 0 * 3 + 2] = 0.5 * (matrix[0 * 6 + 5] + matrix[5 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 1 * 3 + 0] = 0.5 * (matrix[0 * 6 + 3] + matrix[3 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 1 * 3 + 1] = 0.5 * (matrix[0 * 6 + 1] + matrix[1 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 1 * 3 + 2] = 0.5 * (matrix[0 * 6 + 4] + matrix[4 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 2 * 3 + 0] = 0.5 * (matrix[0 * 6 + 5] + matrix[5 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 2 * 3 + 1] = 0.5 * (matrix[0 * 6 + 4] + matrix[4 * 6 + 0]);
	tensor[0 * 27 + 0 * 9 + 2 * 3 + 2] = 0.5 * (matrix[0 * 6 + 2] + matrix[2 * 6 + 0]);
	tensor[0 * 27 + 1 * 9 + 0 * 3 + 0] = 0.5 * (matrix[3 * 6 + 0] + matrix[0 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 0 * 3 + 1] = matrix[3 * 6 + 3];
	tensor[0 * 27 + 1 * 9 + 0 * 3 + 2] = 0.5 * (matrix[3 * 6 + 5] + matrix[5 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 1 * 3 + 0] = 0.5 * (matrix[3 * 6 + 3] + matrix[3 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 1 * 3 + 1] = 0.5 * (matrix[3 * 6 + 1] + matrix[1 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 1 * 3 + 2] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 2 * 3 + 0] = 0.5 * (matrix[3 * 6 + 5] + matrix[5 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 2 * 3 + 1] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[0 * 27 + 1 * 9 + 2 * 3 + 2] = 0.5 * (matrix[3 * 6 + 2] + matrix[2 * 6 + 3]);
	tensor[0 * 27 + 2 * 9 + 0 * 3 + 0] = 0.5 * (matrix[5 * 6 + 0] + matrix[0 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 0 * 3 + 1] = 0.5 * (matrix[5 * 6 + 3] + matrix[3 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 0 * 3 + 2] = matrix[5 * 6 + 5];
	tensor[0 * 27 + 2 * 9 + 1 * 3 + 0] = 0.5 * (matrix[5 * 6 + 3] + matrix[3 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 1 * 3 + 1] = 0.5 * (matrix[5 * 6 + 1] + matrix[1 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 1 * 3 + 2] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 2 * 3 + 0] = 0.5 * (matrix[5 * 6 + 5] + matrix[5 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 2 * 3 + 1] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[0 * 27 + 2 * 9 + 2 * 3 + 2] = 0.5 * (matrix[5 * 6 + 2] + matrix[2 * 6 + 5]);
	tensor[1 * 27 + 0 * 9 + 0 * 3 + 0] = 0.5 * (matrix[3 * 6 + 0] + matrix[0 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 0 * 3 + 1] = 0.5 * (matrix[3 * 6 + 3] + matrix[3 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 0 * 3 + 2] = 0.5 * (matrix[3 * 6 + 5] + matrix[5 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 1 * 3 + 0] = matrix[3 * 6 + 3];
	tensor[1 * 27 + 0 * 9 + 1 * 3 + 1] = 0.5 * (matrix[3 * 6 + 1] + matrix[1 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 1 * 3 + 2] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 2 * 3 + 0] = 0.5 * (matrix[3 * 6 + 5] + matrix[5 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 2 * 3 + 1] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[1 * 27 + 0 * 9 + 2 * 3 + 2] = 0.5 * (matrix[3 * 6 + 2] + matrix[2 * 6 + 3]);
	tensor[1 * 27 + 1 * 9 + 0 * 3 + 0] = 0.5 * (matrix[0 * 6 + 1] + matrix[1 * 6 + 0]);
	tensor[1 * 27 + 1 * 9 + 0 * 3 + 1] = 0.5 * (matrix[3 * 6 + 1] + matrix[1 * 6 + 3]);
	tensor[1 * 27 + 1 * 9 + 0 * 3 + 2] = 0.5 * (matrix[5 * 6 + 1] + matrix[1 * 6 + 5]);
	tensor[1 * 27 + 1 * 9 + 1 * 3 + 0] = 0.5 * (matrix[3 * 6 + 1] + matrix[1 * 6 + 3]);
	tensor[1 * 27 + 1 * 9 + 1 * 3 + 1] = matrix[1 * 6 + 1];
	tensor[1 * 27 + 1 * 9 + 1 * 3 + 2] = 0.5 * (matrix[1 * 6 + 4] + matrix[4 * 6 + 1]);
	tensor[1 * 27 + 1 * 9 + 2 * 3 + 0] = 0.5 * (matrix[5 * 6 + 1] + matrix[1 * 6 + 5]);
	tensor[1 * 27 + 1 * 9 + 2 * 3 + 1] = 0.5 * (matrix[1 * 6 + 4] + matrix[4 * 6 + 1]);
	tensor[1 * 27 + 1 * 9 + 2 * 3 + 2] = 0.5 * (matrix[1 * 6 + 2] + matrix[2 * 6 + 1]);
	tensor[1 * 27 + 2 * 9 + 0 * 3 + 0] = 0.5 * (matrix[0 * 6 + 4] + matrix[4 * 6 + 0]);
	tensor[1 * 27 + 2 * 9 + 0 * 3 + 1] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[1 * 27 + 2 * 9 + 0 * 3 + 2] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[1 * 27 + 2 * 9 + 1 * 3 + 0] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[1 * 27 + 2 * 9 + 1 * 3 + 1] = 0.5 * (matrix[4 * 6 + 1] + matrix[1 * 6 + 4]);
	tensor[1 * 27 + 2 * 9 + 1 * 3 + 2] = matrix[4 * 6 + 4];
	tensor[1 * 27 + 2 * 9 + 2 * 3 + 0] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[1 * 27 + 2 * 9 + 2 * 3 + 1] = 0.5 * (matrix[4 * 6 + 4] + matrix[4 * 6 + 4]);
	tensor[1 * 27 + 2 * 9 + 2 * 3 + 2] = 0.5 * (matrix[4 * 6 + 2] + matrix[2 * 6 + 4]);
	tensor[2 * 27 + 0 * 9 + 0 * 3 + 0] = 0.5 * (matrix[5 * 6 + 0] + matrix[0 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 0 * 3 + 1] = 0.5 * (matrix[5 * 6 + 3] + matrix[3 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 0 * 3 + 2] = 0.5 * (matrix[5 * 6 + 5] + matrix[5 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 1 * 3 + 0] = 0.5 * (matrix[5 * 6 + 3] + matrix[3 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 1 * 3 + 1] = 0.5 * (matrix[5 * 6 + 1] + matrix[1 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 1 * 3 + 2] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 2 * 3 + 0] = matrix[5 * 6 + 5];
	tensor[2 * 27 + 0 * 9 + 2 * 3 + 1] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[2 * 27 + 0 * 9 + 2 * 3 + 2] = 0.5 * (matrix[5 * 6 + 2] + matrix[2 * 6 + 5]);
	tensor[2 * 27 + 1 * 9 + 0 * 3 + 0] = 0.5 * (matrix[0 * 6 + 4] + matrix[4 * 6 + 0]);
	tensor[2 * 27 + 1 * 9 + 0 * 3 + 1] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[2 * 27 + 1 * 9 + 0 * 3 + 2] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[2 * 27 + 1 * 9 + 1 * 3 + 0] = 0.5 * (matrix[3 * 6 + 4] + matrix[4 * 6 + 3]);
	tensor[2 * 27 + 1 * 9 + 1 * 3 + 1] = 0.5 * (matrix[4 * 6 + 1] + matrix[1 * 6 + 4]);
	tensor[2 * 27 + 1 * 9 + 1 * 3 + 2] = 0.5 * (matrix[4 * 6 + 4] + matrix[4 * 6 + 4]);
	tensor[2 * 27 + 1 * 9 + 2 * 3 + 0] = 0.5 * (matrix[5 * 6 + 4] + matrix[4 * 6 + 5]);
	tensor[2 * 27 + 1 * 9 + 2 * 3 + 1] = matrix[4 * 6 + 4];
	tensor[2 * 27 + 1 * 9 + 2 * 3 + 2] = 0.5 * (matrix[4 * 6 + 2] + matrix[2 * 6 + 4]);
	tensor[2 * 27 + 2 * 9 + 0 * 3 + 0] = 0.5 * (matrix[0 * 6 + 2] + matrix[2 * 6 + 0]);
	tensor[2 * 27 + 2 * 9 + 0 * 3 + 1] = 0.5 * (matrix[3 * 6 + 2] + matrix[2 * 6 + 3]);
	tensor[2 * 27 + 2 * 9 + 0 * 3 + 2] = 0.5 * (matrix[5 * 6 + 2] + matrix[2 * 6 + 5]);
	tensor[2 * 27 + 2 * 9 + 1 * 3 + 0] = 0.5 * (matrix[3 * 6 + 2] + matrix[2 * 6 + 3]);
	tensor[2 * 27 + 2 * 9 + 1 * 3 + 1] = 0.5 * (matrix[1 * 6 + 2] + matrix[2 * 6 + 1]);
	tensor[2 * 27 + 2 * 9 + 1 * 3 + 2] = 0.5 * (matrix[4 * 6 + 2] + matrix[2 * 6 + 4]);
	tensor[2 * 27 + 2 * 9 + 2 * 3 + 0] = 0.5 * (matrix[5 * 6 + 2] + matrix[2 * 6 + 5]);
	tensor[2 * 27 + 2 * 9 + 2 * 3 + 1] = 0.5 * (matrix[4 * 6 + 2] + matrix[2 * 6 + 4]);
	tensor[2 * 27 + 2 * 9 + 2 * 3 + 2] = matrix[2 * 6 + 2];
}


void convertSymmetric4thOrderTensorToMatrix(double *matrix, double *tensor)
{
	matrix[0 * 6 + 0] = tensor[0 * 27 + 0 * 9 + 0 * 3 + 0];
	matrix[0 * 6 + 1] = 0.5 * (tensor[0 * 27 + 0 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 0 * 3 + 0]);
	matrix[0 * 6 + 2] = 0.5 * (tensor[0 * 27 + 0 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 0 * 3 + 0]);
	matrix[0 * 6 + 3] = 0.25 * (tensor[0 * 27 + 0 * 9 + 0 * 3 + 1] + tensor[0 * 27 + 0 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 0 * 3 + 0] + tensor[1 * 27 + 0 * 9 + 0 * 3 + 0]);
	matrix[0 * 6 + 4] = 0.25 * (tensor[0 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 0 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 0 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 0]);
	matrix[0 * 6 + 5] = 0.25 * (tensor[0 * 27 + 0 * 9 + 2 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 2 * 9 + 0 * 3 + 0]);
	matrix[1 * 6 + 0] = 0.5 * (tensor[1 * 27 + 1 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 1 * 3 + 1]);
	matrix[1 * 6 + 1] = tensor[1 * 27 + 1 * 9 + 1 * 3 + 1];
	matrix[1 * 6 + 2] = 0.5 * (tensor[1 * 27 + 1 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 1 * 3 + 1]);
	matrix[1 * 6 + 3] = 0.25 * (tensor[1 * 27 + 1 * 9 + 0 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 0 * 9 + 1 * 3 + 1]);
	matrix[1 * 6 + 4] = 0.25 * (tensor[1 * 27 + 1 * 9 + 1 * 3 + 2] + tensor[1 * 27 + 1 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 1 * 3 + 1] + tensor[2 * 27 + 1 * 9 + 1 * 3 + 1]);
	matrix[1 * 6 + 5] = 0.25 * (tensor[1 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 1 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 1 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 1]);
	matrix[2 * 6 + 0] = 0.5 * (tensor[2 * 27 + 2 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 2 * 3 + 2]);
	matrix[2 * 6 + 1] = 0.5 * (tensor[2 * 27 + 2 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 2 * 3 + 2]);
	matrix[2 * 6 + 2] = tensor[2 * 27 + 2 * 9 + 2 * 3 + 2];
	matrix[2 * 6 + 3] = 0.25 * (tensor[2 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 2 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 2 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 2]);
	matrix[2 * 6 + 4] = 0.25 * (tensor[2 * 27 + 2 * 9 + 1 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 1 * 9 + 2 * 3 + 2]);
	matrix[2 * 6 + 5] = 0.25 * (tensor[2 * 27 + 2 * 9 + 2 * 3 + 0] + tensor[2 * 27 + 2 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 2 * 3 + 2] + tensor[0 * 27 + 2 * 9 + 2 * 3 + 2]);
	matrix[3 * 6 + 0] = 0.25 * (tensor[0 * 27 + 1 * 9 + 0 * 3 + 0] + tensor[1 * 27 + 0 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 0 * 3 + 1] + tensor[0 * 27 + 0 * 9 + 1 * 3 + 0]);
	matrix[3 * 6 + 1] = 0.25 * (tensor[0 * 27 + 1 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 0 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 0 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 1 * 3 + 0]);
	matrix[3 * 6 + 2] = 0.25 * (tensor[0 * 27 + 1 * 9 + 2 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 2 * 9 + 1 * 3 + 0]);
	matrix[3 * 6 + 3] = 0.25 * (tensor[0 * 27 + 1 * 9 + 0 * 3 + 1] + tensor[0 * 27 + 1 * 9 + 1 * 3 + 0] + tensor[1 * 27 + 0 * 9 + 0 * 3 + 1] + tensor[1 * 27 + 0 * 9 + 1 * 3 + 0]);
	matrix[3 * 6 + 4] = 0.125 * (tensor[0 * 27 + 1 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 1 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 1 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 1 * 9 + 1 * 3 + 0]);
	matrix[3 * 6 + 5] = 0.125 * (tensor[0 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 0 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 0 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 0 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 0]);
	matrix[4 * 6 + 0] = 0.25 * (tensor[1 * 27 + 2 * 9 + 0 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 0 * 9 + 2 * 3 + 1]);
	matrix[4 * 6 + 1] = 0.25 * (tensor[1 * 27 + 2 * 9 + 1 * 3 + 1] + tensor[2 * 27 + 1 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 1 * 3 + 2] + tensor[1 * 27 + 1 * 9 + 2 * 3 + 1]);
	matrix[4 * 6 + 2] = 0.25 * (tensor[1 * 27 + 2 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 1 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 1 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 2 * 3 + 1]);
	matrix[4 * 6 + 3] = 0.125 * (tensor[1 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 1 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 1 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 1 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 1]);
	matrix[4 * 6 + 4] = 0.25 * (tensor[1 * 27 + 2 * 9 + 1 * 3 + 2] + tensor[1 * 27 + 2 * 9 + 2 * 3 + 1] + tensor[2 * 27 + 1 * 9 + 1 * 3 + 2] + tensor[2 * 27 + 1 * 9 + 2 * 3 + 1]);
	matrix[4 * 6 + 5] = 0.125 * (tensor[1 * 27 + 2 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 2 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 2 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 2 * 9 + 2 * 3 + 1]);
	matrix[5 * 6 + 0] = 0.25 * (tensor[2 * 27 + 0 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 2 * 9 + 0 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 2 * 3 + 0] + tensor[0 * 27 + 0 * 9 + 0 * 3 + 2]);
	matrix[5 * 6 + 1] = 0.25 * (tensor[2 * 27 + 0 * 9 + 1 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 1] + tensor[1 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 1 * 9 + 0 * 3 + 2]);
	matrix[5 * 6 + 2] = 0.25 * (tensor[2 * 27 + 0 * 9 + 2 * 3 + 2] + tensor[0 * 27 + 2 * 9 + 2 * 3 + 2] + tensor[2 * 27 + 2 * 9 + 2 * 3 + 0] + tensor[2 * 27 + 2 * 9 + 0 * 3 + 2]);
	matrix[5 * 6 + 3] = 0.125 * (tensor[2 * 27 + 0 * 9 + 0 * 3 + 1] + tensor[2 * 27 + 0 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 2 * 9 + 0 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[0 * 27 + 1 * 9 + 0 * 3 + 2] + tensor[1 * 27 + 0 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 0 * 9 + 0 * 3 + 2]);
	matrix[5 * 6 + 4] = 0.125 * (tensor[2 * 27 + 0 * 9 + 1 * 3 + 2] + tensor[2 * 27 + 0 * 9 + 2 * 3 + 1] + tensor[0 * 27 + 2 * 9 + 1 * 3 + 2] + tensor[0 * 27 + 2 * 9 + 2 * 3 + 1] + tensor[1 * 27 + 2 * 9 + 2 * 3 + 0] + tensor[1 * 27 + 2 * 9 + 0 * 3 + 2] + tensor[2 * 27 + 1 * 9 + 2 * 3 + 0] + tensor[2 * 27 + 1 * 9 + 0 * 3 + 2]);
	matrix[5 * 6 + 5] = 0.25 * (tensor[2 * 27 + 0 * 9 + 2 * 3 + 0] + tensor[2 * 27 + 0 * 9 + 0 * 3 + 2] + tensor[0 * 27 + 2 * 9 + 2 * 3 + 0] + tensor[0 * 27 + 2 * 9 + 0 * 3 + 2]);
}


// blending function
void get_blending_data(information *info)
{
	ifstream file(info->opt_files[4]);
	if (!file.is_open())
	{
		cout << "Error: blending data file is not found." << endl;
		exit(1);
	}

	string line;
	bool coef_read = false;
	while (getline(file, line))
	{
		// skip comment lines
		if (line.empty() || line[0] == '#')
			continue;

		istringstream stream(line);
		if (!coef_read)
		{
			double coef;
			if (stream >> coef)
			{
				coef_read = true;
				info->curvature_radius_coef = coef;
				continue;
			}
			else
			{
				cerr << "Error: Invalid curvature radius coefficient format" << endl;
				exit(1);
			}
		}

		// blending faces
		blending_faces bf;
		if (stream >> bf.p >> bf.dir >> bf.face >> bf.area_ratio)
		{
			bf.p += info->Total_Patch_to_mesh[1];
			info->bf.emplace_back(bf);
			continue;
		}

		// Reset stream for fictitious area parsing
        stream.clear();
        stream.seekg(0);

		// fictitious area
		fictitious_area fa;
		if (stream >> fa.patch)
		{
			for (int i = 0; i < info->DIMENSION; i++)
			{
				double value[2];
				if (stream >> value[0] >> value[1])
				{
					fa.area_start.emplace_back(value[0]);
					fa.area_end.emplace_back(value[1]);
				}
			}
			info->fa.emplace_back(fa);
			continue;
		}

		cerr << "Error: Invalid blending data file format" << endl;
		exit(1);
	}

	file.close();
}


void get_auxiliary(information *info, char *option_file, char *global_file)
{
	info->aux = std::make_unique<information>();
	information *aux = info->aux.get();

	Get_Option(option_file, aux->opt_files);
	aux->c.get_data(aux->opt_files[0].c_str(), aux);

	// if empty auxiliary file, return
	if (aux->opt_files[5].empty())
	{
		cout << "No auxiliary file is provided." << endl;
		exit(1);
	}
	string auxiliary_file = aux->opt_files[5];
	Allocation(0, aux);
	Get_Input_1(0, global_file, aux);
	Get_Input_1(1, auxiliary_file.c_str(), aux);
	Allocation(1, aux);
	Get_Input_2(0, global_file, aux);
	Get_Input_2(1, auxiliary_file.c_str(), aux);
	Allocation(2, aux);
	Make_INC(aux);
	Allocation(3, aux);
	searchOverlappingEle(aux);
	Allocation(4, aux);

	aux->D = info->D;
	aux->Index_Dof = info->Index_Dof;
	aux->K_Whole_Val = info->K_Whole_Val;
	aux->K_Whole_Ptr = info->K_Whole_Ptr;
	aux->K_Whole_Col = info->K_Whole_Col;
	aux->full_K = info->full_K;
	aux->rhs_vec = info->rhs_vec;

	get_material(aux);
	get_load_surface(aux);
	get_blending_data(aux);

	// set gauss points
	Make_Gauss_points(false, aux);

	#if 0
	if (aux->c.MODE_EX == 0)
	{
		int n[2] = {0};
		n[0] = pow_int(info->c.NG, info->DIMENSION), n[1] = pow_int(info->c.NG_EXTEND, info->DIMENSION);

		// reset gauss points in info
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
			if (aux->eoi[i].size() > 0 && info->eoi[i].size() == 0)
				info->gp[i].setVar(true, false, info->Element_patch[i], i, n[0], n[1], info);

		// set gauss points in aux
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < aux->Total_Element_to_mesh[1]; i++)
			if (info->eoi[i].size() > 0 && aux->eoi[i].size() == 0)
				aux->gp[i].setVar(true, false, aux->Element_patch[i], i, n[0], n[1], aux);
	}
	#endif

	Allocation_init_NL(aux);

	aux->internal_force = info->internal_force;
	aux->rhs_vec_initial = info->rhs_vec_initial;
}


void setBlendingArea(information *info)
{
	cout << "Set blending area" << endl;
	information *aux = info->aux.get();
	aux->indices.resize(aux->Total_Patch_to_mesh[Total_mesh]);
	aux->corner_case.resize(aux->Total_Patch_to_mesh[Total_mesh]);

	if (aux->DIMENSION == 2)
	{
		cerr << "Error: Blending function is not supported in 2D" << endl;
		exit(1);
	}

	// search patch
	for (size_t i = 0; i < aux->bf.size(); i++)
		aux->indices[aux->bf[i].p].emplace_back(i);

	for (int i = 0; i < aux->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		if (aux->indices[i].size() == 0)
			continue;

		for (int j = 0; j < pow_int(2, aux->DIMENSION); j++)
		{
			bf_corner_case current_cc;
			int count = 0;

			short int bit[MAX_DIMENSION] = {0};
			bit[0] = (j & 1) ? 1 : 0;
			bit[1] = (j & 2) ? 1 : 0;
			bit[2] = (j & 4) ? 1 : 0;

			for (int k = 0; k < aux->DIMENSION; k++)
			{
				for (size_t l = 0; l < aux->indices[i].size(); l++)
				{
					int index = aux->indices[i][l];
					blending_faces &bf = aux->bf[index];
					if (bf.dir == k && bf.face == bit[k])
					{
						current_cc.axis_combination.emplace_back(k);
						current_cc.face_combination.emplace_back(bf.face);
						int knot_index_offset = aux->Total_Knot_to_patch_dim[i * aux->DIMENSION + k];
						double para_patch_start = aux->Position_Knots[knot_index_offset];
						double para_patch_end = aux->Position_Knots[knot_index_offset + aux->No_Control_point[i * aux->DIMENSION + k] + aux->Order[i * aux->DIMENSION + k]];
						double bf_range = bf.area_ratio * (para_patch_end - para_patch_start);
						if (bf.face == 0)
						{
							current_cc.area_start.emplace_back(para_patch_start);
							current_cc.area_end.emplace_back(para_patch_start + bf_range * (aux->curvature_radius_coef));
						}
						else if (bf.face == 1)
						{
							current_cc.area_start.emplace_back(para_patch_end - bf_range * (aux->curvature_radius_coef));
							current_cc.area_end.emplace_back(para_patch_end);
						}
						current_cc.ellipse_inner_radius.emplace_back(bf_range * aux->curvature_radius_coef);
						current_cc.ellipse_outer_radius.emplace_back(bf_range * (aux->curvature_radius_coef));
						current_cc.ellipse_center.emplace_back(para_patch_start + bf_range * (aux->curvature_radius_coef));
						count++;
					}
				}
			}

			if (count < 2)
				continue;
			else if (count == 2)
				current_cc.dim_type = 2;
			else if (count == 3)
				current_cc.dim_type = 3;

			aux->corner_case[i].emplace_back(current_cc);
		}
	}

	// sort corner case
	for (int i = 0; i < aux->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		if (aux->corner_case[i].size() == 0)
			continue;

		std::sort(aux->corner_case[i].begin(), aux->corner_case[i].end(), [](const bf_corner_case &a, const bf_corner_case &b) -> bool
		{
			return a.dim_type > b.dim_type;
		});
	}

	// output corner case
	for (int i = 0; i < aux->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		if (aux->corner_case[i].size() == 0)
			continue;

		cout << "Patch: " << i << endl;
		for (size_t j = 0; j < aux->corner_case[i].size(); j++)
		{
			bf_corner_case &cc = aux->corner_case[i][j];
			cout << "Dim type: " << cc.dim_type << endl;
			for (size_t k = 0; k < cc.axis_combination.size(); k++)
			{
				cout << "Axis: " << cc.axis_combination[k] << " Face: " << cc.face_combination[k] << " Start: " << cc.area_start[k] << " End: " << cc.area_end[k] << " Inner radius: " << cc.ellipse_inner_radius[k] << " Outer radius: " << cc.ellipse_outer_radius[k] << " Center: " << cc.ellipse_center[k] << endl;
			}
		}
	}

	// set blending function
	gp_switch(true, info);
	gp_switch(true, aux);

	// global and local patch
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		info->gp[i].setBlendingFunction(info);

	// global and auxiliary patch
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
		aux->gp[i].setBlendingFunction(aux);
}


// blending function using axiliary patch
double blendingFunction(const int e, const int p, double *para, information *target)
{
	// IGA
	if (Total_mesh == 1)
		return 1.0;

	int target_p = p;

	// global and local patch
	bool isGlobalPatch = (target_p < target->Total_Patch_on_mesh[0]) ? true : false;
	bool isGlobalAndLocalPatch = false;
	if (target->aux != nullptr)
	{
		isGlobalAndLocalPatch = true;

		// local patch
		if (!isGlobalPatch)
			return 1.0;
	}

	double (*func_ptr)(double, double *, double, int) = nullptr;
	#if 0
	// cos function
	func_ptr = cosFunction;
	#endif
	#if 1
	// cubic function
	func_ptr = cubicFunction;
	#endif
	#if 0
	// linear function
	func_ptr= linearFunction;
	#endif

	// check fictitious area
	if(isGlobalAndLocalPatch && target->fa.size() > 0 && fictitiousAreaCheck(e, target_p, para, target))
		return target->c.FICT_COEF;

	double out_knot[MAX_DIMENSION] = {0.0};
	if (isGlobalPatch)
	{
		// calc global coord
		double coord[MAX_DIMENSION] = {0.0};
		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), para, e, target);
		for (int i = 0; i < target->No_Control_point_ON_ELEMENT[target_p]; i++)
			for (int j = 0; j < target->DIMENSION; j++)
				coord[j] += R[i] * target->Node_Coordinate[target->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (target->DIMENSION + 1) + j];

		// search loc patch and calc loc knot
		bool isFoundInside = false;
		for (int i = target->Total_Patch_to_mesh[1]; i < target->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int itr_n = calc_patch_parameter_coord(coord, i, out_knot, target);
			if (itr_n != ERROR)
			{
				isFoundInside = true;

				// substitute target_p
				target_p = i;
				break;
			}
		}
		if (!isFoundInside)
			return 1.0;
		
		if (isGlobalAndLocalPatch)
		{
			if (isFoundInside)
				return 0.0;
			else if (!isFoundInside)
				return 1.0;
		}
	}

	// corner case
	#if 1
	if (target->corner_case[target_p].size() > 0)
	{
		vector<double> knot(target->DIMENSION);
		if (isGlobalPatch)
			for (int i = 0; i < target->DIMENSION; i++)
				knot[i] = out_knot[i];
		if (!isGlobalPatch)
			for (int i = 0; i < target->DIMENSION; i++)
			{
				int knot_index_offset = target->Total_Knot_to_patch_dim[target_p * target->DIMENSION + i];
				double para_element_start = target->Position_Knots[knot_index_offset + target->Order[target_p * target->DIMENSION + i] + target->ENC[e * target->DIMENSION + i]];
				double para_element_end = target->Position_Knots[knot_index_offset + target->Order[target_p * target->DIMENSION + i] + target->ENC[e * target->DIMENSION + i] + 1];
				knot[i] = 0.5 * (para[i] + 1.0) * (para_element_end - para_element_start) + para_element_start;
			}

		// search corner case
		for (size_t i = 0; i < target->corner_case[target_p].size(); i++)
		{
			bf_corner_case &cc = target->corner_case[target_p][i];
			bool isInside = true;
			for (int j = 0; j < cc.dim_type; j++)
			{
				int axis = cc.axis_combination[j];
				if (!(fabs(knot[axis] - cc.area_start[j]) < MERGE_ERROR || (cc.area_start[j] < knot[axis] && knot[axis] < cc.area_end[j]) || fabs(knot[axis] - cc.area_end[j]) < MERGE_ERROR))
				{
					isInside = false;
					break;
				}
			}

			if (isInside)
			{
				vector<double> modified_knot(target->DIMENSION);
				for (int j = 0; j < cc.dim_type; j++)
				{
					int axis = cc.axis_combination[j];

					if (cc.face_combination[j] == 0)
						modified_knot[j] = knot[axis];
					else if (cc.face_combination[j] == 1)
						modified_knot[j] = cc.area_end[j] - knot[axis];
				}

				// check singular point
				bool isSingular = true;
				for (int j = 0; j < cc.dim_type; j++)
				{
					if (!(fabs(modified_knot[j]) < MERGE_ERROR))
					{
						isSingular = false;
						break;
					}
				}
				if (isSingular)
					return 1.0;

				// modify reverse knot
				for (int j = 0; j < cc.dim_type; j++)
				{
					modified_knot[j] = cc.ellipse_outer_radius[j] - modified_knot[j];

					if (fabs(modified_knot[j] - cc.ellipse_outer_radius[j]) < MERGE_ERROR)
						modified_knot[j] = cc.ellipse_outer_radius[j];
					else if (fabs(modified_knot[j]) < MERGE_ERROR)
						modified_knot[j] = 0.0;
				}


				// vector<double> ellipse_intersection_coord_1(target->DIMENSION);
				vector<double> ellipse_intersection_coord_2(target->DIMENSION);
				// findEllipseIntersection(ellipse_intersection_coord_1, cc.ellipse_center, cc.ellipse_inner_radius, modified_knot, cc.dim_type);
				findEllipseIntersection(ellipse_intersection_coord_2, cc.ellipse_center, cc.ellipse_outer_radius, modified_knot, cc.dim_type);

				// calc line length
				double OA = 0.0, /*OB1 = 0.0,*/ OB2 = 0.0;
				for (int j = 0; j < cc.dim_type; j++)
				{
					OA += pow(modified_knot[j] - cc.ellipse_center[j], 2);
					// OB1 += pow(ellipse_intersection_coord_1[j] - cc.ellipse_center[j], 2);
					OB2 += pow(ellipse_intersection_coord_2[j] - cc.ellipse_center[j], 2);
				}
				OA = sqrt(OA);
				// OB1 = sqrt(OB1);
				OB2 = sqrt(OB2);

				// bool outOfBlendingArea1 = (OA < OB1) ? true : false;
				bool outOfBlendingArea2 = (OA > OB2) ? true : false;
				// if (outOfBlendingArea1)
				// {
				// 	if (isGlobalPatch)
				// 		return 0.0;
				// 	else if (!isGlobalPatch)
				// 		return 1.0;
				// }
				if (outOfBlendingArea2)
				{
					if (isGlobalPatch)
						return 1.0;
					else if (!isGlobalPatch)
						return 0.0;
				}

				#if 0
				vector<double> axis_intersection_coord(target->DIMENSION);
				findAxisIntersection(axis_intersection_coord, cc.ellipse_center, modified_knot, cc.dim_type);
				#endif

				// calc line length
				double B1B2 = 0.0, AB2 = 0.0;
				for (int j = 0; j < cc.dim_type; j++)
				{
					// B1B2 += pow(ellipse_intersection_coord_1[j] - ellipse_intersection_coord_2[j], 2);
					B1B2 += pow(cc.ellipse_outer_radius[j] - ellipse_intersection_coord_2[j], 2);
					AB2 += pow(ellipse_intersection_coord_2[j] - modified_knot[j], 2);
				}
				B1B2 = sqrt(B1B2);
				AB2 = sqrt(AB2);
				
				double postition = AB2 / B1B2;
				int face = 0;
				double blending_range = 1.0;
				double blending_range_array[2] = {0.0, 1.0};

				double temp_result = func_ptr(blending_range, blending_range_array, postition, face);

				if (isGlobalPatch)
					return temp_result;
				else if (!isGlobalPatch)
					return 1.0 - temp_result;
			}
		}

		// not found
		if (isGlobalPatch)
			return 1.0;
		else if (!isGlobalPatch)
			return 0.0;
	}
	#endif

	bool isFound = false;
	double result = 0.0;
	for (size_t i = 0; i < target->indices[target_p].size(); i++)
	{
		int index = target->indices[target_p][i];
		blending_faces &bf = target->bf[index];

		int knot_index_offset = target->Total_Knot_to_patch_dim[target_p * target->DIMENSION + bf.dir];
		double para_patch_start = target->Position_Knots[knot_index_offset];
		double para_patch_end = target->Position_Knots[knot_index_offset + target->No_Control_point[target_p * target->DIMENSION + bf.dir] + target->Order[target_p * target->DIMENSION + bf.dir]];
		double blending_range = bf.area_ratio * (para_patch_end - para_patch_start);
		double blending_range_start = (bf.face == 0) ? para_patch_start : para_patch_end - blending_range;
		double blending_range_end = (bf.face == 0) ? para_patch_start + blending_range : para_patch_end;
		double blending_range_array[2] = {blending_range_start, blending_range_end};

		// set knot
		double knot = 0.0;
		if (isGlobalPatch)
		{
			knot = out_knot[bf.dir];
		}
		if (!isGlobalPatch)
		{
			double para_element_start = target->Position_Knots[knot_index_offset + target->Order[target_p * target->DIMENSION + bf.dir] + target->ENC[e * target->DIMENSION + bf.dir]];
			double para_element_end = target->Position_Knots[knot_index_offset + target->Order[target_p * target->DIMENSION + bf.dir] + target->ENC[e * target->DIMENSION + bf.dir] + 1];
			knot = 0.5 * (para[bf.dir] + 1.0) * (para_element_end - para_element_start) + para_element_start;
		}

		// check blending area
		if (fabs(blending_range_start - knot) < MERGE_ERROR || (blending_range_start < knot && knot < blending_range_end) || fabs(blending_range_end - knot) < MERGE_ERROR)
		{
			isFound = true;
			result = func_ptr(blending_range, blending_range_array, knot, bf.face);
			break;
		}
	}

	// found
	if (isFound)
	{
		if (isGlobalPatch)
			return 1.0 - result;
		else if (!isGlobalPatch)
			return result;
	}

	// not found
	if (!isFound)
	{
		if (isGlobalPatch)
			return 1.0;
		else if (!isGlobalPatch)
			return 0.0;
	}

	return -1.0;
}


#if 0
double normalizedShapeFunction(int patch, int dim, double blending_range, double *blending_range_array, double knot, int face, information *info)
{
	thread_local int max_support_1D = 2 * MAX_ORDER + 2;
	int offset = (face == 0) ? 0 : info->No_Control_point[patch * info->DIMENSION + dim] - info->Order[patch * info->DIMENSION + dim] - 1;
	int p = info->Order[patch * info->DIMENSION + dim];
	int knot_n = info->No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];
	thread_local double *Shape = (double *)malloc(sizeof(double) * MAX_ORDER * max_support_1D);

	double result = 0.0;

	double modified_knot;
	if (face == 0)
	{
		double knot_range = knot_ptr[offset + p + 1] - knot_ptr[offset + p];
		modified_knot = blending_range_array[0] + knot_range * (knot - blending_range_array[0]) / blending_range;
	}
	else
	{
		double knot_range = knot_ptr[offset + p + 1] - knot_ptr[offset + p];
		modified_knot = blending_range_array[1] - knot_range * (blending_range_array[1] - knot) / blending_range;
	}

	int support = p + 1;
	for (int i = 0; i < support + p; i++)
	{
		if (fabs(knot_ptr[offset + i] - knot_ptr[offset + i + 1]) < MERGE_ERROR)
			Shape[i] = 0.0;

		else if (knot_ptr[offset + i] <= modified_knot && modified_knot < knot_ptr[offset + i + 1])
			Shape[i] = 1.0;

		else if (fabs(knot_ptr[offset + i + 1] - knot_ptr[knot_n - 1]) < MERGE_ERROR && knot_ptr[offset + i] <= modified_knot && modified_knot <= knot_ptr[offset + i + 1])
			Shape[i] = 1.0;

		else
			Shape[i] = 0.0;
	}

	double left_term, right_term;
	for (int i = 1; i <= p; i++)
	{
		int new_support = support + p - i;
		for (int j = 0; j < new_support; j++)
		{
			left_term = 0.0;
			right_term = 0.0;

			// left
			if ((fabs(modified_knot - knot_ptr[offset + j]) * Shape[(i - 1) * max_support_1D + j]) < MERGE_ERROR && knot_ptr[offset + j + i] - knot_ptr[offset + j] == 0)
				left_term = 0.0;
			else
				left_term = (modified_knot - knot_ptr[offset + j]) / (knot_ptr[offset + j + i] - knot_ptr[offset + j]) * Shape[(i - 1) * max_support_1D + j];

			// right
			if ((fabs(knot_ptr[offset + j + i + 1] - modified_knot) * Shape[(i - 1) * max_support_1D + j + 1]) < MERGE_ERROR && knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1] == 0)
				right_term = 0.0;
			else
				right_term = (knot_ptr[offset + j + i + 1] - modified_knot) / (knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1]) * Shape[(i - 1) * max_support_1D + j + 1];

			Shape[i * max_support_1D + j] = left_term + right_term;
		}

		if (i == p)
		{
			int index = (face == 0) ? 0 : support - 1;
			result = Shape[i * max_support_1D + index];
		}
	}

	result = (1.0 - result);
	if (fabs(result) < MERGE_ERROR)
		result = 0.0;
	if (fabs(result - 1.0) < MERGE_ERROR)
		result = 1.0;

	return result;
}
#endif


double cosFunction(double blending_range, double *blending_range_array, double knot, int face)
{
	// y = 0.5 * cos(pi * x) + 0.5
	double modified_x;
	if (face == 0)
	{
		modified_x = blending_range_array[0] + (knot - blending_range_array[0]) / blending_range;
	}
	else
	{
		modified_x = blending_range_array[1] - (blending_range_array[1] - knot) / blending_range;
		modified_x = 1.0 - modified_x;
	}
	double result = 0.5 * cos(M_PI * modified_x) + 0.5;
	return result;
}


double cubicFunction(double blending_range, double *blending_range_array, double knot, int face)
{
	// y = 2x^3 - 3x^2 + 1
	double modified_x;
	if (face == 0)
	{
		modified_x = blending_range_array[0] + (knot - blending_range_array[0]) / blending_range;
	}
	else
	{
		modified_x = blending_range_array[1] - (blending_range_array[1] - knot) / blending_range;
		modified_x = 1.0 - modified_x;
	}
	double result = 2.0 * modified_x * modified_x * modified_x - 3.0 * modified_x * modified_x + 1.0;
	return result;
}


double linearFunction(double blending_range, double *blending_range_array, double knot, int face)
{
	// y = -x + 1
	double modified_x;
	if (face == 0)
	{
		modified_x = blending_range_array[0] + (knot - blending_range_array[0]) / blending_range;
	}
	else
	{
		modified_x = blending_range_array[1] - (blending_range_array[1] - knot) / blending_range;
		modified_x = 1.0 - modified_x;
	}
	double result = -modified_x + 1.0;
	return result;
}


void findEllipseIntersection(vector<double> &result, const vector<double> &o, const vector<double> &r, const vector<double> &a, const int dim)
{
	vector<double> m(dim);
	for (int i = 0; i < dim; ++i)
		m[i] = o[i] - a[i];

	double a_coeff = 0.0, b_coeff = 0.0, c_coeff = -1.0;
	for (int i = 0; i < dim; ++i)
	{
		double m_square = m[i] * m[i];
		double r_square = r[i] * r[i];
		a_coeff += m_square / r_square;
		b_coeff -= 2 * m_square / r_square;
		c_coeff += m_square / r_square;
	}

	double discriminant = b_coeff * b_coeff - 4 * a_coeff * c_coeff;
	if (discriminant < 0)
	{
		cout << "Discriminant: " << discriminant << endl;
		throw std::runtime_error("No intersection with the ellipsoid/ellipse.");
	}

	double t = (-b_coeff - std::sqrt(discriminant)) / (2 * a_coeff);
	for (int i = 0; i < dim; ++i)
		result[i] = a[i] + t * m[i];
}


void findAxisIntersection(std::vector<double> &result, const std::vector<double> &o, const std::vector<double> &a, const int dim)
{
	std::vector<double> m(dim);
	for (int i = 0; i < dim; ++i)
		m[i] = o[i] - a[i];

	std::vector<double> t(dim, std::numeric_limits<double>::max());
	for (int i = 0; i < dim; ++i)
	{
		if (std::abs(m[i]) > MERGE_ERROR)
			t[i] = a[i] / m[i];
	}

	double min_t = std::numeric_limits<double>::max();
	int min_idx = -1;

	for (int i = 0; i < dim; ++i)
	{
		if (fabs(t[i]) < MERGE_ERROR)
			t[i] = 0.0;

		if (t[i] >= 0 && t[i] < min_t)
		{
			bool valid = true;
			std::vector<double> test_point = a;
			for (int j = 0; j < dim; ++j)
			{
				test_point[j] += t[i] * m[j];
				// 各次元に対して交点が座標軸平面内にあるかを確認
				if (test_point[j] < 0)
				{
					valid = false;
					break;
				}
			}
			if (valid)
			{
				min_t = t[i];
				min_idx = i;
			}
		}
	}

	if (min_t == std::numeric_limits<double>::max() || min_idx == -1)
		throw std::runtime_error("No valid intersection with the axis planes.");

	result = a;
	for (int i = 0; i < dim; ++i)
		result[i] -= min_t * m[i];
	result[min_idx] = 0.0;
}


void debugBlendingFunctionGlo(information *info, int itr)
{
	string filename = "blending_glo";
	filename += to_string(itr) + ".vtk";

	const int nelem = info->Total_Element_to_mesh[1];
	const int subdivision[3] = {1, 1, 1};
	const int numNodes[3] = {subdivision[0] + 1, subdivision[1] + 1, subdivision[2] + 1};
	const int nodesPerElement = numNodes[0] * numNodes[1] * numNodes[2];

	// パラメトリック座標の分割点を計算
	std::vector<std::vector<double>> subdivCoords(3);
	for (int dir = 0; dir < 3; ++dir)
	{
		subdivCoords[dir].resize(numNodes[dir]);
		for (int i = 0; i < numNodes[dir]; ++i)
		{
			subdivCoords[dir][i] = -1.0 + (2.0 * i) / (double)subdivision[dir];
		}
	}

	std::ofstream file(filename, std::ios::binary);

	// ヘッダー部分
	file << "# vtk DataFile Version 2.0\n";
	file << "Subdivided hexahedral mesh data\n";
	file << "BINARY\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	// 点の座標を出力
	const int totalPoints = nelem * nodesPerElement;
	file << "POINTS " << totalPoints << " float\n";

	// 各点の座標を計算
	std::vector<float> coordinates(totalPoints * 3);
	#pragma omp parallel for schedule(dynamic) collapse(3)
	for (int elem = 0; elem < nelem; ++elem)
		for (int k = 0; k < numNodes[2]; ++k)
			for (int j = 0; j < numNodes[1]; ++j)
				for (int i = 0; i < numNodes[0]; ++i)
				{
					double para[3] = {subdivCoords[0][i],
									  subdivCoords[1][j],
									  subdivCoords[2][k]};
					double xyz[3] = {0.0};
					vector<double> R(MAX_NO_CP_ON_ELEMENT);
					shape_and_dshape(R.data(), para, elem, info);
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[elem]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							int cp = info->Controlpoint_of_Element[elem * MAX_NO_CP_ON_ELEMENT + l];
							xyz[m] += R[l] * (info->Node_Coordinate[cp * (info->DIMENSION + 1) + m] + info->disp_overlay[cp * info->DIMENSION + m]);
						}
					int index = (elem * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i) * 3;
					coordinates[index] = static_cast<float>(xyz[0]);
					coordinates[index + 1] = static_cast<float>(xyz[1]);
					coordinates[index + 2] = static_cast<float>(xyz[2]);
				}

	// バイナリ形式で座標を書き込み
	for (auto &coord : coordinates)
	{
		swapEndian(coord); // エンディアン変換
	}
	file.write(reinterpret_cast<char *>(coordinates.data()), coordinates.size() * sizeof(float));

	// セルの接続情報を出力
	const int numSubElements = (subdivision[0]) * (subdivision[1]) * (subdivision[2]) * nelem;
	file << "\nCELLS " << numSubElements << " " << (numSubElements * 9) << "\n";

	std::vector<int32_t> connectivity(9); // セル当たりの接続情報
	for (int elem = 0; elem < nelem; ++elem)
	{
		const int elemOffset = elem * nodesPerElement;
		for (int k = 0; k < subdivision[2]; ++k)
			for (int j = 0; j < subdivision[1]; ++j)
				for (int i = 0; i < subdivision[0]; ++i)
				{
					// 1つのサブ要素の8つの節点番号を計算
					connectivity[0] = 8; // 頂点数
					connectivity[1] = elemOffset + i + j * numNodes[0] + k * numNodes[0] * numNodes[1];
					connectivity[2] = connectivity[1] + 1;
					connectivity[3] = connectivity[1] + numNodes[0] + 1;
					connectivity[4] = connectivity[1] + numNodes[0];
					connectivity[5] = connectivity[1] + numNodes[0] * numNodes[1];
					connectivity[6] = connectivity[5] + 1;
					connectivity[7] = connectivity[5] + numNodes[0] + 1;
					connectivity[8] = connectivity[5] + numNodes[0];

					// エンディアン変換して書き込み
					for (auto &value : connectivity)
					{
						swapEndian(value);
					}
					file.write(reinterpret_cast<char *>(connectivity.data()), connectivity.size() * sizeof(int32_t));
				}
	}

	// セルタイプの出力
	file << "\nCELL_TYPES " << numSubElements << "\n";
	int32_t cellType = 12; // 六面体要素
	swapEndian(cellType);
	for (int i = 0; i < numSubElements; ++i)
	{
		file.write(reinterpret_cast<char *>(&cellType), sizeof(int32_t));
	}

	// スカラー値の出力
	file << "\nPOINT_DATA " << totalPoints << "\n";
	file << "SCALARS value float 1\n";
	file << "LOOKUP_TABLE default\n";

	// 各点でのスカラー値を計算
	information *aux = info->aux.get();
	std::vector<float> values(totalPoints);
	#pragma omp parallel for schedule(dynamic) collapse(3)
	for (int elem = 0; elem < nelem; ++elem)
		for (int k = 0; k < numNodes[2]; ++k)
			for (int j = 0; j < numNodes[1]; ++j)
				for (int i = 0; i < numNodes[0]; ++i)
				{
					double para[3] = {subdivCoords[0][i], subdivCoords[1][j], subdivCoords[2][k]};
					double value = std::min(blendingFunction(elem, info->Element_patch[elem], para, info), blendingFunction(elem, aux->Element_patch[elem], para, aux));
					int index = elem * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i;
					values[index] = static_cast<float>(value);
				}

	// バイナリ形式でスカラー値を書き込み
	for (auto &value : values)
	{
		swapEndian(value);
	}
	file.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(float));

	file.close();
}


void debugBlendingFunctionLoc(information *info, int itr)
{
	string filename = "blending_loc";
	filename += to_string(itr) + ".vtk";

	const int nelem = info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1];
	const int subdivision[3] = {1, 1, 1};
	const int numNodes[3] = {subdivision[0] + 1, subdivision[1] + 1, subdivision[2] + 1};
	const int nodesPerElement = numNodes[0] * numNodes[1] * numNodes[2];

	// パラメトリック座標の分割点を計算
	std::vector<std::vector<double>> subdivCoords(3);
	for (int dir = 0; dir < 3; ++dir)
	{
		subdivCoords[dir].resize(numNodes[dir]);
		for (int i = 0; i < numNodes[dir]; ++i)
		{
			subdivCoords[dir][i] = -1.0 + (2.0 * i) / (double)subdivision[dir];
		}
	}

	std::ofstream file(filename, std::ios::binary);

	// ヘッダー部分
	file << "# vtk DataFile Version 2.0\n";
	file << "Subdivided hexahedral mesh data\n";
	file << "BINARY\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	// 点の座標を出力
	const int totalPoints = nelem * nodesPerElement;
	file << "POINTS " << totalPoints << " float\n";

	// 各点の座標を計算
	std::vector<float> coordinates(totalPoints * 3);
	#pragma omp parallel for schedule(dynamic) collapse(3)
	for (int elem = info->Total_Element_to_mesh[1]; elem < info->Total_Element_to_mesh[Total_mesh]; ++elem)
		for (int k = 0; k < numNodes[2]; ++k)
			for (int j = 0; j < numNodes[1]; ++j)
				for (int i = 0; i < numNodes[0]; ++i)
				{
					double para[3] = {subdivCoords[0][i],
									  subdivCoords[1][j],
									  subdivCoords[2][k]};
					double xyz[3] = {0.0};
					vector<double> R(MAX_NO_CP_ON_ELEMENT);
					shape_and_dshape(R.data(), para, elem, info);
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[elem]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							int cp = info->Controlpoint_of_Element[elem * MAX_NO_CP_ON_ELEMENT + l];
							xyz[m] += R[l] * (info->Node_Coordinate[cp * (info->DIMENSION + 1) + m] + info->disp_overlay[cp * info->DIMENSION + m]);
						}
					int index = ((elem - info->Total_Element_to_mesh[1]) * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i) * 3;
					coordinates[index] = static_cast<float>(xyz[0]);
					coordinates[index + 1] = static_cast<float>(xyz[1]);
					coordinates[index + 2] = static_cast<float>(xyz[2]);
				}

	// バイナリ形式で座標を書き込み
	for (auto &coord : coordinates)
	{
		swapEndian(coord); // エンディアン変換
	}
	file.write(reinterpret_cast<char *>(coordinates.data()), coordinates.size() * sizeof(float));

	// セルの接続情報を出力
	const int numSubElements = (subdivision[0]) * (subdivision[1]) * (subdivision[2]) * nelem;
	file << "\nCELLS " << numSubElements << " " << (numSubElements * 9) << "\n";

	std::vector<int32_t> connectivity(9); // セル当たりの接続情報
	for (int elem = 0; elem < info->Total_Element_to_mesh[Total_mesh] - info->Total_Element_to_mesh[1]; ++elem)
	{
		const int elemOffset = elem * nodesPerElement;
		for (int k = 0; k < subdivision[2]; ++k)
			for (int j = 0; j < subdivision[1]; ++j)
				for (int i = 0; i < subdivision[0]; ++i)
				{
					// 1つのサブ要素の8つの節点番号を計算
					connectivity[0] = 8; // 頂点数
					connectivity[1] = elemOffset + i + j * numNodes[0] + k * numNodes[0] * numNodes[1];
					connectivity[2] = connectivity[1] + 1;
					connectivity[3] = connectivity[1] + numNodes[0] + 1;
					connectivity[4] = connectivity[1] + numNodes[0];
					connectivity[5] = connectivity[1] + numNodes[0] * numNodes[1];
					connectivity[6] = connectivity[5] + 1;
					connectivity[7] = connectivity[5] + numNodes[0] + 1;
					connectivity[8] = connectivity[5] + numNodes[0];

					// エンディアン変換して書き込み
					for (auto &value : connectivity)
					{
						swapEndian(value);
					}
					file.write(reinterpret_cast<char *>(connectivity.data()), connectivity.size() * sizeof(int32_t));
				}
	}

	// セルタイプの出力
	file << "\nCELL_TYPES " << numSubElements << "\n";
	int32_t cellType = 12; // 六面体要素
	swapEndian(cellType);
	for (int i = 0; i < numSubElements; ++i)
	{
		file.write(reinterpret_cast<char *>(&cellType), sizeof(int32_t));
	}

	// スカラー値の出力
	file << "\nPOINT_DATA " << totalPoints << "\n";
	file << "SCALARS value float 1\n";
	file << "LOOKUP_TABLE default\n";

	// 各点でのスカラー値を計算
	std::vector<float> values(totalPoints);
	#pragma omp parallel for schedule(dynamic) collapse(3)
		for (int elem = info->Total_Element_to_mesh[1]; elem < info->Total_Element_to_mesh[Total_mesh]; ++elem)
			for (int k = 0; k < numNodes[2]; ++k)
				for (int j = 0; j < numNodes[1]; ++j)
					for (int i = 0; i < numNodes[0]; ++i)
					{
						double para[3] = {subdivCoords[0][i], subdivCoords[1][j], subdivCoords[2][k]};
						double value = blendingFunction(elem, info->Element_patch[elem], para, info);
						int index = (elem - info->Total_Element_to_mesh[1]) * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i;
						values[index] = static_cast<float>(value);
					}

	// バイナリ形式でスカラー値を書き込み
	for (auto &value : values)
	{
		swapEndian(value);
	}
	file.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(float));

	file.close();
}


void debugBlendingFunctionAux(information *info, int itr)
{
	information *aux = info->aux.get();

	string filename = "blending_aux";
	filename += to_string(itr) + ".vtk";

	const int nelem = aux->Total_Element_to_mesh[Total_mesh] - aux->Total_Element_to_mesh[1];
	const int subdivision[3] = {1, 1, 1};
	const int numNodes[3] = {subdivision[0] + 1, subdivision[1] + 1, subdivision[2] + 1};
	const int nodesPerElement = numNodes[0] * numNodes[1] * numNodes[2];

	// パラメトリック座標の分割点を計算
	std::vector<std::vector<double>> subdivCoords(3);
	for (int dir = 0; dir < 3; ++dir)
	{
		subdivCoords[dir].resize(numNodes[dir]);
		for (int i = 0; i < numNodes[dir]; ++i)
		{
			subdivCoords[dir][i] = -1.0 + (2.0 * i) / (double)subdivision[dir];
		}
	}

	std::ofstream file(filename, std::ios::binary);

	// ヘッダー部分
	file << "# vtk DataFile Version 2.0\n";
	file << "Subdivided hexahedral mesh data\n";
	file << "BINARY\n";
	file << "DATASET UNSTRUCTURED_GRID\n";

	// 点の座標を出力
	const int totalPoints = nelem * nodesPerElement;
	file << "POINTS " << totalPoints << " float\n";

	// 各点の座標を計算
	std::vector<float> coordinates(totalPoints * 3);
	#pragma omp parallel for schedule(dynamic) collapse(3)
	for (int elem = aux->Total_Element_to_mesh[1]; elem < aux->Total_Element_to_mesh[Total_mesh]; ++elem)
		for (int k = 0; k < numNodes[2]; ++k)
			for (int j = 0; j < numNodes[1]; ++j)
				for (int i = 0; i < numNodes[0]; ++i)
				{
					double para[3] = {subdivCoords[0][i],
									  subdivCoords[1][j],
									  subdivCoords[2][k]};
					double xyz[3] = {0.0};
					vector<double> R(MAX_NO_CP_ON_ELEMENT);
					shape_and_dshape(R.data(), para, elem, aux);
					for (int l = 0; l < aux->No_Control_point_ON_ELEMENT[aux->Element_patch[elem]]; l++)
						for (int m = 0; m < aux->DIMENSION; m++)
						{
							int cp = aux->Controlpoint_of_Element[elem * MAX_NO_CP_ON_ELEMENT + l];
							xyz[m] += R[l] * (aux->Node_Coordinate[cp * (aux->DIMENSION + 1) + m] + aux->disp_overlay[cp * aux->DIMENSION + m]);
						}
					int index = ((elem - aux->Total_Element_to_mesh[1]) * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i) * 3;
					coordinates[index] = static_cast<float>(xyz[0]);
					coordinates[index + 1] = static_cast<float>(xyz[1]);
					coordinates[index + 2] = static_cast<float>(xyz[2]);
				}

	// バイナリ形式で座標を書き込み
	for (auto &coord : coordinates)
	{
		swapEndian(coord); // エンディアン変換
	}
	file.write(reinterpret_cast<char *>(coordinates.data()), coordinates.size() * sizeof(float));

	// セルの接続情報を出力
	const int numSubElements = (subdivision[0]) * (subdivision[1]) * (subdivision[2]) * nelem;
	file << "\nCELLS " << numSubElements << " " << (numSubElements * 9) << "\n";

	std::vector<int32_t> connectivity(9); // セル当たりの接続情報
	for (int elem = 0; elem < aux->Total_Element_to_mesh[Total_mesh] - aux->Total_Element_to_mesh[1]; ++elem)
	{
		const int elemOffset = elem * nodesPerElement;
		for (int k = 0; k < subdivision[2]; ++k)
			for (int j = 0; j < subdivision[1]; ++j)
				for (int i = 0; i < subdivision[0]; ++i)
				{
					// 1つのサブ要素の8つの節点番号を計算
					connectivity[0] = 8; // 頂点数
					connectivity[1] = elemOffset + i + j * numNodes[0] + k * numNodes[0] * numNodes[1];
					connectivity[2] = connectivity[1] + 1;
					connectivity[3] = connectivity[1] + numNodes[0] + 1;
					connectivity[4] = connectivity[1] + numNodes[0];
					connectivity[5] = connectivity[1] + numNodes[0] * numNodes[1];
					connectivity[6] = connectivity[5] + 1;
					connectivity[7] = connectivity[5] + numNodes[0] + 1;
					connectivity[8] = connectivity[5] + numNodes[0];

					// エンディアン変換して書き込み
					for (auto &value : connectivity)
					{
						swapEndian(value);
					}
					file.write(reinterpret_cast<char *>(connectivity.data()), connectivity.size() * sizeof(int32_t));
				}
	}

	// セルタイプの出力
	file << "\nCELL_TYPES " << numSubElements << "\n";
	int32_t cellType = 12; // 六面体要素
	swapEndian(cellType);
	for (int i = 0; i < numSubElements; ++i)
	{
		file.write(reinterpret_cast<char *>(&cellType), sizeof(int32_t));
	}

	// スカラー値の出力
	file << "\nPOINT_DATA " << totalPoints << "\n";
	file << "SCALARS value float 1\n";
	file << "LOOKUP_TABLE default\n";

	// 各点でのスカラー値を計算
	std::vector<float> values(totalPoints);
	#pragma omp parallel for schedule(dynamic) collapse(3)
		for (int elem = aux->Total_Element_to_mesh[1]; elem < aux->Total_Element_to_mesh[Total_mesh]; ++elem)
			for (int k = 0; k < numNodes[2]; ++k)
				for (int j = 0; j < numNodes[1]; ++j)
					for (int i = 0; i < numNodes[0]; ++i)
					{
						double para[3] = {subdivCoords[0][i], subdivCoords[1][j], subdivCoords[2][k]};
						double value = blendingFunction(elem, aux->Element_patch[elem], para, aux);
						int index = (elem - aux->Total_Element_to_mesh[1]) * nodesPerElement + k * numNodes[1] * numNodes[0] + j * numNodes[0] + i;
						values[index] = static_cast<float>(value);
					}

	// バイナリ形式でスカラー値を書き込み
	for (auto &value : values)
	{
		swapEndian(value);
	}
	file.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(float));

	file.close();
}


template <typename T>
void swapEndian(T &val)
{
	char *bytes = reinterpret_cast<char *>(&val);
	for (size_t i = 0; i < sizeof(T) / 2; ++i)
	{
		std::swap(bytes[i], bytes[sizeof(T) - 1 - i]);
	}
}


bool fictitiousAreaCheck(int e, int p, double *para, information *info)
{
	// if local patch
	if (info->Total_Patch_to_mesh[1] <= p)
		return false;
	
	// check overlap
	bool isOverlap = false;
	double coord[MAX_DIMENSION] = {0.0};
	double temp_para[MAX_DIMENSION] = {0.0};
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	shape_and_dshape(R.data(), para, e, info);
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[p]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			coord[j] += R[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	int itr_n = 0;
	for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		itr_n = calc_patch_parameter_coord(coord, i, temp_para, info);
		if (itr_n != ERROR)
		{
			isOverlap = true;
			break;
		}
	}
	if (isOverlap)
		return false;

	// check fictitious area
	for (size_t i = 0; i < info->fa.size(); i++)
	{
		fictitious_area &fa = info->fa[i];

		// check patch
		if (fa.patch != p)
			continue;
		
		double knot[MAX_DIMENSION];
		double knot_start[MAX_DIMENSION];
		double knot_end[MAX_DIMENSION];
		double knot_diff[MAX_DIMENSION];
		for (int j = 0; j < info->DIMENSION; j++)
		{
			knot_start[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + j] + info->Order[p * info->DIMENSION + j] + info->ENC[e * info->DIMENSION + j]];
			knot_end[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + j] + info->Order[p * info->DIMENSION + j] + info->ENC[e * info->DIMENSION + j] + 1];
			knot_diff[j] = knot_end[j] - knot_start[j];
			knot[j] = knot_start[j] + 0.5 * knot_diff[j] * (para[j * info->DIMENSION + j] + 1.0);
		}

		bool isInside = true;
		for (int j = 0; j < info->DIMENSION; j++)
		{
			if (fabs(fa.area_start[j] - knot[j]) < MERGE_ERROR || (fa.area_start[j] <= knot[j] && knot[j] <= fa.area_end[j]) || fabs(fa.area_end[j] - knot[j]) < MERGE_ERROR)
				isInside &= true;
			else
			{
				isInside &= false;
				break;
			}
		}
		if (isInside)
			return true;
	}
	return false;
}


// B-bar method using L2 projection
void get_decreased_order_data(information *info, char *option_file)
{
	info->dec_order= std::make_unique<information>();
	information *dec_order= info->dec_order.get();

	dec_order->isDecreasedOrderPatch = true;

	Get_Option(option_file, dec_order->opt_files);
	dec_order->c.get_data(dec_order->opt_files[0].c_str(), dec_order);

	vector<string> filename;
	if (!dec_order->opt_files[6].empty())
		filename.emplace_back(dec_order->opt_files[6]);
	if (!dec_order->opt_files[7].empty())
		filename.emplace_back(dec_order->opt_files[7]);

	Allocation(0, dec_order);
	for (int i = 0; i < static_cast<int>(filename.size()); i++)
		Get_Input_1(i, filename[i].c_str(), dec_order);
	Allocation(1, dec_order);
	for (int i = 0; i < static_cast<int>(filename.size()); i++)
		Get_Input_2(i, filename[i].c_str(), dec_order);
	Allocation(2, dec_order);
	Make_INC(dec_order);
	Allocation(3, dec_order);
	searchOverlappingEle(dec_order);
	Allocation(4, dec_order);

	dec_order->D = info->D;
	dec_order->Index_Dof = info->Index_Dof;
	dec_order->K_Whole_Val = info->K_Whole_Val;
	dec_order->K_Whole_Ptr = info->K_Whole_Ptr;
	dec_order->K_Whole_Col = info->K_Whole_Col;
	dec_order->full_K = info->full_K;
	dec_order->rhs_vec = info->rhs_vec;

	get_material(dec_order);

	// set gauss points
	Make_Gauss_points(false, dec_order);

	Allocation_init_NL(dec_order);

	dec_order->internal_force = info->internal_force;
	dec_order->rhs_vec_initial = info->rhs_vec_initial;

	// make Matrix for B-bar
	dec_order->L2_M.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	dec_order->L2_M_ptr.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	dec_order->L2_M_col.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	
	vector<vector<set<int>>> L2_ptr(dec_order->Total_Patch_to_mesh[Total_mesh]);
	vector<vector<std::unique_ptr<std::mutex>>> mtx(dec_order->Total_Patch_to_mesh[Total_mesh]);

	// initialize L2_M_ptr and L2_ptr
	#pragma omp parallel for
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		dec_order->L2_M_ptr[i].resize(dec_order->No_Control_point_in_patch[i] + 1, 0);
		L2_ptr[i].resize(dec_order->No_Control_point_in_patch[i]);
		mtx[i].resize(dec_order->No_Control_point_in_patch[i]);
		for (size_t j = 0; j < mtx[i].size(); j++)
			mtx[i][j] = std::make_unique<std::mutex>();
	}

	// make L2_M ptr
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dec_order->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int p = dec_order->Element_patch[i];
		for (int j = 0; j < dec_order->No_Control_point_ON_ELEMENT[p]; j++)
		{
			int row = dec_order->Controlpoint_of_Element_in_patch[i * MAX_NO_CP_ON_ELEMENT + j];
			for (int k = 0; k < dec_order->No_Control_point_ON_ELEMENT[p]; k++)
			{
				int col = dec_order->Controlpoint_of_Element_in_patch[i * MAX_NO_CP_ON_ELEMENT + k];
				if (col >= row)
				{
					std::lock_guard<std::mutex> lock(*mtx[p][row]);
					L2_ptr[p][row].insert(col);
				}
			}
		}
	}

	// substitute L2_ptr to L2_M_ptr
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		dec_order->L2_M_ptr[i][0] = 0;
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
			dec_order->L2_M_ptr[i][j + 1] = dec_order->L2_M_ptr[i][j] + static_cast<int>(L2_ptr[i][j].size());
	}

	// initialize L2_M and L2_M_col
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		dec_order->L2_M[i].resize(dec_order->L2_M_ptr[i][dec_order->No_Control_point_in_patch[i]], 0.0);
		dec_order->L2_M_col[i].resize(dec_order->L2_M_ptr[i][dec_order->No_Control_point_in_patch[i]], 0);
	}

	// make L2_M col
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
		{
			int count = 0;
			for (auto it = L2_ptr[i][j].begin(); it != L2_ptr[i][j].end(); ++it)
			{
				dec_order->L2_M_col[i][dec_order->L2_M_ptr[i][j] + count] = *it;
				count++;
			}
		}
	}

	// initialize L2_G and L2_M_inv_G
	#if 0
	dec_order->L2_G.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	dec_order->L2_M_inv_G.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		dec_order->L2_G[i].resize(dec_order->No_Control_point_in_patch[i]);
		dec_order->L2_M_inv_G[i].resize(dec_order->No_Control_point_in_patch[i]);
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
		{
			dec_order->L2_G[i][j].resize(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0);
			dec_order->L2_M_inv_G[i][j].resize(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0);
		}
	}
	#else
	dec_order->L2_G.resize(dec_order->Total_Patch_to_mesh[Total_mesh]);
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
		dec_order->L2_G[i].resize(info->No_Control_point_in_patch[i] * info->DIMENSION);
	#endif

	// calculate L2_M and L2_G values at each element
	auto Calc_L2_M_ele_L2_G_ele = [](int e, vector<double> &L2_M_ele, vector<vector<double>> &L2_G_ele, information *dec_order, information *info) {
		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
		int p = dec_order->Element_patch[e];

		// gauss point
		for (int i = 0; i < info->gp[e].n(); i++)
		{
			// make R
			shape_and_dshape(R.data(), &info->gp[e].para()[i * info->DIMENSION], e, dec_order);

			// calc L2_G
			double J = Make_B_component(e, &info->gp[e].para()[i * info->DIMENSION], info->gp[e].isOverlay()[i], info->gp[e].opp_ele()[i], &info->gp[e].opp_para_tilde()[i * info->DIMENSION], b.data(), info);
			double coef = info->gp[e].w()[i] * J;
			for (int j = 0; j < dec_order->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						L2_G_ele[j][k * info->DIMENSION + l] += R[j] * b[l * MAX_NO_CP_ON_ELEMENT + k] * coef;

			// calc L2_M
			for (int j = 0; j < dec_order->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < dec_order->No_Control_point_ON_ELEMENT[p]; k++)
					L2_M_ele[j * dec_order->No_Control_point_ON_ELEMENT[p] + k] += R[j] * R[k] * coef;
		}
	};

	// init L2_G
	for (int p = 0; p < dec_order->Total_Patch_to_mesh[Total_mesh]; p++)
		for (size_t col = 0; col < dec_order->L2_G[p].size(); col++)
			dec_order->L2_G[p][col].init_thread_local();

	// make L2_M and L2_G
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dec_order->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int p = dec_order->Element_patch[i];
		vector<double> L2_M_ele(dec_order->No_Control_point_ON_ELEMENT[p] * dec_order->No_Control_point_ON_ELEMENT[p], 0.0);
		vector<vector<double>> L2_G_ele(dec_order->No_Control_point_ON_ELEMENT[p], vector<double>(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0));
		Calc_L2_M_ele_L2_G_ele(i, L2_M_ele, L2_G_ele, dec_order, info);
		for (int j = 0; j < dec_order->No_Control_point_ON_ELEMENT[p]; j++)
		{
			int row = dec_order->Controlpoint_of_Element_in_patch[i * MAX_NO_CP_ON_ELEMENT + j];

			// L2_G
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
			{
				for (int l = 0; l < info->DIMENSION; l++)
				{
					int col_loc = k * info->DIMENSION + l;
					int col = info->Controlpoint_of_Element_in_patch[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
					double val = L2_G_ele[j][col_loc];
					dec_order->L2_G[p][col].insert(row, val);
				}
			}

			// L2_M
			for (int k = 0; k < dec_order->No_Control_point_ON_ELEMENT[p]; k++)
			{
				int col = dec_order->Controlpoint_of_Element_in_patch[i * MAX_NO_CP_ON_ELEMENT + k];
				if (col >= row)
				{
					int icount = L2_M_RowCol_to_icount(p, row, col, dec_order);
					double val = L2_M_ele[j * dec_order->No_Control_point_ON_ELEMENT[p] + k];
					#pragma omp atomic
					dec_order->L2_M[p][icount] += val;
				}
			}
		}
	}

	// merge L2_G entries
	#pragma omp parallel for schedule(dynamic)
	for (int p = 0; p < dec_order->Total_Patch_to_mesh[Total_mesh]; p++)
		for (size_t col = 0; col < dec_order->L2_G[p].size(); col++)
			dec_order->L2_G[p][col].compress();

	// calc [L2_M]^(-1){L2_G} and make stiffness matrix
	cout << endl << "L2 projection: Start making 'L2_M_inv_G'" << endl;
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
		intel_PARDISO_for_L2(i, dec_order->No_Control_point_in_patch[i], info, dec_order);
	cout << "Finish preprocess of L2 projection" << endl << endl;

	#if 0
	// calc {L2_G}^T[L2_M]^(-1){L2_G} and search non-zero pattern
	vector<vector<set<int>>> L2_M_inv_G_T_L2_M_inv_G_ptr_set(dec_order->Total_Patch_to_mesh[Total_mesh]);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		L2_M_inv_G_T_L2_M_inv_G_ptr_set[i].resize(dec_order->No_Control_point_in_patch[i] + 1, 0);
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
		{
			for (int k = 0; k < dec_order->No_Control_point_in_patch[i]; k++)
			{
				double val = 0.0;
				for (int l = 0; l < MAX_NO_CP_ON_ELEMENT * info->DIMENSION; l++)
					val += dec_order->L2_M_inv_G[i][j][l] * dec_order->L2_G[i][k][l];
				if (fabs(val) > MERGE_ERROR)
				{
					#pragma omp critical
					L2_M_inv_G_T_L2_M_inv_G_ptr_set[i][j].insert(k);
				}
			}
		}
	}
	// substitute L2_M_inv_G_T_L2_M_inv_G_ptr_set to L2_M_inv_G_T_L2_M_inv_G_ptr
	vector<vector<int>> L2_M_inv_G_T_L2_M_inv_G_ptr(dec_order->Total_Patch_to_mesh[Total_mesh]);
	#pragma omp parallel for
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		L2_M_inv_G_T_L2_M_inv_G_ptr[i].resize(dec_order->No_Control_point_in_patch[i] + 1, 0);
		L2_M_inv_G_T_L2_M_inv_G_ptr[i][0] = 0;
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
			L2_M_inv_G_T_L2_M_inv_G_ptr[i][j + 1] = L2_M_inv_G_T_L2_M_inv_G_ptr[i][j] + static_cast<int>(L2_M_inv_G_T_L2_M_inv_G_ptr_set[i].count(j));
	}
	vector<vector<double>> L2_M_inv_G_T_L2_M_inv_G_col(dec_order->Total_Patch_to_mesh[Total_mesh]);
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		L2_M_inv_G_T_L2_M_inv_G_col[i].resize(L2_M_inv_G_T_L2_M_inv_G_ptr[i][dec_order->No_Control_point_in_patch[i]], 0);
		#pragma omp parallel for
		for (int j = 0; j < dec_order->No_Control_point_in_patch[i]; j++)
		{
			int count = 0;
			for (auto it = L2_M_inv_G_T_L2_M_inv_G_ptr_set[i].begin(); it != L2_M_inv_G_T_L2_M_inv_G_ptr_set[i].end(); ++it)
			{
				if (*it == j)
				{
					L2_M_inv_G_T_L2_M_inv_G_col[i][L2_M_inv_G_T_L2_M_inv_G_ptr[i][j] + count] = *it;
					count++;
				}
			}
		}
	}

	// debug L2_M_inv_G_T_L2_M_inv_G_ptr
	for (int i = 0; i < dec_order->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		cout << "Patch " << i << " L2_M_inv_G_T_L2_M_inv_G_ptr: ";
		for (size_t j = 0; j < L2_M_inv_G_T_L2_M_inv_G_ptr[i].size(); j++)
			cout << L2_M_inv_G_T_L2_M_inv_G_ptr[i][j] << " ";
		cout << endl;
	}
	exit(0);
	#endif

	// check elements
	if (info->Total_Element_to_mesh[Total_mesh] != dec_order->Total_Element_to_mesh[Total_mesh])
	{
		cerr << "Error: The number of elements in the decreased order model does not match." << endl;
		exit(1);
	}
	else
	{
		cout << "\033[33m";
		cout << "The elements of the decreased order model have been successfully checked." << endl;
		cout << "\033[0m";
	}
}


int L2_M_RowCol_to_icount(int &p, int &row, int &col, information *dec_order)
{
	// 二分探索
	int min_id = dec_order->L2_M_ptr[p][row], max_id = dec_order->L2_M_ptr[p][row + 1] - 1, mid_id;
	while (min_id <= max_id)
	{
		mid_id = (min_id + max_id) / 2;

		if (dec_order->L2_M_col[p][mid_id] == col)
			return mid_id;
		else if (dec_order->L2_M_col[p][mid_id] < col)
			min_id = mid_id + 1;
		else
			max_id = mid_id - 1;
	}

	// error
	cerr << "Error: L2_M_RowCol_to_icount did not find the column." << endl;
	exit(1);

	return 0;
}


int L2_debug_RowCol_to_icount(vector<int> &M_ptr_debug, vector<int> &M_col_debug, int &row, int &col)
{
	// 二分探索
	int min_id = M_ptr_debug[row], max_id = M_ptr_debug[row + 1] - 1, mid_id;
	while (min_id <= max_id)
	{
		mid_id = (min_id + max_id) / 2;

		if (M_col_debug[mid_id] == col)
			return mid_id;
		else if (M_col_debug[mid_id] < col)
			min_id = mid_id + 1;
		else
			max_id = mid_id - 1;
	}

	// error
	cerr << "Error: L2_debug_RowCol_to_icount did not find the column." << endl;
	exit(1);

	return 0;
}


// Update Control Point
int Update_Control_Point_glo_old(information *info)
{
	// グローバルパッチ上のコントロールポイントに対応するノットベクトル, 要素を算出
	int CP_glo_n = info->Total_Control_Point_to_mesh[1];
	static int *ele_at_cp = (int *)malloc(sizeof(int) * CP_glo_n);
	static double *knot_vec_at_cp = (double *)malloc(sizeof(double) * CP_glo_n * info->DIMENSION);
	int temp_max = 0;
	for (int i = 0; i < info->Total_Patch_to_mesh[1]; i++)
	{
		int temp = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			temp *= info->No_Control_point[i * info->DIMENSION + j];
		if (temp_max < temp)
			temp_max = temp;
	}
	static double *knot_array = (double *)malloc(sizeof(double) * temp_max * info->DIMENSION);

	static int counter = 0;
	if (counter == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Patch_to_mesh[1]; i++)
		{
			// Calc
			Calc_MI3D(info->DIMENSION, &info->Order[i * info->DIMENSION], &info->No_Control_point[i * info->DIMENSION], knot_array);

			// Update
			int cp_to_patch = info->Total_Control_Point_to_patch[i];
			for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
			{
				int connectivity = info->Patch_Control_point[cp_to_patch + j];
				ele_at_cp[connectivity] = ele_check(i, &knot_array[j * info->DIMENSION], info);
				tilde_coord(&knot_vec_at_cp[connectivity * info->DIMENSION], &knot_array[j * info->DIMENSION], i, ele_at_cp[connectivity], info);
			}
		}
		counter++;
	}

	// cout << "glo" << endl;
	// for (int i = 0; i < CP_glo_n; i++)
	// {
	// 	cout << i << " ";
	// 	cout << ele_at_cp[i] << " ";
	// 	for (int j = 0; j < info->DIMENSION; j++)
	// 	{
	// 		printf("%.7e ", knot_vec_at_cp[i * info->DIMENSION + j]);
	// 	}
	// 	cout << endl;
	// }
	// for (int i = 0; i < info->Total_Patch_to_mesh[1]; i++)
	// {
	// 	for (int j = 0; j < info->DIMENSION; j++)
	// 	{
	// 		for (int k = 0; k < info->No_Control_point_in_patch[i]; k++)

	// 	}
	// }

	// newton raphson
	static const int max_itr = info->c.UPDATE_CONTROL_POINT_MAX_ITR;
	vector<double> delta_sol(CP_glo_n * info->DIMENSION, 0.0);
	// double *delta_sol = (double *)calloc(CP_glo_n * info->DIMENSION, sizeof(double));
	static double *right_vec = (double *)malloc(sizeof(double) * CP_glo_n * info->DIMENSION);
	static double *true_f_coord = (double *)malloc(sizeof(double) * CP_glo_n * info->DIMENSION);
	static double *temp_f = (double *)malloc(sizeof(double) * CP_glo_n * info->DIMENSION);

	static double *glo_coordinate = (double *)malloc(sizeof(double) * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION);
	size_t M_size = CP_glo_n * info->DIMENSION;
	static double *M_on_point = (double *)malloc(sizeof(double) * M_size * M_size);
	// static double *M_small = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
	// static double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);

	static double *temp_glo_coordinate = (double *)malloc(sizeof(double) * CP_glo_n * info->DIMENSION);
	double temp_min = 0.0;

	// make true_f
	#pragma omp parallel for
	for (int i = 0; i < CP_glo_n; i++)
	{
		int e = ele_at_cp[i];

		int element_loc = 0;
		double temp_para_loc[MAX_DIMENSION];
		double temp_point_loc[MAX_DIMENSION];
		double temp_coord[MAX_DIMENSION];

		// init
		for (int j = 0; j < info->DIMENSION; j++)
		{
			true_f_coord[i * info->DIMENSION + j] = 0.0;
			temp_coord[j] = 0.0;
		}

		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), &knot_vec_at_cp[i * info->DIMENSION], e, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
				double d = info->disp[id] + info->disp_increment[id];
				double coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				temp_coord[k] += R[j] * coord;
				true_f_coord[i * info->DIMENSION + k] += R[j] * (coord + d);
			}

		int status = 0;

		// overlay
		int itr_n = 0, loc_patch = 0;
		for (int j = info->Total_Patch_to_mesh[1]; j < info->Total_Patch_to_mesh[Total_mesh]; j++)
		{
			itr_n = calc_patch_parameter_coord(temp_coord, j, temp_para_loc, info);
			loc_patch = j;
			if (itr_n != ERROR)
			{
				status = 1;
				break;
			}
		}

		if (status)
		{
			element_loc = ele_check(loc_patch, temp_para_loc, info);
			tilde_coord(temp_point_loc, temp_para_loc, loc_patch, element_loc, info);
		}

		if (status)
		{
			vector<double> R_loc(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int id = info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
					double d_loc = info->disp[id] + info->disp_increment[id];
					true_f_coord[i * info->DIMENSION + k] += R_loc[j] * d_loc;
				}
		}
	}

	// init M_on_point
	#pragma omp parallel for
	for (size_t i = 0; i < M_size * M_size; i++)
	{
		M_on_point[i] = 0.0;
	}

	// newton raphson loop
	int check_convergence = 0;
	int count_NR = 0;
	for (int i = 0; i < max_itr; i++)
	{
		// update glo_coordinate
		if (i == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_glo_n; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					// glo_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[j * (info->DIMENSION + 1) + k];
					glo_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[j * (info->DIMENSION + 1) + k] + info->disp_overlay[j * info->DIMENSION + k] + info->disp_overlay_increment[j * info->DIMENSION + k];
				}
		}
		else
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				glo_coordinate[j] += delta_sol[j];
		}

		// make right_vec
		#pragma omp parallel for
		for (int j = 0; j < CP_glo_n; j++)
		{
			int e = ele_at_cp[j];

			// init
			for (int k = 0; k < info->DIMENSION; k++)
				temp_f[j * info->DIMENSION + k] = 0.0;

			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &knot_vec_at_cp[j * info->DIMENSION], e, info);
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
				{
					int ptr = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
					temp_f[j * info->DIMENSION + l] += R[k] * glo_coordinate[ptr];

					// make M_on_point
					if (i == 0)
					{
						size_t index = (j * info->DIMENSION + l) * M_size + ptr;
						double val = R[k];

						#pragma omp atomic
						M_on_point[index] += val;
					}
				}
		}

		#pragma omp parallel for
		for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
			right_vec[j] = true_f_coord[j] - temp_f[j];

		// solve
		// GaussianElimination(delta_sol, right_vec, M_on_point, M_size);
		// temp_solver(delta_sol, right_vec, M_on_point, M_size, M_size);
		sparseLU_update_CP(delta_sol.data(), right_vec, M_on_point, M_size, M_size);

		// check
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
		{
			temp1 += pow(true_f_coord[j] - temp_f[j], 2);
			temp2 += pow(true_f_coord[j], 2);
			temp3 += pow(delta_sol[j], 2);
		}
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);
		temp3 = sqrt(temp3);
		double temp_norm = temp1 / temp2;
		printf("\titr_norm: %d %.6e\n", i, temp_norm);
		printf("\tdelta_sol_norm: %d %.6e\n", i, temp3);
		if (temp_norm < MERGE_ERROR)
		{
			check_convergence = 1;
			// update
			for (int i = 0; i < CP_glo_n * info->DIMENSION; i++)
				glo_coordinate[i] += delta_sol[i];
			for (int i = 0; i < CP_glo_n; i++)
				for (int j = 0; j < info->DIMENSION; j++)
					info->disp_overlay_increment[i * info->DIMENSION + j] = glo_coordinate[i * info->DIMENSION + j] - info->Node_Coordinate[i * (info->DIMENSION + 1) + j] - info->disp_overlay[i * info->DIMENSION + j];
			break;
		}

		if (i == 0)
			temp_min = temp_norm;
		else if (temp_norm < temp_min)
		{
			count_NR = 0;
			temp_min = temp_norm;
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				temp_glo_coordinate[j] = glo_coordinate[j] + delta_sol[j];
		}
		else
			count_NR++;

		if (temp_norm > temp_min * 1.0e+2 && count_NR > 10)
		{
			// update
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				glo_coordinate[j] = temp_glo_coordinate[j];
			break;
		}
	}

	if (check_convergence == 0)
		return 1;
	else if (check_convergence == 1)
		return 0;
	return 0;
}


int Update_Control_Point_loc_old(information *info)
{
	// ローカルパッチ上のコントロールポイントに対応するノットベクトル, 要素を算出
	int CP_loc_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	static int *ele_at_cp = (int *)malloc(sizeof(int) * CP_loc_n);
	static double *knot_vec_at_cp = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	int temp_max = 0;
	for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		int temp = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			temp *= info->No_Control_point[i * info->DIMENSION + j];
		if (temp_max < temp)
			temp_max = temp;
	}
	static double *knot_array = (double *)malloc(sizeof(double) * temp_max * info->DIMENSION);

	static int counter = 0;
	if (counter == 0)
	{
		for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			// Calc
			Calc_MI3D(info->DIMENSION, &info->Order[i * info->DIMENSION], &info->No_Control_point[i * info->DIMENSION], knot_array);

			// Update
			int cp_to_patch = info->Total_Control_Point_to_patch[i];
			for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
			{
				int connectivity = info->Patch_Control_point[cp_to_patch + j] - info->Total_Control_Point_to_mesh[1];
				ele_at_cp[connectivity] = ele_check(i, &knot_array[j * info->DIMENSION], info);
				tilde_coord(&knot_vec_at_cp[connectivity * info->DIMENSION], &knot_array[j * info->DIMENSION], i, ele_at_cp[connectivity], info);
			}
		}
		counter++;
	}

	// cout << "loc" << endl;
	// for (int i = 0; i < CP_loc_n; i++)
	// {
	// 	cout << i << " ";
	// 	cout << ele_at_cp[i] << " ";
	// 	for (int j = 0; j < info->DIMENSION; j++)
	// 	{
	// 		printf("%.7e ", knot_vec_at_cp[i * info->DIMENSION + j]);
	// 	}
	// 	cout << endl;
	// }

	// newton raphson
	static const int max_itr = info->c.UPDATE_CONTROL_POINT_MAX_ITR;
	vector<double> delta_sol(CP_loc_n * info->DIMENSION, 0.0);
	// double *delta_sol = (double *)calloc(CP_loc_n * info->DIMENSION, sizeof(double));
	static double *right_vec = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	static double *true_f_coord = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	static double *temp_f = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);

	static double *loc_coordinate = (double *)malloc(sizeof(double) * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION);
	size_t M_size = CP_loc_n * info->DIMENSION;
	static double *M_on_point = (double *)malloc(sizeof(double) * M_size * M_size);
	// static double *M_small = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
	// static double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);

	static double *temp_loc_coordinate = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	double temp_min = 0.0;

	// make true_f
	#pragma omp parallel for
	for (int i = 0; i < CP_loc_n; i++)
	{
		int e = ele_at_cp[i];

		int element_glo = 0;
		double temp_para_glo[MAX_DIMENSION];
		double temp_point_glo[MAX_DIMENSION];
		double temp_coord[MAX_DIMENSION];

		// init
		for (int j = 0; j < info->DIMENSION; j++)
		{
			true_f_coord[i * info->DIMENSION + j] = 0.0;
			temp_coord[j] = 0.0;
		}

		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), &knot_vec_at_cp[i * info->DIMENSION], e, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
				double d = info->disp[id] + info->disp_increment[id];
				double coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				temp_coord[k] += R[j] * coord;
				true_f_coord[i * info->DIMENSION + k] += R[j] * (coord + d);
			}

		// overlay
		int itr_n = 0, glo_patch = 0;
		for (int j = 0; j < info->Total_Patch_to_mesh[1]; j++)
		{
			itr_n = calc_patch_parameter_coord(temp_coord, j, temp_para_glo, info);
			glo_patch = j;
			if (itr_n != ERROR)
				break;
		}

		element_glo = ele_check(glo_patch, temp_para_glo, info);
		tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

		vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
				double d_glo = info->disp[id] + info->disp_increment[id];
				true_f_coord[i * info->DIMENSION + k] += R_glo[j] * d_glo;
			}
	}

	// init M_on_point
	#pragma omp parallel for
	for (size_t i = 0; i < M_size * M_size; i++)
	{
		M_on_point[i] = 0.0;
	}

	// newton raphson loop
	int check_convergence = 0;
	int count_NR = 0;
	for (int i = 0; i < max_itr; i++)
	{
		// update loc_coordinate
		if (i == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_loc_n; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					// loc_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[(j + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + k];
					int index = (j + info->Total_Control_Point_to_mesh[1]);
					loc_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[index * (info->DIMENSION + 1) + k] + info->disp_overlay[index * info->DIMENSION + k] + info->disp_overlay_increment[index * info->DIMENSION + k];
				}
		}
		else
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] += delta_sol[j];
		}

		// make right_vec
		#pragma omp parallel for
		for (int j = 0; j < CP_loc_n; j++)
		{
			int e = ele_at_cp[j];

			// init
			for (int k = 0; k < info->DIMENSION; k++)
				temp_f[j * info->DIMENSION + k] = 0.0;

			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &knot_vec_at_cp[j * info->DIMENSION], e, info);
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
				{
					int ptr = (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + l;
					temp_f[j * info->DIMENSION + l] += R[k] * loc_coordinate[ptr];

					// make M_on_point
					if (i == 0)
					{
						size_t index = (j * info->DIMENSION + l) * M_size + ptr;
						double val = R[k];

						#pragma omp atomic
						M_on_point[index] += val;
					}
				}
		}

		#pragma omp parallel for
		for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
			right_vec[j] = true_f_coord[j] - temp_f[j];

		// solve
		// GaussianElimination(delta_sol, right_vec, M_on_point, M_size);
		// temp_solver(delta_sol, right_vec, M_on_point, M_size, M_size);
		sparseLU_update_CP(delta_sol.data(), right_vec, M_on_point, M_size, M_size);

		// check
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
		{
			temp1 += pow(true_f_coord[j] - temp_f[j], 2);
			temp2 += pow(true_f_coord[j], 2);
			temp3 += pow(delta_sol[j], 2);
		}
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);
		temp3 = sqrt(temp3);
		double temp_norm = temp1 / temp2;
		printf("\titr_norm: %d %.6e\n", i, temp_norm);
		printf("\tdelta_sol_norm: %d %.6e\n", i, temp3);
		if (temp_norm < MERGE_ERROR)
		{
			check_convergence = 1;
			// update
			for (int i = 0; i < CP_loc_n * info->DIMENSION; i++)
				loc_coordinate[i] += delta_sol[i];
			for (int i = 0; i < CP_loc_n; i++)
				for (int j = 0; j < info->DIMENSION; j++)
					info->disp_overlay_increment[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j] = loc_coordinate[i * info->DIMENSION + j] - info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + j] - info->disp_overlay[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j];
			break;
		}

		if (i == 0)
			temp_min = temp_norm;
		if (temp_norm < temp_min)
		{
			count_NR = 0;
			temp_min = temp_norm;
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				temp_loc_coordinate[j] = loc_coordinate[j] + delta_sol[j];
		}
		else
			count_NR++;

		if (temp_norm > temp_min * 1.0e+2 && count_NR > 10)
		{
			// update
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] = temp_loc_coordinate[j];
			break;
		}
	}

	if (check_convergence == 0)
		return 1;
	else if (check_convergence == 1)
		return 0;
	return 0;
}


void Update_Control_Point_substitute_glo_aux(information *info, information *aux)
{
	for (int i = 0; i < info->Total_Control_Point_to_mesh[1] * info->DIMENSION; i++)
		aux->disp_overlay_increment[i] = info->disp_overlay_increment[i];
}


int Update_Control_Point_aux_old(information *info)
{
	// ローカルパッチ上のコントロールポイントに対応するノットベクトル, 要素を算出
	int CP_loc_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	static int *ele_at_cp = (int *)malloc(sizeof(int) * CP_loc_n);
	static double *knot_vec_at_cp = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	int temp_max = 0;
	for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		int temp = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			temp *= info->No_Control_point[i * info->DIMENSION + j];
		if (temp_max < temp)
			temp_max = temp;
	}
	static double *knot_array = (double *)malloc(sizeof(double) * temp_max * info->DIMENSION);

	static int counter = 0;
	if (counter == 0)
	{
		for (int i = info->Total_Patch_to_mesh[1]; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			// Calc
			Calc_MI3D(info->DIMENSION, &info->Order[i * info->DIMENSION], &info->No_Control_point[i * info->DIMENSION], knot_array);

			// Update
			int cp_to_patch = info->Total_Control_Point_to_patch[i];
			for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
			{
				int connectivity = info->Patch_Control_point[cp_to_patch + j] - info->Total_Control_Point_to_mesh[1];
				ele_at_cp[connectivity] = ele_check(i, &knot_array[j * info->DIMENSION], info);
				tilde_coord(&knot_vec_at_cp[connectivity * info->DIMENSION], &knot_array[j * info->DIMENSION], i, ele_at_cp[connectivity], info);
			}
		}
		counter++;
	}

	// cout << "loc" << endl;
	// for (int i = 0; i < CP_loc_n; i++)
	// {
	// 	cout << i << " ";
	// 	cout << ele_at_cp[i] << " ";
	// 	for (int j = 0; j < info->DIMENSION; j++)
	// 	{
	// 		printf("%.7e ", knot_vec_at_cp[i * info->DIMENSION + j]);
	// 	}
	// 	cout << endl;
	// }

	// newton raphson
	static const int max_itr = info->c.UPDATE_CONTROL_POINT_MAX_ITR;
	vector<double> delta_sol(CP_loc_n * info->DIMENSION, 0.0);
	// double *delta_sol = (double *)calloc(CP_loc_n * info->DIMENSION, sizeof(double));
	static double *right_vec = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	static double *true_f_coord = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	static double *temp_f = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);

	static double *loc_coordinate = (double *)malloc(sizeof(double) * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION);
	size_t M_size = CP_loc_n * info->DIMENSION;
	static double *M_on_point = (double *)malloc(sizeof(double) * M_size * M_size);
	// static double *M_small = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);
	// static double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);

	static double *temp_loc_coordinate = (double *)malloc(sizeof(double) * CP_loc_n * info->DIMENSION);
	double temp_min = 0.0;

	// make true_f
	#pragma omp parallel for
	for (int i = 0; i < CP_loc_n; i++)
	{
		int e = ele_at_cp[i];

		int element_glo = 0;
		double temp_para_glo[MAX_DIMENSION];
		double temp_point_glo[MAX_DIMENSION];
		double temp_coord[MAX_DIMENSION];

		// init
		for (int j = 0; j < info->DIMENSION; j++)
		{
			true_f_coord[i * info->DIMENSION + j] = 0.0;
			temp_coord[j] = 0.0;
		}

		vector<double> R(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R.data(), &knot_vec_at_cp[i * info->DIMENSION], e, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
				double d = info->disp[id] + info->disp_increment[id];
				double coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				temp_coord[k] += R[j] * coord;
				true_f_coord[i * info->DIMENSION + k] += R[j] * (coord + d);
			}

		// overlay
		int itr_n = 0, glo_patch = 0;
		for (int j = 0; j < info->Total_Patch_to_mesh[1]; j++)
		{
			itr_n = calc_patch_parameter_coord(temp_coord, j, temp_para_glo, info);
			glo_patch = j;
			if (itr_n != ERROR)
				break;
		}

		element_glo = ele_check(glo_patch, temp_para_glo, info);
		tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

		vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
				double d_glo = info->disp[id] + info->disp_increment[id];
				true_f_coord[i * info->DIMENSION + k] += R_glo[j] * d_glo;
			}
	}

	// init M_on_point
	#pragma omp parallel for
	for (size_t i = 0; i < M_size * M_size; i++)
	{
		M_on_point[i] = 0.0;
	}

	// newton raphson loop
	int check_convergence = 0;
	int count_NR = 0;
	for (int i = 0; i < max_itr; i++)
	{
		// update loc_coordinate
		if (i == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_loc_n; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					// loc_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[(j + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + k];
					int index = (j + info->Total_Control_Point_to_mesh[1]);
					loc_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[index * (info->DIMENSION + 1) + k] + info->disp_overlay[index * info->DIMENSION + k] + info->disp_overlay_increment[index * info->DIMENSION + k];
				}
		}
		else
		{
			#pragma omp parallel for
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] += delta_sol[j];
		}

		// make right_vec
		#pragma omp parallel for
		for (int j = 0; j < CP_loc_n; j++)
		{
			int e = ele_at_cp[j];

			// init
			for (int k = 0; k < info->DIMENSION; k++)
				temp_f[j * info->DIMENSION + k] = 0.0;

			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &knot_vec_at_cp[j * info->DIMENSION], e, info);
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				for (int l = 0; l < info->DIMENSION; l++)
				{
					int ptr = (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + l;
					temp_f[j * info->DIMENSION + l] += R[k] * loc_coordinate[ptr];

					// make M_on_point
					if (i == 0)
					{
						size_t index = (j * info->DIMENSION + l) * M_size + ptr;
						double val = R[k];

						#pragma omp atomic
						M_on_point[index] += val;
					}
				}
		}

		#pragma omp parallel for
		for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
			right_vec[j] = true_f_coord[j] - temp_f[j];

		// solve
		// GaussianElimination(delta_sol, right_vec, M_on_point, M_size);
		// temp_solver(delta_sol, right_vec, M_on_point, M_size, M_size);
		sparseLU_update_CP(delta_sol.data(), right_vec, M_on_point, M_size, M_size);

		// check
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
		{
			temp1 += pow(true_f_coord[j] - temp_f[j], 2);
			temp2 += pow(true_f_coord[j], 2);
			temp3 += pow(delta_sol[j], 2);
		}
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);
		temp3 = sqrt(temp3);
		double temp_norm = temp1 / temp2;
		printf("\titr_norm: %d %.6e\n", i, temp_norm);
		printf("\tdelta_sol_norm: %d %.6e\n", i, temp3);
		if (temp_norm < MERGE_ERROR)
		{
			check_convergence = 1;
			// update
			for (int i = 0; i < CP_loc_n * info->DIMENSION; i++)
				loc_coordinate[i] += delta_sol[i];
			for (int i = 0; i < CP_loc_n; i++)
				for (int j = 0; j < info->DIMENSION; j++)
					info->disp_overlay_increment[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j] = loc_coordinate[i * info->DIMENSION + j] - info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + j] - info->disp_overlay[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j];
			break;
		}

		if (i == 0)
			temp_min = temp_norm;
		if (temp_norm < temp_min)
		{
			count_NR = 0;
			temp_min = temp_norm;
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				temp_loc_coordinate[j] = loc_coordinate[j] + delta_sol[j];
		}
		else
			count_NR++;

		if (temp_norm > temp_min * 1.0e+2 && count_NR > 10)
		{
			// update
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] = temp_loc_coordinate[j];
			break;
		}
	}

	if (check_convergence == 0)
		return 1;
	else if (check_convergence == 1)
		return 0;
	return 0;
}


int Update_Control_Point_glo(information *info)
{
	// gp_switch(false, info);

	int CP_glo_n = info->Total_Control_Point_to_mesh[1];
	int M_row = info->Total_Element_on_mesh[0] * info->gp[0].n() * info->DIMENSION;
	int M_col = CP_glo_n * info->DIMENSION;

	// newton raphson
	static const int max_itr = info->c.UPDATE_CONTROL_POINT_MAX_ITR;
	double *delta_sol = (double *)calloc(M_col, sizeof(double));
	static double *right_vec = (double *)malloc(sizeof(double) * M_row);
	static double *true_f_coord = (double *)malloc(sizeof(double) * M_row);
	static double *temp_f = (double *)malloc(sizeof(double) * M_row);
	static double *glo_coordinate = (double *)malloc(sizeof(double) * M_col);
	double *temp_glo_coordinate = (double *)malloc(sizeof(double) * M_col);
	static double *M_on_point = (double *)malloc((long)sizeof(double) * M_row * M_col);

	double temp_min = 0.0;

	// make true_f
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Element_on_mesh[0]; i++)
	{
		int e = i + info->Total_Element_to_mesh[0];

		for (int n = 0; n < info->gp[e].n(); n++)
		{
			int element_loc = 0;
			double temp_para_loc[MAX_DIMENSION];
			double temp_point_loc[MAX_DIMENSION];
			double temp_coord[MAX_DIMENSION];

			// init
			for (int j = 0; j < info->DIMENSION; j++)
			{
				true_f_coord[(i * info->gp[e].n() + n) * info->DIMENSION + j] = 0.0;
				temp_coord[j] = 0.0;
			}

			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &info->gp[e].para()[n * info->DIMENSION], e, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
					double d = info->disp[id] + info->disp_increment[id];
					double coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
					temp_coord[k] += R[j] * coord;

					long index = (i * info->gp[e].n() + n) * info->DIMENSION + k;
					double val = R[j] * (coord + d);

					#pragma omp atomic
					true_f_coord[index] += val;
				}

			int status = 0;

			// overlay
			int itr_n = 0, loc_patch = 0;
			for (int j = info->Total_Patch_to_mesh[1]; j < info->Total_Patch_to_mesh[Total_mesh]; j++)
			{
				itr_n = calc_patch_parameter_coord(temp_coord, j, temp_para_loc, info);
				loc_patch = j;
				if (itr_n != ERROR)
				{
					status = 1;
					break;
				}
			}

			if (status)
			{
				element_loc = ele_check(loc_patch, temp_para_loc, info);
				tilde_coord(temp_point_loc, temp_para_loc, loc_patch, element_loc, info);
			}

			if (status)
			{
				vector<double> R_loc(MAX_NO_CP_ON_ELEMENT);
				shape_and_dshape(R_loc.data(), temp_point_loc, element_loc, info);
				for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_loc]]; j++)
					for (int k = 0; k < info->DIMENSION; k++)
					{
						int id = info->Controlpoint_of_Element[element_loc * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
						double d_loc = info->disp[id] + info->disp_increment[id];

						long index = (i * info->gp[e].n() + n) * info->DIMENSION + k;
						double val = R_loc[j] * d_loc;

						#pragma omp atomic
						true_f_coord[index] += val;
					}
			}
		}
	}

	// init M_on_point
	for (long i = 0; i < (long)M_row * M_col; i++)
	{
		M_on_point[i] = 0.0;
	}

	// newton raphson loop
	int check_convergence = 0;
	int count_NR = 0;
	for (int i = 0; i < max_itr; i++)
	{
		// update glo_coordinate
		if (i == 0)
			for (int j = 0; j < CP_glo_n; j++)
				for (int k = 0; k < info->DIMENSION; k++)
					glo_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[j * (info->DIMENSION + 1) + k];
		else
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				glo_coordinate[j] += delta_sol[j];

		// make right_vec
		#pragma omp parallel for
		for (int j = 0; j < info->Total_Element_on_mesh[0]; j++)
		{
			int e = j + info->Total_Element_to_mesh[0];

			for (int n = 0; n < info->gp[e].n(); n++)
			{
				// init
				for (int k = 0; k < info->DIMENSION; k++)
					temp_f[(j * info->gp[e].n() + n) * info->DIMENSION + k] = 0.0;

				vector<double> R(MAX_NO_CP_ON_ELEMENT);
				shape_and_dshape(R.data(), &info->gp[e].para()[n * info->DIMENSION], e, info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int ptr = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;

						long index_1 = (j * info->gp[e].n() + n) * info->DIMENSION + l;
						double val_1 = R[k] * glo_coordinate[ptr];

						#pragma omp atomic
						temp_f[index_1] += val_1;

						// make M_on_point
						if (i == 0)
						{
							long index_2 = index_1 * M_col + ptr;
							double val_2 = R[k];

							#pragma omp atomic
							M_on_point[index_2] += val_2;
						}
					}
			}
		}
		for (int j = 0; j < M_row; j++)
			right_vec[j] = true_f_coord[j] - temp_f[j];

		// solve
		temp_solver(delta_sol, right_vec, M_on_point, M_row, M_col);

		// check
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		for (int j = 0; j < M_row; j++)
		{
			temp1 += pow(true_f_coord[j] - temp_f[j], 2);
			temp2 += pow(true_f_coord[j], 2);
		}
		for (int j = 0; j < M_col; j++)
			temp3 += pow(delta_sol[j], 2);
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);
		temp3 = sqrt(temp3);
		double temp_norm = temp1 / temp2;
		printf("\titr_norm: %d %.6e\n", i, temp_norm);
		printf("\tdelta_sol_norm: %d %.6e\n", i, temp3);
		if (temp_norm < MERGE_ERROR)
		{
			check_convergence = 1;
			// update
			for (int i = 0; i < CP_glo_n * info->DIMENSION; i++)
				glo_coordinate[i] += delta_sol[i];
			for (int i = 0; i < CP_glo_n; i++)
				for (int j = 0; j < info->DIMENSION; j++)
					info->disp_overlay_increment[i * info->DIMENSION + j] = glo_coordinate[i * info->DIMENSION + j] - info->Node_Coordinate[i * (info->DIMENSION + 1) + j] - info->disp_overlay[i * info->DIMENSION + j];
			break;
		}

		if (i == 0)
			temp_min = temp_norm;
		else if (temp_norm < temp_min)
		{
			count_NR = 0;
			temp_min = temp_norm;
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				temp_glo_coordinate[j] = glo_coordinate[j] + delta_sol[j];
		}
		else
			count_NR++;

		if (temp_norm > temp_min * 1.0e+2 && count_NR > 10)
		{
			// update
			for (int j = 0; j < CP_glo_n * info->DIMENSION; j++)
				glo_coordinate[j] = temp_glo_coordinate[j];
			break;
		}
	}

	if (check_convergence == 0)
		return 1;
	else if (check_convergence == 1)
		return 0;
	return 0;
}


int Update_Control_Point_loc(information *info)
{
	// gp_switch(false, info);

	int CP_loc_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	int M_row = info->Total_Element_on_mesh[1] * info->gp[0].n() * info->DIMENSION;
	int M_col = CP_loc_n * info->DIMENSION;

	// newton raphson
	static const int max_itr = info->c.UPDATE_CONTROL_POINT_MAX_ITR;
	double *delta_sol = (double *)calloc(M_col, sizeof(double));
	static double *right_vec = (double *)malloc(sizeof(double) * M_row);
	static double *true_f_coord = (double *)malloc(sizeof(double) * M_row);
	static double *temp_f = (double *)malloc(sizeof(double) * M_row);
	static double *loc_coordinate = (double *)malloc(sizeof(double) * M_col);
	double *temp_loc_coordinate = (double *)malloc(sizeof(double) * M_col);
	static double *M_on_point = (double *)malloc(sizeof(double) * M_row * M_col);

	double temp_min = 0.0;

	// make true_f
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Element_on_mesh[1]; i++)
	{
		int e = i + info->Total_Element_to_mesh[1];

		for (int n = 0; n < info->gp[e].n(); n++)
		{
			int element_glo = 0;
			double temp_para_glo[MAX_DIMENSION];
			double temp_point_glo[MAX_DIMENSION];
			double temp_coord[MAX_DIMENSION];

			// init
			for (int j = 0; j < info->DIMENSION; j++)
			{
				true_f_coord[(i * info->gp[e].n() + n) * info->DIMENSION + j] = 0.0;
				temp_coord[j] = 0.0;
			}

			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R.data(), &info->gp[e].para()[n * info->DIMENSION], e, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
					double d = info->disp[id] + info->disp_increment[id];
					double coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
					temp_coord[k] += R[j] * coord;

					long index = (i * info->gp[e].n() + n) * info->DIMENSION + k;
					double val = R[j] * (coord + d);

					#pragma omp atomic
					true_f_coord[index] += val;
				}

			// overlay
			int itr_n = 0, glo_patch = 0;
			for (int j = 0; j < info->Total_Patch_to_mesh[1]; j++)
			{
				itr_n = calc_patch_parameter_coord(temp_coord, j, temp_para_glo, info);
				glo_patch = j;
				if (itr_n != ERROR)
					break;
			}

			element_glo = ele_check(glo_patch, temp_para_glo, info);
			tilde_coord(temp_point_glo, temp_para_glo, glo_patch, element_glo, info);

			vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
			shape_and_dshape(R_glo.data(), temp_point_glo, element_glo, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int id = info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k;
					double d_glo = info->disp[id] + info->disp_increment[id];

					long index = (i * info->gp[e].n() + n) * info->DIMENSION + k;
					double val = R_glo[j] * d_glo;

					#pragma omp atomic
					true_f_coord[index] += val;
				}
		}
	}

	// init M_on_point
	for (long i = 0; i < (long)M_row * M_col; i++)
	{
		M_on_point[i] = 0.0;
	}

	// newton raphson loop
	int check_convergence = 0;
	int count_NR = 0;
	for (int i = 0; i < max_itr; i++)
	{
		// update loc_coordinate
		if (i == 0)
			for (int j = 0; j < CP_loc_n; j++)
				for (int k = 0; k < info->DIMENSION; k++)
					loc_coordinate[j * info->DIMENSION + k] = info->Node_Coordinate[(j + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + k];
		else
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] += delta_sol[j];

		// make right_vec
		#pragma omp parallel for
		for (int j = 0; j < info->Total_Element_on_mesh[1]; j++)
		{
			int e = j + info->Total_Element_to_mesh[1];

			for (int n = 0; n < info->gp[e].n(); n++)
			{
				// init
				for (int k = 0; k < info->DIMENSION; k++)
					temp_f[(j * info->gp[e].n() + n) * info->DIMENSION + k] = 0.0;

				vector<double> R(MAX_NO_CP_ON_ELEMENT);
				shape_and_dshape(R.data(), &info->gp[e].para()[n * info->DIMENSION], e, info);
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int ptr = (info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] - info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + l;
						
						long index_1 = (j * info->gp[e].n() + n) * info->DIMENSION + l;
						double val_1 = R[k] * loc_coordinate[ptr];

						#pragma omp atomic
						temp_f[index_1] += val_1;

						// make M_on_point
						if (i == 0)
						{
							long index_2 = index_1 * M_col + ptr;
							double val_2 = R[k];

							#pragma omp atomic
							M_on_point[index_2] += val_2;
						}
					}
			}
		}
		for (int j = 0; j < M_row; j++)
			right_vec[j] = true_f_coord[j] - temp_f[j];

		// solve
		temp_solver(delta_sol, right_vec, M_on_point, M_row, M_col);

		// check
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		for (int j = 0; j < M_row; j++)
		{
			temp1 += pow(true_f_coord[j] - temp_f[j], 2);
			temp2 += pow(true_f_coord[j], 2);
		}
		for (int j = 0; j < M_col; j++)
			temp3 += pow(delta_sol[j], 2);
		temp1 = sqrt(temp1);
		temp2 = sqrt(temp2);
		temp3 = sqrt(temp3);
		double temp_norm = temp1 / temp2;
		printf("\titr_norm: %d %.6e\n", i, temp_norm);
		printf("\tdelta_sol_norm: %d %.6e\n", i, temp3);
		if (temp_norm < MERGE_ERROR)
		{
			check_convergence = 1;
			// update
			for (int i = 0; i < CP_loc_n * info->DIMENSION; i++)
				loc_coordinate[i] += delta_sol[i];
			for (int i = 0; i < CP_loc_n; i++)
				for (int j = 0; j < info->DIMENSION; j++)
					info->disp_overlay_increment[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j] = loc_coordinate[i * info->DIMENSION + j] - info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[1]) * (info->DIMENSION + 1) + j] - info->disp_overlay[(i + info->Total_Control_Point_to_mesh[1]) * info->DIMENSION + j];
			break;
		}

		if (i == 0)
			temp_min = temp_norm;
		if (temp_norm < temp_min)
		{
			count_NR = 0;
			temp_min = temp_norm;
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				temp_loc_coordinate[j] = loc_coordinate[j] + delta_sol[j];
		}
		else
			count_NR++;

		if (temp_norm > temp_min * 1.0e+2 && count_NR > 10)
		{
			// update
			for (int j = 0; j < CP_loc_n * info->DIMENSION; j++)
				loc_coordinate[j] = temp_loc_coordinate[j];
			break;
		}
	}

	if (check_convergence == 0)
		return 1;
	else if (check_convergence == 1)
		return 0;
	return 0;
}


inline void Compute_u_at_gauss_point(int e, int gp_id, const vector<double> &R, vector<double> &u, information *info)
{
	int p = info->Element_patch[e];

	// make J and R
	for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
		for (int k = 0; k < info->DIMENSION; k++)
		{
			int cp = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j];
			u[k] += R[j] * (info->disp[cp * info->DIMENSION + k] + info->disp_increment[cp * info->DIMENSION + k]);
		}
	
	// overlay
	bool overlay = info->gp[e].isOverlay()[gp_id];
	if (overlay)
	{
		int overlay_ele = info->gp[e].opp_ele()[gp_id];
		int overlay_patch = info->gp[e].opp_patch()[gp_id];
		static thread_local vector<double> R_overlay(MAX_NO_CP_ON_ELEMENT);

		shape_and_dshape(R_overlay.data(), &info->gp[e].opp_para_tilde()[gp_id * info->DIMENSION], overlay_ele, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[overlay_patch]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int cp = info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + j];
				u[k] += R_overlay[j] * (info->disp[cp * info->DIMENSION + k] + info->disp_increment[cp * info->DIMENSION + k]);
			}
	}
}


inline void Calc_M_G_ele(int e, bool calc_M, vector<vector<double>> &M_ele, vector<vector<double>> &G_ele, information *info)
{
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR_dummy(MAX_NO_CP_ON_ELEMENT * MAX_DIMENSION);
	vector<double> b_dummy(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	int p = info->Element_patch[e];

	// gauss point
	for (int i = 0; i < info->gp[e].n(); i++)
	{
		// make J and R
		double J = Make_B_component(e, &info->gp[e].para()[i * info->DIMENSION], info->gp[e].isOverlay()[i], info->gp[e].opp_ele()[i], &info->gp[e].opp_para_tilde()[i * info->DIMENSION], b_dummy.data(), R, dR_dummy, info);
		vector<double> u(info->DIMENSION, 0.0);
		Compute_u_at_gauss_point(e, i, R, u, info);

		double coef = info->gp[e].w()[i] * J;

		// calc G_ele
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
				G_ele[j][k] += R[j] * u[k] * coef;

		// calc M_ele
		if (calc_M)
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
					M_ele[j][k] += R[j] * R[k] * coef;
	}
}


void Update_Control_Point_glo_using_L2_projection(information *info)
{
	// init
	static vector<double> M;
	static vector<int> M_ptr(info->Total_Control_Point_to_mesh[1] + 1, 0);
	static vector<int> M_col;
	static vector<vector<double>> G(info->Total_Control_Point_to_mesh[1], vector<double>(info->DIMENSION, 0.0));
	static vector<vector<double>> sol(info->Total_Control_Point_to_mesh[1], vector<double>(info->DIMENSION, 0.0));

	static bool isInitialized = false;
	if (!isInitialized)
	{
		vector<set<int>> ptr_debug(info->Total_Control_Point_to_mesh[1]);
		vector<std::mutex> mtx_debug(info->Total_Control_Point_to_mesh[1]);
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
		{
			int p = info->Element_patch[i];
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			{
				int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
				{
					int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k];
					if (col >= row)
					{
						std::lock_guard<std::mutex> lock(mtx_debug[row]);
						ptr_debug[row].insert(col);
					}
				}
			}
		}
		M_ptr[0] = 0;
		for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
			M_ptr[i + 1] = M_ptr[i] + static_cast<int>(ptr_debug[i].size());
		M.resize(M_ptr[info->Total_Control_Point_to_mesh[1]], 0.0);
		M_col.resize(M_ptr[info->Total_Control_Point_to_mesh[1]], 0);
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
		{
			int count = 0;
			for (auto it = ptr_debug[i].begin(); it != ptr_debug[i].end(); ++it)
			{
				M_col[M_ptr[i] + count] = *it;
				count++;
			}
		}
	}

	bool clac_M = (!isInitialized) ? true : false;
	if (!clac_M)
	{
		for (auto &row : G)
			std::fill(row.begin(), row.end(), 0.0);
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[1]; i++)
	{
		int p = info->Element_patch[i];

		vector<vector<double>> M_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->No_Control_point_ON_ELEMENT[p], 0.0));
		vector<vector<double>> G_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->DIMENSION, 0.0));
		Calc_M_G_ele(i, clac_M, M_ele, G_ele, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
		{
			int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];

			// G
			for (int k = 0; k < info->DIMENSION; k++)
			{
				double G_val = G_ele[j][k];
				double *G_ptr_local = &G[row][k];
				#pragma omp atomic
				*G_ptr_local += G_val;
			}

			// M
			if (!clac_M) continue;
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
			{
				int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k];
				if (col >= row)
				{
					int icount = L2_debug_RowCol_to_icount(M_ptr, M_col, row, col);
					double M_val = M_ele[j][k];
					double *M_ptr_local = &M[icount];
					#pragma omp atomic
					*M_ptr_local += M_val;
				}
			}
		}
	}

	// solve
	intel_PARDISO_for_L2_shape_update(M, M_ptr, M_col, G, sol, info->Total_Control_Point_to_mesh[1]);
	cout << "Finish preprocess of L2 projection" << endl << endl;

	// update disp_overlay_increment
	for (int i = 0; i < info->Total_Control_Point_to_mesh[1]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			info->disp_overlay_increment[i * info->DIMENSION + j] = sol[i][j] - info->disp_overlay[i * info->DIMENSION + j];

	if (!isInitialized)
		isInitialized = true;
}


void Update_Control_Point_loc_using_L2_projection(information *info)
{
	// init
	static int local_cp_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	static int local_cp_offset = info->Total_Control_Point_to_mesh[1];
	static vector<double> M;
	static vector<int> M_ptr(local_cp_n + 1, 0);
	static vector<int> M_col;
	static vector<vector<double>> G(local_cp_n, vector<double>(info->DIMENSION, 0.0));
	static vector<vector<double>> sol(local_cp_n, vector<double>(info->DIMENSION, 0.0));

	static bool isInitialized = false;
	if (!isInitialized)
	{
		vector<set<int>> ptr_debug(local_cp_n);
		vector<std::mutex> mtx_debug(local_cp_n);
		#pragma omp parallel for
		for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			{
				int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] - local_cp_offset;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
				{
					int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] - local_cp_offset;
					if (col >= row)
					{
						std::lock_guard<std::mutex> lock(mtx_debug[row]);
						ptr_debug[row].insert(col);
					}
				}
			}
		}
		M_ptr[0] = 0;
		for (int i = 0; i < local_cp_n; i++)
			M_ptr[i + 1] = M_ptr[i] + static_cast<int>(ptr_debug[i].size());
		M.resize(M_ptr[local_cp_n], 0.0);
		M_col.resize(M_ptr[local_cp_n], 0);
		#pragma omp parallel for
		for (int i = 0; i < local_cp_n; i++)
		{
			int count = 0;
			for (auto it = ptr_debug[i].begin(); it != ptr_debug[i].end(); ++it)
			{
				M_col[M_ptr[i] + count] = *it;
				count++;
			}
		}
	}

	bool clac_M = (!isInitialized) ? true : false;
	if (!clac_M)
	{
		for (auto &row : G)
			std::fill(row.begin(), row.end(), 0.0);
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int p = info->Element_patch[i];

		vector<vector<double>> M_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->No_Control_point_ON_ELEMENT[p], 0.0));
		vector<vector<double>> G_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->DIMENSION, 0.0));
		Calc_M_G_ele(i, clac_M, M_ele, G_ele, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
		{
			int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] - local_cp_offset;
			
			// G
			for (int k = 0; k < info->DIMENSION; k++)
			{
				double G_val = G_ele[j][k];
				double *G_ptr_local = &G[row][k];
				#pragma omp atomic
				*G_ptr_local += G_val;
			}
			
			// M
			if (!clac_M) continue;
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
			{
				int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] - local_cp_offset;
				if (col >= row)
				{
					int icount = L2_debug_RowCol_to_icount(M_ptr, M_col, row, col);
					double M_val = M_ele[j][k];
					double *M_ptr_local = &M[icount];
					#pragma omp atomic
					*M_ptr_local += M_val;
				}
			}
		}
	}

	// solve
	intel_PARDISO_for_L2_shape_update(M, M_ptr, M_col, G, sol, local_cp_n);
	cout << "Finish preprocess of L2 projection" << endl << endl;

	// update disp_overlay_increment
	for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		int offset_i = i - local_cp_offset;
		for (int j = 0; j < info->DIMENSION; j++)
			info->disp_overlay_increment[i * info->DIMENSION + j] = sol[offset_i][j] - info->disp_overlay[i * info->DIMENSION + j];
	}

	if (!isInitialized)
		isInitialized = true;
}


void Update_Control_Point_aux_using_L2_projection(information *info)
{
	// init
	static int local_cp_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	static int local_cp_offset = info->Total_Control_Point_to_mesh[1];
	static vector<double> M;
	static vector<int> M_ptr(local_cp_n + 1, 0);
	static vector<int> M_col;
	static vector<vector<double>> G(local_cp_n, vector<double>(info->DIMENSION, 0.0));
	static vector<vector<double>> sol(local_cp_n, vector<double>(info->DIMENSION, 0.0));

	static bool isInitialized = false;
	if (!isInitialized)
	{
		vector<set<int>> ptr_debug(local_cp_n);
		vector<std::mutex> mtx_debug(local_cp_n);
		#pragma omp parallel for
		for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
			{
				int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] - local_cp_offset;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
				{
					int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] - local_cp_offset;
					if (col >= row)
					{
						std::lock_guard<std::mutex> lock(mtx_debug[row]);
						ptr_debug[row].insert(col);
					}
				}
			}
		}
		M_ptr[0] = 0;
		for (int i = 0; i < local_cp_n; i++)
			M_ptr[i + 1] = M_ptr[i] + static_cast<int>(ptr_debug[i].size());
		M.resize(M_ptr[local_cp_n], 0.0);
		M_col.resize(M_ptr[local_cp_n], 0);
		#pragma omp parallel for
		for (int i = 0; i < local_cp_n; i++)
		{
			int count = 0;
			for (auto it = ptr_debug[i].begin(); it != ptr_debug[i].end(); ++it)
			{
				M_col[M_ptr[i] + count] = *it;
				count++;
			}
		}
	}

	bool clac_M = (!isInitialized) ? true : false;
	if (!clac_M)
	{
		for (auto &row : G)
			std::fill(row.begin(), row.end(), 0.0);
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		int p = info->Element_patch[i];

		vector<vector<double>> M_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->No_Control_point_ON_ELEMENT[p], 0.0));
		vector<vector<double>> G_ele(info->No_Control_point_ON_ELEMENT[p], vector<double>(info->DIMENSION, 0.0));
		Calc_M_G_ele(i, clac_M, M_ele, G_ele, info);
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
		{
			int row = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] - local_cp_offset;
			
			// G
			for (int k = 0; k < info->DIMENSION; k++)
			{
				double G_val = G_ele[j][k];
				double *G_ptr_local = &G[row][k];
				#pragma omp atomic
				*G_ptr_local += G_val;
			}
			
			// M
			if (!clac_M) continue;
			for (int k = 0; k < info->No_Control_point_ON_ELEMENT[p]; k++)
			{
				int col = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] - local_cp_offset;
				if (col >= row)
				{
					int icount = L2_debug_RowCol_to_icount(M_ptr, M_col, row, col);
					double M_val = M_ele[j][k];
					double *M_ptr_local = &M[icount];
					#pragma omp atomic
					*M_ptr_local += M_val;
				}
			}
		}
	}

	// solve
	intel_PARDISO_for_L2_shape_update(M, M_ptr, M_col, G, sol, local_cp_n);
	cout << "Finish preprocess of L2 projection" << endl << endl;

	// update disp_overlay_increment
	for (int i = info->Total_Control_Point_to_mesh[1]; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		int offset_i = i - local_cp_offset;
		for (int j = 0; j < info->DIMENSION; j++)
			info->disp_overlay_increment[i * info->DIMENSION + j] = sol[offset_i][j] - info->disp_overlay[i * info->DIMENSION + j];
	}

	if (!isInitialized)
		isInitialized = true;
}


// gauss point for SIGA large deformation
void Update_Gauss_points(information *info)
{
	// temporaly substitute
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			info->Node_Coordinate[i * (info->DIMENSION + 1) + j] += info->disp_overlay[i * info->DIMENSION + j] + info->disp_overlay_increment[i * info->DIMENSION + j];

	// reset eoi
	info->eoi.clear();
	info->eoi.resize(info->Total_Element_to_mesh[Total_mesh]);
	searchOverlappingEle(info);

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		info->gp[i].reset(info);
	gp_switch(true, info);

	// restore
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		for (int j = 0; j < info->DIMENSION; j++)
			info->Node_Coordinate[i * (info->DIMENSION + 1) + j] -= (info->disp_overlay[i * info->DIMENSION + j] + info->disp_overlay_increment[i * info->DIMENSION + j]);

	// auxiliary patch
	if (info->c.BLENDING == 1)
	{
		information *aux = info->aux.get();

		// temporaly substitute
		for (int i = 0; i < aux->Total_Control_Point_to_mesh[Total_mesh]; i++)
			for (int j = 0; j < aux->DIMENSION; j++)
				aux->Node_Coordinate[i * (aux->DIMENSION + 1) + j] += aux->disp_overlay[i * aux->DIMENSION + j] + aux->disp_overlay_increment[i * aux->DIMENSION + j];

		// reset eoi
		aux->eoi.clear();
		aux->eoi.resize(aux->Total_Element_to_mesh[Total_mesh]);
		searchOverlappingEle(aux);

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < aux->Total_Element_to_mesh[Total_mesh]; i++)
			aux->gp[i].reset(aux);
		gp_switch(true, aux);

		// restore
		for (int i = 0; i < aux->Total_Control_Point_to_mesh[Total_mesh]; i++)
			for (int j = 0; j < aux->DIMENSION; j++)
				aux->Node_Coordinate[i * (aux->DIMENSION + 1) + j] -= (aux->disp_overlay[i * aux->DIMENSION + j] + aux->disp_overlay_increment[i * aux->DIMENSION + j]);
	}
}


vector<debug_J_vtk> debug_J;
vector<int> e_index;


bool get_info_J(vector<J_point_info> &jp_list, vector<J_info> &J_list, information *info, vector<double> &circle_r)
{
	#if 0 // 縮退
	circle_r[0] = 10.0;
	circle_r[1] = 10.0;
	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, temp_p + 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 2, 1, 0, 1, 2);

	J_list.emplace_back(temp_p + 0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 1, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 2, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 3, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 4, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 5, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 6, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 7, 1, 0, 1, 10 /*sub_n*/);

	#elif 0 // ローカルのみ
	circle_r[0] = 10.0;
	circle_r[1] = 10.0;
	jp_list.emplace_back(info, 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, 2, 1, 0, 1, 2);

	J_list.emplace_back(0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(1, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(2, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(3, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(4, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(5, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(6, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(7, 1, 0, 1, 10 /*sub_n*/);

	#elif 1 // 特異パッチ
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

	#elif 0 // omar

	#define OMAR
	// #define CALC_BY_INTERPOLATION

	// #define INTERACTION_INTEGRAL_METHOD

	circle_r[0] = 2.0;
	// circle_r[1] = 2.0;
	// circle_r[1] = 1.2;
	circle_r[1] = 0.4;
	jp_list.emplace_back(info, 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, 1, 1, 0, 1, 2);

	J_list.emplace_back(0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(1, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(2, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(3, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(4, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(5, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(6, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(7, 1, 0, 1, 10 /*sub_n*/);

	#elif 0 // hamano SIGA

	#define HAMANO
	// #define INTERACTION_INTEGRAL_METHOD

	circle_r[0] = 10.0;
	circle_r[1] = 10.0;

	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, temp_p +  0, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p +  4, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p +  8, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 12, 1, 0, 1, 2);

	J_list.emplace_back(temp_p +  0, 1, 0, 0, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  1, 1, 0, 0, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  2, 1, 0, 0, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  3, 1, 0, 0, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  4, 1, 0, 1, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  5, 1, 0, 1, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  6, 1, 0, 1, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  7, 1, 0, 1, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  8, 1, 0, 2, 30 /*sub_n*/);
	J_list.emplace_back(temp_p +  9, 1, 0, 2, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 10, 1, 0, 2, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 11, 1, 0, 2, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 12, 1, 0, 3, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 13, 1, 0, 3, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 14, 1, 0, 3, 30 /*sub_n*/);
	J_list.emplace_back(temp_p + 15, 1, 0, 3, 30 /*sub_n*/);

	// jp_list.emplace_back(info, /*patch*/, /*crack_dir*/, /*r_dir*/, /*crack_dir*/, /*theta_dir*/);
	// J_list.emplace_back(/*patch*/, /*crack_dir*/, /*r_dir*/, /*normal_info_num(jp_list_num)*/, /*sub_n*/);

	#elif 0 // APCOM IGA

	#define SQUARE_PATCH
	#define OMAR
	// #define CALC_BY_INTERPOLATION

	circle_r[0] = 2.0;
	circle_r[1] = 2.0;
	// circle_r[1] = 1.6;
	// circle_r[1] = 1.2;
	// circle_r[1] = 0.8;
	// circle_r[1] = 0.4;

	jp_list.emplace_back(info, 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, 2, 1, 0, 1, 2);

	J_list.emplace_back(0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(1, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(2, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(3, 1, 0, 1, 10 /*sub_n*/);

	J_list[0].quadrant = 1; // patch 0 -> first quadrant
	J_list[1].quadrant = 2; // patch 1 -> second quadrant
	J_list[2].quadrant = 1; // patch 2 -> first quadrant
	J_list[3].quadrant = 2; // patch 3 -> second quadrant

	#elif 0 // APCOM SIGA

	#define SQUARE_PATCH
	#define OMAR
	// #define CALC_BY_INTERPOLATION

	circle_r[0] = 2.0;
	// circle_r[1] = 2.0;
	// circle_r[1] = 1.6;
	// circle_r[1] = 1.2;
	// circle_r[1] = 0.8;
	circle_r[1] = 0.4;

	int temp_p = info->Total_Patch_to_mesh[1];
	jp_list.emplace_back(info, temp_p + 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 2, 1, 0, 1, 2);

	J_list.emplace_back(temp_p + 0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 1, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 2, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 3, 1, 0, 1, 10 /*sub_n*/);

	J_list[0].quadrant = 1; // patch 0 -> first quadrant
	J_list[1].quadrant = 2; // patch 1 -> second quadrant
	J_list[2].quadrant = 1; // patch 2 -> first quadrant
	J_list[3].quadrant = 2; // patch 3 -> second quadrant

	#elif 0 // tokutome SIGA

	#define OMAR

	circle_r[0] = 2.0;
	circle_r[1] = 2.0;
	// circle_r[1] = 1.2;
	// circle_r[1] = 0.4;

	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, temp_p + 0, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 1, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 2, 1, 0, 1, 2);
	jp_list.emplace_back(info, temp_p + 3, 1, 0, 1, 2);

	J_list.emplace_back(temp_p +  0, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  4, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  8, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 12, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 44, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 48, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 52, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 56, 1, 0, 0, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  1, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  5, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  9, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 13, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 45, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 49, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 53, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 57, 1, 0, 1, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  2, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  6, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 10, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 14, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 46, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 50, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 54, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 58, 1, 0, 2, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  3, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p +  7, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 11, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 15, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 47, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 51, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 55, 1, 0, 3, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 59, 1, 0, 3, 10 /*sub_n*/);

	// J_list.emplace_back(temp_p +  0, 1, 0, 0, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  4, 1, 0, 0, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  8, 1, 0, 0, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 12, 1, 0, 0, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  1, 1, 0, 1, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  5, 1, 0, 1, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  9, 1, 0, 1, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 13, 1, 0, 1, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  2, 1, 0, 2, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  6, 1, 0, 2, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 10, 1, 0, 2, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 14, 1, 0, 2, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  3, 1, 0, 3, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p +  7, 1, 0, 3, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 11, 1, 0, 3, 10 /*sub_n*/);
	// J_list.emplace_back(temp_p + 15, 1, 0, 3, 10 /*sub_n*/);

	#elif 0 // IGA square integral

	#define SQUARE_PATCH

	jp_list.emplace_back(info, 1 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 2 /*crack_dir*/, 1 /*theta_dir*/);

	J_list.emplace_back(0 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 20 /*sub_n*/);
	J_list.emplace_back(1 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 20 /*sub_n*/);

	J_list[0].quadrant = 2; // patch 0 -> second quadrant
	J_list[1].quadrant = 1; // patch 1 -> first quadrant

	#elif 0 // SIGA square integral

	#define SQUARE_PATCH

	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, temp_p + 1 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 2 /*crack_dir*/, 1 /*theta_dir*/);

	J_list.emplace_back(temp_p + 0 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(temp_p + 1 /*patch*/, 2 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);

	J_list[0].quadrant = 2; // patch 0 -> second quadrant
	J_list[1].quadrant = 1; // patch 1 -> first quadrant

	#elif 0 // IGA square 4tenmage
	// 溶接協会 四点曲げ試験 IGA

	#define SQUARE_PATCH
	#define YONTENMAGE

	jp_list.emplace_back(info, 2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);
	jp_list.emplace_back(info, 3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);

	J_list.emplace_back(2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(4 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(5 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);

	J_list[0].quadrant = 1; // patch 2
	J_list[1].quadrant = 2; // patch 4
	J_list[2].quadrant = 1; // patch 3
	J_list[3].quadrant = 2; // patch 5

	#elif 0 // SIGA square 4tenmage
	// 溶接協会 四点曲げ試験 SIGA

	#define SQUARE_PATCH
	#define YONTENMAGE

	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, temp_p + 3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);
	jp_list.emplace_back(info, temp_p + 4 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);

	J_list.emplace_back(temp_p + 3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 1 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 4 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);
	J_list.emplace_back(temp_p + 2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*normal_info_num(jp_list_num)*/, 10 /*sub_n*/);

	J_list[0].quadrant = 1; // patch 3
	J_list[1].quadrant = 2; // patch 1
	J_list[2].quadrant = 1; // patch 4
	J_list[3].quadrant = 2; // patch 2

	#elif 0 // oishi IGA

	// #define SQUARE_PATCH

	jp_list.emplace_back(info, 0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);

	J_list.emplace_back(0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(1 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);

	#elif 0 // oishi SIGA

	// #define SQUARE_PATCH

	int temp_p = info->Total_Patch_to_mesh[1];

	jp_list.emplace_back(info, 0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 1 /*crack_dir*/, 2 /*theta_dir*/);

	J_list.emplace_back(temp_p + 0 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(temp_p + 1 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(temp_p + 2 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);
	J_list.emplace_back(temp_p + 3 /*patch*/, 1 /*crack_dir*/, 0 /*r_dir*/, 0 /*normal_info_num(jp_list_num)*/, 1 /*sub_n*/);

	#endif

	return true;
}


// calc J integral
void calc_J(information *info)
{
	cout << "start calc_J" << endl;

	#if 0
	return;
	#endif

	if (info->c.ANALYSIS_MODE_1 == 1)
	{
		cout << "skip calc_J in finite deformation analysis" << endl;
		return;
	}

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
			printf("J_list[%d]: %d %d %d\n", i, J_list[i].p, J_list[i].crack_dir, J_list[i].q_dir);
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
		int pow_ng = pow_int(info->c.NG, info->DIMENSION);
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
				make_virtual_crack_extension_area(jp, J_list[i], j, e_on_curve, info);
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
	for (size_t i = 0; i < jp_list.size(); i++)
	{
		// output K
		for (size_t j = 0; j < jp_list[i].J_val.size(); j++)
		{	
			double J = jp_list[i].J_val[j];
			J *= 2.0; // symmetric
			double E_prime = E / (1.0 - nu * nu);
			double K = sqrt(E_prime * J);
			printf("i: %zu j: %zu (θ,K,J): %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, K, J);

			#if defined(OMAR)
				#if !defined(INTERACTION_INTEGRAL_METHOD)
				// double E_prime = E / (1.0 - nu * nu);
				// J *= 2.0; // symmetric
				// double K = sqrt(E_prime * J);
				// printf("i: %zu j: %zu x: %le y: %le (θ,K/K_ref): %.15e %.15e\n", i, j, jp_list[i].coord[j][0], jp_list[i].coord[j][1], atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, K / 112.837916709552);
				// printf("i: %zu j: %zu (θ,K/K_ref,J): %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, K / 159.576912160574, J);
				// printf("i: %zu j: %zu (θ,K,J): %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, K, J);
				#else
				double E_prime_1 = E / (1.0 - nu * nu);
				double factor_1 = E_prime_1 * 0.5; // plane strain
				double factor_2 = E / (2.0 * (1.0 + nu)); // mu
				double K1   = factor_1 * jp_list[i].IIM_K_val[j][0];
				double K2  = factor_1 * jp_list[i].IIM_K_val[j][1];
				double K3 = factor_2 * jp_list[i].IIM_K_val[j][2];

				if (i == 0 && j == 0)
					printf("i j θ J K_I(IIM) K_II(IIM) K_III(IIM)\n");
				printf("%zu %zu %.15e %.15e %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, J, K1, K2, K3);
				#endif
			#elif defined(HAMANO)
				#if defined(INTERACTION_INTEGRAL_METHOD)

				double E_prime = E / (1.0 - nu * nu);
				double factor_1 = E_prime * 0.5; // plane strain
				double factor_2 = E / (2.0 * (1.0 + nu)); // mu
				double K1 = factor_1 * jp_list[i].IIM_K_val[j][0];
				double K2 = factor_1 * jp_list[i].IIM_K_val[j][1];
				double K3 = factor_2 * jp_list[i].IIM_K_val[j][2];

				if (i == 0 && j == 0)
					printf("i j θ J K_I(IIM) K_II(IIM) K_III(IIM)\n");
				printf("%zu %zu %.15e %.15e %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, J, K1, K2, K3);
				#else
				double E_prime = E / (1.0 - nu * nu);
				J *= 2.0; // symmetric
				double K1 = sqrt(E_prime * J);
				if (i == 0 && j == 0)
					printf("i j θ K1 J\n");
				printf("%zu %zu %.15e %.15e %.15e\n", i, j, atan2(jp_list[i].coord[j][1] / circle_r[1], jp_list[i].coord[j][0] / circle_r[0]) * 180.0 / M_PI, K1, J);
				#endif
			#elif defined(YONTENMAGE)
				J *= 2.0; // symmetric
				static bool isInitialized = false;
				static double total_arc_length = 0.0;
				if (!isInitialized)
				{
					for (size_t k = 0; k < jp_list.size(); k++)
						for (size_t l = 0; l < jp_list[k].delta_arc_length.size(); l++)
							total_arc_length += jp_list[k].delta_arc_length[l];
					isInitialized = true;
				}
				static double sum_arc_length = 0.0;
				if (i == 0 && j == 0)
					sum_arc_length = 0.0;
				double current_l = jp_list[i].delta_arc_length_at_half_para[j] + sum_arc_length;
				double normalized_current_l = current_l / total_arc_length;
				if (i == 0 && j == 0)
					printf("i j arc_length normalized_arc_length J\n");
				printf("%zu %zu %.15e %.15e %.15e\n", i, j, current_l, normalized_current_l, J);
				sum_arc_length += jp_list[i].delta_arc_length[j];
			#endif
		}
	}

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


void make_virtual_crack_extension_area(J_point_info &jp, J_info &J, size_t line, int e_on_curve, information *info)
{
	vector<double> &J_integral_area = jp.J_integral_area;
	vector<vector<double>> &normal = jp.normal;
	vector<vector<double>> &coord = jp.coord;

	double delta = 1.0 / static_cast<double>(J.sub_n);

	for (int i = 0; i < J.sub_n; i++)
	{
		// calc para
		vector<double> para_c(info->DIMENSION, 0.0);
		vector<double> para_s(info->DIMENSION, 0.0);
		vector<double> para_e(info->DIMENSION, 0.0);
		vector<double> para_n(info->DIMENSION, 0.0);
		for (int m = 0; m < info->DIMENSION; m++)
		{
			// para
			if (m == J.crack_dir)
			{
				double normalized_para_c = (static_cast<double>(i) + 0.5) * delta;
				double normalized_para_s = (static_cast<double>(i) + 0.0) * delta;
				double normalized_para_e = (static_cast<double>(i) + 1.0) * delta;
				double normalized_para_n = (static_cast<double>(i) + 0.0) * delta;
				para_c[m] = normalized_para_c * 2.0 - 1.0;
				para_s[m] = normalized_para_s * 2.0 - 1.0;
				para_e[m] = normalized_para_e * 2.0 - 1.0;
				para_n[m] = normalized_para_n * 2.0 - 1.0;
			}
			else
			{
				para_c[m] = -1.0;
				para_s[m] = -1.0;
				para_e[m] = -1.0;
				para_n[m] = -1.0;
			}
		}
		para_n[J.r_dir] = 1.0;

		vector<double> coord_c(info->DIMENSION, 0.0);
		vector<double> coord_s(info->DIMENSION, 0.0);
		vector<double> coord_e(info->DIMENSION, 0.0);
		vector<double> coord_n(info->DIMENSION, 0.0);
		vector<double> R_c(MAX_NO_CP_ON_ELEMENT);
		vector<double> R_s(MAX_NO_CP_ON_ELEMENT);
		vector<double> R_e(MAX_NO_CP_ON_ELEMENT);
		vector<double> R_n(MAX_NO_CP_ON_ELEMENT);
		shape_and_dshape(R_c.data(), para_c.data(), e_on_curve, info);
		shape_and_dshape(R_s.data(), para_s.data(), e_on_curve, info);
		shape_and_dshape(R_e.data(), para_e.data(), e_on_curve, info);
		shape_and_dshape(R_n.data(), para_n.data(), e_on_curve, info);
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
	D_elastic_smart(D_e.data(), info);

	// init J_val
	for (size_t i = 0; i < jp_list.size(); i++)
		for (size_t j = 0; j < jp_list[i].J_val.size(); j++)
			jp_list[i].J_val[j] = 0.0;
	
	// init IIM_K_val
	#if defined(INTERACTION_INTEGRAL_METHOD)
	for (size_t i = 0; i < jp_list.size(); i++)
		for (size_t j = 0; j < jp_list[i].IIM_K_val.size(); j++)
			for (size_t k = 0; k < jp_list[i].IIM_K_val[j].size(); k++)
				jp_list[i].IIM_K_val[j][k] = 0.0;
	#endif

	// #pragma omp parallel
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

		// #pragma omp for schedule(dynamic) nowait
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

				// bool debug_point = (
				// 	current_sub_ele.e == J_list[0].e_list[0][0] &&
				// 	current_sub_ele.line == 0 &&
				// 	current_sub_ele.slice == 0 &&
				// 	current_sub_ele.face == 0 &&
				// 	j == 0
				// );

				// calc para
				vector<double> disp_grad(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> delta_total_strain(D_MATRIX_SIZE, 0.0);
				vector<double> strain_trial(D_MATRIX_SIZE, 0.0);
				vector<double> bl(D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
				vector<double> bnl(info->DIMENSION * info->DIMENSION * MAX_KIEL_SIZE, 0.0);
				vector<double> u(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
				vector<double> du(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
				double *para_ptr = current_sub_ele.para_tilde.data() + j * info->DIMENSION;

				Make_B_Linear_small_strain(e, para_ptr, bl.data(), info);
				Make_B_Nonlinear_small_strain(e, para_ptr, bnl.data(), info);

				// if (debug_point)
				// {
				// 	for (int k = 0; k < info->DIMENSION; k++)
				// 		for (int l = 0; l < info->DIMENSION; l++)
				// 		{
				// 			int offset = k * info->DIMENSION + l;
				// 			for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]] * info->DIMENSION; m++)
				// 				cout << "grad" << offset << " " << bnl[offset * MAX_KIEL_SIZE + m] << endl;
				// 		}
				// 	exit(0);
				// }

				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int id = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
						u[k * info->DIMENSION + l] = info->disp[id];
						du[k * info->DIMENSION + l] = info->disp[id] - info->disp_previous[id];
					}
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
						for (int m = 0; m < info->DIMENSION; m++)
						{
							delta_total_strain[k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * du[l * info->DIMENSION + m];
							strain_trial[k] += bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * u[l * info->DIMENSION + m];
						}
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
					{
						int offset = k * info->DIMENSION + l;
						for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]] * info->DIMENSION; m++)
							disp_grad[k * info->DIMENSION + l] += u[m] * bnl[offset * MAX_KIEL_SIZE + m];
					}

				// overlay
				if (Total_mesh > 1 && current_sub_ele.overlay[j])
				{
					// calc opp_para
					double *opp_para_ptr = current_sub_ele.opp_para_tilde.data() + j * info->DIMENSION;
					int overlay_ele = current_sub_ele.overlay_ele[j];
					vector<double> opp_bl(D_MATRIX_SIZE * info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);
					vector<double> opp_bnl(info->DIMENSION * info->DIMENSION * MAX_KIEL_SIZE, 0.0);
					vector<double> opp_u(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
					vector<double> opp_du(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
					Make_B_Linear_small_strain(overlay_ele, opp_para_ptr, opp_bl.data(), info);
					Make_B_Nonlinear_small_strain(overlay_ele, opp_para_ptr, opp_bnl.data(), info);
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int id = info->Controlpoint_of_Element[overlay_ele * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l;
							opp_u[k * info->DIMENSION + l] = info->disp[id];
							opp_du[k * info->DIMENSION + l] = info->disp[id] - info->disp_previous[id];
						}
		
					for (int k = 0; k < D_MATRIX_SIZE; k++)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
							{
								delta_total_strain[k] += opp_bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * opp_du[l * info->DIMENSION + m];
								strain_trial[k] += opp_bl[k * MAX_KIEL_SIZE + l * info->DIMENSION + m] * opp_u[l * info->DIMENSION + m];
							}
					for (int k = 0; k < info->DIMENSION; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int offset = k * info->DIMENSION + l;
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[overlay_ele]] * info->DIMENSION; m++)
								disp_grad[k * info->DIMENSION + l] += opp_u[m] * opp_bnl[offset * MAX_KIEL_SIZE + m];
						}
				}

				// if (debug_point)
				// {
				// 	cout << "disp_grad: ";
				// 	for (int k = 0; k < info->DIMENSION * info->DIMENSION; k++)
				// 		cout << disp_grad[k] << " ";
				// 	cout << endl;
				// }

				// Calculate trial stresses: {sigma}^trial = [D] * {epsilon}^trial (small deformation)
				vector<double> stress_trial(D_MATRIX_SIZE, 0.0);
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					for (int l = 0; l < D_MATRIX_SIZE; l++)
						stress_trial[k] += D_e[k * D_MATRIX_SIZE + l] * strain_trial[l];
				double equivalent_stress_trial = calc_equivalent_stress(stress_trial.data());

				// update previous
				current_sub_ele.previous_elastic_strain[j] = current_sub_ele.elastic_strain[j];
				current_sub_ele.previous_stress[j] = current_sub_ele.stress[j];

				// update current trial
				current_sub_ele.elastic_strain[j] = strain_trial;
				current_sub_ele.stress[j] = stress_trial;

				// elastoplastic (if elastic, do nothing)
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

				// W
				for (int k = 0; k < D_MATRIX_SIZE; k++)
					current_sub_ele.previous_W[j] += 0.5 * (current_sub_ele.stress[j][k] + current_sub_ele.previous_stress[j][k]) * delta_total_strain[k];
				double W = current_sub_ele.previous_W[j];

				// (∂q/∂x_i) = ∂q/∂ξ * ∂ξ/∂x_i
				vector<double> q_grad(info->DIMENSION, 0.0);
				vector<double> q_xitildei(info->DIMENSION, 0.0);
				vector<double> a_inv(info->DIMENSION * info->DIMENSION, 0.0);
	
				vector<double> xi_xitilde(info->DIMENSION * info->DIMENSION, 0.0);
				for (int k = 0; k < info->DIMENSION; k++)
					xi_xitilde[k] = dShapeFunc_from_paren(k, e, info);

				#ifndef SQUARE_PATCH
				calc_q_gradient(e, para_ptr, current_sub_ele.slice, current_sub_ele.face, current_J_list, q_xitildei, info);
				#else
				calc_q_gradient_square_patch(e, para_ptr, current_sub_ele.slice, current_sub_ele.face, current_J_list, q_xitildei, info);
				#endif

				double jac = Make_a_q(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, delta, true);
				Make_a_q(e, para_ptr, a_inv.data(), info, current_J_list.crack_dir, 0.0, false);

				vector<double> coef_vec(info->DIMENSION, 1.0);
				coef_vec[current_J_list.crack_dir] = delta;
	
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
						// 変更前
						// J_int[k] += (stress_disp_grad[k * info->DIMENSION + l] - W_I[k * info->DIMENSION + l]) * q_grad[l] * coef;

						// 変更後（正しい？）
						J_int[k] += (stress_disp_grad[l * info->DIMENSION + k] - W_I[k * info->DIMENSION + l]) * q_grad[l] * coef;
					
						// if (debug_point)
						// {
						// 	cout << "J_int[" << k << "] = " << J_int[k] << ", stress_disp_grad[" << l << "][" << k << "] = " << stress_disp_grad[l * info->DIMENSION + k] << ", W_I[" << k << "][" << l << "] = " << W_I[k * info->DIMENSION + l] << ", q_grad[" << l << "] = " << q_grad[l] << ", coef = " << coef << endl;
						// 	fflush(stdout);
						// }
					
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
							// IIM_int[l] += (stress_aux_disp_grad[l * info->DIMENSION + m] + aux_stress_disp_grad[l * info->DIMENSION + m] - aux_W_I[l * info->DIMENSION + m]) * q_grad[m] * coef;
						}

					// add to IIM_K_val
					for (int l = 0; l < info->DIMENSION; l++)
					{
						double val = IIM_int[l] * normal[index_J][l];
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
	const int gp_1d = info->c.NG;
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


// interaction integral method
void calc_polar_coord_for_auxiliary_fields(double &r, double &theta, const vector<vector<double>> &Q, const double *crack_tip_coord, const double *current_coord)
{
	// vector from crack tip to current coord
	vector<double> rc_vector(3, 0.0);
	for (int i = 0; i < 3; i++)
		rc_vector[i] = current_coord[i] - crack_tip_coord[i];

	// transform to crack tip local coord
	vector<double> rc_local(3, 0.0);
	for (int i = 0; i < 3; i++)
	{
		// Dot product with i-th column of Q (basis vector i)
		// rc_local[0] = rc_vector · binormal (Q[*][0])
		// rc_local[1] = rc_vector · normal (Q[*][1])
		// rc_local[2] = rc_vector · tangent (Q[*][2])
		rc_local[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			rc_local[i] += Q[j][i] * rc_vector[j];
			// rc_local[i] += Q[i][j] * rc_vector[j];
		}
	}

	// calc r (only in crack plane: crack propagation direction and normal direction)
	// rc_local[0]: crack propagation direction (local x)
	// rc_local[1]: normal direction (local y)
	// rc_local[2]: tangent direction (local z)
	r = sqrt(rc_local[0] * rc_local[0] + rc_local[1] * rc_local[1]);

	// calc theta (angle from crack propagation direction to normal direction)
	theta = atan2(rc_local[1], rc_local[0]); // atan2(local y, local x)
}


void calc_auxiliary_fields(vector<double> &aux_disp_grad, vector<double> &aux_strain_tensor, vector<double> &aux_stress_tensor, const double *r, const double *theta, const vector<vector<double>> &Q, int mode, vector<double> &D_e)
{
	const double rr = *r;
	const double th = *theta;

	const double sqrt_r = sqrt(rr / (2.0 * M_PI));
	const double inv_sqrt_r = 1.0 / sqrt(2.0 * M_PI * rr);

	// plane strain
	const double kappa = 3.0 - 4.0 * nu;
	const double inv_2mu = (1.0 + nu) / E;

	const double c1 = cos(th / 2.0);
	const double s1 = sin(th / 2.0);
	const double c3 = cos(3.0 * th / 2.0);
	const double s3 = sin(3.0 * th / 2.0);

	// clear
	std::fill(aux_disp_grad.begin(), aux_disp_grad.end(), 0.0);
	std::fill(aux_strain_tensor.begin(), aux_strain_tensor.end(), 0.0);
	std::fill(aux_stress_tensor.begin(), aux_stress_tensor.end(), 0.0);

	// =========================
	// auxiliary displacement gradient
	// =========================
	// Compute local derivatives in crack plane (x': binormal, y': normal)
	double cos_th = cos(th);
	double sin_th = sin(th);
	vector<double> aux_disp_grad_local(9, 0.0);
	for (int j = 0; j < 3; ++j)
	{
		double drdx = (j == 0) ? cos_th : (j == 1) ? sin_th : 0.0;
		double dtdx = (j == 0) ? (-sin_th / rr) : (j == 1) ? (cos_th / rr) : 0.0;

		if (mode == 0)
		{
			// ---------- Mode I ----------
			aux_disp_grad_local[0 * 3 + j] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa - 1.0) * c1 - c3) * drdx + sqrt_r * (-(kappa - 1.0) * 0.5 * s1 + 1.5 * s3) * dtdx);
			aux_disp_grad_local[1 * 3 + j] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa + 1.0) * s1 - s3) * drdx + sqrt_r * ((kappa + 1.0) * 0.5 * c1 - 1.5 * c3) * dtdx);
			// aux_disp_grad_local[j * 3 + 0] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa - 1.0) * c1 - c3) * drdx + sqrt_r * (-(kappa - 1.0) * 0.5 * s1 + 1.5 * s3) * dtdx);
			// aux_disp_grad_local[j * 3 + 1] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa + 1.0) * s1 - s3) * drdx + sqrt_r * ((kappa + 1.0) * 0.5 * c1 - 1.5 * c3) * dtdx);
			// du3/dx'j = 0.0
		}
		else if (mode == 1)
		{
			// ---------- Mode II ----------
			aux_disp_grad_local[0 * 3 + j] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa + 1.0) * s1 + s3) * drdx + sqrt_r * ((kappa + 1.0) * 0.5 * c1 + 1.5 * c3) * dtdx);
			aux_disp_grad_local[1 * 3 + j] = inv_2mu * (0.5 * inv_sqrt_r * (-(kappa - 1.0) * c1 - c3) * drdx + sqrt_r * ((kappa - 1.0) * 0.5 * s1 + 1.5 * s3) * dtdx);
			// aux_disp_grad_local[j * 3 + 0] = inv_2mu * (0.5 * inv_sqrt_r * ((kappa + 1.0) * s1 + s3) * drdx + sqrt_r * ((kappa + 1.0) * 0.5 * c1 + 1.5 * c3) * dtdx);
			// aux_disp_grad_local[j * 3 + 1] = inv_2mu * (0.5 * inv_sqrt_r * (-(kappa - 1.0) * c1 - c3) * drdx + sqrt_r * ((kappa - 1.0) * 0.5 * s1 + 1.5 * s3) * dtdx);
			// du3/dx'j = 0.0
		}
		else if (mode == 2)
		{
			// ---------- Mode III ----------
			aux_disp_grad_local[2 * 3 + j] = inv_2mu * (0.5 * inv_sqrt_r * s1 * drdx + sqrt_r * 0.5 * c1 * dtdx);
			// aux_disp_grad_local[j * 3 + 2] = inv_2mu * (0.5 * inv_sqrt_r * s1 * drdx + sqrt_r * 0.5 * c1 * dtdx);
		}
	}

	// Rotate local tensor to global: A_global = Q * A_local * Q^T (Q columns are basis vectors)
	aux_disp_grad.assign(aux_disp_grad_local.begin(), aux_disp_grad_local.end());
	vector<vector<double>> Q_tf(3, vector<double>(3, 0.0));
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
		{
			// Q_tf[i][j] = Q[i][j];
			Q_tf[i][j] = Q[j][i];
		}
	tensor_transformation(aux_disp_grad.data(), Q_tf);

	// =========================
	// auxiliary strain tensor
	// =========================
	aux_strain_tensor[0] = aux_disp_grad[0 * 3 + 0]; // epsilon_11
	aux_strain_tensor[1] = 0.5 * (aux_disp_grad[0 * 3 + 1] + aux_disp_grad[1 * 3 + 0]); // epsilon_12
	aux_strain_tensor[2] = 0.5 * (aux_disp_grad[0 * 3 + 2] + aux_disp_grad[2 * 3 + 0]); // epsilon_13
	aux_strain_tensor[3] = 0.5 * (aux_disp_grad[1 * 3 + 0] + aux_disp_grad[0 * 3 + 1]); // epsilon_21
	aux_strain_tensor[4] = aux_disp_grad[1 * 3 + 1]; // epsilon_22
	aux_strain_tensor[5] = 0.5 * (aux_disp_grad[1 * 3 + 2] + aux_disp_grad[2 * 3 + 1]); // epsilon_23
	aux_strain_tensor[6] = 0.5 * (aux_disp_grad[2 * 3 + 0] + aux_disp_grad[0 * 3 + 2]); // epsilon_31
	aux_strain_tensor[7] = 0.5 * (aux_disp_grad[2 * 3 + 1] + aux_disp_grad[1 * 3 + 2]); // epsilon_32
	aux_strain_tensor[8] = aux_disp_grad[2 * 3 + 2]; // epsilon_33

	// temporary voigt notation for strain and stress tensor
	static thread_local vector<double> aux_strain_voigt(D_MATRIX_SIZE, 0.0);
	static thread_local vector<double> aux_stress_voigt(D_MATRIX_SIZE, 0.0);
	std::fill(aux_strain_voigt.begin(), aux_strain_voigt.end(), 0.0);
	std::fill(aux_stress_voigt.begin(), aux_stress_voigt.end(), 0.0);
	aux_strain_voigt[0] = aux_strain_tensor[0]; // epsilon_11
	aux_strain_voigt[1] = aux_strain_tensor[4]; // epsilon_22
	aux_strain_voigt[2] = aux_strain_tensor[8]; // epsilon_33
	aux_strain_voigt[3] = aux_strain_tensor[1] + aux_strain_tensor[3]; // epsilon_12
	aux_strain_voigt[4] = aux_strain_tensor[5] + aux_strain_tensor[7]; // epsilon_23
	aux_strain_voigt[5] = aux_strain_tensor[6] + aux_strain_tensor[2]; // epsilon_13

	// =========================
	// auxiliary stress tensor
	// =========================
	for (int i = 0; i < D_MATRIX_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			aux_stress_voigt[i] += D_e[i * D_MATRIX_SIZE + j] * aux_strain_voigt[j];

	// convert voigt notation to tensor
	aux_stress_tensor[0] = aux_stress_voigt[0]; // sigma_11
	aux_stress_tensor[1] = aux_stress_voigt[3]; // sigma_12
	aux_stress_tensor[2] = aux_stress_voigt[5]; // sigma_13
	aux_stress_tensor[3] = aux_stress_voigt[3]; // sigma_21
	aux_stress_tensor[4] = aux_stress_voigt[1]; // sigma_22
	aux_stress_tensor[5] = aux_stress_voigt[4]; // sigma_23
	aux_stress_tensor[6] = aux_stress_voigt[5]; // sigma_31
	aux_stress_tensor[7] = aux_stress_voigt[4]; // sigma_32
	aux_stress_tensor[8] = aux_stress_voigt[2]; // sigma_33
}

#define NOTRANSPOSE 0 // 1: [T] * [A] * [T]^T, 0: [T]^T * [A] * [T]
void tensor_transformation(double *target_tensor, vector<vector<double>> &tf_tensor)
{
	vector<double> temp_tensor(9, 0.0);
	
	#if NOTRANSPOSE
	// [A'] = [T] * [A] * [T]^T
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp_tensor[i * 3 + j] = 0.0;
			for (int k = 0; k < 3; k++)
			{
				temp_tensor[i * 3 + j] += tf_tensor[i][k] * target_tensor[k * 3 + j];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			target_tensor[i * 3 + j] = 0.0;
			for (int k = 0; k < 3; k++)
			{
				target_tensor[i * 3 + j] += temp_tensor[i * 3 + k] * tf_tensor[j][k];
			}
		}
	}
	#else
	// [A'] = [T]^T * [A] * [T]
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp_tensor[i * 3 + j] = 0.0;
			for (int k = 0; k < 3; k++)
			{
				temp_tensor[i * 3 + j] += tf_tensor[k][i] * target_tensor[k * 3 + j];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			target_tensor[i * 3 + j] = 0.0;
			for (int k = 0; k < 3; k++)
			{
				target_tensor[i * 3 + j] += temp_tensor[i * 3 + k] * tf_tensor[k][j];
			}
		}
	}
	#endif
}


void vector_transformation(double *target_vector, vector<vector<double>> &tf_vector)
{
	vector<double> temp_vector(3, 0.0);
	
	#if NOTRANSPOSE
	// [v'] = [T] * [v]
	for (int i = 0; i < 3; i++)
	{
		temp_vector[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			temp_vector[i] += tf_vector[i][j] * target_vector[j];
		}
	}

	for (int i = 0; i < 3; i++)
		target_vector[i] = temp_vector[i];
	#else
	// [v'] = [T]^T * [v]
	for (int i = 0; i < 3; i++)
	{
		temp_vector[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			temp_vector[i] += tf_vector[j][i] * target_vector[j];
		}
	}
	for (int i = 0; i < 3; i++)
		target_vector[i] = temp_vector[i];
	#endif
}


// Nonlinear tool
void identify_tensor(int rank, double *A)
{
	for (int i = 0; i < rank; i++)
		for (int j = 0; j < rank; j++)
		{
			if (i == j)
				A[i * rank + j] = 1.0;
			else
				A[i * rank + j] = 0.0;
		}
}


double tensor_norm(int rank, double *A)
{
	double norm = 0.0;
	for (int i = 0; i < rank; i++)
		for (int j = 0; j < rank; j++)
			norm += A[i * rank + j] * A[i * rank + j];
	return sqrt(norm);
}


void tensor_exponential(int rank, double *A_in, double *A_out)
{
	constexpr int taylor_term_max = 1000;
	constexpr double taylor_term_tolerance = MATH_TOLERANCE;
	thread_local double *power_tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);
	thread_local double *previous_power_tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);

	// identify tensor
	identify_tensor(rank, power_tensor);
	identify_tensor(rank, A_out);

	double nn = 1.0, factorial = 1.0;
	for (int n = 1; n < taylor_term_max; n++)
	{
		// substitution
		for (int i = 0; i < rank * rank; i++)
			previous_power_tensor[i] = power_tensor[i];

		// make power_tensor
		for (int i = 0; i < rank; i++)
			for (int j = 0; j < rank; j++)
			{
				double a = 0.0;
				for (int k = 0; k < rank; k++)
					a += previous_power_tensor[i * rank + k] * A_in[k * rank + j];
				power_tensor[i * rank + j] = a;
			}

		for (int i = 0; i < rank; i++)
			for (int j = 0; j < rank; j++)
				A_out[i * rank + j] += power_tensor[i * rank + j] / factorial;
		
		// norm
		if (tensor_norm(rank, power_tensor) / factorial <= taylor_term_tolerance)
			break;

		nn += 1.0;
		factorial *= nn;
	}
}


void tensor_logarithm(int rank, double *A_in, double *A_out)
{
	constexpr int taylor_term_max = 1000;
	constexpr double taylor_term_tolerance = MATH_TOLERANCE;
	thread_local double *tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);
	thread_local double *square_tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);
	thread_local double *power_tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);
	thread_local double *previous_power_tensor = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION);

	// make tensor
	for (int i = 0; i < rank * rank; i++)
		tensor[i] = A_in[i];
	for (int i = 0; i < rank; i++)
		tensor[i * rank + i] += 1.0;
	if (rank == 2)
		InverseMatrix_2x2(tensor);
	else if (rank == 3)
		InverseMatrix_3x3(tensor);
	for (int i = 0; i < rank * rank; i++)
		tensor[i] *= -2.0;
	for (int i = 0; i < rank; i++)
		tensor[i * rank + i] += 1.0;

	// substitution
	for (int i = 0; i < rank * rank; i++)
		power_tensor[i] = tensor[i];
	
	// make square_tensor
	for (int i = 0; i < rank; i++)
		for (int j = 0; j < rank; j++)
		{
			double a = 0.0;
			for (int k = 0; k < rank; k++)
				a += tensor[i * rank + k] * tensor[k * rank + j];
			square_tensor[i * rank + j] = a;
		}

	// init
	for (int i = 0; i < rank * rank; i++)
		A_out[i] = 0.0;

	int n;
	double nn = 0.0;
	for (n = 0; n < taylor_term_max; n++)
	{
		for (int i = 0; i < rank; i++)
			for (int j = 0; j < rank; j++)
				A_out[i * rank + j] += power_tensor[i * rank + j] / (2.0 * nn + 1.0);
		
		// norm
		if (tensor_norm(rank, power_tensor) / (2.0 * nn + 1.0) <= taylor_term_tolerance)
			break;

		// substitution
		for (int i = 0; i < rank * rank; i++)
			previous_power_tensor[i] = power_tensor[i];

		// make power_tensor
		for (int i = 0; i < rank; i++)
			for (int j = 0; j < rank; j++)
			{
				double a = 0.0;
				for (int k = 0; k < rank; k++)
					a += previous_power_tensor[i * rank + k] * square_tensor[k * rank + j];
				power_tensor[i * rank + j] = a;
			}

		nn++;
	}

	if (n == taylor_term_max)
		for (int i = 0; i < rank * rank; i++)
			A_out[i] = nan("");
	else
		for (int i = 0; i < rank * rank; i++)
			A_out[i] *= 2.0;
}


double calc_inverse(const double value)
{
	return 1.0 / value;
}


void tensor_logarithm_derivative(int rank, double *A_in, double *A_out)
{
	isotropic_tensor_function_derivative(rank, A_in, A_out, log, calc_inverse);
}


void calc_3x3_square(double *A_in, double *A_out)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			double temp = 0.0;
			for (int k = 0; k < 3; k++)
				temp += A_in[i * 3 + k] * A_in[k * 3 + j];
			A_out[i * 3 + j] = temp;
		}
}


double calc_3x3_trace(double *A_in)
{
	return A_in[0] + A_in[4] + A_in[8];
}


double calc_3x3_determinant(double *A_in)
{
	return A_in[0] * A_in[4] * A_in[8] + A_in[1] * A_in[5] * A_in[6] + A_in[2] * A_in[3] * A_in[7] - A_in[0] * A_in[5] * A_in[7] - A_in[1] * A_in[3] * A_in[8] - A_in[2] * A_in[4] * A_in[6];
}


void swap_reals(double *value1, double *value2)
{
	double temp = *value1;
	*value1  = *value2;
	*value2  = temp;
}


void sort_vector(double *vector, const int num)
{
	for (int i = 0; i < num - 1; i++)
		for (int j = num - 1; j > i; j--)
			if (vector[j - 1] > vector[j])
				swap_reals(&vector[j - 1], &vector[j]);
}


void reverse_vector(double *vector, const int num)
{
	for (int i = 0; i < num / 2; i++)
		swap_reals(&vector[i], &vector[num - 1 - i]);
}


double calc_square(const double value)
{
	return value * value;
}


double calc_cube(const double value)
{
	return value * value * value;
}


void calc_3x3_eigenvalues(double *eigenvalues, double *A_in)
{
	const double discriminant_tolerance = MATH_TOLERANCE;
	double square_tensor[9];

	// Calculate square of tensor
	calc_3x3_square(A_in, square_tensor);

	// Calculate invariants
	double i1 = calc_3x3_trace(A_in);
	double i2 = 0.5 * (i1 * i1 - calc_3x3_trace(square_tensor));
	double i3 = calc_3x3_determinant(A_in);

	// Calculate Q and R
	double q  = (i1 * i1 - 3.0 * i2) / 9.0;
	double sqrt_q3 = sqrt(q * q * q);
	double r  = (-2.0 * i1 * i1 * i1 + 9.0 * i1 * i2 - 27.0 * i3) / 54.0;

	// Calculate eigenvalues
	if (q < 0.0 || fabs(q) <= discriminant_tolerance * (i1 * i1 / 9.0))
	{
		const double b = i1 / 3.0;

		eigenvalues[0] = b;
		eigenvalues[1] = b;
		eigenvalues[2] = b;
	}
	else if (r - sqrt_q3 > 0 || fabs(r - sqrt_q3) <= discriminant_tolerance * fabs(sqrt_q3))
	{
		const double a = sqrt(q);
		const double b = i1 / 3.0;

		eigenvalues[0] = -2.0 * a + b;
		eigenvalues[1] = a + b;
		eigenvalues[2] = a + b;
	}
	else if (r + sqrt_q3 < 0 || fabs(r + sqrt_q3) <= discriminant_tolerance * fabs(sqrt_q3))
	{
		const double a = sqrt(q);
		const double b = i1 / 3.0;

		eigenvalues[0] = -a + b;
		eigenvalues[1] = 2.0 * a + b;
		eigenvalues[2] = -a + b;
	}
	else
	{
		const double theta = acos(r / sqrt_q3) / 3.0;
		const double cos_theta = cos(theta);
		const double sin_theta = sin(theta);
		const double a = sqrt(q);
		const double b = i1 / 3.0;

		eigenvalues[0] = -2.0 * a * cos_theta + b;
		eigenvalues[1] = a * (cos_theta + sqrt(3.0) * sin_theta) + b;
		eigenvalues[2] = a * (cos_theta - sqrt(3.0) * sin_theta) + b;
	}

	// Sort and reverse eigenvalues
	sort_vector(eigenvalues, 3);
	reverse_vector(eigenvalues, 3);
}


void calc_3x3_eigenprojections(double *eigenprojections, double *eigenvalues, double *A_in)
{
	double square_tensor[9];
	double identity_tensor[9];
	identify_tensor(3, identity_tensor);
 
	// Calculate square of tensor
	calc_3x3_square(A_in, square_tensor);

	// Calculate eigenprojections
	for (int i = 0; i < 3; i++)
	{
		const double i1 = calc_3x3_trace(A_in);
		const double i3 = calc_3x3_determinant(A_in);
		const double coefficient = eigenvalues[i] / (2.0 * calc_cube(eigenvalues[i]) - i1 * calc_square(eigenvalues[i]) + i3);

		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				eigenprojections[i * 9 + j * 3 + k] = coefficient * (square_tensor[j * 3 + k] - (i1 - eigenvalues[i]) * A_in[j * 3 + k] + i3 / eigenvalues[i] * identity_tensor[j * 3 + k]);
	}
}


void isotropic_tensor_function_derivative(int rank, double *A_in, double *A_out, double (*function)(const double variable), double (*function_derivative)(const double variable))
{
	const double eigenvalue_tolerance = MATH_TOLERANCE;

	double new_eigenvalues[3];
	double new_eigenvalue_derivatives[3];
	double eigenvalues[3];
	double eigenprojections[27];
	double identity_tensor[9];
	identify_tensor(rank, identity_tensor);

	// Calculate eigenvalues and eigenprojections
	calc_3x3_eigenvalues(eigenvalues, A_in);
	calc_3x3_eigenprojections(eigenprojections, eigenvalues, A_in);

	// Calculate new eigenvalues
	for (int i = 0; i < 3; i++)
		new_eigenvalues[i] = function(eigenvalues[i]);
	for (int i = 0; i < 3; i++)
		new_eigenvalue_derivatives[i] = function_derivative(eigenvalues[i]);

	// Calculate eigenvalue tolerance
	double tolerance = eigenvalue_tolerance * (fabs(eigenvalues[0]) + fabs(eigenvalues[1]) + fabs(eigenvalues[2])) / 3.0;

	// Calculate new tensor when x_1 == x_2 == x_3
	if (fabs(eigenvalues[0] - eigenvalues[1]) <= tolerance && fabs(eigenvalues[1] - eigenvalues[2]) <= tolerance && fabs(eigenvalues[2] - eigenvalues[0]) <= tolerance)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						A_out[i * 27 + j * 9 + k * 3 + l] = new_eigenvalue_derivatives[0] * 0.5 * (identity_tensor[i * 3 + k] * identity_tensor[j * 3 + l] + identity_tensor[i * 3 + l] * identity_tensor[j * 3 + k]);
		return;
	}

	// Calculate new tensor when x_a != x_b == x_c
	for (int m = 0; m < 3; m++)
		if (fabs(eigenvalues[(m + 1) % 3] - eigenvalues[(m + 2) % 3]) <= tolerance)
		{
			const double s1 = (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3]) / calc_square(eigenvalues[m] - eigenvalues[(m + 2) % 3]) - new_eigenvalue_derivatives[(m + 2) % 3] / (eigenvalues[m] - eigenvalues[(m + 2) % 3]);
			const double s2 = 2.0 * eigenvalues[(m + 2) % 3] * (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3]) / calc_square(eigenvalues[m] - eigenvalues[(m + 2) % 3]) - (eigenvalues[m] + eigenvalues[(m + 2) % 3]) / (eigenvalues[m] - eigenvalues[(m + 2) % 3]) * new_eigenvalue_derivatives[(m + 2) % 3];
			const double s3 = 2.0 * (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3]) / calc_cube(eigenvalues[m] - eigenvalues[(m + 2) % 3]) - (new_eigenvalue_derivatives[m] + new_eigenvalue_derivatives[(m + 2) % 3]) / calc_square(eigenvalues[m] - eigenvalues[(m + 2) % 3]);
			const double s4 = eigenvalues[(m + 2) % 3] * s3;
			const double s5 = s4;
			const double s6 = calc_square(eigenvalues[(m + 2) % 3]) * s3;

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						for (int l = 0; l < 3; l++)
							A_out[i * 27 + j * 9 + k * 3 + l] = s1 * 0.5 * (identity_tensor[i * 3 + k] * A_in[l * 3 + j] + identity_tensor[i * 3 + l] * A_in[k * 3 + j] + identity_tensor[j * 3 + l] * A_in[i * 3 + k] + identity_tensor[k * 3 + j] * A_in[i * 3 + l]) - s2 * 0.5 * (identity_tensor[i * 3 + k] * identity_tensor[j * 3 + l] + identity_tensor[i * 3 + l] * identity_tensor[j * 3 + k]) - s3 * A_in[i * 3 + j] * A_in[k * 3 + l] + s4 * A_in[i * 3 + j] * identity_tensor[k * 3 + l] + s5 * identity_tensor[i * 3 + j] * A_in[k * 3 + l] - s6 * identity_tensor[i * 3 + j] * identity_tensor[k * 3 + l];
			return;
		}

	// Calculate new tensor when x_1 != x_2 != x_3
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
				{
					double temp = 0.0;
					for (int m = 0; m < 3; m++)
					{
						const double s1 = new_eigenvalues[m] / ((eigenvalues[m] - eigenvalues[(m + 1) % 3]) * (eigenvalues[m] - eigenvalues[(m + 2) % 3]));
						const double s2 = s1 * (eigenvalues[(m + 1) % 3] + eigenvalues[(m + 2) % 3]);
						const double s3 = s1 * ((eigenvalues[m] - eigenvalues[(m + 1) % 3]) + (eigenvalues[m] - eigenvalues[(m + 2) % 3]));
						const double s4 = s1 * (eigenvalues[(m + 1) % 3] - eigenvalues[(m + 2) % 3]);

						temp += s1 * 0.5 * (identity_tensor[i * 3 + k] * A_in[l * 3 + j] + identity_tensor[i * 3 + l] * A_in[k * 3 + j] + identity_tensor[j * 3 + l] * A_in[i * 3 + k] + identity_tensor[k * 3 + j] * A_in[i * 3 + l]) - s2 * 0.5 * (identity_tensor[i * 3 + k] * identity_tensor[j * 3 + l] + identity_tensor[i * 3 + l] * identity_tensor[j * 3 + k]) - s3 * eigenprojections[m * 9 + i * 3 + j] * eigenprojections[m * 9 + k * 3 + l] - s4 * (eigenprojections[((m + 1) % 3) * 9 + i * 3 + j] * eigenprojections[((m + 1) % 3) * 9 + k * 3 + l] - eigenprojections[((m + 2) % 3) * 9 + i * 3 + j] * eigenprojections[((m + 2) % 3) * 9 + k * 3 + l]) + new_eigenvalue_derivatives[m] * eigenprojections[m * 9 + i * 3 + j] * eigenprojections[m * 9 + k * 3 + l];
					}
					A_out[i * 27 + j * 9 + k * 3 + l] = temp;
				}
}

