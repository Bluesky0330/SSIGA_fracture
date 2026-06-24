// header
#include "_header.hpp"
#include "_main.hpp"
#include <H5Cpp.h>

using namespace std;

int main(int argc, char **argv)
{
	cout << "----- SS-IGA Analysis Program -----" << endl << endl;

	// Avoid HDF5 global atexit finalization crash observed after analysis end.
	try
	{
		H5::H5Library::dontAtExit();
	}
	catch (...)
	{
		// Continue even if the runtime does not support this hook.
	}

	// omp_set_dynamic(1);
	kmp_set_defaults("KMP_BLOCKTIME=0");
	kmp_set_blocktime(0);

	information info;
	int alloc_count = 0, glo_var_count = 0;
	Total_mesh = argc - 2;

	cout << Total_mesh << " input files detected." << endl << endl;

	printf("start Get Option Data\n\n");

	// get option file
	Get_Option(argv[argc - 1], info.opt_files);

	// set option data
	Set_Option_Data(&info);

	// get constant data
	info.c.get_data(info.opt_files[0].c_str(), &info);

	// 処理時間計測
	double time;
	chrono::system_clock::time_point start[2], end[2];
	start[0] = chrono::system_clock::now();
	start[1] = chrono::system_clock::now();

	// check arguments, IGA: input file = 1, S-IGA: input file > 1
	if (Total_mesh <= 0)
		printf("Argument is missing\n\n");
	else if (Total_mesh == 1)
		printf("IGA carried out.(No local mesh)\n\n");
	else if (Total_mesh >= 2)
		printf("SS-IGA carried out.(%d local meshes)\n\n", Total_mesh - 1);

	printf("start Get Input Data\n\n");

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 1st time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_1(i, argv[i + 1], &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 2nd time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_2(i, argv[i + 1], &info);

	// MAX value
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// INC 等の作成
	printf("\nstart Make INC\n");
	if (Total_mesh >= 2)
		Make_INC_for_local_geometric(&info);
	Make_INC(&info);
	
	// set D_MATRIX_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	printf("start Preprocessing\n");
	if (Total_mesh == 1) // IGA
	{
		Make_Gauss_points(true, &info);
	}
	else if (Total_mesh >= 2) // SS-IGA
	{
		searchOverlappingEle(&info);
		Make_Gauss_points(false, &info);
	}

	if (info.c.FRACTURE_MODE == 0 && info.c.CALCLATE_DISPLACEMENT == 0)
	{
		printf("error: No valid analysis mode selected.\n");
		exit(0);
	}
	else if (info.c.CALCLATE_DISPLACEMENT == 1)
	{
		// MAX_K_WHOLE_SIZE
		Global_var(glo_var_count++, &info);

		// memory allocation
		Allocation(alloc_count++, &info);

		// check geometry only output
		bool isGeometryOnly = (info.c.GEOMETRY_ONLY_OUTPUT == 1) ? true : false;
		if (isGeometryOnly && info.c.OUTPUT_GLOBAL_PARAMETERS == 1)
		{
			printf("start Make_connectivity\n");
			Allocation(7, &info);
			Allocation(8, &info);
			Make_connectivity(&info);

			printf("Geometry only output mode in global natural coordinates\n");
			output_global_parameters(&info);
			printf("Finish geometry only output\n");
			return 0;
		}
		if (isGeometryOnly)
		{
			// make connectivity
			printf("start Make_connectivity\n");
			Allocation(7, &info);
			Allocation(8, &info);
			Make_connectivity(&info);

			printf("Geometry only output mode\n\n");
			output_for_paraview_timestep(&info, isGeometryOnly, 0.0);
			printf("Finish geometry only output\n");
			return 0;
		}

		// make K matrix structure
		Make_D_Matrix(&info);
		Make_Index_Dof(&info);
		Make_K_Whole_Ptr_Col(&info, 0);

		// memory allocation
		Allocation(alloc_count++, &info);

		// make K matrix structure
		printf("start Make_K_Whole_Ptr_Col\n");
		Make_K_Whole_Ptr_Col(&info, 1);

		// memory allocation
		Allocation(alloc_count++, &info);

		// make connectivity
		printf("start Make_connectivity\n");
		Allocation(alloc_count++, &info);
		Allocation(alloc_count++, &info);
		Make_connectivity(&info);

		printf("start output_for_paraview_timestep\n");
		output_for_paraview_timestep(&info, isGeometryOnly, 0.0); // start time = 0.0

		// make f vector
		printf("start Make_F_Vec\n");
		Make_F_Vec(&info);
		Make_F_Vec_disp_const(&info);
		Add_Equivalent_Nodal_Force_to_F_Vec(&info);

		// make K matrix
		printf("start calc K\n\n");
		Calc_K(&info);

		end[1] = chrono::system_clock::now();
		time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
		printf("\tPreprocess time:\t%.3f[s]\n\n", time);

		start[1] = chrono::system_clock::now();
		
		printf("start pardiso\n\n");
		intel_PARDISO(info.sol_vec, info.rhs_vec, K_Whole_Size, &info);

		#if 0
		printf("start GMRES solver\n\n");
		GMRES_Solver(K_Whole_Size, 1e-6, &info);
		#endif
		
		// printf("start PCG solver\n\n");
		// PCG_Solver(1000, 1e-6, &info);

		end[1] = chrono::system_clock::now();
		time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
		printf("\tsolver time:\t%.3f[s]\n\n", time);

		start[1] = chrono::system_clock::now();
		printf("start postprocess\n\n");
		Make_Displacement(&info);
		Output_Displacement(&info);
	
		// for output
		if (info.c.CALC_ON_ELE_VERTEX == 1 || info.c.CALC_ON_GP == 1)
			Allocation(9, &info);
		if (info.c.CALC_ON_ELE_VERTEX == 1)
		{
			printf("start calc on element vertex\n\n");
			Calc_on_Element_Vertex(&info);
		}
		if (info.c.CALC_ON_GP == 1)
		{
			printf("start calc on gauss point\n\n");
			Calc_on_Gauss_Point(&info);
		}
		if (info.c.OUTPUT_PARAVIEW == 1)
		{
			printf("start output_for_paraview_timestep\n\n");
			output_for_paraview_timestep(&info, isGeometryOnly, 1.0); // end time = 1.0
		}
	}
	else if (info.c.CALCLATE_DISPLACEMENT == 0)
	{
		printf("Displacement calculation is skipped.\n\n");
		printf("start reading displacement file.\n");

		// Prepare sizes and elastic matrix for fracture/J evaluation path.
		Global_var(glo_var_count++, &info);
		Allocation(4, &info);
		Make_D_Matrix(&info);

		Allocation(7, &info);
		Substitute_Displacement(&info);
	}
	
	if (info.c.FRACTURE_MODE == 1)
	{
		printf("start Fracture analysis\n\n");
		Fracture_analysis(&info);
	}

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\tPostprocess time:\t%.3f[s]\n\n", time);

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\tAll analysis time:\t%.3f[s]\n", time);

	return 0;
}