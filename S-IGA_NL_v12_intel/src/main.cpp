// header
#include "_header.hpp"
#include "_main.hpp"

using namespace std;

int main(int argc, char **argv)
{
	// omp_set_dynamic(1);
	kmp_set_defaults("KMP_BLOCKTIME=0");
	kmp_set_blocktime(0);

	information info;
	int alloc_count = 0, glo_var_count = 0;
	Total_mesh = argc - 2;

	// get option file data
	Get_Option(argv[argc - 1], info.opt_files);

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
		printf("S-IGA carried out.(%d local meshes)\n\n", Total_mesh - 1);

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

	printf("\n\nFinish Get Input data\n\n");

	// MAX value
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// INC 等の作成
	Make_INC(&info);
	printf("\nFinish Make INC\n\n");

	// set D_MATRIX_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	if (Total_mesh == 1) // IGA
	{
		Make_Gauss_points(true, &info);
	}
	else if (Total_mesh >= 2) // S-IGA
	{
		searchOverlappingEle(&info);
		Make_Gauss_points(false, &info);
	}
	printf("\nFinish Preprocessing\n\n");

	// MAX_K_WHOLE_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make K matrix
	Make_D_Matrix(&info);
	Make_Index_Dof(&info);
	Make_K_Whole_Ptr_Col(&info, 0);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make K matrix
	Make_K_Whole_Ptr_Col(&info, 1);
	printf("\nFinish Make_K_Whole_Ptr_Col\n\n");

	// memory allocation
	Allocation(alloc_count++, &info);

	// make f vector
	Make_F_Vec(&info);
	Add_Equivalent_Nodal_Force_to_F_Vec(&info);
	printf("\nFinish Make_F_Vec\n\n");

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\nPreprocess time: %.3f[s]\n\n", time);

	start[1] = chrono::system_clock::now();
	if (info.c.ANALYSIS_MODE_0 == 0)
		printf("\nStart Linear analysis\n\n");
	else 
		printf("\nStart Non-linear analysis\n\n");

	// Nonlinear loop
	NL_loop(&info, argv[argc - 1], argv[1]);

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\nNon-linear analysis time: %.3f[s]\n\n", time);

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\nAll analysis time: %.3f[s]\n\n", time);

	return 0;
}