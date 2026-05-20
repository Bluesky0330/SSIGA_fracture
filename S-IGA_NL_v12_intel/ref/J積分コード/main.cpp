/****************************************************************
 * S-IGA :)
 *
 * ローカルメッシュ同士は原則被りなしと仮定
 * 要素の重なりの判定はガウス点で行っている
 * ガウス積分の積分点数：NG(一部NG_EXTEND)
 * inputファイルを複数読み込む, 1つ目がグローバルパッチ, 2つ目以降がローカルパッチ，最後の引数は変位場のファイル
 *
 * 
 ****************************************************************/

#include <iostream>
#include <fstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// header
#include "_header.hpp"
#include "_main.hpp"
#include "_IIM2D.hpp"
#include "_IIM3D.hpp"
#include "_extension.hpp"

using namespace std;

int main(int argc, char **argv)
{
	information info;
	int alloc_count = 0, glo_var_count = 0;

	// 処理時間計測
	double time;
	chrono::system_clock::time_point start[2], end[2];
	start[0] = chrono::system_clock::now();
	start[1] = chrono::system_clock::now();

	// モデル数のための引数の調整
	cout << "argc: " << argc << endl;
	if (SIGA_DATA_RECEPTION == 0)       // 変位場ファイルのみ除算
		argc--;
	else if (SIGA_DATA_RECEPTION == 1)  // 変位場と要素の重なり情報ファイルを除算
		argc -= 2;

	argc--;             // 実行ファイル分除算する
	Total_mesh = argc;  // 全メッシュ数, グローバルメッシュ: 1, ローカルメッシュ: 2 ~

	// check arguments, IGA: input file = 1, S-IGA: input file > 1
	if (Total_mesh <= 0)
		printf("\nArgument is missing\n\n");
	else if (Total_mesh == 1)
		printf("\nIGA carried out.(No local mesh)\n\n");
	else if (Total_mesh >= 2)
		printf("\nS-IGA carried out.(%d local meshes)\n\n", Total_mesh - 1);

	// memory allocation 0
	Allocation(alloc_count++, &info);  

	// Read file 1st time
	for (int i = 0; i < Total_mesh; i++)
	{
		Get_Input_1(i, argv, &info);
	}

	// memory allocation 1
	Allocation(alloc_count++, &info);

	// Read file 2nd time
	for (int i = 0; i < Total_mesh; i++)
	{
		Get_Input_2(i, argv, &info);
	}
	printf("\n\nFinish Get Input data\n\n");

	// MAX value
	Global_var(glo_var_count++, &info);

	// memory allocation 2
	Allocation(alloc_count++, &info);

	// INC 等の作成
	Make_INC(&info);

	printf("\nFinish Make INC\n\n");

	// set D_MATRIX_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation 3
	Allocation(alloc_count++, &info);

	if (Total_mesh == 1) // IGA
	{
		Preprocessing_IGA(&info);  
		printf("\nFinish Preprocessing IGA\n\n");
	}
	else if (Total_mesh >= 2) // S-IGA
	{   
		if (SIGA_DATA_RECEPTION == 0)
		{
			searchOverlappingEle(&info);  // 重なり要素の検索
			printf("\nFinish search overlapping element\n\n");
			Preprocessing_S_IGA(&info);
		} 
		else if (SIGA_DATA_RECEPTION == 1)
		{
			Substitute_EOI(&info, argv);  // 変位場の代入
			printf("\nFinish Substitute eoi\n\n");  
			Preprocessing_S_IGA(&info);
		}
	}
	printf("\nFinish check Overlapping SIGA\n\n");

	

	// MAX_K_WHOLE_SIZE
	Global_var(glo_var_count++, &info);
	printf("\nFinish Make global second var\n\n");

	// memory allocation 4
	Allocation(alloc_count++, &info);

	// make D matrix
	Make_D_Matrix(&info);  // 弾性マトリクスの作成

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\nPreprocess time: %.3f[s]\n\n", time);

	start[1] = chrono::system_clock::now();

	Make_gauss_array(0, &info);

	// memory allocation 5
	Allocation(alloc_count++, &info);
	
	// substitute displacement
	Substitute_Displacement(&info, argv);
	printf("\nFinish Substitute Displacement\n\n");

	// 破壊力学解析
	printf("Start fracture analysis\n\n");

	chrono::system_clock::time_point start_frac, end_frac;
	start_frac = chrono::system_clock::now();

	if (FRACTURE_ANALYSIS != 1)
	{
		// 2次元破壊力学解析
		if (info.DIMENSION == 2) 
		{
			unique_ptr<CrackDataSetting> crackData = make_unique<CrackDataSetting>(info.Total_Control_Point_to_mesh[Total_mesh], info.DIMENSION);  
			unique_ptr<IIM_2D> IIM = make_unique<IIM_2D>();                                                                                       

			// read and arrange the crack data
			cout << "Start Read crack data" << endl << endl;
			IIM->ReadCrackData(crackData.get(), &info);        // 2次元き裂情報の取得     
			IIM->MakeUnitBasisVector(crackData.get(), &info);  // 座標変換テンソルの計算

			// calculate Stress Intensity Factors (SIFs) values
			IIM->IIMComputation(crackData.get(), &info);       // 応力拡大係数計算


		// 3次元破壊力学解析
		} else if (info.DIMENSION == 3) 
		{
			unique_ptr<MathFormula> math_ = make_unique<MathFormula>();
			// unique_ptr<CrackDataSetting_3D> crackData = make_unique<CrackDataSetting_3D>();
			unique_ptr<CrackDataSetting_3D> crackData = make_unique<CrackDataSetting_3D>(info.DIMENSION);

			// read and arrange the crack data
			crackData->readCrackData();

			//set crack tip patch element
			crackData->setCrackTipPatchElement(&info);

			// calculate Stress Intensity Factors (SIFs) values
			// unique_ptr<IIM_3D> IIM3D = make_unique<IIM_3D>(info.DIMENSION,info.Total_Patch_to_mesh[Total_mesh],math_.get(),crackData.get());
			unique_ptr<IIM_3D> IIM3D = make_unique<IIM_3D>(info.DIMENSION, info.Total_Patch_to_mesh[Total_mesh], crackData->getdomain_types(), math_.get(), crackData.get());
			IIM3D->computeIIM3D(crackData.get(), &info);
		}
		end_frac = chrono::system_clock::now();
		time = (double)(chrono::duration_cast<chrono::milliseconds>(end_frac - start_frac).count()) / 1000.0;
		printf("\nfracture mechanics analysis time: %.3f[s]\n\n", time);
	}


	// き裂進展解析
	if(FRACTURE_ANALYSIS != 0)
	{	
		printf("start crack extension analysis\n\n");
		chrono::system_clock::time_point start_crack_extension, end_crack_extension;
		start_crack_extension = chrono::system_clock::now();

		// 2D
		if (info.DIMENSION == 2)
		{	
			printf("Error: 2D's crack extension system is not supported\n");
		}
		// 3D
		else if (info.DIMENSION == 3) 
		{
			// き裂進展解析の実行
			unique_ptr<MathFormula> math_ = make_unique<MathFormula>();
			unique_ptr<CrackDataSetting_3D> crackData = make_unique<CrackDataSetting_3D>(info.DIMENSION);

			// き裂進展解析情報の読みとり
			crackData->readCrackData();

			//set crack tip patch element
			crackData->setCrackTipPatchElement(&info);
			
			// csvファイルからのSIFの読みとり
			unique_ptr<IIM_3D> IIM3D = make_unique<IIM_3D>(info.DIMENSION, info.Total_Patch_to_mesh[Total_mesh], crackData->getdomain_types(), math_.get(), crackData.get());   // ここ修正
			printf("start read SIF data");
			IIM3D->readSIFsFromCSV("SIFs.csv", crackData.get());

			// き裂前縁情報の読みとり
			printf("start read crack front data");
			IIM3D->readCrackFrontData(crackData.get(), &info);
			
			// き裂進展解析の実行
			FatigueCrackIO::crackExtensionAnalysis(&info, crackData.get(), IIM3D.get());
		}
		end_crack_extension = chrono::system_clock::now();
		time = (double)(chrono::duration_cast<chrono::milliseconds>(end_crack_extension - start_crack_extension).count()) / 1000.0;
		printf("\nfracture mechanics analysis time: %.3f[s]\n\n", time);
	}

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\nAll analysis time: %.3f[s]\n\n", time);

	return 0;
}
