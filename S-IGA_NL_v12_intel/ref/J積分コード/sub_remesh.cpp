// header
#include "_header.hpp"
#include "_sub.hpp"
#include "_MI3D.hpp"
#include "_extension.hpp"

using namespace std;

#define HARDCODED 1

void remesh(information *info)
{   
    // input関数にする
    int offset = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;              // 前パッチまでの自由度数
    vector<int> patch_remesh;                                                                  // リメッシュするパッチリスト
    vector<int> direction;                                                                     // リメッシュする方向
    vector<int> division_ele_n;                                                                // リメッシュ後の要素数
    vector<int> patch_around_crack;                                                            // き裂周りのパッチリスト
    vector<double> u_vec(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, 0);  // リメッシュ後の変位ベクトル

    #if HARDCODED
        int glo_patch_n = info->Total_Patch_to_mesh[Total_mesh - 1];  // ローカルパッチ番号の初期値

        // set patch_remesh リメッシングするパッチID
        patch_remesh.emplace_back(glo_patch_n + 8);
        patch_remesh.emplace_back(glo_patch_n + 9);
        patch_remesh.emplace_back(glo_patch_n + 12);
        patch_remesh.emplace_back(glo_patch_n + 13);

        // set direction (判定する方向) 0: ξ, 1: η, 2: ζ
        direction.emplace_back(0);
        direction.emplace_back(0);
        direction.emplace_back(0);
        direction.emplace_back(0);

        // set リメッシュ後の要素数 
        division_ele_n.emplace_back(1);
        division_ele_n.emplace_back(1);
        division_ele_n.emplace_back(1);
        division_ele_n.emplace_back(1);

        // set patch_around_crack 
        patch_around_crack.emplace_back(glo_patch_n + 0);
        patch_around_crack.emplace_back(glo_patch_n + 1);
        patch_around_crack.emplace_back(glo_patch_n + 2);
        patch_around_crack.emplace_back(glo_patch_n + 3);
        patch_around_crack.emplace_back(glo_patch_n + 4);
        patch_around_crack.emplace_back(glo_patch_n + 5);
        patch_around_crack.emplace_back(glo_patch_n + 6);
        patch_around_crack.emplace_back(glo_patch_n + 7);
        patch_around_crack.emplace_back(glo_patch_n + 14);
        patch_around_crack.emplace_back(glo_patch_n + 15);
        patch_around_crack.emplace_back(glo_patch_n + 16);
        patch_around_crack.emplace_back(glo_patch_n + 17);
        // patch_around_crack.emplace_back(glo_patch_n + 18);
        // patch_around_crack.emplace_back(glo_patch_n + 19);
        // patch_around_crack.emplace_back(glo_patch_n + 20);
        // patch_around_crack.emplace_back(glo_patch_n + 21);


        // set u_vec
        vector<bool> isSet(info->Total_Control_Point_to_mesh[Total_mesh], false);
        for (size_t i = 0; i < patch_around_crack.size(); i++)
        {
            int p = patch_around_crack[i];
            int cp_num_offset = info->Total_Control_Point_to_patch[p];

            for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
            {
                int connectivity_num = info->Patch_Control_point[cp_num_offset + j];

                if (isSet[connectivity_num])
                    continue;

                for (int k = 0; k < info->DIMENSION; k++)
                {
                    double bias = (k == 2) ? 1.0 : 1.2;  //き裂進展量をバイアスを用いた仮の設定, u_vecには最終的には計算された情報を入れる
                    u_vec[connectivity_num * info->DIMENSION + k] = bias * info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + k];
                }
                isSet[connectivity_num] = true;
            }
        }
        
        for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
        {
            // skip patch_around_crack
            if (find(patch_around_crack.begin(), patch_around_crack.end(), i) != patch_around_crack.end())
                continue;

            for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
            {
                if (isSet[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]])
                    continue;
                int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j];
                for (int k = 0; k < info->DIMENSION; k++)
                    u_vec[cp_num * info->DIMENSION + k] = info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + k];
            }
        }
    #endif

    // update coord
    vector<vector<double>> new_coord(info->Total_Patch_to_mesh[Total_mesh]);
    vector<vector<double>> new_w(info->Total_Patch_to_mesh[Total_mesh]);
    for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
    {
        auto it = find(patch_remesh.begin(), patch_remesh.end(), i);  // リメッシュするパッチかどうか
        int index = distance(patch_remesh.begin(), it);               // リメッシュするパッチのindex
        int p = patch_remesh[index];  
        int dir = direction[index];

        for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
        {
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j];  // 制御点番号を取得
            // 端部の制御点でない場合はスキップ
            if (it != patch_remesh.end() && !(info->INC[p * offset + cp_num * info->DIMENSION + dir] == 0 || info->INC[p * offset + cp_num * info->DIMENSION + dir] == info->No_Control_point[p * info->DIMENSION + dir] - 1))
                continue;
            
            for (int k = 0; k < info->DIMENSION; k++)
                new_coord[i].emplace_back(u_vec[cp_num * info->DIMENSION + k]);
            new_w[i].emplace_back(info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + info->DIMENSION]);
        }
    }

    // remesh
    for (size_t i = 0; i < patch_remesh.size(); i++)
    {
        int p = patch_remesh[i];
    
        // set MAX_INCREASED_CP
        int MAX_INCREASED_CP_total = 1;
        for (int j = 0; j < info->DIMENSION; j++)
        {
            int temp = info->No_Control_point[p * info->DIMENSION + j];
            if (direction[i] == j)
                temp = division_ele_n[i] + info->Order[p * info->DIMENSION + j];
            MAX_INCREASED_CP_total *= temp;
        }

        // set MAX_INCREASED_KV
        int MAX_INCREASED_KV_1D = 0;
        for (int j = 0; j < info->DIMENSION; j++)
        {
            int temp = info->No_knot[p * info->DIMENSION + j];
            if (direction[i] == j)
            {
                temp += division_ele_n[i] + 2 * info->Order[p * info->DIMENSION + j] + 1;
            }
            if (MAX_INCREASED_KV_1D < temp)
                MAX_INCREASED_KV_1D = temp;
        }

        // set info_global
        info_global glo;
        glo.DIMENSION = info->DIMENSION;
        glo.Total_Control_Point = new_coord[p].size() / info->DIMENSION;
        (void)glo.mode;

        vector<double> Weight(MAX_INCREASED_CP_total);
        vector<double> temp_Weight(MAX_INCREASED_CP_total);
        glo.Weight = Weight.data();
        int count = 0;
        for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
        {
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[p] + j];
            if (!(info->INC[p * offset + cp_num * info->DIMENSION + direction[i]] == 0 || info->INC[p * offset + cp_num * info->DIMENSION + direction[i]] == info->No_Control_point[p * info->DIMENSION + direction[i]] - 1))
                continue;

            glo.Weight[count++] = info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + info->DIMENSION];
        }
        glo.temp_Weight = temp_Weight.data();

        // set info_each_DIMENSION
        vector<info_each_DIMENSION> nurbs;
        vector<vector<double>> CP(info->DIMENSION, vector<double>(MAX_INCREASED_CP_total));
        vector<vector<double>> temp_CP(info->DIMENSION, vector<double>(MAX_INCREASED_CP_total));
        vector<vector<double>> KV(info->DIMENSION, vector<double>(MAX_INCREASED_KV_1D));
        vector<vector<double>> temp_KV(info->DIMENSION, vector<double>(MAX_INCREASED_KV_1D));
        for (int j = 0; j < info->DIMENSION; j++)
        {
            info_each_DIMENSION temp;
            temp.CP_n = info->No_Control_point[p * info->DIMENSION + j];
            if (j == direction[i])
                temp.CP_n = 2;
            temp.Order = info->Order[p * info->DIMENSION + j];
            if (j == direction[i])
                temp.Order = 1;
            temp.OE_n = 0;
            if (j == direction[i])
                temp.OE_n = info->Order[p * info->DIMENSION + j] - 1;
            temp.knot_n = info->No_knot[p * info->DIMENSION + j];
            if (j == direction[i])
                temp.knot_n = 4;
            if (j == direction[i])
                temp.KI_cp_n = division_ele_n[i] + info->Order[p * info->DIMENSION + j];
            else
                temp.KI_cp_n = 0;
            temp.KI_non_uniform_n = 0;
            temp.CP = CP[j].data();
            for (size_t k = 0; k < new_coord[p].size() / info->DIMENSION; k++)
                temp.CP[k] = new_coord[p][k * info->DIMENSION + j];
            temp.KV = KV[j].data();
            for (int k = 0; k < info->No_knot[p * info->DIMENSION + j]; k++)
                temp.KV[k] = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + j] + k];
            if (j == direction[i])
            {
                temp.KV[0] = 0.0; temp.KV[1] = 0.0; temp.KV[2] = 1.0; temp.KV[3] = 1.0;
            }
            temp.temp_CP = temp_CP[j].data();
            temp.temp_KV = temp_KV[j].data();
            (void)temp.insert_knot;
            nurbs.emplace_back(temp);
        }
        
        // knot removal and knot insertion
        for (int j = 0; j < info->DIMENSION; j++)
        {
            if (j != direction[i])
                continue;

            // order elevation
            OE(j, &glo, &nurbs[j]);

            // knot insertion
            KI_cp(j, &glo, &nurbs[j]);
        }

        // update new coord
        new_coord[p].resize(glo.Total_Control_Point * info->DIMENSION);
        for (int j = 0; j < glo.Total_Control_Point; j++)
            for (int k = 0; k < info->DIMENSION; k++)
                new_coord[p][j * info->DIMENSION + k] = nurbs[k].CP[j];
    }

    // output as vtk file
    ofstream vtkFile("output.vtk");
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "Remesh output" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET POLYDATA" << endl;

    // calc total control point
    int total_cp = 0;
    for (int i = 1; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
    {
        if (find(patch_remesh.begin(), patch_remesh.end(), i) != patch_remesh.end())
            total_cp += new_coord[i].size() / info->DIMENSION;
        else
            total_cp += info->No_Control_point_in_patch[i];
        cout << i << ": ";
        cout << new_coord[i].size() / info->DIMENSION << endl;
    }

    // Write points
    vtkFile << "POINTS " << total_cp << " float" << endl;

    for (int i = 1; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
    {
        bool isModified = (find(patch_remesh.begin(), patch_remesh.end(), i) != patch_remesh.end());

        // All patch
        int cp_n = isModified ? new_coord[i].size() / info->DIMENSION : info->No_Control_point_in_patch[i];
        // Only Modified patch
        // int cp_n = isModified ? new_coord[i].size() / info->DIMENSION : 0;
        for (int j = 0; j < cp_n; j++)
        {
            for (int k = 0; k < info->DIMENSION; k++)
            {
                vtkFile << new_coord[i][j * info->DIMENSION + k] << " ";
            }
            if (info->DIMENSION == 2)
            {
                vtkFile << "0.0"; // Add z-coordinate for 2D points
            }
            vtkFile << endl;
        }
    }

    // Write vertices
    vtkFile << "VERTICES " << total_cp << " " << total_cp * 2 << endl;
    for (int i = 0; i < total_cp; i++)
    {
        vtkFile << "1 " << i << endl;
    }

    vtkFile.close();

    cout << "finnish remesh" << endl;



    cout << "remesh" << endl;
    // exit(0);

}