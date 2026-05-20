
/*
補助場計算に使用している変位の厳密解が平面ひずみ条件のみに適用できるものなのでデバッグの必要あり．
*/

#include "_header.hpp"
#include "_IIM2D.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;

void IIM_2D::IIMComputation(CrackDataSetting *crackData, information *info)
{
    int i, j, k, e;
	double M1 = 0.0, M2 = 0.0;
    double onlyAuxM1 = 0.0, onlyAuxM2 = 0.0, JValue = 0.0;
    int non_qzero_element = 0;

    Make_gauss_array(0, info);

    for (int re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++) {
        e = info->real_element[re];

        // x'- y'座標における各制御点の仮想き裂進展ベクトルの作成(x-y → x'-y')
        vector<double> virtualCrackExtension(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
        for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++) {
            for (j = 0; j < info->DIMENSION; j++) {
                for (k = 0; k < info->DIMENSION; k++) 
                {   
                    virtualCrackExtension[i * info->DIMENSION + j] += crackData->virtualCrackExtension_cp[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + k] * crackData->T[k][j];
                }
            }
        }

        // for debug
        // for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++) {
        //     for (j = 0; j < info->DIMENSION; j++) {
        //         if (fabs(virtualCrackExtension[i * info->DIMENSION + j]) > 0) 
        //         {
        //             printf("virtualCrackExtension[%d] = %.15e\n", i, virtualCrackExtension[i * info->DIMENSION + j]);
        //         }
        //     }
        // }

        // for debug
            // for (j = 0; j < info->DIMENSION; j++) 
            // {
            //     if (fabs(virtualCrackExtension[i * info->DIMENSION + j]) > 0) 
            //     {
            //         printf("virtualCrackExtension[%d] = %.15e\n", i, virtualCrackExtension[i * info->DIMENSION + j]);
            //     }
            // }


        // virtualCrackExtension がゼロの要素数を数える
        int zeroFlag = 0;
        for (i = 0; i < info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++) {
            if (fabs(virtualCrackExtension[i]) < MERGE_ERROR) {
                zeroFlag++;
            }
        }

        // VirtualCrackExtensionゼロでない要素が存在すれば処理を行う
        if (zeroFlag < info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]) 
        {   
            non_qzero_element++;
            printf("VirtualCrackExtension is not zero element : %d\n", e);

            for (int N = 0; N < GP_ON_ELEMENT; N++) {
                // initialize the physical quantities at gauss point
                std::unique_ptr<PQSetting> act = std::make_unique<PQSetting>(info->DIMENSION);
                std::unique_ptr<PQSetting> actGloInLoc = std::make_unique<PQSetting>(info->DIMENSION);
                std::unique_ptr<PQSetting> actOverlay = std::make_unique<PQSetting>(info->DIMENSION);
                std::unique_ptr<PQSetting> aux1 = std::make_unique<PQSetting>(info->DIMENSION);
                std::unique_ptr<PQSetting> aux2 = std::make_unique<PQSetting>(info->DIMENSION);

                // calculate the physical quantity at gauss point
                PhysicalQuantity phy;
                phy.MakeDispGradientAtGP(act.get(), e, N, info);  // ガウス点における変位勾配
                phy.MakeStrainAtGP(act.get(), e, N, info);        // ガウス点におけるひずみ
                phy.MakeStressAtGP(act.get(), info);              // ガウス点における応力
                phy.MakeStrainEnergyDensityAtGP(act.get());       // ガウス点におけるひずみエネルギ密度

                // calculate the global patch physical quantity in local coordinate
                phy.CheckGlobalInLocalCoordinates(actGloInLoc.get(), e, N, info);  // ローカルパッチに重なるグローバル要素の判定とグローバル領域におけるローカル座標ガウス点の物理量（変位勾配，ひずみ）を計算
                phy.MakeStressAtGP(actGloInLoc.get(), info);                       // D マトリクスを用いて該当グローバルパッチの応力を計算

                // overlay the global and local patch
                phy.Overlay(actOverlay->Disp_grad, act->Disp_grad, actGloInLoc->Disp_grad, info->DIMENSION * info->DIMENSION);  // ローカルパッチとグローバルパッチにおける変位勾配を重ねる
                phy.Overlay(actOverlay->Strain, act->Strain, actGloInLoc->Strain, D_MATRIX_SIZE);                               // ローカルパッチとグローバルパッチにおけるひずみを重ねる
                phy.Overlay(actOverlay->Stress, act->Stress, actGloInLoc->Stress, D_MATRIX_SIZE);                               // ローカルパッチとグローバルパッチにおける応力を重ねる
                phy.MakeStrainEnergyDensityAtGP(actOverlay.get());                                                              // ローカルパッチとグローバルパッチにおけるひずみエネルギ密度を計算

                // x'-y'座標における q 関数勾配の計算
                vector<double> D_grad(info->DIMENSION * info->DIMENSION * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
                phy.Make_Disp_Gradient_Matrix(e, D_grad, &info->Gxi[N * info->DIMENSION], info);
                vector<double> qFuncGradient(info->DIMENSION * info->DIMENSION, 0.0);
                for (i = 0; i < info->DIMENSION; i++) {
                    for (j = 0; j < info->DIMENSION; j++) {
                        for (k = 0; k < info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++) {
                            qFuncGradient[i * info->DIMENSION + j] += D_grad[(i * info->DIMENSION + j) * MAX_KIEL_SIZE + k] * virtualCrackExtension[k];  // 
                        }
                    }
                }

                // debug (q gradient)
                // for (i = 0; i < info->DIMENSION; i++) {
                //     for (j = 0; j < info->DIMENSION; j++) {
                //         if (fabs(qFuncGradient[i * info->DIMENSION + j]) > 0) {
                //             printf("qFuncGradient[%d][%d] = %.15e\n", i, j, qFuncGradient[i * info->DIMENSION + j]);
                //         }
                //     }
                // }

                //  debug (q_gradiaen_x)
                // for (i = 0; i < info->DIMENSION; i++) {
                //     if (fabs(qFuncGradient[i]) > 0) {
                //         printf("qFuncGradient[%d] = %.15e\n", 0, qFuncGradient[0]);
                //     }
                // }

                //  debug (q_gradiaen_y)
                // for (i = 0; i < info->DIMENSION; i++) {
                //     if (fabs(qFuncGradient[i]) > 0) {
                //         printf("qFuncGradient[%d] = %.15e\n", 1, qFuncGradient[1]);
                //     }
                // }

                // calculate the auxiliary field
                Auxiliary aux;
                bool isMode1 = true;
                aux.MakeAuxiliary(aux1.get(), actOverlay.get(), e, N, crackData->crackTipCoordinates, crackData->T, isMode1, info);  // モード１補助場の物理量計算
                isMode1 = false;
                aux.MakeAuxiliary(aux2.get(), actOverlay.get(), e, N, crackData->crackTipCoordinates, crackData->T, isMode1, info);  // モード２補助場の物理量計算

                // calculate EMT for Interaction Integral method
                M1 += CalculateEMT(actOverlay.get(), aux1.get(), qFuncGradient, info) * info->w[N] * info->Jac[e * GP_ON_ELEMENT + N];  // モード1実場
                M2 += CalculateEMT(actOverlay.get(), aux2.get(), qFuncGradient, info) * info->w[N] * info->Jac[e * GP_ON_ELEMENT + N];  // モード2実場

                // calculate only auxiliary EMT for J integral method
                bool isAuxiliary = true;
                onlyAuxM1 += CalculateEMT(aux1.get(), qFuncGradient, isAuxiliary, info) * info->w[N] * info->Jac[e * GP_ON_ELEMENT + N];  // モード1補助場
                onlyAuxM2 += CalculateEMT(aux2.get(), qFuncGradient, isAuxiliary, info) * info->w[N] * info->Jac[e * GP_ON_ELEMENT + N];  // モード2補助場

                // calculate EMT for J integral method
                isAuxiliary = false;
                JValue += CalculateEMT(actOverlay.get(), qFuncGradient, isAuxiliary, info) * info->w[N] * info->Jac[e * GP_ON_ELEMENT + N];  // J積分
            }
        }
    }
    printf("Total non q-zero element : %d\n", non_qzero_element);

    // calculate SIFs
    double SIFMode1 = 0.0, SIFMode2 = 0.0, SIFauxMode1 = 0.0, SIFauxMode2 = 0.0, SIF = 0.0;
    bool isIIM = true;

    // debug
    // printf("\n**** M1, M2, onlyAuxM1, onlyAuxM2, JValue ****\n");
    // printf("M1 = %.15e\n", M1);
    // printf("M2 = %.15e\n", M2);
    // printf("onlyAuxM1 = %.15e\n", onlyAuxM1);
    // printf("onlyAuxM2 = %.15e\n", onlyAuxM2);
    // printf("JValue = %.15e\n\n", JValue);

    // モデル形式の設定（フルモデルか1/4モデルか）
    double model_size_para = 0;
    if (MODEL_SIZE_PARA_2D == 0) 
    {   
        printf("\nAnalysis Model : Full Model\n");
        model_size_para = 2.0;
    }
    if (MODEL_SIZE_PARA_2D == 1)
    {   
        printf("\nAnalysis Model : 1/4 Model\n");
        model_size_para = 1.0;
    }

    // SIF計算
    SIFMode1 = CalculateSIF(crackData, &M1, model_size_para, isIIM);            // モード1 SIF計算
    SIFMode2 = CalculateSIF(crackData, &M2, model_size_para, isIIM);            // モード2 SIF計算
    isIIM = false;
    SIFauxMode1 = CalculateSIF(crackData, &onlyAuxM1, model_size_para, isIIM);  // モード1補助場 SIF計算
    SIFauxMode2 = CalculateSIF(crackData, &onlyAuxM2, model_size_para, isIIM);  // モード2補助場 SIF計算
    SIF = CalculateSIF(crackData, &JValue, model_size_para, isIIM);             // J積分 SIF計算

    printf("\n**** J Integral value ****\n");
    printf("SIF = %.15e\n", SIF);
	printf("****                                      ***\n\n");

    printf("**** J Integral value only auxiliary field ****\n");
    printf("SIF_onlyAux_mode1 = %.15e\n", SIFauxMode1);
    printf("SIF_onlyAux_mode2 = %.15e\n", SIFauxMode2);
	printf("****                                      ***\n\n");

    printf("**** Interaction Integral Value ****\n");
    printf("SIF mode1 = %.15e\n", SIFMode1);
    printf("SIF mode2 = %.15e\n", SIFMode2);    
	printf("****                                      ***\n\n");
}

void IIM_2D::GaussianElimination(vector<double> &A, int num)
{
    int i, j, k;
    int pivot_row = 0;
    double tmp;
    vector<double> C(num * (num * 2), 0.0);

    // make enlarged coefficient matrix
    for (i = 0; i < num; i++) {
        for (j = 0; j < num; j++) {
            C[i * (num * 2) + j] = A[i * num + j];
            C[i * (num * 2) + (i + num)] = 1.0;
        }
    }

    // forward elimination
    for (i = 0; i < num; i++) {
        // search for maximum value
        double big = 0.0;
        for (j = i; j < num; j++) {
            if (fabs(C[j * (num * 2) + i]) > big) {
                big = fabs(C[j * (num * 2) + i]);
                pivot_row = j;
            } 
        }
        if (big == 0.0) {
            printf("Error(Gaussian Elimination):underspecified\n");
            exit(1);
        }
        // switching row
        if (i != pivot_row) {
            for (j = num * 2 - 1; j >= 0; j--) {
                tmp = C[i * (num * 2) + j];
                C[i * (num * 2) + j] = C[pivot_row * (num * 2) + j];
                C[pivot_row * (num * 2) + j] = tmp;
            }
        }
        // diagonal component = 1
        for (j = num * 2 - 1; j >= 0; j--) {
            C[i * (num * 2) + j] = C[i * (num * 2) + j] / C[i * (num * 2) + i];
        }
        // pivot column = 0
        for (j = i + 1; j < num; j++) {
            for (k = num * 2 - 1; k >= 0; k--) {
                C[j * (num * 2) + k] -= C[j * (num * 2) + i] * C[i * (num * 2) + k];
            }
        }
    }
    
    // back substitution
    for (i = num - 1; i >= 0; i--) {
        for (j = i - 1; j >= 0; j--) {
            for (k = num * 2; k >= i; k--) {
                C[j * (num * 2) + k] -= C[i * (num * 2) + k] * C[j * (num * 2) + i];
            }
        }
    }
    
    // make convert matrix
    for (i = 0; i < num; i++) {
        for (j = 0; j < num; j++) {
            A[i * num + j] = C[i * (num * 2) + (j + num)];
        }
    }
}

void IIM_2D::AppendDataToCSV(string fileName, double data)
{
	// open file
	std::ofstream csvFile;
	csvFile.open(fileName, std::ios::app);	// append mode

	if (!csvFile.is_open()) {	// if file is not open
		cout << "Make new csv file" << endl;
		csvFile.open(fileName, std::ios::out);
		if (!csvFile.is_open()) {
			cout << "Error: cannot open file" << endl;
			exit(1);
		}
	}

	// append data to csv file
	csvFile << std::setprecision(16) << data << endl;

	// close file
	csvFile.close();
}


void IIM_2D::ReadCrackData(CrackDataSetting *crackData, information *info)
{
    int i, j;
    int cpNum, cpID;

    // read a crack data file
    std::ifstream crackInputFile("J_int.txt");
    if (!crackInputFile.is_open()) {
        cout << "Error opening file" << endl;
        exit(1);
    }

    // crack origin coordinates
    for (i = 0; i < info->DIMENSION; i++) {
        crackInputFile >> crackData->crackOriginCoordinates[i];
        printf("crackOriginCoordinates[%d] = %.15e\n", i, crackData->crackOriginCoordinates[i]);
    }
    cout << endl;

    // crack tip coordinates
    for (i = 0; i < info->DIMENSION; i++) {
        crackInputFile >> crackData->crackTipCoordinates[i];
        printf("crackTipCoordinates[%d] = %.15e\n", i, crackData->crackTipCoordinates[i]);
    }
    cout << endl;

    // virtual crack extension value
    crackInputFile >> crackData->delta;
    printf("delta = %.15e\n\n", crackData->delta);

    // q function
    crackInputFile >> cpNum;
    printf("Number of cp given q-function : %d\n\n", cpNum);
    for (i = 0; i < cpNum; i++) {
        crackInputFile >> cpID;
        for (j = 0; j < info->DIMENSION; j++) {
            crackInputFile >> crackData->virtualCrackExtension_cp[(cpID + info->Total_Control_Point_to_mesh[Total_mesh - 1]) * info->DIMENSION + j];
        }
    }
    crackInputFile.close();
}

void IIM_2D::MakeUnitBasisVector(CrackDataSetting *crackData, information *info)
{
    int i;

    vector<double> unitBasisLocal(info->DIMENSION, 0.0);
    double r_tipRoot = 0.0;
    for (i = 0; i < info->DIMENSION; i++) {
        r_tipRoot += (crackData->crackTipCoordinates[i] - crackData->crackOriginCoordinates[i])
                 * (crackData->crackTipCoordinates[i] - crackData->crackOriginCoordinates[i]);
    }
    double r_tip = sqrt(r_tipRoot);  // r_tip = (x_crack_length^2 + y_crack_length^2)^(1/2)

    // local basis vector (x'-y' coordinates)
    for (i = 0; i < info->DIMENSION; i++) {
        unitBasisLocal[i] = (crackData->crackTipCoordinates[i] - crackData->crackOriginCoordinates[i]) / r_tip;  // (cosθ, sinθ)
    } 

    if (info->DIMENSION == 2) {
        // T matrix (2x2)
        crackData->T[0][0] = unitBasisLocal[0];   crackData->T[0][1] = unitBasisLocal[1];  // (cosθ, sinθ)
        crackData->T[1][0] = -unitBasisLocal[1];  crackData->T[1][1] = unitBasisLocal[0];  // (-sinθ, cosθ)
    }
}

// 補助場もしくはJ積分のための∫_A (σ_ij * ∂u_j/∂x_i - W * δ_1i) * ∂q/∂x_i dA の計算
double IIM_2D::CalculateEMT(PQSetting *act, vector<double> &qFuncGradient, bool isAuxiliary, information *info)
{
    double StrainEnergyDensity;
    vector<double> EMT(info->DIMENSION * info->DIMENSION);
    double value = 0.0;
    
    // true: only auxiliary field, false: global field (J-integral)
    if (isAuxiliary) {
        StrainEnergyDensity = act->StrainEnergyDensity_aux;
    } else {
        StrainEnergyDensity = act->StrainEnergyDensity;
    }

    // (σ_ij * ∂u_j/∂x_i - W * δ_1i)
    EMT[0] = StrainEnergyDensity - (act->Stress[0] * act->Disp_grad[0] + act->Stress[2] * act->Disp_grad[2]);  // 
    EMT[1] =                     - (act->Stress[0] * act->Disp_grad[1] + act->Stress[2] * act->Disp_grad[3]);  // 
    EMT[2] =                     - (act->Stress[2] * act->Disp_grad[0] + act->Stress[1] * act->Disp_grad[2]);  // 
    EMT[3] = StrainEnergyDensity - (act->Stress[2] * act->Disp_grad[1] + act->Stress[1] * act->Disp_grad[3]);  // 

    value -= EMT[0] * qFuncGradient[0] + EMT[2] * qFuncGradient[1];  // 
    value -= EMT[1] * qFuncGradient[2] + EMT[3] * qFuncGradient[3];  // 

    return value;
}

// 相互積分のための積分計算
double IIM_2D::CalculateEMT(PQSetting *act, PQSetting *aux, vector<double> &qFuncGradient, information *info)
{
    vector<double> EMT(info->DIMENSION * info->DIMENSION);
    double value = 0.0;
    
    // 
    EMT[0] = aux->StrainEnergyDensity - (act->Stress[0] * aux->Disp_grad[0] + aux->Stress[0] * act->Disp_grad[0] + act->Stress[2] * aux->Disp_grad[2] + aux->Stress[2] * act->Disp_grad[2]);  // 
    EMT[1] =                          - (act->Stress[0] * aux->Disp_grad[1] + aux->Stress[0] * act->Disp_grad[1] + act->Stress[2] * aux->Disp_grad[3] + aux->Stress[2] * act->Disp_grad[3]);  //
    EMT[2] =                          - (act->Stress[2] * aux->Disp_grad[0] + aux->Stress[2] * act->Disp_grad[0] + act->Stress[1] * aux->Disp_grad[2] + aux->Stress[1] * act->Disp_grad[2]);  // 
    EMT[3] = aux->StrainEnergyDensity - (act->Stress[2] * aux->Disp_grad[1] + aux->Stress[2] * act->Disp_grad[1] + act->Stress[1] * aux->Disp_grad[3] + aux->Stress[1] * act->Disp_grad[3]);  // 

    value -= EMT[0] * qFuncGradient[0] + EMT[2] * qFuncGradient[1];
    value -= EMT[1] * qFuncGradient[2] + EMT[3] * qFuncGradient[3];

    return value;
}

double IIM_2D::CalculateSIF(CrackDataSetting *crackData, double *value, double modelSize, bool isIIM)
{
    double SIF;

    if (DM == 0) {
        // true: Interaction Integral, false: J Integral
        if(isIIM) {
            SIF = *value / crackData->delta / (modelSize * (1.0 / E));  // K = J*E/ 2δa
        } else {
            SIF = sqrt(modelSize * *value * E / crackData->delta);  // K = √(J*E/δa)
        }
    } else if (DM == 1) {
        if (isIIM) {
            SIF = *value / crackData->delta / (modelSize * ((1.0 - nu * nu) / E)); // K = J*E/ (1 - ν^2)δa
        } else {
            SIF = sqrt(modelSize * *value * E / (1.0 - nu * nu) / crackData->delta);  // K = √J*E/δa(1-ν^2)
        }
    }

    return SIF;
}

void Auxiliary::MakeAuxiliary(PQSetting *aux, PQSetting *actOverlay, int e, int GPNum, vector<double> &crackTipCoordinates, double T[MAX_DIMENSION][MAX_DIMENSION], bool isMode1, information *info)
{
    int i, j;
    double mu = E / (2.0 * (1.0 + nu));
    vector<double> crackTipCoords_local(info->DIMENSION, 0.0);
    vector<double> data_result_shape_local(info->DIMENSION, 0.0);

    Make_D_Matrix(info);

    vector<double> Dinv(D_MATRIX_SIZE * D_MATRIX_SIZE, 0.0);  // D_vin : D 逆行列
    for (i = 0; i < D_MATRIX_SIZE; i++) {
        for (j = 0; j < D_MATRIX_SIZE; j++) {
            Dinv[i * D_MATRIX_SIZE + j] = info->D[i * D_MATRIX_SIZE + j];
        }
    }

    // calculate the inverse matrix
    GaussianElimination(Dinv, D_MATRIX_SIZE);  // D の逆行列をDinvとして計算

    // calculate physical coordinates at Gaussian points for each element
    vector<double> data_result_shape(info->DIMENSION, 0.0);
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; i++) {
        double R = Shape_func(i, &info->Gxi[GPNum * info->DIMENSION], e, info);  // ガウス点での基底関数値
        for (j = 0; j < info->DIMENSION; j++) {
            data_result_shape[j] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];  // ガウス点物理座標取得
        }
    }

    // Transform physical coordinates from x-y to x'-y' coordinates
    for (i = 0; i < info->DIMENSION; i++) {
        for (j = 0; j < info->DIMENSION; j++) {
            crackTipCoords_local[i] += T[i][j] * crackTipCoordinates[j];
            data_result_shape_local[i] += T[i][j] * data_result_shape[j];
        }
    }

    // Transform x'-y' coordinates at Gaussian points into Polar coordinates (r, θ) at the crack tip
    double r_root = 0.0, r = 0.0, rad = 0.0;
    for (i = 0; i < info->DIMENSION; i++) {
        r_root += (data_result_shape_local[i] - crackTipCoords_local[i]) * (data_result_shape_local[i] - crackTipCoords_local[i]);  // ガウス点でのx'-y'座標からき裂先端座標までの距離の2乗 r^2
    }
    r = sqrt(r_root);  // x'-y'におけるガウス点からき裂先端座標までの距離 r
    rad = atan2(data_result_shape_local[1] - crackTipCoords_local[1], data_result_shape_local[0] - crackTipCoords_local[0]);  //  ガウス点からき裂先端までのラジアン θ

    // calculate the physical quantities in auxiliary field (true: mode I, false: mode II)
    if (isMode1) {
        MakeDispGradient_aux1(aux, r, rad, mu);  // き裂先端を原点とした極座標系でのモード１補助場の変位勾配
        MakeStress_aux1(aux, r, rad);            // き裂先端を原点とした極座標系でのモード１補助場の応力
    } else {
        MakeDispGradient_aux2(aux, r, rad, mu);  // き裂先端を原点とした極座標系でのモード２補助場の変位勾配
        MakeStress_aux2(aux, r, rad);            // き裂先端を原点とした極座標系でのモード２補助場の応力
    }

    // transform from x'-y' to x-y coordinates
    vector<double> temp_dispGradient(info->DIMENSION * info->DIMENSION, 0.0);
    CoordTransform(T, aux->Disp_grad, temp_dispGradient, info);
	aux->Disp_grad[0] = temp_dispGradient[0];  // 
	aux->Disp_grad[1] = temp_dispGradient[1];  // 
	aux->Disp_grad[2] = temp_dispGradient[2];  // 
	aux->Disp_grad[3] = temp_dispGradient[3];  // 

    vector<double> temp_stress_local(info->DIMENSION * info->DIMENSION, 0.0);
    vector<double> temp_stress(info->DIMENSION * info->DIMENSION, 0.0);
	temp_stress_local[0] = aux->Stress[0];  // σ'_xx
	temp_stress_local[1] = aux->Stress[2];  // τ'_xy
	temp_stress_local[2] = aux->Stress[2];  // τ'_yx
	temp_stress_local[3] = aux->Stress[1];  // σ'_yy
    CoordTransform(T, temp_stress_local, temp_stress, info);
	aux->Stress[0] = temp_stress[0];  // σ_xx
	aux->Stress[2] = temp_stress[1];  // τ_xy
	aux->Stress[1] = temp_stress[3];  // σ_yy

    for (i = 0; i < D_MATRIX_SIZE; i++) {
        for (j = 0; j < D_MATRIX_SIZE; j++) {
            aux->Strain[i] += Dinv[i * D_MATRIX_SIZE + j] * aux->Stress[j];  // x-y 座標系での補助場のひずみ
        }
    }

    for (i = 0; i < D_MATRIX_SIZE; i++) {
        aux->StrainEnergyDensity += actOverlay->Stress[i] * aux->Strain[i];       
        aux->StrainEnergyDensity_aux += (aux->Stress[i] * aux->Strain[i]) / 2.0;  // x-y 座標系での補助場のひずみエネルギー密度
    }
}

// 物理量をx'-y'座標系からx-y座標系に変換する
void Auxiliary::CoordTransform(double T[MAX_DIMENSION][MAX_DIMENSION],  vector<double> &S, vector<double> &result, information *info)
{
    int i, j, k;

    if (info->DIMENSION == 2) {
        vector<double> temp(info->DIMENSION * info->DIMENSION, 0.0);
        for (i = 0; i < info->DIMENSION; i++) {
            for (j = 0; j < info->DIMENSION; j++) {
                for (k = 0; k < info->DIMENSION; k++) {
                    temp[i * info->DIMENSION + j] += T[k][i] * S[k * info->DIMENSION + j];  
                }
            }
        }
        for (i = 0; i < info->DIMENSION; i++) {
            for (j = 0; j < info->DIMENSION; j++) {
                for (k = 0; k < info->DIMENSION; k++) {
                    result[i * info->DIMENSION + j] += temp[i * info->DIMENSION + k] * T[k][j];  // [result] = ([T]^T[S])[T] 
                }
            }
        }
    }
}

void Auxiliary::MakeDispGradient_aux1(PQSetting *aux1, double r, double rad, double mu)
{   
    // if (DM == 0) 
    // {
        // double kappa = 3.0 - 4.0 * nu;
    // }
    // else 
    // {
        // double kappa = (3.0 - nu) / (1 + nu);
    // }
    
    // Irwin の厳密化を偏微分（平面ひずみ条件のみ）
    aux1->Disp_grad[0] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (1.0 - 2.0 * nu + sin(rad / 2.0) * sin(rad / 2.0) - 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
	aux1->Disp_grad[1] = sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * (1.0 - 2.0 * nu + sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * cos(rad / 2.0) * cos(rad / 2.0) * cos(rad)) / (2.0 * mu);
	aux1->Disp_grad[2] = -sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * (2.0 - 2.0 * nu - cos(rad / 2.0) * cos(rad / 2.0) + 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
	aux1->Disp_grad[3] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (2.0 - 2.0 * nu - cos(rad / 2.0) * cos(rad / 2.0) + 2.0 * sin(rad / 2.0) * sin(rad / 2.0) * cos(rad)) / (2.0 * mu);
}

void Auxiliary::MakeDispGradient_aux2(PQSetting *aux2, double r, double rad, double mu)
{
	aux2->Disp_grad[0] = -sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * (2.0 - 2.0 * nu + cos(rad / 2.0) * cos(rad / 2.0) - 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
	aux2->Disp_grad[1] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (2.0 - 2.0 * nu + cos(rad / 2.0) * cos(rad / 2.0) - 2.0 * sin(rad / 2.0) * sin(rad / 2.0) *  cos(rad)) / (2.0 * mu);
	aux2->Disp_grad[2] = -sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (1.0 - 2.0 * nu - sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
	aux2->Disp_grad[3] = -sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * (1.0 - 2.0 * nu - sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * cos(rad / 2.0) * cos(rad / 2.0) * cos(rad)) / (2.0 * mu);
}

void Auxiliary::MakeStress_aux1(PQSetting *aux1, double r, double rad)
{
	aux1->Stress[0] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (1.0 - sin(rad / 2.0) * sin(3.0 * rad / 2.0));
	aux1->Stress[1] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (1.0 + sin(rad / 2.0) * sin(3.0 * rad / 2.0));
	aux1->Stress[2] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * sin(rad / 2.0) * cos(3.0 * rad / 2.0);
}

void Auxiliary::MakeStress_aux2(PQSetting *aux2, double r, double rad)
{
	aux2->Stress[0] = -sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * (2.0 + cos(rad / 2.0) * cos(3.0 * rad / 2.0));
	aux2->Stress[1] = sqrt(1.0 / (2.0 * M_PI * r)) * sin(rad / 2.0) * cos(rad / 2.0) * cos(3.0 * rad / 2.0);
	aux2->Stress[2] = sqrt(1.0 / (2.0 * M_PI * r)) * cos(rad / 2.0) * (1.0 - sin(rad / 2.0) * sin(3.0 * rad / 2.0));
}

void PhysicalQuantity::CheckGlobalInLocalCoordinates(PQSetting *globalInLoc, int elementID, int GPNum, information *info)
{
    if (info->Element_mesh[elementID] > 0)  // ローカルパッチのみ処理の対処
    {
        // eoiを使用して重複要素をチェック
        const auto& overlapping_elements = info->eoi[elementID];
        
        if (!overlapping_elements.empty())  // ローカルメッシュに重なるグローバルパッチがあれば処理
        {
            for (int overlapping_element : overlapping_elements)  // ローカル要素に重なるグローバル要素分だけ処理
            {
                // ローカル要素ガウス点座標におけるグローバル要素の変位勾配とひずみを計算
                MakeDispGradientAndStrain(globalInLoc, elementID, overlapping_element, GPNum, info);
            }
        }
    }
}



void PhysicalQuantity::MakeDispGradientAndStrain(PQSetting *globalInLoc, int localElemID, int globalElemID, int GPNum, information *info)
{
    int i, j, k;
    int BDBJFlag = 0;

    vector<double> U(MAX_KIEL_SIZE);
    vector<double> B_grad_glo (info->DIMENSION * info->DIMENSION * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
    vector<double> B_glo(D_MATRIX_SIZE * MAX_KIEL_SIZE);
    double G_Gxi[MAX_POW_NG_EXTEND][MAX_DIMENSION];

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[globalElemID]]; i++) {
		for (j = 0; j < info->DIMENSION; j++) {
			U[i * info->DIMENSION + j] = info->Displacement[info->Controlpoint_of_Element[globalElemID * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j];
		}
	}

	// determine if Local Gaussian points are included in the Global element
	// calculate physical coordinates of Local element at Gaussian point
	double data_result_shape[2] = {0.0};
	for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[localElemID]]; j++) {
		double R_shape_func = Shape_func(j, &info->Gxi[GPNum * info->DIMENSION], localElemID, info);
		for (k = 0; k < info->DIMENSION; k++) {
			data_result_shape[k] += R_shape_func * info->Node_Coordinate[info->Controlpoint_of_Element[localElemID * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
		}
	}

	// calculate parameter space coordinates of Local element at Gaussian points on the Global patch
	double para[MAX_DIMENSION] = {0.0};
	int patch_n = 0;
	for (j = 0; j < info->Total_Patch_on_mesh[0]; j++) {
		int iteration = Calc_xi_eta(data_result_shape[0], data_result_shape[1], 
                                    &info->Position_Knots[info->Total_Knot_to_patch_dim[j * info->DIMENSION + 0]], &info->Position_Knots[info->Total_Knot_to_patch_dim[j * info->DIMENSION + 1]],
                                    info->No_Control_point[j * info->DIMENSION + 0], info->No_Control_point[j * info->DIMENSION + 1], 
                                    info->Order[j * info->DIMENSION + 0], info->Order[j * info->DIMENSION + 1],
                                    &para[0], &para[1], j, info);
		patch_n = j;

		if (iteration == ERROR) {
			printf("Not converged\n");
			exit(1);
		}
	}

	//要素内外判定
    if (info->DIMENSION == 2) {
        if (para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElemID * info->DIMENSION + 0] + 0] &&
            para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElemID * info->DIMENSION + 0] + 1] &&
            para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElemID * info->DIMENSION + 1] + 0] &&
            para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElemID * info->DIMENSION + 1] + 1])	//要素内であるとき
        {
            BDBJFlag = 1;

            //親要素座標の算出
            G_Gxi[GPNum][0] = -1.0 + 2.0 * (para[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElemID * info->DIMENSION + 0] + 0])
                             / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElemID * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElemID * info->DIMENSION + 0] + 0]);
            G_Gxi[GPNum][1] = -1.0 + 2.0 * (para[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElemID * info->DIMENSION + 1] + 0])
                             / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElemID * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElemID * info->DIMENSION + 1] + 0]);
        }
        else	// 要素外であるとき
        {
            BDBJFlag = 0;
        }
    }

	// 要素内であるとき、次を計算
	if (BDBJFlag) {
		// 重なるグローバル要素のDとBマトリックス
		Make_Disp_Gradient_Matrix(globalElemID, B_grad_glo, G_Gxi[GPNum], info);
        Make_B_Matrix_anypoint(globalElemID, B_glo, G_Gxi[GPNum], info);

        // 重なるグローバル要素の変位勾配を計算
		for (j = 0; j < info->DIMENSION * info->DIMENSION; j++) {
			for (k = 0; k < info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[globalElemID]]; k++) {
				globalInLoc->Disp_grad[j] += B_grad_glo[j * MAX_KIEL_SIZE + k] * U[k];
			}
		}

        // 重なるグローバル要素のひずみを計算
		for (j = 0; j < D_MATRIX_SIZE; j++) {
			int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[globalElemID]] * info->DIMENSION;
			for (k = 0; k < KIEL_SIZE; k++) {
				globalInLoc->Strain[j] += B_glo[j * MAX_KIEL_SIZE + k] * U[k];
			}
		}
	}
}

void PhysicalQuantity::Make_B_Matrix_anypoint(int elementID, vector<double> &B, double *localCoords, information *info)
{
	int i, j, k;
	static vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);

	if (info->DIMENSION == 2) {
        double a_2x2[2][2] = {0.0};
		for (i = 0; i < info->DIMENSION; i++) {
			for (j = 0; j < info->DIMENSION; j++) {
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; k++)
					a_2x2[i][j] += dShape_func(k, j, localCoords, elementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];  
			}
		}

		InverseMatrix_2x2(a_2x2);

		for (i = 0; i < info->DIMENSION; i++) {
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; j++)	{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++) {
                    b[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k][i] * dShape_func(j, k, localCoords, elementID, info);
                }
			}
		}

		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++)	
        {
			B[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];  // ∂Ni/∂x
			B[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;

			B[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			B[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i]; // ∂Ni/∂y

			B[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i]; // ∂Ni/∂y
			B[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i]; // ∂Ni/∂x
		}
        
    } else if (info->DIMENSION == 3) {
        printf("Error: 3D is not supported\n");
        exit(1);
    }
}

void PhysicalQuantity::Make_Disp_Gradient_Matrix(int elementID, vector<double> &D_grad, double *localCoords, information *info)
{
    int i, j, k;
    vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);

    if (info->DIMENSION == 2) {
        double a[2][2];  // ヤコビ行列の成分
        for (i = 0; i < info->DIMENSION; i++) {
            for (j = 0; j < info->DIMENSION; j++) {
                a[i][j] = 0.0;
                for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; k++) {
                    a[i][j] += dShape_func(k, j, localCoords, elementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];  // (∂x/∂xi, ∂y/∂xi), (∂x/∂eta, ∂y/∂eta)
                }
            }
        }

        InverseMatrix_2x2(a);  // ヤコビ行列の逆行列を計算

        for (i = 0; i < info->DIMENSION; i++) {
            for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; j++) {
                b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
                for (k = 0; k < info->DIMENSION; k++) {
                    b[i * MAX_NO_CP_ON_ELEMENT + j] += a[k][i] * dShape_func(j, k, localCoords, elementID, info);
                }
            }
        }

        for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) 
        {
            D_grad[0 * MAX_KIEL_SIZE + (2 * i)] = b[0 * MAX_NO_CP_ON_ELEMENT + i];  // ∂Ni/∂x
            D_grad[0 * MAX_KIEL_SIZE + (2 * i + 1)] = 0.0;

            D_grad[1 * MAX_KIEL_SIZE + (2 * i)] = b[1 * MAX_NO_CP_ON_ELEMENT + i];  // ∂Ni/∂y
            D_grad[1 * MAX_KIEL_SIZE + (2 * i + 1)] = 0.0;

            D_grad[2 * MAX_KIEL_SIZE + (2 * i)] = 0.0;                              
            D_grad[2 * MAX_KIEL_SIZE + (2 * i + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + i];  // ∂Ni/∂y

            D_grad[3 * MAX_KIEL_SIZE + (2 * i)] = 0.0;
            D_grad[3 * MAX_KIEL_SIZE + (2 * i + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + i];  // ∂Ni/∂x
        }
    }
}

void PhysicalQuantity::MakeDispGradientAtGP(PQSetting *act, int elementID, int GPNum, information *info)
{
    int i, j;

    vector<double> U(MAX_KIEL_SIZE);
    vector<double> D_grad(info->DIMENSION * info->DIMENSION * MAX_KIEL_SIZE, 0.0);

    // get the d matrix and the displacement matrix of each element
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) {
        for (j = 0; j < info->DIMENSION; j++) {
            U[i * info->DIMENSION + j] = info->Displacement[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j];  
        }
    }

    // make the displacement gradient
    Make_Disp_Gradient_Matrix(elementID, D_grad, &info->Gxi[GPNum * info->DIMENSION], info);
    for (i = 0; i < info->DIMENSION * info->DIMENSION; i++) {
        for (j = 0; j < info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; j++) {
            act->Disp_grad[i] += D_grad[i * MAX_KIEL_SIZE + j] * U[j];  // ∂u(x)/∂x = ∑ ∂Ni(x)/∂x * ui, ∂u(y)/∂y = ∑ ∂Ni(y)/∂y * ui
        }
    }
}

void PhysicalQuantity::MakeStrainAtGP(PQSetting *act, int elementID, int GPNum, information *info)
{
	int i, j;

	vector<double> U(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	vector<double> B_local(D_MATRIX_SIZE * MAX_KIEL_SIZE);

	//get the B matrix and the displacement of each element
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) {
		for (j = 0; j < info->DIMENSION; j++) {
			U[i * info->DIMENSION + j] = info->Displacement[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j];
		}
	}

	// make strains
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]] * info->DIMENSION;
	Make_B_Matrix_anypoint(elementID, B_local, &info->Gxi[GPNum * info->DIMENSION], info);
	for (i = 0; i < D_MATRIX_SIZE; i++) {
		for (j = 0; j < KIEL_SIZE; j++) {
			act->Strain[i] += B_local[i * MAX_KIEL_SIZE + j] * U[j];  // {ε_x, ε_y, γ_xy} = {Σ ∂Ni/∂x * uxi, Σ ∂Ni/∂y * uyi, (Σ ∂Ni/∂y * uxi + Σ ∂Ni/∂x * uyi)}
		}
	}
}

void PhysicalQuantity::MakeStressAtGP(PQSetting *act, information *info)
{
    for (int i = 0; i < D_MATRIX_SIZE; i++) {
        for (int j = 0; j < D_MATRIX_SIZE; j++) {
            act->Stress[i] += info->D[i * D_MATRIX_SIZE + j] * act->Strain[j]; // {σ}=[D]{ε}
        }
    }
}

void PhysicalQuantity::MakeStrainEnergyDensityAtGP(PQSetting *act)
{
	for (int i = 0; i < D_MATRIX_SIZE; i++) {
		act->StrainEnergyDensity += (act->Stress[i] * act->Strain[i]) / 2.0; // W = 1/2 * ∑ σ_ij * ε_ij
	}
}

void PhysicalQuantity::Overlay(vector<double> &value1, vector<double> &value2, vector<double> &value3, int num)
{
    for (int i = 0; i < num; i++) {
        value1[i] += value2[i] + value3[i]; 
    }
}