/*
ファイル名)
破壊力学解析用のインプットファイル          J_int.txt
応力拡大係数を格納するファイル              SIF.csv
き裂進展前後の前縁点と制御点座標のファイル   Crack_Front_Points_and_CP.csv
繰り返し荷重サイクル数を格納するファイル     Nk.csv

memo)
SIFs.csvの各モードのSIFはdomain_typesと同一の数だけあることを前提とした実装になっている．

き裂進展方向はき裂前縁ベクトルをそのまま格納．キンク等は考慮していない．
最小二乗法近似による制御点座標の計算はパッチごとに行っているが，全体のジオメトリに対して行う方が適切であった．
最小二乗法における形状関数の重みは１としている．
最小二乗法計算における逆行列の計算にEigen のLU分解法を使用．

リメッシング工程はブロック体における埋め込みき裂 or 表面き裂問題を仮定している．
また，き裂前縁方向に45°の角度を持つパッチによるモデリングを仮定し，90°のモデルと180°のモデルに対応している．
45 * 2 = 90°のモデルは未デバッグ．
ローカルパッチ外郭も，き裂前縁の制御点移動量を基に拡張を行っている．（0°と90°でのき裂前縁方向のインデックス）

リメッシング後も要素数とコネクティビティ等が同一であることを仮定した実装．

現在は応力拡大係数振幅が
*/

#include "_header.hpp"
#include "_sub.hpp"
#include "_IIM3D.hpp"
#include "_extension.hpp"
#include "_MI3D.hpp"

#if 1
	// windows, cygwin
	#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
        #include <Eigen/Core>
		#include <Eigen/Dense>
		#include <Eigen/Sparse>
	#endif

	// linux
	#if defined(__linux__)
        #include <eigen3/Eigen/Core>
		#include <eigen3/Eigen/Dense>
		#include <eigen3/Eigen/Sparse>
	#endif
	using namespace Eigen;
#endif


// き裂進展解析
void FatigueCrackIO::crackExtensionAnalysis(information *info, CrackDataSetting_3D *crackData, IIM_3D *IIM3D) 
{   
	// き裂進量計算に関する情報の読みとり
	unique_ptr<ParisLawCrack> paris = make_unique<ParisLawCrack>(*crackData, info, crackData->getdomain_types());
	paris->readParameters("J_int.txt");
	std::cout << "finish read paris law parameters" << std::endl << std::endl;

	// 応力拡大係数振幅
	// unique_ptr<StressIntensityFactorRange> deltaK = make_unique<StressIntensityFactorRange>(*crackData, MODE, DOMAIN_TYPES, 0.0);
	unique_ptr<StressIntensityFactorRange> deltaK = make_unique<StressIntensityFactorRange>(*crackData, MODE, crackData->getdomain_types());
	deltaK->calcDeltaK(IIM3D->getMintegralSIF(), *crackData, *paris);
    std::cout << "finish calc delta K" << std::endl << std::endl;

    // き裂進展速度の計算
    paris->calcCrackGrowthRate(*crackData, *deltaK.get());
    std::cout << "finish calc crack growth rate" << std::endl << std::endl;

    // き裂進展方向の計算
    unique_ptr<CrackDirection> crackExtensionVector = make_unique<CrackDirection>(*crackData, info, crackData->getdomain_types());
    crackExtensionVector->calcCrackGrowthDirection(*crackData, *IIM3D, *paris.get(), info);
    std::cout << "finish calc crack growth direction" << std::endl << std::endl;

    // 荷重繰り返しサイクル数の計算
    paris->calcCrackGrowthCycle(*crackData);
    std::cout << "finish calc crack growth cycle" << std::endl << std::endl;

    // き裂進展ベクトルの計算 
    paris->calcCrackGrowthRateComponent(*crackData, *crackExtensionVector, info);  // き裂進展速度成分の計算
    crackExtensionVector->calcCrackGrowthVector(*crackData, *paris.get(), info);   // き裂進展ベクトルの計算
    std::cout << "finish calc crack growth vector" << std::endl << std::endl;

    // 進展後き裂前縁点の計算
    unique_ptr<CrackFrontManager> crackFrontCoord = make_unique<CrackFrontManager>(*crackData, info, crackData->getdomain_types());
    crackFrontCoord->getbeforeCrackFrontPoints(*IIM3D);                          // 進展前き裂前縁の幾何情報の取得
    crackFrontCoord->calcaferCrackFrontPoints(*crackExtensionVector.get(), *crackData, info);      // 進展後き裂前縁の幾何情報の計算
    std::cout << "finish calc crack front points" << std::endl << std::endl;

    // 進展前後き裂前縁離散点の出力
    unique_ptr<FatigueCrackIO> CrackIO = make_unique<FatigueCrackIO>(0);
    CrackIO->read_select_domain("J_int.txt", *crackData);                    // 進展に基づく積分領域の選択
    CrackIO->writeCrackFrontData(*crackFrontCoord.get(), info);  // き裂前縁離散点の出力
    std::cout << "finish write crack front data" << std::endl << std::endl;

    // き裂前縁をB-スプライン補間により生成
    unique_ptr<LeastSquaresApproximation> calc_cp_LSA = make_unique<LeastSquaresApproximation>(*crackData, info);
    int selected_domain = CrackIO->getSelectDomain();
    bool use_total = true;
    calc_cp_LSA->preLeastSquaresApproximation(*crackData, *IIM3D, info, use_total);                                    // 
    calc_cp_LSA->DoLeastSquaresApproximation(*crackData, *crackFrontCoord, *IIM3D, info, selected_domain, use_total);  // 最小二乗法による制御点座標の算出
    CrackIO->writeCrackFrontCPData(*calc_cp_LSA, info);                                                     // 進展前後制御点座標のcsvファイル出力
    std::cout << "finish calc control point coordinate" << std::endl << std::endl;

    // リメッシング
    unique_ptr<Remeshing> remesh = make_unique<Remeshing>(*calc_cp_LSA, info);
    remesh->DoRemeshing(*crackData, *calc_cp_LSA, info);                     
    std::cout << "finish remesh" << std::endl << std::endl;
}


// パリス則パラメータの読み込み
void ParisLawCrack::readParameters(const std::string& filename) 
{
    std::ifstream file(filename);
    if (!file) 
    {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(1);
    }

    // 積分パラメータをスキップ
    int domainWidth, JintegralPoints, model_size_para_3D, domain_types, domain_area_divisions, patchType;
    int num_edge_crack_patch, num_singular_patch, patchID;

    file >> domainWidth >> JintegralPoints >> model_size_para_3D >> domain_types >> domain_area_divisions >>patchType;
    file >> num_edge_crack_patch >> num_singular_patch;

    for (int i = 0; i < num_edge_crack_patch * num_singular_patch; i++)
    {
        file >> patchID;
    }

    // き裂モードを読みとり (0: mode I, 1: mixed mode)
    int mode;
    file >> mode;
    setCrackMode(mode);

    // Paris則パラメータの読み込み
    double C_val, m_val;
    file >> C_val >> m_val;
    setC(C_val);
    setM(m_val);
	// printf("C: %.5f, m: %.5f\n", C_val, m_val);

    // 最大き裂進展長さの読み込み
    double max_growth;
    file >> max_growth;
    setMaxGrowth(max_growth);
	// printf("max crack growth length Δa: %.5f\n", max_growth);

    // 繰り返し荷重サイクル数 N の計算 (0:進展前のN-a曲線の傾きを使用，1:進展前後のN-a曲線の傾きを平均化)
    // int selectNk;
    // file >> selectNk;
    // setselectedNk(selectNk);

    // き裂進展量に対して荷重移動平均を計算するか
    int weight_move_ave, num_move_ave;
    file >> weight_move_ave >> num_move_ave;
    setWeightedMovingAverage(weight_move_ave);
    setNumberMovingAverage(num_move_ave);
    printf("weighted moving average: %d\n", weighted_moving_average);
    printf("weighted moving average: %d\n", weighted_moving_average);

    // error handring
    if (getCrackMode() != 0 && getCrackMode() != 1) 
    {
        std::cerr << "Error: Invalid crack mode" << std::endl;
        exit(1);
    }
    if (weight_move_ave != 0 && weight_move_ave != 1 && weight_move_ave != 2) 
    {
        std::cerr << "Error: Invalid weighted moving average" << std::endl;
        exit(1);
    }
    
    // 読み込んだ値の表示
	// std::cout << std::fixed << std::setprecision(10);
    std::cout << "Read analysis parameters:" << std::endl;
    std::cout << "crack mode" << getCrackMode() << std::endl;
    std::cout << "Paris law constant C: " << getC() << std::endl;
    std::cout << "Paris law exponent m: " << getM() << std::endl;
	std::cout << "max crack growth length Δa:" << getMaxGrowth() << std::endl;
    // std::cout << "selected Nk: " << selected_Nk << std::endl;
    std::cout << "weighted moving average: " << weighted_moving_average << std::endl;
    std::cout << "number of moving average: " << number_moving_average << std::endl;
    std::cout << std::endl;

    file.close();
}
//SIF振幅を計算
void StressIntensityFactorRange::calcDeltaK(const std::vector<double>& K_max, const CrackDataSetting_3D& crackData, const ParisLawCrack& paris)  
{
    // 引数の配列サイズチェック
	std::cout << "K_max size: " << K_max.size() << std::endl;      
	std::cout << "delta_K size: " << delta_K.size() << std::endl; 

    // SIF配列数の確認
    if ((K_max.size() != delta_K.size()) || K_max.size() != static_cast<std::size_t>(crackData.getJintegralPoints() * MODE * crackData.getdomain_types())) 
    {
        std::cerr << "Error: Size mismatch in calcDeltaK" << std::endl;
        return;
    }

    // モード別のSIF振幅を計算
    for (size_t i = 0; i < K_max.size(); i++) 
    {
        delta_K[i] = K_max[i] - 0.0;
        // setDeltaK(delta_K);
    }
    // for dubug
    #if 1
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        std::cout << "\n=== Calculation Point " << point << " ===" << std::endl;
        for (int mode = 0; mode < MODE; mode++)
        {
            std::cout << "\nMode " << mode + 1 << ":" << std::endl;
            for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
            {
                int index = point * (MODE * crackData.getdomain_types()) + mode * crackData.getdomain_types() + domain;
                std::cout << "Domain Size " << domain + 1 << ": "
                          << "ΔK" << mode + 1 << " = "
                          << std::scientific << std::setprecision(16) 
                        //   << delta_K[index] 
                          << getDeltaK()[index] << std::endl;
            }
        }
    }
    #endif

    // ΔK_bar = ΔKI + |ΔKII| + |ΔKIII| の計算
    calcDeltaKbar(crackData);

    // 等価SIF振幅の計算 (ΔKeq = (1/2)*√(ΔKI^2 + 4(α1ΔKII)^2 + 4(α2ΔKIII)^2)
    calcDeltaKeq(crackData, paris);
}
//ΔK_barを計算
void StressIntensityFactorRange::calcDeltaKbar(const CrackDataSetting_3D& crackData) 
{
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            // モード I, II, III の応力拡大係数を取得
            double deltaKI = delta_K[point * (MODE * crackData.getdomain_types()) + 0 * crackData.getdomain_types() + domain];
            double deltaKII = delta_K[point * (MODE * crackData.getdomain_types()) + 1 * crackData.getdomain_types() + domain];
            double deltaKIII = delta_K[point * (MODE * crackData.getdomain_types()) + 2 * crackData.getdomain_types() + domain];

            // ΔK_bar = ΔKI + |ΔKII| + |ΔKIII| の計算
            int index = point * crackData.getdomain_types() + domain;
            delta_K_bar[index] = deltaKI + std::abs(deltaKII) + std::abs(deltaKIII);
        }
    }

    // for debug
    #if 1
    std::cout << std::endl << "=== ΔK_bar ===" << std::endl;
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        std::cout << "\n=== Calculation Point " << point << " ===" << std::endl;
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            std::cout << "Domain Size " << domain + 1 << ": ΔK_bar = " 
                      << std::scientific << std::setprecision(6) 
                      << getDeltaKbar()[index] << std::endl;
        }
    }
    #endif
}
//等価応力拡大係数振幅ΔK_eqを計算
void StressIntensityFactorRange::calcDeltaKeq(const CrackDataSetting_3D& crackData, const ParisLawCrack& paris) 
{
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            // モード I, II, III の応力拡大係数を取得
            double deltaKI = delta_K[point * (MODE * crackData.getdomain_types()) + 0 * crackData.getdomain_types() + domain];
            double deltaKII = delta_K[point * (MODE * crackData.getdomain_types()) + 1 * crackData.getdomain_types() + domain];
            double deltaKIII = delta_K[point * (MODE * crackData.getdomain_types()) + 2 * crackData.getdomain_types() + domain];

            // ΔK_eq
            int index = point * crackData.getdomain_types() + domain;
            if (paris.getCrackMode() == 0) 
            {
                delta_K_eq[index] = deltaKI;  // mode I
            }
            else if (paris.getCrackMode() == 1)
            {
                delta_K_eq[index] = (deltaKI + sqrt(pow(deltaKI, 2) + 4 * pow(alfa_1 * deltaKII, 2) + 4 * pow(alfa_2 * deltaKIII, 2))) / 2;  // mixed mode (Richard)
            }
        }
    }

    // for debug
    #if 1
    std::cout << "=== ΔK_eq ===" << std::endl;
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        std::cout << "\n=== Calculation Point " << point << " ===" << std::endl;
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            std::cout << "Domain Size " << domain + 1 << ": ΔK_eq = " 
                      << std::scientific << std::setprecision(6) 
                      << getDeltaKeq()[index] << std::endl;
        }
    }
    #endif
}

//き裂進展速度の計算と重み付き移動平均処理
void ParisLawCrack::calcCrackGrowthRate(const CrackDataSetting_3D& crackData, const StressIntensityFactorRange& deltaK) 
{
    const auto& delta_K_eq = deltaK.getDeltaKeq();  // 等価応力拡大係数を取得

    // 各計算点、各領域サイズでの進展速度を計算
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            da_dN[index] = C * pow(delta_K_eq[index], m);  // Paris則: da/dN = C(ΔK_eq)^m
        }
    }

    // for debug
    #if 1
    std::cout << std::endl << "=== da/dN ===" << std::endl;
    for (int point = 0; point < crackData.getJintegralPoints(); point++) 
    {
        std::cout << "\n=== Calculation Point " << point << " ===" << std::endl;
        for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            std::cout << "Domain Size " << domain + 1 
                        << ": da/dN = " << std::scientific << std::setprecision(16) 
                        << getCrackGrowthRate()[index] << " [mm/cycle]" << std::endl;
        }
    }
    #endif

    // き裂進展速度に対して重み付き移動平均化処理
    if (weighted_moving_average == 1 || weighted_moving_average == 2) 
    {   
        cout << "Do Weighted Moving Average" << endl;
        calcWeightedMovingAverageCrackGrowthRate(crackData);
    }
    else 
    {
        cout << "Do not Weighted Moving Average" << endl;
    }

    // 最大進展速度の計算
    calcMaxCrackGrowthRate(crackData);
}
#if 0
//重み付き移動平均
void ParisLawCrack::calcWeightedMovingAverageCrackGrowthRate(const CrackDataSetting_3D& crackData)
{
    std::vector<double> temp_da_dN = da_dN;  // 平滑化前のda/dNを一時保存

    // データ数の確認_き裂先端を基準に対称
    if (number_moving_average % 2 == 0) 
    {
        std::cerr << "Error: number of data points must be odd" << std::endl;
        exit(1);
    }

    //const int half_window = (number_moving_average - 1) / 2;  // 片側のデータ数

    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            double sum = 0.0;
            double weight_sum = 0.0;
            int index = point * crackData.getdomain_types() + domain;
            int scope = 5;
            int half_scope = (scope - 1) / 2;  //中心点を除いた半分の計算幅
            int start = -1;
            int end = -1;

            // 利用可能なデータの範囲を決定
            int N = crackData.getJintegralPoints();
            int max_left = std::min(half_scope, point);
            int max_right = std::min(half_scope, N - 1 - point);
            int actual_half = std::min(max_left, max_right);

            start = point - actual_half;
            end   = point + actual_half;

            //int start = std::max(0, point - half_window);
            //int end = std::min(crackData.getJintegralPoints() - 1, point + half_window);
            std::cout << "start: " << start << ", end: " << end << std::endl;


            // 移動平均の種類に応じた計算
            if (weighted_moving_average == 1)  // 単純移動平均
            {
                for (int i = start; i <= end; i++) 
                {
                    sum += temp_da_dN[i * crackData.getdomain_types() + domain];
                }
                da_dN[index] = sum / (end - start + 1);
            }
            else if (weighted_moving_average == 2)  // 重み付き移動平均
            {
                for (int i = start; i <= end; i++) 
                {
                    // 中心からの距離に応じた重みを計算
                    double sigma = (double) (end - start + 1) / 3.0;
                    int distance = std::abs(i - point);
                    double weight = std::exp(-(distance * distance) / (2.0 * sigma * sigma)); //ガウス型 w_i= exp(-(i-μ)^2/2σ^2)
                    //int weight = number_moving_average - distance;

                    sum += temp_da_dN[i * crackData.getdomain_types() + domain] * weight;
                    weight_sum += weight;
                }
                da_dN[index] = sum / weight_sum;
            }
        }
    }

    #if 1
    // 平滑化の種類に応じたデバッグ出力
    std::cout << std::endl << "=== ";
    if (weighted_moving_average == 1)
        std::cout << "Moving Average";
    else if (weighted_moving_average == 2)
        std::cout << "Weighted Moving Average";
    std::cout << " da/dN ===" << std::endl << std::endl;

    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        std::cout << "\n=== Domain Size " << domain + 1 << " ===" << std::endl;
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            std::cout << "Point " << point 
                     << ": Original = " << std::scientific << std::setprecision(6) 
                     << temp_da_dN[index] 
                     << ", Smoothed = " << da_dN[index] << " [mm/cycle]" << std::endl;
        }
    }
    #endif
}
#endif

void ParisLawCrack::calcWeightedMovingAverageCrackGrowthRate(const CrackDataSetting_3D& crackData)
{
    std::vector<double> temp_da_dN = da_dN;  // 平滑化前のda/dNを一時保存

    // データ数の確認_き裂先端を基準に対称
    if (number_moving_average % 2 == 0) 
    {
        std::cerr << "Error: number of data points must be odd" << std::endl;
        exit(1);
    }

    const int half_window = (number_moving_average - 1) / 2;  // 片側のデータ数

    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            double sum = 0.0;
            double weight_sum = 0.0;
            int index = point * crackData.getdomain_types() + domain;

            // 利用可能なデータの範囲を決定
            int start = std::max(0, point - half_window);
            int end = std::min(crackData.getJintegralPoints() - 1, point + half_window);
            std::cout << "start: " << start << ", end: " << end << std::endl;


            // 移動平均の種類に応じた計算
            if (weighted_moving_average == 1)  // 単純移動平均
            {
                for (int i = start; i <= end; i++) 
                {
                    sum += temp_da_dN[i * crackData.getdomain_types() + domain];
                }
                da_dN[index] = sum / (end - start + 1);
            }
            else if (weighted_moving_average == 2)  // 重み付き移動平均
            {
                for (int i = start; i <= end; i++) 
                {
                    // 中心からの距離に応じた重みを計算
                    int distance = std::abs(i - point);
                    int weight = number_moving_average - distance;

                    sum += temp_da_dN[i * crackData.getdomain_types() + domain] * weight;
                    weight_sum += weight;
                }
                da_dN[index] = sum / weight_sum;
            }
        }
    }

    #if 1
    // 平滑化の種類に応じたデバッグ出力
    std::cout << std::endl << "=== ";
    if (weighted_moving_average == 1)
        std::cout << "Moving Average";
    else if (weighted_moving_average == 2)
        std::cout << "Weighted Moving Average";
    std::cout << " da/dN ===" << std::endl << std::endl;

    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        std::cout << "\n=== Domain Size " << domain + 1 << " ===" << std::endl;
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            std::cout << "Point " << point 
                     << ": Original = " << std::scientific << std::setprecision(6) 
                     << temp_da_dN[index] 
                     << ", Smoothed = " << da_dN[index] << " [mm/cycle]" << std::endl;
        }
    }
    #endif
}

void ParisLawCrack::calcMaxCrackGrowthRate(const CrackDataSetting_3D& crackData)
{
    std::cout << std::endl <<"=== Maximum crack growth rate ===" << std::endl <<std::endl;
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        double max_growth_rate = 0.0;
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            int index = point * crackData.getdomain_types() + domain;
            if (da_dN[index] > max_growth_rate) 
            {
                max_growth_rate = da_dN[index];
            }
        }
        da_dN_max[domain] = max_growth_rate;
        // setMaxCrackGrowthRate(da_dN_max);
        // for debug
        #if 1
        // cout << "Domain Size " << domain + 1 << ": " << da_dN_max[domain] << " [mm/cycle]" << endl;
        cout << "Domian Size " << domain + 1 << ": (da/dN)_max = " << getMaxCrackGrowthRate()[domain] << " [mm/cycle]" << endl;
        #endif
    }
} 


void ParisLawCrack::calcCrackGrowthCycle(const CrackDataSetting_3D& crackData) 
{   
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        Nk[domain] =  delta_a_max / getMaxCrackGrowthRate()[domain];  // Nk = Δa / (da/dN)max
    }

    // 荷重サイクル数に進展前後の平均値を使用する場合 Nk = (Nk + Nk_after) / 2.0
    if (CALC_LOAD_CYCLES == 1) 
    {   
        std::vector<double> Nk_after(crackData.getdomain_types(), 0.0);  // 進展後のサイクル数を格納
        readAfterCrackGrowthCycle(crackData, Nk_after);                  // 進展後のサイクル数を読み取り
        calcAverageCrackGrowthCycle(crackData, Nk_after);                // 平均サイクル数を計算
    }

    // 荷重サイクル数ファイルを出力
    writeCrackGrowthCycle(crackData);

    // for debug
    #if 0
    std::cout << std::endl << "=== Crack Growth Cycle ===" << std::endl << std::endl;
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        std::cout << "Domain Size " << domain + 1 << ": " << getCrackGrowthCycle()[domain] << " [cycle]" << endl;
    #endif
}


void ParisLawCrack::readAfterCrackGrowthCycle(const CrackDataSetting_3D& crackData, std::vector<double>& Nk_after)
{
    std::ifstream file("Nk_after.csv");
    if (!file) {
        std::cerr << "Error: Cannot open file: Nk_after.csv" << std::endl;
        exit(1);
    }

    std::string line;
    // ヘッダー行をスキップ
    std::getline(file, line);

    // データ行を読み込み
    std::getline(file, line);
    std::stringstream ss(line);
    std::string cell;

    // domain列を読み込み
    std::getline(ss, cell, ',');
    // if (cell != "Nk_after") 
    // {
    //     std::cerr << "Error: Invalid data type in after_Nk.csv. Expected 'after_Nk' but got '" << cell << "'" << std::endl;
    //     exit(1);
    // }

    // 各積分領域のNk値を読み込み
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        std::getline(ss, cell, ',');
        Nk_after[domain] = std::stod(cell);
    }

    file.close();

    // for debug
    #if 0
    std::cout << std::endl << "=== After Crack Growth Cycle ===" << std::endl << std::endl;
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        std::cout << "Domain Size " << domain + 1 << ": " << Nk_after[domain] << " [cycle]" << std::endl;
    #endif
}


void ParisLawCrack::calcAverageCrackGrowthCycle(const CrackDataSetting_3D& crackData, std::vector<double>& Nk_after)
{
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        Nk[domain] = (Nk[domain] + Nk_after[domain]) / 2.0;
    }

    // for debug
    #if 0
    std::cout << std::endl << "=== Average Crack Growth Cycle ===" << std::endl << std::endl;
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
        std::cout << "Domain Size " << domain + 1 << ": " << getCrackGrowthCycle()[domain] << " [cycle]" << std::endl;
    #endif
}


void ParisLawCrack::writeCrackGrowthCycle(const CrackDataSetting_3D& crackData) 
{
    std::ofstream file("Nk.csv");
    if (!file) {
        std::cerr << "Error: Cannot open file: Nk.csv" << std::endl;
        exit(1);
    }

    // ヘッダーの書き込み
    file << "domain,";
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        file << "all";
        if (domain > 0) file << " - " << domain;
        if (domain < crackData.getdomain_types() - 1) file << ",";
    }
    file << std::endl;

    // データの書き込み
    file << "Nk_" << (CALC_LOAD_CYCLES == 1 ? "average" : "before") << ",";
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        file << std::scientific << std::setprecision(16) 
             << Nk[domain];
        if (domain < crackData.getdomain_types() - 1) file << ",";
    }
    file << std::endl;

    file.close();
}


void ParisLawCrack::calcCrackGrowthRateComponent(const CrackDataSetting_3D& crackData, const CrackDirection& crackExtensionVector, information *info) 
{   
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            for (int comp = 0; comp < info->DIMENSION; comp++) 
            {
                int growthIndex = domain * crackData.getJintegralPoints() * info->DIMENSION
                                  + point * info->DIMENSION + comp;
                int rateIndex = point * crackData.getdomain_types() + domain;
                da_dN_component[growthIndex] = getCrackGrowthRate()[rateIndex] * crackExtensionVector.getCrackGrowthDirection()[point * info->DIMENSION + comp];
            }
        }
    }

    // for debug
    #if 0
    printf("\n=== Crack Growth Rate Component ===\n\n");
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        for (int point = 0; point < crackData.getJintegralPoints(); point++) 
        {
            for (int comp = 0; comp < info->DIMENSION; comp++) 
            {
                int growthIndex = domain * crackData.getJintegralPoints() * info->DIMENSION
                                    + point * info->DIMENSION + comp;
                std::cout << "Domain Size " << domain + 1
                            << ", Jcp = " << point
                            << ", Component = " << comp
                            << ", Value = " << getda_dN_component()[growthIndex] << " [mm/cycle]" << std::endl;
            }
        }
    }
    #endif
}


void CrackDirection::calcCrackGrowthDirection(const CrackDataSetting_3D& crackData, const IIM_3D& IIM3D, const ParisLawCrack&, information *info) 
{   
    // き裂進展方向の取得
    const auto& tempcrackFrontVector = IIM3D.getJcpCrackFrontVector();
    std::vector<double> newDirection(crackData.getJintegralPoints() * info->DIMENSION, 0.0);

    // for debug
    #if 0
    std::cout << "\nTotal size of crackFrontVector: " << tempcrackFrontVector.size() << std::endl;
    std::cout << "Expected size: " << crackData.getJintegralPoints() * info->DIMENSION * info->DIMENSION << std::endl;

    std::cout << "\n=== Raw Vector Data ===\n";
    for (size_t i = 0; i < tempcrackFrontVector.size(); ++i) 
    {
        std::cout << "Index " << i << ": " << tempcrackFrontVector[i] << std::endl;
    }
    #endif

    // normalベクトル（き裂進展方向）の抽出
    int normal_index = 2;  
    for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {   
            int src_index = Jcp * info->DIMENSION * info->DIMENSION 
                           + normal_index * info->DIMENSION + j;
            int dest_index = Jcp * info->DIMENSION + j;
            newDirection[dest_index] = tempcrackFrontVector[src_index];
        }
    }
    setCrackGrowthDirection(newDirection);

    // for debug
    #if 0
    std::cout << "\n=== Crack Growth Direction ===" << std::endl;
    for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) 
    {
        for (int j = 0; j < info->DIMENSION; j++)
        {
            int index = Jcp * info->DIMENSION + j;
            std::cout << "Jcp = " << Jcp
                    << ", Component = " << j 
                    << ", Value = " << getCrackGrowthDirection()[index] << std::endl;
        }
    }
    #endif
}


void CrackDirection::calcCrackGrowthVector(const CrackDataSetting_3D& crackData, const ParisLawCrack& paris, information *info)
{   
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) // domainループ
    {
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) // J-integral points
        {
            for (int comp = 0; comp < info->DIMENSION; comp++) // ベクトル成分 (x, y, z)
            {
                // インデックス計算
                int growthIndex = domain * crackData.getJintegralPoints() * info->DIMENSION
                                + Jcp * info->DIMENSION + comp;
                // き裂進展ベクトルの計算
                CrackGrowthVector[growthIndex] = paris.getCrackGrowthCycle()[domain] * paris.getda_dN_component()[growthIndex];
            }
        }
    }

    // for debug
    #if 1
    printf("\n=== Crack Growth Extension Vector ===\n\n");
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) // domainループ
    {
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) // J-integral points
        {
            for (int comp = 0; comp < info->DIMENSION; comp++) // ベクトル成分 (x, y, z)
            {
                // インデックス計算
                int growthIndex = domain * crackData.getJintegralPoints() * info->DIMENSION
                                + Jcp * info->DIMENSION + comp;
                // き裂進展ベクトルの計算
                std::cout << "Domain Size " << domain + 1 
                          << ", Jcp = " << Jcp 
                          << ", Component = " << comp 
                          << ", Value = " << getCrackGrowthVector()[growthIndex] << " [mm]" << std::endl;
            }
        }
    }
    #endif
}



// 進展前き裂前縁座標を取得
void CrackFrontManager::getbeforeCrackFrontPoints(const IIM_3D& IIM3D)
{
    before_crack_front_points = IIM3D.getCalculationPointCoords();
    setbeforeCrackFrontPoints(before_crack_front_points);

    // for debug
    #if 0
    std::cout << std::endl << "=== Before Crack Front Points ===" << std::endl;
    for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) 
    {
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {
            int index = Jcp * info->DIMENSION + dim;
            // std::cout << "Jcp: " << Jcp << ", Jcenter in physical coordinate[" << dim << "]: " << before_crack_front_points[index] << std::endl;
            std::cout << "Jcp: " << Jcp << ", Jcenter in physical coordinate[" << dim << "]: " << getbeforeCrackFrontPoints()[index] << std::endl;
        }
    }
    #endif
}

void CrackFrontManager::calcaferCrackFrontPoints(const CrackDirection& crackExtensionVector, const CrackDataSetting_3D& crackData, information *info)
{
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {    
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) 
        {
            for (int dim = 0; dim < info->DIMENSION; dim++) 
            {
                // ドメインを考慮したインデックス計算
                int after_index = domain * (crackData.getJintegralPoints() * info->DIMENSION) + Jcp * info->DIMENSION + dim;
                after_crack_front_points[after_index] = before_crack_front_points[Jcp * info->DIMENSION + dim] + crackExtensionVector.getCrackGrowthVector()[after_index];
            }
        }
    }

    // for debug
    #if 0
    std::cout << std::endl << "=== After Crack Front Points ===" << std::endl;
    for (int domain = 0; domain < crackData.getdomain_types(); domain++) 
    {
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) 
        {
            for (int dim = 0; dim < info->DIMENSION; dim++) 
            {
                // ドメインを考慮したインデックス計算
                int after_index = domain * (crackData.getJintegralPoints() * info->DIMENSION) + Jcp * info->DIMENSION + dim;
                std::cout << "Domain Size " << crackData.getdomain_types() << ", Jcp: " << Jcp << ", Jcenter in physical coordinate[" << dim << "]: " << after_crack_front_points[after_index] << std::endl;
            }
        }
    }
    #endif
}


void FatigueCrackIO::read_select_domain(const std::string& filename, const CrackDataSetting_3D& crackData)
{
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(1);
    }

    // 積分パラメータをスキップ
    int domainWidth, JintegralPoints, model_size_para_3D, domain_types, domain_area_divisions, patchType;
    int num_edge_crack_patch, num_singular_patch, patchID;

    file >> domainWidth >> JintegralPoints >> model_size_para_3D >> domain_types >> domain_area_divisions >>patchType;
    file >> num_edge_crack_patch >> num_singular_patch;

    for (int i = 0; i < num_edge_crack_patch * num_singular_patch; i++)
    {
        file >> patchID;
    }

    // き裂モードをスキップ
    int mode;
    file >> mode;

    // Paris則パラメータのスキップ
    double C_val, m_val;
    file >> C_val >> m_val;

    // 最大き裂進展長さのスキップ
    double max_growth;
    file >> max_growth;

    // 繰り返し荷重サイクル数の計算選択をスキップ
    // int selectNk;
    // file >> selectNk;

    // き裂進展量に対して荷重移動平均を計算するかをスキップ
    int weight_move_ave, num_move_ave;
    file >> weight_move_ave >> num_move_ave;

    // 積分領域の選択を読み込み
    int domainval;
    file >> domainval;

    // 選択した積分領域が有効かチェック
    if (domainval > crackData.getdomain_types()) 
    {
        std::cerr << "Error: Invalid domain size" << std::endl;
        exit(1);
    }
    setSelectDomain(domainval);
    
    # if 0
    std::cout << std::endl <<"Selected domain size: " << domainval << std::endl;
    std::cout << std::endl << "Selected domain size: " << getSelectDomain() << std::endl;
    std::cout << std::endl;
    #endif

    file.close();
}


void FatigueCrackIO::writeCrackFrontData(const CrackFrontManager& crackFrontCoord, information *info)
{
    std::cout << std::endl << "=== Write Crack Front Data ===" << std::endl << std::endl;
    
    // 選択されたドメインサイズの確認
    int selected_domain = getSelectDomain();
    std::cout << "Selected domain size: " << selected_domain << std::endl;

    // ファイルを開く（存在しない場合は新規作成）
    std::ofstream outFile;
    outFile.open("Crack_Front_Points_and_CP.csv", std::ios::out | std::ios::trunc);
    if (!outFile) 
    {
        std::cerr << "Error: Cannot create file: Crack_Front_Points_and_CP.csv" << std::endl;
        return;
    }

    // ヘッダーの書き込み
    outFile << "Jcp,x_before, y_before, z_before, x_after, y_after, z_after" << std::endl;

    // 進展前と進展後の座標データを取得
    const auto& before_points = crackFrontCoord.getbeforeCrackFrontPoints();
    const auto& after_points = crackFrontCoord.getafterCrackFrontPoints();
    cout << "before_points size: " << before_points.size() << endl;
    cout << "after_points size: " << after_points.size() << endl;

    // 各計算点でのき裂前縁座標を書き込み
    int JintegralPoints = before_points.size() / info->DIMENSION; // き裂前縁点数
    cout << "JintegralPoints: " << JintegralPoints << endl;
    for (int Jcp = 0; Jcp < JintegralPoints; Jcp++)
    {
        // Jcp番号
        outFile << Jcp;
        outFile << std::fixed << std::setprecision(16);

        // 進展前の座標 (step0)
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {
            outFile << "," << before_points[Jcp * info->DIMENSION + dim];
        }

        // 進展後の座標 (step1) - 選択されたドメインサイズの結果
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {   
            int after_index = selected_domain * (JintegralPoints * info->DIMENSION) + Jcp * info->DIMENSION + dim;
            // int after_index = Jcp * info->DIMENSION + dim;
            outFile << "," << after_points[after_index];
        }

        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Successfully wrote crack front points to Crack_Front_Points.csv" << std::endl;
}


void FatigueCrackIO::writeCrackFrontCPData(const LeastSquaresApproximation& calc_cp_LSA, information *info)
{
    std::cout << std::endl << "=== Write Control Point Data to CSV ===" << std::endl << std::endl;
    
    const std::string filename = "Crack_Front_Points_and_CP.csv";
    
    // 既存ファイルの内容を読み込む
    std::vector<std::string> existing_lines;
    std::ifstream inFile(filename);
    if (!inFile) 
    {
        std::cerr << "Error: Cannot open file for reading: " << filename << std::endl;
        exit(1);
    }
    
    std::string line;
    while (std::getline(inFile, line)) 
    {
        existing_lines.push_back(line);
    }
    inFile.close();

    // 新しいファイルを書き出し
    std::ofstream outFile(filename);
    if (!outFile) 
    {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        exit(1);
    }

    // 既存のJcpデータを書き出し
    for (const auto& line : existing_lines) 
    {
        outFile << line << std::endl;
    }
    outFile << std::endl;

    

    // 制御点データの取得
    const auto& cp_list = calc_cp_LSA.getcp_list();
    const auto& cp_coordinates = calc_cp_LSA.getControlPointCoordinate();

    // 制御点データの書き出し
    outFile << "CP_ID,Bx_before, By_before,Bz_before,Bx_after,By_after,Bz_after" << std::endl;

    std::cout << "Total control points written: " << cp_list.size() << std::endl;

    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        // 進展前の制御点座標
        outFile << std::scientific << std::setprecision(16);
        outFile << cp_list[i] << ",";
        for (int dir = 0; dir < info->DIMENSION; dir++) 
        {
            outFile << info->Node_Coordinate[cp_list[i] * (info->DIMENSION + 1) + dir];
            outFile << ",";
        }

        // 進展後の制御点座標
        for (int dir = 0; dir < info->DIMENSION; dir++) 
        {
            outFile << cp_coordinates[i * info->DIMENSION + dir];
            if (dir < info->DIMENSION - 1) 
            {
                outFile << ",";
            }
        }
        outFile << std::endl;
    }

    outFile.close();
    
    std::cout << "Successfully wrote control point data to CSV file" << std::endl;
    std::cout << "Total control points written: " << cp_list.size() << std::endl;
    for (size_t i = 0; i < cp_list.size(); i++) {
        std::cout << "[CSV DEBUG] cp_list[" << i << "] = " << cp_list[i] << std::endl;
    }
    
}

// パッチごとの処理 or 全体での処理を選択して 前処理
void LeastSquaresApproximation::preLeastSquaresApproximation(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information *info, bool use_total)
{
    if (use_total) {
        std::cout << "All patches - パッチ全体で最小二乗法を実行" << std::endl;
        preLSA_Global(crackData, IIM3D, info);
    } else {
        std::cout << "Patchwise - パッチごとに最小二乗法を実行" << std::endl;
        preLSA_Patchwise(crackData, IIM3D, info);
    }
}
// パッチごとの処理 or 全体での処理を選択して 実行
void LeastSquaresApproximation::DoLeastSquaresApproximation(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, IIM_3D& IIM3D, information *info, int selected_domain, bool use_total)
{
    if (use_total) 
    {
        DoLeastSquaresApproximation_Global(crackData, crackFrontCoord, IIM3D, info, selected_domain);
    }
    else
    {
        DoLeastSquaresApproximation_Patchwise(crackData, crackFrontCoord, IIM3D, info, selected_domain);
    }
}

// 全体での処理
// 前処理
void LeastSquaresApproximation::preLSA_Global(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information* info)
{
     // 制御点リスト、インデックスマップ、J積分点情報リストを初期化
    cp_list.clear();
    Jcp_info_list.clear();
    cp_global_index_map.clear();

    // ローカルパッチ番号の初期値と、制御点インクリメントアクセス用のオフセットを取得
    int local_init_patch = info->Total_Patch_to_mesh[Total_mesh -1];
    int offset = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;

     // 最小二乗法を行うパッチリストの取得
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {   
        int index = edge_patch * crackData.getnum_singular_patch();
        LSA_patchList[edge_patch] = crackData.getpatchList()[index] + local_init_patch;
    }
    #if 1
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {
        std::cout << "LSA Patch List: " << LSA_patchList[edge_patch] << std::endl;
    }
    #endif

    // パラメータ空間の方向（ξ, η, ζ）と自然座標の一時保存ベクトルを定義
    std::vector<int> direction = {0, 1, 2};
    std::vector<double> temp_Natural_Coord(info->DIMENSION, 0.0);

    //コントロールポイント数とJ積分点数
    int total_cp = 0;
    int total_Jcp = 0;

    // き裂前縁を分割するパッチ数でループ
    for(int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++)
    {
        // パッチIDを取得(全体⇒メッシュ内ローカル番号に変換)
        int patchID = crackData.getpatchList()[edge_patch * crackData.getnum_singular_patch()] + local_init_patch; 

        cp_size[edge_patch] = 0;

        // 各コントロールポイントについて調べる
        for (int i = 0; i < info->No_Control_point_in_patch[patchID]; i++) 
        {
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[patchID]+ i]; //全制御点に対する番号

            // ξ=0かつζ=0を満たすコントロールポイントかチェック
            int inc_xi  = info->INC[patchID * offset + cp_num * info->DIMENSION + direction[0]];
            int inc_zeta = info->INC[patchID * offset + cp_num * info->DIMENSION + direction[2]];
            bool is_crack_front_cp = (inc_xi == 0 && inc_zeta == 0);

            //き裂前縁のコントロールポイントかつ重複していないなら、cp_list(制御点IDのリスト(順番付き))とマップ(制御点ID → 列番号のmap)に追加
            if(is_crack_front_cp)
            {
                cp_list.push_back(cp_num);
                cp_size[edge_patch]++;

                if(cp_global_index_map.count(cp_num) == 0)
                {
                    cp_global_index_map[cp_num] = total_cp++;
                }
            }

        }

        // J積分点について調べる
        for(int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++)
        {
            // η方向の自然座標を取得
            temp_Natural_Coord[1] = IIM3D.getcrackEdgeJcpNatural()[Jcp];
            bool element_found = false;
            int elementID = -1;
            double eta_total =0.0;

            // パッチに属する要素の中で、このJ積分点が含まれる要素を探す
            for(int e = 0; e < info->Total_Element_to_mesh[Total_mesh]; e++)
            {
                //この要素が調べているパッチにあるか確認
                if(info->Element_patch[e] != patchID) continue;

                // [0, 1]空間で与えられた座標点が特定の要素内に存在するかどうかを判定する関数
                if (IIM3D.isCoordinateInElement(temp_Natural_Coord, patchID, e, info)) {
                    elementID = e;
                    element_found = true;
                    eta_total = (Jcp + 1.0) / ( (crackData.getJintegralPoints() + 1.0));
                    cout << "eta_total = " << eta_total << endl;
                    break;
                }
            }

            //要素が見つかれば最後の処理
            if (!element_found) continue;

            // このJ点が確かにこのパッチのこの要素に属していれば登録
            if (elementID == IIM3D.getcrackEdgeJcpElementID()[Jcp]) {
                JcpInfo info_item = {
                    IIM3D.getcrackEdgeJcpNatural()[Jcp],
                    patchID,
                    elementID,
                    eta_total
                };
                Jcp_info_list.push_back(info_item);
                total_Jcp++;
                Jcp_size[edge_patch]++;
            }
        }
    }

    // 全体のJ積分点数が制御点数より少ないと、最小二乗法が解けないためエラー
    if (total_Jcp < total_cp) {
        std::cerr << "Error: J-integral points (" << total_Jcp << ") are fewer than control points (" << total_cp << ") in global mode." << std::endl;
        exit(1);
    }

    // 各edgeパッチがcp_listのどこで始まるか
    patch_cp_start_index.clear();
    int running_index = 0;
    for (size_t edge_patch_num = 0; edge_patch_num < cp_size.size(); edge_patch_num++) 
    {
        patch_cp_start_index[edge_patch_num] = running_index;
        running_index += cp_size[edge_patch_num];
    }

    // for debug
    #if 1
    std::cout << "\n=== cp_list -> cp_global_index_map 対応確認 ===\n";
for (size_t i = 0; i < cp_list.size(); ++i) {
    int cp_id = cp_list[i];
    auto it = cp_global_index_map.find(cp_id);
    if (it != cp_global_index_map.end()) {
        std::cout << "cp_list[" << i << "] = " << cp_id 
                  << ", global_index = " << it->second << std::endl;
    } else {
        std::cout << "cp_list[" << i << "] = " << cp_id 
                  << " → [ERROR] not found in cp_global_index_map!" << std::endl;
    }
}

    #endif
}
// 実行
void LeastSquaresApproximation::DoLeastSquaresApproximation_Global(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, IIM_3D& IIM_3D, information *info, int selected_domain)
{
    std::vector<double> temp_Natural_Coord(info->DIMENSION, 0.0);
    std::vector<int> direction = {0 ,1, 2};
    static int max_support_1D = 2 * MAX_ORDER + 2;

    int total_cp = cp_global_index_map.size();
    int total_Jcp = Jcp_info_list.size();

    cout << "total_Jcp = " << total_Jcp << endl;
    cout << "total_cp = " << total_cp << endl;

    #if 1
    std::cout << "cp_list.size() " << total_cp << "Jcp_info_list.size() " << total_Jcp << endl;
    #endif

    //Nマトリクス, 
    std::vector<double> N_matrix(total_Jcp * total_cp, 0.0);
    std::vector<double> NTN_inv_NT_matrix(total_cp * total_Jcp, 0.0);
    int cp_list_offset = 0;
    std::vector<double> knot_vector = generate_Open_Knot_Vector(total_cp, info->Order[4]);

    for(int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++)
    {
        int patchID =LSA_patchList[edge_patch];
        cout<< "patchID = " << patchID << endl;
        int offset_row = 0;
        for(int ep = 0; ep < edge_patch; ep++) offset_row += getJcp_size(ep);
        cout << "cp_list_offset = "<< cp_list_offset << endl;
        cout << "offset_row = "<< offset_row << endl;

        for(int j = 0; j < getJcp_size(edge_patch); j++)
        {
            // jcpにJcp_info_listを入れてく
            //const auto& jcp = Jcp_info_list[j];
            // j番目のJ積分点のη座標
            temp_Natural_Coord[1] = IIM_3D.getcrackEdgeJcpNatural()[offset_row + j];
            int elementID = IIM_3D.getcrackEdgeJcpElementID()[offset_row + j];
            int offset_ENC = info->ENC[elementID * info->DIMENSION + direction[1]];

            //std::cout << "Jcp = " << j << endl;
            //std::cout << "jcp.eta = " << temp_Natural_Coord[1] << endl;
            //std::cout << "jcp.eta_total = "<< jcp.eta_total << endl;


            //基底関数マトリクス、微分マトリクス
            std::vector<double> N_row(MAX_ORDER * MAX_CP, 0.0);
            std::vector<double> dN_row(MAX_CP, 0.0);

            int range_eta = info->Order[patchID * info->DIMENSION + direction[1]] + 1;
            shape_function_1D(temp_Natural_Coord[1], direction[1], elementID, N_row.data(), dN_row.data(), 0, info);
            //calcShapeFunction1D_Global(jcp.eta_total, knot_vector, info->Order[4], N_row);



            for(int i = 0; i < range_eta; i++) 
            {
                int cp_id =cp_list[cp_list_offset + j];
                int k = direction[1] * MAX_ORDER + info->Order[patchID * info->DIMENSION + direction[1]];
                double val = N_row[k * max_support_1D + i];
                if (val == 0.0) continue; // 非ゼロの基底関数だけ処理
                int col = cp_global_index_map[cp_id];
                int row = offset_row + j;
                //N_matrix[row * total_cp + col] += val;
                N_matrix[j * total_cp + ((getcp_size(edge_patch) - 1) * edge_patch + i)] += val;
            }

            #if 0 
        for(int i = 0; i < order_eta + 1; i++)
        {

            // J積分点がいる要素のη方向i番目のコントロールポイントの全体でのindex
            int global_cp_id = info->Controlpoint_of_Element[jcp.elementID * MAX_NO_CP_ON_ELEMENT + i];

            //for(int k = 0; k < (order_eta + 1) * (order_eta + 1); k++)std::cout << info->Controlpoint_of_Element[jcp.elementID * MAX_NO_CP_ON_ELEMENT + k] << endl;

            int col =  cp_global_index_map[global_cp_id];
            int row = j;

            //cout << row << "行 " << col << "列 " << endl;

            // j番目のJ積分点のη方向i番目のコントロールポイントに対応する形状関数を配列に入れる
            N_matrix[row * total_cp + col] += N_row[col];
            //std::cout << "cp_id = " << global_cp_id << ", col = " << col << ", N = " << Shape[shape_index] << std::endl;

             // for debug
            
            std::cout << "=== Shape Function Debug ===" << std::endl;
                for(int i = 0; i < total_Jcp; i++)
                {
                    for(int j = 0; j < total_cp; j++)
                    {
                        std::cout << N_matrix[i * total_cp + j];
                    }
                    std::cout << endl;
                }
        } 
            #endif

            #if 1 
                 // for debug
                std::cout << "=== Shape Function Debug ===" << std::endl;
                    for(int i = 0; i < total_Jcp; i++)
                    {
                        for(int j = 0; j < total_cp; j++)
                        {
                            std::cout << N_matrix[i * total_cp + j] << " ";
                        }
                        std::cout << endl;
                    }
            #endif
        }
                    cp_list_offset += cp_size[edge_patch];
    }

// (N^T N)^(-1) N^T の計算
calcNTN_inverse_NT_Global(N_matrix, NTN_inv_NT_matrix);

// ((N^T N)^(-1) * N^T) * Jcp_coordinate の計算
calcControlPointCoordinate_Global(crackFrontCoord, NTN_inv_NT_matrix, info, selected_domain);

//補正
correctControlPointCoordinate_Global(info);

// ここで複製
copyGlobalCPResultToAll();

// 重複点は全体処理では生じないので不要
// deleteDuplicateCP(crackData, info);

}

std::vector<double> LeastSquaresApproximation::generate_Open_Knot_Vector(int num_cp, int degree)
{
    int num_knots = num_cp + degree + 1; // ノット数
    std::vector<double> knot_vector(num_knots, 0.0);

    int num_internal_knots = num_knots - 2 * (degree + 1); // 中間ノット数

    for(int i =  0; i < num_internal_knots ; i++)   
    {
        knot_vector[degree + 1 + i] = (double) (i + 1) / (num_internal_knots + 1);
    }


    // 明示的に最初と最後のノットを固定（p+1 回繰り返す）
    for (int i = 0; i < degree + 1; i++)
    {
        knot_vector[i] = 0.0;
        knot_vector[num_knots - 1 - i] = 1.0;
    }

    for(int j = 0; j < num_knots; j++)
    {
        std::cout << knot_vector[j] << endl;
    }
    return knot_vector;

}

void LeastSquaresApproximation::calcShapeFunction1D_Global(double eta, const std::vector<double>& knots, int degree, std::vector<double>& N)
{
    int n = knots.size() - degree - 1;
    int i, p;
    cout << knots.size() << endl;
    for(i = 0; i < knots.size() - 1; i++)
    {
        if(knots[i] == knots[i + 1])
        {
            N[i * MAX_ORDER + 0] = 0.0;
        }
        else if(knots[i] != knots[i + 1] && knots[i] <= eta && eta < knots[i + 1])
        {
            N[i * MAX_ORDER + 0] = 1.0;
        }
        else
        {
            N[i * MAX_ORDER + 0] = 0.0;
        }
        //cout << "N[0] = "<< N[i * MAX_ORDER + 0] << endl;
    }

    double left_term, right_term;
    for(p = 1; p <= degree; p++)
    {
        for(i = 0; i < n; i++)
        {
            left_term = 0.0;
            right_term = 0.0;

            if((eta - knots[i]) * N[i * MAX_ORDER + p - 1] == 0 && knots[i + p] - knots[i] == 0)
                left_term = 0.0;
            else
                left_term = (eta - knots[i]) / (knots[i + p] - knots[i]) * N[i * MAX_ORDER + p - 1];
            
            if((knots[i + p + 1] - eta) * N[(i + 1) * MAX_ORDER + p - 1] == 0 && knots[i + p + 1] - knots[i + 1] == 0)
                right_term = 0.0;
            else
                right_term = (knots[i + p + 1] - eta) / (knots[i + p + 1] - knots[i + 1]) * N[(i + 1) * MAX_ORDER + p - 1];
        
                N[i * MAX_ORDER + p] = left_term + right_term;
            } 
    }
    #if 0
    for(i = 0; i < knots.size(); i++) cout<<"N [1] = "<< N[i * MAX_ORDER + 1] << endl;
    for(i = 0; i < knots.size(); i++) cout<<"N [2] = "<< N[i * MAX_ORDER + 2] << endl;
    #endif
}


void LeastSquaresApproximation::calcNTN_inverse_NT_Global(std::vector<double>& N_matrix, std::vector<double>& NTN_inverse_NT_matrix)
{
    int total_cp = cp_global_index_map.size();
    int total_Jcp = Jcp_info_list.size();

    std::vector<double> NTN_matrix(total_cp * total_cp, 0.0);
    std::vector<double> NTN_inverse_matrix(total_cp * total_cp, 0.0);
    
    // NTNの計算  N_matrixは J*C, N^T_matrixは C*J, NTN_matrixは C*C
    for(int i = 0; i < total_cp; i++)
    {
        for(int j = 0; j < total_cp; j ++)
        {
            for(int k = 0; k < total_Jcp; k++)
            {
                // NTN_matrix[i][j] = Σ_k N[k][i] * N[k][j]
                NTN_matrix[i * total_cp + j] += N_matrix[k * total_cp + i] * N_matrix[k * total_cp + j];
            }
        }
    }

    #if 0
    std::cout << "NTN_matrix" << endl;
    for(int i = 0; i < total_cp; i++)
    {
        for(int j = 0; j < total_cp; j ++)
        {
                // NTN_matrix[i][j] = Σ_k N[k][i] * N[k][j]
                std::cout << NTN_matrix[i * total_cp + j] << " ";
        }
        std::cout << endl;
    }
    #endif
    
    //NTN_inverse_NT_matrixまで計算
    calcInverseMatrix_byEigen_Global(NTN_matrix, NTN_inverse_matrix, total_cp, total_Jcp);

    // (N^T)^(-1) * N^T 行列の計算
    for (int i = 0; i < total_cp; i++) 
    {
        for (int j = 0; j < total_Jcp; j++) 
        {
            for (int k = 0; k < total_cp; k++) 
            {
                NTN_inverse_NT_matrix[i * total_Jcp + j] += NTN_inverse_matrix[i * total_cp + k] * N_matrix[j * total_cp + k];
            }
        }
    }

    #if 0
    std::cout << "NTN_inverse_NT_matrix" << endl;
    for(int i = 0; i < total_cp; i++)
    {
        for(int j = 0; j < total_Jcp; j ++)
        {
                // NTN_matrix[i][j] = Σ_k N[k][i] * N[k][j]
                std::cout << NTN_inverse_NT_matrix[i * total_Jcp + j] << " ";
        }
        std::cout << endl;
    }
    #endif
}

void LeastSquaresApproximation::calcInverseMatrix_byEigen_Global(std::vector<double>& NTN_matrix, std::vector<double>& NTN_inverse_matrix, int total_cp, int total_Jcp)
{
    using namespace Eigen;

  // NTN行列をEigen形式に変換
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> NTN = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
(NTN_matrix.data(), total_cp, total_cp);

// LU分解による逆行列計算
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Identity = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(total_cp, total_cp);
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> NTN_inverse = NTN.lu().solve(Identity);

// 結果を標準のvectorに戻す
Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(NTN_inverse_matrix.data(), total_cp, total_cp) = NTN_inverse;

// for debug
#if 0
std::cout << "\n=== NTN Inverse Matrix ===\n";
for (int i = 0; i < cp_size[edge_patch]; i++)
{
for (int j = 0; j < cp_size[edge_patch]; j++)
{
std::cout << NTN_inverse_matrix[i * cp_size[edge_patch] + j] << " ";
}
std::cout << std::endl;
}
#endif

// NTN * NTN^(-1)が単位行列になることを確認
#if 0
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result = NTN * NTN_inverse;
std::cout << "\n=== Verification (NTN * NTN^(-1)) ===\n" << result << std::endl;
#endif
}

//後でEigen使う形式に変える
// 重複しない各制御点の座標を計算
void LeastSquaresApproximation::calcControlPointCoordinate_Global(const CrackFrontManager& crackFrontCoord, vector<double>& NTN_inverse_NT_matrix, information *info, int selected_domain)
{
    int total_cp = cp_global_index_map.size();
    int total_Jcp = Jcp_info_list.size();
    int dim = info->DIMENSION;
    
    // 各制御点（行）について座標を計算
    for (int cp = 0; cp < total_cp; cp++) {
        for (int d = 0; d < dim - 1; d++) {
            double temp_coord_value = 0.0;

            // 最小二乗近似で計算したマトリクスとき裂前縁座標を掛け算して、コントロールポイント座標を計算
            for (int j = 0; j < total_Jcp; ++j) {
                int j_index = selected_domain * (total_Jcp * dim) + j * dim + d;
                double matrix = NTN_inverse_NT_matrix[cp * total_Jcp + j];
                double cp_coord = crackFrontCoord.getafterCrackFrontPoints()[j_index];
                temp_coord_value += matrix * cp_coord;
            }

            // 結果を格納（制御点cp_idx、方向d）
            crack_front_cp_coordinate[cp * dim + d] = temp_coord_value;

            #if 0
            std::cout << temp_coord_value << endl;
            #endif
        }
        cout << "crack_front_cp_coordinate[cp * dim + 0] = " << crack_front_cp_coordinate[cp * dim + 0] << endl;
        cout << "crack_front_cp_coordinate[cp * dim + 1] = " << crack_front_cp_coordinate[cp * dim + 1] << endl;
        cout << "crack_front_cp_coordinate[cp * dim + 2] = " << crack_front_cp_coordinate[cp * dim + 2] << endl;

    }
}

void LeastSquaresApproximation::correctControlPointCoordinate_Global(information *info)
{
    const int dim = info->DIMENSION;
    
    // 補正を行う制御点のリストを取得
    std::vector<int> correction_patch_set, edge_parameter_set, fixed_direction_set;
    ReadCPSettingLSA("J_int.txt", correction_patch_set, edge_parameter_set, fixed_direction_set);
    
    for(size_t i = 0; i < correction_patch_set.size(); i++)
    {
        int patch = correction_patch_set[i];
        int is_end = edge_parameter_set[i];
        int dir = fixed_direction_set[i];

        cout << "patch = " << patch << "is_end = "<< is_end << "dir = "<< dir << endl;

         // 制御点が存在しないパッチはスキップ
        if (cp_size[patch] == 0) 
        continue;

        //指定されたパッチ内の始点 0 or 終点の制御点IDを取得 cp_size[patch] -1
        int local_id = (is_end == 0) ? 0 : cp_size[patch] - 1;

        // preLSAで作った配列を使ってmapのアドレスを作成 該当パッチの始まる番号 + (始点 or 終点)
        size_t cp_index_in_list = patch_cp_start_index.at(patch) + local_id;

        // 範囲外アクセス防止
        if (cp_index_in_list >= cp_list.size()) {
            std::cerr << "Warning: cp_index_in_list " << cp_index_in_list << " out of range." << std::endl;
            continue;
        }

        int cp_id = cp_list[cp_index_in_list];
        

        //mapに該当するコントロールポイントがあるか確認,無ければ警告
        if(cp_global_index_map.count(cp_id) == 0)
        {
            std::cerr << "Warning: cp_id " << cp_id << " not found in cp_global_index_map" << std::endl;
            continue;
        }

        int col = cp_global_index_map[cp_id];
        cout << "col = " << col << endl;

        //補正(指定の方向に沿って進展した)コントロールポイントを代入
        cout << "crack_front_cp_coordinate[col * dim + dir] = " << crack_front_cp_coordinate[col * dim + dir] << endl;
        cout << "info->Node_Coordinate[cp_id  * (dim + 1) + dir] = " << info->Node_Coordinate[cp_id  * (dim + 1) + dir] << endl;
        crack_front_cp_coordinate[col * dim + dir] = info->Node_Coordinate[cp_id  * (dim + 1) + dir];
    }
}

void LeastSquaresApproximation::copyGlobalCPResultToAll()
{
    const int dim = 3;  // = DIMENSION
    std::vector<double> full_coordinates(cp_list.size() * dim, 0.0);
    std::vector<int> new_cp_list;

    for (size_t i = 0; i < cp_list.size(); ++i) 
    {
        int cp_id = cp_list[i];
        // この制御点IDが初めて出現する場合のみ処理
        if (std::find(new_cp_list.begin(), new_cp_list.end(), cp_id) == new_cp_list.end()) 
        {
            new_cp_list.push_back(cp_id);

            if (cp_global_index_map.count(cp_id) == 0) {
                std::cerr << "[ERROR] cp_id " << cp_id << " not found in cp_global_index_map!" << std::endl;
                continue;
            }
    
            int global_idx = cp_global_index_map[cp_id];
            std::cout << "[copy] cp_id = " << cp_id << ", global_idx = " << global_idx << std::endl;
            for (int d = 0; d < dim - 1; d++) {
                full_coordinates[global_idx * dim + d] = crack_front_cp_coordinate[global_idx * dim + d];
                cout << "full_coordinates = " << full_coordinates[global_idx * dim + d] << endl;
            }
        }

        
    }

    cp_list = std::move(new_cp_list);
    crack_front_cp_coordinate = std::move(full_coordinates);

    //For debug
    #if 1
    for(int cp = 0; cp < cp_list.size(); cp++)
    {
    cout << "cp_list" << cp << endl;
    cout << "crack_front_cp_coordinate[cp * dim + 0] = " << crack_front_cp_coordinate[cp * dim + 0] << endl;
    cout << "crack_front_cp_coordinate[cp * dim + 1] = " << crack_front_cp_coordinate[cp * dim + 1] << endl;
    cout << "crack_front_cp_coordinate[cp * dim + 2] = " << crack_front_cp_coordinate[cp * dim + 2] << endl;
    }
    #endif

}

#if 0
void LeastSquaresApproximation::preLeastSquaresApproximation(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information *info)
{
    int local_init_patch = info->Total_Patch_to_mesh[Total_mesh - 1];              // ローカルパッチ番号の初期値
    int offset = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;  // 前パッチまでの自由度数
    std::vector<int> direction = {0, 1, 2};                                        // パラメータ空間における方向(ξ, η, ζ)
    std::vector<double> temp_Natural_Coord(info->DIMENSION, 0.0);                  // 自然座標の一時保存用
    
    // 最小二乗法を行うパッチリストの取得
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {   
        int index = edge_patch * crackData.getnum_singular_patch();
        LSA_patchList[edge_patch] = crackData.getpatchList()[index] + local_init_patch;
    }
    #if 0
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {
        std::cout << "LSA Patch List: " << LSA_patchList[edge_patch] << std::endl;
    }
    #endif


    // 各パッチの最小二乗法を行う制御点のリストを取得
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {   
        int patchID = LSA_patchList[edge_patch];
        for (int i = 0; i < info->No_Control_point_in_patch[patchID]; i++) 
        {
            //  ξ, ζ両方向で0である制御点（き裂前縁制御点）のみを処理
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[patchID] + i];
            if (info->INC[patchID * offset + cp_num * info->DIMENSION + direction[0]] == 0 && 
                info->INC[patchID * offset + cp_num * info->DIMENSION + direction[2]] == 0) 
            {
                cp_list.push_back(cp_num);
                cp_size[edge_patch]++;  // パッチごとの制御点数をカウント
            }
        }
        
        // J 計算点が属するパッチの判別と各パッチのLSAを行う制御点数をカウント
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++ )
        {   
            // 要素の特定に使用
            bool element_found = false;
            int elementID = -1;
            temp_Natural_Coord[1] = IIM3D.getcrackEdgeJcpNatural()[Jcp]; // η方向の自然座標を取得
            
            // J 計算点が含まれる要素を探索
            for (int e = 0; e < info->Total_Element_to_mesh[Total_mesh]; e++) 
            {   
                // パッチ内の要素を探索
                if (info->Element_patch[e] != patchID) continue;
                
                // J計算点を含む要素の探索
                if (IIM3D.isCoordinateInElement(temp_Natural_Coord, patchID, e, info))
                {
                    elementID = e;
                    element_found = true;      
                    break;
                }
            }
            // J 計算点を持つ要素が見つからない場合はエラー
            if (!element_found) 
            {
                std::cerr << "No element found for J-integral point " << Jcp << " in patch " << patchID << std::endl;
                exit(1);
            }
            // J 計算点が属するパッチの判別
            if (elementID == IIM3D.getcrackEdgeJcpElementID()[Jcp])
            {   
                Jcp_size[edge_patch]++;
            }
        }
        // 各き裂前縁パッチのき裂前縁離散点が制御点より少ない場合は最小二乗法計算を行えない為エラー
        if (Jcp_size[edge_patch] < cp_size[edge_patch]) 
        {
            std::cerr << "Error: crack front points are less than control points in patch " << patchID << std::endl;
            exit(1);
        }
    }
    // for debug
    #if 0
    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "Control Point: " << cp_list[i] << std::endl;
    }
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++)
    {
        std::cout << "Patch " << LSA_patchList[edge_patch] << " has " << cp_size[edge_patch] << " control points" << std::endl;
    }
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++)
    {
        std::cout << "Patch " << LSA_patchList[edge_patch] << " has " << getJcp_size(edge_patch) << " J-integral points" << std::endl;
    }
    #endif
}
#endif

void LeastSquaresApproximation::preLSA_Patchwise(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information* info) 
{
    int local_init_patch = info->Total_Patch_to_mesh[Total_mesh - 1];
    int offset = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;
    std::vector<int> direction = {0, 1, 2};
    std::vector<double> temp_Natural_Coord(info->DIMENSION, 0.0);
    

    //最小二乗法を行うパッチリストの取得
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) {
        int index = edge_patch * crackData.getnum_singular_patch();
        LSA_patchList[edge_patch] = crackData.getpatchList()[index] + local_init_patch;
    }

    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) {
        int patchID = LSA_patchList[edge_patch];
        int local_idx = 0;

        int cp_start_index = cp_list.size(); // 現在の位置
        patch_cp_start_index[edge_patch] = cp_start_index; //各edge_patchが始まる位置

        std::cout << cp_start_index << std::endl;


        for (int i = 0; i < info->No_Control_point_in_patch[patchID]; i++) 
        {
            // ξ, ζ 両方向で0である制御点 (き裂前縁制御点) のみ処理
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[patchID] + i];
            if (info->INC[patchID * offset + cp_num * info->DIMENSION + direction[0]] == 0 &&
                info->INC[patchID * offset + cp_num * info->DIMENSION + direction[2]] == 0) 
            {
                cp_local_index_map[cp_num] = local_idx++; // 制御点ID → パッチ内のローカル列番号 
                cp_list.push_back(cp_num); //cp_listに追加
                cp_global_index_map[cp_num] = cp_list.size() - 1; // グローバルインデックスを登録
                cp_size[edge_patch]++; // パッチごとの制御点数をカウント
                std::cout << "[DEBUG][edge_patch " << edge_patch << "] added CP: " << cp_num << std::endl;

            }
        
        }

        std::cout << "1cp_list.size() = " << cp_list.size() << std::endl;


        // J 計算点が属するパッチの判別と各パッチでLSAを行う制御点数のカウント
        for (int Jcp = 0; Jcp < crackData.getJintegralPoints(); Jcp++) {
            bool element_found = false;
            int elementID = -1;
            double eta_total;
            temp_Natural_Coord[1] = IIM3D.getcrackEdgeJcpNatural()[Jcp]; // Jcpの [0,1] 座標

            for (int e = 0; e < info->Total_Element_to_mesh[Total_mesh]; e++) 
            {
                if (info->Element_patch[e] != patchID) continue; // パッチが一致したら続行
                // Jcpが要素に入ってたらtrueに
                if (IIM3D.isCoordinateInElement(temp_Natural_Coord, patchID, e, info)) {
                    elementID = e;
                    element_found = true;
                    eta_total = (temp_Natural_Coord[1] + (double)edge_patch) / (double)crackData.getnum_edge_crack_patch();
                    break;
                }
            }

            // 見つからなかったらエラー
            if (!element_found) {
                std::cerr << "No element found for J-integral point " << Jcp << " in patch " << patchID << std::endl;
                exit(1);
            }

            // edge_patchごとの J 積分点の数を記録
            if (elementID == IIM3D.getcrackEdgeJcpElementID()[Jcp]) 
            {
                Jcp_size[edge_patch]++;
            }
        }

        // J 積分点の数 が き裂前縁離散点 の数を下回ったらエラー
        if (Jcp_size[edge_patch] < cp_size[edge_patch]) 
        {
            std::cerr << "Error: crack front points are less than control points in patch " << patchID << std::endl;
            exit(1);
        }
    }
}

void LeastSquaresApproximation::DoLeastSquaresApproximation_Patchwise(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, IIM_3D& IIM3D, information *info, int selected_domain)
{   
    static int max_support_1D = 2 * MAX_ORDER + 2;                                       // 1次元の最大サポート数
    std::vector<double> temp_Natural_Coord(info->DIMENSION, 0.0);                        // 自然座標の一時保存用
    std::vector<int> direction = {0, 1, 2};                                              // パラメータ空間における方向(ξ, η, ζ)
    int Jcp_flag = 0;                                                                    // J 計算点のフラグ

    // き裂前縁パッチごとの処理
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {   
        std::vector<double> N_matrix(getJcp_size(edge_patch) * getcp_size(edge_patch), 0.0);         // N 行列 (cp * Jcp)
        std::vector<double> NTN_inverse_NT_matrix(cp_size[edge_patch] * Jcp_size[edge_patch], 0.0);  // (N^T)^(-1) * N^T 行列 (Jcp * cp) 
        int patchID = LSA_patchList[edge_patch];

        // 
        for (int Jcp = 0; Jcp < getJcp_size(edge_patch); Jcp++ )
        {   
            temp_Natural_Coord[1] = IIM3D.getcrackEdgeJcpNatural()[Jcp_flag];
            int elementID = IIM3D.getcrackEdgeJcpElementID()[Jcp_flag];   
            int offset_ENC = info->ENC[elementID * info->DIMENSION + direction[1]];
            std::cout << "patch ID: " << patchID << ", Jcp: " << Jcp << ", Element ID: " << elementID << ", Jcp Natural Coordinate: " << temp_Natural_Coord[1] <<std::endl;

            // 形状関数の計算
            std::vector<double> Shape(MAX_ORDER * MAX_CP, 0.0);
            std::vector<double> dShape(MAX_CP, 0.0);    
            shape_function_1D(temp_Natural_Coord[1], direction[1], elementID, Shape.data(), dShape.data(), 0, info);  // η方向の形状関数の計算
                

            // N マトリクスに格納
            for (int i = 0; i < info->Order[patchID * info->DIMENSION + direction[1]] + 1; i++)
            {
               N_matrix[Jcp * getcp_size(edge_patch) + offset_ENC + i] += Shape[(direction[1] * MAX_ORDER + info->Order[patchID * info->DIMENSION + direction[1]]) * max_support_1D + i];
                /*
                int local_cp_num = offset_ENC + i;
                int global_cp_id = info->Patch_Control_point[info->Total_Control_Point_to_patch[patchID] + local_cp_num];
                int col = cp_local_index_map[global_cp_id];

                double Nval = Shape[(direction[1] * MAX_ORDER + info->Order[patchID * info->DIMENSION + direction[1]]) * max_support_1D + i];
                N_matrix[Jcp * getcp_size(edge_patch) + col] += Nval; //  = ⇒ += に変更
                */
            }

            // for debug
            #if 1
            if (!Jcp_flag) std::cout << std::endl << "\n=== Shape Function Values ===" << std::endl << std::endl;
            std::cout << "Jcp " << Jcp_flag << ": ";
            for (int i = 0; i < getcp_size(edge_patch); i++) 
            {
                std::cout << N_matrix[Jcp * getcp_size(edge_patch) + i] << " ";
            }
            std::cout << std::endl;
            #endif
            Jcp_flag++;
        }
        cout << endl;

        // (N^T)^(-1) * N^T の計算
        calcNTN_inverse_Patchwise(N_matrix, NTN_inverse_NT_matrix, edge_patch);

        // ((N^T)^(-1) * N^T) * Jcp_coordinate の計算
        calcControlPointCoordinate_Patchwise(crackData, crackFrontCoord, NTN_inverse_NT_matrix, edge_patch, info, selected_domain);
    }

    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "cp_list[" << i << "] = " << cp_list[i] << " : ";
        for (int d = 0; d < info->DIMENSION; d++) 
        {
            std::cout << std::scientific << std::setprecision(6)
                      << crack_front_cp_coordinate[i * info->DIMENSION + d];
            if (d < info->DIMENSION - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }

    // パッチ境界上の制御点座標の補正
    correctControlPointCoordinate_Patchwise(info);

    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "cp_list[" << i << "] = " << cp_list[i] << " : ";
        for (int d = 0; d < info->DIMENSION; d++) 
        {
            std::cout << std::scientific << std::setprecision(6)
                      << crack_front_cp_coordinate[i * info->DIMENSION + d];
            if (d < info->DIMENSION - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }

    // 重複制御点IDの制御点座標データの平均化と削除
    deleteDuplicateCP(crackData, info);
}

void LeastSquaresApproximation::calcNTN_inverse_Patchwise(std::vector<double>& N_matrix, std::vector<double>& NTN_inverse_NT_matrix,int edge_patch)
{
    std::vector<double> NTN_matrix(cp_size[edge_patch] * cp_size[edge_patch], 0.0);              // (N^T) * N 行列 (cp * cp)
    std::vector<double> NTN_inverse_matrix(cp_size[edge_patch] * cp_size[edge_patch], 0.0);      // (N^T)^(-1) 行列 (cp * cp)
    
    // for debug
    #if 0
    std::cout << "=== Shape Function Values in caclNTN_inverse_NT function===" << std::endl << std::endl;
    for (int Jcp = 0; Jcp < getJcp_size(edge_patch); Jcp++)
    {
        for (int i = 0; i < getcp_size(edge_patch); i++) 
        {
            std::cout << N_matrix[Jcp * getcp_size(edge_patch) + i] << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // (N^T) * N 行列の計算
    for (int i = 0; i < cp_size[edge_patch]; i++) 
    {
        for (int j = 0; j < cp_size[edge_patch]; j++) 
        {
            for (int k = 0; k < Jcp_size[edge_patch]; k++) 
            {
                NTN_matrix[i * cp_size[edge_patch] + j] +=  N_matrix[k * cp_size[edge_patch] + i] * N_matrix[k * cp_size[edge_patch] + j];
            }
        }
    }
    
    // for debug
    #if 0
    std::cout << std::endl << "=== (N^T)*N Matrix ===" << std::endl << std::endl;
    for (int i = 0; i < cp_size[edge_patch]; i++) 
    {
        for (int j = 0; j < cp_size[edge_patch]; j++) 
        {
            std::cout << NTN_matrix[i * cp_size[edge_patch] + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    #endif

    // N^TNの逆行列を計算
    calcInverseMatrix_byEigen_Patchwise(NTN_matrix, NTN_inverse_matrix, edge_patch);

    // (N^T N)^(-1) * N^T 行列の計算
    for (int i = 0; i < cp_size[edge_patch]; i++) 
    {
        for (int j = 0; j < Jcp_size[edge_patch]; j++) 
        {
            for (int k = 0; k < cp_size[edge_patch]; k++) 
            {
                NTN_inverse_NT_matrix[i * Jcp_size[edge_patch] + j] += NTN_inverse_matrix[i * cp_size[edge_patch] + k] * N_matrix[j * cp_size[edge_patch] + k];
            }
        }
    }

    // for debug
    #if 1
    std::cout << std::endl << "=== (N^T)^(-1) * N^T Matrix ===" << std::endl << std::endl;
    for (int i = 0; i < cp_size[edge_patch]; i++) 
    {
        for (int j = 0; j < Jcp_size[edge_patch]; j++) 
        {
            std::cout << NTN_inverse_NT_matrix[i * Jcp_size[edge_patch] + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    #endif
}

void LeastSquaresApproximation::calcControlPointCoordinate_Patchwise(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, vector<double>& NTN_inverse_NT_matrix, int edge_patch, information *info, int selected_domain)
{
    int dim = info->DIMENSION;
    int num_cp = cp_size[edge_patch];
    int num_jcp = Jcp_size[edge_patch];

    std::vector<double> temp_crack_front_points(num_jcp * dim, 0.0);

    // 対応するJ積分点の座標を取得
    int Jcp_offset = 0;
    for (int i = 0; i < edge_patch; i++) Jcp_offset += Jcp_size[i];

    for (int j = 0; j < num_jcp; j++) {
        for (int d = 0; d < dim; d++) {
            int idx = j * dim + d;
            int global_jcp_idx = Jcp_offset + j;
            int global_idx = selected_domain * (crackData.getJintegralPoints() * dim) + global_jcp_idx * dim + d;
            temp_crack_front_points[idx] = crackFrontCoord.getafterCrackFrontPoints()[global_idx];
        }
    }

    // 制御点座標の計算（NTN_inv_NT * Jcp）
    for (int i = 0; i < num_cp; i++) {
        for (int d = 0; d < dim; d++) {
            double value = 0.0;
            for (int j = 0; j < num_jcp; j++) {
                double N = NTN_inverse_NT_matrix[i * num_jcp + j];
                value += N * temp_crack_front_points[j * dim + d];
            }

            // cp_listのインデックスを取得（edge_patchにおけるi番目）
            int cp_id = cp_list[patch_cp_start_index[edge_patch] + i];
            if (cp_global_index_map.count(cp_id) == 0) {
                std::cerr << "[ERROR] cp_id " << cp_id << " not found in cp_global_index_map!" << std::endl;
                continue;
            }
            

            int global_cp_idx = patch_cp_start_index[edge_patch] + i;
            crack_front_cp_coordinate[global_cp_idx * dim + d] = value;
        }
    }
    
    #if 0
    for(int i = 0; i < num_cp; i++){
        int cp_id = cp_list[patch_cp_start_index[edge_patch] + i];
        int global_cp_idx = cp_global_index_map[cp_id];
        std::cout << "CP " << i << " "<< crack_front_cp_coordinate[global_cp_idx * dim + 0] <<" "<< crack_front_cp_coordinate[global_cp_idx * dim + 1] << " "<< crack_front_cp_coordinate[global_cp_idx * dim + 2] << std::endl;
    }
    #endif
}

void LeastSquaresApproximation::correctControlPointCoordinate_Patchwise(information *info)
{
    const int dim = info->DIMENSION;

    
    std::vector<int> correction_patch_set, edge_parameter_set, fixed_direction_set;
    ReadCPSettingLSA("J_int.txt", correction_patch_set, edge_parameter_set, fixed_direction_set);
    
    for (size_t i = 0; i < correction_patch_set.size(); i++)
    {
        int patch = correction_patch_set[i];
        int is_end = edge_parameter_set[i];
        int dir = fixed_direction_set[i];

        // パッチ先頭インデックス + 始点(0) or 終点(cp_size-1)
        int local_id = (is_end == 0) ? 0 : cp_size[patch] - 1;
        size_t cp_index_in_list = patch_cp_start_index[patch] + local_id;
        int cp_id = cp_list[cp_index_in_list];

        #if 1
        std::cout << cp_id << " "<< cp_index_in_list <<" " << crack_front_cp_coordinate[cp_index_in_list * dim + dir]<<" " << info->Node_Coordinate[cp_id * (dim + 1) + dir] << endl;
        #endif

        // 指定方向だけ補正
        crack_front_cp_coordinate[cp_index_in_list * dim + dir] =
            info->Node_Coordinate[cp_id * (dim + 1) + dir];

    }

    // debug出力（必要なら有効化）
    #if 1
    std::cout << "\n=== Control Points After Patchwise Correction ===" << std::endl;
    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "Control Point " << i << " (ID: " << cp_list[i] << "): ";
        for (int d = 0; d < dim; d++) {
            std::cout << crack_front_cp_coordinate[i * dim + d] << (d < dim - 1 ? ", " : "\n");
        }
    }
    #endif
}

void LeastSquaresApproximation::calcInverseMatrix_byEigen_Patchwise(std::vector<double>& NTN_matrix, std::vector<double>& NTN_inverse_matrix, int edge_patch)
{
    //edge_patchのCP数
    int n = cp_size[edge_patch];

    // NTN行列をEigen形式に変換
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> NTN = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
                                                                (NTN_matrix.data(), n, n);
    
    // LU分解による逆行列計算
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Identity = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(n, n);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> NTN_inverse = NTN.lu().solve(Identity);

    // 結果を標準のvectorに戻す
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(NTN_inverse_matrix.data(), n, n) = NTN_inverse;

    // for debug
    #if 0
    std::cout << "\n=== NTN Inverse Matrix ===\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << NTN_inverse_matrix[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // NTN * NTN^(-1)が単位行列になることを確認
    #if 0
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result = NTN * NTN_inverse;
    std::cout << "\n=== Verification (NTN * NTN^(-1)) ===\n" << result << std::endl;
    #endif
}

void LeastSquaresApproximation::deleteDuplicateCP(const CrackDataSetting_3D& crackData, information *info)
{   
    std::map<int, std::vector<double>> unique_cp_coords;  // 制御点ID -> 座標値のリスト
    std::map<int, int> cp_count;                          // 制御点ID -> 出現回数
    int cp_list_offset = 0;                               // 

    // 各パッチの制御点を処理
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {
        // パッチ内の各制御点を処理
        for (int i = 0; i < cp_size[edge_patch]; i++) 
        {
            int cp_id = cp_list[cp_list_offset + i];
            
            // 座標値の取得
            std::vector<double> coord(info->DIMENSION);
            for (int dim = 0; dim < info->DIMENSION; dim++) 
            {
                coord[dim] = crack_front_cp_coordinate[(cp_list_offset + i) * info->DIMENSION + dim];
            }

            // 重複制御点の処理
            if (unique_cp_coords.find(cp_id) == unique_cp_coords.end()) 
            {
                // 新しい制御点の場合
                unique_cp_coords[cp_id] = coord;
                cp_count[cp_id] = 1;
            } else {
                // 既存の制御点の場合、座標を加算
                for (int dim = 0; dim < info->DIMENSION; dim++) {
                    unique_cp_coords[cp_id][dim] += coord[dim];
                }
                cp_count[cp_id]++;
            }
        }
        cp_list_offset += cp_size[edge_patch];
    }

    // 重複を除いた新しいcp_listと座標値の作成
    std::vector<int> new_cp_list;
    std::vector<double> new_coordinates;

    // 元の順序を維持しながら重複を除去
    cp_list_offset = 0;
    for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch(); edge_patch++) 
    {
        for (int i = 0; i < cp_size[edge_patch]; i++) 
        {
            int cp_id = cp_list[cp_list_offset + i];
            
            // この制御点IDが初めて出現する場合のみ処理
            if (std::find(new_cp_list.begin(), new_cp_list.end(), cp_id) == new_cp_list.end()) 
            {
                new_cp_list.push_back(cp_id);
                
                // 平均座標値の計算と追加
                const auto& coord = unique_cp_coords[cp_id];
                for (int dim = 0; dim < info->DIMENSION; dim++) {
                    new_coordinates.push_back(coord[dim] / cp_count[cp_id]);
                }
            }
        }
        cp_list_offset += cp_size[edge_patch];
    }

    // メンバ変数の更新
    cp_list = std::move(new_cp_list);
    crack_front_cp_coordinate = std::move(new_coordinates);

    // for debug
    #if 0
    std::cout << "\n=== Final Control Points After Duplicate Removal ===" << std::endl;
    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "Control Point " << i << ", ID: " << cp_list[i] << ": ";
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {
            std::cout << crack_front_cp_coordinate[i * info->DIMENSION + dim];
            if (dim < info->DIMENSION - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << endl;
    #endif
}

void LeastSquaresApproximation::ReadCPSettingLSA(const std::string& filename, vector<int> &correction_patch_set, vector<int> &edge_parameter_set, vector<int> &direction_set)
{
    cout << "=== Read Control Point Setting ===" << endl;

    // ファイルの読み込み
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(1);
    }

    // 積分パラメータをスキップ
    int domainWidth, JintegralPoints, model_size_para_3D, domain_types, domain_area_divisions, patchType;
    int num_edge_crack_patch, num_singular_patch, patchID;

    file >> domainWidth >> JintegralPoints >> model_size_para_3D >> domain_types >> domain_area_divisions >>patchType;
    file >> num_edge_crack_patch >> num_singular_patch;

    for (int i = 0; i < num_edge_crack_patch * num_singular_patch; i++)
    {
        file >> patchID;
    }

    // き裂モードをスキップ
    int mode;
    file >> mode;

    // Paris則パラメータのスキップ
    double C_val, m_val;
    file >> C_val >> m_val;

    // 最大き裂進展長さのスキップ
    double max_growth;
    file >> max_growth;

    // 繰り返し荷重サイクル数の計算選択をスキップ
    // int selectNk;
    // file >> selectNk;

    // き裂進展量に対して荷重移動平均を計算するかをスキップ
    int weight_move_ave, num_move_ave;
    file >> weight_move_ave >> num_move_ave;

    // 積分領域の選択をスキップ
    int domainval;
    file >> domainval;
    
    // 制御点の補正を行う制御点と方向のセット数を取得
    int correction_cp_num;
    file >> correction_cp_num;

    // 制御点の補正を行う制御点と方向のセットを取得
    for (int i = 0; i < correction_cp_num; i++)
    {
        int corection_patch, edge_parameter, fixed_direction;
        file >> corection_patch >> edge_parameter >> fixed_direction;
        correction_patch_set.push_back(corection_patch);
        edge_parameter_set.push_back(edge_parameter);
        direction_set.push_back(fixed_direction);
    }
    // for debug and error handring
    #if 1
    std::cout << std::endl << "Correction CP Num: " << correction_cp_num << std::endl;
    for (int i = 0; i < correction_cp_num; i++)
    {   
        if (edge_parameter_set[i] != 0 && edge_parameter_set[i] != 1) 
        {
            std::cerr << "Error: edge parameter must be 0 or 1" << std::endl;
            exit(1);
        }
        if (direction_set[i] != 0 && direction_set[i] != 1 && direction_set[i] != 2) 
        {
            std::cerr << "Error: direction must be 0, 1 or 2" << std::endl;
            exit(1);
        }
        std::cout << "Correction Patch: " << correction_patch_set[i] << 
                     ", edge parameter: " << edge_parameter_set[i] <<
                     ", fixed Direction: " << direction_set[i] << std::endl;
    }
    #endif
    std::cout << std::endl;

    file.close();
}


void Remeshing::DoRemeshing(const CrackDataSetting_3D& crackData, LeastSquaresApproximation calc_cp_LSA, information *info)
{   
    int offset = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;              // 前パッチまでの自由度数
    int glo_patch_n = info->Total_Patch_to_mesh[Total_mesh - 1];                               // ローカルパッチ番号の初期値
    vector<double> u_vec(info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION, 0);  // リメッシュ後の各パッチの端部制御点座標
    vector<vector<double>> new_coord(info->Total_Patch_to_mesh[Total_mesh]);                   // リメッシング後の制御点座標
    vector<vector<double>> new_w(info->Total_Patch_to_mesh[Total_mesh]);                       // リメッシング後の重み

    // リメッシングの設定を読み込み
    readRemeshPatch("J_int.txt");

    // き裂前縁パッチ, 後方パッチ，リメッシングパッチをローカルパッチ番号に変換
    for (size_t i = 0; i < patch_around_crack.size(); i++) 
    {
        patch_around_crack[i] += glo_patch_n;
    }
    for (size_t i = 0; i < patch_behind_crack.size(); i++) 
    {
        patch_behind_crack[i] += glo_patch_n;
    }
    for (size_t i = 0; i < patch_remesh.size(); i++) 
    {
        patch_remesh[i] += glo_patch_n;
    }

    // for debug
    std::cout << "=== patch_around_crack offset global ===" << std::endl;
    for (size_t i = 0; i < patch_around_crack.size(); i++)
    {
        std::cout << patch_around_crack[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "=== patch_behind_crack offset global ===" << std::endl;
    for (size_t i = 0; i < patch_behind_crack.size(); i++)
    {
        std::cout << patch_behind_crack[i] << std::endl;
    }
    std::cout << "=== patch_remeshoffset global  ===" << std::endl;
    for (size_t i = 0; i < patch_remesh.size(); i++)
    {
        std::cout << patch_remesh[i] << std::endl;
    }
    std::cout << std::endl;

    // き裂前縁の各制御点の移動量を計算
    calcCPmoveVector(calc_cp_LSA, info);

    // パッチのジオメトリの再定義
    RedefinePatchGeometry(crackData, calc_cp_LSA, info, u_vec, patch_around_crack, offset, new_coord, new_w);
    std::cout << "finish RedefinePatchGeometry" << std::endl;

    // パッチ内部の制御点座標の計算
    for (size_t i = 0; i < patch_remesh.size(); i++)
    {
        int p = patch_remesh[i];
    
        // set MAX_INCREASED_CP
        int MAX_INCREASED_CP_total = 1;
        for (int j = 0; j < info->DIMENSION; j++)
        {
            int temp = info->No_Control_point[p * info->DIMENSION + j];
            MAX_INCREASED_CP_total *= temp;
        }

        // set MAX_INCREASED_KV
        int MAX_INCREASED_KV_1D = 0;
        for (int j = 0; j < info->DIMENSION; j++)
        {
            int temp = info->No_knot[p * info->DIMENSION + j];
            if (remesh_direction[i] == j)
            {
                temp += division_ele_n[i] + 2 * info->Order[p * info->DIMENSION + j] + 1;
                // temp += (division_cp_n[i] - info->Order[p * info->DIMENSION + j]) + 2 * info->Order[p * info->DIMENSION + j] + 1;
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
            if (!(info->INC[p * offset + cp_num * info->DIMENSION + remesh_direction[i]] == 0 || info->INC[p * offset + cp_num * info->DIMENSION + remesh_direction[i]] 
                    == info->No_Control_point[p * info->DIMENSION + remesh_direction[i]] - 1))
                continue;

            // glo.Weight[count++] = info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + info->DIMENSION];
            glo.Weight[count++] = 1.0000000000000000e+00;
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
            if (j == remesh_direction[i])
                temp.CP_n = 2;
            temp.Order = info->Order[p * info->DIMENSION + j];
            if (j == remesh_direction[i])
                temp.Order = 1;
            temp.OE_n = 0;
            if (j == remesh_direction[i])
                temp.OE_n = info->Order[p * info->DIMENSION + j] - 1;
            temp.knot_n = info->No_knot[p * info->DIMENSION + j];
            if (j == remesh_direction[i])
                temp.knot_n = 4;
            if (j == remesh_direction[i])
                temp.KI_cp_n = division_ele_n[i] + info->Order[p * info->DIMENSION + j];
                // temp.KI_cp_n = division_cp_n[i];
            else
                temp.KI_cp_n = 0;
            temp.KI_non_uniform_n = 0;
            temp.CP = CP[j].data();
            for (size_t k = 0; k < new_coord[p].size() / info->DIMENSION; k++)
                temp.CP[k] = new_coord[p][k * info->DIMENSION + j];
            temp.KV = KV[j].data();
            for (int k = 0; k < info->No_knot[p * info->DIMENSION + j]; k++)
                temp.KV[k] = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + j] + k];
            if (j == remesh_direction[i])
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
            if (j != remesh_direction[i])
                continue;

            // order elevation
            OE(j, &glo, &nurbs[j]);

            // knot insertion
            KI_cp(j, &glo, &nurbs[j]);
        }
        // cout << "finish order elevation and knot insertion" << endl;

        // update new coord
        new_coord[p].resize(glo.Total_Control_Point * info->DIMENSION);
        for (int j = 0; j < glo.Total_Control_Point; j++)
            for (int k = 0; k < info->DIMENSION; k++)
                new_coord[p][j * info->DIMENSION + k] = nurbs[k].CP[j];


        // cout << "finish update new coord" << endl;        
    }
    cout << "finish remeshing" << endl;

    // for debug
    #if 0
    std::cout << "\n=== Remeshing Data ===" << std::endl;
    for (size_t i = 0; i < new_coord.size(); i++) 
    {
        std::cout << "Patch " << i << ": " << std::endl;
        for (size_t j = 0; j < new_coord[i].size(); j += info->DIMENSION) 
        {
            std::cout << "Control Point " << j / info->DIMENSION << ": ";
            for (int k = 0; k < info->DIMENSION; k++) 
            {
                std::cout << new_coord[i][j + k];
                if (k < info->DIMENSION - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }
    #endif

    // out put input file for next step analysis
    makeInputfileNextStep(info, new_coord);
    std::cout << "finish makeInputfileNextStep" << std::endl;

    // output as vtk file 
    outputRemeshingData(info, new_coord);
    std::cout << "finish outputRemeshingData" << std::endl;

    // output as some patch vtk file 
    // outputEachPatchData(info);
}


void Remeshing::readRemeshPatch(const std::string& filename)
{   
    cout << "=== Read Remesh Patch Setting ===" << endl;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(1);
    }

    // 積分パラメータをスキップ
    int domainWidth, JintegralPoints, model_size_para_3D, domain_types, domain_area_divisions, patchType;
    int num_edge_crack_patch, num_singular_patch, patchID;

    file >> domainWidth >> JintegralPoints >> model_size_para_3D >> domain_types >> domain_area_divisions >>patchType;
    file >> num_edge_crack_patch >> num_singular_patch;

    for (int i = 0; i < num_edge_crack_patch * num_singular_patch; i++)
    {
        file >> patchID;
    }

    // き裂モードをスキップ
    int mode;
    file >> mode;

    // Paris則パラメータのスキップ
    double C_val, m_val;
    file >> C_val >> m_val;

    // 最大き裂進展長さのスキップ
    double max_growth;
    file >> max_growth;

    // 繰り返し荷重サイクル数の計算選択をスキップ
    // int selectNk;
    // file >> selectNk;

    // き裂進展量に対して荷重移動平均を計算するかをスキップ
    int weight_move_ave, num_move_ave;
    file >> weight_move_ave >> num_move_ave;

    // 積分領域の選択をスキップ
    int domainval;
    file >> domainval;
    
    // 制御点の補正を行う制御点と方向のセット数をスキップ
    int correction_cp_num;
    file >> correction_cp_num;

    // 制御点の補正を行う制御点と方向のセットをスキップ
    for (int i = 0; i < correction_cp_num; i++)
    {
        int corection_patch, edge_parameter, fixed_direction;
        file >> corection_patch >> edge_parameter >> fixed_direction;
    }

    // リメッシング関連のデータを読み込み
    std::string line;
    int temp_value;

    // き裂前縁パッチの読み込み
    std::getline(file >> std::ws, line);
    std::istringstream iss_around(line);
    while (iss_around >> temp_value) {
        patch_around_crack.emplace_back(temp_value);
    }
    for (size_t around = 0; around < patch_around_crack.size(); around++)
    {
        std::cout << "patch_around_crack[" << around << "]: " << patch_around_crack[around] << std::endl;
    }

    // き裂前縁後方パッチの読み込み
    std::getline(file >> std::ws, line);
    std::istringstream iss_behind(line);
    while (iss_behind >> temp_value) {
        patch_behind_crack.emplace_back(temp_value);
    }
    for (size_t behind = 0; behind < patch_behind_crack.size(); behind++)
    {
        std::cout << "patch_behind_crack[" << behind << "]: " << patch_behind_crack[behind] << std::endl;
    }

    // リメッシュパッチの読み込み
    std::getline(file >> std::ws, line);
    std::istringstream iss_remesh(line);
    while (iss_remesh >> temp_value) {
        patch_remesh.emplace_back(temp_value);
    }
    for (size_t remesh = 0; remesh < patch_remesh.size(); remesh++)
    {
        std::cout << "patch_remesh[" << remesh << "]: " << patch_remesh[remesh] << std::endl;
    }

    // リメッシュ方向の読み込み
    std::getline(file >> std::ws, line);
    std::istringstream iss_direction(line);
    while (iss_direction >> temp_value) 
    {
        remesh_direction.emplace_back(temp_value);
    }
    for (size_t direction = 0; direction < remesh_direction.size(); direction++)
    {
        std::cout << "remesh_direction[" << direction << "]: " << remesh_direction[direction] << std::endl;
    }

    // 要素数の読み込み
    std::getline(file >> std::ws, line);
    std::istringstream iss_division(line);
    while (iss_division >> temp_value) {
        division_ele_n.emplace_back(temp_value);
    }
    for (size_t division = 0; division < division_ele_n.size(); division++)
    {
        std::cout << "division_ele_n[" << division << "]: " << division_ele_n[division] << std::endl;
    }

    file.close();
}


void Remeshing::calcCPmoveVector(LeastSquaresApproximation calc_cp_LSA, information *info) 
{
    const auto& cp_list = calc_cp_LSA.getcp_list();
    const auto& after_coordinates = calc_cp_LSA.getControlPointCoordinate();

    // 各制御点の移動量を計算
    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        int cp_id = cp_list[i];
        // std::cout << "Control Point " << cp_id << ": ";

        std::cout << "[copy] cp_id = " << cp_id  << std::endl;
        
        // 各座標成分(x,y,z)の移動量を計算
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {
            double before_coord = info->Node_Coordinate[cp_id * (info->DIMENSION + 1) + dim];
            double after_coord = after_coordinates[i * info->DIMENSION + dim];  
            cout << "before_coord = " << before_coord << "after_coord = " << after_coord << endl;

            // 移動量を計算 
            cp_move_vector[i * info->DIMENSION + dim] = after_coord - before_coord;
        }
    }
    // for debug
    #if 0
    std::cout << "\n=== Control Point move vector ===" << std::endl;
    for (size_t i = 0; i < cp_list.size(); i++) 
    {
        std::cout << "Control Point " << cp_list[i] << ": ";
        for (int dim = 0; dim < info->DIMENSION; dim++) 
        {
            printf("%.20e", cp_move_vector[i * info->DIMENSION + dim]);
            if (dim < info->DIMENSION - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    #endif
}


void Remeshing::RedefinePatchGeometry(const CrackDataSetting_3D& crackData, LeastSquaresApproximation calc_cp_LSA, information *info, vector<double>& u_vec, vector<int>& patch_around_crack, int offset, vector<vector<double>>& new_coord, vector<vector<double>>& new_w)
{
    vector<bool> isSet(info->Total_Control_Point_to_mesh[Total_mesh], false); // 制御点が設定されたかどうかのフラグ

    // patch_around_crack の各パッチ境界に対して対応する移動量を与える
    for (size_t i = 0; i < patch_around_crack.size(); i++)
    {
        int p = patch_around_crack[i]; //現在のパッチ番号
        int cp_num_offset = info->Total_Control_Point_to_patch[p]; //パッチ p の先頭番号
        int crack_edge_dir = 1;  // き裂前縁方向のインデックス

        // 各パッチの制御点のき裂前縁方向の始点のインデックスを取得
        int segment_size = patch_around_crack.size() / crackData.getnum_edge_crack_patch();
        int segment_index = i / segment_size; //何番目のセグメントか
        int offset_move_vec = 0;
        if (segment_index > 0)
        {
            for (int j = 0; j < segment_index; j++) 
            {
                offset_move_vec += calc_cp_LSA.getcp_size(j) - 1; //
            }
        }
        std::cout << "offset_move_vec: " << offset_move_vec << std::endl;

        for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
        {
            int connectivity_num = info->Patch_Control_point[cp_num_offset + j]; //制御点のグローバル番号

            if (isSet[connectivity_num]) 
                continue;

            // 制御点座標に移動量を与える
            for (int dim = 0; dim < info->DIMENSION; dim++)
            {
                // 上面と下面でインデックスが異なるため分岐
                int index = 0;
                // 上面 または下面の特異パッチの場合
                // if (MODEL_SIZE_PARA_3D == 1 || ( MODEL_SIZE_PARA_3D == 0 && ((p - info->Total_Patch_to_mesh[1]) < (info->Total_Patch_on_mesh[1] / 2 ) + crackData.getpatchList().size() / 2)  ))
                // if (MODEL_SIZE_PARA_3D == 1 || (MODEL_SIZE_PARA_3D == 0 && ((p - info->Total_Patch_to_mesh[1]) < (info->Total_Patch_on_mesh[1] / 2) + static_cast<int>(crackData.getpatchList().size() / 2))))
                if (crackData.getmodel_size_para_3D() == 1 || (crackData.getmodel_size_para_3D() == 0 && ((p - info->Total_Patch_to_mesh[1]) < (info->Total_Patch_on_mesh[1] / 2) + static_cast<int>(crackData.getpatchList().size() / 2))))
                {
                    index = info->INC[p * offset + connectivity_num * info->DIMENSION + crack_edge_dir];
                }
                // フルモデルであれば下面
                else
                {
                    int last_connectivity_num = info->Patch_Control_point[cp_num_offset + info->No_Control_point_in_patch[p] - 1];
                    index = info->INC[p * offset + last_connectivity_num * info->DIMENSION + crack_edge_dir] - info->INC[p * offset + connectivity_num * info->DIMENSION + crack_edge_dir];
                }
                
                double u = cp_move_vector[(offset_move_vec + index) * info->DIMENSION + dim];
                // std::cout << "u: " << u << std::endl;
                u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] + u;
                #if 0
                std::cout << "[DEBUG] p = " << p 
                << ", cp_num = " << connectivity_num 
                << ", segment_index = " << segment_index 
                << ", index = " << index 
                << ", offset_move_vec = " << offset_move_vec 
                << ", final_index = " << (offset_move_vec + index) 
                << ", dim = " << dim 
                << ", u = " << cp_move_vector[(offset_move_vec + index) * info->DIMENSION + dim]
                << ", u_vec = " << u_vec[connectivity_num * info->DIMENSION + dim]
                << std::endl;
                #endif

            }
            isSet[connectivity_num] = true;
        }
    }

    // patch_behind_crack の各パッチ境界に対して対応する移動量を与える
    #if 0
    vector<double> max(2, -1.0e+10);  // x, y 方向の最大座標
    vector<double> min(2, 1.0e+10);   // x, y 方向の最小座標

    // き裂前縁の後方のパッチの座標の最大値と最小値を取得
    for (size_t i = 0; i < patch_behind_crack.size(); i++)
    {
        int p = patch_behind_crack[i];
        if (i < 2) 
        {
            for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
            {
                int connectivity_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[p] + j];
                if (max[1] < info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1])
                    max[1] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1];
                if (min[1] > info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1])
                    min[1] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1];
            }
        }
        else
        {
            for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
            {
                int connectivity_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[p] + j];
                if (max[0] < info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0])
                    max[0] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0];
                if (min[0] > info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0])
                    min[0] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0];
            }
        }
    }
    #endif

    vector<double> max(2, -1.0e+10);  // x, y 方向の最大座標
    vector<double> min(2, 1.0e+10);   // x, y 方向の最小座標

    // き裂前縁の後方のパッチの座標の最大値と最小値を取得
    for (size_t i = 0; i < patch_behind_crack.size(); i++)
    {
        int p = patch_behind_crack[i];

        // き裂前縁方向に0インデックスのパッチ群からy方向の最大値と最小値を取得
        if (i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch())  
        {   
            for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
            {
                int connectivity_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[p] + j];
                if (max[1] < info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1])
                    max[1] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1];
                if (min[1] > info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1])
                    min[1] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 1];
            }
        }
        // き裂前縁方向に1インデックスのパッチ群からx方向の最大値と最小値を取得
        else if (i >= patch_behind_crack.size() / crackData.getnum_edge_crack_patch() && i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 2)  
        {   
            for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
            {
                int connectivity_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[p] + j];
                if (max[0] < info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0])
                    max[0] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0];
                if (min[0] > info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0])
                    min[0] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + 0];
            }
        }
    }
    // for debug
    #if 0
    std::cout << "max[0]: " << max[0] << ", min[0]: " << min[0] << std::endl;
    std::cout << "max[1]: " << max[1] << ", min[1]: " << min[1] << std::endl;
    std::cout << std::endl;
    #endif

    // き裂前縁の後方のパッチの制御点に対して移動量を与える
    for (size_t i = 0; i < patch_behind_crack.size(); i++)
    {
        int p = patch_behind_crack[i];
        int cp_num_offset = info->Total_Control_Point_to_patch[p];
        int new_surface_dir = 0;  // 更新するパラメータ空間の面に垂直方向のインデックス

        // パッチ内の各制御点について処理
        for (int j = 0; j < info->No_Control_point_in_patch[p]; j++)
        {
            int connectivity_num = info->Patch_Control_point[cp_num_offset + j];
            if (!(info->INC[p * offset + connectivity_num * info->DIMENSION + new_surface_dir] == info->No_Control_point[p * info->DIMENSION + new_surface_dir] - 1))
                continue;

            if (isSet[connectivity_num]) 
                continue;

            // 制御点移動に使用するき裂前縁方向のインデックスを取得（90°モデルなら始点と終点，180°モデルなら始点と90°位置と終点）
            std::vector<int> move_vec_index;  
            move_vec_index.emplace_back(0);  // 始点のインデックス

            // き裂前縁方向に45°のパッチが2つ分で90°のモデルの場合
            if (crackData.getnum_edge_crack_patch() == 2) 
            {   
                int temp_move_vec_indexd = calc_cp_LSA.getcp_list().size() - 1;
                move_vec_index.emplace_back(temp_move_vec_indexd);  // 終点のインデックス
            }
            // き裂前縁方向に45°のパッチが4つ分で180°のモデルの場合
            else if (crackData.getnum_edge_crack_patch() == 4) 
            {   
                int temp_move_vec_index = 0;
                for (int edge_patch = 0; edge_patch < crackData.getnum_edge_crack_patch() / 2; edge_patch++) 
                {
                    temp_move_vec_index += calc_cp_LSA.getcp_size(edge_patch);
                    temp_move_vec_index = temp_move_vec_index - 1;
                }
                move_vec_index.emplace_back(temp_move_vec_index);                  // 90°位置のインデックス
                move_vec_index.emplace_back(calc_cp_LSA.getcp_list().size() - 1);  // 終点のインデックス
            }
            #if 0
            for (int k = 0; k < move_vec_index.size(); k++)
            {
                std::cout << "move_vec_index[" << k << "]: " << move_vec_index[k] << std::endl;
            }
            #endif

            // 各パッチに対して同一のき裂前縁方向のインデックスに制御点移動量を与える
            if (i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch())  // 0~45°のパッチ
            {   
                for (int dim = 0; dim < info->DIMENSION; dim++)
                {
                    if (dim == 0)
                    {
                        double u = cp_move_vector[move_vec_index[0] * info->DIMENSION + dim];
                        // std::cout << "u_11: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] + u;
                    }
                    else if (dim == 1)
                    {
                        double l = max[1] - min[1];
                        double u = (l + cp_move_vector[move_vec_index[1] * info->DIMENSION + dim]) * ((info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] - min[1]) / l);
                        // std::cout << "u_12: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = min[1] + u;
                    }
                    else if (dim == 2)
                    {
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim];
                    }
                }
            }
            else if (i >= patch_behind_crack.size() / crackData.getnum_edge_crack_patch() && i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 2) // 45° ~ 90°のパッチ
            {   
                for (int dim = 0; dim < info->DIMENSION; dim++)
                {
                    if (dim == 0)
                    {
                        double l = max[0] - min[0];
                        double u = (l + cp_move_vector[move_vec_index[0] * info->DIMENSION + dim]) * ((info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] - min[0]) / l);
                        // std::cout << "u_21: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = min[0] + u;
                    }
                    else if (dim == 1)
                    {
                        double u = cp_move_vector[move_vec_index[1] * info->DIMENSION + dim];
                        // std::cout << "u_22: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] + u;
                    }
                    else if (dim == 2)
                    {
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim];
                    }
                }
            }
            else if (i >= patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 2 && i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 3)  // 90° ~ 135°のパッチ
            {   
                for (int dim = 0; dim < info->DIMENSION; dim++)
                {
                    if (dim == 0)
                    {
                        double l = max[0] - min[0];
                        double u = (l - cp_move_vector[move_vec_index[2] * info->DIMENSION + dim]) * ((info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] - min[0]) / l);
                        // double u = (l + cp_move_vector[move_vec_index[0] * info->DIMENSION + dim]) * ((info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] - min[0]) / l);
                        // std::cout << "u_31: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = min[0] + u;
                    }
                    else if (dim == 1)
                    {
                        double u = cp_move_vector[move_vec_index[1] * info->DIMENSION + dim];
                        // std::cout << "u_32: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] + u;
                    }
                    else if (dim == 2)
                    {
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim];
                    }
                }
            }
           else if (i >= patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 3 && i < patch_behind_crack.size() / crackData.getnum_edge_crack_patch() * 4)  // 135° ~ 180°のパッチ
            {   
                for (int dim = 0; dim < info->DIMENSION; dim++)
                {
                    if (dim == 0)
                    {
                        double u = cp_move_vector[move_vec_index[2] * info->DIMENSION + dim];
                        // double u = cp_move_vector[move_vec_index[0] * info->DIMENSION + dim];
                        // std::cout << "u_41: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] + u;
                    }
                    else if (dim == 1)
                    {
                        double l = max[1] - min[1];
                        double u = (l + cp_move_vector[move_vec_index[1] * info->DIMENSION + dim]) * ((info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim] - min[1]) / l);
                        // std::cout << "u_42: " << u << std::endl;
                        u_vec[connectivity_num * info->DIMENSION + dim] = min[1] + u;
                    }
                    else if (dim == 2)
                    {
                        u_vec[connectivity_num * info->DIMENSION + dim] = info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + dim];
                    }
                }
            }
            isSet[connectivity_num] = true;
        }
    }

    // patch_around_crack以外に対してジオメトリを再定義
    cout << "=== Redefine Patch Geometry ===" << endl;
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
    
    // update coord
    std::cout << "=== Update Control Point Coordinate ===" << std::endl;
    for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
    {
        auto it = find(patch_remesh.begin(), patch_remesh.end(), i);  
        int index = distance(patch_remesh.begin(), it);               
        int p = patch_remesh[index];  
        int dir = remesh_direction[index];

        for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
        {
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j];

            // 端部の制御点でない場合はスキップ
            if (it != patch_remesh.end() && !(info->INC[p * offset + cp_num * info->DIMENSION + dir] == 0 || info->INC[p * offset + cp_num * info->DIMENSION + dir] == info->No_Control_point[p * info->DIMENSION + dir] - 1))
                continue;
            
            for (int k = 0; k < info->DIMENSION; k++)
            {   
                new_coord[i].emplace_back(u_vec[cp_num * info->DIMENSION + k]);
            }
            new_w[i].emplace_back(info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + info->DIMENSION]);  
        }
    }

    // for debug
    #if 0
    std::cout << "\n=== New Control Point Coordinate ===" << std::endl;
    // 各制御点ごとに改行して表示
    // 制御点IDvx y z 
    for (size_t i = 0; i < new_coord.size(); i++) 
    {
        std::cout << "Patch " << i << ": " << std::endl;
        for (size_t j = 0; j < new_coord[i].size(); j += info->DIMENSION) 
        {
            std::cout << "Control Point " << j / info->DIMENSION << ": ";
            for (int k = 0; k < info->DIMENSION; k++) 
            {
                std::cout << new_coord[i][j + k];
                if (k < info->DIMENSION - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }
    #endif

    // output new_coord size
    for (size_t i = 0; i < new_coord.size(); i++)
    {
        std::cout << "Patch " << i << ": " << new_coord[i].size() / info->DIMENSION << std::endl;
    }
}


void Remeshing::makeInputfileNextStep(information *info, vector<vector<double>>& new_coord)
{
    int local_init_patch = info->Total_Patch_to_mesh[1];       // グローバルパッチ数
    int local_init_cp = info->Total_Control_Point_to_mesh[1];  // グローバル制御点数   
    
    // ファイルの書き出し
    std::ofstream file("input_loc_next_step.txt");
    if (!file) 
    {
        std::cerr << "Error: Cannot open file: input_next_step.txt" << std::endl;
        exit(1);
    }

    // 空間の次元数
    file << info->DIMENSION << std::endl;
    file << std::endl;

    // ヤング率 ポアソン比
    file << E << " " << nu << std::endl;
    file << std::endl;

    // ローカルのパッチ数
    file << info->Total_Patch_to_mesh[Total_mesh] - local_init_patch << std::endl;
    file << std::endl;

    // ローカルの制御点数
    file << info->Total_Control_Point_to_mesh[Total_mesh] - local_init_cp << std::endl;
    file << std::endl;

    // パッチ数分の各方向の次数
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {
            file << info->Order[i * info->DIMENSION + j] << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    // パッチ数分の各方向のノットベクトル数
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {
            file << info->No_knot[i * info->DIMENSION + j] << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    // パッチ数分の各方向の制御点数
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {
            file << info->No_Control_point[i * info->DIMENSION + j] << " ";
        }
        file << std::endl;
    }
    file << std::endl;
    std::cout << "finish output number of control point" << std::endl;

    // 各パッチのコネクティビティ
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        for (int j = 0; j < info->No_Control_point_in_patch[i]; j++) 
        {
            file << info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j] - local_init_cp<< " ";
        }
        file << std::endl;
    }
    file << std::endl;
    std::cout << "finish output connectivity" << std::endl;

    // 変位指定の制御点数 集中荷重の制御点数 分布荷重の制御点数
    file << info->Total_Constraint_to_mesh[Total_mesh] - info->Total_Constraint_to_mesh[1] << " "
         << info->Total_Load_to_mesh[Total_mesh] - info->Total_Load_to_mesh[1] << " "
         << info->Total_DistributeForce_to_mesh[Total_mesh] - info->Total_DistributeForce_to_mesh[1] << std::endl;
    file << std::endl;
    std::cout << "finish output number of constraint and force" << std::endl;    

    // 各パッチの各方向のノットベクトル
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {
            for (int k = 0; k < info->No_knot[i * info->DIMENSION + j]; k++) 
            {
                file << std::scientific << std::setprecision(16) << 
                     info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + j] + k] << " ";
            }
            file << std::endl;
        }
    }
    file << std::endl;
    std::cout << "finish output knot vector" << std::endl;

    vector<bool> isSet(info->Total_Control_Point_to_mesh[Total_mesh], false);
    int temp_count = 0;
    for (int i = local_init_patch; i < info->Total_Patch_to_mesh[Total_mesh]; i++) 
    {
        // for (int j = 0; j < info->No_Control_point_in_patch[i]; j++)
        for (size_t j = 0; j < new_coord[i].size() / info->DIMENSION; j++)
        {
            int cp_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j];
            if (isSet[cp_num]) 
                continue;

            temp_count++;

            file << cp_num - local_init_cp << " ";
            for (int k = 0; k < info->DIMENSION; k++) 
            {
                file << std::scientific << std::setprecision(16) << new_coord[i][j * info->DIMENSION + k] << " ";
            }
            // file << std::scientific << std::setprecision(16) << info->Node_Coordinate[cp_num * (info->DIMENSION + 1) + info->DIMENSION] << std::endl;
            file << std::scientific << std::setprecision(16) << 1.0000000000000000 << std::endl;
            isSet[cp_num] = true;
        }
    }
    file << std::endl;
    

    // 変位制御する制御点番号 方向0 or 1 or 2(整数) 変位量(倍精度)
    int local_constraint_count = info->Total_Constraint_to_mesh[Total_mesh] - info->Total_Constraint_to_mesh[1];
    for (int i = 0; i < local_constraint_count; i++)
    {   
        int idx = info->Total_Constraint_to_mesh[1] + i;
        file << info->Constraint_Node_Dir[2 * idx + 0] - local_init_cp << " "
            << info->Constraint_Node_Dir[2 * idx + 1] << " "
            << info->Value_of_Constraint[idx] << std::endl;
    }
    file << std::endl;

    std::cout << "finish input file for next step" << std::endl;
}


void Remeshing::outputRemeshingData(information *info, vector<vector<double>>& new_coord)
{
    ofstream vtkFile("remesh_cp.vtk");
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
}


void Remeshing::outputEachPatchData(information *info)
{
    for (int i = 1; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
    {
        if (find(patch_remesh.begin(), patch_remesh.end(), i) == patch_remesh.end())
            continue;

        string filename = to_string(i) + "_cp" + ".vtk";
        ofstream vtkFile(filename);
        vtkFile << "# vtk DataFile Version 3.0" << endl;
        vtkFile << "Remesh output" << endl;
        vtkFile << "ASCII" << endl;
        vtkFile << "DATASET POLYDATA" << endl;

        // calc total control point
        int total_cp = info->No_Control_point_in_patch[i];

        // Write points
        vtkFile << "POINTS " << total_cp << " float" << endl;

        for (int j = 0; j < total_cp; j++)
        {
            int connectivity_num = info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j];
            for (int k = 0; k < info->DIMENSION; k++)
            {
                vtkFile << info->Node_Coordinate[connectivity_num * (info->DIMENSION + 1) + k] << " ";
            }
            vtkFile << endl;
        }

        // Write vertices
        vtkFile << "VERTICES " << total_cp << " " << total_cp * 2 << endl;
        for (int j = 0; j < total_cp; j++)
        {
            vtkFile << "1 " << j << endl;
        }

        vtkFile.close();
    }
}