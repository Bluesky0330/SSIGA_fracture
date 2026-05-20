#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <unordered_set>
#include <map>

using std::cout;
using std::endl;
using std::string;
using std::vector;


// 前方宣言
class CrackDirection;
class LeastSquaresApproximation;
class ParisLawCrack;


// 応力拡大係数振幅の計算クラス
class StressIntensityFactorRange {
private:
    std::vector<double> delta_K;          // 各モードのSIF振幅 [ΔKI, ΔKII, ΔKIII]
    std::vector<double> delta_K_bar;      // ΔK_bar = ΔKI + |ΔKII| + |ΔKIII|
    std::vector<double> delta_K_eq;       // 等価SIF振幅
    double alfa_1 = 1.115, alfa_2 = 1.0;  // ΔK_eq の計算に使用する係数(Richard の式)
    
public:
    // constructon
    StressIntensityFactorRange(const CrackDataSetting_3D& crackData, const int mode, const int domain_types)
    : delta_K(mode * domain_types *crackData.getJintegralPoints(), 0.0),
      delta_K_bar(mode * domain_types * crackData.getJintegralPoints(), 0.0),
      delta_K_eq(domain_types * crackData.getJintegralPoints(), 0.0) {}

    // member fucntion
    void calcDeltaK(const std::vector<double>& K_max, const CrackDataSetting_3D& crackData, const ParisLawCrack& paris);   // モード別SIF振幅の計算
    void calcDeltaKbar(const CrackDataSetting_3D& crackData);                                                              // ΔK_bar = ΔKI + |ΔKII| + |ΔKIII| の計算
    void calcDeltaKeq(const CrackDataSetting_3D& crackData, const ParisLawCrack& paris);                                   // 等価SIF振幅の計算

    // getter
    const std::vector<double>& getDeltaK() const { return delta_K; }
    const std::vector<double>& getDeltaKbar() const { return delta_K_bar; }
    const std::vector<double>& getDeltaKeq() const { return delta_K_eq; }
    
    // setter
    void setDeltaK(const std::vector<double>& dK) { delta_K = dK; }
};


// き裂進展速度，き裂進展サイクル数，き裂進展量の計算クラス
class ParisLawCrack {
private:
    int crack_mode;                       // き裂モード (0: mode I, 1: mixed mode)
    double C;                             // Paris則定数C
    double m;                             // Paris則定数m
    double delta_a_max;                   // 最大き裂進展長さ
    int weighted_moving_average;          // 重み付き移動平均化の有無 (0: なし, 1: 移動平均, 2: 重み付き移動平均)
    int number_moving_average;            // 移動平均を行うデータ数
    std::vector<double> da_dN;            // き裂進展速度
    std::vector<double> da_dN_component;  // き裂進展速度の成分ごと
    std::vector<double> da_dN_max;        // 最大き裂進展速度
    // int selected_Nk;                      // 繰り返し荷重サイクル数の選択 (0:進展前N-a曲線の傾き, 1:進展前後N-a曲線の傾きの平均)
    std::vector<double> Nk;               // き裂進展サイクル数

public:
    // constructor
    ParisLawCrack(const CrackDataSetting_3D& crackData, information *info, const int domain_types)
        : crack_mode(0), 
          C(0.0), 
          m(0.0), 
          delta_a_max(0.0), 
          weighted_moving_average(0),
          number_moving_average(0), 
          da_dN(crackData.getJintegralPoints() * domain_types, 0.0),
          da_dN_component(crackData.getJintegralPoints() * domain_types * info->DIMENSION, 0.0),
          da_dN_max(domain_types, 0.0),
        //   selecteds_Nk(0),
          Nk(domain_types, 0.0){}

    // member fucntion
    void readParameters(const std::string& filename);                                                                                         // ファイルからパラメータを読み込む関数
    void calcCrackGrowthRate(const CrackDataSetting_3D& crackData, const StressIntensityFactorRange& deltaK);                                 // き裂進展速度の計算 da/dN = C(ΔK)^m
    void calcWeightedMovingAverageCrackGrowthRate(const CrackDataSetting_3D& crackData);                                                      // き裂進展速度の重み付き移動平均化
    void calcCrackGrowthRateComponent(const CrackDataSetting_3D& crackData, const CrackDirection& crackExtensionVector, information *info);   // き裂進展速度の計算(da/dN)の成分ごと
    void calcMaxCrackGrowthRate(const CrackDataSetting_3D& crackData);                                                                        // 最大き裂進展速度(da/dN)maxの計算
    void calcCrackGrowthCycle(const CrackDataSetting_3D& crackData);                                                                          // き裂進展サイクル数Nkの計算
    void readAfterCrackGrowthCycle(const CrackDataSetting_3D& crackData, std::vector<double>& Nk_after);                                      // 進展後のき裂進展サイクル数の読み込み
    void calcAverageCrackGrowthCycle(const CrackDataSetting_3D& crackData, std::vector<double>& Nk_after);                                    // き裂進展サイクル数の進展前後の平均化
    void writeCrackGrowthCycle(const CrackDataSetting_3D& crackData);                                                                       // き裂進展サイクル数の書き込み

    // getter
    int getCrackMode() const { return crack_mode; }
    double getC() const { return C; }
    double getM() const { return m; }
    double getMaxGrowth() const { return delta_a_max; }
    const std::vector<double>& getCrackGrowthRate() const { return da_dN; }
    const std::vector<double>& getMaxCrackGrowthRate() const { return da_dN_max; }
    const std::vector<double>& getCrackGrowthCycle() const { return Nk; }
    const std::vector<double>& getda_dN_component() const { return da_dN_component; }

    // setter
    void setCrackMode(int val) { crack_mode = val; }    
    void setC(double val) { C = val; }
    void setM(double val) { m = val; }
    void setWeightedMovingAverage(int val) { weighted_moving_average = val; }
    void setNumberMovingAverage(int val) { number_moving_average = val; }
    void setMaxGrowth(double val) { delta_a_max = val; }
};


// き裂進展方向,ベクトルの計算クラス
class CrackDirection {
private:
    double kinkAngle;                          // キンク角
    double twistAngle;                         // ねじれ角
    std::vector<double> CrackGrowthDirection;  // き裂進展方向
    std::vector<double> CrackGrowthVector;     // き裂進展ベクトル

public:
    // constructor
    CrackDirection(const CrackDataSetting_3D& crackData, information *info, const int domain_types) 
    : CrackGrowthDirection(crackData.getJintegralPoints() * info->DIMENSION, 0.0),  
      CrackGrowthVector(crackData.getJintegralPoints()  * domain_types * info->DIMENSION, 0.0){}  

    // member fucntion
    // キンク角の計算
    // ねじれ角の計算
    void calcCrackGrowthDirection(const CrackDataSetting_3D& crackData, const IIM_3D& IIM3D, const ParisLawCrack& paris, information *info);  // き裂進展方向の計算(キンク角は未考慮，き裂前縁ベクトルをそのまま格納)
    void calcCrackGrowthVector(const CrackDataSetting_3D& crackData, const ParisLawCrack& paris, information *info);                          // き裂進展ベクトルの計算

    // getter
    const std::vector<double>& getCrackGrowthDirection() const { return CrackGrowthDirection; }
    const std::vector<double>& getCrackGrowthVector() const { return CrackGrowthVector; }

    // setter
    void setCrackGrowthDirection(const std::vector<double>& value) { CrackGrowthDirection = value; }
};


// き裂前縁上の離散点を管理するクラス
class CrackFrontManager {
private:
    std::vector<double> before_crack_front_points;    // 現在のき裂前縁点群
    std::vector<double> after_crack_front_points;     // 現在のき裂前縁点群

public:
    // constructor
    CrackFrontManager(const CrackDataSetting_3D& crackData, information *info, const int domain_types)
        : before_crack_front_points(crackData.getJintegralPoints() * info->DIMENSION, 0.0),
          after_crack_front_points(domain_types * crackData.getJintegralPoints() * info->DIMENSION, 0.0) {}

    // member function
    void getbeforeCrackFrontPoints(const IIM_3D& IIM3D);                                                                                 // 進展前き裂前縁の幾何情報の取得
    void calcaferCrackFrontPoints(const CrackDirection& crackExtensionVector, const CrackDataSetting_3D& crackData, information *info);  // 進展後き裂前縁の幾何情報の計算

    // getter
    const std::vector<double>& getbeforeCrackFrontPoints() const { return before_crack_front_points; }
    const std::vector<double>& getafterCrackFrontPoints() const { return after_crack_front_points; }

    // setter
    void setbeforeCrackFrontPoints(const std::vector<double>& value) { before_crack_front_points = value; }
    void setafterCrackFrontPoints(const std::vector<double>& value) { after_crack_front_points = value; }
};


// き裂前縁離散点から最小二乗法により制御点座標を計算するクラス
class LeastSquaresApproximation {
private:
    std::vector<double> crack_front_cp_coordinate;  // き裂前縁離散点から導かれる制御点座標 (cp * dim)
    std::vector<int> LSA_patchList;                 // 最小二乗法を行うためのパッチリスト
    std::vector<int> cp_list;                       // 最小二乗法計算に用いる全き裂前縁パッチの制御点IDのセット
    std::vector<int> cp_size;                       // 最小二乗法計算に用いる各パッチの制御点数
    std::vector<int> Jcp_size;                      // 最小二乗法計算に用いる各パッチのき裂前縁点数
    int selected_domain;                            // き裂前縁データ取得時の積分領域の選択を反映
    bool use_total;

     // 🔽 追加（全体一括処理用）
    std::map<int, int> cp_global_index_map;  // 制御点ID → 行列列番号への対応
    struct JcpInfo {
        double eta;
        double eta_total;
        int patchID;
        int elementID;
    };
    std::vector<JcpInfo> Jcp_info_list;    // 全J積分点の情報（η座標・パッチ・要素）
    std::unordered_map<int, int> patch_cp_start_index; //各パッチが cp_list の中で始まるインデックスを記録するマップ
    std::unordered_map<int, int> cp_local_index_map;  // 全体の制御点ID → パッチ内のローカル列番号


public:
    //constructor
    LeastSquaresApproximation(const CrackDataSetting_3D& crackData, information *info)
        : crack_front_cp_coordinate(crackData.getJintegralPoints() * info->DIMENSION, 0.0),  
          LSA_patchList(crackData.getnum_edge_crack_patch(), 0),
          cp_list(),
          cp_size(crackData.getnum_edge_crack_patch(), 0),
          Jcp_size(crackData.getnum_edge_crack_patch(), 0),
          selected_domain(0)
          {} 
        
     // 🔽 追加インターフェース
    void preLeastSquaresApproximation(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information* info, bool use_global);
    void preLSA_Patchwise(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information* info);
    void preLSA_Global(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information* info);
    void DoLeastSquaresApproximation_Global(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, IIM_3D& IIM_3D, information *info, int selected_domain);  
    void calcNTN_inverse_NT_Global(std::vector<double>& N_matrix, std::vector<double>& NTN_inverse_NT_matrix);
    void calcControlPointCoordinate_Global(const CrackFrontManager& crackFrontCoord, vector<double>& NTN_inverse_NT_matrix, information *info, int selected_domain);
    void correctControlPointCoordinate_Global(information *info);
    void DoLeastSquaresApproximation(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, IIM_3D& IIM3D, information *info, int selected_domain, bool use_total); 
    void calcInverseMatrix_byEigen_Global(std::vector<double>& N_matrix, std::vector<double>& NTN_inverse_NT_matrix, int total_cp, int total_Jcp);
    void setUseTotal(bool flag) { use_total = flag; };
    void copyGlobalCPResultToAll(); // 引数なし
    std::vector<double> generate_Open_Knot_Vector(int num_cp, int degree);
    void calcShapeFunction1D_Global(double eta, const std::vector<double>& knots, int degree, std::vector<double>& N);



    // member function
    //void preLeastSquaresApproximation(const CrackDataSetting_3D& crackData, IIM_3D& IIM3D, information *info);                       // 最小二乗法による制御点座標の計算のための準備
    void DoLeastSquaresApproximation_Patchwise(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord,             
                                     IIM_3D& IIM3D, information *info, int selected_domain);                                         // 最小二乗法による制御点座標
    void calcNTN_inverse_Patchwise(std::vector<double>& N_Matrix, std::vector<double>& NTN_inverse_NT_matrix, int edge_patch);              // (N^T N)^(-1)N^Tの計算 
    void calcInverseMatrix_byEigen_Patchwise(std::vector<double>& NTN_matrix, std::vector<double>& NTN_inverse_matrix, int edge_patch);        // N^T N の逆行列の計算                                 
    void calcControlPointCoordinate_Patchwise(const CrackDataSetting_3D& crackData, const CrackFrontManager& crackFrontCoord, 
                                    vector<double>& NTN_inverse_NT_matrix, int edge_patch, information *info, int selected_domain);  // (N^T N)^(-1)N^T * Jcp 計算
    void deleteDuplicateCP(const CrackDataSetting_3D& crackData, information *info);                                                 // 重複する制御点IDの制御点座標の平均化と削除
    void correctControlPointCoordinate_Patchwise(information *info);                                                                           // パッチ境界上の制御点座標の補正(元の座標からずれるべきでない制御点座標の補正)
    void ReadCPSettingLSA(const std::string& filename, vector<int> &correction_patch_set, 
                          vector<int> &edge_parameter_set, vector<int> &direction_set);                                              // 制御点座標の補正のための制御点情報の読み込み

    // getter
    const std::vector<double>& getControlPointCoordinate() const { return crack_front_cp_coordinate; }
    const std::vector<int>& getLSA_patchList() const { return LSA_patchList; }
    const std::vector<int>& getcp_list() const { return cp_list; }
    const std::vector<int>& getcp_size() const { return cp_size; }
    int getcp_size(int edge_patch) const { return cp_size[edge_patch]; }
    int getJcp_size(int edge_patch) const { return Jcp_size[edge_patch]; }
    const std::map<int, int>& getcp_global_index_map() const {return cp_global_index_map;}

    // setter
    void setControlPointCoordinate(const std::vector<double>& value) { crack_front_cp_coordinate = value; }
};


// き裂進展解析の入出力クラス
class FatigueCrackIO {
private:
    int select_domain;  // 進展に基づく積分領域の選択

public:
    // constructor
    FatigueCrackIO(int domainval)
    : select_domain(domainval){};  

    // member function
    static void crackExtensionAnalysis(information *info, CrackDataSetting_3D *crackData, IIM_3D *IIM3D);  // き裂進展解析の実行
    void read_select_domain(const std::string& filename, const CrackDataSetting_3D& crackData);                                                  // 進展に基づく積分領域の選択
    void writeCrackFrontData(const CrackFrontManager& crackFrontCoord, information *info);                 // き裂前縁離散点の出力  ステップ情報も追加したい
    void writeCrackFrontCPData(const LeastSquaresApproximation& calc_cp_LSA, information *info);                  // 進展前後制御点座標の出力
    

    // getter
    int getSelectDomain() const { return select_domain; }

    // setter
    void setSelectDomain(int domainval) { select_domain = domainval; }
};


// リメッシングクラス
class Remeshing {
private:
    vector<int> patch_around_crack;  // き裂前縁パッチ
    vector<int> patch_behind_crack;  // き裂前縁パッチの後方のパッチ
    vector<int> patch_remesh;        // リメッシングを行うパッチ
    vector<int> remesh_direction;    // リメッシュ方向
    vector<int> division_ele_n;      // リメッシュ後の制御点数
    vector<double> cp_move_vector;   // 各制御点の移動率

public:
    // constructor
    Remeshing(const LeastSquaresApproximation& calc_cp_LSA, information *info) 
    : patch_around_crack(),
      patch_behind_crack(),
      patch_remesh(),
      remesh_direction(), 
      division_ele_n(),
      cp_move_vector(calc_cp_LSA.getcp_list().size() * info->DIMENSION, 0.0) 
    {}

    // member function
    void readRemeshPatch(const std::string& filename);                                                                 // リメッシュデータの読み込み
    void DoRemeshing(const CrackDataSetting_3D& crackData, LeastSquaresApproximation calc_cp_LSA, information *info);  // リメッシュ処理の実行
    void calcCPmoveVector(LeastSquaresApproximation calc_cp_LSA, information *info);                                   // き裂前縁を構成する各制御点の移動率の計算
    void RedefinePatchGeometry(const CrackDataSetting_3D& crackData, LeastSquaresApproximation calc_cp_LSA,
                               information *info, vector<double>& u_vec, vector<int>& patch_around_crack,
                               int offset, vector<vector<double>>& new_coord, vector<vector<double>>& new_w);          // パッチのジオメトリの再定義
    void outputRemeshingData(information *info, vector<vector<double>>& new_coord);                                    // リメッシング後の可視化データの出力
    void outputEachPatchData(information *info);                                                                       // リメッシング後の各パッチの可視化データの出力
    void makeInputfileNextStep(information *info, vector<vector<double>>& new_coord);                                  // リメッシング後のインプットファイルの出力
};