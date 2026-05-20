#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <memory>

using std::cout;
using std::endl;
using std::string;
using std::vector;

// external variables
extern double E;                  // Young's modulus
extern double nu;                 // Poisson's ratio
extern int Total_mesh;            // total number of meshes
extern int MAX_NO_CP_ON_ELEMENT;  // maximum number of control points (support basis function) on element 
extern int MAX_KIEL_SIZE;         // element stiffness matrix
extern int D_MATRIX_SIZE;         // dispalacement matrix size
extern int N_STRAIN;              // 
extern int N_STRESS;              // 
extern int GP_ON_ELEMENT;         // number of Gauss points on element


namespace {
    // const int MODEL_SIZE_PARA_3D = 0;            // き裂面の積分領域を設定 (0:上下両面, 1:どちらかのみ)
    const int MODE = 3;                          // number of fracture mode (1:mode 1 only, 2:mode 1 and 2, 3:all mode)
    // const int DOMAIN_TYPES = 5;                  // number of domain size type(1:only Whole singular area, 2:Whole singular & Whole - 1, 3:Whole & Whole -1 & Whole - 2, 3:Whole & Whole -1 & Whole - 2 & Whole - 3, ...) ξ方向要素数分だけ選択可能
    const int DOMAIN_SIDES = 2;                  // 積分領域として設定する幅(η)方向要素数, 計算点を挟んで左右に1つずつ領域を確保するため 2 とする
    // const int DOMAIN_AREA_DIVIDES = 10;          // ΔAの計算に置けるき裂前縁方向のき裂長さの計算に使用する領域分割数 
    const int MOVE_DOMAIN = 1;                   // 積分領域がパッチ境界をまたぐ場合の領域移動 (0 : 移動なし, 1 : 移動あり)
    const int OpenMP = 0;                        // OpenMP の使用 (0 : 使用しない, 1 : 使用する)
}


// class
// 物理量を格納するクラス
class PQSetting_3D {
private:
    vector<double> dispGrad1;
    vector<double> dispGrad2;
    vector<double> dispGrad3;
    vector<double> strain;
    vector<double> stress;
    double strainEnergyDensity = 0.0;
    double strainEnergyDensityAux = 0.0;
public:
    // constructor
    PQSetting_3D(int dimension) :
    dispGrad1(dimension, 0.0), 
    dispGrad2(dimension, 0.0),
    dispGrad3(dimension, 0.0),
    strain(D_MATRIX_SIZE, 0.0), 
    stress(D_MATRIX_SIZE, 0.0) {};

    // getter
    const vector<double> &getDispGrad1() const { return dispGrad1; }
    const vector<double> &getDispGrad2() const { return dispGrad2; }
    const vector<double> &getDispGrad3() const { return dispGrad3; }
    const vector<double> &getStrain() const { return strain; }
    const vector<double> &getStress() const { return stress; }
    double getStrainEnergyDensity() const { return strainEnergyDensity; }
    double getStrainEnergyDensityAux() const { return strainEnergyDensityAux; }

    // setter
    void setDispGrad1(const vector<double> &dispGrad1) { PQSetting_3D::dispGrad1 = dispGrad1; }
    void setDispGrad2(const vector<double> &dispGrad2) { PQSetting_3D::dispGrad2 = dispGrad2; }
    void setDispGrad3(const vector<double> &dispGrad3) { PQSetting_3D::dispGrad3 = dispGrad3; }
    void setStrain(const vector<double> &strain) { PQSetting_3D::strain = strain; }
    void setStress(const vector<double> &stress) { PQSetting_3D::stress = stress; }
    void setStrainEnergyDensity(double strainEnergyDensity) { PQSetting_3D::strainEnergyDensity = strainEnergyDensity; }
    void setStrainEnergyDensityAux(double strainEnergyDensityAux) { PQSetting_3D::strainEnergyDensityAux = strainEnergyDensityAux; }
};


// き裂に関するデータを設定するクラス
class CrackDataSetting_3D {
private:
    vector<int> integralDomainDivisions;           // ξ, η, ζ 方向特異パッチ要素数 (η は1とする)  
    // vector<int> integralDomainDivisions{5, 1, 3};   // ξ, η, ζ 方向特異パッチ要素数 (η は1とする)  singular_suna
    // vector<int> integralDomainDivisions{8, 1, 8};  // ξ, η, ζ 方向特異パッチ要素数 (η は1とする)  singular_oishi
    // vector<int> integralDomainDivisions{2, 1, 2};  // ξ, η, ζ 方向特異パッチ要素数 (η は1とする)  debug
    int domainWidth;                                  // 1 / domainWidth がき裂前縁(η)方向積分領域幅
    int JintegralPoints;                              // き裂前縁に沿って SIF を計算する点の数を決定
    double integralDomainNaturalWidth;                // η方向積分領域幅
    int model_size_para_3D;                           // き裂面の積分領域を設定 (0:上下両面, 1:どちらかのみ)
    int domain_types;                                 // number of domain size type(1:only Whole singular area, 2:Whole singular & Whole - 1, 3:Whole & Whole -1 & Whole - 2, 3:Whole & Whole -1 & Whole - 2 & Whole - 3, ...) ξ方向要素数分だけ選択可能
    int domain_area_divisions;                        // ΔAの計算に置けるき裂前縁方向のき裂長さの計算に使用する領域分割数
    int patchType;                                    // き裂先端のパッチの種類 ( 0 : 縮退・特異パッチを使用したき裂先端パッチを極座標形式とした場合，1 : 通常パッチを使用 )
    int local_init_patch;                             // Total mesh - 1 の全グローバルパッチ数（ローカル初期パッチIDの特定に用いる）
    int num_edge_crack_patch;                         // き裂前縁(η)方向のパッチ数
    int num_singular_patch;                           // 一回のSIF計算に関与する特異パッチ数
    vector<int> patchList;                            // 計算領域として使用するパッチID を格納 

public:
    // constructor
    // CrackDataSetting_3D(){};

    CrackDataSetting_3D(int dimension)
    : integralDomainDivisions(dimension, 1){};

    // member function
    void readCrackData();                             // 破壊力学解析用のインプットファイルの読みとり
    void setCrackTipPatchElement(information *info);  // き裂先端のパッチの要素数の設定

    // getter
    const vector<int> &getIntegralDomainDivisions() const { return integralDomainDivisions; }
    int getJintegralPoints() const { return JintegralPoints; }
    double getIntegralDomainNaturalWidth() const { return integralDomainNaturalWidth; }
    int getmodel_size_para_3D() const { return model_size_para_3D; }
    int getdomain_types() const { return domain_types; }
    int getdomain_area_divisions() const { return domain_area_divisions; }
    int getpatchType() const { return patchType; }
    int getnum_edge_crack_patch() const { return num_edge_crack_patch; }
    int getnum_singular_patch() const { return num_singular_patch; }
    const std::vector<int>& getpatchList() const { return patchList; }
    

    //setter
    void setJintegralPoints(const int JintegralPoints) { CrackDataSetting_3D::JintegralPoints = JintegralPoints; }
    void setIntegralDomainNaturalWidth(const double integralDomainNaturalWidth) { CrackDataSetting_3D::integralDomainNaturalWidth = integralDomainNaturalWidth; }
    // void setlocal_init_patch(const int local_init_patch) { CrackDataSetting_3D::local_init_patch = local_init_patch; }
    void setpatchType(const int patchType) { CrackDataSetting_3D::patchType = patchType; }
    void setnum_edge_crack_patch(const int num_edge_crack_patch) { CrackDataSetting_3D::num_edge_crack_patch = num_edge_crack_patch; }
    void setnum_singular_patch(const int num_singular_patch) { CrackDataSetting_3D::num_singular_patch = num_singular_patch; }
    void setpatchList(const vector<int>& patchList) { this->patchList = patchList; }
};


// 座標変換や行列演算などの数学的な操作を行うメソッドを含むクラス
class MathFormula {
public:
    // member function
    vector<double> convertNaturalCoordinatesToTilde(const int patchID, const int elementID, const vector<double> &naturalCoords, information *info);               //
    vector<double> crossProduct3x3(const vector<double> &a, const vector<double> &b, information *info);                                                           //
    vector<double> tensorTransformation(const vector<double> &unitLocal, const vector<double> &unitGlobal, const vector<double> &inputTensor, information *info);  //
    vector<double> overlay(const vector<double> &value1, const vector<double> &value2, const int num);                                                             //
};


// 補助場の計算に関するメソッドを含むクラス
class AuxiliaryField {
public:
    // member function
    vector<double> makeDisplacementGradientAux1(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info);
    vector<double> makeDisplacementGradientAux2(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info);
    vector<double> makeDisplacementGradientAux3(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info);
    void setDisplacementGradientAux(const vector<double> &tempDispGradient, PQSetting_3D *aux, information *info);
    void makeStrainAux(const vector<double> &tempDispGradient, PQSetting_3D *aux, information *info);
};


// q 関数の勾配を計算するクラス
class QGradientFunction {
private:
    vector<double> qGradient1;
    vector<double> qGradient2;
    vector<double> qGradient3;

public:
    // constructor
    QGradientFunction(int dimension, const CrackDataSetting_3D& crackData) 
        : qGradient1(crackData.getIntegralDomainDivisions()[0] * dimension, 0.0), 
          qGradient2(crackData.getIntegralDomainDivisions()[0] * dimension, 0.0), 
          qGradient3(crackData.getIntegralDomainDivisions()[0] * dimension, 0.0) {};

    // getter
    const vector<double> &getQGradient1() const { return qGradient1; }
    const vector<double> &getQGradient2() const { return qGradient2; }
    const vector<double> &getQGradient3() const { return qGradient3; }

    // setter
    void setQGradient1(const int index, const double value) { qGradient1[index] = value; }
    void setQGradient2(const int index, const double value) { qGradient2[index] = value; }
    void setQGradient3(const int index, const double value) { qGradient3[index] = value; }

    // member function
    void makeQGradient(int patchID, double crackJintegralNatural, const vector<double> &Jac, const vector<double> &dXiTilde_dXiPrime, const vector<double> &patchGaussPointCoords, CrackDataSetting_3D *crackData, information *info);
};


// J積分および相互積分法（Interaction Integral Method）の3次元実装を行うクラス
class IIM_3D {
private:
    vector<vector<int>> elementID_in_patch;   // 各パッチに属する要素ID
    vector<int> elements_in_patch;            // 各パッチに含まれる要素の数
    vector<double> crackFrontVector;          // 0: binormal, 1: tangent, 2: normal (Jcpごとに上書き)
    vector<double> JcpcrackFrontVector;       // 0: binormal, 1: tangent, 2: normal (Jcp全体の情報)
    vector<double> crackEdgeJcpNatural;       // き裂前縁のJ積分計算点の自然座標(η のみ)
    vector<int> crackEdgeJcpElementID;        // き裂前縁のJ積分計算点の要素ID
    vector<double> calculationPointCoords;    // き裂前縁のJ積分計算点の物理座標(x, y, z)
    //vector<double> crackFrontDegrees;       // き裂前縁角度
    vector<double> JintegralScalar;           // J積分値を格納
    vector<double> JintegralSIF;              // J積分によるSIFを格納
    vector<double> MintegralScalar;           // IIM積分値を格納
    vector<double> MintegralSIF;              // IIMによるSIFを格納
    MathFormula *math_;
public:
    // constructor
    // IIM_3D(int dimension, int local_init_patch, MathFormula *math, CrackDataSetting_3D *crackData)
    IIM_3D(int dimension, int local_init_patch, int domain_types, MathFormula *math, CrackDataSetting_3D *crackData)
    : elementID_in_patch(local_init_patch, vector<int>(0)), 
      elements_in_patch(local_init_patch, 0), 
      crackFrontVector(dimension * dimension, 0.0),
      JcpcrackFrontVector(dimension * dimension * crackData->getJintegralPoints(), 0.0),
      crackEdgeJcpNatural(),
      crackEdgeJcpElementID(),
      calculationPointCoords(dimension * crackData->getJintegralPoints(), 0.0),
      // crackFrontDegrees(crackData->getJintegralPoints(), 0.0),
    //   JintegralScalar(DOMAIN_TYPES * crackData->getJintegralPoints(), 0.0),
      JintegralScalar(domain_types * crackData->getJintegralPoints(), 0.0),
    //   JintegralSIF(DOMAIN_TYPES * crackData->getJintegralPoints(), 0.0),
      JintegralSIF(domain_types * crackData->getJintegralPoints(), 0.0),
    //   MintegralScalar(MODE * DOMAIN_TYPES * crackData->getJintegralPoints(), 0.0),  
      MintegralScalar(MODE * domain_types * crackData->getJintegralPoints(), 0.0),  
    //   MintegralSIF(MODE * DOMAIN_TYPES * crackData->getJintegralPoints(), 0.0),
      MintegralSIF(MODE * domain_types * crackData->getJintegralPoints(), 0.0),
      math_(math) {};

    // getter
    const vector<double>& getJcpCrackFrontVector() const { return JcpcrackFrontVector; }
    const vector<double>& getcrackEdgeJcpNatural() const { return crackEdgeJcpNatural; }
    const vector<int>& getcrackEdgeJcpElementID() const { return crackEdgeJcpElementID; }
    const vector<double>& getCalculationPointCoords() const { return calculationPointCoords; }
    // const vector<double>& getcrackFrontDegrees() const { return crackFrontVector; }
    const vector<double>& getJintegralScalar() const { return JintegralScalar; }
    const vector<double>& getJintegralSIF() const { return JintegralSIF; }
    const vector<double>& getMintegralScalar() const { return MintegralScalar; }
    const vector<double>& getMintegralSIF() const { return MintegralSIF; }
    

    // setter
    void setJcpCrackFrontVector(const int index, const double value) { JcpcrackFrontVector[index] = value; }
    void setCalculationPointCoords(const int index, const double value) { calculationPointCoords[index] = value; }
    // void setcrackFrontDegrees(const int index, const double value) { crackFrontVector[index] = value; }
    void setJintegralScalar(const int index, const double value) { JintegralScalar[index] = value; }
    void setJintegralSIF(const int index, const double value) { JintegralSIF[index] = value; }
    void setMintegralScalar(const int index, const double value) { MintegralScalar[index] = value; }
    void setMintegralSIF(const int index, const double value) { MintegralSIF[index] = value; }


    // member function
    void computeIIM3D(CrackDataSetting_3D *crackData, information *info);                                                                                // SIF計算のメイン関数
    double makeCrackFrontVector(double &crackJintegralNatural, int Jcp,vector<int> &domainPatchID,                                                       
                                vector<double> &Jcenter, CrackDataSetting_3D *crackData, information *info);                                             // SIF計算点におけるき裂前縁情報の計算
    void readCrackFrontData(CrackDataSetting_3D *crackData, information *info);                                                                          // き裂前縁データの読み込み
    void makeJintegralDomain(const int patchID, const double crackJintegralNatural, int &integralElementID,                                              
                             vector<double> &integralPointNatural, vector<double> &integralPointNaturalTilde, vector<double> &integralPointCoords,
                             double &integralJac, QGradientFunction *qGrad, CrackDataSetting_3D *crackData, information *info);                          // 積分領域の設定
    void computeActualField(const int gaussID, const int integralElementID, vector<double> &integralPointNaturalTilde, 
                            PQSetting_3D *actOverlay, information *info);                                                                                // 実場における物理量の計算
    void computeAuxiliaryField(const int mode, vector<double> &JintegralCenter, vector<double> &patchDivPointCoords,                                     
                               AuxiliaryField *auxField, PQSetting_3D *actOverlay, PQSetting_3D *aux, information *info);                                // 補助場における物理量の計算
    bool isCoordinateInElement(const vector<double> &center, const int patchID, const int elementID, information *info);                                 // 任意の点がどの要素内にあるかを判定
    void makeBgradientMatrix(const int elementID, vector<double> &bGradient1, vector<double> &bGradient2, vector<double> &bGradient3,                    
                             vector<double> &localCoords, information *info);                                                                            // Bマトリックスの作成
    void makeOrthogonalBasisVector(const vector<double> &vector1, const vector<double> &vector2, information *info);
    void storeJcpCrackFrontVector(information *info, int Jcp);                                                                                           // 各計算点におけるき裂前縁ベクトルを格納
    void makeDisplacementGradientAtGP(const int elementID, vector<double> &localCoords, PQSetting_3D *act, information *info);                           // ガウス点における変位勾配の計算
    void makeStrainAtGP(const int elementID, vector<double> &localCoords, PQSetting_3D *act, information *info);                                         // ガウス点におけるひずみの計算
    void makeStress(PQSetting_3D *act, information *info);                                                                                               // 応力の計算
    void makeStrainEnergyDensityAtGP(PQSetting_3D *act);                                                                                                 // ガウス点におけるひずみエネルギー密度の計算
    void checkGlobalInLocalCoordinates(const int elementID, const int gaussID, PQSetting_3D *globalInLoc, information *info);                            // グローバル座標系を局所座標系に変換
    void makeDispGradientAndStrain(const int localElemID, const int globalElemID, const int gaussID, PQSetting_3D *globalInLoc, information *info);      // 変位勾配とひずみの計算
    void writeSIFsToCSV(const string &fileName, const int Jcp, const vector<double> &calculatePointCoords, const double degree, 
                        information *info, CrackDataSetting_3D *crackData);                                                                              // SIFをsSCVファイルに書き込み
    // void readSIFsFromCSV(const std::string &fileName, information *info, CrackDataSetting_3D *crackData);                                                // SIFをCSVファイルから読み込み
    void readSIFsFromCSV(const std::string &fileName, CrackDataSetting_3D *crackData);                                                                   // SIFをCSVファイルから読み込み
};
