#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <memory>

using std::cout;
using std::endl;
using std::vector;
using std::string;

// external variables
extern double E;                  // Young's modulus
extern double nu;                 // Poisson's ratio
extern int Total_mesh;            // total number of meshes
extern int MAX_NO_CP_ON_ELEMENT;  // maximum number of control points (support basis function) on element 
extern int MAX_KIEL_SIZE;         // element stiffness matrix
extern int D_MATRIX_SIZE;         // dispalcement matrix size
extern int N_STRAIN;              // 
extern int N_STRESS;              // 
extern int GP_ON_ELEMENT;         // 要素内のガウス点数

namespace {
    const int MODEL_SIZE_PARA_2D = 0;         // 0 : き裂上下面パッチを積分領域とする, 1 : き裂上面パッチのみを積分領域とする
    const int PARALLEL_CALC_FRACTURE_2D = 0;  // 2D破壊力学計算部分の並列計算の有無 (0 : 並列化なし, 1 : OpenMP)
}

// class
class PQSetting {
public:
	vector<double> Disp_grad;
	vector<double> Strain;
	vector<double> Stress;
	double StrainEnergyDensity;
	double StrainEnergyDensity_aux;
    // constructor
    PQSetting(int dimension)
        : Disp_grad(dimension *dimension, 0.0)
        , Strain(D_MATRIX_SIZE, 0.0)
        , Stress(D_MATRIX_SIZE, 0.0)
        , StrainEnergyDensity(0.0)
        , StrainEnergyDensity_aux(0.0) {};
    // destructor
    ~PQSetting() {};
};

class CrackDataSetting {
public:
    vector<double> crackOriginCoordinates;
    vector<double> crackTipCoordinates;
    vector<double> virtualCrackExtension_cp;
    double delta;
    double T[MAX_DIMENSION][MAX_DIMENSION];
    // constructor
    CrackDataSetting(int totalCPToMesh, int dimension)
        : crackOriginCoordinates(dimension, 0.0)
        , crackTipCoordinates(dimension, 0.0)
        , virtualCrackExtension_cp(totalCPToMesh * dimension, 0.0) {};
    // destructor
    ~CrackDataSetting() {};
};

class IIM_2D {
    public:
        void IIMComputation(CrackDataSetting *crackData, information *info);
        void GaussianElimination(vector<double> &A, int num);
        void AppendDataToCSV(string fileName, double data);
        void ReadCrackData(CrackDataSetting *crackData, information *info);
        void MakeUnitBasisVector(CrackDataSetting *crackData, information *info);
        double CalculateEMT(PQSetting *act, vector<double> &qFuncGradient, bool isAuxiliary, information *info);
        double CalculateEMT(PQSetting *act, PQSetting *aux, vector<double> &qFuncGradient, information *info);
        double CalculateSIF(CrackDataSetting *crackData, double *value, double modelSize, bool isIIM);
};

class Auxiliary : public IIM_2D {
    public:
        void MakeAuxiliary(PQSetting *aux, PQSetting *actOverlay, int e, int GPNum, vector<double> &crackTipCoordinates,double T[MAX_DIMENSION][MAX_DIMENSION], bool isMode1, information *info);
        void CoordTransform(double T[MAX_DIMENSION][MAX_DIMENSION], vector<double> &S, vector<double> &result, information *info);
        void MakeDispGradient_aux1(PQSetting *aux1, double r, double rad, double mu);
        void MakeDispGradient_aux2(PQSetting *aux2, double r, double rad, double mu);
        void MakeStress_aux1(PQSetting *aux1, double r, double rad);
        void MakeStress_aux2(PQSetting *aux2, double r, double rad);
};

class PhysicalQuantity {
    public:
        void CheckGlobalInLocalCoordinates(PQSetting *globalInLoc, int elementID, int GPNum, information *info);
        void MakeDispGradientAndStrain(PQSetting *globalInLoc, int localElemID, int globalElemID, int GP_num, information *info);
        void Make_B_Matrix_anypoint(int elementID, vector<double> &B, double *localCoords, information *info);
        void Make_Disp_Gradient_Matrix(int elementID, vector<double> &Disp_grad, double *localCoords, information *info);
        void MakeDispGradientAtGP(PQSetting *act, int elementID, int GPNum, information *info);
        void MakeStrainAtGP(PQSetting *act, int elementID, int GPNum, information *info);
        void MakeStressAtGP(PQSetting *act, information *info);
        void MakeStrainEnergyDensityAtGP(PQSetting *act);
        void Overlay(vector<double> &value1, vector<double> &value2, vector<double> &value3, int num);
};
