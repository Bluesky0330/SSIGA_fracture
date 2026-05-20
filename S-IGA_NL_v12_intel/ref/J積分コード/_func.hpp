#ifndef _FUNC
#define _FUNC

#include <vector>
#include "_class.hpp"

using namespace std;


// memory allocation
void Allocation(const int num, information *info);
void Global_var(const int num, information *info);
// Get input file data
void Get_Input_1(int tm, char **argv, information *info);
void Get_Input_2(int tm, char **argv, information *info);
void Make_INC(information *info);
// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info);
void Setting_Dist_Load_3D(int mesh_n, int iPatch, int iCoord, int jCoord, double val_Coord, double *iRange_Coord, double *jRange_Coord, int type_load, double *iCoeff_Dist_Load, double *jCoeff_Dist_Load, information *info);
int SearchForElement_2D(int mesh_n, int iPatch, int iX, int iY, information *info);
int SearchForElement_3D(int mesh_n, int iPatch, int iX, int iY, int iZ, information *info);
// for IGA
void Preprocessing_IGA(information *info);
// for S_IGA
void Preprocessing_S_IGA(information *info);
int ele_check(int patch_n, double *para_coord, information *info);
void tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info);
void sort(int total, int *element_n_point);

// Preprocessing
void Preprocessing(int m, int e, information *info);
void Make_Gauss_Coordinate(int m, int e, information *info);
void Make_dSF(int m, int e, double *dsf, information *info);
void Make_Jac(int m, int e, double *a_matrix, double *dsf, information *info);
void Make_B_Matrix(int m, int e, double *a_matrix, double *dsf, information *info);
void Make_gauss_array(int select_GP, information *info);
void Gauss_point(information *info, int temp_GP_1D, double *w_ptr, double *point_ptr, double *w_1D_ptr, double *point_1D_ptr);
// K matrix
void Make_D_Matrix(information *info);
void Make_K_Whole_Ptr_Col(information *info, int mode_select);
double Make_B_Matrix_anypoint(int El_No, double *B, double *Local_coord, information *info);
double Make_Jac_anypoint(int El_No, double *Local_coord, information *info);
//Gaussian Elimination
bool GaussianElimination2(double *sol, const double *diff, double *J, int dimension);
// tool
double InverseMatrix_2x2(double M[2][2]);
double InverseMatrix_2x2(double *M);
double InverseMatrix_3x3(double M[3][3]);
double InverseMatrix_3x3(double *M);
long pow_int(int val, int n);
int duplicate(int total_n, int *array);
// Shape Function
double Shape_func(int I_No, double *Local_coord, int El_No, information *info);
void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info);
void ShapeFunction1D(double *Position_Data_param, int j, int e, double *Shape, double *dShape, int select_deriv, information *info);
double dShape_func(int I_No, int xez, double *Local_coord, int El_No, information *info);
void NURBS_deriv_2D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, information *info);
void NURBS_deriv_3D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, double *dShape_func3, information *info);
double dShapeFunc_from_paren(int j, int e, information *info);
void shape_function_1D(double knot, int dim, int e, double *Shape, double *dShape, int select_deriv, information *info);
void shape_and_dshape(double *R, double *dR, double *Local_coord, int e, information *info);
int calc_patch_parameter_coord(double *physical_coord, int patch_num, double *out_para_coord, information *info);
// Newton-Raphson
double &BasisFunc(double *knot_vec, const int &knot_index, const int &order, const double &xi, double *output, double *d_output);
double &rBasisFunc(double *knot_vec, const int &knot_index, const int &order, const double &xi, double *output, double *d_output);
double &lBasisFunc(double *knot_vec, const int &knot_index, const int &cntl_p_n, const int &order, const double &xi, double *output, double *d_output);
double NURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                     double *cntl_px, double *cntl_py,
                     int cntl_p_n_xi, int cntl_p_n_eta,
                     double *weight, int order_xi, int order_eta,
                     double xi, double eta,
                     double *output_x, double *output_y,
                     double *output_dxi_x, double *output_deta_x,
                     double *output_dxi_y, double *output_deta_y);
double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y);
double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y);
double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y);
double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y);
double rrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double lrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double rlrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double rrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double llrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double lrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double rllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
double lllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z);
int Calc_xi_eta(double px, double py,
                double *input_knot_vec_xi, double *input_knot_vec_eta,
                int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
                double *output_xi, double *output_eta, int global_patch_num, information *info);
int Calc_xi_eta_zeta(double px, double py, double pz,
                     double *input_knot_vec_xi, double *input_knot_vec_eta, double *input_knot_vec_zeta,
                     int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                     int order_xi, int order_eta, int order_zeta,
                     double *output_xi, double *output_eta, double *output_zeta, int global_patch_num, information *info);
// Postprocessing
void Substitute_Displacement(information *info, char **argv);
void Substitute_EOI(information *info, char **argv);

// paraview
void Search_ele_point_2D(const int &xi, const int &eta, const int &e_x_max, const int &e_y_max, int *p_x, int *p_y, int *e_x, int *e_y, int *point);
void Search_ele_point_3D(const int &xi, const int &eta, const int &zeta, const int &e_x_max, const int &e_y_max, const int &e_z_max, int *p_x, int *p_y, int *p_z, int *e_x, int *e_y, int *e_z, int *point);

// search overlapping element
void searchOverlappingEle(information *info);
int setStartEle(const int e, information *info);
bool checkInner(target_domain &td_base, target_domain &td, information *info);
bool checkInnerCenter(target_domain &td_base, int e_center, information *info);
bool checkInnerCenterEachPatch(target_domain &td_base, target_domain &td, information *info);
bool checkPosition(target_domain &td_base, double *check_para, information *info);
void addInnerEle(vector<int> &eoi_e, vector<int> &enoi_e, target_domain &td_base, information *info);
template <typename T>
void sortUniqueErase(std::vector<T> &vec);
void physical_coord(const int e, double *para, double *out_coord, information *info);
void physical_coord_R(const int e, double *para, double *R, double *out_coord, information *info);
void makeChildPosition(int *child_position, int i, int *parent_position, information *info);
bool search2D(target_domain &td_base, target_domain &td, information *info);
void setEdge2D(int edge_num, int *target_axis, double *para, target_domain &temp);
bool search3D(target_domain &td_base, target_domain &td, information *info);
void setEdge3D(int edge_or_face, int edge_or_face_num, int *target_axis, double *para, target_domain &temp);
void elementConnectivity(information *info, vector<vector<int>> &ecn);

#endif