#ifndef _FUNC
#define _FUNC

#include <vector>
#include <string>
#include "_class.hpp"

using namespace std;

// option file
void Get_Option(const char *s, vector<string> &opt_files);
void Set_Option_Data(information *info);
// memory allocation
void Allocation(const int num, information *info);
void Global_var(const int num, information *info);
// Get input file data
void Get_Input_1(int tm, const char *filename, information *info);
void Get_Input_2(int tm, const char *filename, information *info);
void Make_INC(information *info);
void Make_INC_for_local_geometric(information *info);
void setDistLoad(int current_mesh, int patch, int coord, double target_knot, double *range, int type_load, double *dist_load_coeff, information *info);
void setDistLoad(int current_mesh, int patch, int coord_i, int coord_j, double target_knot, double *range_i, double *range_j, int type_load, double *dist_load_coeff_i, double *dist_load_coeff_j, information *info);
void setDistLoad_infinite_plate_with_hole(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info);
int SearchForElement_2D(int mesh_n, int iPatch, int iX, int iY, information *info);
// search element
void trans_ele_patch_coord(double* xi_patch, const double* xi_elem, int patch_num, int ele_num, information* info);
void geo_trans_ele_patch_coord(double* xi_patch, const double* xi_elem, int patch_num, int ele_num, information* info);
int trans_local_para_to_local_geo_para(int ele_disp, double *para_disp, double *para_geo, information* info);
int trans_local_para_to_global_para(int ele_loc, double *para_loc, double *para_glo, information* info);
int ele_check(int patch_n, double *para_coord, information *info);
int geo_ele_check(int patch_n, double *para_coord, information *info);
void tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info);
void geo_tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info);
void sort(int total, int *element_n_point);
// gauss points
void Make_gauss_array(information *info);
void Gauss_point(information *info, int temp_GP_1D, double *w_ptr, double *point_ptr, double *w_1D_ptr, double *point_1D_ptr);
void Make_Gauss_points(bool isSinglePatch, information *info);
void gp_switch(bool flag, information *info);
// calcurate B
double Make_B_component(int ele, double *para, double *out_B_component, information *info);
double Make_B_component_for_SSIGA(int ele, double *para, double *out_B_component, information *info);
double Make_updated_jacobian(int ele, double *para, information *info);
double Make_B_Linear(int ele, double *para, double *out_BL, information *info);
// K matrix
void Make_D_Matrix(information *info);
void Make_Index_Dof(information *info);
void Make_K_Whole_Ptr_Col(information *info, int mode_select);
void Calc_K(information *info);
void Calc_K_linear_EL(int El_No, double *K_EL, information *info);
void Calc_coupled_K_linear_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info);
void BDBJ(int KIEL_SIZE, double *B, double *D, double J, double *K_EL);
void coupled_BDBJ(int KIEL_SIZE_glo, int KIEL_SIZE_loc, double *B, double *BG, double *D, double J, double *K_EL);
// F vector
void Make_F_Vec(information *info);
void Make_F_Vec_disp_const(information *info);
void Add_Equivalent_Nodal_Force_to_F_Vec(information *info);
// PCG solver
void PCG_Solver(int max_itetarion, double eps, information *info);
void Make_M(int M_mode, double *M, int *M_Ptr, int *M_Col, double *M_diag, int ndof, information *info);
void CG(const int &ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *M_diag, double *right_vec, double *gg, double *dd, double *pp, double *temp_r);
void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int &ndof);
double inner_product(int ndof, double *vec1, double *vec2);
int RowCol_to_icount(int row, int col, information *info);
// GMRES solver
void GMRES(int length, double eps, information *info);
// Gaussian Elimination
void GaussianElimination(double *sol, double *r, double *A, int size);
bool GaussianElimination2(double *sol, const double *diff, double *J, int dimension);
void swap_matrix(double *A, double *r, int row1, int row2, int size, double *r_for_swap);
// LU decomposition
void LU(double *sol, double *r, double *A, int size);
void intel_PARDISO(double *sol, double *rhs, int size, information *info);
void temp_solver(double *sol, double *r, double *M, int row_n, int col_n);
void overdetermined_system(double *sol,double *rhs, double *A, int m, int n);
// tool
double InverseMatrix_2x2(double *M);
double InverseMatrix_2x2(double M[2][2]);
double InverseMatrix_3x3(double *M);
double InverseMatrix_3x3(double M[3][3]);
long pow_int(int val, int n);
// Shape Function
void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info);
void Geo_ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info);
double dShapeFunc_from_paren(int j, int e, information *info);
double Geo_dShapeFunc_from_paren(int j, int e, information *info);
void shape_function_1D(double knot, int dim, int e, double *Shape, information *info);
void shape_function_1D(double knot, int dim, int e, double *Shape, double *dShape, information *info);
void Geo_Shape_Function_1D(double knot, int dim, int e, double *Shape, information *info);
void Geo_Shape_Function_1D(double knot, int dim, int e, double *Shape, double *dShape, information *info);
void shape_and_dshape(double *R, double *Local_coord, int e, information *info);
void shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info);
void geo_shape_and_dshape(double *R, double *Local_coord, int e, information *info);
void geo_shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info);
void Bspline_shape_and_dshape(double *R, double *Local_coord, int e, information *info);
void Bspline_shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info);
int calc_patch_parameter_coord(double *physical_coord, int patch_num, double *out_para_coord, information *info);
int calc_global_patch_parameter_coord(double *local_para_coord, int patch_num, double *out_para_coord, information *info);
int calc_local_patch_parameter_coord(double *global_para_coord, int patch_num, double *out_para_coord, information *info);
bool check_nonfinite(double *vec, int size);
// Postprocessing
void Make_Displacement(information *info);
void Substitute_Displacement(information *info, char **argv);
void Output_Displacement(information *info);
// output
void Calc_on_Element_Vertex(information *info);
void Calc_on_Gauss_Point(information *info);
// output SVG K matrix
void K_output_svg(information *info);
// paraview
void output_global_parameters(information *info);
void output_for_paraview_timestep(information *info, bool isGeometryOnly, double factor);
void Make_connectivity(information *info);
void Search_ele_point_2D(const int &xi, const int &eta, const int &e_x_max, const int &e_y_max, int *p_x, int *p_y, int *e_x, int *e_y, int *point);
void Search_ele_point_3D(const int &xi, const int &eta, const int &zeta, const int &e_x_max, const int &e_y_max, const int &e_z_max, int *p_x, int *p_y, int *p_z, int *e_x, int *e_y, int *e_z, int *point);
void Init_viewer_info(information *info);
void Make_info_for_viewer_by_shape_func(information *info);
void Make_info_for_viewer_by_shape_func_at_global_parameter_space(information *info);
void Make_boundary_line(information *info);
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
void geo_parameter_coord(const int e, double *para, double *out_coord, information *info);
void physical_coord_R(const int e, double *para, double *R, double *out_coord, information *info);
void geo_parameter_coord_R(const int e_local, double *R_local, double *out_coord, information *info);
void makeChildPosition(int *child_position, int i, int *parent_position, information *info);
bool search2D(target_domain &td_base, target_domain &td, information *info);
void setEdge2D(int edge_num, int *target_axis, double *para, target_domain &temp);
bool search3D(target_domain &td_base, target_domain &td, information *info);
void setEdge3D(int edge_or_face, int edge_or_face_num, int *target_axis, double *para, target_domain &temp);
void elementConnectivity(information *info, vector<vector<int>> &ecn);
// fracture analysis
void Init_increment_field(information *info);
void Fracture_analysis(information *info);
bool get_info_J(vector<J_point_info> &jp_list, vector<J_info> &J_list, information *info, vector<double> &circle_r);
void calc_J(information *info);
void writeVTKFile(const vector<debug_J_vtk>& debug_J, const string& filename, int dim);
void J_sub_element(vector<J_info> &J_list, vector<J_point_info> &jp_list, vector<sub_ele_J> &sub_ele, information *info);
void normal_vector(J_info &J, double *para, int e, J_point_info &nl, vector<vector<double>> &normal, information *info);
bool calc_normalized_cross(double *cross, double *a, double *b);
void vector_normalize(double *vec, int dim);
void calc_q_gradient(int e, double *tilde_coord, int line, int face, J_info &J_list, vector<double> &q_grad, information *info);
void calc_q_gradient_square_patch(int e, double *tilde_coord, int line, int face, J_info &J_list, vector<double> &q_grad, information *info);
void make_virtual_crack_extension_area(J_point_info &jp, J_info &J, size_t line, int e_on_curve, bool islocal, information *info);
double calc_virtual_crack_extension_area(J_info &J, int e, const int face_n, int sub_num, information *info);
double calc_q_1D(double normalized_para_1D, int face);
double calc_q(int e, vector<double> &tilde_coord, int line, int face, J_info &J_list, information *info);
// stress
double calc_equivalent_stress(double *stress);
// make a (for make q gradient)
double Make_a_q(int ele, double *para, double *a_inv_T, information *info, int crack_dir, double delta, bool coef_flag);
double Make_a_q_for_SSIGA(int ele, double *para, double *a_inv_T, information *info, int crack_dir, double delta, bool coef_flag);
void Make_a(int ele, double *para, double *a, information *info);
void Make_a_tilde_prime(int ele, double *para, double *a, information *info, int crack_dir, double delta);

#endif