#ifndef _MI3D
#define _MI3D

class info_global {
public:
    int DIMENSION;
    int Total_Control_Point;
    int mode;
    double *Weight;
    double *temp_Weight;
};

class info_each_DIMENSION {
public:
    int CP_n;
    int Order;
    int OE_n;
    int knot_n;
    int KI_cp_n;
    int KI_non_uniform_n;
    double *CP;
    double *KV;
    double *temp_CP;
    double *temp_KV;
    double *insert_knot;
};

class line_coordinate {
public:
    double *line;
    double *new_line;
};

// Main
void Calc_MI3D(const int DIM, const int *Order, const int *CP_n, double *knot_array);
// set default patch
void Set_default_patch_1(const int *Order, const int *CP_n, info_global *info_glo, info_each_DIMENSION *info);
void Set_default_patch_2(info_global *info_glo, info_each_DIMENSION *info);
// Knot Insertion
void KI_non_uniform(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_global *info_glo, info_each_DIMENSION *info);
void KI_calc_knot_1D(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_each_DIMENSION *info);
void KI_calc_T_1D(int insert_axis, int insert_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM);
void KI_cp(int insert_axis, info_global *info_glo, info_each_DIMENSION *info);
void KI_cp_not_open_knot_vec(int insert_axis, info_global *info_glo, info_each_DIMENSION *info);
// Order Elevation
void OE(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info);
void OR(int reduction_axis, int OR_n, info_global *info_glo, info_each_DIMENSION *info);
void Calc_insert_knot_in_OE(int elevation_axis, int *insert_knot_n, double *insert_knot, int *removal_knot_n, double *removal_knot, info_each_DIMENSION *info);
void Calc_Bezier(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM);
void Bezier_Order_Elevation(int elevation_axis, int Bezier_line, int *counter, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM);
// Knot Removal
void KR_non_uniform(int removal_axis, int removal_knot_n, double *removal_knot, info_global *info_glo, info_each_DIMENSION *info);
void KR_calc_knot_1D(int removal_axis, int removal_knot_n, double *removal_knot, info_each_DIMENSION *info);
void KR_calc_Tinv_1D(int removal_axis, int removal_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM);

#endif