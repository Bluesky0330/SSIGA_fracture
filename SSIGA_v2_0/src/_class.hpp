#ifndef _CLASS
#define _CLASS

#include <iostream>
#include <vector>
#include <string>
#include <any>
#include <memory>
#include <set>
#include "_constant.hpp"

class constant;
class information;
class target_domain;
class max_min;
class gauss_points;
class switch_gauss_points;
class average_points;
class distance_norm;
class J_info;
class J_point_info;

using namespace std;

class alignas(64) constant {
public:
	vector<any> data_vec;

	// general
	int ANALYSIS_MODE;
	int GEOMETRY_ONLY_OUTPUT;
	int OUTPUT_GLOBAL_PARAMETERS;
	int OUTPUT_PARAVIEW;
	int CALC_ON_GP;
	int CALC_ON_ELE_VERTEX;

	// gaussian quadrature
	int USE_EXTENDED_QUADRATURE;
	int NUM_GAUSS_POINTS;
	int NUM_GAUSS_POINTS_EXTENDED;
	
	// fracture analysis settings
	int FRACTURE_MODE;
	int CALCLATE_DISPLACEMENT;
	int INTEGRAL_DOMAIN_TYPE;
	int CRACK_TYPE;
	
	void get_data(const char *s, information *info);
};


class alignas(64) information {
public:
	constant c;
	vector<string> opt_files;

	int DIMENSION;

	int *Total_Knot_to_mesh;
	int *Total_Patch_on_mesh;
	int *Total_Patch_to_mesh;
	int *Total_Control_Point_on_mesh;
	int *Total_Control_Point_to_mesh;
	int *Total_Element_on_mesh;
	int *Total_Element_to_mesh;
	int *real_Total_Element_on_mesh;
	int *real_Total_Element_to_mesh;
	int *Total_Load_to_mesh;
	int *Total_Constraint_to_mesh;
	int *Total_DistributeForce_to_mesh;

	int *Order;
	int *No_knot;
	int *Total_Control_Point_to_patch;
	int *Total_Knot_to_patch_dim;
	double *Position_Knots;
	int *No_Control_point;
	int *No_Control_point_in_patch;
	int *Patch_Control_point;
	int *No_Control_point_ON_ELEMENT;
	double *Node_Coordinate;
	double *Control_Coord_x;
	double *Control_Coord_y;
	double *Control_Coord_z;
	double *Control_Weight;
	int *Constraint_Node_Dir;
	double *Value_of_Constraint;
	int *Load_Node_Dir;
	double *Value_of_Load;
	int *iPatch_array;
	int *iCoord_array;
	int *jCoord_array;
	int *type_load_array;
	double *val_Coord_array;
	double *Range_Coord_array;
	double *Coeff_Dist_Load_array;

	int *INC;
	int *Controlpoint_of_Element;
	int *Element_patch;
	int *Element_mesh;
	int *line_No_real_element;
	int *line_No_Total_element;
	double *difference;
	int *Total_element_all_ID;
	int *ENC;
	int *real_element_line;
	int *real_element;
	int *real_El_No_on_mesh;
	double *Equivalent_Nodal_Force;

	vector<vector<int>> eoi;

	vector<set<int>> ptr;

	double *D;
	int *Index_Dof;
	int *K_Whole_Ptr;
	int *K_Whole_Col;
	double *K_Whole_Val;

	double *sol_vec;
	double *rhs_vec;

	double *Displacement;

	int *Connectivity;
	int *Connectivity_ele;
	int *Connectivity_point;
	double *Connectivity_coord;
	int *Patch_check;
	int *Patch_array;
	int *Face_Edge_info;
	int *Opponent_patch_num;
	double *Edge_coord;

	double *disp_at_connectivity;
	double *strain_at_connectivity;
	double *stress_at_connectivity;
	double *vm_at_connectivity;
	double *hs_at_connectivity;

	// output
	double *Displacement_at_GP;
	double *PhysicalCoordinate_at_GP;
	double *Strain_at_GP;
	double *Stress_at_GP;
	double *Displacement_at_ele_vertex;
	double *PhysicalCoordinate_at_ele_vertex;
	double *Strain_at_ele_vertex;
	double *Stress_at_ele_vertex;
	double *ReactionForce;

	// gauss point array
	bool isSetGaussArray = false;
	double *gauss_w;
	double *gauss_point;
	double *gauss_w_ex;
	double *gauss_point_ex;
	double *gauss_w_1D;
	double *gauss_point_1D;
	double *gauss_w_1D_ex;
	double *gauss_point_1D_ex;
	vector<switch_gauss_points> gp;

	// element connectivity
	int *e_Patch_check;
	int *e_Face_Edge_info;
	int *e_Opponent_patch_num;
	int *adjacent_element;

	// for Elastoplasticity
	double *delta_F;
	double *F;
	double *residual_vec;
	double *internal_force;
	double *external_force;
	double *previous_external_force;
	double *rhs_vec_initial;
	double *forced_disp_T;
	double *delta_u;
	double *total_delta_u;

	double *parent_para_loc;
	int *opp_ele_totaldisp_loc;
	double *opp_parent_para_totaldisp_loc;

	// deformed output for paraview
	double *deformed_Connectivity_coord;
	double *deformed_Edge_coord;

	double *disp_increment;
	double *disp;
	double *disp_previous;

	double *deformation_gradient;
	double *current_deformation_gradient;
	double *elastic_strain;
	double *elastic_strain_trial;
	double *current_elastic_strain;
	double *back_stress;
	double *current_back_stress;
	double *stress;
	double *current_stress;
	double *equivalent_stress;
	double *yield_stress;
	double *current_yield_stress;
	double *equivalent_plastic_strain;
	double *equivalent_plastic_strain_increment;

	double kinematic_hardening_fraction;
	int ss_curve_num;
	double *ss_curve_stress;
	double *ss_curve_plastic_strain;

	// for SSIGA
	int Global_local_patch; 					// ローカルパッチが属するグローバルパッチの番号
	int Geo_real_Total_Element_on_mesh;
	int Geo_real_Total_Element_to_mesh = 0;
	int Geo_Total_patch_on_mesh;
	int Geo_Total_patch_to_mesh = 0;			// ジオメトリ表現のためのメッシュは一個だけ
	int Geo_Total_Element_on_mesh;
	int Geo_Total_Element_to_mesh = 0;			// ジオメトリ表現のためのメッシュは一個だけ
	int Geo_Total_Control_Point_on_mesh;
	int Geo_Total_Control_Point_to_mesh = 0;	// ジオメトリ表現のためのメッシュは一個だけ
	int Geo_Total_Knot_on_mesh;
	int Geo_Total_Knot_to_mesh = 0;				// ジオメトリ表現のためのメッシュは一個だけ
	int *Geo_Order;
	int *Geo_No_knot;
	int *Geo_No_Control_point;
	int *Geo_No_Control_point_in_patch;
	int *Geo_Total_Control_Point_to_patch;
	int *Geo_Patch_Control_point;
	double *Geo_Position_Knots;
	int *Geo_No_Control_point_ON_ELEMENT;
	int *Geo_Total_Knot_to_patch_dim;
	int *Geo_Element_patch;
	double *Geo_Node_Coordinate;
	double *Control_Coord_x_geo;
	double *Control_Coord_y_geo;
	double *Control_Coord_z_geo;
	double *Control_Weight_geo;
	int *Geo_INC;
	int *Geo_Controlpoint_of_Element;
	int *Geo_Element_mesh;
	int *Geo_line_No_real_element;
	int *Geo_line_No_Total_element;
	double *Geo_difference;
	int *Geo_Total_element_all_ID;
	int *Geo_ENC;
	int *Geo_real_element_line;
	int *Geo_real_element;
	int *Geo_real_El_No_on_mesh;
	int Geo_real_El_No_to_mesh = 0;
	double *Geo_Equivalent_Nodal_Force;

	// paraview
	vector<vector<int>> Connectivity_all_ele;
	vector<vector<int>> ecn;

	// load disp curve
	int load_disp_element_n;
	int *load_disp_element;
	int *load_disp_element_direction;
	int *load_disp;

	// disp overlay increment
	double *disp_overlay;
	double *disp_overlay_increment;

	// crack pair
	int *crack_pair_n_to_mesh;
	int *crack_pair;
};


class alignas(64) target_domain {
public:
	int e;
	double para_start[MAX_DIMENSION];
	double para_end[MAX_DIMENSION];

	target_domain(int input_e) : e(input_e) {}
	~target_domain() = default;

	void setPara(int octree_level, int *octree_position, information *info);

	bool operator<(const target_domain& other) const {
		return e < other.e;
	}

	bool operator==(const target_domain& other) const {
		return e == other.e;
	}
};


class alignas(64) max_min {
public:
	int patch_num;
	double min[MAX_DIMENSION];
	double max[MAX_DIMENSION];
	max_min(int input_patch_num) : patch_num(input_patch_num) {}
	~max_min() = default;
};


class alignas(64) gauss_points {
public:
	int patch;
	int ele;
	int n;
	double *para;
	double *w;

	vector<double> coord;
	vector<bool> isOverlay;
	vector<int> opp_patch;
	vector<int> opp_ele;
	vector<double> opp_para;
	vector<double> opp_para_tilde;

	vector<double> deformation_gradient;
	vector<double> current_deformation_gradient;
	vector<double> elastic_strain;
	vector<double> elastic_strain_trial;
	vector<double> current_elastic_strain;
	vector<double> back_stress;
	vector<double> current_back_stress;
	vector<double> stress;
	vector<double> current_stress;
	vector<double> equivalent_stress;
	vector<double> yield_stress;
	vector<double> current_yield_stress;
	vector<double> equivalent_plastic_strain;
	vector<double> equivalent_plastic_strain_increment;

	vector<double> strain_energy_density;
	vector<double> delta_strain_energy_density;

	vector<double> previous_elastic_strain;
	vector<double> previous_stress;

	void set(int input_patch, int input_ele, int input_n, bool opp_flag, information *info);
	void reset(information *info);
	// void setBlendingFunction(information *info);

private:
	void makeCoord(bool islocal, information *info);
	void makeOpp(information *info);
};


class alignas(64) switch_gauss_points {
private:
	bool switch_flag;
	gauss_points *gp;
	gauss_points gp0;
	gauss_points gp1;

	// template accesser
	template<typename T>
	T& access(T gauss_points::* member) {
		return gp->*member;
	}

public:
	switch_gauss_points() : gp0(), gp1() {}

	void setVar(bool opp_flag, int input_patch, int input_ele, int input_n0, information *info)
	{
		switch_flag = false;
		gp0.set(input_patch, input_ele, input_n0, opp_flag, info);
		gp = &gp0;
	}

	void setVar(bool opp_flag, int input_patch, int input_ele, int input_n0, int input_n1, information *info)
	{
		switch_flag = true;
		gp0.set(input_patch, input_ele, input_n0, opp_flag, info);
		gp1.set(input_patch, input_ele, input_n1, opp_flag, info);
		gp = &gp0;
	}

	void switch_gp(bool flag)
	{
		if (switch_flag)
			gp = flag ? &gp1 : &gp0;
	}

	bool get_switch_flag() const { return switch_flag; }

	void set(int input_patch, int input_ele, int input_n, bool opp_flag, information *info) const { gp->set(input_patch, input_ele, input_n, opp_flag, info); }

	// access wrappers
	int& n() { return access(&gauss_points::n); }
	double* para() { return access(&gauss_points::para); }
	double* w() { return access(&gauss_points::w); }

	vector<double>& coord() { return access(&gauss_points::coord); }
	vector<bool>& isOverlay() { return access(&gauss_points::isOverlay); }
	vector<int>& opp_patch() { return access(&gauss_points::opp_patch); }
	vector<int>& opp_ele() { return access(&gauss_points::opp_ele); }
	vector<double>& opp_para() { return access(&gauss_points::opp_para); }
	vector<double>& opp_para_tilde() { return access(&gauss_points::opp_para_tilde); }

	vector<double>& deformation_gradient() { return access(&gauss_points::deformation_gradient); }
	vector<double>& current_deformation_gradient() { return access(&gauss_points::current_deformation_gradient); }
	vector<double>& elastic_strain() { return access(&gauss_points::elastic_strain); }
	vector<double>& elastic_strain_trial() { return access(&gauss_points::elastic_strain_trial); }
	vector<double>& current_elastic_strain() { return access(&gauss_points::current_elastic_strain); }
	vector<double>& back_stress() { return access(&gauss_points::back_stress); }
	vector<double>& current_back_stress() { return access(&gauss_points::current_back_stress); }
	vector<double>& stress() { return access(&gauss_points::stress); }
	vector<double>& current_stress() { return access(&gauss_points::current_stress); }
	vector<double>& equivalent_stress() { return access(&gauss_points::equivalent_stress); }
	vector<double>& yield_stress() { return access(&gauss_points::yield_stress); }
	vector<double>& current_yield_stress() { return access(&gauss_points::current_yield_stress); }
	vector<double>& equivalent_plastic_strain() { return access(&gauss_points::equivalent_plastic_strain); }
	vector<double>& equivalent_plastic_strain_increment() { return access(&gauss_points::equivalent_plastic_strain_increment); }

	vector<double>& strain_energy_density() { return access(&gauss_points::strain_energy_density); }
	vector<double>& delta_strain_energy_density() { return access(&gauss_points::delta_strain_energy_density); }

	vector<double>& previous_elastic_strain() { return access(&gauss_points::previous_elastic_strain); }
	vector<double>& previous_stress() { return access(&gauss_points::previous_stress); }

};


class alignas(64) average_points {
public:
	int n;
	vector<double> w;
	vector<double> target_e;
	vector<double> target_gp_num;
};


class alignas(64) distance_norm {
public:
	double r;
	int e;
	int gp_num;

	distance_norm() : r(0.0), e(0), gp_num(0) {}

	bool operator<(const distance_norm &other) const {
		return r < other.r;
	}
};


class alignas(64) J_info {
public:
	int p;
	int crack_dir;
	int r_dir;
	int normal_info_num;
	int sub_n;
	vector<vector<int>> e_list;

	int quadrant;

	J_info(int patch, int dir1, int dir2, int normal_info_num1, int sub_n1) : p(patch), crack_dir(dir1), r_dir(dir2), normal_info_num(normal_info_num1), sub_n(sub_n1) {}
};


class alignas(64) J_point_info {
public:
	int p;
	vector<int> ele_in_normal_patch;
	int crack_dir;
	int r_dir;
	int tangent_dir[2];
	vector<vector<double>> coord;
	vector<vector<double>> normal;
	vector<double> J_val;
	vector<double> J_integral_area;

	vector<double> delta_arc_length;
	vector<double> delta_arc_length_at_half_para;

	// interaction integral method
	vector<vector<double>> IIM_K_val;
	vector<vector<vector<double>>> Q;

	J_point_info(information *info, int patch, int dir1, int dir2, int dir3, int dir4);
	void init(J_info &J, information *info);
};


// debug J
class alignas(64) debug_J_vtk {
public:
	vector<double> coord;

	double q;
	vector<double> q_grad;
	vector<double> normal;

	debug_J_vtk(int dim) : coord(dim, 0.0), q(0.0), q_grad(dim, 0.0), normal(dim, 0.0) {}
};


// subelement for J
class alignas(64) sub_ele_J {
public:
	int gp_num_1d;
	int gp_num;

	int e;
	int p;

	int sub_n;
	int face_n;

	int line;
	int slice;
	int face;

	int J_list_index;
	int jp_list_index;

	int crack_dir;

	vector<double> para_tilde;
	double *gp;
	double *w;

	// for overlay
	vector<bool> overlay;
	vector<int> overlay_ele;
	vector<double> coord;
	vector<double> opp_para_tilde;

	// internal state variables
	vector<vector<double>> elastic_strain;
	vector<vector<double>> stress;
	vector<vector<double>> previous_elastic_strain;
	vector<vector<double>> previous_stress;
	vector<double> yield_stress;
	vector<double> equivalent_plastic_strain;
	vector<double> previous_W;

	sub_ele_J(int gp_num_1d, int gp_num, int e, int p, int sub_n, int face_n, int line, int slice, int face, int J_list_index, int jp_list_index, int crack_dir, double *sub_gp_para, double *sub_gp_weight) : gp_num_1d(gp_num_1d), gp_num(gp_num), e(e), p(p), sub_n(sub_n), face_n(face_n), line(line), slice(slice), face(face), J_list_index(J_list_index), jp_list_index(jp_list_index), crack_dir(crack_dir), gp(sub_gp_para), w(sub_gp_weight) {}
	void init(information *info);
};

#endif