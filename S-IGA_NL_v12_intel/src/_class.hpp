#ifndef _CLASS
#define _CLASS

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <any>
#include <memory>
#include <set>
#include <omp.h>
#include "_constant.hpp"

class constant;
class information;
class target_domain;
class max_min;
class octree_node;
class octree_root;
class subcell;
class gauss_points;
class switch_gauss_points;
class average_points;
class distance_norm;
class bf_corner_case;
class blending_faces;
class fictitious_area;
class J_info;
class J_point_info;
class LoadStepNode;
class sparse_vector;

using namespace std;

class alignas(64) constant {
public:
	vector<any> data_vec;

	// general
	int ANALYSIS_MODE_0;
	int ANALYSIS_MODE_1;
	int BASE_PATCH_SIGA;
	int SOLVER;
	int M_MODE;
	int CALC_ON_GP;
	int CALC_ON_ELE_VERTEX;
	int OUTPUT_SVG;
	int OUTPUT_DEFORMED;
	int PARAVIEW_GLO_MODE;
	int PARAVIEW_CRACK_REPRESENTATION;
	int BIN_MODE;
	int MODE_EX;
	int MAX_OCTREE_N;
	int SKIP_J;
	int DM;
	int NG;
	int NG_EXTEND;
	int BLENDING;

	// B-bar method
	int B_BAR;

	// numerical parameters
	double EPS;

	// nonlinear
	double TOTAL_SECONDS;
	int STEP_N;
	int CUTBACK;
	int NEWTON_RAPHSON_MAX_ITR;
	double NEWTON_RAPHSON_EPS;
	int UPDATE_CONTROL_POINT_MAX_ITR;

	// fictitous
	double FICT_COEF;

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
	double *full_K;

	double *sol_vec;
	double *rhs_vec;

	double *Displacement;
	double *Displacement_at_GP;
	double *Strain_at_GP;
	double *Stress_at_GP;
	double *Displacement_at_ele_vertex;
	double *Strain_at_ele_vertex;
	double *Stress_at_ele_vertex;
	double *ReactionForce;

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

	double *glo_coordinate;
	double *loc_coordinate;

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

	// octree
	vector<vector<subcell>> octree_subcell;

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

	// paraview
	vector<vector<int>> Connectivity_all_ele;
	vector<vector<int>> ecn;
	double *elastic_strain_at_connectivity;
	double *equivalent_plastic_strain_at_connectivity;
	double *yield_stress_at_connectivity;
	double *back_stress_at_connectivity;

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

	// blending
	double curvature_radius_coef;
	vector<vector<int>> indices;
	vector<vector<bf_corner_case>> corner_case;
	vector<blending_faces> bf;

	// fictitious area
	vector<fictitious_area> fa;

	// auxiliary
	std::unique_ptr<information> aux = nullptr;

	// decreased order
	std::unique_ptr<information> dec_order = nullptr;
	bool isDecreasedOrderPatch = false;
	vector<vector<double>> L2_M;
	vector<vector<int>> L2_M_ptr;
	vector<vector<int>> L2_M_col;
	vector<vector<sparse_vector>> L2_G;
	// vector<vector<vector<double>>> L2_M_inv_G; // control variable for dilatational strain obtained by L2 projection
	vector<int> Controlpoint_of_Element_in_patch;

	// B-bar (old)
	double *b_average;
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


class alignas(64) octree_node {
public:
	bool is_leaf;
	int level;
	int octree_position[MAX_DIMENSION];
	double para_start[MAX_DIMENSION];
	double para_end[MAX_DIMENSION];
	vector<octree_node *> children;
	vector<int> overlapping_ele;

	octree_node(int *input_octree_position, int input_level, information *info);
	~octree_node()
	{
		for (auto child : children)
			delete child;
	}

	void recieve_para(target_domain &input_td, information *info);
};


class alignas(64) octree_root {
public:
	octree_node *root;
	octree_root() : root(nullptr) {}
	~octree_root() { delete root; }
};


class alignas(64) subcell {
public:
	int e;
	int patch;
	int level;
	double subcell_jac_1D;
	double para_start[MAX_DIMENSION];
	double para_end[MAX_DIMENSION];
	vector<int> overlapping_ele;

	subcell(int input_e, octree_node &input_octree_node, information *info);
	~subcell() = default;

	void get_tilde_para(double *tilde_para, double *child_tilde_para, information *info);
	void get_physical_coord(double *tilde_para, double *child_tilde_para, double *coord, information *info);
};


class alignas(64) gauss_points {
public:
	int patch;
	int ele;
	int n;
	double *para;
	double *w;

	vector<double> blending_function;
	vector<bool> isSkip;

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
	void setBlendingFunction(information *info);

private:
	void makeCoord(information *info);
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

	void setVar(bool input_switch_flag, bool opp_flag, int input_patch, int input_ele, int input_n0, information *info)
	{
		switch_flag = input_switch_flag;
		gp0.set(input_patch, input_ele, input_n0, opp_flag, info);
		gp = &gp0;
	}

	void setVar(bool input_switch_flag, bool opp_flag, int input_patch, int input_ele, int input_n0, int input_n1, information *info)
	{
		switch_flag = input_switch_flag;
		gp0.set(input_patch, input_ele, input_n0, opp_flag, info);
		gp1.set(input_patch, input_ele, input_n1, opp_flag, info);
		gp = &gp0;
	}

	void resetVar(information *info)
	{
		gp0.reset(info);

		if (switch_flag)
			gp1.reset(info);
	}

	void switch_gp(bool flag)
	{
		if (switch_flag)
			gp = flag ? &gp1 : &gp0;
	}

	bool get_switch_flag() const { return switch_flag; }

	void set(int input_patch, int input_ele, int input_n, bool opp_flag, information *info) const { gp->set(input_patch, input_ele, input_n, opp_flag, info); }
	void reset(information *info) const { gp->reset(info); }
	void setBlendingFunction(information *info) const { gp->setBlendingFunction(info); }

	// access wrappers
	int& n() { return access(&gauss_points::n); }
	double* para() { return access(&gauss_points::para); }
	double* w() { return access(&gauss_points::w); }

	vector<double>& blending_function() { return access(&gauss_points::blending_function); }
	vector<bool>& isSkip() { return access(&gauss_points::isSkip); }

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


class alignas(64) bf_corner_case {
public:
	int dim_type;
	vector<int> axis_combination;
	vector<int> face_combination;
	vector<double> area_start;
	vector<double> area_end;

	// ellipse
	vector<double> ellipse_inner_radius;
	vector<double> ellipse_outer_radius;
	vector<double> ellipse_center;
};


class alignas(64) blending_faces {
public:
	int p;
	int dir;
	int face;
	double area_ratio;
};


class alignas(64) fictitious_area {
public:
	int patch;
	vector<double> area_start;
	vector<double> area_end;
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


class alignas(64) LoadStepNode {
public:
	double start;
	double end;
	bool converged = false;
	int step_id = -1; // 元のstep番号（推定用）
	std::vector<std::shared_ptr<LoadStepNode>> children;

	LoadStepNode(double s, double e, int id) : start(s), end(e), step_id(id) {}

	bool is_leaf() const { return children.empty(); }
	double size() const { return end - start; }

	void subdivide()
	{
		double mid = 0.5 * (start + end);
		children.push_back(std::make_shared<LoadStepNode>(start, mid, step_id));
		children.push_back(std::make_shared<LoadStepNode>(mid, end, step_id));
	}
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


// for B-bar method
class alignas(64) sparse_vector {
public:
	sparse_vector() = default;

	void init_thread_local()
	{
		int nt = omp_get_max_threads();
		tl_buffers_.resize(nt);
	}

	void insert(int col, double val)
	{
		if (val != 0.0)
		{
			tl_buffers_[omp_get_thread_num()].push_back({col, val});
		}
	}

	void compress()
	{
		// merge all local buffers
		for (auto &buf : tl_buffers_)
		{
			if (!buf.empty())
			{
				data_.insert(data_.end(), buf.begin(), buf.end());
				buf.clear();
			}
		}

		if (data_.empty())
			return;

		std::sort(data_.begin(), data_.end(),
				  [](auto &a, auto &b)
				  { return a.first < b.first; });

		cols_.clear();
		vals_.clear();
		cols_.reserve(data_.size());
		vals_.reserve(data_.size());

		int cur_c = data_[0].first;
		double cur_v = data_[0].second;

		for (size_t i = 1; i < data_.size(); i++)
		{
			auto [c, v] = data_[i];
			if (c == cur_c)
			{
				cur_v += v;
			}
			else
			{
				cols_.push_back(cur_c);
				vals_.push_back(cur_v);
				cur_c = c;
				cur_v = v;
			}
		}
		cols_.push_back(cur_c);
		vals_.push_back(cur_v);

		data_.clear();
		data_.shrink_to_fit();
	}

	double get(int col) const noexcept
	{
		auto it = std::lower_bound(cols_.begin(), cols_.end(), col);
		if (it != cols_.end() && *it == col)
		{
			return vals_[it - cols_.begin()];
		}
		return 0.0;
	}

	size_t size() const noexcept { return cols_.size(); }

	int col_at(size_t index) const { return cols_[index]; }
	double val_at(size_t index) const { return vals_[index]; }

private:
	// instance-owned thread-local buffers
	std::vector<std::vector<std::pair<int, double>>> tl_buffers_;

	std::vector<std::pair<int, double>> data_;
	std::vector<int> cols_;
	std::vector<double> vals_;
};

#endif