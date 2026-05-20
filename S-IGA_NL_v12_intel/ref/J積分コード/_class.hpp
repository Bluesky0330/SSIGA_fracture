#ifndef _CLASS
#define _CLASS

#include <vector>
#include "_constant.hpp"

class information;
class target_domain;
class max_min;
class octree_node;
class octree_root;
class subcell;

using namespace std;

class information {
public:
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
	double *a_matrix;
	double *dSF;
	double *dSF_ex;
	double *Gauss_Coordinate;
	double *Gauss_Coordinate_ex;
	double *Jac;
	double *Jac_ex;
	double *B_Matrix;
	double *B_Matrix_ex;
	double *Loc_parameter_on_Glo;
	double *Loc_parameter_on_Glo_ex;
	int *opp_patch_num;
	int *opp_patch_num_ex;

	double *D;
	int *Node_To_Node;
	int *Total_Control_Point_To_Node;
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

	double *deformed_Connectivity_coord;
	double *deformed_Edge_coord;

	double *glo_coordinate;
	double *loc_coordinate;

	// gauss point array
	double *gauss_w;
	double *gauss_point;
	double *gauss_w_ex;
	double *gauss_point_ex;
	double *gauss_w_1D;
	double *gauss_point_1D;
	double *gauss_w_1D_ex;
	double *gauss_point_1D_ex;
	double *Gxi_1D;
	double *w_1D;
	double *Gxi;
	double *w;

	// octree
	vector<vector<subcell>> octree_subcell;

	// element connectivity
	int *e_Patch_check;
	int *e_Face_Edge_info;
	int *e_Opponent_patch_num;
	int *adjacent_element;
};


class target_domain {
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


class max_min {
public:
	int patch_num;
	double min[MAX_DIMENSION];
	double max[MAX_DIMENSION];
	max_min(int input_patch_num) : patch_num(input_patch_num) {}
	~max_min() = default;
};


class octree_node {
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

class octree_root {
public:
	octree_node *root;
	octree_root() : root(nullptr) {}
	~octree_root() { delete root; }
};


class subcell {
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


#endif