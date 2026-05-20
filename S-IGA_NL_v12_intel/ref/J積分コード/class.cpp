// header
#include "_class.hpp"
#include "_func.hpp"
#include "_sub.hpp"

void target_domain::setPara(int octree_level, int *octree_position, information *info)
{
	int division = pow_int(2, octree_level);

	int p = info->Element_patch[e];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int offset = info->ENC[e * info->DIMENSION + i];
		double s = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + i] + info->Order[p * info->DIMENSION + i] + offset];
		double e = info->Position_Knots[info->Total_Knot_to_patch_dim[p * info->DIMENSION + i] + info->Order[p * info->DIMENSION + i] + offset + 1];
		double id = (double)octree_position[i];
		para_start[i] = (s + (e - s) * id / (double)division);
		para_end[i] = (s + (e - s) * (id + 1.0) / (double)division);
	}
}


octree_node::octree_node(int *input_octree_position, int input_level, information *info) : is_leaf(true), level(input_level)
{
	for (int i = 0; i < info->DIMENSION; i++)
		octree_position[i] = input_octree_position[i];
}


void octree_node::recieve_para(target_domain &input_td, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		para_start[i] = input_td.para_start[i];
		para_end[i] = input_td.para_end[i];
	}
}


subcell::subcell(int input_e, octree_node &input_octree_node, information *info) : e(input_e), patch(info->Element_patch[input_e]), level(input_octree_node.level)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		para_start[i] = input_octree_node.para_start[i];
		para_end[i] = input_octree_node.para_end[i];
	}
	for (size_t i = 0; i < input_octree_node.overlapping_ele.size(); i++)
		overlapping_ele.emplace_back(input_octree_node.overlapping_ele[i]);

	subcell_jac_1D = 1.0 / static_cast<double>(pow_int(2, level));
}


void subcell::get_tilde_para(double *tilde_para, double *child_tilde_para, information *info)
{
	double para[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		para[i] = para_start[i] + 0.5 * (child_tilde_para[i] + 1.0) * (para_end[i] - para_start[i]);
	tilde_coord(tilde_para, para, patch, e, info);
}


void subcell::get_physical_coord(double *tilde_para, double *child_tilde_para, double *coord, information *info)
{
	get_tilde_para(tilde_para, child_tilde_para, info);

	for (int i = 0; i < info->DIMENSION; i++)
		coord[i] = 0.0;

	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
	{
		double R = Shape_func(i, tilde_para, e, info);
		for (int j = 0; j < info->DIMENSION; j++)
			coord[j] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
	}
}