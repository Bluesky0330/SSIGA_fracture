// header
#include "_header.hpp"
#include "_sub.hpp"

using namespace std;

// tool
double InverseMatrix_2x2(double Matrix[2][2])
{
	int i, j;
	double a[2][2];
	double det = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			a[i][j] = Matrix[i][j];

	Matrix[0][0] = a[1][1] / det;
	Matrix[0][1] = a[0][1] * (-1) / det;
	Matrix[1][0] = a[1][0] * (-1) / det;
	Matrix[1][1] = a[0][0] / det;

	return det;
}


double InverseMatrix_2x2(double *Matrix)
{
	int i;
	double a[4];
	double det = Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2];
	double one_over_det = 1.0 / det;

	for (i = 0; i < 4; i++)
		a[i] = Matrix[i];

	Matrix[0] =  a[3] * one_over_det;
	Matrix[1] = -a[1] * one_over_det;
	Matrix[2] = -a[2] * one_over_det;
	Matrix[3] =  a[0] * one_over_det;

	return det;
}


double InverseMatrix_3x3(double Matrix[3][3])
{
	int i, j;
	double a[3][3];
	double det = Matrix[0][0] * Matrix[1][1] * Matrix[2][2] + Matrix[1][0] * Matrix[2][1] * Matrix[0][2] + Matrix[2][0] * Matrix[0][1] * Matrix[1][2] - Matrix[0][0] * Matrix[2][1] * Matrix[1][2] - Matrix[2][0] * Matrix[1][1] * Matrix[0][2] - Matrix[1][0] * Matrix[0][1] * Matrix[2][2];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			a[i][j] = Matrix[i][j];

	Matrix[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	Matrix[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	Matrix[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	Matrix[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	Matrix[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	Matrix[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	Matrix[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	Matrix[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	Matrix[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;

	return det;
}


double InverseMatrix_3x3(double *Matrix)
{
	int i;
	double a[9];
	double det = Matrix[0] * Matrix[4] * Matrix[8] + Matrix[3] * Matrix[7] * Matrix[2] + Matrix[6] * Matrix[1] * Matrix[5] - Matrix[0] * Matrix[7] * Matrix[5] - Matrix[6] * Matrix[4] * Matrix[2] - Matrix[3] * Matrix[1] * Matrix[8];
	double one_over_det = 1.0 / det;

	for (i = 0; i < 9; i++)
		a[i] = Matrix[i];

	Matrix[0] = (a[4] * a[8] - a[5] * a[7]) * one_over_det;
	Matrix[1] = (a[2] * a[7] - a[1] * a[8]) * one_over_det;
	Matrix[2] = (a[1] * a[5] - a[2] * a[4]) * one_over_det;
	Matrix[3] = (a[5] * a[6] - a[3] * a[8]) * one_over_det;
	Matrix[4] = (a[0] * a[8] - a[2] * a[6]) * one_over_det;
	Matrix[5] = (a[2] * a[3] - a[0] * a[5]) * one_over_det;
	Matrix[6] = (a[3] * a[7] - a[4] * a[6]) * one_over_det;
	Matrix[7] = (a[1] * a[6] - a[0] * a[7]) * one_over_det;
	Matrix[8] = (a[0] * a[4] - a[1] * a[3]) * one_over_det;

	return det;
}


long pow_int(int val, int n)
{
	long result = 1;
	for (int i = 0; i < n; i++)
		result *= val;
	return result;
}


double Determinant_2x2(const double *Matrix)
{
	return Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2];
}


double Determinant_3x3(const double *Matrix)
{
	return Matrix[0] * Matrix[4] * Matrix[8] + Matrix[3] * Matrix[7] * Matrix[2] + Matrix[6] * Matrix[1] * Matrix[5] - Matrix[0] * Matrix[7] * Matrix[5] - Matrix[6] * Matrix[4] * Matrix[2] - Matrix[3] * Matrix[1] * Matrix[8];
}


// Shape Function
void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	Position_Data_param[j] = ((info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) * Local_coord[j] + (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] + info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i])) / 2.0;
}


// ローカル形状表現のための要素パラメータ空間からの形状関数
void Geo_ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info)
{
	int i = info->Geo_INC[info->Geo_Element_patch[e] * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	Position_Data_param[j] = ((info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i]) * Local_coord[j] + (info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i + 1] + info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i])) / 2.0;
}


double dShapeFunc_from_paren(int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	double dPosition_Data_param = (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) / 2.0;
	return dPosition_Data_param;
}


// ローカル形状表現のための要素パラメータ空間からの形状関数の微分
double Geo_dShapeFunc_from_paren(int j, int e, information *info)
{
	int i = info->Geo_INC[info->Geo_Element_patch[e] * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	double dPosition_Data_param = (info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[info->Geo_Element_patch[e] * info->DIMENSION + j] + i]) / 2.0;
	return dPosition_Data_param;
}


void shape_function_1D(double knot, int dim, int e, double *Shape, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	int patch = info->Element_patch[e];
	int offset = info->ENC[e * info->DIMENSION + dim];
	int p = info->Order[patch * info->DIMENSION + dim];
	int knot_n = info->No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];

	int support = p + 1;
	for (int i = 0; i < support + p; i++)
	{
		if (fabs(knot_ptr[offset + i] - knot_ptr[offset + i + 1]) < MERGE_ERROR)
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;

		else if (knot_ptr[offset + i] <= knot && knot < knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else if (fabs(knot_ptr[offset + i + 1] - knot_ptr[knot_n - 1]) < MERGE_ERROR && knot_ptr[offset + i] <= knot && knot <= knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;
	}

	double left_term, right_term;
	for (int i = 1; i <= p; i++)
	{
		int new_support = support + p - i;
		for (int j = 0; j < new_support; j++)
		{
			left_term = 0.0;
			right_term = 0.0;

			// left
			if ((knot - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j] == 0 && knot_ptr[offset + j + i] - knot_ptr[offset + j] == 0)
				left_term = 0.0;
			else
				left_term = (knot - knot_ptr[offset + j]) / (knot_ptr[offset + j + i] - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j];

			// right
			if ((knot_ptr[offset + j + i + 1] - knot) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1] == 0 && knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1] == 0)
				right_term = 0.0;
			else
				right_term = (knot_ptr[offset + j + i + 1] - knot) / (knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1];

			Shape[(dim * MAX_ORDER + i) * max_support_1D + j] = left_term + right_term;
		}
	}
}


void shape_function_1D(double knot, int dim, int e, double *Shape, double *dShape, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	int patch = info->Element_patch[e];
	int offset = info->ENC[e * info->DIMENSION + dim];
	int p = info->Order[patch * info->DIMENSION + dim];
	int knot_n = info->No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];

	int support = p + 1;
	for (int i = 0; i < support + p; i++)
	{
		if (fabs(knot_ptr[offset + i] - knot_ptr[offset + i + 1]) < MERGE_ERROR)
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;

		else if (knot_ptr[offset + i] <= knot && knot < knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else if (fabs(knot_ptr[offset + i + 1] - knot_ptr[knot_n - 1]) < MERGE_ERROR && knot_ptr[offset + i] <= knot && knot <= knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;
	}

	double left_term, right_term;
	for (int i = 1; i <= p; i++)
	{
		int new_support = support + p - i;
		for (int j = 0; j < new_support; j++)
		{
			left_term = 0.0;
			right_term = 0.0;

			// left
			if ((knot - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j] == 0 && knot_ptr[offset + j + i] - knot_ptr[offset + j] == 0)
				left_term = 0.0;
			else
				left_term = (knot - knot_ptr[offset + j]) / (knot_ptr[offset + j + i] - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j];

			// right
			if ((knot_ptr[offset + j + i + 1] - knot) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1] == 0 && knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1] == 0)
				right_term = 0.0;
			else
				right_term = (knot_ptr[offset + j + i + 1] - knot) / (knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1];

			Shape[(dim * MAX_ORDER + i) * max_support_1D + j] = left_term + right_term;
		}
	}

	// derivative
	double dleft_term, dright_term;
	for (int i = 0; i < support; i++)
	{
		dleft_term = 0.0;
		dright_term = 0.0;

		// left
		if (p * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i] == 0 && knot_ptr[offset + i + p] - knot_ptr[offset + i] == 0)
			dleft_term = 0.0;
		else
			dleft_term = p / (knot_ptr[offset + i + p] - knot_ptr[offset + i]) * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i];

		// right
		if (p * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i + 1] == 0 && knot_ptr[offset + i + p + 1] - knot_ptr[offset + i + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = p / (knot_ptr[offset + i + p + 1] - knot_ptr[offset + i + 1]) * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i + 1];

		dShape[dim * max_support_1D + i] = dleft_term - dright_term;
	}
}


// ローカル形状表現のための1方向形状関数
void Geo_Shape_Function_1D(double knot, int dim, int e, double *Shape, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	int patch = info->Geo_Element_patch[e];
	int offset = info->Geo_ENC[e * info->DIMENSION + dim];
	int p = info->Geo_Order[patch * info->DIMENSION + dim];
	int knot_n = info->Geo_No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];

	int support = p + 1;
	for (int i = 0; i < support + p; i++)
	{
		if (fabs(knot_ptr[offset + i] - knot_ptr[offset + i + 1]) < MERGE_ERROR)
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;

		else if (knot_ptr[offset + i] <= knot && knot < knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else if (fabs(knot_ptr[offset + i + 1] - knot_ptr[knot_n - 1]) < MERGE_ERROR && knot_ptr[offset + i] <= knot && knot <= knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;
	}

	double left_term, right_term;
	for (int i = 1; i <= p; i++)
	{
		int new_support = support + p - i;
		for (int j = 0; j < new_support; j++)
		{
			left_term = 0.0;
			right_term = 0.0;

			// left
			if ((knot - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j] == 0 && knot_ptr[offset + j + i] - knot_ptr[offset + j] == 0)
				left_term = 0.0;
			else
				left_term = (knot - knot_ptr[offset + j]) / (knot_ptr[offset + j + i] - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j];

			// right
			if ((knot_ptr[offset + j + i + 1] - knot) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1] == 0 && knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1] == 0)
				right_term = 0.0;
			else
				right_term = (knot_ptr[offset + j + i + 1] - knot) / (knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1];

			Shape[(dim * MAX_ORDER + i) * max_support_1D + j] = left_term + right_term;
		}
	}
}


void Geo_Shape_Function_1D(double knot, int dim, int e, double *Shape, double *dShape, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	int patch = info->Geo_Element_patch[e];
	int offset = info->Geo_ENC[e * info->DIMENSION + dim];
	int p = info->Geo_Order[patch * info->DIMENSION + dim];
	int knot_n = info->Geo_No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];
	int support = p + 1;
	for (int i = 0; i < support + p; i++)
	{
		if (fabs(knot_ptr[offset + i] - knot_ptr[offset + i + 1]) < MERGE_ERROR)
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;

		else if (knot_ptr[offset + i] <= knot && knot < knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else if (fabs(knot_ptr[offset + i + 1] - knot_ptr[knot_n - 1]) < MERGE_ERROR && knot_ptr[offset + i] <= knot && knot <= knot_ptr[offset + i + 1])
			Shape[dim * MAX_ORDER * max_support_1D + i] = 1.0;

		else
			Shape[dim * MAX_ORDER * max_support_1D + i] = 0.0;
	}

	double left_term, right_term;
	for (int i = 1; i <= p; i++)
	{
		int new_support = support + p - i;
		for (int j = 0; j < new_support; j++)
		{
			left_term = 0.0;
			right_term = 0.0;

			// left
			if ((knot - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j] == 0 && knot_ptr[offset + j + i] - knot_ptr[offset + j] == 0)
				left_term = 0.0;
			else
				left_term = (knot - knot_ptr[offset + j]) / (knot_ptr[offset + j + i] - knot_ptr[offset + j]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j];

			// right
			if ((knot_ptr[offset + j + i + 1] - knot) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1] == 0 && knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1] == 0)
				right_term = 0.0;
			else
				right_term = (knot_ptr[offset + j + i + 1] - knot) / (knot_ptr[offset + j + i + 1] - knot_ptr[offset + j + 1]) * Shape[(dim * MAX_ORDER + i - 1) * max_support_1D + j + 1];

			Shape[(dim * MAX_ORDER + i) * max_support_1D + j] = left_term + right_term;
		}
	}

	// derivative
	double dleft_term, dright_term;
	for (int i = 0; i < support; i++)
	{
		dleft_term = 0.0;
		dright_term = 0.0;

		// left
		if (p * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i] == 0 && knot_ptr[offset + i + p] - knot_ptr[offset + i] == 0)
			dleft_term = 0.0;
		else
			dleft_term = p / (knot_ptr[offset + i + p] - knot_ptr[offset + i]) * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i];

		// right
		if (p * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i + 1] == 0 && knot_ptr[offset + i + p + 1] - knot_ptr[offset + i + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = p / (knot_ptr[offset + i + p + 1] - knot_ptr[offset + i + 1]) * Shape[(dim * MAX_ORDER + p - 1) * max_support_1D + i + 1];

		dShape[dim * max_support_1D + i] = dleft_term - dright_term;
	}
}


void shape_and_dshape(double *R, double *Local_coord, int e, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);

	int patch = info->Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		shape_function_1D(Position_Data_param[i], i, e, Shape.data(), info);
		basis_num_offset[i] = info->ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;
}


void shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);
	static thread_local vector<double> dShape(info->DIMENSION * max_support_1D);

	int patch = info->Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		shape_function_1D(Position_Data_param[i], i, e, Shape.data(), dShape.data(), info);
		basis_num_offset[i] = info->ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;

	// make dR
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;
	double dWeight_func3 = 0.0;
	int basis_num_dim[MAX_DIMENSION];
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * w_coord;
		}
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func3 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] * w_coord;
		}

		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 2] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] - dWeight_func3 * shape_func[i]) / (weight_func * weight_func);
		}
	}

	if (calc_dR_coef_flag)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			double coef = dShapeFunc_from_paren(i, e, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
				dR[j * info->DIMENSION + i] *= coef;
		}
	}
}


// ローカル形状表現用の形状関数とその微分
void geo_shape_and_dshape(double *R, double *Local_coord, int e, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);

	int patch = info->Geo_Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Geo_Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		Geo_ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		Geo_Shape_Function_1D(Position_Data_param[i], i, e, Shape.data(), info);
		basis_num_offset[i] = info->Geo_ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;
}


void geo_shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);
	static thread_local vector<double> dShape(info->DIMENSION * max_support_1D);

	int patch = info->Geo_Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Geo_Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		Geo_ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		Geo_Shape_Function_1D(Position_Data_param[i], i, e, Shape.data(), dShape.data(), info);
		basis_num_offset[i] = info->Geo_ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;

	// make dR
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;
	double dWeight_func3 = 0.0;
	int basis_num_dim[MAX_DIMENSION];
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * w_coord;
		}
		for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func3 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] * w_coord;
		}

		for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 2] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] - dWeight_func3 * shape_func[i]) / (weight_func * weight_func);
		}
	}

	if (calc_dR_coef_flag)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			double coef = Geo_dShapeFunc_from_paren(i, e, info);
			for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[patch]; j++)
				dR[j * info->DIMENSION + i] *= coef;
		}
	}
}


// Bspline基底関数
void Bspline_shape_and_dshape(double *R, double *Local_coord, int e, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);

	int patch = info->Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		shape_function_1D(Position_Data_param[i], i, e, Shape.data(), info);
		basis_num_offset[i] = info->ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * 1.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * 1.0 / weight_func;
}


void Bspline_shape_and_dshape(double *R, double *dR, double *Local_coord, int e, bool calc_dR_coef_flag, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	static thread_local vector<double> Position_Data_param(MAX_DIMENSION);
	static thread_local vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);
	static thread_local vector<double> dShape(info->DIMENSION * max_support_1D);

	int patch = info->Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		ShapeFunc_from_paren(Position_Data_param.data(), Local_coord, i, e, info);
		shape_function_1D(Position_Data_param[i], i, e, Shape.data(), dShape.data(), info);
		basis_num_offset[i] = info->ENC[e * info->DIMENSION + i];
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
		{
			int basis_num = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + i] - basis_num_offset[i];
			shape_func[j] *= Shape[(i * MAX_ORDER + order[i]) * max_support_1D + basis_num];
		}
	}

	// make R
	double weight_func = 0.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		weight_func += shape_func[i] * 1.0;
	for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		R[i] = shape_func[i] * 1.0 / weight_func;

	// make dR
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;
	double dWeight_func3 = 0.0;
	int basis_num_dim[MAX_DIMENSION];
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = 1.0;
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * w_coord;
		}
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = 1.0;
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = 1.0;
			dWeight_func1 += dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func2 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] * w_coord;
			dWeight_func3 += Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] * w_coord;
		}

		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
				basis_num_dim[j] = info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] - basis_num_offset[j];

			double w_coord = 1.0;
			dR[i * info->DIMENSION + 0] = w_coord * (weight_func * dShape[0 * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 1] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * dShape[1 * max_support_1D + basis_num_dim[1]] * Shape[(2 * MAX_ORDER + order[2]) * max_support_1D + basis_num_dim[2]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
			dR[i * info->DIMENSION + 2] = w_coord * (weight_func * Shape[(0 * MAX_ORDER + order[0]) * max_support_1D + basis_num_dim[0]] * Shape[(1 * MAX_ORDER + order[1]) * max_support_1D + basis_num_dim[1]] * dShape[2 * max_support_1D + basis_num_dim[2]] - dWeight_func3 * shape_func[i]) / (weight_func * weight_func);
		}
	}

	if (calc_dR_coef_flag)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			double coef = dShapeFunc_from_paren(i, e, info);
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[patch]; j++)
				dR[j * info->DIMENSION + i] *= coef;
		}
	}
}


int calc_patch_parameter_coord(double *physical_coord, int patch_num, double *out_para_coord, information *info)
{
	// constexpr int max_itr = 60;
	constexpr int max_itr = 20;
	constexpr double tol = MERGE_ERROR;

	static thread_local vector<double> delta_current_coord(MAX_DIMENSION);
	static thread_local vector<double> local_coord(MAX_DIMENSION);
	static thread_local vector<double> M(MAX_DIMENSION * MAX_DIMENSION);
	static thread_local vector<double> R(MAX_NO_CP_ON_ELEMENT);
	static thread_local vector<double> dR(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);

	// init
	for (int i = 0; i < info->DIMENSION; i++)
		out_para_coord[i] = 0.5 * (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]]);

	for (int itr = 0; itr < max_itr; itr++)
	{
		// search element
		if (check_nonfinite(out_para_coord, info->DIMENSION))
			return ERROR;
		int e = ele_check(patch_num, out_para_coord, info);
		tilde_coord(local_coord.data(), out_para_coord, patch_num, e, info);

		// calc shape function and its derivative
		shape_and_dshape(R.data(), dR.data(), local_coord.data(), e, false, info);

		// delta_current_coord x - x_current
		for (int i = 0; i < info->DIMENSION; i++)
			delta_current_coord[i] = physical_coord[i];
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[patch_num]; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				delta_current_coord[j] -= R[i] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];

		// check convergence
		double residual = 0.0;
		for (int i = 0; i < info->DIMENSION; i++)
			residual += delta_current_coord[i] * delta_current_coord[i];
		residual = sqrt(residual);
		if (residual < tol)
			return itr;
		if (!isfinite(residual))
			return ERROR;

		// make matrix
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				M[i * info->DIMENSION + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[patch_num]; k++)
					M[i * info->DIMENSION + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}

		double det = 0.0;
		if (info->DIMENSION == 2)
			det = Determinant_2x2(M.data());
		else if (info->DIMENSION == 3)
			det = Determinant_3x3(M.data());
		if (!isfinite(det) || fabs(det) < tol)
			continue;

		// calc inverse
		if (info->DIMENSION == 2)
			InverseMatrix_2x2(M.data());
		else if (info->DIMENSION == 3)
			InverseMatrix_3x3(M.data());

		// solve and update
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				out_para_coord[i] += M[i * info->DIMENSION + j] * delta_current_coord[j];

		// modify result
		for (int i = 0; i < info->DIMENSION; i++)
		{
			if (out_para_coord[i] > info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]];
			if (out_para_coord[i] < info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]];
		}
	}

	return ERROR;
}


// ローカルパラメータ座標から、グローバルのパラメータ座標を得る関数
int calc_global_patch_parameter_coord(double *local_para_coord, int patch_num, double *out_para_coord, information *info)
{
	constexpr int max_itr = 20;
	constexpr double tol = MERGE_ERROR;

	vector<double> diff(MAX_DIMENSION);
	vector<double> M(MAX_DIMENSION * MAX_DIMENSION);

	// make M 
	for (int i = 0; i < info->DIMENSION; i++)
	    for (int j = 0; j < info->DIMENSION; j++)
	        M[i * info->DIMENSION + j] = (i == j) ? 1.0 : 0.0;

	// init
	for (int i = 0; i < info->DIMENSION; i++)
		out_para_coord[i] = 0.5 * (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]]);

	for (int itr = 0; itr < max_itr; itr++)
	{
		// delta_current_coord x - x_current
		for (int i = 0; i < info->DIMENSION; i++)
			diff[i] = local_para_coord[i] - out_para_coord[i];

		// check convergence
		double residual = 0.0;
		for (int i = 0; i < info->DIMENSION; i++)
			residual += diff[i] * diff[i];
		residual = sqrt(residual);


		if (residual < tol)
			return itr;

		// calc inverse
		if (info->DIMENSION == 2)
			InverseMatrix_2x2(M.data());
		else if (info->DIMENSION == 3)
			InverseMatrix_3x3(M.data());

		// solve and update
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				out_para_coord[i] += M[i * info->DIMENSION + j] * diff[j];

		// modify result
		for (int i = 0; i < info->DIMENSION; i++)
		{
			if (out_para_coord[i] > info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]];
			if (out_para_coord[i] < info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]];
		}
	}

	return ERROR;
}


// グローバルパラメータ座標から、ローカルのパラメータ座標を得る関数
int calc_local_patch_parameter_coord(double *global_para_coord, int patch_num, double *out_para_coord, information *info)
{
	constexpr int max_itr = 20;
	constexpr double tol = MERGE_ERROR;

	vector<double> delta_current_coord(MAX_DIMENSION);
	vector<double> local_coord(MAX_DIMENSION);
	vector<double> M(MAX_DIMENSION * MAX_DIMENSION);
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);

	// init
	int geo_patch_num = patch_num - info->Total_Patch_on_mesh[0];
	for (int i = 0; i < info->DIMENSION; i++)
		out_para_coord[i] = 0.5 * (info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i] + info->Geo_No_Control_point[geo_patch_num * info->DIMENSION + i] + info->Geo_Order[geo_patch_num * info->DIMENSION + i]] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i]]);

	for (int itr = 0; itr < max_itr; itr++)
	{
		// search element
		int geo_e = geo_ele_check(geo_patch_num, out_para_coord, info);
		geo_tilde_coord(local_coord.data(), out_para_coord, geo_patch_num, geo_e, info);

		#if 0
		// 特異点チェック：頂点座標（±1.0）での特異点回避
		double modified_local_coord[MAX_DIMENSION];
		bool needs_modification = false;
		
		for (int dim = 0; dim < info->DIMENSION; dim++) {
			if (fabs(local_coord[dim] - 1.0) < 1e-15) {
				modified_local_coord[dim] = 1.0 - 1e-15;  // わずかに内側にシフト
				needs_modification = true;
			} else if (fabs(local_coord[dim] + 1.0) < 1e-15) {
				modified_local_coord[dim] = -1.0 + 1e-15;  // わずかに内側にシフト
				needs_modification = true;
			} else {
				modified_local_coord[dim] = local_coord[dim];
			}
		}

		// 修正が必要な場合は修正座標を使用
		if (needs_modification) {
			for (int dim = 0; dim < info->DIMENSION; dim++)
				local_coord[dim] = modified_local_coord[dim];
		}
		#endif

		// calc shape function and its derivative
		geo_shape_and_dshape(R.data(), dR.data(), local_coord.data(), geo_e, false, info);

		// delta_current_coord x - x_current
		for (int i = 0; i < info->DIMENSION; i++)
			delta_current_coord[i] = global_para_coord[i];
		for (int i = 0; i < info->Geo_No_Control_point_ON_ELEMENT[geo_patch_num]; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				delta_current_coord[j] -= R[i] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_e * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];

		// check convergence
		double residual = 0.0;
		for (int i = 0; i < info->DIMENSION; i++)
			residual += delta_current_coord[i] * delta_current_coord[i];
		residual = sqrt(residual);
		if (residual < tol)
			return itr;

		// make matrix
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				M[i * info->DIMENSION + j] = 0.0;
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[geo_patch_num]; k++)
					M[i * info->DIMENSION + j] += dR[k * info->DIMENSION + j] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[geo_e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}

		double det = 0.0;
		if (info->DIMENSION == 2)
			det = Determinant_2x2(M.data());
		else if (info->DIMENSION == 3)
			det = Determinant_3x3(M.data());
		if (!isfinite(det) || fabs(det) < tol)
			continue;

		// calc inverse
		if (info->DIMENSION == 2)
			InverseMatrix_2x2(M.data());
		else if (info->DIMENSION == 3)
			InverseMatrix_3x3(M.data());

		// solve and update
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
				out_para_coord[i] += M[i * info->DIMENSION + j] * delta_current_coord[j];

		// modify result
		for (int i = 0; i < info->DIMENSION; i++)
		{
			if (out_para_coord[i] > info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i] + info->Geo_No_Control_point[geo_patch_num * info->DIMENSION + i] + info->Geo_Order[geo_patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i] + info->Geo_No_Control_point[geo_patch_num * info->DIMENSION + i] + info->Geo_Order[geo_patch_num * info->DIMENSION + i]];
			if (out_para_coord[i] < info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i]])
				out_para_coord[i] = info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[geo_patch_num * info->DIMENSION + i]];
		}
	}

	return ERROR;
}


bool check_nonfinite(double *vec, int size)
{
	for (int i = 0; i < size; i++)
	{
		if (!isfinite(vec[i]))
			return true;
	}
	return false;
}