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


int duplicate(int total_n, int *array)
{
	int i, j;

	j = 0;
	j++;
	for (i = 1; i < total_n; i++)
	{
		if (array[i] != array[i - 1])
		{
			array[j] = array[i];
			j++;
		}
	}

	return j;
}


// Shape Function
double Shape_func(int I_No, double *Local_coord, int El_No, information *info)
{
	int i, j;
	double R, weight_func = 0.0;

	thread_local double *Position_Data_param = (double *)malloc(sizeof(double) * MAX_DIMENSION);
	thread_local double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	thread_local double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER);
	thread_local double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, 0, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

	if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		R = shape_func[I_No] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + I_No] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;
	else
		R = ERROR;

	return R;
}


void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	Position_Data_param[j] = ((info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) * Local_coord[j] + (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] + info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i])) / 2.0;
}


void ShapeFunction1D(double *Position_Data_param, int j, int e, double *Shape, double *dShape, int select_deriv, information *info)
{
	int ii, p;
	for (ii = 0; ii < info->No_knot[info->Element_patch[e] * info->DIMENSION + j]; ii++)
	{
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
		else if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] != info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] < info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] != info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + (info->No_knot[info->Element_patch[e] * info->DIMENSION + j] - 1)] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] <= info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
	}

	double left_term, right_term;
	for (p = 1; p <= info->Order[info->Element_patch[e] * info->DIMENSION + j]; p++)
	{
		for (ii = 0; ii < info->No_knot[info->Element_patch[e] * info->DIMENSION + j]; ii++)
		{
			left_term = 0.0;
			right_term = 0.0;

			if ((Position_Data_param[j] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == 0)
				left_term = 0.0;
			else
				left_term = (Position_Data_param[j] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1];

			if ((info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == 0)
				right_term = 0.0;
			else
				right_term = (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1];

			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = left_term + right_term;
		}
	}

	if (select_deriv == 0)
		return;
	else
	{
		double dleft_term, dright_term;
		for (ii = 0; ii < info->No_Control_point[info->Element_patch[e] * info->DIMENSION + j] + 1; ii++)
		{
			dleft_term = 0.0;
			dright_term = 0.0;

			if (info->Order[info->Element_patch[e] * info->DIMENSION + j] * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j]] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == 0)
				dleft_term = 0.0;
			else
				dleft_term = info->Order[info->Element_patch[e] * info->DIMENSION + j] / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j]] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1];

			if (info->Order[info->Element_patch[e] * info->DIMENSION + j] * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == 0)
				dright_term = 0.0;
			else
				dright_term = info->Order[info->Element_patch[e] * info->DIMENSION + j] / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1];

			dShape[j * info->Total_Control_Point_to_mesh[Total_mesh] + ii] = dleft_term - dright_term;
		}
	}
}


double dShape_func(int I_No, int xez, double *Local_coord, int El_No, information *info)
{
	thread_local double *dShape_func1 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	thread_local double *dShape_func2 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	thread_local double *dShape_func3 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);

	double dR = 0.0;

	if (info->DIMENSION == 2)
	{
		NURBS_deriv_2D(Local_coord, El_No, dShape_func1, dShape_func2, info);

		if (xez != 0 && xez != 1)
			dR = ERROR;
		else if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		{
			if (xez == 0)
				dR = dShape_func1[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			else if (xez == 1)
				dR = dShape_func2[I_No] * dShapeFunc_from_paren(xez, El_No, info);
		}
		else
			dR = ERROR;
	}
	else if (info->DIMENSION == 3)
	{
		NURBS_deriv_3D(Local_coord, El_No, dShape_func1, dShape_func2, dShape_func3, info);

		if (xez != 0 && xez != 1 && xez != 2)
			dR = ERROR;
		else if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		{
			if (xez == 0)
				dR = dShape_func1[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			else if (xez == 1)
				dR = dShape_func2[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			else if (xez == 2)
				dR = dShape_func3[I_No] * dShapeFunc_from_paren(xez, El_No, info);
		}
		else
			dR = ERROR;
	}

	return dR;
}


void NURBS_deriv_2D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, information *info)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;

	thread_local double *Position_Data_param = (double *)malloc(sizeof(double) * MAX_DIMENSION);
	thread_local double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	thread_local double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER);
	thread_local double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, 1, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
	}

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]] * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
		dWeight_func2 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]] * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dShape_func1[i] = info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] * (weight_func * dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]] * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
		dShape_func2[i] = info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] * (weight_func * Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]] * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
	}
}


void NURBS_deriv_3D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, double *dShape_func3, information *info)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;
	double dWeight_func3 = 0.0;

	thread_local double *Position_Data_param = (double *)malloc(sizeof(double) * MAX_DIMENSION);
	thread_local double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	thread_local double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER);
	thread_local double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, 1, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
	}

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]]
					   * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
					   * Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

		dWeight_func2 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
					   * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]]
					   * Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

		dWeight_func3 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
					   * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
					   * dShape[2 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}

	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dShape_func1[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]]
			* Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
			* Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
			- dWeight_func1	* shape_func[i])
			/ (weight_func * weight_func);

		dShape_func2[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
			* dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]]
			* Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
			- dWeight_func2	* shape_func[i])
			/ (weight_func * weight_func);

		dShape_func3[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
			* Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
			* dShape[2 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2]]
			- dWeight_func3 * shape_func[i])
			/ (weight_func * weight_func);
	}
}


double dShapeFunc_from_paren(int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	double dPosition_Data_param = (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) / 2.0;
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


void shape_and_dshape(double *R, double *Local_coord, int e, information *info)
{
	int max_support_1D = 2 * MAX_ORDER + 2;
	vector<double> Position_Data_param(MAX_DIMENSION);
	vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);

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
	vector<double> Position_Data_param(MAX_DIMENSION);
	vector<double> shape_func(MAX_NO_CP_ON_ELEMENT);
	vector<double> Shape(info->DIMENSION * MAX_ORDER * max_support_1D);
	vector<double> dShape(info->DIMENSION * max_support_1D);

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


int calc_patch_parameter_coord(double *physical_coord, int patch_num, double *out_para_coord, information *info)
{
	// constexpr int max_itr = 60;
	constexpr int max_itr = 20;
	constexpr double tol = MERGE_ERROR;

	vector<double> delta_current_coord(MAX_DIMENSION);
	vector<double> local_coord(MAX_DIMENSION);
	vector<double> M(MAX_DIMENSION * MAX_DIMENSION);
	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT);

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
	// Fast Math最適化対策: residualの異常値チェック
    if (!is_value_finite_robust(residual))
        return ERROR;
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				M[i * info->DIMENSION + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[patch_num]; k++)
					M[i * info->DIMENSION + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}

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


bool check_nonfinite(double *vec, int size)
{
	for (int i = 0; i < size; i++)
	{
		// Fast Math最適化対策: vec[i]の異常値チェック
		volatile double val_vol = vec[i];
		if (!isfinite(val_vol) || vec[i] != vec[i] || (vec[i] > 0 && vec[i] > 1.0e100) || (vec[i] < 0 && vec[i] < -1.0e100))
			return true;
	}
	return false;
}

// Fast Math最適化下でも堅牢な単一値の有限性チェック
// 戻り値: true = 有限値, false = NaN/Inf/異常値
bool is_value_finite_robust(double value)
{
	// Fast Math最適化対策: volatile宣言で最適化を抑制
	volatile double val_vol = value;
	
	// 複数の判定方法を組み合わせて堅牢性を確保
	// 1. isfinite (Fast Mathでは機能しない可能性)
	// 2. 自己不等号 (NaN != NaN は常にtrue)
	// 3. 範囲チェック (Infを補完的に検出)
	if (!isfinite(val_vol) || value != value || 
	    (value > 0 && value > 1.0e100) || (value < 0 && value < -1.0e100))
		return false;
	
	return true;
}