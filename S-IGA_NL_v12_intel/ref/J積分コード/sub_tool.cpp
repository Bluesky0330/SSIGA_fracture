#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

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

	static double *Position_Data_param_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *shape_func_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *Shape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_max_threads());
	static double *dShape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_max_threads());
	double *Position_Data_param = &Position_Data_param_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *shape_func = &shape_func_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *Shape = &Shape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_thread_num()];
	double *dShape = &dShape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_thread_num()];

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
	static double *dShape_func1_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *dShape_func2_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *dShape_func3_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	double *dShape_func1 = &dShape_func1_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *dShape_func2 = &dShape_func2_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *dShape_func3 = &dShape_func3_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
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

	static double *Position_Data_param_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *shape_func_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *Shape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_max_threads());
	static double *dShape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_max_threads());
	double *Position_Data_param = &Position_Data_param_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *shape_func = &shape_func_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *Shape = &Shape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_thread_num()];
	double *dShape = &dShape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_thread_num()];

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

	static double *Position_Data_param_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *shape_func_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *Shape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_max_threads());
	static double *dShape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_max_threads());
	double *Position_Data_param = &Position_Data_param_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *shape_func = &shape_func_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *Shape = &Shape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER * omp_get_thread_num()];
	double *dShape = &dShape_parallel[info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * omp_get_thread_num()];

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


// Newton-Raphson法, from NURBSviewer
double &BasisFunc(double *knot_vec, const int &knot_index, const int &order, const double &xi, double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	static double *temp_basis = (double *)malloc(sizeof(double) * MAX_ORDER * MAX_ORDER);
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && knot_vec[knot_index + order + 1] >= xi)
	{
		for (j = 0; j <= order; j++)
		{
			if ((knot_vec[knot_index + j] <= xi) && (xi <= knot_vec[knot_index + j + 1]))
				temp_basis[j * MAX_ORDER + 0] = 1.0;
			else
				temp_basis[j * MAX_ORDER + 0] = 0.0;
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j * MAX_ORDER + p - 1];
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[(j + 1) * MAX_ORDER + p - 1];
					temp_basis[j * MAX_ORDER + p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0 * MAX_ORDER + order - 1];
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1 * MAX_ORDER + order - 1];
		}
		(*output) = temp_basis[0 * MAX_ORDER + order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double &rBasisFunc(double *knot_vec, const int &knot_index, const int &order, const double &xi, double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	static double *temp_basis = (double *)malloc(sizeof(double) * MAX_ORDER * MAX_ORDER);
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1])
	{
		if (knot_index == 0)
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[j] <= xi) && (xi <= knot_vec[j + 1]))
					temp_basis[j * MAX_ORDER + 0] = 1.0;
				else
					temp_basis[j * MAX_ORDER + 0] = 0.0;
			}
		}
		else
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[knot_index + j] < xi) && (xi <= knot_vec[knot_index + j + 1]))
					temp_basis[j * MAX_ORDER + 0] = 1.0;
				else
					temp_basis[j * MAX_ORDER + 0] = 0.0;
			}
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j * MAX_ORDER + p - 1];
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[(j + 1) * MAX_ORDER + p - 1];
					temp_basis[j * MAX_ORDER + p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;

			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0 * MAX_ORDER + order - 1];
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1 * MAX_ORDER + order - 1];
		}
		(*output) = temp_basis[0 * MAX_ORDER + order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double &lBasisFunc(double *knot_vec, const int &knot_index, const int &cntl_p_n, const int &order, const double &xi, double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	static double *temp_basis = (double *)malloc(sizeof(double) * MAX_ORDER * MAX_ORDER);
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1])
	{
		if (knot_index == cntl_p_n - 1)
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[cntl_p_n - 1 + j] <= xi) && (xi <= knot_vec[cntl_p_n + j]))
					temp_basis[j * MAX_ORDER + 0] = 1.0;
				else
					temp_basis[j * MAX_ORDER + 0] = 0.0;
			}
		}
		else
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[knot_index + j] <= xi) && (xi < knot_vec[knot_index + j + 1]))
					temp_basis[j * MAX_ORDER + 0] = 1.0;
				else
					temp_basis[j * MAX_ORDER + 0] = 0.0;
			}
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j * MAX_ORDER + p - 1];
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[(j + 1) * MAX_ORDER + p - 1];
					temp_basis[j * MAX_ORDER + p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0 * MAX_ORDER + order - 1];
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1 * MAX_ORDER + order - 1];
		}
		(*output) = temp_basis[0 * MAX_ORDER + order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double NURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					 double *cntl_px, double *cntl_py,
					 int cntl_p_n_xi, int cntl_p_n_eta,
					 double *weight, int order_xi, int order_eta,
					 double xi, double eta,
					 double *output_x, double *output_y,
					 double *output_dxi_x, double *output_deta_x,
					 double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	// int index_max_xi = cntl_p_n_xi; //2020_09_12
	int index_max_xi = cntl_p_n_xi - 1; // 2020_09_12
	int index_min_eta = 0;
	// int index_max_eta = cntl_p_n_eta; //2020_09_12
	int index_max_eta = cntl_p_n_eta - 1; // 2020_09_12

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] > xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi; // 2020_09_12

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] > eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta; // 2020_09_12

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		BasisFunc(input_knot_vec_xi, i, order_xi, xi,
				  &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			BasisFunc(input_knot_vec_eta, j, order_eta, eta,
					  &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double rrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + (j * cntl_p_n_xi) + (k * cntl_p_n_xi * cntl_p_n_eta);

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rlrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double llrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_deta_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_deta_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}
	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}
	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}
	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[temp_index];
				deta_molecule_y += temp3 * cntl_py[temp_index];
				deta_molecule_z += temp3 * cntl_pz[temp_index];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_deta_x)  = weight[temp_index] * (deta_molecule_x  * denominator - molecule_x * deta_denominator)  / (denominator * denominator);
	(*output_deta_y)  = weight[temp_index] * (deta_molecule_y  * denominator - molecule_y * deta_denominator)  / (denominator * denominator);
	(*output_deta_z)  = weight[temp_index] * (deta_molecule_z  * denominator - molecule_z * deta_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


// 算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
int Calc_xi_eta(double px, double py,
				double *input_knot_vec_xi, double *input_knot_vec_eta,
				int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
				double *output_xi, double *output_eta,
				int global_patch_num, information *info)
{
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x, temp_tol_y;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// int repeat = 1000;
	// double tol = 10e-8;
	int repeat = 20;
	double tol = 10e-14;

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	double *ptr_x = &(info->Control_Coord_x[info->Total_Control_Point_to_patch[global_patch_num]]);
    double *ptr_y = &(info->Control_Coord_y[info->Total_Control_Point_to_patch[global_patch_num]]);
    double *ptr_w = &(info->Control_Weight[info->Total_Control_Point_to_patch[global_patch_num]]);

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
					   ptr_x, ptr_y,
					   cntl_p_n_xi, cntl_p_n_eta,
					   ptr_w, order_xi, order_eta,
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
					   ptr_x, ptr_y,
					   cntl_p_n_xi, cntl_p_n_eta,
					   ptr_w, order_xi, order_eta,
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		rlNURBS_surface(input_knot_vec_eta, input_knot_vec_eta,
						ptr_x, ptr_y,
						cntl_p_n_xi, cntl_p_n_eta,
						ptr_w, order_xi, order_eta,
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						ptr_x, ptr_y,
						cntl_p_n_xi, cntl_p_n_eta,
						ptr_w, order_xi, order_eta,
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
	}

	return ERROR;
}


int Calc_xi_eta_zeta(double px, double py, double pz,
				     double *input_knot_vec_xi, double *input_knot_vec_eta, double *input_knot_vec_zeta,
				     int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                     int order_xi, int order_eta, int order_zeta,
				     double *output_xi, double *output_eta, double *output_zeta, int global_patch_num, information *info)
{
	double temp_xi, temp_eta, temp_zeta;
	double temp_x, temp_y, temp_z;
	double temp_matrix[3][3];
	double temp_dxi, temp_deta, temp_dzeta;
    double temp_tol_x, temp_tol_y, temp_tol_z;

	(*output_xi)   = 0.0;
	(*output_eta)  = 0.0;
	(*output_zeta) = 0.0;

	int i, j;
	int repeat  = 5;
	int repeat2 = 5;
	// double tol = 10e-22;

	double tol = 10e-14;
	double coef = 0.0;

	double *ptr_x = &(info->Control_Coord_x[info->Total_Control_Point_to_patch[global_patch_num]]);
    double *ptr_y = &(info->Control_Coord_y[info->Total_Control_Point_to_patch[global_patch_num]]);
    double *ptr_z = &(info->Control_Coord_z[info->Total_Control_Point_to_patch[global_patch_num]]);
    double *ptr_w = &(info->Control_Weight[info->Total_Control_Point_to_patch[global_patch_num]]);

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rrrNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lrrNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rlrNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rrlNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			llrNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lrlNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rllNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lllNURBS_volume(input_knot_vec_xi, input_knot_vec_eta, input_knot_vec_zeta,
							ptr_x, ptr_y, ptr_z,
							cntl_p_n_xi, cntl_p_n_eta, cntl_p_n_zeta,
							ptr_w, order_xi, order_eta, order_zeta,
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	return ERROR;
}


void shape_function_1D(double knot, int dim, int e, double *Shape, double *dShape, int select_deriv, information *info)
{
	static int max_support_1D = 2 * MAX_ORDER + 2;
	int patch = info->Element_patch[e];
	int offset = info->ENC[e * info->DIMENSION + dim];
	int p = info->Order[patch * info->DIMENSION + dim];
	int knot_n = info->No_knot[patch * info->DIMENSION + dim];
	double *knot_ptr = &info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + dim]];
	int support = p + 1;  // 基底関数のサポート範囲

	// 0 次基底関数
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

	// 1 次以上の基底関数
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
	if (select_deriv == 0)
		return;

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


#if 1
void shape_and_dshape(double *R, double *dR, double *Local_coord, int e, information *info)
{
	static int max_support_1D = 2 * MAX_ORDER + 2;
	static double *Position_Data_param_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *shape_func_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *Shape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_ORDER * max_support_1D * omp_get_max_threads());
	static double *dShape_parallel = (double *)malloc(sizeof(double) * info->DIMENSION * max_support_1D * omp_get_max_threads());
	// static double *dShape_func1_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	// static double *dShape_func2_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	// static double *dShape_func3_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());

	double *Position_Data_param = &Position_Data_param_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *shape_func = &shape_func_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *Shape = &Shape_parallel[info->DIMENSION * MAX_ORDER * max_support_1D * omp_get_thread_num()];
	double *dShape = &dShape_parallel[info->DIMENSION * max_support_1D * omp_get_thread_num()];
	// double *dShape_func1 = &dShape_func1_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	// double *dShape_func2 = &dShape_func2_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	// double *dShape_func3 = &dShape_func3_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];

	int patch = info->Element_patch[e];
	int order[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
		order[i] = info->Order[patch * info->DIMENSION + i];

	for (int i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
		shape_func[i] = 1.0;

	int basis_num_offset[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, i, e, info);
		shape_function_1D(Position_Data_param[i], i, e, Shape, dShape, 1, info);
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
}
#endif


int calc_patch_parameter_coord(double *physical_coord, int patch_num, double *out_para_coord, information *info)
{
	static const int max_itr = 30;
	static const double tol = MERGE_ERROR_PARA_COORD;

	static double *delta_current_coord_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *local_coord_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * omp_get_max_threads());
	static double *M_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_DIMENSION * omp_get_max_threads());
	static double *R_parallel = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());
	static double *dR_parallel = (double *)malloc(sizeof(double) * MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_max_threads());

	double *delta_current_coord = &delta_current_coord_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *local_coord = &local_coord_parallel[MAX_DIMENSION * omp_get_thread_num()];
	double *M = &M_parallel[MAX_DIMENSION * MAX_DIMENSION * omp_get_thread_num()];
	double *R = &R_parallel[MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];
	double *dR = &dR_parallel[MAX_DIMENSION * MAX_NO_CP_ON_ELEMENT * omp_get_thread_num()];

	// init
	for (int i = 0; i < info->DIMENSION; i++)
		out_para_coord[i] = 0.5 * (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->No_Control_point[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i]] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i]]);

	for (int itr = 0; itr < max_itr; itr++)
	{
		// search element
		int e = ele_check(patch_num, out_para_coord, info);
		tilde_coord(local_coord, out_para_coord, patch_num, e, info);

		// calc shape function and its derivative
		shape_and_dshape(R, dR, local_coord, e, info);

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

		// make matrix
		for (int i = 0; i < info->DIMENSION; i++)
			for (int j = 0; j < info->DIMENSION; j++)
			{
				M[i * info->DIMENSION + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[patch_num]; k++)
					M[i * info->DIMENSION + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
			}

		// calc inverse
		if (info->DIMENSION == 2)
			InverseMatrix_2x2(M);
		else if (info->DIMENSION == 3)
			InverseMatrix_3x3(M);

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