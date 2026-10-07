// header
#include "_header.hpp"
#include "_sub.hpp"

#if 0
	// windows, cygwin
	#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
		#include <Eigen/Dense>
		#include <Eigen/Sparse>
	#endif

	// linux
	#if defined(__linux__)
		#include <eigen3/Eigen/Dense>
		#include <eigen3/Eigen/Sparse>
	#endif

	using namespace Eigen;
#endif

using namespace std;


// Print progress bar (using \r to overwrite)
// Intermediate updates go to stderr so typical stdout logs keep only final 100% line.
static void Print_Solver_Progress(const char *label, int current, int total)
{
	if (total <= 0)
		total = 1;
	if (current < 0)
		current = 0;
	if (current > total)
		current = total;

	const int BAR_WIDTH = 30;
	int percent = (100 * current) / total;
	int filled = (percent * BAR_WIDTH) / 100;
	FILE *out = (percent >= 100) ? stdout : stderr;

	fprintf(out, "\r[%s]: %3d%% || ", label, percent);
	for (int i = 0; i < BAR_WIDTH; i++)
	{
		if (i < filled)
			fprintf(out, "█");
		else
			fprintf(out, " ");
	}

	if (percent >= 100)
		fprintf(out, "|| (%d/%d iterations)\n", current, total);
	else
		fprintf(out, "|| (%d/%d iterations)                    ", current, total);

	fflush(out);
}

// PCG solver, 前処理付共役勾配法(PCG)により[K]{d}={f}を解く, 対角スケーリングを行ったCG法を用いる
void PCG_Solver(int max_itetarion, double eps, information *info)
{
	int i, j, k;
	int ndof = K_Whole_Size;
	double alpha, beta, e = 0.0;
	int M_mode = 2; // 0: M = K, 1: M = グローバルパッチとローカルパッチのKの対角バンド成分, 2: M = [[K^G, 0], [0, K^L]]

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));

	double *gg = (double *)malloc(sizeof(double) * ndof);
	double *dd = (double *)malloc(sizeof(double) * ndof);
	double *pp = (double *)malloc(sizeof(double) * ndof);

	double *temp_r = (double *)malloc(sizeof(double) * ndof);

	// 初期化
	for (i = 0; i < ndof; i++)
		info->sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)malloc(sizeof(double) * info->K_Whole_Ptr[ndof]);
	int *M_Ptr = (int *)malloc(sizeof(int) * (ndof + 1));
	int *M_Col = (int *)malloc(sizeof(int) * info->K_Whole_Ptr[ndof]);
	double *M_diag = (double *)malloc(sizeof(double) * ndof);
	Make_M(M_mode, M, M_Ptr, M_Col, M_diag, ndof, info);

	// 第0近似解に対する残差の計算
	for (i = 0; i < ndof; i++)
		r[i] = info->rhs_vec[i];

	// p_0 = (LDL^T)^-1 r_0 の計算 <- CG法で M = [[K^G, 0], [0, K^L]] とし, p_0 = (LDL^T)^-1 r_0 = M^-1 r_0
	CG(ndof, p, M, M_Ptr, M_Col, M_diag, r, gg, dd, pp, temp_r);

	for (k = 0; k < max_itetarion; k++)
	{
		// rr0 の計算
		double rr0 = inner_product(ndof, r, p);

		// y = AP の計算
		for (i = 0; i < ndof; i++)
		{
			y[i] = 0.0;
			for (j = 0; j < i; j++)
			{
				int temp1 = RowCol_to_icount(j, i, info); // temp_array_K[i][j] = temp_array_K[j][i]
				if (temp1 != -1)
					y[i] += info->K_Whole_Val[temp1] * p[j];
			}
			for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
				y[i] += info->K_Whole_Val[j] * p[info->K_Whole_Col[j]];
		}

		// alpha = r * r / (P * AP) の計算
		double temp_scaler = inner_product(ndof, p, y);
		alpha = rr0 / temp_scaler;

		// 解x, 残差rの更新
		for (i = 0; i < ndof; i++)
		{
			info->sol_vec[i] += alpha * p[i];
			r[i] -= alpha * y[i];
		}

		// (r * r)_(k+1)の計算
		CG(ndof, r2, M, M_Ptr, M_Col, M_diag, r, gg, dd, pp, temp_r);

		// 収束判定 (CG法と同じ)
		double e1 = inner_product(ndof, p, p), e2 = inner_product(ndof, info->sol_vec, info->sol_vec);
		e = fabs(alpha) * sqrt(e1 / e2);
		if (e < eps)
		{
			k++;
			break;
		}

		// βの計算とPの更新
		beta = inner_product(ndof, y, r2) / temp_scaler;

		for (i = 0; i < ndof; i++)
			p[i] = r2[i] - beta * p[i];

		printf("PCG itr %d\t", k);
		printf("eps %le", e);
	}

	int max_itr_result = k;
	double eps_result = e;

	printf("\nndof = %d\n", ndof);
	printf("itr_result = %d\n", max_itr_result);
	printf("eps_result = %.15e\n", eps_result);

	free(r), free(p), free(y), free(r2);
	free(M), free(M_Ptr), free(M_Col), free(M_diag);
	free(gg), free(dd), free(pp), free(temp_r);
}


void Make_M(int M_mode, double *M, int *M_Ptr, int *M_Col, double *M_diag, int ndof, information *info)
{
	int i, j;
	int ndof_glo = 0, ndof_loc;
	int counter = 0;

	// グローバルパッチのdofを求める
	for (i = 0; i < info->Total_Control_Point_on_mesh[0] * info->DIMENSION; i++)
	{
		if (info->Index_Dof[i] != ERROR)
		{
			ndof_glo++;
		}
	}
	ndof_loc = ndof - ndof_glo;
	printf("ndof		%d\n", ndof);
	printf("ndof_glo	%d\n", ndof_glo);
	printf("ndof_loc	%d\n", ndof_loc);

	if (M_mode == 0)
	{
		// M = K と M_diag を作成
		M_Ptr[0] = 0;
		for (i = 0; i < ndof; i++)
		{
			M_Ptr[i + 1] = M_Ptr[i];
			for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
			{
				M[j] = info->K_Whole_Val[j];
				M_Col[j] = info->K_Whole_Col[j];
				M_Ptr[i + 1]++;
			}
		}
	}
	else if (M_mode == 1)
	{
		// M = グローバルパッチとローカルパッチのKの対角バンド成分 と M_diag を作成
		int band_width[2];
		band_width[0] = info->No_Control_point_in_patch[0] * info->DIMENSION; 								// global patch 0 の自由度数
		band_width[1] = info->No_Control_point_in_patch[info->Total_Patch_on_mesh[0]] * info->DIMENSION;	// local patch 0 の自由度数
		if (ndof_glo < band_width[0])
			band_width[0] = ndof_glo;
		if (ndof_loc < band_width[1])
			band_width[1] = ndof_loc;
		printf("global: band width / ndof_glo = %d / %d\n", band_width[0], ndof_glo);
		printf("local : band width / ndof_loc = %d / %d\n", band_width[1], ndof_loc);

		M_Ptr[0] = 0;
		for (i = 0; i < ndof; i++)
		{
			M_Ptr[i + 1] = M_Ptr[i];
			for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
			{
				if (i < ndof_glo && info->K_Whole_Col[j] < ndof_glo)
				{
					if (info->K_Whole_Col[j] - i == band_width[0])
						break;
					M[counter] = info->K_Whole_Val[j];
					M_Col[counter] = info->K_Whole_Col[j];
					counter++;
					M_Ptr[i + 1]++;
				}
				else if (i >= ndof_glo)
				{
					if (info->K_Whole_Col[j] - i == band_width[1])
						break;
					M[counter] = info->K_Whole_Val[j];
					M_Col[counter] = info->K_Whole_Col[j];
					counter++;
					M_Ptr[i + 1]++;
				}
			}
		}
	}
	else if (M_mode == 2)
	{
		// M = [[K^G, 0], [0, K^L]] と M_diag を作成
		M_Ptr[0] = 0;
		for (i = 0; i < ndof; i++)
		{
			M_Ptr[i + 1] = M_Ptr[i];
			for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
			{
				if (i < ndof_glo && info->K_Whole_Col[j] < ndof_glo)
				{
					M[counter] = info->K_Whole_Val[j];
					M_Col[counter] = info->K_Whole_Col[j];
					counter++;
					M_Ptr[i + 1]++;
				}
				else if (i >= ndof_glo)
				{
					M[counter] = info->K_Whole_Val[j];
					M_Col[counter] = info->K_Whole_Col[j];
					counter++;
					M_Ptr[i + 1]++;
				}
			}
		}
	}

	// diag scaling preprocess
	M_diag[0] = 1.0 / sqrt(M[0]);
	for (i = 1; i < ndof; i++)
		M_diag[i] = 1.0 / sqrt(M[M_Ptr[i]]);

	counter = 0;
	for (i = 0; i < ndof; i++)
	{
		for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
		{
			M[counter] = M[counter] * M_diag[i] * M_diag[M_Col[j]];
			counter++;
		}
	}
}


void CG(const int &ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *M_diag, double *right_vec, double *gg, double *dd, double *pp, double *temp_r)
{
	int i, j, itr;
	int max_itr = ndof;
	double eps = 1.0e-5;
	double rrr3 = 0.0;

	// diag scaling preprocess
	for (i = 0; i < ndof; i++)
		temp_r[i] = right_vec[i] * M_diag[i];

	// CG solver (diag scaling)
	for (i = 0; i < ndof; i++)
		solution_vec[i] = 0.0;

	M_mat_vec_crs(M, M_Ptr, M_Col, dd, solution_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		gg[i] = temp_r[i] - dd[i];
		pp[i] = gg[i];
	}

	for (itr = 0; itr < max_itr; itr++)
	{
		double ppp = inner_product(ndof, gg, gg);
		M_mat_vec_crs(M, M_Ptr, M_Col, dd, pp, ndof);
		double rrr = inner_product(ndof, dd, pp);
		double alpha = ppp / rrr;
		for (j = 0; j < ndof; j++)
		{
			solution_vec[j] += alpha * pp[j];
			gg[j] -= alpha * dd[j];
		}
		double qqq = inner_product(ndof, gg, dd);
		double beta = qqq / rrr;
		for (j = 0; j < ndof; j++)
			pp[j] = gg[j] - beta * pp[j];

		// M_check_conv_CG
		double rrr1 = inner_product(ndof, pp, pp), rrr2 = inner_product(ndof, solution_vec, solution_vec);
		rrr3 = fabs(alpha) * sqrt(rrr1 / rrr2);
		if (rrr3 < eps)
			break;
	}

	// diag scaling postprocess
	for (i = 0; i < ndof; i++)
		solution_vec[i] = solution_vec[i] * M_diag[i];

	printf("\t\tCG itr %d\teps %le\n", itr, rrr3);
}


void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int &ndof)
{
	int i, j, icount = 0;

	for (i = 0; i < ndof; i++)
		vec_result[i] = 0;
	for (i = 0; i < ndof; i++)
	{
		for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
		{
			vec_result[i] += M[icount] * vec[M_Col[j]];
			if (i != M_Col[j])
				vec_result[M_Col[j]] += M[icount] * vec[i];
			icount++;
		}
	}
}


double inner_product(int ndof, double *vec1, double *vec2)
{
	double result = 0.0;
	for (int i = 0; i < ndof; i++)
		result += vec1[i] * vec2[i];

	return result;
}


int RowCol_to_icount(int row, int col, information *info)
{
	// Guard invalid row/col and empty row range.
	if (row < 0 || row >= K_Whole_Size)
		return -1;
	if (col < row)
		return -1;

	long long start = info->K_Whole_Ptr[row];
	long long end = info->K_Whole_Ptr[row + 1];
	if (start < 0 || end < start)
		return -1;
	if (start == end)
		return -1;

	// Quick reject when col is out of this row's sorted column range.
	if (col < info->K_Whole_Col[start] || col > info->K_Whole_Col[end - 1])
		return -1;

	// 二分探索
	long long min_id = start;
	long long max_id = end - 1;
	while (min_id <= max_id)
	{
		long long mid_id = min_id + (max_id - min_id) / 2;

		if (info->K_Whole_Col[mid_id] == col)
			return mid_id;
		else if (info->K_Whole_Col[mid_id] < col)
			min_id = mid_id + 1;
		else
			max_id = mid_id - 1;
	}
	return -1;
}


// GMRES solver
void GMRES(int length, double eps, information *info)
{
	int i, j, k, l, m, n;
	constexpr int max_itr = 70;
	constexpr int itr_inf = 10;
	const double tol = eps;
	long long size = max_itr + 1;
	int itr;
	double norm;
	const int progress_total = itr_inf * max_itr;
	const int progress_interval = max(1, progress_total / 100);
	int progress_current = 0;
	bool progress_completed = false;
	int M_mode = 2; // 0: M = K, 1: M = 各パッチのKの対角ブロック, 2: M = [[K^G, 0], [0, K^L]]

	double *H = (double *)malloc(sizeof(double) * size * size);
	double *V = (double *)malloc(sizeof(double) * size * length);
	double *Omega = (double *)malloc(sizeof(double) * size * size);
	double *QH = (double *)malloc(sizeof(double) * size * size);
	double *QH_temp = (double *)malloc(sizeof(double) * size * size);
	double *P = (double *)malloc(sizeof(double) * size * length);

	double *e1 = (double *)malloc(sizeof(double) * size);
	double *g = (double *)malloc(sizeof(double) * size);
	double *g_temp = (double *)malloc(sizeof(double) * size);
	double *y = (double *)malloc(sizeof(double) * size);
	double *Ax = (double *)malloc(sizeof(double) * length);
	double *v = (double *)malloc(sizeof(double) * length);
	double *r = (double *)malloc(sizeof(double) * length);
	double *v_hat = (double *)malloc(sizeof(double) * length);
	double *Av = (double *)malloc(sizeof(double) * length);
	double *h = (double *)malloc(sizeof(double) * size);
	double *p = (double *)malloc(sizeof(double) * length);

	double *gg = (double *)malloc(sizeof(double) * length);
	double *dd = (double *)malloc(sizeof(double) * length);
	double *pp = (double *)malloc(sizeof(double) * length);
	double *temp_r = (double *)malloc(sizeof(double) * length);

	double *M = (double *)malloc(sizeof(double) * info->K_Whole_Ptr[length]);
	int *M_Ptr = (int *)malloc(sizeof(int) * (length + 1));
	int *M_Col = (int *)malloc(sizeof(int) * info->K_Whole_Ptr[length]);
	double *M_diag = (double *)malloc(sizeof(double) * length);

	// 前処理行列作成
	Make_M(M_mode, M, M_Ptr, M_Col, M_diag, length, info);

	// x0
	for (i = 0; i < length; i++)
		info->sol_vec[i] = 0.0;

	fprintf(stdout, "Solving with GMRES...\n");
	fflush(stdout);

	// loop itr_inf, set x0
	for (n = 0; n < itr_inf; n++)
	{
		// Ax0
		for (i = 0; i < length; i++)
		{
			double dot = 0.0;
			for (j = 0; j < length; j++)
			{
				int temp = 0;
				if (i <= j)
					temp = RowCol_to_icount(i, j, info);
				else if (i > j)
					temp = RowCol_to_icount(j, i, info);

				if (temp != -1)
					dot += info->K_Whole_Val[temp] * info->sol_vec[j];
			}
			Ax[i] = dot;
		}

		// r initialization
		for (i = 0; i < length; i++)
			r[i] = info->rhs_vec[i] - Ax[i];

		// l2 norm
		norm = 0.0;
		for (i = 0; i < length; i++)
			norm += pow(r[i], 2);
		norm = sqrt(norm);

		// e1
		e1[0] = norm;
		for (i = 1; i < size; i++)
			e1[i] = 0.0;

		// v
		for (i = 0; i < length; i++)
			v[i] = r[i] / norm;

		// loop
		for (i = 0; i < max_itr; i++)
		{
			CG(length, p, M, M_Ptr, M_Col, M_diag, v, gg, dd, pp, temp_r);

			// Av
			for (j = 0; j < length; j++)
			{
				double dot = 0.0;
				for (k = 0; k < length; k++)
				{
					int temp = 0;
					if (j <= k)
						temp = RowCol_to_icount(j, k, info);
					else if (j > k)
						temp = RowCol_to_icount(k, j, info);

					if (temp != -1)
						dot += info->K_Whole_Val[temp] * p[k];
				}
				Av[j] = dot;
			}

			// P
			for (j = 0; j < length; j++)
				P[i * length + j] = p[j];

			// V
			for (j = 0; j < length; j++)
				V[i * length + j] = v[j];

			// h
			for (j = 0; j <= i; j++)
			{
				double dot = 0.0;
				for (k = 0; k < length; k++)
					dot += Av[k] * V[j * length + k];
				h[j] = dot;
			}

			// v_hat
			for (j = 0; j < length; j++)
			{
				v_hat[j] = Av[j];
				for (k = 0; k <= i; k++)
					v_hat[j] -= h[k] * V[k * length + j];
			}

			// v_hat l2 norm
			norm = 0.0;
			for (j = 0; j < length; j++)
				norm += pow(v_hat[j], 2);
			norm = sqrt(norm);

			// h
			h[i + 1] = norm;

			// H
			for (j = 0; j <= i + 1; j++)
				H[j * size + i] = h[j];

			// v
			for (j = 0; j < length; j++)
				v[j] = v_hat[j] / norm;

			// givens rotation
			for (j = 0; j <= i; j++)
			{
				double nu, c_i, s_i;
				double H_a = 0, H_b = 0;

				if (j == 0)
					H_a = H[j * size + j], H_b = H[(j + 1) * size + j];
				else
					H_a = QH[j * size + j], H_b = QH[(j + 1) * size + j];

				nu = sqrt(pow(H_a, 2) + pow(H_b, 2));
				c_i = H_a / nu;
				s_i = H_b / nu;

				// Omega
				for (k = 0; k <= i + 1; k++)
					for (l = 0; l <= i + 1; l++)
					{
						if (k == l)
							Omega[k * size + l] = 1.0;
						else
							Omega[k * size + l] = 0.0;
					}

				Omega[j * size + j] = c_i;
				Omega[j * size + j + 1] = s_i;
				Omega[(j + 1) * size + j] = -s_i;
				Omega[(j + 1) * size + j + 1] = c_i;

				// QH
				if (j == 0)
					for (k = 0; k <= i + 1; k++)
						for (l = 0; l <= i; l++)
							QH[k * size + l] = H[k * size + l];
				for (k = 0; k <= i + 1; k++)
					for (l = 0; l <= i; l++)
						QH_temp[k * size + l] = 0.0;
				for (k = 0; k <= i + 1; k++)
					for (m = 0; m <= i + 1; m++)
						for (l = 0; l <= i; l++)
							QH_temp[k * size + l] += Omega[k * size + m] * QH[m * size + l];
				for (k = 0; k <= i + 1; k++)
					for (l = 0; l <= i; l++)
						QH[k * size + l] = QH_temp[k * size + l];

				// g
				if (j == 0)
				{
					for (k = 0; k <= i + 1; k++)
						g[k] = 0.0;
					for (k = 0; k <= i + 1; k++)
						for (l = 0; l <= i + 1; l++)
							g[k] += Omega[k * size + l] * e1[l];
					for (k = 0; k <= i + 1; k++)
						g_temp[k] = g[k];
				}
				else
				{
					for (k = 0; k <= i + 1; k++)
						g_temp[k] = 0.0;
					for (k = 0; k <= i + 1; k++)
						for (l = 0; l <= i + 1; l++)
							g_temp[k] += Omega[k * size + l] * g[l];
					for (k = 0; k <= i + 1; k++)
						g[k] = g_temp[k];
				}
			}

			// r_i norm
			norm = fabs(g[i + 1]);

			progress_current = n * max_itr + (i + 1);
			if (progress_current % progress_interval == 0 || norm <= tol || progress_current == progress_total)
			{
				Print_Solver_Progress("GMRES", progress_current, progress_total);
				if (progress_current == progress_total)
					progress_completed = true;
			}

			itr = i;
			if (norm <= tol)
				break;
		}

		if (norm <= tol)
			break;

		// y
		for (j = itr; j >= 0; j--)
		{
			y[j] = g_temp[j] / QH[j * size + j];
			for (k = 0; k <= j; k++)
				g_temp[k] -= QH[k * size + j] * y[j];
		}

		// update x0
		if (n != itr_inf - 1)
			for (i = 0; i < length; i++)
				for (j = 0; j <= itr; j++)
					info->sol_vec[i] += P[j * length + i] * y[j];
	}

	// y
	for (j = itr; j >= 0; j--)
	{
		y[j] = g_temp[j] / QH[j * size + j];
		for (k = 0; k <= j; k++)
			g_temp[k] -= QH[k * size + j] * y[j];
	}

	// sol_vec
	for (i = 0; i < length; i++)
		for (j = 0; j <= itr; j++)
			info->sol_vec[i] += P[j * length + i] * y[j];

	if (!progress_completed)
		Print_Solver_Progress("GMRES", progress_total, progress_total);
}


// Gaussian Elimination
void GaussianElimination(double *sol, double *r, double *A, int size)
{
	cout << "GaussianElimination is performed." << endl;
	double *r_for_swap = (double *)malloc(sizeof(double) * size);

	// 前進消去
	for (int i = 0; i < size; i++)
	{
		// 同列で絶対値が最大の行とピボット
		int swap_row = i;
		double max_val = 0.0;
		for (int j = i; j < size; j++)
		{
			if (double temp = fabs(A[j * size + i]); max_val < temp)
			{
				max_val = temp;
				swap_row = j;
			}
		}
		swap_matrix(A, r, i, swap_row, size, r_for_swap);

		// 行の対角項の値を1にする
		double a = A[i * size + i];
		A[i * size + i] = 1.0;
		r[i] /= a;
		if (i == size - 1)
			break;
		for (int j = i + 1; j < size; j++)
			A[i * size + j] /= a;

		// 消去
		for (int j = i + 1; j < size; j++)
		{
			double coefficient = A[j * size + i];
			r[j] -= coefficient * r[i];
			for (int k = i; k < size; k++)
				A[j * size + k] -= coefficient * A[i * size + k];
		}
	}

	// 後退代入
	for (int j = size - 1; j >= 0; j--)
	{
		double coefficient = 0.0;
		for (int i = j; i >= 0; i--)
		{
			if (i == j)
			{
				coefficient = r[i];
				sol[i] = coefficient;
			}
			else
				r[i] -= coefficient * A[i * size + j];
		}
	}

	cout << "Finish GaussianElimination" << endl;
}


bool GaussianElimination2(double *sol, const double *diff, double *J, int dimension)
{
	const double EPSILON = std::numeric_limits<double>::epsilon();

	// initialize solution
	for (int i = 0; i < dimension; ++i)
		sol[i] = diff[i];

	// Forward elimination
	for (int i = 0; i < dimension; ++i)
	{
		// Find the pivot
		int pivot = i;
		for (int j = i + 1; j < dimension; ++j)
		{
			if (std::abs(J[j * dimension + i]) > std::abs(J[pivot * dimension + i]))
			{
				pivot = j;
			}
		}

		// Swap rows if necessary
		if (pivot != i)
		{
			for (int k = 0; k < dimension; ++k)
			{
				std::swap(J[i * dimension + k], J[pivot * dimension + k]);
			}
			std::swap(sol[i], sol[pivot]);
		}

		// Check for singular matrix
		if (std::abs(J[i * dimension + i]) < EPSILON)
		{
			return false; // Matrix is singular
		}

		// Eliminate column below pivot
		for (int j = i + 1; j < dimension; ++j)
		{
			double factor = J[j * dimension + i] / J[i * dimension + i];
			for (int k = i; k < dimension; ++k)
			{
				J[j * dimension + k] -= factor * J[i * dimension + k];
			}
			sol[j] -= factor * sol[i];
		}
	}

	// Back substitution
	for (int i = dimension - 1; i >= 0; --i)
	{
		sol[i] /= J[i * dimension + i];
		for (int j = 0; j < i; ++j)
		{
			sol[j] -= J[j * dimension + i] * sol[i];
		}
	}

	return true; // Success
}


void swap_matrix(double *A, double *r, int row1, int row2, int size, double *r_for_swap)
{
	if (row1 == row2)
		return;

	int i;
	double temp;
	for (i = 0; i < size; i++)
		r_for_swap[i] = A[row1 * size + i];
	for (i = 0; i < size; i++)
	{
		temp = A[row2 * size + i];
		A[row2 * size + i] = r_for_swap[i];
		r_for_swap[i] = temp;
	}
	for (i = 0; i < size; i++)
		A[row1 * size + i] = r_for_swap[i];

	temp = r[row2];
	r[row2] = r[row1];
	r[row1] = temp;
}


#if 0
// LU decomposition
void LU(double *sol, double *r, double *A, int size)
{
	// partialPivLu
	MatrixXd U(size, size);
	VectorXd x(size);
	VectorXd b(size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			U(i, j) = A[(long)i * size + j];
		}
		b(i) = r[i];
	}

	x = U.partialPivLu().solve(b);
	// x = U.fullPivLu().solve(b);

	for (int i = 0; i < size; i++)
	{
		sol[i] = x(i);
	}
}


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void sparseLU(double *sol, double *r, int size, information *info)
{
    std::vector<T> coefficients;
    coefficients.reserve(size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // Assuming A is stored in row-major format
			if (i <= j && info->K_Whole_Val[RowCol_to_icount(i, j, info)] != -1)
                coefficients.push_back(T(i, j, info->K_Whole_Val[RowCol_to_icount(i, j, info)]));
			else if (i > j && info->K_Whole_Val[RowCol_to_icount(j, i, info)] != -1)
                coefficients.push_back(T(i, j, info->K_Whole_Val[RowCol_to_icount(j, i, info)]));
        }
    }

    SpMat A_sparse(size, size);
    A_sparse.setFromTriplets(coefficients.begin(), coefficients.end());

    // Perform LU decomposition
    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(A_sparse);
    solver.factorize(A_sparse);

    if (solver.info() != Eigen::Success)
    {
        // Handle failure
        // (e.g., matrix is singular or not positive definite)
        return;
    }

    Eigen::VectorXd b(size);
    for (int i = 0; i < size; i++)
    {
        b(i) = r[i];
    }

    Eigen::VectorXd x = solver.solve(b);

    for (int i = 0; i < size; i++)
    {
        sol[i] = x(i);
    }
}
#endif


#ifndef MKL_ILP64
#define MKL_ILP64
#endif
#include <mkl.h>
void intel_PARDISO(double *sol, double *rhs, int size, information *info)
{
	MKL_INT n = size;
	long long nnz = info->K_Whole_Ptr[n];

	// int *ia = (int *)malloc(sizeof(int) * (n + 1));
	// int *ja = (int *)malloc(sizeof(int) * nnz);
	// double *a = (double *)malloc(sizeof(double) * nnz);
	vector<MKL_INT> ia(n + 1);
	vector<MKL_INT> ja((MKL_INT)nnz);
	vector<double> a((MKL_INT)nnz);
	double *r = rhs;

	#pragma omp parallel for
	for (MKL_INT i = 0; i < n; i++)
	{
		ia[i] = (MKL_INT)info->K_Whole_Ptr[i];
		for (long long j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
		{
			ja[(MKL_INT)j] = info->K_Whole_Col[j];
			a[(MKL_INT)j] = info->K_Whole_Val[j];
		}
	}
	ia[n] = (MKL_INT)nnz;

	// int *perm = (int *)malloc(sizeof(int) * n);

	// Real symmetric matrix
	MKL_INT mtype = -2;

	// Number of right hand sides
	MKL_INT nrhs = 1;

	// Internal solver memory pointer
	void *pt[64];

	// Pardiso control parameters
	MKL_INT iparm[64] = {0};
	MKL_INT maxfct, mnum, phase, error, msglvl;

	// dummy
	MKL_INT idum;
	double ddum;

	// Setup Pardiso control parameters
	error = 0;

	// iparm[0] = 0;
	// iparm[1] = 2;
	// iparm[2] = omp_get_max_threads();
	// iparm[3] = 1;
	// iparm[4] = 0;

	iparm[0] = 1;	/* No solver default */
	iparm[1] = 2;	/* Fill-in reordering from METIS */
	iparm[3] = 0;	/* No iterative-direct algorithm */
	iparm[4] = 0;	/* No user fill-in reducing permutation */
	iparm[5] = 0;	/* Write solution into x */
	iparm[7] = 2;	/* Max numbers of iterative refinement steps */
	iparm[9] = 13;	/* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;	/* Use nonsymmetric permutation and scaling MPS */
	iparm[12] = 0;	/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;	/* Output: Number of perturbed pivots */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0;	/* Output: Numbers of CG Iterations */
	iparm[34] = 1;	/* PARDISO use C-style indexing for ia and ja arrays */

	// Maximum number of numerical factorizations
	maxfct = 1;

	// Which factorization to use
	mnum = 1;

	// Print statistical information
	// msglvl = 0;
	msglvl = 1;

	// Initialize error flag
	error = 0;

	/* ----------------------------------------------------------------*/
	/* .. Initialize the internal solver memory pointer. This is only  */
	/*   necessary for the FIRST call of the PARDISO solver.           */
	/* ----------------------------------------------------------------*/
	for (int i = 0; i < 64; i++)
		pt[i] = 0;

	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase = 11;

	printf("PARDISO: attempting symbolic factorization (nnz=%lld)...\n", nnz);
	fflush(stdout);

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %lld\n", error);
		exit(2);
	}
	// printf("\nReordering completed ... ");
	// printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	// printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase = 22;

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %lld\n", error);
		exit(2);
	}
	// printf("\nFactorization completed ... ");

	// printf("perturbed pivots     = %lld\n", (long long)iparm[13]);
	// printf("positive eigenvalues = %lld\n", (long long)iparm[21]);
	// printf("negative eigenvalues = %lld\n", (long long)iparm[22]);
	// printf("zero eigenvalues     = %lld\n", (long long)(n - iparm[21] - iparm[22]));
	// fflush(stdout);

	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm, &msglvl, r, sol, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %lld", error);
		exit(3);
	}

	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia.data(), ja.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	// phase = 0;
	// PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, r, sol, &error);
}


#include <eigen3/Eigen/Dense>
using namespace Eigen;
void temp_solver(double *sol, double *r, double *M, int row_n, int col_n)
{
	// partialPivLu
	MatrixXd U(row_n, col_n);
	VectorXd x(col_n);
	VectorXd b(row_n);
	for (int i = 0; i < row_n; i++)
	{
		for (int j = 0; j < col_n; j++)
		{
			U(i, j) = M[(long)i * col_n + j];
		}
		b(i) = r[i];
	}

	x = U.colPivHouseholderQr().solve(b);
	// x = U.fullPivHouseholderQr().solve(b);
	// x = U.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	// x = U.partialPivLu().solve(b);
	// x = U.fullPivLu().solve(b);

	for (int i = 0; i < col_n; i++)
	{
		sol[i] = x(i);
	}
}


void temp_solver_v2(double *sol, double *r, double *A, int row_n, int col_n)
{
	// 疑似逆行列の計算
	int i, j, k;

	// double *A = (double *)malloc(sizeof(double) * (row_n) * col_n);
	double *A_T = (double *)malloc(sizeof(double) * col_n * (row_n));
	double *B = (double *)malloc(sizeof(double) * col_n * col_n);
	double *B_inv = (double *)malloc(sizeof(double) * col_n * col_n);
	double *A_pinv = (double *)malloc(sizeof(double) * (row_n) * col_n);
	double *A_pinv_T = (double *)malloc(sizeof(double) * col_n * (row_n));

	// for (i = 0; i < row_n; i++)
	// {
	// 	for (j = 0; j < col_n; j++)
	// 	{
	// 		A[i * col_n + j] = T[n * (row_n) * col_n + i * col_n + j];
	// 	}
	// }

	for (i = 0; i < row_n; i++)
	{
		for (j = 0; j < col_n; j++)
		{
			A_T[j * (row_n) + i] = A[i * col_n + j];
		}
	}

	for (i = 0; i < col_n; i++)
	{
		for (j = 0; j < col_n; j++)
		{
			B[i * col_n + j] = 0.0;
			for (k = 0; k < row_n; k++)
			{
				B[i * col_n + j] += A_T[i * (row_n) + k] * A[k * col_n + j];
			}
		}
	}

	double temp;

	for (i = 0; i < col_n; i++)
	{
		for (j = 0; j < col_n; j++)
		{
			B_inv[i * col_n + j] = 0.0;
		}
	}

	for (i = 0; i < col_n; i++)
	{
		B_inv[i * col_n + i] = 1.0;
	}

	for (k = 0; k < col_n; k++)
	{
		temp = B[k * col_n + k];
		for (i = 0; i < col_n; i++)
		{
			B[k * col_n + i] /= temp;
			B_inv[k * col_n + i] /= temp;
		}

		for (i = 0; i < col_n; i++)
		{
			if (i != k)
			{
				temp = B[i * col_n + k];
				for (j = 0; j < col_n; j++)
				{
					B[i * col_n + j] -= B[k * col_n + j] * temp;
					B_inv[i * col_n + j] -= B_inv[k * col_n + j] * temp;
				}
			}
		}
	}

	for (i = 0; i < row_n; i++)
	{
		for (j = 0; j < col_n; j++)
		{
			A_pinv[i * col_n + j] = 0.0;
			for (k = 0; k < col_n; k++)
			{
				A_pinv[i * col_n + j] += A[i * col_n + k] * B_inv[k * col_n + j];
			}
		}
	}

	for (i = 0; i < row_n; i++)
	{
		for (j = 0; j < col_n; j++)
		{
			A_pinv_T[j * (row_n) + i] = A_pinv[i * col_n + j];
		}
	}

	// inner product
	for (i = 0; i < col_n; i++)
	{
		sol[i] = 0.0;
	}
	for (i = 0; i < col_n; i++)
	{
		for (j = 0; j < row_n; j++)
		{
			sol[i] += A_pinv_T[i * (row_n) + j] * r[j];
		}
	}
}


void overdetermined_system(double *sol,double *rhs, double *A, int m, int n)
{
    cout << "m = " << m << " n = " << n << endl;

    // double Z[] = {1, 2, 3, 4, 5, 6};
    // double Y[] = {6, 7, 9};
    
    // Eigen 行列に既存の行列をコピーする
    Eigen::Map<Eigen::MatrixXd> A_transpose(A, n, m);

    // Eigen ベクトルに既存のベクトルBをコピーする
    Eigen::Map<Eigen::VectorXd> B_eigen(rhs, m);

    // Aの転置行列
    Eigen::MatrixXd A_eigen = A_transpose.transpose();

    // cout << "A_transpose" << endl;
    // cout << A_transpose << endl;
    // cout << endl;

    // cout << "A_eigen" << endl;
    // cout << A_eigen << endl;
    // cout << endl;

    // A_transpose * A の逆行列を計算
    Eigen::MatrixXd A_transpose_A_inv = (A_transpose * A_eigen).inverse();
    
    // cout << "A_inverse" << endl;
    // cout << A_transpose_A_inv << endl;
    // cout << endl;

    // 解を計算
    Eigen::VectorXd x = A_transpose_A_inv * A_transpose * B_eigen;

    // 結果を出力
    // cout << "Solution x:\n" << x << "\n";

	copy(x.data(), x.data() + x.size(), sol);

    return;
}

#include <cfloat>
#include <chrono>
static const int MAX_DENSE_NDOF = 5000;

void Condition_Number_K(information *info)
{
	const int n = K_Whole_Size;
	// 密行列化で必要なメモリ (GB)
	const double mem_gb = 8.0 * (double)n * (double)n / 1e9;
 
	printf("\n============ Condition_Number_K ============\n");
	printf("ndof = %d, dense matrix memory ~ %.2f GB\n", n, mem_gb);
 
	auto t0 = std::chrono::steady_clock::now();
 
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n; i++)
		for (int j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
		{
			const int c = info->K_Whole_Col[j];
			K(i, c) = info->K_Whole_Val[j];
			K(c, i) = info->K_Whole_Val[j];
		}
 
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K, Eigen::EigenvaluesOnly);
	if (es.info() != Eigen::Success)
	{
		printf("固有値計算に失敗しました (Eigen::SelfAdjointEigenSolver)\n");
		return;
	}
 
	const Eigen::VectorXd &ev = es.eigenvalues(); // 昇順
	const double lam_min = ev(0);
	const double lam_max = ev(n - 1);
 
	double abs_min = DBL_MAX, abs_max = 0.0;
	int n_neg = 0;
	for (int i = 0; i < n; i++)
	{
		abs_min = std::min(abs_min, std::fabs(ev(i)));
		abs_max = std::max(abs_max, std::fabs(ev(i)));
		if (ev(i) < 0.0)
			n_neg++;
	}
 
	auto t1 = std::chrono::steady_clock::now();
 
	printf("lambda_min (signed) = % .10e\n", lam_min);
	printf("lambda_max (signed) = % .10e\n", lam_max);
	printf("# negative          = %d\n", n_neg);
	printf("|lambda|_min        = %.10e\n", abs_min);
	printf("|lambda|_max        = %.10e\n", abs_max);
	printf("cond(K) = |lambda|_max / |lambda|_min = %.10e\n", abs_max / abs_min);
	printf("elapsed = %.3f sec\n", std::chrono::duration<double>(t1 - t0).count());
	printf("==============================================\n");
}

// 昇順に並んだ固有値 ev から統計を表示する
static void Eigen_Report(const char *label, const Eigen::VectorXd &ev)
{
	const int n = (int)ev.size();
	if (n == 0)
		return;
 
	const double lam_min = ev(0);
	const double lam_max = ev(n - 1);
	const double abs_max = std::max(std::fabs(lam_min), std::fabs(lam_max));
 
	// 数値ランク判定の閾値 (LAPACK/NumPy と同じ考え方: n * eps * |λ|max)
	const double tol = abs_max * n * DBL_EPSILON;
 
	int n_neg = 0, n_zero = 0;
	double abs_min = abs_max;
	for (int i = 0; i < n; i++)
	{
		if (ev(i) < -tol)
			n_neg++;
		else if (ev(i) <= tol)
			n_zero++;
		abs_min = std::min(abs_min, std::fabs(ev(i)));
	}
 
	printf("\n---- %s (n = %d) ----\n", label, n);
	printf("  lambda_min      = % .6e\n", lam_min);
	printf("  lambda_max      = % .6e\n", lam_max);
	printf("  zero threshold  = %.3e  (n*eps*|lambda|max)\n", tol);
	printf("  # negative      = %d\n", n_neg);
	printf("  # zero (|l|<=tol) = %d\n", n_zero);
 
	if (n_zero > 0)
		printf("  cond (2-norm)   = inf  -> 特異 (rank deficient, 零空間の次元 = %d)\n", n_zero);
	else
		printf("  cond (2-norm)   = %.6e  (= |lambda|max / |lambda|min)\n", abs_max / abs_min);

	if (n_zero == 0)
		printf("  -> 正定値 (SPD)\n");
 
	// 小さい方から最大5個, 大きい方から3個
	printf("  smallest eigenvalues:");
	for (int i = 0; i < std::min(n, 5); i++)
		printf(" % .3e", ev(i));
	printf("\n  largest  eigenvalues:");
	for (int i = std::max(0, n - 3); i < n; i++)
		printf(" % .3e", ev(i));
	printf("\n");
 
	// cond が 1/eps 級のときは固有値の絶対誤差 (≒ eps*|λ|max) に埋もれる
	if (n_zero == 0 && abs_max / abs_min > 1.0 / DBL_EPSILON)
		printf("  !! cond > 1/eps: 最小固有値は倍精度の丸め誤差以下で信頼できません\n");
}

static void Eigen_Solve_And_Report(const char *label, const Eigen::MatrixXd &A)
{
	if (A.rows() == 0)
		return;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Eigen::EigenvaluesOnly);
	if (es.info() != Eigen::Success)
	{
		printf("\n---- %s ---- 固有値計算に失敗しました\n", label);
		return;
	}
	Eigen_Report(label, es.eigenvalues()); // 昇順
}
 
void Check_K_Eigen(information *info)
{
	const int n = K_Whole_Size;
 
	printf("\n================ Check_K_Eigen ================\n");
	if (n > MAX_DENSE_NDOF)
	{
		printf("ndof = %d > MAX_DENSE_NDOF (%d) のためスキップします\n", n, MAX_DENSE_NDOF);
		return;
	}
 
	auto t0 = std::chrono::steady_clock::now();
 
	// グローバルパッチの自由度数 (Make_M と同じ求め方)
	int ng = 0;
	for (int i = 0; i < info->Total_Control_Point_on_mesh[0] * info->DIMENSION; i++)
		if (info->Index_Dof[i] != ERROR)
			ng++;
	const int nl = n - ng;
 
	// 上三角 CRS -> 対称な密行列
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n; i++)
		for (int j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
		{
			const int c = info->K_Whole_Col[j];
			K(i, c) = info->K_Whole_Val[j];
			K(c, i) = info->K_Whole_Val[j];
		}
 
	// (1) K そのもの
	Eigen_Solve_And_Report("K", K);
 
	// (2) 対角スケーリング後 D^-1/2 K D^-1/2 (Make_M の M_diag と同じスケーリング)
	{
		Eigen::VectorXd d(n);
		bool ok = true;
		for (int i = 0; i < n; i++)
		{
			if (K(i, i) <= 0.0)
			{
				printf("\n!! K(%d,%d) = %.3e <= 0 : 対角が非正のためスケーリング不可 (Make_M の sqrt が NaN になります)\n", i, i, K(i, i));
				ok = false;
				break;
			}
			d(i) = 1.0 / std::sqrt(K(i, i));
		}
		if (ok)
		{
			Eigen::MatrixXd Ks = d.asDiagonal() * K * d.asDiagonal();
			Eigen_Solve_And_Report("D^-1/2 K D^-1/2 (diag scaled)", Ks);
		}
	}
 
	// (3) 対角ブロック K^G, K^L
	if (ng > 0 && nl > 0)
	{
		Eigen_Solve_And_Report("K^G (global block)", K.topLeftCorner(ng, ng));
		Eigen_Solve_And_Report("K^L (local block)", K.bottomRightCorner(nl, nl));
 
		// (4) 前処理後 M^-1 K の固有値 (M = blockdiag(K^G, K^L))
		Eigen::MatrixXd M = K;
		M.topRightCorner(ng, nl).setZero();
		M.bottomLeftCorner(nl, ng).setZero();
 
		Eigen::LLT<Eigen::MatrixXd> llt(M);
		if (llt.info() == Eigen::Success)
		{
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K, M, Eigen::EigenvaluesOnly);
			if (ges.info() == Eigen::Success)
				Eigen_Report("M^-1 K (preconditioned, K x = lambda M x)", ges.eigenvalues());
			else
				printf("\n---- M^-1 K ---- 一般化固有値計算に失敗しました\n");
		}
		else
		{
			printf("\n---- M^-1 K ---- M = blockdiag(K^G, K^L) が正定値でないため (Cholesky 失敗) スキップ\n");
		}
	}
 
	auto t1 = std::chrono::steady_clock::now();
	printf("\nCheck_K_Eigen: %.3f sec\n", std::chrono::duration<double>(t1 - t0).count());
	printf("===============================================\n");
}