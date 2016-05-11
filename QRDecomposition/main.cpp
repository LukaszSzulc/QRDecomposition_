#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <float.h>
//#define float double

struct Matrix
{
	int n;
	int m;
	float **matrix;
};
struct Vector
{
	float n;
	float *vector;
};

struct QResult
{
	Matrix aPrim;
	Matrix qMatrix;
};

struct RResult
{
	Matrix rMatrix;
};

void PrintMatrix(Matrix matrix)
{
	for (int i = 0; i < matrix.n; i++)
	{
		for (int j = 0; j < matrix.m; j++)
		{
			printf("%f ", matrix.matrix[i][j]);
		}

		printf("\n");
	}
}

void PrintVector(Vector v)
{
	for (int i = 0; i < v.n; i++)
	{
		printf("%f ", v.vector[i]);
	}

	printf("\n");
}

void FreeVector(Vector v)
{
	free(v.vector);
}

Vector InitializeVector(int n)
{
	Vector result;
	result.n = n;
	result.vector = (float*)malloc(n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		result.vector[i] = 0.0f;
	}

	return result;
}

Matrix Initialize(int n, int m)
{
	Matrix result;
	result.m = m;
	result.n = n;
	result.matrix = (float**)malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++)
	{
		result.matrix[i] = (float*)malloc(m * sizeof(float));
		for (int j = 0; j < m; j++)
		{
			result.matrix[i][j] = 0.0f;
		}
	}

	return result;
}

Matrix InitializeMatrixFromFile(char *path)
{
	Matrix result;
	int row = 0;
	int column = 0;
	FILE *f;
	char line[100];
	f = fopen(path, "r");
	fgets(line, 1024, f);
	char * n_ = strtok(line, " ");
	int n = atoi(n_);
	result.n = n;
	result.m = n;
	result.matrix = (float**)malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++)
	{
		result.matrix[i] = (float*)malloc(n * sizeof(float));
	}

	while (fgets(line, 1024, f) != NULL)
	{
		char *token = strtok(line, " ");
		result.matrix[row][column++] = atof(token);
		while (token)
		{
			token = strtok(NULL, " ");
			if (token == NULL)
			{
				break;
			}

			result.matrix[row][column] = atof(token);
			column++;
		}
		row++;
		column = 0;
	}

	//PrintMatrix(result);
	return result;
}

Vector ExtractColumn(Matrix matrix, int columnToExtract)
{
	Vector v = InitializeVector(matrix.n);
	for (int i = 0; i < matrix.n; i++)
	{
		v.vector[i] = matrix.matrix[i][columnToExtract];
	}

	return v;
}

float CalculateNorm(Vector v)
{
	float norm = 0;
	for (int i = 0; i < v.n; i++)
	{
		norm += v.vector[i]*v.vector[i];
	}

	return sqrt(norm);
}

Vector NormalizeVector(Vector v, float norm)
{
	Vector result = InitializeVector(v.n);
	for (int i = 0; i < result.n; i++)
	{
		result.vector[i] = v.vector[i] / norm;
	}

	return result;
}

void SetVectorInMatrixColumn(Matrix m, Vector v, int column)
{
	for (int i = 0; i < m.n; i++)
	{
		m.matrix[i][column] = v.vector[i];
	}
}

float MultiplyVectorByVector(Vector v1, Vector v2)
{
	float sum = 0;
	for (int i = 0; i < v1.n; i++)
	{
		sum += v1.vector[i] * v2.vector[i];
	}

	return sum;
}

void MultiplyVectorByScalar(Vector v, float scalar)
{
	for (int i = 0; i < v.n; i++)
	{
		v.vector[i] *= scalar;
	}
}

void AddVectorToVector(Vector v1, Vector v2)
{
	for (int i = 0; i < v1.n; i++)
	{
		v1.vector[i] += v2.vector[i];
	}
}

void SubstractVector(Vector v1, Vector v2)
{
	for (int i = 0; i < v1.n; i++)
	{
		v1.vector[i] -= v2.vector[i];
	}
}

QResult CreateQMatrix(Matrix m)
{
	Matrix qMatrix = Initialize(m.n, m.m);
	Matrix aPrim = Initialize(m.n, m.n);
	Vector firstColumn = ExtractColumn(m, 0);
	float norm = CalculateNorm(firstColumn);
	Vector normalizedVector = NormalizeVector(firstColumn, norm);
	SetVectorInMatrixColumn(qMatrix, normalizedVector, 0);
	for (int i = 1; i < qMatrix.n; i++)
	{
		Vector fromAMatrix = ExtractColumn(m, i);
		Vector tmpResult = InitializeVector(qMatrix.n);
		for (int j = 0; j < i; j++)
		{
			Vector fromQMatrix = ExtractColumn(qMatrix, j);
			float multiplyResult = MultiplyVectorByVector(fromAMatrix, fromQMatrix);
			MultiplyVectorByScalar(fromQMatrix, multiplyResult);
			AddVectorToVector(tmpResult, fromQMatrix);
			FreeVector(fromQMatrix);
		}
		SubstractVector(fromAMatrix, tmpResult);
		float _norm = CalculateNorm(fromAMatrix);
		Vector normalized = NormalizeVector(fromAMatrix,_norm);
		SetVectorInMatrixColumn(aPrim, fromAMatrix, i);
		SetVectorInMatrixColumn(qMatrix, normalized, i);
		FreeVector(fromAMatrix);
		FreeVector(tmpResult);
		FreeVector(normalized);
	}

	QResult q;
	q.aPrim = aPrim;
	q.qMatrix = qMatrix;
	return q;
}

RResult CreateRMatrix(QResult qResult, Matrix orginal)
{
	RResult result;
	result.rMatrix = Initialize(qResult.aPrim.n, qResult.aPrim.m);
	for (int i = 0; i < qResult.qMatrix.n; i++)
	{
		for (int j = 0; j < qResult.qMatrix.m; j++)
		{
			if (i == j)
			{
				Vector v = ExtractColumn(qResult.aPrim, i);
				float norm = CalculateNorm(v);
				result.rMatrix.matrix[i][i] = norm;
				FreeVector(v);
			}
			else
			{
				if (i < j)
				{
					Vector a = ExtractColumn(qResult.qMatrix, i);
					Vector b = ExtractColumn(orginal, j);
					float multiplyResult = MultiplyVectorByVector(a, b);
					result.rMatrix.matrix[i][j] = multiplyResult;
					FreeVector(a);
					FreeVector(b);
				}
				else
				{
					result.rMatrix.matrix[i][j] = 0.0f;
				}
			}

		}
	}
	result.rMatrix.matrix[0][0] = CalculateNorm(ExtractColumn(orginal, 0));
	return result;
}

Matrix MultipyMatrixByMatrix(Matrix a, Matrix b)
{
	Matrix result = Initialize(a.n, b.n);
	for (int i = 0; i < a.n; i++)
	{
		for (int j = 0; j < b.m; j++)
		{
			float sum = 0.0f;
			for (int k = 0; k < a.m; k++)
			{
				sum += a.matrix[i][k] * b.matrix[k][j];
			}

			result.matrix[i][j] += sum;
		}
	}

	return result;
}

void FreeMatrix(Matrix m)
{
	for (int i = 0; i < m.n; i++)
	{
		free(m.matrix[i]);
	}

	free(m.matrix);
}

bool ChangesDetected(Vector first, Vector second)
{
	for (int i = 0; i < first.n; i++)
	{
		if (second.vector[i] - first.vector[i] > FLT_EPSILON)
		{
			return false;
		}
	}

	return true;
}

Vector GetDiagonal(Matrix m)
{

	Vector diag = InitializeVector(m.n);
	for (int i = 0; i < m.n; i++)
	{
		diag.vector[i] = m.matrix[i][i];
	}

	return diag;
}

Vector FindEigenValues(Matrix orginal)
{
	Matrix tmp = orginal;
	Vector diag1;
	Vector diag2;
	do
	{
		QResult q = CreateQMatrix(tmp);
		RResult r = CreateRMatrix(q, tmp);
		Matrix result = MultipyMatrixByMatrix(q.qMatrix, r.rMatrix);
		tmp = MultipyMatrixByMatrix(r.rMatrix, q.qMatrix);
		diag1 = GetDiagonal(result);
		diag2 = GetDiagonal(tmp);
		FreeMatrix(q.qMatrix);
		FreeMatrix(r.rMatrix);
	} while (!ChangesDetected(diag1, diag2));
	FreeVector(diag1);
	return diag2;
}

int main(int argc, char* argv[])
{
	char *path = argv[3];
	Matrix initial = InitializeMatrixFromFile(path);
	Vector eig = FindEigenValues(initial);
	PrintVector(eig);
	FreeMatrix(initial);
	FreeVector(eig);
	system("pause");
	return EXIT_SUCCESS;
}