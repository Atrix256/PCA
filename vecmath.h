#pragma once

#include <array>
#include <stdio.h>

// ==================================================== BASE TYPES ====================================================

template <typename T, size_t SIZE>
using BaseVec = std::array<T, SIZE>;

template <typename T, size_t WIDTH, size_t HEIGHT>
using BaseMtx = std::array<BaseVec<T, WIDTH>, HEIGHT>;

template <typename T, size_t SIZE>
using BaseMtxSq = std::array<BaseVec<T, SIZE>, SIZE>;

// Matrices are indexed as [y][x].  Matrix[i] is the i'th row.

// ==================================================== Scalar Math ====================================================

template <typename T>
T BinomialCoefficient(T n, T k)
{
	T numerator = (T)1.0f;
	T denominator = (T)1.0f;

	for (int i = 1; i <= k; ++i)
	{
		numerator *= T(n + 1 - i);
		denominator *= T(i);
	}

	return numerator / denominator;
}

template <typename T>
T Clamp(T value, T themin, T themax)
{
	if (value < themin)
		return themin;

	if (value > themax)
		return themax;

	return value;
}

inline float Lerp(float A, float B, float t)
{
	return (1.0f - t) * A + t * B;
}

// ==================================================== Vector Math ====================================================

template <typename T, size_t SIZE>
T Dot(const BaseVec<T, SIZE>& A, const BaseVec<T, SIZE>& B)
{
	T sum = (T)0.0f;
	for (size_t i = 0; i < SIZE; ++i)
		sum += A[i] * B[i];
	return sum;
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> operator-(const BaseVec<T, SIZE>& A, const BaseVec<T, SIZE>& B)
{
	BaseVec<T, SIZE> ret;
	for (size_t i = 0; i < SIZE; ++i)
		ret[i] = A[i] - B[i];
	return ret;
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> operator+(const BaseVec<T, SIZE>& A, const BaseVec<T, SIZE>& B)
{
	BaseVec<T, SIZE> ret;
	for (size_t i = 0; i < SIZE; ++i)
		ret[i] = A[i] + B[i];
	return ret;
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> operator*(const BaseVec<T, SIZE>& A, float B)
{
	BaseVec<T, SIZE> ret;
	for (size_t i = 0; i < SIZE; ++i)
		ret[i] = A[i] * B;
	return ret;
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> operator*(float A, const BaseVec<T, SIZE>& B)
{
	return B * A;
}

template <typename T, size_t SIZE>
T Length(const BaseVec<T, SIZE>& A)
{
	return (T)sqrt(Dot(A, A));
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> Normalize(const BaseVec<T, SIZE>& A)
{
	T len = Length(A);
	BaseVec<T, SIZE> ret;
	for (size_t i = 0; i < SIZE; ++i)
		ret[i] = A[i] / len;
	return ret;
}

template <typename T, size_t SIZE>
BaseVec<T, SIZE> ProjectVectorAOntoVectorB(const BaseVec<T, SIZE>& A, const BaseVec<T, SIZE>& B)
{
	return B * (Dot(A, B) / Dot(B, B));
}

// ==================================================== Matrix Math ====================================================

template <typename T, size_t WIDTH, size_t HEIGHT>
constexpr size_t Rows(const BaseMtx<T, WIDTH, HEIGHT>& A)
{
	return HEIGHT;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
constexpr size_t Columns(const BaseMtx<T, WIDTH, HEIGHT>& A)
{
	return WIDTH;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
BaseVec<T, HEIGHT> Column(const BaseMtx<T, WIDTH, HEIGHT>& A, size_t index)
{
	BaseVec<T, HEIGHT> ret;
	for (size_t i = 0; i < HEIGHT; ++i)
		ret[i] = A[i][index];
	return ret;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
void SetColumn(BaseMtx<T, WIDTH, HEIGHT>& A, size_t index, const BaseVec<T, HEIGHT>& data)
{
	for (size_t i = 0; i < HEIGHT; ++i)
		A[i][index] = data[i];
}

template <typename T, size_t WIDTH, size_t HEIGHT>
BaseVec<T, WIDTH> Row(const BaseMtx<T, WIDTH, HEIGHT>& A, size_t index)
{
	return A[index];
}

template <typename T, size_t WIDTH, size_t HEIGHT>
void SetRow(BaseMtx<T, WIDTH, HEIGHT>& A, size_t index, BaseVec<T, WIDTH>& data)
{
	A[index] = data;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
BaseMtx<T, HEIGHT, WIDTH> Transpose(const BaseMtx<T, WIDTH, HEIGHT>& mtx)
{
	BaseMtx<T, HEIGHT, WIDTH> ret;
	for (size_t i = 0; i < WIDTH; ++i)
		ret[i] = Column(mtx, i);
	return ret;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
BaseVec<T, HEIGHT> Multiply(const BaseMtx<T, WIDTH, HEIGHT>& mtx, const BaseVec<T, WIDTH>& vec)
{
	BaseVec<T, HEIGHT> ret{};

	for (size_t i = 0; i < HEIGHT; ++i)
		ret[i] = Dot(Row(mtx, i), vec);

	return ret;
}

template <typename T, size_t WIDTH, size_t HEIGHT>
BaseVec<T, WIDTH> Multiply(const BaseVec<T, HEIGHT>& vec, const BaseMtx<T, WIDTH, HEIGHT>& mtx)
{
	return Multiply(Transpose(mtx), vec);
}

template <typename T, size_t AWIDTH, size_t AHEIGHT, size_t BWIDTH>
BaseMtx<T, BWIDTH, AHEIGHT> Multiply(const BaseMtx<T, AWIDTH, AHEIGHT>& A, const BaseMtx<T, BWIDTH, AWIDTH>& B)
{
	BaseMtx<T, BWIDTH, AHEIGHT> ret{};

	for (int im = 0; im < BWIDTH; ++im)
		for (int in = 0; in < AHEIGHT; ++in)
			ret[in][im] = Dot(Row(A, in), Column(B, im));

	return ret;
}

// ==================================================== Complex Functions ====================================================

template <typename T, size_t WIDTH, size_t HEIGHT>
void GaussJordanElimination(BaseMtx<T, WIDTH, HEIGHT>& augmentedMtx)
{
	// augmentedMtx is assumed to be an NxN matrix possibly with some "extra bits" of size (M-N, N) on the right side.

	// make each column in the matrix have only a single row with a value in it, and have that value be 1.
	for (int column = 0; column < HEIGHT; ++column)
	{
		// find the row that has the maximum absolute value for this column
		int maxValueRowIndex = column;
		T maxValue = augmentedMtx[column][column];
		for (int row = column + 1; row < HEIGHT; ++row)
		{
			if (abs(augmentedMtx[row][column]) > abs(maxValue))
			{
				maxValue = augmentedMtx[row][column];
				maxValueRowIndex = row;
			}
		}

		// swap rows if we need to
		if (column != maxValueRowIndex)
		{
			for (size_t ix = 0; ix < WIDTH; ++ix)
				std::swap(augmentedMtx[column][ix], augmentedMtx[maxValueRowIndex][ix]);
		}

		// scale this row by the value
		{
			T scale = augmentedMtx[column][column];
			for (size_t ix = 0; ix < WIDTH; ++ix)
				augmentedMtx[column][ix] /= scale;
		}

		// make only this row have a value in it, by adding or subtracting multiples of it from the other rows
		for (size_t iy = 0; iy < HEIGHT; ++iy)
		{
			if (iy == column)
				continue;

			T scale = augmentedMtx[iy][column];
			if (scale == 0.0f)
				continue;

			for (size_t ix = 0; ix < WIDTH; ++ix)
				augmentedMtx[iy][ix] -= augmentedMtx[column][ix] * scale;
		}
	}
}

// Makes columns of matrix orthonormal
// Uses the modified Gram-Schmidt process which is more numerically stable
// https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
template <typename T, size_t WIDTH, size_t HEIGHT>
BaseMtx<T, WIDTH, HEIGHT> GramSchmidt(const BaseMtx<T, WIDTH, HEIGHT>& mtx)
{
	BaseMtx<T, WIDTH, HEIGHT> ret{};

	// The 0th column is unchanged
	SetColumn(ret, 0, Column(mtx, 0));

	// every other column is adjusted to be orthogonal to the others before it
	for (int column = 1; column < WIDTH; ++column)
	{
		BaseVec<T, HEIGHT> vec = Column(mtx, column);
		for (int i = 0; i < column; ++i)
			vec = vec - ProjectVectorAOntoVectorB(vec, Column(ret, i));
		SetColumn(ret, column, vec);
	}

	// Normalize the column
	for (int column = 0; column < WIDTH; ++column)
		SetColumn(ret, column, Normalize(Column(ret, column)));

	return ret;
}

// https://en.wikipedia.org/wiki/QR_decomposition
template <typename T, size_t SIZE>
void QRDecomposition(const BaseMtxSq<T, SIZE>& mtx, BaseMtxSq<T, SIZE>& Q, BaseMtxSq<T, SIZE>& R)
{
	Q = GramSchmidt(mtx);
	R = Multiply(Transpose(Q), mtx);
}

// Calculates Eigenvalues of a matrix
// Iteratively uses QR decomposition and remakes the matrix as R*Q until convergence.
// QR decomposition uses the Gram-Schmidt process (modified to the more stable version)
template <typename T, size_t SIZE>
BaseVec<T, SIZE> QRAlgorithm(const BaseMtxSq<T, SIZE>& mtx, int maxIterations, float maxError, const char* debugFileName = nullptr)
{
	BaseMtxSq<T, SIZE> m = mtx;
	float maxabserror = 0.0f;
	float avgerror = 0.0f;
	int counterror = 0;

	// Calculate the starting stats
	for (int y = 0; y < SIZE; ++y)
	{
		for (int x = 0; x < SIZE; ++x)
		{
			if (x < y)
			{
				maxabserror = std::max(maxabserror, std::abs(m[y][x]));
				counterror++;
				avgerror = Lerp(avgerror, std::abs(m[y][x]), 1.0f / float(counterror));
			}
		}
	}

	FILE* file = nullptr;
	if (debugFileName)
	{
		fopen_s(&file, debugFileName, "wt");

		if (file)
		{
			fprintf(file, "\"Index\",\"Max\",\"Average\"\n");
			fprintf(file, "\"0\",\"%f\",\"%f\"\n", maxabserror, avgerror);
		}
	}

	for (int i = 0; i < maxIterations; ++i)
	{
		if (maxabserror < maxError)
			break;

		BaseMtxSq<T, SIZE> Q, R;
		QRDecomposition(m, Q, R);
		m = Multiply(R, Q);
		
		// Calculate stats
		maxabserror = 0.0f;
		avgerror = 0.0f;
		counterror = 0;
		for (int y = 0; y < SIZE; ++y)
		{
			for (int x = 0; x < SIZE; ++x)
			{
				if (x < y)
				{
					maxabserror = std::max(maxabserror, std::abs(m[y][x]));
					counterror++;
					avgerror = Lerp(avgerror, std::abs(m[y][x]), 1.0f / float(counterror));
				}
			}
		}

		// Log convergence if file is open
		if (file)
			fprintf(file, "\"%i\",\"%f\",\"%f\"\n", i + 1, maxabserror, avgerror);
	}

	if (file)
		fclose(file);

	// get the eigenvector
	BaseVec<T, SIZE> ret;
	for (int i = 0; i < SIZE; ++i)
		ret[i] = m[i][i];
	return ret;
}

// Get the approximation of the eigenvector associated with a given approximate eigenvalue
// https://en.wikipedia.org/wiki/Inverse_iteration
template <typename T, size_t SIZE>
BaseVec<T, SIZE> InverseIteration(const BaseMtxSq<T, SIZE>& mtx, float eigenValue, int iterationCount)
{
	// first invert the matrix
	BaseMtxSq<T, SIZE> mtxInverse;
	{
		BaseMtx<T, SIZE * 2, SIZE> augmentedMtx;

		for (int row = 0; row < SIZE; ++row)
		{
			for (int column = 0; column < SIZE; ++column)
			{
				augmentedMtx[row][column] = mtx[row][column];
				if (row == column)
					augmentedMtx[row][column] -= eigenValue;

				augmentedMtx[row][column + SIZE] = (row == column) ? 1.0f : 0.0f;
			}
		}

		GaussJordanElimination(augmentedMtx);

		for (int row = 0; row < SIZE; ++row)
			for (int column = 0; column < SIZE; ++column)
				mtxInverse[row][column] = augmentedMtx[row][column + SIZE];
	}

	// iterate
	BaseVec<T, SIZE> eigenVector{};
	eigenVector[0] = 1.0f;
	for (int i = 0; i < iterationCount; ++i)
		eigenVector = Normalize(Multiply(mtxInverse, eigenVector));

	// return value
	return eigenVector;
}
