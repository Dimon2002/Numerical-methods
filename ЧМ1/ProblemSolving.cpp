#include "PrototypeMatrix.h"

#pragma region Разложение матрицы A на LU*

void Matrix::CalcLUStar()
{
	type AccumulatorAl = 0.0,
		AccumulatorAu = 0.0,
		AccumulatorD = 0.0;

	// Частный случай с обработкой дополненных нулей в al и au
	for (size_t i = 1; i < HalfTapeSize; ++i)
	{
		type AccumulatorD = 0.0;

		for (size_t j = HalfTapeSize - i; j < HalfTapeSize; ++j)
		{
			AccumulatorAl = 0.0;
			AccumulatorAu = 0.0;

			for (size_t k = 0; k < j - (HalfTapeSize - i); ++k)
			{
				AccumulatorAl += al[i][HalfTapeSize - i + k] * au[j - (HalfTapeSize - i)][(static_cast<unsigned long long>(2) * HalfTapeSize) - j - i + k];
				AccumulatorAu += al[HalfTapeSize - i + k][i] * au[(static_cast<unsigned long long>(2) * HalfTapeSize) - j - i + k][j - (HalfTapeSize - i)];
			}

			al[i][j] -= AccumulatorAl;
			al[i][j] /= di[j - (HalfTapeSize - i)];
			au[i][j] -= AccumulatorAu;
		}
		// Для диагонали
		for (size_t k = HalfTapeSize - i; k < HalfTapeSize; ++k)
		{
			AccumulatorD += al[i][k] * au[i][k];
		}
		di[i] -= AccumulatorD;
	}

	// Основное разложение
	for (int i = HalfTapeSize; i < n; ++i)
	{
		AccumulatorD = 0.0;

		for (int j = 0; j < HalfTapeSize; ++j)
		{
			AccumulatorAl = 0.0;
			AccumulatorAu = 0.0;

			for (size_t k = 1; k <= j; ++k)
			{
				AccumulatorAl += al[i][k - 1] * au[j - (HalfTapeSize - i)][HalfTapeSize - static_cast<unsigned long long>(j) + k - 1];
				AccumulatorAu += al[j - (HalfTapeSize - i)][HalfTapeSize - static_cast<unsigned long long>(j) + k - 1] * au[i][k - 1];
			}
			al[i][j] -= AccumulatorAl;
			al[i][j] /= di[i - (HalfTapeSize - j)];
			au[i][j] -= AccumulatorAu;
		}

		// Диагональ
		for (size_t j = 0; j < HalfTapeSize; j++)
			AccumulatorD += al[i][j] * au[i][j];
		di[i] -= AccumulatorD;
	}
}

#pragma endregion

#pragma region Решение уравнения методом прямого хода

void Matrix::ForwardSubstitution()
{
	type Accumulator = 0;

	y[0] = b[0];
	for (size_t i = 1; i < HalfTapeSize; ++i)
	{
		Accumulator = 0;
		size_t k = HalfTapeSize - i;
		for (size_t j = 0; j < i; ++j) // Взять столько элементов сколько в строке al
			Accumulator += al[i][k + j] * y[j];
		y[i] = b[i] - Accumulator;
	}
	// Все остальные элементы матрицы al
	for (size_t i = HalfTapeSize; i < n; ++i)
	{
		Accumulator = 0;
		for (size_t j = 0; j < HalfTapeSize; ++j)
			Accumulator += al[i][j] * y[j + i - HalfTapeSize];
		y[i] = b[i] - Accumulator;
	}
}

#pragma endregion

#pragma region Решение уравнения методом обратного хода

void Matrix::BackwardSubstitution()
{
	type Accumulator = 0;

	x[n - 1] = y[n - 1] / di[n - 1];

	for (int i = n - 1; i >= n - HalfTapeSize; --i)
	{
		for (int k = 0; k < n - i; ++k)
			Accumulator += au[i + k][HalfTapeSize - k - 1] * x[i + k];
		int index = i - 1;
		x[index] = (y[index] - Accumulator) / di[index];
		Accumulator = 0;
	}

	for (int i = n - HalfTapeSize; i > 0; --i)
	{
		for (int k = 0; k < HalfTapeSize; ++k)
			Accumulator += au[i + k][HalfTapeSize - k - 1] * x[i + k];
		int index = i - 1;
		x[index] = (y[index] - Accumulator) / di[index];
		Accumulator = 0;
	}
}

#pragma endregion