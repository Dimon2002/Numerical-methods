#include "PrototypeMatrix.h"

#pragma region ���������� ������� A �� LU*

void Matrix::CalcLUStar()
{
	type AccumulatorAl = 0.0,
		AccumulatorAu = 0.0,
		AccumulatorD = 0.0;

	for (int i = 1; i < n; ++i)
	{
		AccumulatorD = 0.0;

		if (i < HalfTapeSize)
		{
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
		}
		else
		{
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
		}

		for (size_t j = 0; j < HalfTapeSize; j++)
			AccumulatorD += al[i][j] * au[i][j];
		di[i] -= AccumulatorD;
	}
}

#pragma endregion

#pragma region ������� ��������� ������� ������� ����

void Matrix::ForwardSubstitution()
{
	type Accumulator = 0;

	y[0] = b[0];
	for (size_t i = 1; i < n; ++i)
	{
		if (i < HalfTapeSize)
		{
			size_t k = HalfTapeSize - i;
			for (size_t j = 0; j < i; ++j) // ����� ������� ��������� ������� � ������ al
				Accumulator += al[i][k + j] * y[j];
			y[i] = b[i] - Accumulator;
		}
		else
		{
			for (size_t j = 0; j < HalfTapeSize; ++j)
				Accumulator += al[i][j] * y[j + i - HalfTapeSize];
			y[i] = b[i] - Accumulator;
		}
		Accumulator = 0;
	}
}

#pragma endregion

#pragma region ������� ��������� ������� ��������� ����

void Matrix::BackwardSubstitution()
{
	type Accumulator = 0;

	x[n - 1] = y[n - 1] / di[n - 1];

	for (int i = n - 1; i > 0; --i)
	{	
		int index = i - 1;
		if (i > n - HalfTapeSize)
		{
			for (int k = 0; k < n - i; ++k)
				Accumulator += au[i + k][HalfTapeSize - k - 1] * x[i + k];
		}
		else
		{
			for (int k = 0; k < HalfTapeSize; ++k)
				Accumulator += au[i + k][HalfTapeSize - k - 1] * x[i + k];
		}
		x[index] = (y[index] - Accumulator) / di[index];
		Accumulator = 0;
	}
}

#pragma endregion