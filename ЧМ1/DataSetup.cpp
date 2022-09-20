#include "PrototypeMatrix.h"

Matrix::Matrix(
	const std::string& PathSize,
	const std::string& PathAu,
	const std::string& PathAl,
	const std::string& PathDi,
	const std::string& PathVector)
{
	InputSize(PathSize);
	AllocateMemory();
	InputMatrix(PathAu, PathAl, PathDi);
	ExistenceDecomposition();
	InputVector(PathVector);
}

Matrix::~Matrix()
{
	delete[] x;
	delete[] y;
	delete[] b;
	delete[] di;

	for (size_t i = 0; i < n; i++)
	{
		delete[] au[i];
		delete[] al[i];
	}
}

void Matrix::InputSize(const std::string& PathSize) {

	std::ifstream ReaderSize(PathSize + ".txt");
	if (!ReaderSize.is_open()) {
		Errors = "���� SizeMatrix.txt �� ������!";
		return;
	}
	ReaderSize >> n >> HalfTapeSize;
	ReaderSize.close();
	if (HalfTapeSize > n) {
		Errors = "�������� ������� ������� � ����������!";
		return;
	}
}

void Matrix::AllocateMemory()
{
	x = new type[n];
	y = new type[n];
	b = new type[n];

	di = new type[n];
	al = new type * [n];
	au = new type * [n];
	for (size_t i = 0; i < n; i++)
	{
		al[i] = new type[HalfTapeSize];
		au[i] = new type[HalfTapeSize];
	}
}

void Matrix::InputMatrix(
	const std::string& PathAu,
	const std::string& PathAl,
	const std::string& PathDi)
{
	if (Errors != "") return;

	std::ifstream ReaderAu(PathAu + ".txt");
	std::ifstream ReaderAl(PathAl + ".txt");
	std::ifstream ReaderDi(PathDi + ".txt");

	if (!ReaderAu.is_open()) {
		ReaderAu.close();
		Errors = "���� Au.txt �� ������!";
		return;
	}
	if (!ReaderAl.is_open()) {
		ReaderAl.close();
		Errors = "���� Al.txt �� ������!";
		return;
	}
	if (!ReaderDi.is_open()) {
		ReaderDi.close();
		Errors = "���� Di.txt �� ������!";
		return;
	}
	for (size_t i = 0; i < n; i++)
	{
		ReaderDi >> di[i];
		for (size_t j = 0; j < HalfTapeSize; j++)
		{
			ReaderAl >> al[i][j];
			ReaderAu >> au[i][j];
		}
	}
	ReaderAu.close();
	ReaderAl.close();
	ReaderDi.close();
}

void Matrix::InputVector(const std::string& PathVector) {

	if (Errors != "") return;
	std::ifstream ReaderVector(PathVector + ".txt");
	if (!ReaderVector.is_open()) {
		ReaderVector.close();
		Errors = "���� V.txt �� ������!";
		return;
	}
	for (size_t i = 0; i < n; i++)
		ReaderVector >> b[i];
	ReaderVector.close();
}

void Matrix::ExistenceDecomposition()
{
	if (Errors != "") return;
	for (size_t i = 0; i < n; i++)
		if (di[i] == 0)
		{
			Errors = "������� �� ����� LU* ����������";
			return;
		}
}

void Matrix::ShowSolution()
{
	std::cout << std::fixed << std::setprecision(6);
	for (size_t i = 0; i < n; i++)
		std::cout << x[i] << std::endl;
}

void Matrix::ShowDecompositions()
{
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "L:" << std::endl;
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < HalfTapeSize; ++j)
		{
			std::cout << al[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	std::cout << "U:" << std::endl;
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < HalfTapeSize; ++j)
		{
			std::cout << au[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	std::cout << "Ud:" << std::endl;
	for (size_t i = 0; i < n; i++)
	{
		std::cout << di[i] << std::endl;
	}
}