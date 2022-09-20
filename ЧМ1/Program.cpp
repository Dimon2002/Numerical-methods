#include "PrototypeMatrix.h"

const std::string PathSize = "SizeMatrix";
const std::string PathAu = "Au";
const std::string PathAl = "Al";
const std::string PathDi = "Di";
const std::string PathVector = "V";

int main()
{
	setlocale(LC_ALL, "ru");
	Matrix a(PathSize, PathAu, PathAl, PathDi, PathVector);

	const auto status = a.Errors;
	if (status.size() != 0)
		goto stop;

	a.CalcLUStar();
	a.ForwardSubstitution();
	a.BackwardSubstitution();
	a.ShowSolution();
	//a.ShowDecompositions();
	return 0;

stop:
	std::cout << status << std::endl;
}