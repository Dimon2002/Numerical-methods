#pragma once

#include "DataConfiguration.h"
#include "Libraries.h"

class Matrix
{
public:
	std::string Errors = "";

	Matrix(
		const std::string& PathSize,
		const std::string& PathAu,
		const std::string& PathAl,
		const std::string& PathDi,
		const std::string& PathVector);
	~Matrix();

	void CalcLUStar();   // Разложение матрицы A на LU*
	void ForwardSubstitution(); // Решение уравнения методом прямого хода
	void BackwardSubstitution();    // Решение уравнения методом обратного хода
	
	void ShowSolution(); // Вывод результата решения СЛАУ
	void ShowDecompositions(); // Вывод Lu* Разложения

private:
	int n;
	int HalfTapeSize;

	type** al, ** au, * di; // Хранение матрицы
	type* x, * y, *b; // Резултат, Вектор прямого хода, Вектор правой части

	void InputSize(const std::string& PathSize);
	void InputMatrix(
		const std::string& PathAu,
		const std::string& PathAl,
		const std::string& PathDi);
	void InputVector(const std::string& PathVector);

	void AllocateMemory(); // Выделение памяти для хранения данных
	void ExistenceDecomposition(); // Проверка на существование разложения
};