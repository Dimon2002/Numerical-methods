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

	void CalcLUStar();   // ���������� ������� A �� LU*
	void ForwardSubstitution(); // ������� ��������� ������� ������� ����
	void BackwardSubstitution();    // ������� ��������� ������� ��������� ����

	void ShowSolution(); // ����� ���������� ������� ����
	void ShowDecompositions(); // ����� Lu* ����������

private:
	int n;
	int HalfTapeSize;

	type** al, ** au, * di; // �������� �������
	type* x, * y, * b; // ��������, ������ ������� ����, ������ ������ �����

	void InputSize(const std::string& PathSize);
	void InputMatrix(
		const std::string& PathAu,
		const std::string& PathAl,
		const std::string& PathDi);
	void InputVector(const std::string& PathVector);

	void AllocateMemory(); // ��������� ������ ��� �������� ������
	void ExistenceDecomposition(); // �������� �� ������������� ����������
};