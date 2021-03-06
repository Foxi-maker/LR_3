#pragma once
#include "Header.h"


void FUN::Show(double** matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << std::setw(12) << matrix[i][j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void FUN::MultMatrix(double** a, double** b, int n)
{
	double** c = new double*[n];
	for (int i = 0; i < n; i++)
		c[i] = new double[n];

	for (int indexRow = 0; indexRow < n; indexRow++)
	{
		for (int indexCol = 0; indexCol < n; indexCol++)
		{
			double temp = 0;
			for (int index = 0; index < n; index++)
			{
				temp += a[indexRow][index] * b[index][indexCol];
			}
			c[indexRow][indexCol] = temp;
		}
	}

	for (int indexRow = 0; indexRow < n; indexRow++)
		for (int indexCol = 0; indexCol < n; indexCol++)
			a[indexRow][indexCol] = c[indexRow][indexCol];

	for (int i = 0; i < n; i++)
	{
		delete[](c[i]);
	}
	delete[](c);
}

bool FUN::StopCondition(double ** a, int n, std::string flag)
{
	double max = 0;

	if (flag == "HF")
	{
		for (int i = 1; i < n; i++)
			if (fabs(a[i][i - 1]) > max)
				max = fabs(a[i][i - 1]);
	}
	else
	{
		for (int i = 1; i < n; i++)
			for (int j = 0; j < i; j++)
				if (fabs(a[i][j]) > max)
					max = fabs(a[i][j]);
	}


	if (max < eps)
		return true;
	else
		return false;
}

bool FUN::ConditionShifts(double** a, int n, std::string flag)
{
	double max = 0;
	if (flag == "HF")
	{
		max = fabs(a[n - 1][n - 2]);
	}
	else
	{
		for (int i = 0; i < n - 1; i++)
			if (fabs(a[n - 1][i]) > max)
				max = fabs(a[n - 1][i]);
	}

	if (max < eps)
		return true;
	else
		return false;
}

void FUN::HesForm(double ** a, int n)
{
	double c, s, znam;
	//k ���� �� �������
	for (int k = 1; k < n - 1; k++)
	{
		//i ���� �� ��������
		for (int l = k + 1; l < n; l++)
		{
			znam = (sqrt(a[k][k - 1] * a[k][k - 1] + a[l][k - 1] * a[l][k - 1]));
			c = a[k][k - 1] / znam;
			s = a[l][k - 1] / znam;

			double temp;

			//��������� ����� �� ������� ��������
			for (int columnIterator = k; columnIterator < n + 1; columnIterator++)
			{
				temp = a[k][columnIterator - 1];
				a[k][columnIterator - 1] = c * a[k][columnIterator - 1] + s * a[l][columnIterator - 1];
				a[l][columnIterator - 1] = -s * temp + c * a[l][columnIterator - 1];
			}
			//��������� ������ �� ����������������� ������� ��������
			for (int columnIterator = k; columnIterator < n + 1; columnIterator++)
			{
				temp = a[columnIterator - 1][k];
				a[columnIterator - 1][k] = c * a[columnIterator - 1][k] + s * a[columnIterator - 1][l];
				a[columnIterator - 1][l] = -s * temp + c * a[columnIterator - 1][l];
			}
		}
	}
}

double FUN::NormInf(double* x, int n)
{
	int indexMax = 0;
	for (int index = 0; index < n; index++)
	{
		if (fabs(x[indexMax]) < fabs(x[index]))
		{
			indexMax = index;
		}
	}
	return abs(x[indexMax]);
}