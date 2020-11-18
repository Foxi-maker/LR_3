#pragma once

namespace FUN
{
	void Show(double**, int);

	void MultMatrix(double**, double**, int);

	bool StopCondition(double**, int, std::string = "SF");

	bool ConditionShifts(double**, int, std::string = "SF");

	void HesForm(double**, int);

	//void Difference(double**, int, double*, double*, int);

	double NormInf(double*,int);
}