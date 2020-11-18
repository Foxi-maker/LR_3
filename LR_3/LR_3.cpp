#include "Header.h"

std::fstream Matrix::stream;

int main()
{

	for (const auto& f : fileNames)
	{
		Matrix elem(f);

		elem.Show();

		elem.OutStreamOpen(f);

		elem.OutToFile();

		elem.QR_EigVal();
		elem.QR_EigVal_shifts();
		elem.QR_EigVal_HF();

		elem.OutStreamClose();

		//elem.~Matrix();
	}
	
}
