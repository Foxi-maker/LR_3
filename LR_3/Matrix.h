#pragma once
class Matrix
{
	double** A;
	int matrixSize;

	//int operations = 0;
	static std::fstream stream;

	double* eigval;
public:
	Matrix();
	Matrix(const std::string& f);

	void Show() const;

	void OutStreamOpen(const std::string&);
	void OutStreamClose();
	void OutToFile()const;
	void OutToFile(const std::string&) const;
	void OutToFile(const std::string&, const double*)const;
	void OutToFile(const std::string&, double**)const;
	void OutToFile(const std::string&, double)const;

	void HesFormOut()const;

	void QR_EigVal(std::string = "SF");
	void QR_EigVal_shifts(std::string = "SF");

	void ReverseIterations();
	void QRforSLAE(double**, double**, int);

	void RI_Rayleigh();

	void QRMethod(double**, double**,int);

	~Matrix();
};

