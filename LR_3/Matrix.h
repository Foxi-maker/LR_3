#pragma once
class Matrix
{
	double** A;
	int matrixSize;

	static std::fstream stream;
public:
	Matrix();
	Matrix(const std::string& f);

	void Show() const;

	void OutStreamOpen(const std::string&);
	void OutStreamClose();
	void OutToFile()const;
	void OutToFile(const std::string&) const;
	void OutToFile(const std::string&, const double*)const;
	void OutToFile(const std::string&, double)const;

	void QR_EigVal(std::string);
	void QR_EigVal_shifts(std::string);
	void QR_EigVal_HF();
	void QR_EigVal_HF_shifts();

	void QRMethod(double**, double**,int);

	~Matrix();
};

