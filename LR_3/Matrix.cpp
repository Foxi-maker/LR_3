#include "Header.h"


Matrix::Matrix()
{
	A = nullptr;
	eigval = nullptr;
}

Matrix::Matrix(const std::string& f)
{
	stream.open("Tests/" + f, std::ios_base::in);

	if (!stream.is_open())
	{
		std::cerr << "The file was not opened!\n";
		exit(1);
	}

	stream >> matrixSize;

	A = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		A[i] = new double[matrixSize];

	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			stream >> A[i][j];
		}
		stream.get();
	}
	stream.close();

	eigval = new double[matrixSize];
}

void Matrix::Show() const
{
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			std::cout << std::setw(5) << A[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void Matrix::OutStreamOpen(const std::string& f)
{
	stream.open("Results/" + f, std::ios_base::out);
	if (!stream.is_open())
	{
		std::cout << "File " << f << " was not opened for writing!";
		exit(1);
	}
	OutToFile("eps = ", eps);
}

void Matrix::OutStreamClose()
{
	if (stream.is_open())
		stream.close();
}

//Запись в файл *proj dir*\Results
void Matrix::OutToFile()const
{
	stream << matrixSize << "\n";
	stream.width(10);
	for (int i = 0; i < matrixSize; ++i)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			stream << std::setw(8) << std::left << A[i][j] << " ";
		}
		stream << "\n";
	}
	stream << "\n";
}

void Matrix::OutToFile(const std::string& str) const
{
	stream << str;
}

void Matrix::OutToFile(const std::string& str, const double* x)const
{
	stream << str;
	for (int index = 0; index < matrixSize; index++)
		stream << std::setw(12) << std::left << x[index] << " ";
	stream << "\n";
}

void Matrix::OutToFile(const std::string& str, double** matrix)const
{
	stream << str << "\n";
	for (int i = 0; i < matrixSize; ++i)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			stream << std::setw(12) << std::left << matrix[i][j] << " ";
		}
		stream << "\n";
	}
	stream << "\n";
}

void Matrix::OutToFile(const std::string& str, double num)const
{
	stream << str << num << "\n";
}

void Matrix::HesFormOut() const
{
	double** R = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		R[i] = new double[matrixSize];

	for (int i = 0; i < matrixSize; ++i)
		for (int j = 0; j < matrixSize; j++)
			R[i][j] = A[i][j];

	FUN::HesForm(R, matrixSize);

	OutToFile("Hessenberg decomposition:", R);

	//FUN::Show(R, matrixSize);

	for (int i = 0; i < matrixSize; i++)
	{
		delete[](R[i]);
	}
	delete[](R);
}

void Matrix::QR_EigVal(std::string flag)
{
	int operations = 0;
	double** R = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		R[i] = new double[matrixSize];

	for (int i = 0; i < matrixSize; ++i)
		for (int j = 0; j < matrixSize; j++)
			R[i][j] = A[i][j];

	double** Q = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		Q[i] = new double[matrixSize];

	if (flag == "HF")
	{
		OutToFile("QR (Hessenberg form)\n");
		FUN::HesForm(R, matrixSize);
		//число вращений 
		operations += matrixSize * matrixSize - 3 * matrixSize + 2;

	}
	else
	{
		OutToFile("QR\n");
	}


	int numIter = 0;
	while (!FUN::StopCondition(R, matrixSize, flag))
	{
		numIter++;

		QRMethod(Q, R, matrixSize);

		operations += matrixSize * matrixSize - 3 * matrixSize + 2;

		//результат умножения записывается в матрицу R
		FUN::MultMatrix(R, Q, matrixSize);
		operations += matrixSize * matrixSize * matrixSize;

		//std::cout << "A:\n";
		//FUN::Show(R, matrixSize);

	}

	for (int i = 0; i < matrixSize; i++)
		eigval[i] = R[i][i];

	//Вывод в файл числа итераций и собственных чисел
	OutToFile("Number of iterations: ", numIter);
	OutToFile("Number of spins + multiplications: ", operations);
	OutToFile("Eigenvalues: ", eigval);
	OutToFile("\n");

	for (int i = 0; i < matrixSize; i++)
	{
		delete[](R[i]);
		delete[](Q[i]);
	}
	delete[](Q);
	delete[](R);
	
}

void Matrix::QR_EigVal_shifts(std::string flag)
{
	int operations = 0;

	double** R = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		R[i] = new double[matrixSize];

	for (int i = 0; i < matrixSize; ++i)
		for (int j = 0; j < matrixSize; j++)
			R[i][j] = A[i][j];

	double** Q = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		Q[i] = new double[matrixSize];

	if (flag == "HF")
	{
		OutToFile("QR with shifts (Hessenberg form)\n");
		FUN::HesForm(R, matrixSize);

		//std::cout << "Hessenberg form:\n";
		//FUN::Show(R, matrixSize);

		operations += matrixSize * matrixSize - 3 * matrixSize + 2;
	}
	else
	{
		OutToFile("QR with shifts\n");
	}

	int numIter = 0;
	int dynamicSize = matrixSize;

	double sigma;

	while (dynamicSize > 1)
	{
		sigma = R[dynamicSize - 1][dynamicSize - 1];

		for (int j = 0; j < dynamicSize; j++)
			R[j][j] -= sigma;

		while (!FUN::ConditionShifts(R, dynamicSize, flag))
		{
		//	std::cout << "As:\n";
		//	FUN::Show(R, dynamicSize);

			QRMethod(Q, R, dynamicSize);
			operations += dynamicSize * dynamicSize - 3 * dynamicSize + 2;

			//результат умножения записывается в матрицу R
			FUN::MultMatrix(R, Q, dynamicSize);
			operations += dynamicSize * dynamicSize * dynamicSize;

			//std::cout << "A:\n";
			//FUN::Show(R, dynamicSize);

			//numIter++;
			//std::cout << "numIter: " << numIter << "\n";
		}

		for (int i = 0; i < dynamicSize; i++)
			R[i][i] += sigma;

		dynamicSize--;
	}


	for (int i = 0; i < matrixSize; i++)
		eigval[i] = R[i][i];

	//Вывод в файл числа итераций и собственных чисел
	OutToFile("Number of iterations: ", numIter);
	OutToFile("Number of spins + multiplications: ", operations);
	OutToFile("Eigenvalues: ", eigval);
	OutToFile("\n");

	for (int i = 0; i < matrixSize; i++)
	{
		delete[](R[i]);
		delete[](Q[i]);
	}
	delete[](Q);
	delete[](R);
}

void Matrix::QRMethod(double** Q, double** R, int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				Q[i][j] = 1.;
			else
				Q[i][j] = 0.;
		}
	double c, s;

	for (int j = 0; j < n - 1; j++)
	{
		//Проверка на нулевой элемент на главной диагонали?
		//

		for (int i = 1 + j; i < n; i++)
		{
			double temp = sqrt(R[j][j] * R[j][j] + R[i][j] * R[i][j]);
			c = R[j][j] / temp;
			s = R[i][j] / temp;

			//Умножение на матрицу поворота
			for (int index = 0; index < n; index++)
			{
				temp = Q[j][index];
				Q[j][index] = c * Q[j][index] + s * Q[i][index];
				Q[i][index] = -s * temp + c * Q[i][index];
			}

			for (int columnIterator = j; columnIterator < n; columnIterator++)
			{
				temp = R[j][columnIterator];
				R[j][columnIterator] = c * R[j][columnIterator] + s * R[i][columnIterator];
				R[i][columnIterator] = -s * temp + c * R[i][columnIterator];
			}
		}
	}

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			std::swap(Q[i][j], Q[j][i]);
		}
	}
}

void Matrix::ReverseIterations()
{
	OutToFile("Reverse iterations\n");
	double** x = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		x[i] = new double[matrixSize];

	for (int i = 0; i < matrixSize; ++i)
		for (int j = 0; j < matrixSize; j++)
			if (i == j)
				x[i][j] = 1.;
			else
				x[i][j] = 0.;
	
	double** B = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)
		B[i] = new double[matrixSize];

	double* x_K = new double[matrixSize];
	double* dif = new double[matrixSize];

	//double* numIter = new double*[matrixSize];
	int numIter;

	for (int index=0;index<matrixSize;index++)
	{
		numIter = 0;

		for (int i = 0; i < matrixSize; ++i)
		{
			x_K[i] = 0.;

			for (int j = 0; j < matrixSize; j++)
				if (i == j)
					B[i][j] = A[i][j] - eigval[index];
				else
					B[i][j] = A[i][j];
		}
			
		
		do{
			QRforSLAE(B, x, index);

			//Нормировка ||*||_2
			double Norm_x = 0.;
			for (int i = 0; i < matrixSize; ++i)
				Norm_x += x[index][i]* x[index][i];

			Norm_x = std::sqrt(Norm_x);

			for (int i = 0; i < matrixSize; ++i)
				x[index][i] /= Norm_x;

			numIter++;

			for (int i = 0; i < matrixSize; i++)
				dif[i] = fabs(x[index][i]) - fabs(x_K[i]);

			FUN::Show(x, matrixSize);
			for (int i = 0; i < matrixSize; i++)
				x_K[i] = x[index][i];
			std::cout << "FUN::NormInf(dif, matrixSize): " << FUN::NormInf(dif, matrixSize) << "\n";
		}while (FUN::NormInf(dif, matrixSize) > eps);

		stream << index + 1 <<" eigenvector: ";
		OutToFile("", x[index]);
		OutToFile("Number of iterations: ", numIter);
	}
	OutToFile("\n");

	for (int i = 0; i < matrixSize; i++)
	{
		delete[](x[i]);
		delete[](B[i]);
	}
	delete[](x);
	delete[](B);
	delete[](x_K);
	delete[](dif);
}

void Matrix::QRforSLAE(double** A, double** e, int k)
{
	double c, s;
	for (int j = 0; j < matrixSize - 1; j++)
	{
		if (!A[j][j])
		{
			int index;
			for (index = j + 1; !A[index][j]; index++) {}

			std::swap(A[index], A[j]);
		}
		for (int i = 1 + j; i < matrixSize; i++)
		{
			double temp = sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]);
			c = A[j][j] / temp;
			s = A[i][j] / temp;

			temp = e[k][j];
			e[k][j] = c * e[k][j] + s * e[k][i];
			e[k][i] = -s * temp + c * e[k][i];
			//Умножение на матрицу поворота
			for (int columnIterator = j; columnIterator < matrixSize; columnIterator++)
			{
				temp = A[j][columnIterator];
				A[j][columnIterator] = c * A[j][columnIterator] + s * A[i][columnIterator];
				A[i][columnIterator] = -s * temp + c * A[i][columnIterator];
			}
		}
	}
	//FUN::Show(A,matrixSize);

	for (int columIterator = matrixSize - 1; columIterator >= 0; columIterator--)
	{
		for (int rowIterator = matrixSize - 1; rowIterator > columIterator; rowIterator--)
		{
			e[k][columIterator] -= A[columIterator][rowIterator] * e[k][rowIterator];
		}
		e[k][columIterator] /= A[columIterator][columIterator];
	}
}

Matrix::~Matrix()
{
	if (A)
	{
		for (int i = 0; i < matrixSize; i++)
		{
			delete[](A[i]);
		}
		delete[](A);
		delete[](eigval);
	}
	if (stream.is_open())
		stream.close();
}