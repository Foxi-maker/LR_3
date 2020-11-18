#include "Header.h"


Matrix::Matrix()
{
	A = nullptr;
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
	double* eigval = new double[matrixSize];

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
	delete[](eigval);
}

void Matrix::QR_EigVal_shifts(std::string flag)
{
	int operations = 0;
	double sigma = A[matrixSize - 1][matrixSize - 1];

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


		operations += matrixSize * matrixSize - 3 * matrixSize + 2;
	}
	else
	{
		OutToFile("QR with shifts\n");
	}

	for (int j = 0; j < matrixSize; j++)
		R[j][j] -= sigma;

	int numIter = 0;
	int dynamicSize = matrixSize;

	while (dynamicSize > 1)
	{
		//TODO: изменить структуру
		if (FUN::ConditionShifts(R, dynamicSize, flag))
		{
			for (int i = 0; i < dynamicSize; i++)
				R[i][i] += sigma;

			dynamicSize--;

			if (dynamicSize > 1)
			{
				sigma = R[dynamicSize - 1][dynamicSize - 1];

				for (int i = 0; i < dynamicSize; i++)
					R[i][i] -= sigma;
			}
		}
		else
		{
			QRMethod(Q, R, dynamicSize);
			operations += dynamicSize * dynamicSize - 3 * dynamicSize + 2;

			//результат умножения записывается в матрицу R
			FUN::MultMatrix(R, Q, dynamicSize);
			operations += dynamicSize * dynamicSize * dynamicSize;
		}
		numIter++;
	}
	double* eigval = new double[matrixSize];

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
	delete[](eigval);
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
		//Проверка на нулевой элемент на главной диагонали
		/*if (!R[j][j])
		{
			int index;
			for (index = j + 1; !R[index][j]; index++) {}

			std::swap(R[index], R[j]);
		}*/

		for (int i = 1 + j; i < n; i++)
		{
			c = R[j][j] / (sqrt(R[j][j] * R[j][j] + R[i][j] * R[i][j]));
			s = R[i][j] / (sqrt(R[j][j] * R[j][j] + R[i][j] * R[i][j]));

			double temp;

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

Matrix::~Matrix()
{
	if (A)
	{
		for (int i = 0; i < matrixSize; i++)
		{
			delete[](A[i]);
		}
		delete[](A);
	}
	if (stream.is_open())
		stream.close();
}