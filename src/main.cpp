#define NUMCPP_NO_USE_BOOST

#include <NumCpp.hpp>
#include <Eigen>
//#include <numeric/ublas/vector.hpp>
//#include <numeric/ublas/matrix.hpp>
#include "KMatrix.h"

int main()
{
	// 行列 (1 2 3)
	// 　　 (2 4 6)
	// 　　 (3 6 9)

	//clock_t start = clock();    // スタート時間
	//boost::numeric::ublas::matrix<int> a(150, 10);
	//for (int i = 0; i < a.size1() ; ++i)
	//	for (int j = 0; j < a.size2(); ++j)
	//		a(i, j) = (i + 1) * (j + 1);
	//boost::numeric::ublas::matrix<int> a_ = prod(a, trans(a));
	//clock_t end = clock();     // 終了時間
	//std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	clock_t start = clock();    // スタート時間
	nc::NdArray<int> b(150, 10);
	for (int i = 0; i < b.numRows(); ++i)
	{
		for (int j = 0; j < b.numCols(); ++j)
		{
			b(i, j) = (i + 1) * (j + 1);
		}
	}
	nc::NdArray<int> b_ = matmul(b, b.transpose());
	clock_t end = clock();     // 終了時間
	std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	start = clock();    // スタート時間
	Eigen::MatrixXd c(150, 10);
	for (int i = 0; i < c.rows(); ++i)
		for (int j = 0; j < c.cols(); ++j)
			c(i, j) = (i + 1) * (j + 1);
	Eigen::MatrixXd c_ = c * c.transpose();
	end = clock();     // 終了時間
	std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
	for (int i = 0; i < c_.rows(); ++i)
	{
		for (int j = 0; j < c_.cols(); ++j)
		{
			std::cout << c_(i, j) << ",";
		}
		std::cout << std::endl;
	}

	start = clock();    // スタート時間
	KMatrix<int> a(150, 10);
	for (int i = 0; i < a.row(); ++i)
	{
		for (int j = 0; j < a.col(); ++j)
		{
			a(i, j) = (i + 1) * (j + 1);
		}
	}
	KMatrix<int> a_ = a * a.transpose();
	end = clock();     // 終了時間
	std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
	for (int i = 0; i < a_.row(); ++i)
	{
		for (int j = 0; j < a_.col(); ++j)
		{
			std::cout << a_(i, j) << ",";
		}
		std::cout << std::endl;
	}

	return 0;
}