//#define NUMCPP_NO_USE_BOOST
//#define VALUE_CHECK
#define SPEED_TEST

#ifdef SPEED_TEST
#include <queue>
#include <vector>
#endif

#include <iostream>
#include <chrono>
#include <fstream>
//#include <NumCpp.hpp>
#include <Eigen>
#include "KMatrix.h"

int main()
{
#ifdef SPEED_TEST
	std::ofstream speedlog("speed_log.csv");
	enum FUNCTION_NAME
	{
		name,
		constructor = 1,
		//size,
		//row,
		//col,
		//data,
		rowdata,
		coldata,
		amax,
		amin,
		argmax,
		argmin,
		transpose,
		zeros,
		reserve,
		resize,
		diag,
		ones,
		abs,
		multiply,
		mean,
		exp,
		//arrayaccess,
		matmul,
		//equal,
		flip,
		reverse,
		MAX
	};
	std::vector<std::string> FUNCTION_LIST;
	{
		FUNCTION_LIST.push_back("[nanoseconds]");
		FUNCTION_LIST.push_back("constructor");
		//FUNCTION_LIST.push_back("size");
		//FUNCTION_LIST.push_back("row");
		//FUNCTION_LIST.push_back("col");
		//FUNCTION_LIST.push_back("data");
		FUNCTION_LIST.push_back("rowdata");
		FUNCTION_LIST.push_back("coldata");
		FUNCTION_LIST.push_back("amax");
		FUNCTION_LIST.push_back("amin");
		FUNCTION_LIST.push_back("argmax");
		FUNCTION_LIST.push_back("argmin");
		FUNCTION_LIST.push_back("transpose");
		FUNCTION_LIST.push_back("zeros");
		FUNCTION_LIST.push_back("reserve");
		FUNCTION_LIST.push_back("resize");
		FUNCTION_LIST.push_back("diag");
		FUNCTION_LIST.push_back("ones");
		FUNCTION_LIST.push_back("abs");
		FUNCTION_LIST.push_back("multiply");
		FUNCTION_LIST.push_back("mean");
		FUNCTION_LIST.push_back("exp");
		//FUNCTION_LIST.push_back("arrayaccess");
		FUNCTION_LIST.push_back("matmul");
		//FUNCTION_LIST.push_back("equal");
		FUNCTION_LIST.push_back("flip");
		FUNCTION_LIST.push_back("reverse");
	};
	for (unsigned int i = 0; i < FUNCTION_LIST.size(); ++i)
	{
		if (i != FUNCTION_LIST.size() - 1)
		{
			speedlog << FUNCTION_LIST[i] << ",";
		}
		else
		{
			speedlog << FUNCTION_LIST[i] << std::endl;
		}
	}
	std::vector<double> eigen_speed(MAX);
	std::vector<std::string> matrix_speed(MAX);
	matrix_speed[name] = "KMatrix";
#endif
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
	std::chrono::system_clock::time_point start, end;

	// =======================================================================================
	// NumCpp
	// =======================================================================================
	//start = std::chrono::system_clock::now(); // 計測開始時間
	//nc::NdArray<int> b(150, 10);
	//for (int i = 0; i < b.numRows(); ++i)
	//{
	//	for (int j = 0; j < b.numCols(); ++j)
	//	{
	//		b(i, j) = (i + 1) * (j + 1);
	//	}
	//}
	//nc::NdArray<int> b_ = matmul(b, b.transpose());
	//end = std::chrono::system_clock::now();  // 計測終了時間
	//double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	//std::cout << "duration = " << elapsed << "nanosec.\n";

	// =======================================================================================
	// Eigen
	// =======================================================================================
	start = std::chrono::system_clock::now(); // 計測開始時間
	Eigen::MatrixXd c(5, 10);
	for (int i = 0; i < c.rows(); ++i)
		for (int j = 0; j < c.cols(); ++j)
			c(i, j) = (i + 1) * (j + 1);
	end = std::chrono::system_clock::now();  // 計測終了時間
	double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	eigen_speed[constructor] = elapsed;
#ifdef VALUE_CHECK
	std::cout << "Eigen::constructor" << std::endl;
	for (int i = 0; i < c.rows(); ++i)
	{
		for (int j = 0; j < c.cols(); ++j)
		{
			std::cout << c(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	Eigen::MatrixXd c_trans = c.transpose();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	eigen_speed[transpose] = elapsed;
#ifdef VALUE_CHECK
	std::cout << "Eigen::transpose" << std::endl;
	for (int i = 0; i < c_trans.rows(); ++i)
	{
		for (int j = 0; j < c_trans.cols(); ++j)
		{
			std::cout << c_trans(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	Eigen::MatrixXd c_ = c * c_trans;
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	eigen_speed[matmul] = elapsed;
#ifdef VALUE_CHECK
	std::cout << "Eigen::matmul" << std::endl;
	for (int i = 0; i < c_.rows(); ++i)
	{
		for (int j = 0; j < c_.cols(); ++j)
		{
			std::cout << c_(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	// =======================================================================================
	// KMatrix
	// =======================================================================================
	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a(5, 10);
	for (int i = 0; i < a.row(); ++i)
	{
		for (int j = 0; j < a.col(); ++j)
		{
			a(i, j) = (i + 1) * (j + 1);
		}
	}
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[constructor] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::constructor" << std::endl;
	for (int i = 0; i < a.row(); ++i)
	{
		for (int j = 0; j < a.col(); ++j)
		{
			std::cout << a(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_flip = a.flip();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[flip] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::flip" << std::endl;
	for (int i = 0; i < a_flip.row(); ++i)
	{
		for (int j = 0; j < a_flip.col(); ++j)
		{
			std::cout << a_flip(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_trans = a.transpose();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[transpose] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::transpose" << std::endl;
	for (int i = 0; i < a_trans.row(); ++i)
	{
		for (int j = 0; j < a_trans.col(); ++j)
		{
			std::cout << a_trans(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_ = a * a_trans;
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[matmul] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::matmul" << std::endl;
	for (int i = 0; i < a_.row(); ++i)
	{
		for (int j = 0; j < a_.col(); ++j)
		{
			std::cout << a_(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_multiply = a.multiply(a);
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[multiply] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::multiply" << std::endl;
	for (int i = 0; i < a_multiply.row(); ++i)
	{
		for (int j = 0; j < a_multiply.col(); ++j)
		{
			std::cout << a_multiply(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_rowdata = a_.rowdata(1);
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[rowdata] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::rowdata" << std::endl;
	for (int i = 0; i < a_rowdata.row(); ++i)
	{
		for (int j = 0; j < a_rowdata.col(); ++j)
		{
			std::cout << a_rowdata(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_reverse = a_.rowdata(1).reverse();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[reverse] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::reverse" << std::endl;
	for (int i = 0; i < a_reverse.row(); ++i)
	{
		for (int j = 0; j < a_reverse.col(); ++j)
		{
			std::cout << a_reverse(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_coldata = a_.coldata(2);
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[coldata] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::coldata" << std::endl;
	for (int i = 0; i < a_coldata.row(); ++i)
	{
		for (int j = 0; j < a_coldata.col(); ++j)
		{
			std::cout << a_coldata(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	int a_amax = a_.rowdata(2).amax();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[amax] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::amax" << std::endl;
	std::cout << a_amax << std::endl;
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	int a_amin = a_.rowdata(2).amin();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[amin] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::amin" << std::endl;
	std::cout << a_amin << std::endl;
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	int a_argmax = a_.coldata(1).argmax();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[argmax] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::argmax" << std::endl;
	std::cout << a_argmax << std::endl;
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	int a_argmin = a_.coldata(1).argmin();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[argmin] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::argmin" << std::endl;
	std::cout << a_argmin << std::endl;
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	int a_argmean = a_.coldata(2).mean();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[mean] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::mean" << std::endl;
	std::cout << a_argmean << std::endl;
#endif

	KMatrix<int> a2(5, 10);
	for (int i = 0; i < a2.row(); ++i)
	{
		for (int j = 0; j < a2.col(); ++j)
		{
			a2(i, j) = std::pow(-1 * i, j);
		}
	}
#ifdef VALUE_CHECK
	std::cout << "KMatrix::a2" << std::endl;
	for (int i = 0; i < a2.row(); ++i)
	{
		for (int j = 0; j < a2.col(); ++j)
		{
			std::cout <<
				a2(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif
	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_abs = a2.abs();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[abs] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::abs" << std::endl;
	for (int i = 0; i < a_abs.row(); ++i)
	{
		for (int j = 0; j < a_abs.col(); ++j)
		{
			std::cout <<
				a_abs(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	KMatrix<double> a3(5, 10);
	for (int i = 0; i < a3.row(); ++i)
	{
		for (int j = 0; j < a3.col(); ++j)
		{
			a3(i, j) = i + j;
		}
	}
#ifdef VALUE_CHECK
	std::cout << "KMatrix::a3" << std::endl;
	for (int i = 0; i < a3.row(); ++i)
	{
		for (int j = 0; j < a3.col(); ++j)
		{
			std::cout <<
				a3(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif
	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<double> a_exp = a3.exp();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[exp] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::exp" << std::endl;
	for (int i = 0; i < a_exp.row(); ++i)
	{
		for (int j = 0; j < a_exp.col(); ++j)
		{
			std::cout <<
				a_exp(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_coldiag = a_.coldata(1).diag();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[diag] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::diag" << std::endl;
	for (int i = 0; i < a_coldiag.row(); ++i)
	{
		for (int j = 0; j < a_coldiag.col(); ++j)
		{
			std::cout <<
				a_coldiag(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	KMatrix<int> a_ones(5, 7);
	a_ones.ones();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[ones] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::ones" << std::endl;
	for (int i = 0; i < a_ones.row(); ++i)
	{
		for (int j = 0; j < a_ones.col(); ++j)
		{
			std::cout <<
				a_ones(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	a.resize(25, 2);
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[resize] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::resize" << std::endl;
	for (int i = 0; i < a.row(); ++i)
	{
		for (int j = 0; j < a.col(); ++j)
		{
			std::cout << a(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	a.zeros();
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[zeros] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::zeros" << std::endl;
	for (int i = 0; i < a.row(); ++i)
	{
		for (int j = 0; j < a.col(); ++j)
		{
			std::cout << 
				a(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	start = std::chrono::system_clock::now(); // 計測開始時間
	a_.reserve(3, 5);
	end = std::chrono::system_clock::now();  // 計測終了時間
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); //処理に要した時間をnano秒に変換
	std::cout << "duration = " << elapsed << "nanosec.\n";
	matrix_speed[reserve] = std::to_string(elapsed);
#ifdef VALUE_CHECK
	std::cout << "KMatrix::reserve" << std::endl;
	for (int i = 0; i < a_.row(); ++i)
	{
		for (int j = 0; j < a_.col(); ++j)
		{
			std::cout << a_(i, j) << ",";
		}
		std::cout << std::endl;
	}
#endif

	for (unsigned int i = 0; i < MAX; ++i)
	{
		if (i != MAX - 1)
		{
			speedlog << matrix_speed[i] << ",";
		}
		else
		{
			speedlog << matrix_speed[i] << std::endl;
		}
	}
	speedlog.close();

	return 0;
}