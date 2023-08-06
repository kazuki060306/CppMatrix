

/*
 * 行列演算クラス
 * ToDo; 逆行列生成
 */
template<typename arraytype>
class KMatrix
{
private:

	int m_Rows;			// 行数
	int m_Cols;			// 列数
	arraytype* m_data;	// 配列データ

public:

	/*
	 * コンストラクタ
	 */
	KMatrix(int row, int col)
		: m_Rows(row)
		, m_Cols(col)
	{
		m_data = new arraytype[row * col];
	}
	 
	 
	 /*
	 * 行数の取得 
	 */
	const int row()
	{
		return this->m_Rows;
	}

	/*
	 * 列数の取得
	 */
	const int col()
	{
		return this->m_Cols;
	}

	/*
	 * 配列データの取得
	 */
	arraytype* data()
	{
		return this->m_data;
	}

	/*
	 * 転置行列の取得
	 */
	KMatrix<arraytype>& transpose()
	{
		KMatrix<arraytype> transmat(m_Cols, m_Rows);
		for (int i = 0; i < m_Rows; ++i)
		{
			for (int j = 0; j < m_Cols; ++j)
			{
				transmat(j, i) = operator()(i, j);
			}
		}
		return transmat;
	}

	/*
	 * ゼロベクトルの取得
	 */
	KMatrix<arraytype>& zeros(int row, int col)
	{
		KMatrix<arraytype> zaromat(row, col);
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
			{
				zaromat(j, i) = arraytype(0);
			}
		}
		return zaromat;
	}

	/*
	 * 対角行列の取得
	 */
	KMatrix<arraytype>& diag()
	{
		if (m_Rows != 1)
		{
			// 1次元のみ対応
		}

		KMatrix<arraytype> diagmat = zeros(m_Cols, m_Cols);
		for (int i = 0; i < m_Cols; ++i)
		{
			diagmat(i, i) = operator()(0, i);
		}
		return diagmat;
	}

	/*
	 * ()演算子:行列アクセス
	 */
	arraytype& operator()(int row, int col)
	{
		return this->m_data[row * m_Cols + col];
	}

	/*
	 * *演算子:行列積
	 */
	KMatrix<arraytype>& operator*(const KMatrix<arraytype>& b)
	{
		if (m_Rows != b.m_Cols)
		{
			// 計算不可
		}

		KMatrix<arraytype> matmul(m_Rows, b.m_Cols);
		for (int k = 0; k < m_Rows; ++k)
		{
			for (int i = 0; i < b.m_Cols; ++i)
			{
				for (int j = 0; j < b.m_Rows; ++j)
				{
					matmul(k, i) += operator()(k, j) * b(j, i);
				}
			}
		}

		return matmul;
	}

	/*
	 * =演算子:代入
	 */
	KMatrix<arraytype>& operator=(const KMatrix<arraytype>& b)
	{
		m_Rows = b.m_Rows;
		m_Cols = b.m_Cols;
		m_data = b.m_data;
	}

};