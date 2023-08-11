#ifndef KMATRIX_H
#define KMATRIX_H

#include <algorithm>
#include <iterator>

/*
 * 行列演算クラス
 * ToDo; 逆行列生成
 */
template<typename arraytype>
class KMatrix
{
private:

	unsigned int m_Rows;	// 行数
	unsigned int m_Cols;	// 列数
	unsigned int m_Size;	// 配列サイズ
	arraytype* m_Data;		// 配列データ

public:

	enum AXIS_KIND
	{
		ROW = 0,
		COL,
		MAX
	};

	/*
	 * コンストラクタ
	 */
	KMatrix()
		: m_Rows()
		, m_Cols()
		, m_Size()
		, m_Data()
	{
	}

	/*
	 * コンストラクタ
	 */
	KMatrix(unsigned int row, unsigned int col)
		: m_Rows(row)
		, m_Cols(col)
		, m_Size(row * col)
	{
		m_Data = new arraytype[m_Size];
		memset(m_Data, 0, sizeof(arraytype) * m_Size);
	}
	 
	/*
	 * デストラクタ
	 */
	~KMatrix()
	{
		if (m_Data != nullptr)
		{
			//delete[] m_Data;
		}
	}

	/*
	 * 行数を取得
	 */
	const unsigned int size() const
	{
		return this->m_Size;
	}

	/*
	 * 行数を取得 
	 */
	const unsigned int row() const
	{
		return this->m_Rows;
	}

	/*
	 * 列数を取得
	 */
	const unsigned int col() const
	{
		return this->m_Cols;
	}

	/*
	 * 配列データを取得
	 */
	arraytype* data() const
	{
		return this->m_Data;
	}

	/*
	 * 指定行のデータを取得
	 */
	KMatrix<arraytype> rowdata(unsigned int row) const
	{

		KMatrix<arraytype> rowdata(1, this->m_Cols);
		memcpy(rowdata.m_Data, this->m_Data + (row * this->m_Cols), sizeof(arraytype) * this->m_Cols);
		return rowdata;
	}

	/*
	 * 指定列のデータを取得
	 */
	KMatrix<arraytype> coldata(unsigned int col) const
	{
		KMatrix<arraytype> coldata(this->m_Rows, 1);
		memcpy(coldata.m_Data, this->transpose().m_Data + (col * this->m_Rows), sizeof(arraytype) * this->m_Rows);
		return coldata;
	}

	/*
	 * 配列データの最大値を取得
	 */
	arraytype amax() const
	{
		// 1次元のみ対応
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::amax error : data is not 1 demension.");
		}

		arraytype maxval = this->m_Data[0];
		for (unsigned int i = 1; i < this->m_Size; ++i)
		{
			if (maxval < this->m_Data[i])
			{
				maxval = this->m_Data[i];
			}
		}
		return maxval;
	}

	/*
	 * 配列データの最小値を取得
	 */
	arraytype amin() const
	{
		// 1次元のみ対応
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::diag error : data is not 1 demension.");
		}

		arraytype minval = this->m_Data[0];
		for (unsigned int i = 1; i < this->m_Size; ++i)
		{
			if (minval > this->m_Data[i])
			{
				minval = this->m_Data[i];
			}
		}
		return minval;
	}

	/*
	 * 配列データの最大値のindexを取得
	 */
	unsigned int argmax() const
	{
		// 1次元のみ対応
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::diag error : data is not 1 demension.");
		}
		arraytype maxval = this->m_Data[0];
		unsigned int maxelement = 0;
		for (unsigned int i = 1; i < this->m_Size; ++i)
		{
			if (maxval < this->m_Data[i])
			{
				maxelement = i;
			}
		}
		return maxelement;
	}

	/*
	 * 配列データの最小値のindexを取得
	 */
	unsigned int argmin() const
	{
		// 1次元のみ対応
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::diag error : data is not 1 demension.");
		}

		arraytype minval = this->m_Data[0];
		unsigned int minelement = 0;
		for (unsigned int i = 1; i < this->m_Size; ++i)
		{
			if (minval > this->m_Data[i])
			{
				minelement = i;
			}
		}
		return minelement;
	}

	/*
	 * 転置行列を取得
	 */
	KMatrix<arraytype> transpose() const
	{
		KMatrix<arraytype> transmat(m_Cols, m_Rows);

		// 1次元の場合
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			transmat.m_Cols = this->m_Rows;
			transmat.m_Rows = this->m_Cols;
			transmat.m_Data = this->m_Data;
		}
		else
		{
			for (unsigned int i = 0; i < this->m_Rows; ++i)
			{
				for (unsigned int j = 0; j < this->m_Cols; ++j)
				{
					transmat(j, i) = operator()(i, j);
				}
			}
		}
		return transmat;
	}

	/*
	 * ゼロベクトルを取得
	 */
	void zeros()
	{
		memset(this->m_Data, 0, sizeof(arraytype) * this->m_Size);
	}

	/*
	 * サイズの確保
	 */
	void reserve(unsigned int row, unsigned int col)
	{
		if (this->m_Data != nullptr)
		{
			delete[] this->m_Data;
		}
		this->m_Rows = row;
		this->m_Cols = col;
		this->m_Size = row * col;
		this->m_Data = new arraytype[m_Size];
		memset(m_Data, 0, sizeof(arraytype) * this->m_Size);
	}

	/*
	 * 行・列の変更
	 */
	void resize(unsigned int row, unsigned int col)
	{
		// 変更不可
		if (this->m_Size != row * col)
		{
			throw std::runtime_error("KMatrix::resize error : data size is not same arg (row * col).");
		}
		this->m_Rows = row;
		this->m_Cols = col;
	}

	/*
	 * 対角行列を取得
	 */
	KMatrix<arraytype>& diag() const
	{
		// 1次元のみ対応
 		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::diag error : data is not 1 demension.");
		}

		// 対角上に配列を詰める
		KMatrix<arraytype> diagmat(this->m_Size, this->m_Size);
		for (unsigned int i = 0; i < this->m_Size; ++i)
		{
			diagmat(i, i) = this->m_Data[i];
		}
		return diagmat;
	}

	/*
	 * 全要素1の行列を取得
	 */
	void ones()
	{
		// 対角上に配列を詰める
		for (unsigned int i = 0; i < this->m_Size; ++i)
		{
			this->m_Data[i] = arraytype(1);
		}
	}

	/*
	 * 行列の絶対値を取得
	 */
	KMatrix<arraytype>& abs() const
	{
		KMatrix<arraytype> absmat(this->m_Rows, this->m_Cols);
		memcpy(absmat.m_Data, this->m_Data, sizeof(arraytype) * this->m_Size);
		for (unsigned int i = 0; i < this->m_Size; ++i)
		{
			if (0 > absmat.m_Data[i])
			{
				absmat.m_Data[i] = std::abs(absmat.m_Data[i]);
			}
		}
		return absmat;
	}

	/*
	 * 要素毎の積を取得
	 */
	KMatrix<arraytype>& multiply(const KMatrix<arraytype>& b) const
	{
		// 計算不可
		if (this->m_Cols != b.m_Cols || this->m_Rows != b.m_Rows)
		{
			throw std::runtime_error("KMatrix::operator* error ; size is not same.");
		}

		KMatrix<arraytype> multiply(this->m_Rows, this->m_Cols);
		for (unsigned int i = 0; i < this->m_Size; ++i)
		{
			multiply.m_Data[i] = this->m_Data[i] * b.m_Data[i];
		}

		return multiply;
	}

	/*
	 * 平均値を取得
	 */
	arraytype mean() const
	{
		// 1次元のみ対応
		if (this->m_Rows != 1 && this->m_Cols != 1)
		{
			throw std::runtime_error("KMatrix::diag error : data is not 1 demension.");
		}

		arraytype meanval = 0;
		for (unsigned int i = 0; i < this->m_Size; ++i)
		{
			meanval += this->m_Data[i];
		}

		return meanval /= this->m_Size;
	}

	/*
	 * ()演算子:行列アクセス
	 */
	arraytype& operator()(unsigned int row, unsigned int col) const
	{
		return this->m_Data[row * this->m_Cols + col];
	}

	/*
	 * *演算子:行列積
	 */
	KMatrix<arraytype>& operator*(const KMatrix<arraytype>& b) const
	{
		// 計算不可
		if (this->m_Cols != b.m_Rows)
		{
			throw std::runtime_error("KMatrix::operator* error ; cols is not same row.");
		}

		// 行列積計算
		KMatrix<arraytype> matmul(this->m_Rows, b.m_Cols);
		for (unsigned int k = 0; k < this->m_Rows; ++k)
		{
			for (unsigned int i = 0; i < b.m_Cols; ++i)
			{
				for (unsigned int j = 0; j < b.m_Rows; ++j)
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
		this->m_Rows = b.m_Rows;
		this->m_Cols = b.m_Cols;
		this->m_Size = b.m_Size;
		memcpy(this->m_Data, b.m_Data, sizeof(arraytype) * this->m_Size);
	}

};
#endif	// KMATRIX_H
