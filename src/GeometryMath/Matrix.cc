#include <GeometryMath/Matrix.h>

using namespace FiniteElement;

template<class T, int rows, int cols>
Matrix<T, rows, cols>::Matrix()
    : d_num_rows(rows), d_num_columns(cols)
{
}

template<class T, int rows, int cols>
Matrix<T, rows, cols>::Matrix(const T& initialValue)
 : d_num_rows(rows), d_num_columns(cols)
{
   d_data.resize(d_num_rows*d_num_columns, initialValue);
}


template<class T, int rows, int cols>
const T& Matrix<T, rows, cols>::get(int row, int column) const
{
  return d_data[row][column];
}


template<class T, int rows, int cols>
T& Matrix<T, rows, cols>::get(int row, int column)
{
  return d_data[row][column];
}

template<class T, int rows, int cols>
void Matrix<T, rows, cols>::set(int row, int column, const T& value)
{
  d_data[row][column] = value;
}


template<class T, int rows, int cols>
Matrix<T, rows, cols> Matrix<T, rows, cols>::inverse(const Matrix<T, rows, cols>& mat)
{
  return inv(mat);
}












