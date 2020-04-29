#ifndef __QS_MATRIX_CPP
#define __QS_MATRIX_CPP

#include "matrix.h"

/* Code for the Matrix class found at https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File/, further modified to create the linear 
solver functions */

// Parameter Constructor                                                                                                                                                      
template<typename T>
QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _initial) {
    mat.resize(_rows);
    for (unsigned i = 0; i < mat.size(); i++) {
        mat[i].resize(_cols, _initial);
    }
    rows = _rows;
    cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) {
    mat = rhs.mat;
    rows = rhs.get_rows();
    cols = rhs.get_cols();
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
QSMatrix<T>::~QSMatrix() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) {
    if (&rhs == this)
        return *this;

    unsigned new_rows = rhs.get_rows();
    unsigned new_cols = rhs.get_cols();

    mat.resize(new_rows);
    for (unsigned i = 0; i < mat.size(); i++) {
        mat[i].resize(new_cols);
    }

    for (unsigned i = 0; i < new_rows; i++) {
        for (unsigned j = 0; j < new_cols; j++) {
            mat[i][j] = rhs(i, j);
        }
    }
    rows = new_rows;
    cols = new_cols;

    return *this;
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const QSMatrix<T>& rhs) {
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] + rhs(i, j);
        }
    }

    return result;
}

// Cumulative addition of this matrix and another                                                                                                                             
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            this->mat[i][j] += rhs(i, j);
        }
    }

    return *this;
}

// Subtraction of this matrix and another                                                                                                                                     
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] - rhs(i, j);
        }
    }

    return result;
}

// Cumulative subtraction of this matrix and another                                                                                                                          
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator-=(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            this->mat[i][j] -= rhs(i, j);
        }
    }

    return *this;
}

// Left multiplication of this matrix and another. This function will only multiply matrices of the same size                                                                                                                             
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            for (unsigned k = 0; k < rows; k++) {
                result(i, j) += this->mat[i][k] * rhs(k, j);
            }
        }
    }

    return result;
}

// Cumulative left multiplication of this matrix and another                                                                                                                  
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix<T>& rhs) {
    QSMatrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}

// Adjusted Matrix multiplication function that will allow matrices of different sizes to be mutliplied (as long as the mutliplcation is valid)                                                                                                                               
template<typename T>
QSMatrix<T> multiply(QSMatrix<T> a, QSMatrix<T> b) {
    unsigned rows1 = a.get_rows();
    unsigned rows2 = b.get_rows();
    unsigned cols = b.get_cols();
    QSMatrix<double> result(rows1, cols, 0.0);

    for (unsigned i = 0; i < rows1; i++) {
        for (unsigned j = 0; j < cols; j++) {
            for (unsigned k = 0; k < rows2; k++) {
                result(i, j) += a(i,k) * b(k, j);
            }
        }
    }

    return result;
}

// Calculate a transpose of this matrix. Adjusted formula to allow                                                                                                                                      
template<typename T>
QSMatrix<T> transpose(QSMatrix<T> a) {
    unsigned rows = a.get_cols();
    unsigned cols = a.get_rows();

    QSMatrix<double> result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = a(j,i);
        }
    }

    return result;
}

// Matrix/scalar addition                                                                                                                                                     
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] + rhs;
        }
    }

    return result;
}

// Matrix/scalar subtraction                                                                                                                                                  
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] - rhs;
        }
    }

    return result;
}

// Matrix/scalar multiplication                                                                                                                                               
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] * rhs;
        }
    }

    return result;
}

// Matrix/scalar division                                                                                                                                                     
template<typename T>
QSMatrix<T> QSMatrix<T>::operator/(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result(i, j) = this->mat[i][j] / rhs;
        }
    }

    return result;
}

// Multiply a matrix with a vector                                                                                                                                            
template<typename T>
std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs) {
    std::vector<T> result(rhs.size(), 0.0);

    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            result[i] = this->mat[i][j] * rhs[j];
        }
    }

    return result;
}

// Obtain a vector of the diagonal elements                                                                                                                                   
template<typename T>
std::vector<T> QSMatrix<T>::diag_vec() {
    std::vector<T> result(rows, 0.0);

    for (unsigned i = 0; i < rows; i++) {
        result[i] = this->mat[i][i];
    }

    return result;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) {
    return this->mat[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) const {
    return this->mat[row][col];
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned QSMatrix<T>::get_rows() const {
    return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned QSMatrix<T>::get_cols() const {
    return this->cols;
}

// Forward Substitution: Adapted from pseudocode table 2.1, page 40 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> forward_subst(const QSMatrix<T> l, const QSMatrix<T> b) {
    int matrix_size = l.get_rows() - 1;
    QSMatrix<double> results(b.get_rows(), b.get_cols(), 0.0);

    results(0, 0) = b(0, 0) / l(0, 0); // initial value for first x vector entry

    for (int i = 1; i < matrix_size+1; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum = sum + l(i, j) * results(j, 0);
        }
        results(i, 0) = (b(i, 0) - sum) / l(i, i);
    }
        
    return results;
}

// Backward Substitution: Adapted from pseudocode table 2.3, page 45 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> backward_subst(const QSMatrix<T> u, const QSMatrix<T> b) {
    int matrix_size = u.get_rows() - 1;
    QSMatrix<double> results(b.get_rows(), b.get_cols(), 0.0);

    results(matrix_size, 0) = b(matrix_size, 0) / u(matrix_size, matrix_size); // initial value for final x vector entry

    for (int i = matrix_size - 1; i > -1; i--) { // from matrix size n - 1 to 0
        double sum = 0;
        for (int j = i + 1; j < matrix_size + 1; j++) { // from i to matrix size 
            sum = sum + u(i, j) * results(j, 0);
        }
        results(i, 0) = (b(i, 0) - sum) / u(i, i);
    }

    return results;
}

// Cholesky: Adapted from pseudocode table 6.1, page 170 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> cholesky(QSMatrix<T> choleskyA) {
    int matrix_size = choleskyA.get_rows() - 1;
    QSMatrix<double> results(choleskyA.get_rows(), choleskyA.get_cols(), 0.0);

    for (int i = 0; i < matrix_size; i++) {
        results(i, i) = pow(choleskyA(i, i), 0.5);
        for (int k = i + 1; k < matrix_size + 1; k++) {
            results(i, k) = choleskyA(i, k) / results(i, i);
        }
        for (int j = i + 1; j < matrix_size + 1; j++) {
            for (int k = j; k < matrix_size + 1; k++) {
                choleskyA(j, k) = choleskyA(j, k) - results(i, j) * results(i, k);
            }
        }
    }
    results(matrix_size, matrix_size) = pow(choleskyA(matrix_size, matrix_size), 0.5);

    return results;
}

// Linear Solve Cholesky: Adapted from pseudocode table 6.2. page 174 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> linear_solve_cholesky(const QSMatrix<T> choleskyA, const QSMatrix<T> b) {
    
    QSMatrix<double> choleskyU = cholesky(choleskyA);
    QSMatrix<double> choleskyU_transpose = transpose(choleskyU);
    QSMatrix<double> y = forward_subst(choleskyU_transpose, b);
    QSMatrix<double> results = backward_subst(choleskyU, y);

    return results;
}

// Minmum Variance Portfolio: Adapted from pseudocode table 9.1, page 255 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> min_var_portfolio(const QSMatrix<T> cov_matrix, QSMatrix<T> mu_bar, double risk_free, double expected_return) {
    QSMatrix<double> x = linear_solve_cholesky(cov_matrix, mu_bar);
    QSMatrix<double> denominator = multiply(transpose(mu_bar), x);
    QSMatrix<double> min_weight = x * (expected_return - risk_free) / (denominator(0, 0));
    
    return min_weight;
}
// Tangency Portfolio: Adapted from pseudocode table 9.2, page 256 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> tangency_portfolio(const QSMatrix<T> cov_matrix, QSMatrix<T> mu_bar, double risk_free) {
    QSMatrix<double> x = linear_solve_cholesky(cov_matrix, mu_bar);
    QSMatrix<double> one_t(1, mu_bar.get_rows(), 1.0); // transponsed vector of 1s with the length of mu
    QSMatrix<double> denominator = multiply(one_t, x);
    QSMatrix<double> tan_weight = x * (1 / (denominator(0, 0)));

    return tan_weight;
}

// Maximum Return Portfolio: Adapted from pseudocode table 9.3, page 257 of A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
template<typename T>
QSMatrix<T> max_ret_portfolio(const QSMatrix<T> cov_matrix, QSMatrix<T> mu_bar, double risk_free, double expected_standard_dev) {
    QSMatrix<double> x = linear_solve_cholesky(cov_matrix, mu_bar);
    QSMatrix<double> denominator = multiply(transpose(mu_bar), x);
    QSMatrix<double> max_weight = x * (expected_standard_dev / pow(denominator(0, 0),0.5));
    return max_weight;
}

#endif