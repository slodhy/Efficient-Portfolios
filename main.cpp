#include "matrix.h"
#include <iostream>

int main(int argc, char** argv) {

    QSMatrix<double> one_by_one(1, 1, 0.0); // dummy variable to access the scalar result from 1 x 1 output from matrix multiplication

    QSMatrix<double> cov_matrix(4, 4, 0.0);

    cov_matrix(0, 0) = 0.09; // Example from page 263, A Linear Algebra Primer for Financial Engineering, Dan Stefanica 2014
    cov_matrix(0, 1) = 0.01;
    cov_matrix(0, 2) = 0.03;
    cov_matrix(0, 3) = -0.015;
    cov_matrix(1, 0) = 0.01;
    cov_matrix(1, 1) = 0.0625;
    cov_matrix(1, 2) = -0.02;
    cov_matrix(1, 3) = -0.01;
    cov_matrix(2, 0) = 0.03;
    cov_matrix(2, 1) = -0.02;
    cov_matrix(2, 2) = 0.1225;
    cov_matrix(2, 3) = 0.02;
    cov_matrix(3, 0) = -0.015;
    cov_matrix(3, 1) = -0.01;
    cov_matrix(3, 2) = 0.02;
    cov_matrix(3, 3) = 0.0576;

    QSMatrix<double> mu(4, 1, 0.0); // There is probably a better way to code this than assuming the 4 x 1 matrix as a column vector

    mu(0, 0) = 0.02;
    mu(1, 0) = 0.0175;
    mu(2, 0) = 0.025;
    mu(3, 0) = 0.015;

    double risk_free = 0.01;
    double expected_return = 0.0225;
    QSMatrix<double> mu_bar = mu - risk_free;
    
    QSMatrix<double> min_var_weight = min_var_portfolio(cov_matrix, mu_bar, risk_free, expected_return); 
    /*  Pseudocode found on page 255 shows a return for all the results of the weight vector, cash weight scalar and min standard deviation scalar, 
    but I had trouble creating a function that would return mutliple values (using tuple or auto) including a Matrix class variable. As a result,
    having the function return the weight vector and calculating the cash weight scalar and min standard deviation scalar outside was my solution*/

    QSMatrix<double> one_t(1, mu.get_rows(), 1.0); // transposed vector of 1s with the length of mu

    one_by_one = multiply(one_t, min_var_weight) * (-1) + 1;
    double min_cash_weight = one_by_one(0, 0);
    one_by_one = multiply(transpose(min_var_weight), multiply(cov_matrix, min_var_weight));
    double min_standard_dev = pow(one_by_one(0, 0),0.5);

    for (int i = 0; i < min_var_weight.get_rows(); i++) {
        for (int j = 0; j < min_var_weight.get_cols(); j++) {
            std::cout << min_var_weight(i, j) << ", ";
        }
        std::cout << std::endl;
    }   
    std::cout << std::endl;
    std::cout << min_cash_weight;
    std::cout << std::endl;
    std::cout << min_standard_dev;
    std::cout << std::endl;

    std::cout << std::endl;
    QSMatrix<double> tangency_weight = tangency_portfolio(cov_matrix, mu_bar, risk_free); 
    /* I created only one function for the tangency portfolio, rather than 2 separate ones for minimum variance and maximum return. The thought
    process behind this is since the tangency portfolio is fully invested (0 cash), the tangency portfolio for both minimum variance and maximum return 
    are the same. Once that portfolio weight vector is calculated, the minimum variance and maximum return portfolio weight vectors can be calculated 
    using the other data points*/
    
    for (int i = 0; i < tangency_weight.get_rows(); i++) {
        for (int j = 0; j < tangency_weight.get_cols(); j++) {
            std::cout << tangency_weight(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    
    one_by_one = multiply(transpose(mu_bar), tangency_weight);
    double min_cash_weight2 = 1 - ((expected_return - risk_free) / one_by_one(0, 0));
    
    QSMatrix<double> min_var_weight2 = tangency_weight * (1 - min_cash_weight2);
    one_by_one = multiply(transpose(min_var_weight2), multiply(cov_matrix, min_var_weight2));
    double min_standard_dev2 = pow(one_by_one(0, 0), 0.5);

    // Test to see of both minimum variance portfolios match
    
    std::cout << std::endl;
    for (int i = 0; i < min_var_weight2.get_rows(); i++) {
        for (int j = 0; j < min_var_weight2.get_cols(); j++) {
            std::cout << min_var_weight2(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << min_cash_weight2;
    std::cout << std::endl;
    std::cout << min_standard_dev2;
    std::cout << std::endl;

    double expected_standard_dev = 0.27;

    QSMatrix<double> max_ret_weight = max_ret_portfolio(cov_matrix, mu_bar, risk_free, expected_standard_dev);
    /* Similar issue as the minimum variance function, could only manage to return just the matrix class*/

    std::cout << std::endl;
    for (int i = 0; i < max_ret_weight.get_rows(); i++) {
        for (int j = 0; j < max_ret_weight.get_cols(); j++) {
            std::cout << max_ret_weight(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    one_by_one = multiply(one_t, max_ret_weight);
    double max_cash_weight = 1 - one_by_one(0, 0);
    one_by_one = multiply(transpose(mu_bar), max_ret_weight);
    double max_return = risk_free + one_by_one(0, 0);
        
    std::cout << std::endl;
    std::cout << max_cash_weight;
    std::cout << std::endl;
    std::cout << max_return;
    std::cout << std::endl;

    one_by_one = multiply(one_t, linear_solve_cholesky(cov_matrix, mu_bar)); 
    
    /* it would probably be more computationally efficient to call the linear solve outside of the portfolio weight return functions, 
    because then you would only need to run it once, rather than mutliple times within each function. I have left it in for now for follow the 
    pseudocode as much as possible*/
    
    double max_cash_weight2 = 0;

    if (one_by_one(0, 0) > 0) {
        one_by_one = multiply(transpose(tangency_weight), multiply(cov_matrix, tangency_weight));
        max_cash_weight2 = 1 - (expected_standard_dev / pow(one_by_one(0, 0),0.5));
    }
    else
    {
        one_by_one = multiply(transpose(tangency_weight), multiply(cov_matrix, tangency_weight));
        max_cash_weight2 = 1 + (expected_standard_dev / pow(one_by_one(0, 0),0.5));
    }
    
    QSMatrix<double> max_ret_weight2 = tangency_weight * (1 - max_cash_weight2);
    one_by_one = multiply(transpose(mu_bar), max_ret_weight2);
    double max_return2 = risk_free + one_by_one(0, 0);
    
    // Test to see of both maximum return portfolios match

    std::cout << std::endl;
    for (int i = 0; i < max_ret_weight2.get_rows(); i++) {
        for (int j = 0; j < max_ret_weight2.get_cols(); j++) {
            std::cout << max_ret_weight2(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << max_cash_weight2;
    std::cout << std::endl;
    std::cout << max_return2;
    std::cout << std::endl;

}