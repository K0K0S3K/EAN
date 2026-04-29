#include <iostream>
#include "../../include/Interval.h"
#include <mpfr.h>
#include "../../include/mpreal.h"

using namespace std;
using namespace interval_arithmetic;

struct complex_interval
{
    mpreal real;
    mpreal imag;
};

struct newton_output
{
    mpreal r;
    mpreal q;

    vector<mpreal> coefficients;
};

void b(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &a_values)
{
    b_values.at(0) = a_values.at(a_values.size() - 1);
    b_values.at(1) = a_values.at(a_values.size() - 2) + b_values.at(0) * r;

    for(int i = 2; i < (int)a_values.size(); i++)
    {
        b_values.at(i) = b_values.at(i - 2) * q + b_values.at(i - 1) * r + a_values.at(a_values.size() - 1 - i);
    }
}

void c(mpreal r, mpreal q, vector<mpreal> &c_values, vector<mpreal> &b_values)
{   
    c_values.at(0) = b_values.at(0);
    c_values.at(1) = b_values.at(1) + c_values.at(0) * r;

    for(int i = 2; i < (int)c_values.size(); i++)
    {
        c_values.at(i) = c_values.at(i - 2) * q + c_values.at(i - 1) * r + b_values.at(i);
    }
}


mpreal A(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &a_values, int degree)
{
    return b_values.at(degree - 3) * q + b_values.at(degree - 2) * r + a_values.at(1);
}

mpreal B(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &a_values, int degree)
{
    return b_values.at(degree - 2) * q + a_values.at(0);
}

mpreal A_(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &c_values, int degree)
{
    if(degree < 3)
        return 0;
    if (degree == 3)
        return b_values.at(0); 
    if (degree == 4)
        return c_values.at(0) * r + b_values.at(1); 
    return c_values.at(degree - 5) * q + c_values.at(degree - 4) * r + b_values.at(degree - 3);
  
}

mpreal B_(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &c_values, int degree)
{
    if(degree < 3)
        return 0;
    if (degree == 3)
        return b_values.at(1);
    return c_values.at(degree - 4) * q + b_values.at(degree - 2);
}

newton_output newton(vector<mpreal> coefficients, int max_iterations, mpreal tolerance)
{
    newton_output output;

    int degree = (int)coefficients.size() - 1;

    mpreal r = 0.5;
    mpreal q = 0.5;

    vector<mpreal> b_values(coefficients.size(), 0);
    vector<mpreal> c_values(coefficients.size(), 0);

    for(int iter = 0; iter < max_iterations; iter++)
    {
        b(r, q, b_values, coefficients);
        c(r, q, c_values, b_values);

        mpreal A_val = A(r, q, b_values, coefficients, degree);
        mpreal B_val = B(r, q, b_values, coefficients, degree);
        mpreal A_prime = A_(r, q, b_values, c_values, degree);
        mpreal B_prime = B_(r, q, b_values, c_values, degree);
        
        mpreal Ar = r*A_prime + B_prime;
        mpreal Aq = A_prime;
        mpreal Br = q*A_prime;
        mpreal Bq = B_prime;

        mpreal denominator = 1/(Ar * Bq - Aq * Br);

        mpreal dr = (A_val * Bq - Aq * B_val) * denominator;
        mpreal dq = (Ar * B_val - A_val * Br) * denominator;

        r = r - dr;
        q = q - dq;

        if (abs(dr) < tolerance && abs(dq) < tolerance) {
            break;  
        }
    }

    output.r = r;
    output.q = q;
    output.coefficients = vector<mpreal>(b_values.begin(),b_values.end() - 2);

    reverse(output.coefficients.begin(), output.coefficients.end());

    return output;

}

vector<complex_interval> first_degree_roots(vector<mpreal> coefficients, mpreal zerodet)
{
    vector<complex_interval> output;

    complex_interval root;
    root.real = -coefficients.at(0) / coefficients.at(1);
    root.imag = 0;
    output.push_back(root);

    for (auto& root : output)
    {
        if (abs(root.real) < zerodet) root.real = 0;
        if (abs(root.imag) < zerodet) root.imag = 0;
    }

    return output;
}

vector<complex_interval> second_degree_roots(vector<mpreal> coefficients, mpreal zerodet)
{
    vector<complex_interval> output;

    mpreal delta = coefficients.at(1) * coefficients.at(1) - 4 * coefficients.at(2) * coefficients.at(0);
    complex_interval r1, r2;
    if (delta >= 0)
    {
        r1.real = (-coefficients.at(1) + sqrt(delta)) / (2 * coefficients.at(2));
        r1.imag = 0;
        r2.real = (-coefficients.at(1) - sqrt(delta)) / (2 * coefficients.at(2));
        r2.imag = 0;
    }
    else
    {
        r1.real = -coefficients.at(1) / (2 * coefficients.at(2));
        r1.imag =  sqrt(-delta)       / (2 * coefficients.at(2));
        r2.real = r1.real;
        r2.imag = -r1.imag;
    }
    output.push_back(r1);
    output.push_back(r2);

    for (auto& root : output)
    {
        if (abs(root.real) < zerodet) root.real = 0;
        if (abs(root.imag) < zerodet) root.imag = 0;
    }

    return output;
}


vector<complex_interval> bairstow_method(int degree, vector<mpreal> coefficients, int max_iterations, mpreal tolerance, mpreal zerodet)
{
    vector<complex_interval> output;

    if (degree == 0)
    {
        return output;
    }

    if (degree == 1)
    {
        return first_degree_roots(coefficients, zerodet);
    }

    if (degree == 2)
    {
        return second_degree_roots(coefficients, zerodet);
    }

    newton_output result = newton(coefficients, max_iterations, tolerance);

    vector<complex_interval> temp = second_degree_roots({result.coefficients.at(0), result.r, result.q}, zerodet);
    
    output.insert(output.end(), temp.begin(), temp.end());

    vector<complex_interval> rest = bairstow_method(static_cast<int>(result.coefficients.size() - 1), result.coefficients, max_iterations, tolerance, zerodet);
    output.insert(output.end(), rest.begin(), rest.end());

    return output;
}
    

int main()
{

    vector<vector<mpreal>> tests = {{1,0,1}};

    int i = 1;

    for(const auto& coefficients : tests)
    {
        cout << "Test " << i++ << " :" << endl;
        vector<complex_interval> roots = bairstow_method(coefficients.size() - 1, coefficients, 500, 1e-20, 1e-10);
        for (const auto& root : roots)
        {
            cout << "Root: " << root.real << " + " << root.imag << "i" << endl;
        }
        cout << endl;
    }

    

}


