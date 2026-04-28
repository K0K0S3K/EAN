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

    for(int i = 2; i < a_values.size(); i++)
    {
        b_values.at(i) = b_values.at(i - 2) * q + b_values.at(i - 1) * r + a_values.at(a_values.size() - 1 - i);
    }
}

void c(mpreal r, mpreal q, vector<mpreal> &c_values, vector<mpreal> &b_values)
{   
    c_values.at(0) = b_values.at(0);
    c_values.at(1) = b_values.at(1) + c_values.at(0) * r;

    for(int i = 2; i < c_values.size(); i++)
    {
        c_values.at(i) = c_values.at(i - 2) * q + c_values.at(i - 1) * r + b_values.at(i);
    }
}


mpreal A(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &a_values, int degree)
{
    return b_values.at(degree - 3) * q + b_values.at(degree - 2) * r + a_values.at(degree - 2);
}

mpreal B(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &a_values, int degree)
{
    return b_values.at(degree - 2) * q + a_values.at(degree - 1);
}

mpreal A_(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &c_values, int degree)
{
    return c_values.at(degree - 5) * q + c_values.at(degree - 4) * r + b_values.at(degree - 3);
}

mpreal B_(mpreal r, mpreal q, vector<mpreal> &b_values, vector<mpreal> &c_values, int degree)
{
    return c_values.at(degree - 4) * q + b_values.at(degree - 2);
}

newton_output newton(vector<mpreal> coefficients, int max_iterations, mpreal tolerance)
{
    newton_output output;

    mpreal r = 0.5;
    mpreal q = 0.5;

    vector<mpreal> b_values(coefficients.size(), 0);
    vector<mpreal> c_values(coefficients.size(), 0);

    for(int iter = 0; iter < max_iterations; iter++)
    {
        b(r, q, b_values, coefficients);
        c(r, q, c_values, b_values);

        mpreal A_val = A(r, q, b_values, coefficients, coefficients.size());
        mpreal B_val = B(r, q, b_values, coefficients, coefficients.size());
        mpreal A_prime = A_(r, q, b_values, c_values, coefficients.size());
        mpreal B_prime = B_(r, q, b_values, c_values, coefficients.size());
        
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

    return output;

}

vector<complex_interval> bairstow_method(int degree, vector<mpreal> coefficients, int max_iterations, mpreal tolerance, mpreal zerodet)
{
    vector<complex_interval> output;    


    return output;
}
    

int main()
{

    vector<mpreal> coefficients = {4, 10, 10, 5, 1}; // a0 to wyraz wolny an to wyraz przy najwyzszej potedze

    mpreal tolerance = 1e-20;

    newton_output output = newton(coefficients, 100, tolerance);

    cout << "R: " << output.r << endl;
    cout << "Q: " << output.q << endl;
    
    cout << "B values: ";
    for(const auto& val : output.coefficients)
    {
        cout << val << " ";
    }
    cout << endl;


}


