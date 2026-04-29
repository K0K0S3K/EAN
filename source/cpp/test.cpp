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

    /*
    cout << "degree: " << degree << endl;
    cout << "coefficients.size(): " << coefficients.size() << endl;
    */

    mpreal r = 0.5;
    mpreal q = 0.5;

    vector<mpreal> b_values(coefficients.size(), 0);
    vector<mpreal> c_values(coefficients.size(), 0);

    for(int iter = 0; iter < max_iterations; iter++)
    {
        b(r, q, b_values, coefficients);
        c(r, q, c_values, b_values);

        /*
        cout << "iter: " << iter << endl;
        cout << "degree-3: " << degree-3 << ", degree-2: " << degree-2 << endl;
        cout << "degree-5: " << degree-5 << ", degree-4: " << degree-4 << endl;
        */

        mpreal A_val = A(r, q, b_values, coefficients, degree);
        mpreal B_val = B(r, q, b_values, coefficients, degree);
        mpreal A_prime = A_(r, q, b_values, c_values, degree);
        mpreal B_prime = B_(r, q, b_values, c_values, degree);
        
        mpreal Ar = r*A_prime + B_prime;
        mpreal Aq = A_prime;
        mpreal Br = q*A_prime;
        mpreal Bq = B_prime;

        /*
        cout << "b: ";
        for(auto& v : b_values) cout << v << " ";
        cout << endl;

        cout << "c: ";
        for(auto& v : c_values) cout << v << " ";
        cout << endl;

        cout << "A_val=" << A_val << " B_val=" << B_val << endl;
        cout << "A_prime=" << A_prime << " B_prime=" << B_prime << endl;
        cout << "det=" << (Ar*Bq - Aq*Br) << endl;
        */

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


vector<complex_interval> bairstow_method(int degree, vector<mpreal> coefficients, int max_iterations, mpreal tolerance, mpreal zerodet)
{
    vector<complex_interval> output;

    if (degree == 0)
    {
        return output;
    }

    if (degree == 1)
    {
        complex_interval root;
        root.real = -coefficients.at(0) / coefficients.at(1);
        root.imag = 0;
        output.push_back(root);
        return output;
    }

    if (degree == 2)
    {
        mpreal delta = coefficients.at(1) * coefficients.at(1)
                     - 4 * coefficients.at(2) * coefficients.at(0);
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
        return output;
    }

    newton_output result = newton(coefficients, max_iterations, tolerance);

    mpreal delta = result.r * result.r + 4 * result.q;
    complex_interval r1, r2;
    if (delta >= 0)
    {
        r1.real = (result.r + sqrt(delta)) / 2;
        r1.imag = 0;
        r2.real = (result.r - sqrt(delta)) / 2;
        r2.imag = 0;
    }
    else
    {
        r1.real =  result.r / 2;
        r1.imag =  sqrt(-delta) / 2;
        r2.real =  result.r / 2;
        r2.imag = -sqrt(-delta) / 2;
    }
    output.push_back(r1);
    output.push_back(r2);

    int new_degree = (int)result.coefficients.size() - 1;
    vector<complex_interval> rest = bairstow_method(new_degree, result.coefficients, max_iterations, tolerance, zerodet);
    output.insert(output.end(), rest.begin(), rest.end());

    for (auto& root : output)
    {
        if (abs(root.real) < zerodet) root.real = 0;
        if (abs(root.imag) < zerodet) root.imag = 0;
    }

    return output;
}
    

int main()
{

    vector<vector<mpreal>> tests = {{-6,11,-6,1},{24,-50,35,-10,1},{4,0,5,0,1},{8,-12,14,15,7,-3,1},{-8, 12, -6, 1},{2},{4, 10, 10, 5, 1}};

    int i = 1;

    for(const auto& coefficients : tests)
    {
        cout << "Test " << i++ << " :" << endl;
        vector<complex_interval> roots = bairstow_method(coefficients.size() - 1, coefficients, 500, 1e-20, 1e-20);
        for (const auto& root : roots)
        {
            cout << "Root: " << root.real << " + " << root.imag << "i" << endl;
        }
        cout << endl;
    }

    

}


