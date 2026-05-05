#include <iostream>
#include "../../include/Interval.h"
#include <mpfr.h>
#include "../../include/mpreal.h"
#include <type_traits>

using namespace std;
using namespace interval_arithmetic;

template<typename T>
struct complex_interval
{
    T real;
    T imag;
};

template<typename T>
struct newton_output
{
    T r;
    T q;

    vector<T> coefficients;
};

//funkcja obliczająca wartości B, najprawdopodobniej nie ma w niej błędów
template<typename T>
void b(T r, T q, vector<T> &b_values, vector<T> &a_values)
{
    b_values.at(0) = a_values.at(a_values.size() - 1);
    b_values.at(1) = a_values.at(a_values.size() - 2) + b_values.at(0) * r;

    for(int i = 2; i < (int)a_values.size(); i++)
    {
        b_values.at(i) = b_values.at(i - 2) * q + b_values.at(i - 1) * r + a_values.at(a_values.size() - 1 - i);
    }
}

//funkcja obliczająca wartości C, najprawdopodobniej nie ma w niej błędów
template<typename T>
void c(T r, T q, vector<T> &c_values, vector<T> &b_values)
{   
    c_values.at(0) = b_values.at(0);
    c_values.at(1) = b_values.at(1) + c_values.at(0) * r;

    for(int i = 2; i < (int)c_values.size(); i++)
    {
        c_values.at(i) = c_values.at(i - 2) * q + c_values.at(i - 1) * r + b_values.at(i);
    }
}

//funkcja obliczająca wartości A, najprawdopodobniej nie ma w niej błędów
template<typename T>
T A(T r, T q, vector<T> &b_values, vector<T> &a_values, int degree)
{
    return b_values.at(degree - 3) * q + b_values.at(degree - 2) * r + a_values.at(1);
}

//funkcja obliczająca wartości B, najprawdopodobniej nie ma w niej błędów
template<typename T>
T B(T r, T q, vector<T> &b_values, vector<T> &a_values, int degree)
{
    return b_values.at(degree - 2) * q + a_values.at(0);
}

//funkcja obliczająca wartości A', najprawdopodobniej nie ma w niej błędów
template<typename T>
T A_(T r, T q, vector<T> &b_values, vector<T> &c_values, int degree)
{
    if(degree < 3)
    {
        if constexpr(is_same_v<T, mpreal> || is_same_v<T, long double>)
        {
            return 0;
        }
        else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
        {
            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                return IntRead<mpreal>("0");
            }
            else if constexpr(is_same_v<T, Interval<long double>>)
            {
                return IntRead<long double>("0");
            }
        }
    }
    if (degree == 3)
        return b_values.at(0); 
    if (degree == 4)
        return c_values.at(0) * r + b_values.at(1); 
    return c_values.at(degree - 5) * q + c_values.at(degree - 4) * r + b_values.at(degree - 3);
  
}


//funkcja obliczająca wartości B', najprawdopodobniej nie ma w niej błędów
template<typename T>
T B_(T r, T q, vector<T> &b_values, vector<T> &c_values, int degree)
{
    if(degree < 3)
    {
        if constexpr(is_same_v<T, mpreal> || is_same_v<T, long double>)
        {
            return 0;
        }
        else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
        {
            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                return IntRead<mpreal>("0");
            }
            else if constexpr(is_same_v<T, Interval<long double>>)
            {
                return IntRead<long double>("0");
            }
        }
    }
    if (degree == 3)
        return b_values.at(1);
    return c_values.at(degree - 4) * q + b_values.at(degree - 2);
}


//funkcja implementująca metodę Newtona dla dwóch zmiennych, trzeba się jej przyjrzeć
//może być problem związany z warunkiem zatrzymania
template<typename T>
newton_output<T> newton(vector<T> coefficients, int max_iterations, long double tolerance)
{
    newton_output<T> output;

    int degree = (int)coefficients.size() - 1;

    T r, q;

    if constexpr (is_same_v<T, mpreal> || is_same_v<T, long double>)
    {
        r = 0.5;
        q = 0.5;
    }
    else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
    {
        if constexpr(is_same_v<T, Interval<mpreal>>)
        {
            r = IntRead<mpreal>("0.5");
            q = IntRead<mpreal>("0.5");
        }
        else
        {
            r = IntRead<long double>("0.5");
            q = IntRead<long double>("0.5");
        }
    }

    vector<T> b_values(coefficients.size());
    vector<T> c_values(coefficients.size());

    for(int iter = 0; iter < max_iterations; iter++)
    {
        b(r, q, b_values, coefficients);
        c(r, q, c_values, b_values);

        T A_val = A(r, q, b_values, coefficients, degree);
        T B_val = B(r, q, b_values, coefficients, degree);
        T A_prime = A_(r, q, b_values, c_values, degree);
        T B_prime = B_(r, q, b_values, c_values, degree);
        
        T Ar = r*A_prime + B_prime;
        T Aq = A_prime;
        T Br = q*A_prime;
        T Bq = B_prime;

        T denominator = Ar * Bq - Aq * Br;

        bool contains_zero;
        if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>) {
            contains_zero = (denominator.a <= 0 && denominator.b >= 0);
        } else {
            contains_zero = (denominator == 0);
        }

        if(contains_zero)
        {
            if constexpr(is_same_v<T, mpreal> || is_same_v<T, long double>)
            {
                r += 0.000000001;
                q -= 0.000000001;
            }
            else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
            {
                if constexpr(is_same_v<T, Interval<mpreal>>)
                {
                    r = r + IntRead<mpreal>("0.000000001");
                    q = q - IntRead<mpreal>("0.000000001");
                }
                else
                {
                    r = r + IntRead<long double>("0.000000001");
                    q = q - IntRead<long double>("0.000000001");
                }
            }
            continue;
        }

        if constexpr (is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
        {
            if constexpr (is_same_v<T, Interval<mpreal>>)
            {
                denominator = IntRead<mpreal>("1") / denominator;
            }
            else
            {
                denominator = IntRead<long double>("1") / denominator;
            }
        }
        else
        {
            denominator = 1 / denominator;
        }

        T dr = (A_val * Bq - Aq * B_val) * denominator;
        T dq = (Ar * B_val - A_val * Br) * denominator;

        r = r - dr;
        q = q - dq;

        if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
        {
            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                if (abs(dr.a) < tolerance && abs(dr.b) < tolerance && abs(dq.a) < tolerance && abs(dq.b) < tolerance) 
                    break;
            }
            else if constexpr(is_same_v<T, Interval<long double>>)
            {
                if (abs(dr.a) < tolerance && abs(dr.b) < tolerance && abs(dq.a) < tolerance && abs(dq.b) < tolerance) 
                    break;
                
            }
        }
        else
        {
            if (abs(dr) < tolerance && abs(dq) < tolerance) 
                break;
           
        }
    }

    output.r = r;
    output.q = q;
    output.coefficients = vector<T>(b_values.begin(),b_values.end() - 2);

    reverse(output.coefficients.begin(), output.coefficients.end());

    return output;

}

template<typename T>
vector<complex_interval<T>> first_degree_roots(vector<T> coefficients, long double zero)
{
    vector<complex_interval<T>> output;

    complex_interval<T> root;

    if constexpr (is_same_v<T, mpreal> || is_same_v<T, long double>)
    {
        root.real = -1 * coefficients.at(0) / coefficients.at(1);
        root.imag = 0;
        output.push_back(root);

        for (auto& root : output)
        {
            if (abs(root.real) < zero) root.real = 0;
            if (abs(root.imag) < zero) root.imag = 0;
        }
    } 
    else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
    {
        
        root.real = -1 * coefficients.at(0) / coefficients.at(1);

        if constexpr(is_same_v<T, Interval<mpreal>>)
        {
            root.imag = IntRead<mpreal>("0");
            root.real = IntRead<mpreal>("0") - (coefficients.at(0) / coefficients.at(1));
        }
        else
        {
            root.imag = IntRead<long double>("0");
            root.real = IntRead<long double>("0") - (coefficients.at(0) / coefficients.at(1));
        }
            
        output.push_back(root);

        for (auto& root : output)
        {
            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                if (abs(root.real.a) < zero && abs(root.real.b) < zero) root.real = IntRead<mpreal>("0");
                if (abs(root.imag.a) < zero && abs(root.imag.b) < zero) root.imag = IntRead<mpreal>("0");
            }
            else
            {
                if (abs(root.real.a) < zero && abs(root.real.b) < zero) root.real = IntRead<long double>("0");
                if (abs(root.imag.a) < zero && abs(root.imag.b) < zero) root.imag = IntRead<long double>("0");
            }
        }
    }

    return output;
}

template<typename T>
vector<complex_interval<T>> second_degree_roots(vector<T> coefficients, long double zero)
{
    vector<complex_interval<T>> output;

    if constexpr (is_same_v<T, mpreal> || is_same_v<T, long double>)
    {
        T delta = coefficients.at(1) * coefficients.at(1) - 4 * coefficients.at(2) * coefficients.at(0);
        complex_interval<T> r1, r2;
        if (delta >= 0)
        {
            r1.real = (-coefficients.at(1) + std::sqrt(delta)) / (2 * coefficients.at(2));
            r1.imag = 0;
            r2.real = (-coefficients.at(1) - std::sqrt(delta)) / (2 * coefficients.at(2));
            r2.imag = 0;
        }
        else
        {
            r1.real = -coefficients.at(1) / (2 * coefficients.at(2));
            r1.imag =  std::sqrt(-delta)       / (2 * coefficients.at(2));
            r2.real = r1.real;
            r2.imag = -r1.imag;
        }
        output.push_back(r1);
        output.push_back(r2);

        for (auto& root : output)
        {
            if (abs(root.real) < zero) root.real = 0;
            if (abs(root.imag) < zero) root.imag = 0;
        }
    } 
    else if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>)
    {
        T a = coefficients.at(2);
        T b = coefficients.at(1);
        T c = coefficients.at(0);

        int x = 0;

        T delta = ISqr(b,x) - 4 * a * c;

        complex_interval<T> r1, r2;

        if (delta.a >= 0 && delta.b >= 0)
        {
            r1.real = (-1 * b + ISqrt(delta,x)) / (2 * a);

            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                r1.imag = IntRead<mpreal>("0");
                r2.imag = IntRead<mpreal>("0");
            }
            else if constexpr(is_same_v<T, Interval<long double>>)
            {
                r1.imag = IntRead<long double>("0");
                r2.imag = IntRead<long double>("0");
            }
            
            r2.real = (-1 * b - ISqrt(delta,x)) / (2 * a);
            
        }
        else
        {
            r1.real = -1 * b / (2 * a);
            r1.imag = ISqrt((-1 * delta),x) / (2 * a);
            r2.real = r1.real;
            r2.imag = -1 * r1.imag;
        }

        output.push_back(r1);
        output.push_back(r2);

        for (auto& root : output)
        {
            if constexpr(is_same_v<T, Interval<mpreal>>)
            {
                if (abs(root.real.a) < zero && abs(root.real.b) < zero) root.real = IntRead<mpreal>("0");
                if (abs(root.imag.a) < zero && abs(root.imag.b) < zero) root.imag = IntRead<mpreal>("0");
            }
            else if constexpr(is_same_v<T, Interval<long double>>)
            {
                if (abs(root.real.a) < zero && abs(root.real.b) < zero) root.real = IntRead<long double>("0");
                if (abs(root.imag.a) < zero && abs(root.imag.b) < zero) root.imag = IntRead<long double>("0");
            }
        }   
    }
    

    return output;
}


template<typename T>
vector<complex_interval<T>> bairstow_method(int degree, vector<T> coefficients, int max_iterations, long double tolerance, long double zerodet)
{
    vector<complex_interval<T>> output;

    if (degree == 0)
    {
        return output;
    }

    if (degree == 1)
    {
        return first_degree_roots<T>(coefficients, zerodet);
    }

    if (degree == 2)
    {
        return second_degree_roots<T>(coefficients, zerodet);
    }

    newton_output result = newton<T>(coefficients, max_iterations, tolerance);

    T a_coef;

    if constexpr (is_same_v<T, mpreal> || is_same_v<T, long double>)
    {
        a_coef = 1.0;
    }
    else if constexpr (is_same_v<T, Interval<mpreal>>)
    {
        a_coef = IntRead<mpreal>("1.0");
    }
    else if constexpr (is_same_v<T, Interval<long double>>)
    {
        a_coef = IntRead<long double>("1.0");
    }

    vector<complex_interval<T>> temp = second_degree_roots<T>({-1 * result.q, -1 * result.r, a_coef}, zerodet);
    output.insert(output.end(), temp.begin(), temp.end());

    vector<complex_interval<T>> rest = bairstow_method<T>(static_cast<int>(result.coefficients.size() - 1), result.coefficients, max_iterations, tolerance, zerodet);
    output.insert(output.end(), rest.begin(), rest.end());

    return output;
}
    

int main()
{

    Interval<mpreal>::Initialize();
    Interval<mpreal>::SetMode(PINT_MODE);

    /*
    vector<vector<long double>> tests = {{1,0,1}};

    int i = 1;

    for(const auto& coefficients : tests)
    {
        cout << "Test " << i++ << " :" << endl;

        vector<complex_interval<long double>> roots = bairstow_method<long double>(coefficients.size() - 1, coefficients, 500, 1e-20, 1e-10);

        for (const auto& root : roots)
        {
            //cout << "Root: [" << root.real.a << "," << root.real.b << "] + [" << root.imag.a << "," << root.imag.b << "]i" << endl;
            cout << "Root: " << root.real << " + " << root.imag << "i" << endl;
        }
        cout << endl;
    }*/
    

    
    vector<vector<Interval<long double>>> tests = {
    // 1. x^2 + 1 = 0 (Pierwiastki: 0 +/- 1i)
    {IntRead<long double>("1.0"), IntRead<long double>("0.0"), IntRead<long double>("1.0")},

    // 2. x^2 - 5x + 6 = 0 (Pierwiastki rzeczywiste: 2.0 i 3.0)
    {IntRead<long double>("6.0"), IntRead<long double>("-5.0"), IntRead<long double>("1.0")},

    // 3. x^4 - 1 = 0 (Pierwiastki: 1, -1, i, -i)
    {IntRead<long double>("-1.0"), IntRead<long double>("0.0"), IntRead<long double>("0.0"), IntRead<long double>("0.0"), IntRead<long double>("1.0")},

    // 4. x^3 - 2x^2 + x - 2 = 0 (Pierwiastki: 2.0, i, -i)
    {IntRead<long double>("-2.0"), IntRead<long double>("1.0"), IntRead<long double>("-2.0"), IntRead<long double>("1.0")},

    // 5. Wielomian z pierwiastkami wielokrotnymi: (x-1)^2 = x^2 - 2x + 1
    {IntRead<long double>("1.0"), IntRead<long double>("-2.0"), IntRead<long double>("1.0")},

    // 6. x^2 + 2x + 5 = 0 (Pierwiastki: -1 +/- 2i)
    {IntRead<long double>("5.0"), IntRead<long double>("2.0"), IntRead<long double>("1.0")},

    // x^4 + 5x^3 + 10x^2 + 10x + 4 = 0 (Pierwiastki: -1 +/- i, -2 +/- i)
    {IntRead<long double>("4.0"), IntRead<long double>("10.0"), IntRead<long double>("10.0"), IntRead<long double>("5.0"), IntRead<long double>("1.0")}
};
    int i = 1;

    for(const auto& coefficients : tests)
    {
        cout << "Test " << i++ << " :" << endl;

        vector<complex_interval<Interval<long double>>> roots = bairstow_method<Interval<long double>>(coefficients.size() - 1, coefficients, 16, 1e-16,1e-16);

        string reLeft, reRight, imLeft, imRight;
        for (const auto& root : roots)
        {
            // Konwersja części rzeczywistej
            const_cast<Interval<long double>&>(root.real).IEndsToStrings(reLeft, reRight);
            
            // Konwersja części urojonej
            const_cast<Interval<long double>&>(root.imag).IEndsToStrings(imLeft, imRight);

            cout << "Root: [" << reLeft << ", " << reRight << "] + [" 
                << imLeft << ", " << imRight << "]i" << endl;
        }
        cout << endl;
    }


}


