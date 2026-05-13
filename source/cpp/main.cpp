#include <iostream>
#include "../../include/Interval.h"
#include <mpfr.h>
#include <vector>
#include <map>
#include "../../include/mpreal.h"
#include <fstream>
#include <cmath>

using namespace std;
using namespace interval_arithmetic;

std::ofstream dbg("debug.log", std::ios::app);

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
        if constexpr(is_same_v<T, long double>)
        {
            return 0;
        }
        else if constexpr(is_same_v<T, Interval<long double>>)
        {
            return IntRead<long double>("0");
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
        if constexpr(is_same_v<T, long double>)
        {
            return 0;
        }
        else if constexpr(is_same_v<T, Interval<long double>>)
        {
           return IntRead<long double>("0");
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

    T q, r;

    if constexpr (is_same_v<T, long double>)
    {
        r = coefficients.at(coefficients.size() - 2) / coefficients.at(coefficients.size() - 1);
        q = coefficients.at(coefficients.size() - 3) / coefficients.at(coefficients.size() - 1);
    }
    else if constexpr (is_same_v<T, Interval<long double>>)
    {
        long double r_mid = coefficients.at(coefficients.size() - 2).Mid() / coefficients.at(coefficients.size() - 1).Mid();
        long double q_mid = coefficients.at(coefficients.size() - 3).Mid() / coefficients.at(coefficients.size() - 1).Mid();

        r.a = r_mid;
        r.b = r_mid;

        q.a = q_mid;
        q.b = q_mid;
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

        bool contains_zero = false;
        if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>) {
            contains_zero = (denominator.a <= 0 && denominator.b >= 0);
        } else {
            contains_zero = (abs(denominator) < 1e-14); 
        }

        if (contains_zero) {
            if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>) {
                auto new_r = r.Mid() * 1.5 + (rand() % 100) / 100.0;
                auto new_q = q.Mid() * 0.5 - (rand() % 100) / 100.0;
                
                r.a = new_r; 
                r.b = new_r;
                q.a = new_q; 
                q.b = new_q;
            } else {
                r = r * 1.5 + (rand() % 100) / 100.0;
                q = q * 0.5 - (rand() % 100) / 100.0;
            }
            continue; 
        }

        if constexpr (is_same_v<T, Interval<long double>>) {
            denominator = IntRead<long double>("1") / denominator;    
        } else {
            denominator = 1 / denominator;
        }

        T dr = (A_val * Bq - Aq * B_val) * denominator;
        T dq = (Ar * B_val - A_val * Br) * denominator;

        r = r - dr;
        q = q - dq;

        if constexpr (is_same_v<T, Interval<long double>>)
        {
            long double mid_r = r.Mid();
            long double mid_q = q.Mid();
            
            r.a = mid_r;
            r.b = mid_r;
            q.a = mid_q;
            q.b = mid_q;
        }


        if constexpr(is_same_v<T, Interval<long double>>) {
            if (abs(dr.Mid()) < tolerance && abs(dq.Mid()) < tolerance) 
                break;
        } else {
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

    if constexpr (is_same_v<T, long double>)
    {
        root.real = -1 * coefficients.at(0) / coefficients.at(1);
        root.imag = 0;
        output.push_back(root);
    } 
    else if constexpr(is_same_v<T, Interval<long double>>)
    {
        root.imag = IntRead<long double>("0");
        root.real = IntRead<long double>("0") - (coefficients.at(0) / coefficients.at(1));               
        output.push_back(root);
    }

    return output;
}

template<typename T>
vector<complex_interval<T>> second_degree_roots(vector<T> coefficients, long double zero)
{
    vector<complex_interval<T>> output;

    if constexpr (is_same_v<T, long double>)
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
    } 
    else if constexpr(is_same_v<T, Interval<long double>>)
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
            r1.imag = IntRead<long double>("0");
            r2.real = (-1 * b - ISqrt(delta,x)) / (2 * a);
            r2.imag = IntRead<long double>("0");
        }
        else if (delta.b < 0) {
            r1.real = -1 * b / (2 * a);
            r1.imag = ISqrt((-1 * delta), x) / (2 * a);
            r2.real = r1.real;
            r2.imag = -1 * r1.imag;
        }
        else {
            T delta_pos = delta; 
            delta_pos.a = 0; 
            
            T delta_neg = delta; 
            delta_neg.b = 0; 
            delta_neg = -1 * delta_neg;

            r1.real = (-1 * b + ISqrt(delta_pos, x)) / (2 * a);
            r2.real = (-1 * b - ISqrt(delta_pos, x)) / (2 * a);

            r1.imag = ISqrt(delta_neg, x) / (2 * a);
            r2.imag = -1 * r1.imag;
        }

        output.push_back(r1);
        output.push_back(r2); 
    }
    

    return output;
}


template<typename T>
vector<complex_interval<T>> bairstow_method(int degree, vector<T> coefficients, int max_iterations, long double tolerance, long double zero)
{
    vector<complex_interval<T>> output;

    if (degree == 0)
    {
        return output;
    }

    if (degree == 1)
    {
        return first_degree_roots<T>(coefficients, zero);
    }

    if (degree == 2)
    {
        return second_degree_roots<T>(coefficients, zero);
    }

    newton_output result = newton<T>(coefficients, max_iterations, tolerance);

    T a_coef;

    if constexpr (is_same_v<T, long double>)
    {
        a_coef = 1.0;
    }
    else if constexpr (is_same_v<T, Interval<long double>>)
    {
        a_coef = IntRead<long double>("1.0");
    }

    vector<complex_interval<T>> temp = second_degree_roots<T>({-1 * result.q, -1 * result.r, a_coef}, zero);
    output.insert(output.end(), temp.begin(), temp.end());

    vector<complex_interval<T>> rest = bairstow_method<T>(static_cast<int>(result.coefficients.size() - 1), result.coefficients, max_iterations, tolerance, zero);
    output.insert(output.end(), rest.begin(), rest.end());

    

    for(auto &i : output)
    {
        if constexpr (is_same_v<T, long double>)
        {
            if(abs(i.real) < zero)
                i.real = 0;
            if(abs(i.imag) < zero)
                i.imag = 0;
        }
        else if constexpr (is_same_v<T, Interval<long double>>)
        {
            if(abs(i.real.Mid()) < zero && abs(i.real.GetWidth()) < zero * 2)
                i.real = IntRead<long double>("0.0");
            if(abs(i.imag.Mid()) < zerodet && abs(i.imag.GetWidth()) < zero * 2)
                i.imag = IntRead<long double>("0.0");

        }
    }

    return output;
}

enum ArithmeticType
{
    RealNum = 0,
    IntervalForRealNum = 1,
    IntervalForInterval = 2
};

template<typename T>
vector<T> parse_polynom(const string polynom_str, int arithmetic)
{
    std::map<int,T> temp;

    stringstream ss(polynom_str);
    string token;

    int max_degree = 0;

    if constexpr (is_same_v<T,long double>) 
    {
        while (ss >> token)
        {
            size_t x_pos = token.find(',');
            string coeff_str = token.substr(0, x_pos);

            int degree = stoi(token.substr(x_pos + 1));
            Interval<long double> t = IntRead<long double>(coeff_str);


            temp[degree] = t.Mid();
            max_degree = max(max_degree, degree);
        }
        
    }
    else if constexpr (is_same_v<T,Interval<long double>>)
    {
        if(arithmetic == IntervalForRealNum)
        {
            while (ss >> token)
            {
                size_t x_pos = token.find(',');
                string coeff_str = token.substr(0, x_pos);

                int degree = stoi(token.substr(x_pos + 1));
                T coeff = IntRead<long double>(coeff_str);

                temp[degree] = coeff;
                max_degree = max(max_degree, degree);
            }
        }
        else
        {
            while (ss >> token)
            {
                size_t x_pos = token.rfind(',');
                string coeff_str = token.substr(0, x_pos);

                coeff_str = coeff_str.substr(1, coeff_str.length() - 2);
                string right = coeff_str.substr(coeff_str.find(',') + 1);
                string left = coeff_str.substr(0, coeff_str.find(','));
                
                T coeff;

                coeff.a = interval_arithmetic::LeftRead<long double>(left);
                coeff.b = interval_arithmetic::RightRead<long double>(right);

                int degree = stoi(token.substr(x_pos + 1));
                
                temp[degree] = coeff;
                max_degree = max(max_degree, degree);
            }
        }
    }

    vector<T> result(max_degree + 1);

    for(int i = 0; i < result.size(); ++i)
    {
        if constexpr (is_same_v<T,long double>)
        {
            result.at(i) = 0;
        }
        else if constexpr (is_same_v<T,Interval<long double>>)
        {
            result.at(i) = IntRead<long double>("0.0");
        }
    }

    for(const auto [deg, data] : temp)
    {
        result.at(deg) = data;
    }

    return result;
}

string format_IntWidth(long double width, int precision)
{
    long double multiplier = std::pow(10, precision);
    long double rounded_width = roundl(width * multiplier) / multiplier;

    // 2. Reszta Twojej logiki z wykorzystaniem biblioteki
    Interval<long double> x;
    x.a = rounded_width;
    x.b = rounded_width;

    string w, temp;
    // const_cast nie jest tu potrzebny, IEndsToStrings nie jest const, 
    // ale x i tak nie jest stałą w tej funkcji
    x.IEndsToStrings(w, temp);

    // 3. Twoje wycinanie formatu EX
    // temp zawiera teraz zaokrągloną wartość dzięki krokowi nr 1
    w = temp.substr(0, precision + 2) + temp.substr(temp.rfind('E'));

    return w;
}


int main()
{
    dbg = std::ofstream("/home/adam/Programowanie/EAN/source/cpp/debug.log", std::ios::app);
    dbg << "=== NOWY TEST ===" << std::endl;

    int max_iter = 250;
    long double p1 = 1e-14;
    long double p2 = 1e-14;
    const int precision = 2;

    int arithmetic;
    string polynom;

    cout << "RUNNINGX" << endl;

    while (true)
    {
        cin >> arithmetic;
        getline(cin >> ws, polynom);
        dbg << "Otrzymany string: [" << polynom << "]" << std::endl;

        if(polynom.size() == 0)
        {
            cout << "ERROR" << endl;
            continue;
        }

        switch (arithmetic)
        {
            case RealNum: 
            {
                vector<long double> polynominal = parse_polynom<long double>(polynom,arithmetic);
                vector<complex_interval<long double>> result = bairstow_method<long double>(polynominal.size() - 1, polynominal,max_iter,p1,p2);
                
                dbg << "Wektor po sparsowaniu: ";
                for(auto v : polynominal) dbg << v << " ";
                    dbg << std::endl;

                dbg << "Output bairstowa: ";
                for(auto v : result) 
                {
                    dbg << v.real << " + " << v.imag << "i" << endl;
                }
                dbg.flush();

                int i = 1;
                string reLeft, reRight, imLeft, imRight;
                for (const auto& root : result)
                {

                    Interval<long double> real;
                    Interval<long double> imag;

                    real.a = root.real + 0;
                    real.b = root.real + 0;

                    imag.a = root.imag + 0;
                    imag.b = root.imag + 0;

                    bool minus = false;

                    if(imag.Mid() < 0)
                    {
                        imag = imag * -1;
                        minus = true;
                    }

                    const_cast<Interval<long double>&>(real).IEndsToStrings(reLeft, reRight);
                    const_cast<Interval<long double>&>(imag).IEndsToStrings(imLeft, imRight);

                    if(minus)
                    {
                        cout << "Root " << i++  <<  " : " << reLeft << " - " 
                    << imLeft << "i" << " ?";
                    }
                    else
                    {
                        cout << "Root " << i++  <<  " : " << reLeft << " + " 
                    << imLeft << "i" << " ?";
                    }

                    
                }

                cout << endl;
                break;
            }

            

            case IntervalForRealNum:
            {
                vector<Interval<long double>> polynominal = parse_polynom<Interval<long double>>(polynom,arithmetic);
                vector<complex_interval<Interval<long double>>> result = bairstow_method<Interval<long double>>(polynominal.size() - 1, polynominal,max_iter,p1,p2);
                
                int i = 1;
                string reLeft, reRight, imLeft, imRight;
                for (const auto& root : result)
                {
                    const_cast<Interval<long double>&>(root.real).IEndsToStrings(reLeft, reRight);
                    const_cast<Interval<long double>&>(root.imag).IEndsToStrings(imLeft, imRight);

                    string realWidth = format_IntWidth(IntWidth(root.real),precision);
                    string imagWidth = format_IntWidth(IntWidth(root.real),precision);


                    cout << "Root " << i++  <<  " : " << " [" << reLeft << ", " << reRight << "] (w: " << realWidth << ") + [" 
                        << imLeft << ", " << imRight << "]i (w: " << realWidth << ") ?";
                }
                
                cout << endl;
                break;
            }

            case IntervalForInterval:
            {
                vector<Interval<long double>> polynominal = parse_polynom<Interval<long double>>(polynom,arithmetic);
                vector<complex_interval<Interval<long double>>> result = bairstow_method<Interval<long double>>(polynominal.size() - 1, polynominal,max_iter,p1,p2);
                
                int i = 1;
                string reLeft, reRight, imLeft, imRight;
                for (const auto& root : result)
                {
                    string realWidth = format_IntWidth(IntWidth(root.real),precision);
                    string imagWidth = format_IntWidth(IntWidth(root.real),precision);


                    cout << "Root " << i++  <<  " : " << " [" << reLeft << ", " << reRight << "] (w: " << realWidth << ") + [" 
                        << imLeft << ", " << imRight << "]i (w: " << realWidth << ") ?";
                }

                cout << endl;

                break;
            }
        }



        
   
    }


}
