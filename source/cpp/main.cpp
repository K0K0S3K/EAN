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
newton_output<T> newton(vector<T> coefficients, int max_iterations, long double zerodet, long double relative_error, int &st, int &it, int input_type)
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
        if(input_type == 1)
        {
            long double r_mid = coefficients.at(coefficients.size() - 2).Mid() / coefficients.at(coefficients.size() - 1).Mid();
            long double q_mid = coefficients.at(coefficients.size() - 3).Mid() / coefficients.at(coefficients.size() - 1).Mid();

            r.a = r_mid;
            r.b = r_mid;
            q.a = q_mid;
            q.b = q_mid;
        } 
        else
        {
            r = coefficients.at(coefficients.size() - 2) / coefficients.at(coefficients.size() - 1);
            q = coefficients.at(coefficients.size() - 3) / coefficients.at(coefficients.size() - 1);
        }


        
    }

    vector<T> b_values(coefficients.size());
    vector<T> c_values(coefficients.size());
    
    // Inicjalizacja domyślnie na błąd limitu iteracji (odpowiednik st=3)
    bool ISi3 = true; 
    int iter = 0;

    // Zmienne z Pascala do śledzenia historii wielkości kroku (odpowiedniki pq0 i pq1)
    long double prev_dr_mag = 1.0e63;
    long double prev_dq_mag = 1.0e63;

    for(iter = 0; iter < max_iterations; iter++)
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

        // Poprawiony warunek dla dzielenia przez przedział przecinający zero
        bool contains_zero = false;
        if constexpr(is_same_v<T, Interval<long double>>) {
            contains_zero = (denominator.a <= zerodet && denominator.b >= -zerodet);
        } else {
            contains_zero = (abs(denominator) < zerodet); 
        }

        if (contains_zero) {
            st = 2;
            return output;
        }

        if constexpr (is_same_v<T, Interval<long double>>) {
            denominator = IntRead<long double>("1") / denominator;    
        } else {
            denominator = 1 / denominator;
        }

        T dr = (A_val * Bq - Aq * B_val) * denominator;
        T dq = (Ar * B_val - A_val * Br) * denominator;


        // Mierzymy moduł bieżącego kroku
        long double dr_mag, dq_mag;
        
        // Zmienne dla warunków przedziałowych
        bool function_contains_zero = false;
        bool step_contains_zero = false;

        if constexpr(is_same_v<T, Interval<long double>>) {
            // Bezpieczny moduł - bierzemy maksymalne wychylenie przedziału (najgorszy przypadek)
            dr_mag = std::max(std::abs(dr.a), std::abs(dr.b));
            dq_mag = std::max(std::abs(dq.a), std::abs(dq.b));
            
            // Sprawdzamy czy wartości funkcji A(r,q) i B(r,q) przecięły zero
            function_contains_zero = (A_val.a <= relative_error && A_val.b >= relative_error) && 
                                     (B_val.a <= relative_error && B_val.b >= relative_error);
                                     
            // Sprawdzamy czy wektor kroku utonął w szumie numerycznym
            step_contains_zero = (dr.a <= relative_error && dr.b >= relative_error) && (dq.a <= relative_error && dq.b >= relative_error);
        } else {
            dr_mag = std::abs(dr);
            dq_mag = std::abs(dq);
        }

        // --- WARUNKI ZATRZYMANIA ---
        bool stop_condition_met = false;
        
        if constexpr(is_same_v<T, Interval<long double>>) {
            // Przerywamy, gdy:
            // 1. Zeszliśmy poniżej zakładanego błędu (relative_error)
            // 2. LUB trafiliśmy w idealne zero numeryczne (function_contains_zero)
            // 3. LUB krok ma mniejszą wartość niż szum obliczeniowy (step_contains_zero)
            stop_condition_met = (dr_mag <= relative_error && dq_mag <= relative_error) || 
                                 function_contains_zero || 
                                 step_contains_zero;
        } else {
            // Standardowy warunek dla double (zatrzymaj jeśli błąd jest mały I przestał maleć)
            stop_condition_met = (dr_mag <= relative_error && dq_mag <= relative_error) && 
                                 (dr_mag >= prev_dr_mag || dq_mag >= prev_dq_mag);
        }

        if (!stop_condition_met)
        {
            // Aplikacja kroku Newtona
            r = r - dr;
            q = q - dq;
            
            // Zapisanie bieżącego kroku do historii
            prev_dr_mag = dr_mag;
            prev_dq_mag = dq_mag;

            // Spłaszczenie przedziałów
            if constexpr (is_same_v<T, Interval<long double>>)
            {
                long double mid_r = r.Mid();
                long double mid_q = q.Mid();
                
                r.a = mid_r;
                r.b = mid_r;
                q.a = mid_q;
                q.b = mid_q;
            }
        }
        else
        {
            // Pętla przerywana - osiągnęliśmy precyzję, zero maszynowe lub stagnację
            ISi3 = false;
            break;
        }

        
    }

    

    it = iter;

    if(ISi3)
    {
        st = 3;
        return output;
    }

    output.r = r;
    output.q = q;
    output.coefficients = vector<T>(b_values.begin(), b_values.end() - 2);

    reverse(output.coefficients.begin(), output.coefficients.end());

    return output;
}

template<typename T>
vector<complex_interval<T>> first_degree_roots(vector<T> coefficients)
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
vector<complex_interval<T>> second_degree_roots(vector<T> coefficients)
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
vector<complex_interval<T>> bairstow_method(int degree, vector<T> coefficients, int max_iterations, int relative_error,long double zerodet, int &st, int &it, int input_type)
{
    vector<complex_interval<T>> output;

    if (degree == 0)
    {
        return output;
    }

    if (degree == 1)
    {
        return first_degree_roots<T>(coefficients);
    }

    if (degree == 2)
    {
        return second_degree_roots<T>(coefficients);
    }

    newton_output result = newton<T>(coefficients, max_iterations, zerodet,relative_error,st,it,input_type);

    if(st != 0)
    {
        return output;
    }

    T a_coef;

    if constexpr (is_same_v<T, long double>)
    {
        a_coef = 1.0;
    }
    else if constexpr (is_same_v<T, Interval<long double>>)
    {
        a_coef = IntRead<long double>("1.0");
    }

    vector<complex_interval<T>> temp = second_degree_roots<T>({-1 * result.q, -1 * result.r, a_coef});
    output.insert(output.end(), temp.begin(), temp.end());

    vector<complex_interval<T>> rest = bairstow_method<T>(static_cast<int>(result.coefficients.size() - 1), result.coefficients, max_iterations, relative_error, zerodet,st,it,input_type);
    output.insert(output.end(), rest.begin(), rest.end());

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
    Interval<long double> x;
    x.a = width;
    x.b = width;

    string w, temp;

    x.IEndsToStrings(w, temp);

    std::string number_part = w.substr(0, w.rfind('E'));
    std::string ex_part = w.substr(w.rfind('E'));

    double float_data = std::stod(number_part);
    
    std::ostringstream stream;

    stream << std::fixed << std::setprecision(precision) << float_data;

    std::string float_formatted_data = stream.str();
    std::string final_result = float_formatted_data + ex_part;

    return final_result;
}


int main()
{
    dbg = std::ofstream("/home/adam/Programowanie/EAN/source/cpp/debug.log", std::ios::app);
    dbg << "=== NOWY TEST ===" << std::endl;

    int arithmetic;
    string polynom;

    cout << "RUNNINGX" << endl;

    while (true)
    {
        int max_iter = 5;
        long double zerodet = 1e-16;
        long double relative_error = 1e-16;

        int it = 0;
        int st = 0;
        const int precision = 1;



        cin >> arithmetic;
        getline(cin >> ws, polynom);
        cin >> max_iter;
        cin >> relative_error;
        cin >> zerodet;
        dbg << "Otrzymany string: [" << polynom << "]" << std::endl;

        switch (arithmetic)
        {
            case RealNum: 
            {
                vector<long double> polynominal = parse_polynom<long double>(polynom,arithmetic);
                if(polynominal.size() < 2 || max_iter < 1 || relative_error <= 0 || zerodet <= 0)
                {
                    cout << "st = 1" << endl;
                    continue;
                }
                vector<complex_interval<long double>> result = bairstow_method<long double>(polynominal.size() - 1, polynominal,max_iter,relative_error,zerodet,st,it,arithmetic);

                if(st != 0)
                {
                    cout << "st = " << st << endl;
                    break;
                }
                
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
                    << imLeft << "i" << " ? ";
                    }
                    else
                    {
                        cout << "Root " << i++  <<  " : " << reLeft << " + " 
                    << imLeft << "i" << " ? ";
                    }

                    
                }

                cout << "st = 0, it = " << it << "?" <<endl;

                break;
            }

            

            case IntervalForRealNum:
            {
                vector<Interval<long double>> polynominal = parse_polynom<Interval<long double>>(polynom,arithmetic);
                if(polynominal.size() < 2 || max_iter < 1 || relative_error <= 0 || zerodet <= 0)
                {
                    cout << "st = 1" << endl;
                    continue;
                }
                vector<complex_interval<Interval<long double>>> result = bairstow_method<Interval<long double>>(polynominal.size() - 1, polynominal,max_iter,relative_error,zerodet,st,it,arithmetic);
                
                if(st != 0)
                {
                    cout << "st = " << st << endl;
                    break;
                }

                int i = 1;
                string reLeft, reRight, imLeft, imRight;
                for (const auto& root : result)
                {
                    const_cast<Interval<long double>&>(root.real).IEndsToStrings(reLeft, reRight);
                    const_cast<Interval<long double>&>(root.imag).IEndsToStrings(imLeft, imRight);

                    string realWidth = format_IntWidth(IntWidth(root.real),precision);
                    string imagWidth = format_IntWidth(IntWidth(root.imag),precision);                    

                    cout << "Root " << i++  <<  " : " << " [" << reLeft << ", " << reRight << "] (w: " << realWidth << ") + [" 
                        << imLeft << ", " << imRight << "]i (w: " << imagWidth << ") ?";
                }

                cout << "st = 0, it = " << it << "?" <<endl;
                break;
            }

            case IntervalForInterval:
            {
                vector<Interval<long double>> polynominal = parse_polynom<Interval<long double>>(polynom,arithmetic);
                if(polynominal.size() < 2 || max_iter < 1 || relative_error <= 0 || zerodet <= 0)
                {
                    cout << "st = 1" << endl;
                    continue;
                }
                vector<complex_interval<Interval<long double>>> result = bairstow_method<Interval<long double>>(polynominal.size() - 1, polynominal,max_iter,relative_error,zerodet,st,it,arithmetic);
                
                if(st != 0)
                {
                    cout << "st = " << st << endl;
                    break;
                }

                int i = 1;
                string reLeft, reRight, imLeft, imRight;
                for (const auto& root : result)
                {

                    const_cast<Interval<long double>&>(root.real).IEndsToStrings(reLeft, reRight);
                    const_cast<Interval<long double>&>(root.imag).IEndsToStrings(imLeft, imRight);

                    string realWidth = format_IntWidth(IntWidth(root.real),precision);
                    string imagWidth = format_IntWidth(IntWidth(root.imag),precision);


                    cout << "Root " << i++  <<  " : " << " [" << reLeft << ", " << reRight << "] (w: " << realWidth << ") + [" 
                        << imLeft << ", " << imRight << "]i (w: " << imagWidth << ") ?";
                }

                cout << "st = 0, it = " << it <<endl;

                break;
            }
        }



        
   
    }


}
