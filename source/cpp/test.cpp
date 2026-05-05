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
        T r = coefficients.at(coefficients.size() - 2) / coefficients.at(coefficients.size() - 1);
        T q = coefficients.at(coefficients.size() - 3) / coefficients.at(coefficients.size() - 1);
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

        
        // Sprawdzenie, czy wyznacznik zawiera zero
        bool contains_zero = false;
        if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>) {
            contains_zero = (denominator.a <= 0 && denominator.b >= 0);
        } else {
            // Zabezpieczenie przed epsilonem maszynowym dla zwykłych typów
            contains_zero = (abs(denominator) < 1e-14); 
        }

        // Logika ratunkowa
        if (contains_zero) {
            // Jakobian jest osobliwy. Silnie i asymetrycznie zaburzamy punkt startowy,
            // aby wyrzucić Newtona z martwej strefy.
            if constexpr(is_same_v<T, Interval<mpreal>> || is_same_v<T, Interval<long double>>) {
                auto new_r = r.Mid() * 1.5 + 0.1234;
                auto new_q = q.Mid() * 0.5 - 0.4321;
                
                r.a = new_r; 
                r.b = new_r;
                q.a = new_q; 
                q.b = new_q;
            } else {
                r = r * 1.5 + 0.1234;
                q = q * 0.5 - 0.4321;
            }
            // Przerywamy obecne obliczenia dr i dq, wracamy na początek pętli
            continue; 
        }

        // Dopiero gdy jesteśmy bezpieczni, wykonujemy dzielenie
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
            // Czysto ujemna delta (pełne zespolone)
            r1.real = -1 * b / (2 * a);
            r1.imag = ISqrt((-1 * delta), x) / (2 * a);
            r2.real = r1.real;
            r2.imag = -1 * r1.imag;
        }
        else {
            // delta przecina zero (np. a < 0, b >= 0). Tniemy ją na dwie osobne obwiednie.
            T delta_pos = delta; 
            delta_pos.a = 0; // Ucinamy część ujemną dla dziedziny rzeczywistej
            
            T delta_neg = delta; 
            delta_neg.b = 0; // Ucinamy część dodatnią
            delta_neg = -1 * delta_neg; // Robimy z niej dodatnią pod pierwiastek zespolony

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

    if constexpr (is_same_v<T, long double>)
    {
        a_coef = 1.0;
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

    Interval<long double>::Initialize();
    Interval<long double>::SetMode(PINT_MODE);

    vector<vector<Interval<long double>>> tests = {
    // 1. x^3 - 3x^2 + 4x - 2 = 0 
    // Pierwiastki: 1.0, 1.0 +/- i
    // Klasyczny test stopnia 3 z jednym pierwiastkiem rzeczywistym i parą zespoloną.
    // --- STOPIEŃ 3 ---

    // 1. x^3 - 2x^2 - x + 2 = 0
{IntRead<long double>("5.0"), IntRead<long double>("6.0"), IntRead<long double>("2.0")},

    // --- STOPIEŃ 3 ---

    // 2. x^3 - 4x^2 + 9x - 10 = 0
    // Pierwiastki: 2.0, 1.0 +/- 2i
    // Asymetryczny układ: jeden pierwiastek rzeczywisty i "szeroka" para zespolona. 
    // Sprawdza, czy Bairstow dobrze odseparuje czynnik liniowy na samym końcu po wycięciu czynnika kwadratowego.
    {IntRead<long double>("-10.0"), IntRead<long double>("9.0"), IntRead<long double>("-4.0"), IntRead<long double>("1.0")},

    // 3. x^3 + 0.7x^2 - 5.58x + 3.96 = 0
    // Pierwiastki: 1.1, 1.2, -3.0
    // Bardzo dobry test na stabilność. Dwa pierwiastki rzeczywiste są blisko siebie (1.1 i 1.2). 
    // Standardowa metoda Newtona często w takich miejscach oscyluje i ma problem z domknięciem przedziału błędu.
    {IntRead<long double>("3.96"), IntRead<long double>("-5.58"), IntRead<long double>("0.7"), IntRead<long double>("1.0")},

    // --- STOPIEŃ 4 ---

    // 4. x^4 + 2x^3 + 3x^2 + 2x + 2 = 0
    // Pierwiastki: +/- i, -1.0 +/- i
    // Mieszanka czysto urojonej pary (+/- i) z parą mającą przesunięcie rzeczywiste. 
    // Algorytm musi poprawnie namierzyć dwa różne czynniki kwadratowe bez potknięcia się na zerowej części rzeczywistej pierwszej pary.
    {IntRead<long double>("2.0"), IntRead<long double>("2.0"), IntRead<long double>("3.0"), IntRead<long double>("2.0"), IntRead<long double>("1.0")},

    // 5. x^4 + 3x^3 + 3x^2 - 37x - 78 = 0
    // Pierwiastki: 3.0, -2.0, -2.0 +/- 3i
    // Solidny test deflacji. Duże współczynniki na końcu (-78, -37). Startowa heurystyka dla r i q (wyciąganie Mid() 
    // z najwyższych współczynników) musi tutaj dobrze zadziałać, żeby nie wpaść w zbyt długą pętlę Newtona.
    {IntRead<long double>("-78.0"), IntRead<long double>("-37.0"), IntRead<long double>("3.0"), IntRead<long double>("3.0"), IntRead<long double>("1.0")},

    // --- STOPIEŃ 5 ---

    // 6. 2x^5 - 21x^4 + 64x^3 - 31x^2 - 78x + 40 = 0
    // Pierwiastki: -1.0, 0.5, 2.0, 4.0, 5.0
    // 5 odrębnych pierwiastków rzeczywistych o różnym "rozstrzale" (od -1 do 5, w tym ułamek 0.5). 
    // Zmusza to metodę Bairstowa do sfabrykowania dwóch czynników kwadratowych ze zlepku tych pierwiastków, 
    // co przy stopniu nieparzystym i współczynniku kierunkowym równym 2 jest doskonałym sprawdzianem dla całkowitej dokładności programu.
    {IntRead<long double>("40.0"), IntRead<long double>("-78.0"), IntRead<long double>("-31.0"), IntRead<long double>("64.0"), IntRead<long double>("-21.0"), IntRead<long double>("2.0")}
    };
    int i = 1;

    for(const auto& coefficients : tests)
    {
        vector<complex_interval<Interval<long double>>> roots = bairstow_method<Interval<long double>>(coefficients.size() - 1, coefficients,150, 1e-11,1e-11);

        string reLeft, reRight, imLeft, imRight;
        for (const auto& root : roots)
        {
            // Konwersja części rzeczywistej
            const_cast<Interval<long double>&>(root.real).IEndsToStrings(reLeft, reRight);
            
            // Konwersja części urojonej
            const_cast<Interval<long double>&>(root.imag).IEndsToStrings(imLeft, imRight);

            cout << "Root: " << i << " [" << reLeft << ", " << reRight << "] + [" 
                << imLeft << ", " << imRight << "]i" << " new ";
        }
        cout << endl;
    }


}


