#include <iostream>
#include "../../include/Interval.h"
#include <mpfr.h>
#include "../../include/mpreal.h"

using namespace std;
using namespace interval_arithmetic;


int main()
{

    string arithmetic;
    string counting;
    string input_type;
    string polynomial;

    cout << "RUNNING" << endl;

    while (true)
    {
        getline(cin, arithmetic);
        getline(cin, counting);
        getline(cin, input_type);
        getline(cin, polynomial);

        cout << "DZIAŁA" << endl;
    }


}





/*
int main() {
    // 1. Inicjalizacja biblioteki dla typu mpreal (biblioteka MPFR)
    // To ustawia domyślną precyzję zdefiniowaną w Interval.h
    Interval<mpreal>::Initialize();
    Interval<mpreal>::SetMode(PINT_MODE);

    cout << "--- Test Arytmetyki Przedzialowej ---" << endl;

    // 2. Tworzenie przedziałów przez wczytanie napisów (bezpieczniejsze dla precyzji)
    // IntRead tworzy przedział zawierający daną liczbę (z uwzględnieniem zaokrągleń)
    Interval<mpreal> A = IntRead<mpreal>("1.0");
    Interval<mpreal> B = IntRead<mpreal>("2.0");
    
    // Ręczne tworzenie przedziału [0.1, 0.2]
    Interval<mpreal> X(0.1, 0.2);

    // 3. Wykonywanie operacji
    Interval<mpreal> Suma = A + B;
    Interval<mpreal> Iloczyn = X * B;

    // 4. Wyświetlanie wyników
    string left, right;
    
    Suma.IEndsToStrings(left, right);
    cout << "A + B = [" << left << ", " << right << "]" << endl;

    Iloczyn.IEndsToStrings(left, right);
    cout << "X * B = [" << left << ", " << right << "]" << endl;

    // 5. Przykład funkcji matematycznej (np. Sinus)
    Interval<mpreal> S = ISin(X);
    S.IEndsToStrings(left, right);
    cout << "sin(X) = [" << left << ", " << right << "]" << endl;

    return 0;
}*/