/*
 * Interval.h
 *
 * Created on: 27 lut 2021
 * Author: tomhof
 *
 * Dokumentacja i analiza działania:
 * Plik definiuje szablonową bibliotekę do obliczeń w arytmetyce przedziałowej 
 * (klasycznej oraz skierowanej - Kaucher'a). Głównym założeniem biblioteki jest 
 * zagwarantowanie, że rzeczywisty wynik operacji zawsze znajduje się wewnątrz 
 * zwracanego przedziału. Osiąga się to poprzez ścisłą kontrolę kierunku 
 * zaokrągleń zmiennoprzecinkowych (w dół dla lewego krańca, w górę dla prawego).
 * Wspiera typy wbudowane (float, double, long double) oraz typy o dowolnej
 * precyzji (mpreal z biblioteki MPFR).
 */

#ifndef INTERVAL_H_
#define INTERVAL_H_

// Zabezpieczenia dla biblioteki MPFR, aby unikać konfliktów makr.
#ifndef MPFR_USE_NO_MACRO
  #define MPFR_USE_NO_MACRO
#endif
#ifndef MPFR_USE_INTMAX_T
  #define MPFR_USE_INTMAX_T
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fenv.h>     // Krytyczne: dostarcza fesetround() do zmiany kierunku zaokrągleń FPU
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <mpfr.h>     // GNU Multiple Precision Floating-Point Reliably
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <float.h>
#include <typeinfo>
#include "mpreal.h"   // Wrapper C++ dla MPFR

using namespace std;
using namespace mpfr;

namespace interval_arithmetic {

/**
 * @brief Precyzja w bitach dla poszczególnych typów danych.
 * Używane głównie do konfiguracji wewnętrznej biblioteki MPFR.
 */
enum IAPrecision {
	LONGDOUBLE_PREC = 63, DOUBLE_PREC = 53, FLOAT_PREC = 32, MPREAL_PREC = 40
};

/**
 * @brief Liczba cyfr znaczących używana podczas rzutowania przedziału na ciąg znaków (string).
 */
enum IAOutDigits {
	LONGDOUBLE_DIGITS = 17, DOUBLE_DIGITS = 16, FLOAT_DIGITS = 7
};

/**
 * @brief Tryb działania arytmetyki przedziałowej.
 * - PINT_MODE (Proper Interval): Standardowa arytmetyka, przedziały [a, b] gdzie a <= b.
 * - DINT_MODE (Directed Interval): Skierowana arytmetyka Kauchera, dopuszcza przedziały niewłaściwe a > b.
 */
enum IAMode {
	DINT_MODE, PINT_MODE
};

// Deklaracja zapowiadająca (forward declaration) głównej klasy
template<typename T> class Interval;

// --- ZAPOWIEDZI FUNKCJI GLOBALNYCH ---
// Pomocnicze funkcje do wczytywania przedziałów z ciągów znaków
template<typename T> Interval<T> IntRead(const string &sa);
template<typename T> T LeftRead(const string &sa);
template<typename T> T RightRead(const string &sa);

// Obliczanie szerokości przedziału
template<typename T> T DIntWidth(const Interval<T> &x);
template<typename T> T IntWidth(const Interval<T> &x);

// --- WŁAŚCIWA ARYTMETYKA PRZEDZIAŁOWA (I*) ---
template<typename T> Interval<T> IAdd(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> ISub(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> IDiv(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> IMul(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> ISin(const Interval<T> &x);
template<typename T> Interval<T> ICos(const Interval<T> &x);
template<typename T> Interval<T> IExp(const Interval<T> &x);

// --- SKIEROWANA ARYTMETYKA PRZEDZIAŁOWA KAUCHERA (DI*) ---
template<typename T> Interval<T> DIAdd(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> DISub(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> DIDiv(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> DIMul(const Interval<T> &x, const Interval<T> &y);
template<typename T> Interval<T> DISin(const Interval<T> &x);
template<typename T> Interval<T> DICos(const Interval<T> &x);
template<typename T> Interval<T> DIExp(const Interval<T> &x);

// Narzędzia
template<typename T> Interval<T> Hull(const Interval<T> &x, const Interval<T> &y); // Otoczka przedziałowa
template<typename T> Interval<T> IAbs(const Interval<T> &x); // Wartość bezwzględna z przedziału

/**
 * @brief Globalna funkcja do bezpiecznej zmiany trybu zaokrąglania procesora (FPU).
 * Gwarantuje, że wyniki obliczeń zmiennoprzecinkowych odpowiednio zaniżą/zawyżą wynik.
 */
template<typename T> int SetRounding(int rounding);
template<> Interval<mpreal> IntRead(const string &sa);

/**
 * @class Interval
 * @brief Główna klasa reprezentująca przedział matematyczny [a, b].
 */
template<typename T> class Interval {
private:
	static IAPrecision precision; // Aktualnie ustawiona precyzja
	static IAOutDigits outdigits; // Cyfry do wyświetlania

public:
	static IAMode mode; // Tryb działania (PINT / DINT), domyślnie współdzielony dla całej aplikacji
	
	T a; // Lewy kraniec przedziału (infimum / dolne ograniczenie)
	T b; // Prawy kraniec przedziału (supremum / górne ograniczenie)
	
	// Konstruktory i destryktor
	Interval();
	Interval(Interval const &copy);
	Interval(T a, T b);
	virtual ~Interval();
	
	// Operatory przypisania i operatory arytmetyczne
	Interval& operator=(const Interval<T> i);
	Interval operator+(const Interval<T> &i);
	Interval operator-(const Interval<T> &i);
	Interval operator*(const Interval<T> &i);
	Interval operator*(const long double &l);
	Interval operator*(const int &i);
	Interval operator/(const Interval<T> &i);
	
	/**
	 * @brief Projekcja rzutuje przedział niewłaściwy (a > b) na właściwy (a <= b) zamieniając krańce.
	 */
	Interval Projection();
	
	/**
	 * @brief Tworzy przedział przeciwny: [-a, -b]. Zauważ, że zachowuje on kolejność (nie jest to standardowa arytmetyka ujemna [-b, -a]).
	 */
	Interval Opposite();
	
	/**
	 * @brief Operator dualności w arytmetyce Kauchera. Zwraca [b, a]. 
	 * Pozwala na rozwiązywanie niektórych równań przedziałowych.
	 */
	Interval Dual();
	
	/**
	 * @brief Zwraca przedział odwrotny 1/X z rygorystycznym zaokrąglaniem, dobierając najszerszy możliwy przedział w trybie DINT.
	 */
	Interval Inverse();
	
	T Mid(); // Środek przedziału
	T GetWidth(); // Szerokość przedziału (odległość od a do b)
	
	static void Initialize();
	
	// Predefiniowane stałe przedziałowe (gwarantują, że stała leży ściśle wewnątrz wygenerowanego przedziału)
	static Interval<T> ISqr2();
	static Interval<T> ISqr3();
	static Interval<T> IPi();
	
	// Settery i gettery statyczne
	static void SetMode(IAMode m) { mode = m; }
	static IAMode GetMode();
	static void SetPrecision(IAPrecision p);
	static IAPrecision GetPrecision();
	static void SetOutDigits(IAOutDigits o);
	static IAOutDigits GetOutDigits();
	static T GetEpsilon();

	/**
	 * @brief Formatuje oba krańce przedziału do pary stringów.
	 */
	void IEndsToStrings(string &left, string &right);

	// Przyjaźnie ułatwiają dostęp do struktury wenętrznej wolnym funkcjom
	friend T DIntWidth<T>(const Interval &x);
	friend T IntWidth<T>(const Interval &x);
	friend Interval IAdd<T>(const Interval &x, const Interval &y);
	friend Interval ISub<T>(const Interval &x, const Interval &y);
	friend Interval IDiv<T>(const Interval &x, const Interval &y);
	friend Interval IMul<T>(const Interval &x, const Interval &y);
	friend Interval ISin<T>(const Interval<T> &x);
	friend Interval ICos<T>(const Interval<T> &x);
	friend Interval IExp<T>(const Interval<T> &x);
	friend Interval IntRead<T>(const string &sa);
	friend T LeftRead<T>(const string &sa);
	friend T RightRead<T>(const string &sa);

	friend Interval DIAdd<T>(const Interval &x, const Interval &y);
	friend Interval DISub<T>(const Interval &x, const Interval &y);
	friend Interval DIDiv<T>(const Interval &x, const Interval &y);
	friend Interval DIMul<T>(const Interval &x, const Interval &y);
	friend Interval DISin<T>(const Interval &x);
	friend Interval DICos<T>(const Interval &x);
	friend Interval DIExp<T>(const Interval &x);

	friend Interval IAbs<T>(const Interval &x);
	friend Interval Hull<T>(const Interval &x, const Interval &y);

	friend int SetRounding<T>(int rounding);
};

// --- IMPLEMENTACJE KONSTRUKTORÓW I METOD KLASY ---

template<typename T>
inline Interval<T>::~Interval() {
}

template<typename T>
Interval<T>::Interval() {
	this->a = 0;
	this->b = 0;
}

template<typename T>
Interval<T>::Interval(Interval const &copy) {
	this->a = copy.a;
	this->b = copy.b;
}

template<typename T>
inline Interval<T>::Interval(T a, T b) {
	this->a = a;
	this->b = b;
}

template<typename T>
inline IAMode Interval<T>::GetMode() {
	return Interval<T>::mode;
}

/**
 * @brief Ustawia kierunek zaokrągleń dla sprzętowego FPU.
 * Jest to najważniejsza funkcja utrzymująca spójność przedziałową -
 * zapobiega utracie dokładności, która mogłaby spowodować wyrzucenie
 * poprawnego wyniku poza granice przedziału.
 */
template<typename T>
int SetRounding(int rounding) {
	fesetround(rounding); // Oczekuje makr standardu C: FE_DOWNWARD, FE_UPWARD, FE_TONEAREST
	return rounding;
}

template<typename T>
inline Interval<T>& Interval<T>::operator =(Interval<T> i) {
	std::swap(this->a, i.a);
	std::swap(this->b, i.b);
	return *this;
}

template<typename T>
inline void Interval<T>::SetPrecision(IAPrecision p) {
	mpreal::set_default_prec(p);
	Interval<T>::precision = p;
}

template<typename T>
inline IAPrecision Interval<T>::GetPrecision() {
	return Interval<T>::precision;
}

template<typename T>
inline void Interval<T>::SetOutDigits(IAOutDigits o) {
	Interval<T>::outdigits = LONGDOUBLE_DIGITS; // BUG?: Twarde zakodowanie pomimo argumentu 'o'.
}

template<typename T>
inline IAOutDigits Interval<T>::GetOutDigits() {
	return Interval<T>::outdigits;
}

template<typename T>
inline T Interval<T>::GetEpsilon() {
	return std::numeric_limits<T>::epsilon();
}

/**
 * @brief Czyta przedział ze zmiennej tekstowej, uwzględniając zaokrąglenia.
 * Prawidłowo konwertuje wprowadzony tekst na reprezentację zmiennoprzecinkową, 
 * kierując lewy kraniec w dół, a prawy w górę, zabezpieczając się przed błędami 
 * konwersji systemu dziesiętnego na binarny.
 */
template<typename T>
inline Interval<T> IntRead(const string &sa) {
	Interval<T> r;
	mpfr_t rop;
	mpfr_init2(rop, Interval<T>::precision);
	
	// Parsowanie dla lewego krańca - celowo zaniżamy dokładność (RNDD - Round Down)
	mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDD);
	T le = 0.0;
	if (strcmp(typeid(T).name(), typeid(long double).name()) == 0) {
		le = mpfr_get_ld(rop, MPFR_RNDD);
	}
	if (strcmp(typeid(T).name(), typeid(double).name()) == 0) {
		le = mpfr_get_d(rop, MPFR_RNDD);
	}
	if (strcmp(typeid(T).name(), typeid(float).name()) == 0) {
		le = mpfr_get_flt(rop, MPFR_RNDD);
	}

	// Parsowanie dla prawego krańca - celowo zawyżamy (RNDU - Round Up)
	mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDU);
	T re = 0.0;
	if (strcmp(typeid(T).name(), typeid(long double).name()) == 0) {
		re = mpfr_get_ld(rop, MPFR_RNDU);
	}
	if (strcmp(typeid(T).name(), typeid(double).name()) == 0) {
		re = mpfr_get_d(rop, MPFR_RNDU);
	}
	if (strcmp(typeid(T).name(), typeid(float).name()) == 0) {
		re = mpfr_get_flt(rop, MPFR_RNDU);
	}

	SetRounding<T>(FE_TONEAREST); // Powrót do domyślnego trybu procesora
	r.a = le;
	r.b = re;
	return r;
}

// Specjalizacja IntRead dla typu wysokiej precyzji mpreal
template<>
inline Interval<mpreal> IntRead(const string &sa) {
	Interval<mpreal> r;
	mpfr_t rop;
	mpfr_init2(rop, Interval<mpreal>::precision);
	mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDD);
	mpreal le = rop;

	mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDU);
	mpreal re = rop;

	r.a = le;
	r.b = re;
	return r;
}

/**
 * @brief Zamienia krańce a i b na reprezentację znakową.
 * Formatowanie w konwencji notacji naukowej.
 */
template<typename T>
inline void Interval<T>::IEndsToStrings(string &left, string &right) {
	mpfr_t rop;
	mpfr_exp_t exponent;
	mpfr_init2(rop, precision);
	char *str = NULL;
	char *buffer = new char(precision + 3);
	
	// Formatowanie dolnego krańca z zaokrągleniem w dół
	mpfr_set_ld(rop, this->a, MPFR_RNDD);
	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDD);
	str = buffer;

	stringstream ss;
	int prec = std::numeric_limits<T>::digits10;
	ss.setf(std::ios_base::scientific);
	bool minus = (str[0] == '-');
	int splitpoint = minus ? 1 : 0;
	string sign = minus ? "-" : "";

	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	left = ss.str();
	ss.str(std::string()); // Czyszczenie strumienia

	// Formatowanie górnego krańca z zaokrągleniem w górę
	mpfr_set_ld(rop, this->b, MPFR_RNDU);
	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDU);
	str = buffer;
	splitpoint = (str[0] == '-') ? 1 : 0;
	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	right = ss.str();
	ss.clear();
}

template<typename T>
inline Interval<T> Interval<T>::Projection() {
	Interval<T> x(this->a, this->b);
	Interval<T> r;
	r = x;
	if (x.a > x.b) {
		r.a = x.b;
		r.b = x.a;
	}
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::Opposite() {
	Interval<T> x(this->a, this->b);
	Interval<T> r;
	r.a = -x.a;
	r.b = -x.b;
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::Dual() {
	Interval<T> x(this->a, this->b);
	Interval<T> r;
	r.a = x.b;
	r.b = x.a;
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::Inverse() {
	Interval<T> x(this->a, this->b);
	Interval<T> z1, z2;

	// Obliczanie odwrotności z rygorem zaokrągleń
	SetRounding<T>(FE_DOWNWARD);
	z1.a = 1 / x.a;
	z2.b = 1 / x.b;
	
	SetRounding<T>(FE_UPWARD);
	z1.b = 1 / x.b;
	z2.a = 1 / x.a;
	
	SetRounding<T>(FE_TONEAREST);
	
	// W trybie skierowanym sprawdzamy, która wersja jest bezpieczniejsza (szersza)
	if (DIntWidth(z1) >= DIntWidth(z2))
		return z1;
	else
		return z2;
}

template<typename T>
inline T LeftRead(const string &sa) {
	Interval<T> int_number;
	int_number = IntRead<T>(sa);
	return int_number.a;
}

template<typename T>
inline T Interval<T>::GetWidth() {
	Interval<T> x(this->a, this->b);
	switch (mode) {
	case PINT_MODE:
		return IntWidth(x);
	case DINT_MODE:
		return DIntWidth(x);
	default:
		return IntWidth(x);
	}
}

template<typename T>
inline T Interval<T>::Mid() {
	return (this->b + this->a) / 2.0; // Nieuwzględnia rygorystycznych zaokrągleń dla średniej, gdyż zazwyczaj służy ona jedynie do punktów początkowych iteracji
}

// Generatory najczęstszych stałych matematycznych. Zwracają bezpieczny przedział pokrywający prawdziwą wartość stałej.
template<typename T>
inline Interval<T> Interval<T>::ISqr2() {
	string i2;
	Interval<T> r;
	i2 = "1.414213562373095048";
	r.a = LeftRead<T>(i2);
	i2 = "1.414213562373095049";
	r.b = RightRead<T>(i2);
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::ISqr3() {
	string i2;
	Interval<T> r;
	i2 = "1.732050807568877293";
	r.a = LeftRead<T>(i2);
	i2 = "1.732050807568877294";
	r.b = RightRead<T>(i2);
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::IPi() {
	string i2;
	Interval<T> r;
	i2 = "3.141592653589793238";
	r.a = LeftRead<T>(i2);
	i2 = "3.141592653589793239";
	r.b = RightRead<T>(i2);
	return r;
}

template<typename T>
inline void Interval<T>::Initialize() {
	if (strcmp(typeid(T).name(), typeid(long double).name()) == 0) {
		Interval<T>::SetPrecision(LONGDOUBLE_PREC);
		Interval<T>::SetOutDigits(LONGDOUBLE_DIGITS);
	}

	if (strcmp(typeid(T).name(), typeid(double).name()) == 0) {
		Interval<T>::SetPrecision(DOUBLE_PREC);
		Interval<T>::SetOutDigits(DOUBLE_DIGITS);
	}

	if (strcmp(typeid(T).name(), typeid(float).name()) == 0) {
		Interval<T>::SetPrecision(FLOAT_PREC);
		Interval<T>::SetOutDigits(FLOAT_DIGITS);
	}

	if (strcmp(typeid(T).name(), typeid(mpreal).name()) == 0) {
		Interval<T>::SetPrecision(MPREAL_PREC);
		Interval<T>::SetOutDigits(FLOAT_DIGITS);
	}

}

template<typename T>
inline T RightRead(const string &sa) {
	Interval<T> int_number;
	int_number = IntRead<T>(sa);
	return int_number.b;
}

/**
 * @brief Zwraca szerokość w standardowej arytmetyce właściwej (PINT).
 * Upewnia się, że obliczona szerokość nie jest zaniżona przez ucięcie FP.
 */
template<typename T>
T IntWidth(const Interval<T> &x) {
	SetRounding<T>(FE_UPWARD);
	T w = x.b - x.a;
	SetRounding<T>(FE_TONEAREST);
	return w;
}

/**
 * @brief Zwraca szerokość przedziału dla arytmetyki Kauchera (DINT).
 * Ponieważ może być 'a > b', bierze pod uwagę moduł szerokości w obu kierunkach
 * zaokrąglenia i wybiera większą wartość dla zachowania bezpieczeństwa granic.
 */
template<typename T>
T DIntWidth(const Interval<T> &x) {
	long double w1, w2;

	SetRounding<T>(FE_UPWARD);
	w1 = x.b - x.a;
	if (w1 < 0)
		w1 = -w1;
	SetRounding<T>(FE_DOWNWARD);
	w2 = x.b - x.a;
	if (w2 < 0)
		w2 = -w2;
	SetRounding<T>(FE_TONEAREST);
	if (w1 > w2)
		return w1;
	else
		return w2;
}

// ====================================================================================
// ========================== STANDARDOWA ARYTMETYKA PRZEDZIAŁOWA =====================
// ====================================================================================

/**
 * @brief Implementacja przedziałowego Sinusa poprzez rozwinięcie w szereg Taylora (Maclaurina).
 */
template<typename T>
Interval<T> ISin(const Interval<T> &x) {
	bool is_even, finished;
	int k;
	int st = 0;
	Interval<T> d, s, w, w1, x2, tmp;
	Interval<T> izero(0, 0);
	string left, right;
	T eps = 1E-18; // Tolerancja dla testu zbieżności szeregu
	T diff = std::numeric_limits<T>::max();
	
	if (x.a > x.b)
		st = 1; // Błąd: przedział niewłaściwy nie jest tu obsługiwany
	else {
		s = x;
		w = x;
		x2 = IMul(x, x);
		k = 1;
		is_even = true;
		finished = false;
		st = 0;

		// Pętla dodająca kolejne wyrazy szeregu naprzemiennego
		do {
			d.a = (k + 1) * (k + 2);
			d.b = d.a;
			s = IMul(s, IDiv(x2, d)); // s = s * x^2 / ((k+1)(k+2))
			if (is_even)
				w1 = ISub(w, s);
			else
				w1 = IAdd(w, s);

			if ((w.a == 0) && (w.b == 0)) {
				return izero;
			}

			// Sprawdzenie kryterium stopu (różnica względna / bezwzględna pomiędzy iteracjami < epsilon)
			if ((w.a != 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) / abs(w.a) < eps)
						&& (abs(w.b - w1.b) / abs(w.b) < eps))
					finished = true;
				else
					;
			} else if ((w.a == 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) < eps)
						&& (abs(w.b - w1.b) / abs(w.b) < eps))
					finished = true;
				else
					;
			} else if (w.a != 0) {
				if ((abs(w.a - w1.a) / abs(w.a) < eps)
						& (abs(w.b - w1.b) < eps))
					finished = true;
				else if ((abs(w.a - w1.a) < eps) & (abs(w.b - w1.b) < eps))
					finished = true;
			}

			if (finished) {
				// Korekta, aby sin(x) nigdy nie uciekł poza matematyczne ograniczenie obwiedni [-1, 1]
				if (w1.b > 1) {
					w1.b = 1;
					if (w1.a > 1)
						w1.a = 1;
				}
				if (w1.a < -1) {
					w1.a = -1;
					if (w1.b < -1)
						w1.b = -1;
				}
				return w1;
			} else {
				w = w1;
				k = k + 2;
				is_even = !is_even;
			}
		} while (!(finished || (k > INT_MAX / 2)));
	}
	if (!finished)
		st = 2; // Osiągnięto limit iteracji bez zbieżności

	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Przedziałowy Cosinus z rozwinięcia szeregu Maclaurina.
 * Analogicznie jak w ISin().
 */
template<typename T>
Interval<T> ICos(const Interval<T> &x) {
	Interval<T> c, d, w, w1, x2;
	int k, st;
	bool is_even, finished;

	c.a = 1;
	c.b = 1;
	w = c;
	x2 = IMul(x, x);
	k = 1;
	is_even = true;
	finished = false;
	st = 0;

	do {
		d.a = k * (k + 1);
		d.b = d.a;
		c = IMul(c, IDiv(x2, d));
		if (is_even) {
			w1 = ISub(w, c);
		} else {
			w1 = IAdd(w, c);
		}
		
		// Kryterium zbieżności względnej
		if ((w.a != 0) && (w.b != 0)) {
			if (((abs(w.a - w1.a) / abs(w.a)) < 1e-18)
					&& (abs(w.b - w1.b) / abs(w.b) < 1e-18))
				finished = true;
		} else {
			if ((w.a == 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) < 1e-18)
						&& (abs(w.b - w1.b) / abs(w.b) < 1e-18)) {
					finished = true;
				}
			} else if (w.a != 0) {
				if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
						&& (abs(w.b - w1.b) < 1e-18))
					finished = true;
			} else if ((abs(w.a - w1.a) < 1e-18) && (abs(w.b - w1.b) < 1e-18))
				finished = true;

		}
		
		if (finished) {
			// Twarde granice [-1, 1] dla cosinusa
			if (w1.b > 1) {
				w1.b = 1;
				if (w1.a > 1)
					w1.a = 1;
			}
			if (w1.a < -1) {
				w1.a = -1;
				if (w1.b < -1)
					w1.b = -1;
			}
			return w1;
		} else {
			w = w1;
			k = k + 2;
			is_even = !is_even;
		}

	} while (!finished || (k > INT_MAX / 2));

	if (!finished)
		st = 2;

	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Przedziałowa funkcja wykładnicza e^x (szereg Maclaurina).
 */
template<typename T>
Interval<T> IExp(const Interval<T> &x) {
	bool finished;
	int k;
	int st = 0;
	Interval<T> d, e, w, w1;
	string left, right;
	Interval<T> tmp = x;
	T eps = 1E-18; 
	T diff = std::numeric_limits<T>::max();
	
	// Przypadek specjalny: e^0 = 1
	if ((x.a < 0) && (x.b > 0))
		return {1,1}; 
		
	if (x.a > x.b)
		st = 1;
	else {
		e.a = 1;
		e.b = 1;
		w = e;
		k = 1;
		finished = false;
		st = 0;
		do {
			d.a = k;
			d.b = k;
			e = IMul(e, IDiv(x, d));
			w1 = IAdd(w, e);
			
			T oldMid = (w.a + w.b) / 2;
			T newMid = (w1.a + w1.b) / 2;
			T currDiff = abs(oldMid - newMid);
			T tmpDiff = diff - currDiff;
			diff = currDiff;
			
			if ((abs(w.a - w1.a) / abs(w.a) < eps)
					&& (abs(w.b - w1.b) / abs(w.b) < eps)) {
				finished = true;
				return w1;
			} else {
				w = w1;
				k = k + 1;
				if (k > 100000) {
					// Wypisanie informacji pomocniczej w przypadku braku zbieżności
					T wdth = w.GetWidth();

					tmp.IEndsToStrings(left, right);
					cout << "x=[" << left << "," << right << "]" << endl;
					w.IEndsToStrings(left, right);
					cout << "[" << left << "," << right << "]" << endl;
					cout << "      width =  " << std::setprecision(17) << wdth
							<< endl << " diff = " << diff << endl
							<< " tmpDiff = " << tmpDiff << endl << "eps = "
							<< eps << endl;
				}
			}
		} while (!(finished || (k > INT_MAX / 2)));
		if (!finished)
			st = 2;
	}
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Podnoszenie przedziału do kwadratu w bezpieczny sposób.
 * Własność: Jeżeli zero należy do przedziału wejściowego, to po podniesieniu
 * do kwadratu dolną granicą MUSI być 0. W przeciwnym razie jest to kwadrat
 * mniejszego elementu.
 */
template<typename T>
Interval<T> ISqr(const Interval<T> &x, int &st) {
	long double minx, maxx;
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	if (x.a > x.b)
		st = 1;
	else {
		st = 0;
		if ((x.a <= 0) && (x.b >= 0))
			minx = 0; // Przedział zawiera 0. Minimalna wartość kwadratu to 0.
		else if (x.a > 0)
			minx = x.a;
		else
			minx = x.b;
			
		if (abs(x.a) > abs(x.b))
			maxx = abs(x.a);
		else
			maxx = abs(x.b);
			
		SetRounding<T>(FE_DOWNWARD);
		r.a = minx * minx;
		SetRounding<T>(FE_UPWARD);
		r.b = maxx * maxx;
		SetRounding<T>(FE_TONEAREST);
	}
	return r;
}

/**
 * @brief Pierwiastkowanie przedziału.
 * Wymaga ostrożności: pierwiastek jest niezdefiniowany dla ujemnych przedziałów.
 */
template<typename T>
Interval<T> ISqrt(const Interval<T> &x, int &st) {
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	if (x.a > x.b) {
		st = 1;
	} else if (x.a < 0) {
		st = 2; // Błąd matematyczny: pierwiastek z liczby ujemnej
	} else {
		st = 0;
		SetRounding<T>(FE_DOWNWARD);
		r.a = std::sqrt(x.a);
		SetRounding<T>(FE_UPWARD);
		r.b = std::sqrt(x.b);
		SetRounding<T>(FE_TONEAREST);
	}

	return r;
}

/**
 * @brief Przedziałowe dodawanie właściwe: [a+c, b+d] 
 * Oblicza lewą krawędź zaokrąglając w dół, prawą w górę.
 */
template<typename T>
Interval<T> IAdd(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> r;
	SetRounding<T>(FE_DOWNWARD);
	r.a = x.a + y.a;
	SetRounding<T>(FE_UPWARD);
	r.b = x.b + y.b;
	SetRounding<T>(FE_TONEAREST);
	return r;
}

/**
 * @brief Przedziałowe odejmowanie właściwe: [a-d, b-c]
 */
template<typename T>
Interval<T> ISub(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> r;
	SetRounding<T>(FE_DOWNWARD);
	r.a = x.a - y.b; // Infimum to najmniejsza wartość: lewy - prawy
	SetRounding<T>(FE_UPWARD);
	r.b = x.b - y.a; // Supremum to największa wartość: prawy - lewy
	SetRounding<T>(FE_TONEAREST);
	return r;
}

/**
 * @brief Przedziałowe mnożenie właściwe.
 * Ponieważ nie znamy a priori znaków wartości (mogą być dowolne z kombinacji ujemne/dodatnie),
 * obliczamy wszystkie 4 iloczyny skrajne i wybieramy min() dla lewego i max() dla prawego,
 * kontrolując rygorystycznie zaokrąglenia na poziomie FPU.
 */
template<typename T>
Interval<T> IMul(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> r(0, 0);
	T x1y1, x1y2, x2y1;

	// OBLICZANIE DOLNEJ KRAWĘDZI (wyszukiwanie minimum)
	SetRounding<T>(FE_DOWNWARD);
	x1y1 = x.a * y.a;
	x1y2 = x.a * y.b;
	x2y1 = x.b * y.a;
	r.a = x.b * y.b;
	if (x2y1 < r.a) r.a = x2y1;
	if (x1y2 < r.a) r.a = x1y2;
	if (x1y1 < r.a) r.a = x1y1;

	// OBLICZANIE GÓRNEJ KRAWĘDZI (wyszukiwanie maksimum)
	SetRounding<T>(FE_UPWARD);
	x1y1 = x.a * y.a;
	x1y2 = x.a * y.b;
	x2y1 = x.b * y.a;

	r.b = x.b * y.b;
	if (x2y1 > r.b) r.b = x2y1;
	if (x1y2 > r.b) r.b = x1y2;
	if (x1y1 > r.b) r.b = x1y1;
	
	SetRounding<T>(FE_TONEAREST);
	return r;
}

/**
 * @brief Przedziałowe dzielenie właściwe. 
 * Chroni przed sytuacją dzielenia przez przedział obejmujący 0. 
 * Podobnie jak w mnożeniu, oblicza 4 kombinacje z odpowiednim trybem zaokrągleń.
 */
template<typename T>
Interval<T> IDiv(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> r;
	T x1y1, x1y2, x2y1, t;

	if ((y.a <= 0) && (y.b >= 0)) {
		throw runtime_error("Division by an interval containing 0.");
	} else {
		// Minimum z 4 kombinacji
		SetRounding<T>(FE_DOWNWARD);
		x1y1 = x.a / y.a;
		x1y2 = x.a / y.b;
		x2y1 = x.b / y.a;
		r.a = x.b / y.b;
		t = r.a;
		if (x2y1 < t) r.a = x2y1;
		if (x1y2 < t) r.a = x1y2;
		if (x1y1 < t) r.a = x1y1;

		// Maksimum z 4 kombinacji
		SetRounding<T>(FE_UPWARD);
		x1y1 = x.a / y.a;
		x1y2 = x.a / y.b;
		x2y1 = x.b / y.a;

		r.b = x.b / y.b;
		t = r.b;
		if (x2y1 > t) r.b = x2y1;
		if (x1y2 > t) r.b = x1y2;
		if (x1y1 > t) r.b = x1y1;

	}
	SetRounding<T>(FE_TONEAREST);
	return r;
}

// ====================================================================================
// ========================== SKIEROWANA ARYTMETYKA PRZEDZIAŁOWA ======================
// ====================================================================================

/**
 * @brief Dodawanie Kauchera (skierowane).
 * Zabezpiecza sytuację, w której któryś przedział może być niewłaściwy (a > b).
 * Wylicza dwie możliwości i zwraca przedział najszerszy (zachowanie gwarancji zawierania wyniku).
 */
template<typename T>
Interval<T> DIAdd(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> z1, z2;
	if ((x.a <= x.b) && (y.a <= y.b)) {
		return IAdd<T>(x, y); // Jeśli oba są właściwe, użyj zwykłego dodawania
	} else {
		SetRounding<T>(FE_DOWNWARD);
		z1.a = x.a + y.a;
		z2.b = x.b + y.b;
		SetRounding<T>(FE_UPWARD);
		z1.b = x.b + y.b;
		z2.a = x.a + y.a;
		SetRounding<T>(FE_TONEAREST);
		
		if (z1.GetWidth() >= z2.GetWidth())
			return z1;
		else
			return z2;
	}
}

/**
 * @brief Odejmowanie Kauchera (skierowane). 
 * Algorytm generuje obwiednię po uwzględnieniu potencjalnie niewłaściwych przedziałów wejściowych.
 */
template<typename T>
Interval<T> DISub(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> z1, z2;
	if ((x.a <= x.b) && (y.a <= y.b)) {
		return ISub(x, y);
	} else {
		SetRounding<T>(FE_DOWNWARD);
		z1.a = x.a - y.b;
		z2.b = x.b - y.a;
		SetRounding<T>(FE_UPWARD);
		z1.b = x.b - y.a;
		z2.a = x.a - y.b;
		SetRounding<T>(FE_TONEAREST);
		
		if (z1.GetWidth() >= z2.GetWidth())
			return z1;
		else
			return z2;
	}
}

/**
 * @brief Mnożenie Kauchera (skierowane). 
 * Ponieważ operuje na rozszerzonej strukturze algebraicznej, bada poszczególne
 * podzbiory płaszczyzny dwuwymiarowej (np. całkowicie dodatnie, ujemne, przecinające oś zer).
 * W implementacji jest to widoczne jako olbrzymia drabinka if-ów pokrywająca całą tabelę Kauchera.
 */
template<typename T>
Interval<T> DIMul(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> z1, z2, r;
	T z;
	bool xn, xp, yn, yp, zero;

	if ((x.a <= x.b) && (y.a <= y.b))
		r = IMul(x, y);
	else {
		// Logika identyfikacji stref przedziałów
		xn = (x.a < 0) and (x.b < 0);
		xp = (x.a > 0) and (x.b > 0);
		yn = (y.a < 0) and (y.b < 0);
		yp = (y.a > 0) and (y.b > 0);
		zero = false;
		
		// A, B w obszarach ściśle poza zerem
		if ((xn || xp) && (yn || yp))
			if (xp && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.a;
				z2.b = x.b * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.b;
				z2.a = x.a * y.a;
			} else if (xp && yn) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.a;
				z2.b = x.a * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.b;
				z2.a = x.b * y.a;
			} else if (xn && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.b;
				z2.b = x.b * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.a;
				z2.a = x.a * y.b;
			} else {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.b;
				z2.b = x.a * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.a;
				z2.a = x.b * y.b;
			}
		// Złożone tabele przejść Kauchera uwzględniające przejście przez 0
		else if ((xn || xp)
				&& (((y.a <= 0) && (y.b >= 0)) || ((y.a >= 0) && (y.b <= 0))))
			if (xp && (y.a <= y.b)) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.a;
				z2.b = x.b * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.b;
				z2.a = x.b * y.a;
			} else if (xp && (y.a > y.b)) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.a;
				z2.b = x.a * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.b;
				z2.a = x.a * y.a;
			} else if (xn && (y.a <= y.b)) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.b;
				z2.b = x.a * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.a;
				z2.a = x.a * y.b;
			} else {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.b;
				z2.b = x.b * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.a;
				z2.a = x.b * y.b;
			}
		// Złożone tabele przejść c.d.
		else if ((((x.a <= 0) && (x.b >= 0)) || ((x.a >= 0) && (x.b <= 0)))
				&& (yn || yp))
			if ((x.a <= x.b) && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.b;
				z2.b = x.b * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.b;
				z2.a = x.a * y.b;
			} else if ((x.a <= 0) && yn) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.a;
				z2.b = x.a * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.a;
				z2.a = x.b * y.a;
			} else if ((x.a > x.b) && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a * y.a;
				z2.b = x.b * y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b * y.a;
				z2.a = x.a * y.a;
			} else {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b * y.b;
				z2.b = x.a * y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a * y.b;
				z2.a = x.b * y.b;
			}
		// Operacje wokół stref zawierających samo zero
		else if ((x.a >= 0) && (x.b <= 0) && (y.a >= 0) && (y.b <= 0)) {
			SetRounding<T>(FE_DOWNWARD);
			z1.a = x.a * y.a;
			z = x.b * y.b;
			if (z1.a < z)
				z1.a = z;
			z2.b = x.a * y.b;
			z = x.b * y.a;
			if (z < z2.b)
				z2.b = z;
			SetRounding<T>(FE_UPWARD);
			z1.b = x.a * y.b;
			z = x.b * y.a;
			if (z < z1.b)
				z1.b = z;
			z2.a = x.a * y.a;
			z = x.b * y.b;
			if (z2.a < z)
				z2.a = z;
		}
		// A in Z and B in Z- or A in Z- and B in Z
		else
			zero = true;
			
		// Wybór bezpieczniejszego oszacowania z uwzględnieniem Zera 
		if (zero) {
			r.a = 0;
			r.b = 0;
		} else if (z1.GetWidth() >= z2.GetWidth())
			r = z1;
		else
			r = z2;
	}

	SetRounding<T>(FE_TONEAREST);
	return r;
}

/**
 * @brief Dzielenie Kauchera (skierowane). Implementuje podobną wieloprzypadkową mechanikę 
 * tabelaryczną co mnożenie skierowane. 
 */
template<typename T>
Interval<T> DIDiv(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> z1, z2, r;
	bool xn, xp, yn, yp, zero;

	if ((x.a <= x.b) && (y.a <= y.b))
		r = IDiv(x, y); // delegacja do zwykłego jeżeli przedziały są właściwe
	else {
		xn = (x.a < 0) && (x.b < 0);
		xp = (x.a > 0) && (x.b > 0);
		yn = (y.a < 0) && (y.b < 0);
		yp = (y.a > 0) && (y.b > 0);
		zero = false;
		// Drabinka reguł Kauchera dla dzielenia
		if ((xn || xp) && (yn || yp))
			if (xp && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a / y.b;
				z2.b = x.b / y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b / y.a;
				z2.a = x.a / y.b;
			} else if (xp && yn) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b / y.b;
				z2.b = x.a / y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a / y.a;
				z2.a = x.b / y.b;
			} else if (xn && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a / y.a;
				z2.b = x.b / y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b / y.b;
				z2.a = x.a / y.a;
			} else {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b / y.a;
				z2.b = x.a / y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a / y.b;
				z2.a = x.b / y.a;
			}
		// A in T, B in H-T
		else if (((x.a <= 0) && (x.b >= 0))
				|| (((x.a >= 0) && (x.b <= 0)) && (yn || yp)))
			if ((x.a <= x.b) && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a / y.a;
				z2.b = x.b / y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b / y.a;
				z2.a = x.a / y.a;
			} else if ((x.a <= x.b) && yn) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b / y.b;
				z2.b = x.a / y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a / y.b;
				z2.a = x.b / y.b;
			} else if ((x.a > x.b) && yp) {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.a / y.b;
				z2.b = x.b / y.b;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.b / y.b;
				z2.a = x.a / y.b;
			} else {
				SetRounding<T>(FE_DOWNWARD);
				z1.a = x.b / y.a;
				z2.b = x.a / y.a;
				SetRounding<T>(FE_UPWARD);
				z1.b = x.a / y.a;
				z2.a = x.b / y.a;
			}
		else
			zero = true;
			
		// Wyjątek dzielenia w sytuacji kolizji z zerem dla trybu niewłaściwego
		if (zero)
			throw runtime_error("Division by an interval containing 0.");
		else if (z1.GetWidth() >= z2.GetWidth())
			r = z1;
		else
			r = z2;
		SetRounding<T>(FE_TONEAREST);
	}
	return r;
}

/**
 * @brief Skierowany trygonometryczny Sinus bazujący wewnętrznie na zespole
 * operatorów skierowanych (DIMul, DIDiv, DIAdd, DISub).
 */
template<typename T>
Interval<T> DISin(const Interval<T> &x) {
	bool is_even, finished;
	int k;
	int st = 0;
	Interval<T> d, s, w, w1, x2;
	if (x.a > x.b)
		st = 1;
	else {
		s = x;
		w = x;
		x2 = DIMul(x, x);
		k = 1;
		is_even = true;
		finished = false;
		st = 0;

		do {
			d.a = (k + 1) * (k + 2);
			d.b = d.a;
			s = DIMul(s, DIDiv(x2, d));
			if (is_even)
				w1 = DISub(w, s);
			else
				w1 = DIAdd(w, s);
				
			if ((w.a != 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
						&& (abs(w.b - w1.b) / abs(w.b) < 1e-18))
					finished = true;
				else
					;
			} else if ((w.a == 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) < 1e-18)
						&& (abs(w.b - w1.b) / abs(w.b) < 1e-18))
					finished = true;
				else
					;
			}

			else if (w.a != 0) {
				if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
						& (abs(w.b - w1.b) < 1e-18))
					finished = true;
				else if ((abs(w.a - w1.a) < 1e-18) & (abs(w.b - w1.b) < 1e-18))
					finished = true;
			}

			if (finished) {
				if (w1.b > 1) {
					w1.b = 1;
					if (w1.a > 1)
						w1.a = 1;
				}
				if (w1.a < -1) {
					w1.a = -1;
					if (w1.b < -1)
						w1.b = -1;
				}
				return w1;
			} else {
				w = w1;
				k = k + 2;
				is_even = !is_even;
			}
		} while (!(finished || (k > INT_MAX / 2)));
	}
	if (!finished)
		st = 2;

	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Skierowany trygonometryczny Cosinus bazujący na operatorach Kauchera.
 */
template<typename T>
Interval<T> DICos(const Interval<T> &x) {
	bool is_even, finished;
	int k, st;
	Interval<T> d, c, w, w1, x2;
	if (x.a > x.b)
		st = 1;
	else {
		c.a = 1;
		c.b = 1;
		w = c;
		x2 = DIMul(x, x);
		k = 1;
		is_even = true;
		finished = false;
		st = 0;

		do {
			d.a = k * (k + 1);
			d.b = d.a;
			c = DIMul(c, DIDiv(x2, d));
			if (is_even)
				w1 = DISub(w, c);
			else
				w1 = DIAdd(w, c);

			if ((w.a != 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
						&& (abs(w.b - w1.b) / abs(w.b) < 1e-18))
					finished = true;
				else
					;
			} else if ((w.a == 0) && (w.b != 0)) {
				if ((abs(w.a - w1.a) < 1e-18)
						&& (abs(w.b - w1.b) / abs(w.b) < 1e-18))
					finished = true;
				else
					;
			}

			else if (w.a != 0) {
				if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
						& (abs(w.b - w1.b) < 1e-18))
					finished = true;
				else if ((abs(w.a - w1.a) < 1e-18) & (abs(w.b - w1.b) < 1e-18))
					finished = true;
			}

			if (finished) {
				if (w1.b > 1) {
					w1.b = 1;
					if (w1.a > 1)
						w1.a = 1;
				}
				if (w1.a < -1) {
					w1.a = -1;
					if (w1.b < -1)
						w1.b = -1;
				}
				return w1;
			} else {
				w = w1;
				k = k + 2;
				is_even = !is_even;
			}
		} while (!(finished || (k > INT_MAX / 2)));
	}
	if (!finished)
		st = 2;

	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Skierowana funkcja wykładnicza oparta o uogólnione operatory DIMul, DIAdd etc.
 */
template<typename T>
Interval<T> DIExp(const Interval<T> &x) {
	bool finished;
	int k;
	int st = 0;
	Interval<T> d, e, w, w1;
	if (x.a > x.b)
		st = 1;
	else {
		e.a = 1;
		e.b = 1;
		w = e;
		k = 1;
		finished = false;
		st = 0;
		do {
			d.a = k;
			d.b = k;
			e = IMul(e, DIDiv(x, d));
			w1 = DIAdd(w, e);
			if ((abs(w.a - w1.a) / abs(w.a) < 1e-18)
					&& (abs(w.b - w1.b) / abs(w.b) < 1e-18)) {
				finished = true;
				return w1;
			} else {
				w = w1;
				k = k + 1;
			}
		} while (!(finished || (k > INT_MAX / 2)));
		if (!finished)
			st = 2;
	}
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	return r;
}

/**
 * @brief Skierowany kwadrat uogólniający standardowy proces.
 */
template<typename T>
Interval<T> DISqr(const Interval<T> &x) {
	long double minx, maxx;
	int st = 0;
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	if (x.a > x.b)
		st = 1;
	else {
		st = 0;
		if ((x.a <= 0) && (x.b >= 0))
			minx = 0;
		else if (x.a > 0)
			minx = x.a;
		else
			minx = x.b;
		if (abs(x.a) > abs(x.b))
			maxx = abs(x.a);
		else
			maxx = abs(x.b);
		SetRounding<T>(FE_DOWNWARD);
		r.a = minx * minx;
		SetRounding<T>(FE_UPWARD);
		r.b = maxx * maxx;
		SetRounding<T>(FE_TONEAREST);
	}
	return r;
}

// ====================================================================================
// ==================== PRZECIĄŻENIA OPERATORÓW DLA KLASY INTERVAL ====================
// ====================================================================================

/**
 * @brief Operatory dodawania, w zależności od globalnego trybu wywołują IAdd (właściwe) lub DIAdd (Kauchera).
 */
template<typename T>
inline Interval<T> Interval<T>::operator +(const Interval<T> &y) {
	Interval<T> x(this->a, this->b);
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = IAdd<T>(x, y);
		break;
	case DINT_MODE:
		r = DIAdd<T>(x, y);
		break;
	default:
		r = IAdd<T>(x, y);
		break;
	}

	return r;
}

template<typename T>
inline Interval<T> operator +(Interval<T> x, const Interval<T> &y) {
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IAdd<T>(x, y);
	case DINT_MODE:
		return DIAdd<T>(x, y);
	default:
		return IAdd<T>(x, y);
	}
}

/**
 * @brief Operatory odejmowania, również zależne od globalnego trybu arytmetyki.
 */
template<typename T>
inline Interval<T> Interval<T>::operator -(const Interval<T> &y) {
	Interval<T> x(this->a, this->b);
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = ISub<T>(x, y);
		break;
	case DINT_MODE:
		r = DISub<T>(x, y);
		break;
	default:
		r = ISub<T>(x, y);
		break;
	}

	return r;
}

template<typename T>
inline Interval<T> operator -(Interval<T> x, const Interval<T> &y) {
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return ISub<T>(x, y);
	case DINT_MODE:
		return DISub<T>(x, y);
	default:
		return ISub<T>(x, y);
	}
}

/**
 * @brief Operatory mnożenia przedziału przez przedział i przez skalary.
 */
template<typename T>
inline Interval<T> Interval<T>::operator *(const Interval<T> &y) {
	Interval<T> x(this->a, this->b);
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = IMul<T>(x, y);
		break;
	case DINT_MODE:
		r = DIMul<T>(x, y);
		break;
	default:
		r = IMul<T>(x, y);
		break;
	}

	return r;
}

template<typename T>
inline Interval<T> operator *(int i, const Interval<T> &y) {
	Interval<T> x = { i, i }; // Niejawna konwersja skalara (punktu) na przedział zdegradowany do punktu
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IMul<T>(x, y);
	case DINT_MODE:
		return DIMul<T>(x, y);
	default:
		return IMul<T>(x, y);
	}
}

template<typename T>
inline Interval<T> operator *(const Interval<T> &y, int i) {
	Interval<T> x = { i, i };
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IMul<T>(x, y);
	case DINT_MODE:
		return DIMul<T>(x, y);
	default:
		return IMul<T>(x, y);
	}
}

template<typename T>
inline Interval<T> operator *(long double i, const Interval<T> &y) {
	Interval<T> x = { i, i };
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IMul<T>(x, y);
	case DINT_MODE:
		return DIMul<T>(x, y);
	default:
		return IMul<T>(x, y);
	}
}

template<typename T>
inline Interval<T> operator *(const Interval<T> &y, long double i) {
	Interval<T> x = { i, i };
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IMul<T>(x, y);
	case DINT_MODE:
		return DIMul<T>(x, y);
	default:
		return IMul<T>(x, y);
	}
}

template<typename T>
inline Interval<T> operator *(Interval<T> x, const Interval<T> &y) {
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IMul<T>(x, y);
	case DINT_MODE:
		return DIMul<T>(x, y);
	default:
		return IMul<T>(x, y);
	}
}

template<typename T>
inline Interval<T> Interval<T>::operator *(const long double &l) {
	Interval<T> x(this->a, this->b);
	Interval<T> y(static_cast<T>(l), static_cast<T>(l));
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = IMul<T>(x, y);
		break;
	case DINT_MODE:
		r = DIMul<T>(x, y);
		break;
	default:
		r = IMul<T>(x, y);
		break;
	}

	return r;
}

template<typename T>
inline Interval<T> Interval<T>::operator *(const int &i) {
	Interval<T> x(this->a, this->b);
	Interval<T> y(static_cast<T>(i), static_cast<T>(i));
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = IMul<T>(x, y);
		break;
	case DINT_MODE:
		r = DIMul<T>(x, y);
		break;
	default:
		r = IMul<T>(x, y);
		break;
	}

	return r;
}

/**
 * @brief Operator dzielenia, analogicznie z uwzględnieniem trybu w całej klasie.
 */
template<typename T>
inline Interval<T> Interval<T>::operator /(const Interval<T> &y) {
	Interval<T> x(this->a, this->b);
	Interval<T> r = { 0, 0 };
	switch (mode) {
	case PINT_MODE:
		r = IDiv<T>(x, y);
		break;
	case DINT_MODE:
		r = DIDiv<T>(x, y);
		break;
	default:
		r = IDiv<T>(x, y);
		break;
	}

	return r;
}

template<typename T>
inline Interval<T> operator /(Interval<T> x, const Interval<T> &y) {
	switch (Interval<T>::mode) {
	case PINT_MODE:
		return IDiv<T>(x, y);
	case DINT_MODE:
		return DIDiv<T>(x, y);
	default:
		return IDiv<T>(x, y);
	}
}

/**
 * @brief Otoczka wypukła dwóch przedziałów (Hull). 
 * Zwraca najmniejszy pojedynczy przedział, który w całości pokrywa parametry x i y.
 */
template<typename T>
Interval<T> Hull(const Interval<T> &x, const Interval<T> &y) {
	Interval<T> r = { 0, 0 };
	// Szukanie absolutnego minimum lewych krawędzi
	r.a = min(x.a, x.b);
	r.a = min(r.a, y.a);
	r.a = min(r.a, y.b);

	// Szukanie absolutnego maksimum prawych krawędzi
	r.b = max(x.a, x.b);
	r.b = max(r.b, y.a);
	r.b = max(r.b, y.b);

	return r;
}

/**
 * @brief Specjalizacja funkcji formatującej końce przedziału dla wieloprecyzyjnego typu mpreal.
 * Wewnętrznie bazuje bezpośrednio na wskaźnikach i buforach biblioteki MPFR C API.
 */
template<>
inline void Interval<mpreal>::IEndsToStrings(string &left, string &right) {
	mpfr_t rop;
	mpfr_exp_t exponent;
	mpfr_init2(rop, precision);
	char *str = NULL;
	char *buffer = new char(precision + 3);
	mpfr_set(rop, this->a.mpfr_ptr(), MPFR_RNDD); // Pobranie wewnętrznego wskaźnika i zaokrąglenie RNDD
	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDD);
	str = buffer;

	stringstream ss;
	int prec = std::numeric_limits<mpreal>::digits10();
	ss.setf(std::ios_base::scientific);
	bool minus = (str[0] == '-');
	int splitpoint = minus ? 1 : 0;
	string sign = minus ? "-" : "";

	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	left = ss.str();
	ss.str(std::string());

	mpfr_set(rop, this->b.mpfr_ptr(), MPFR_RNDU); // Pobranie wewnętrznego wskaźnika i zaokrąglenie RNDU
	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDU);
	str = buffer;
	splitpoint = (str[0] == '-') ? 1 : 0;
	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	right = ss.str();
	ss.clear();
}

/**
 * @brief Szerokość DINT dedykowana dla mpreal, manipulująca kierunkiem 
 * MPFR RNDU/RNDD wewnątrz własnych klas (mpreal API) zamiast bazowego FPU.
 */
template<>
inline mpreal DIntWidth<mpreal>(const Interval<mpreal> &x) {
	mpreal w1, w2;

	mpreal::set_default_rnd(MPFR_RNDU);
	w1 = x.b - x.a;
	if (w1 < 0)
		w1 = -w1;
	mpreal::set_default_rnd(MPFR_RNDD);
	w2 = x.b - x.a;
	if (w2 < 0)
		w2 = -w2;
	mpreal::set_default_rnd(MPFR_RNDN); // Powrót do TONEAREST wg MPFR
	if (w1 > w2)
		return w1;
	else
		return w2;
}

/**
 * @brief Wartość bezwzględna całego przedziału.
 * Dokonuje translacji tak by zawsze |x| >= 0, zamieniając przy tym ewentualnie krawędzie
 * przedziału jeśli byłyby nieustawione poprawnie pod względem nowej orientacji rosnącej.
 */
template<typename T>
Interval<T> IAbs(const Interval<T> &x) {
	T tmp = 0;
	Interval<T> r = { 0, 0 };
	r.a = abs(x.a);
	r.b = abs(x.b);
	if (r.b < r.a) {
		tmp = r.a;
		r.a = r.b;
		r.b = tmp;
	}

	return r;
}

// ====================================================================================
// ==================== JAWNE INSTANCJONOWANIE SZABLONÓW ==============================
// ====================================================================================
// Niezbędne, aby kompilator poprawnie rozwinął kod dla pożądanych typów liczbowych.

template long double DIntWidth(const Interval<long double> &x);

template class Interval<long double> ;
template class Interval<double> ;
template class Interval<float> ;
template class Interval<mpreal> ;

// Inicjalizacja stałych statycznych
template<typename T> IAMode Interval<T>::mode = PINT_MODE;
template<typename T> IAOutDigits Interval<T>::outdigits = LONGDOUBLE_DIGITS;

template<> IAPrecision Interval<mpreal>::precision = MPREAL_PREC;
template<typename T> IAPrecision Interval<T>::precision = LONGDOUBLE_PREC;

//-------------------------------------------------------------------------------------

} /* namespace interval_arithmetic */

#endif /* INTERVAL_H_ */
