markdown_content = r"""# Obliczanie pierwiastków zespolonych wielomianu o współczynnikach rzeczywistych metodą Bairstowa

**Autor:** Adam Kwiatkowski  

---

## 1 Zastosowanie
Procedura `bairstow_method` oblicza pierwiastki zespolone wielomianu (1) o współczynnikach rzeczywistych w zmiennopozycyjnej arytmetyce przedziałowej

$$p(x) = a_nx^n + a_{n-1}x^{n-1} + ... + a_1x + a_0 \tag{1}$$

## 2 Opis metody
Do znalezienia pierwiastków zespolonych wielomianu (1) stosuje się metodę Bairstowa, w której unika się rachunków na liczbach zespolonych. Metoda ta opiera się na twierdzeniu, że pierwiastki wielomianu

$$x^2 - rx - q \tag{2}$$

gdzie $r$ i $q$ oznaczają współczynniki rzeczywiste, są pierwiastkami wielomianu (1), gdy wielomian (1) dzieli się bez reszty przez wielomian (2). W ogólności wykonując dzielenie wielomianu (1) przez wielomian (2) mamy

$$p(x) = p_1(x)(x^2 - rx - q) + Ax + B \tag{3}$$

przy czym stopień wielomianu $p_1(x)$ jest nie większy niż $n - 2$, a wyrażenie $Ax + B$ jest resztą z dzielenia. Współczynniki $A$ i $B$ zależą od $r$ i $q$:
$A = A(r, q), B = B(r, q)$

Reszta z dzielenia będzie równa 0, jeżeli
$A(r, q) = 0, B(r, q) = 0$

Powyższe równania są układem dwu równań nieliniowych o dwu niewiadomych $r$ i $q$. Do ich rozwiązania możemy zastosować proces iteracyjny Newtona:

$$
\begin{bmatrix} r_{j+1} \\ q_{j+1} \end{bmatrix} = \begin{bmatrix} r_j \\ q_j \end{bmatrix} - \begin{bmatrix} \frac{\partial A}{\partial r} & \frac{\partial A}{\partial q} \\ \frac{\partial B}{\partial r} & \frac{\partial B}{\partial q} \end{bmatrix}^{-1}_{r=r_j, q=q_j} \begin{bmatrix} A(r_j, q_j) \\ B(r_j, q_j) \end{bmatrix}
$$

który można zapisać w postaci układu równań liniowych:

$$
\begin{bmatrix} \frac{\partial A}{\partial r} & \frac{\partial A}{\partial q} \\ \frac{\partial B}{\partial r} & \frac{\partial B}{\partial q} \end{bmatrix}_{r=r_j, q=q_j} \begin{bmatrix} r_{j+1} \\ q_{j+1} \end{bmatrix} = \begin{bmatrix} \frac{\partial A}{\partial r} & \frac{\partial A}{\partial q} \\ \frac{\partial B}{\partial r} & \frac{\partial B}{\partial q} \end{bmatrix}_{r=r_j, q=q_j} \begin{bmatrix} r_j \\ q_j \end{bmatrix} - \begin{bmatrix} A(r_j, q_j) \\ B(r_j, q_j) \end{bmatrix} \tag{4}
$$

Wartości $A(r, q)$ i $B(r, q)$ można otrzymać za pomocą następujących wzorów rekurencyjnych:
$$b_0 = a_n, b_1 = b_0r + a_{n-1},$$
$$b_i = b_{i-2}q + b_{i-1}r + a_{n-i}, i = 2, 3, \dots , n - 2,$$
$$A(r, q) = b_{n-3}q + b_{n-2}r + a_1, B(r, q) = b_{n-2}q + a_0.$$

Pochodne cząstkowe występujące w równaniach (4) wyznacza się na podstawie wzorów
$$\frac{\partial A(r, q)}{\partial r} = r\bar{A}(r, q) + \bar{B}(r, q),$$
$$\frac{\partial A(r, q)}{\partial q} = \bar{A}(r, q),$$
$$\frac{\partial B(r, q)}{\partial r} = q\bar{A}(r, q),$$
$$\frac{\partial B(r, q)}{\partial q} = \bar{B}(r, q),$$

gdzie wartości $\bar{A}(r, q)$ i $\bar{B}(r, q)$ mogą być obliczone na podstawie podobnych wzorów rekurencyjnych, a mianowicie
$$c_0 = b_0, c_1 = c_0r + b_1,$$
$$c_i = c_{i-2}q + c_{i-1}r + b_i, i = 2, 3, \dots , n - 4,$$
$$\bar{A}(r, q) = c_{n-5}q + c_{n-4}r + b_{n-3}, \bar{B}(r, q) = c_{n-4}q + b_{n-2}.$$

## 3 Wywołanie funkcji
```cpp
vector<T> bairstow_method(
    degree,
    coefficients,
    max_iterations,
    relative_error,
    double zerodet,
    st,
    it,
    input_type
);
