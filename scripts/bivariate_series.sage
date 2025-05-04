def H_k_m_n(k, m, n):
    """
    Function to generate the series H in https://ia.cr/2024/992
    To change the precision, change the parameter default_prec
    """
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=15)
    termX = (1+X)^(n-k)
    termXY1 = (1+X*Y)^k
    termXY2 = (1+X^2 * Y^2)^m
    H = (termX - (termXY1 * termX) / termXY2)/Y
    return H

def G_k_m_n(k, m, n):
    """
    Function to generate the series G in https://ia.cr/2024/992
    To change the precision, change the parameter default_prec
    """
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=15)
    termX = (1+X)^(n-k)
    termXY1 = (1+X)^k
    termXY2 = (1+X^2)^m
    H_1 = (termX - (termXY1 * termX) / termXY2)
    G = -(Y*H_k_m_n(k, m, n) - H_1)/((1-X)*(1-Y))
    return G

def J_k_m_n(k, m, n):
    """
    Function to generate the series J in https://ia.cr/2024/992
    To change the precision, change the parameter default_prec
    """
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=32)
    term1 = 1 / ((1 - X)*(1 - Y))
    termX = (1 + X)^(n-k)
    termXY1 = (1 + X*Y)^k
    termXY2 = (1 + X^2 * Y^2)^m
    grosTerme = (termX * termXY1) / termXY2
    J = term1 * (grosTerme - ((1 + X)^n / (1 + X^2)^m) - ((1 + Y)^k / (1 + Y^2)^m))
    return J

def nb_monomials(d1, d2, k, n):
    term1 = binomial(k, d2 + 1)
    if d2 > d1 - 1:
        return 0
    term2 = binomial((n - k), d1 - d2 - 1)
    return term1 * term2

def parametres_admissibles(k, m, n):
    parametres = G_k_m_n(k, m, n).monomial_coefficients()
    parametres_admissibles = {}
    for (d1, d2), value in parametres.items():
        if value > 0:
            if d2 <= d1:
                parametres_admissibles[(d1, d2)] = value
    return parametres_admissibles

def parametres_admissibles_crossbred(k, m, n):
    """
    Returns the admissible parameters for the crossbred
    algorithm as in https://ia.cr/2024/992
    """
    parametres = J_k_m_n(k, m, n).monomial_coefficients()
    parametres_admissibles = {}
    for (d1, d2), value in parametres.items():
        if value > 0:
            if d2 <= d1:
                parametres_admissibles[(d1, d2)] = value
    return parametres_admissibles

def try_parameters(m, n, k_min, k_max, fn):
    #file header
    with open(fn, "w") as f:
        f.write("d1,d2,k,m,n,nouveaux polynômes pre\n")

        for k in range(k_min, k_max+1):
            print(f"k: {k}")
            p_admi = parametres_admissibles(k, m, n)
            for (d1, d2), value in p_admi.items():
                f.write(f"{d1},{d2},{k},{m},{n},{value}\n")

def try_parameters_crossbred(m, n, k_min, k_max, fn):
    """
    Tries multiple parameters with different k_min <= k <= k_max
    and writes the result in fn
    """
    #file header
    with open(fn, "w") as f:
        f.write("d1,d2,nb_new_poly,nb_cols,complexity_pre\n")
        f.write(f"#n = {n} m = {m}\n")
        for k in range(k_min, k_max+1):
            print(f"k: {k}")
            p_admi = parametres_admissibles_crossbred(k, m, n)
            f.write(f"\n\nexhaustive search over {n - k} bits:\n")
            for (d1, d2), value in p_admi.items():
                nb_cols = nb_monomials(d1, d2, k, n)
                complexity_pre = nb_cols.nbits() + 2
                f.write(f"{d1},{d2},{value},{nb_cols},{complexity_pre}\n")

def parametres_crossbred_large(min_n, max_n):
    for n in range(min_n, max_n + 1, 2):
        print(f"n: {n}")
        m = n + 1
        k_min = n - 90
        k_max = n - 30
        fn = f"parametres_admissibles_crossbred/{m}_{n}.csv"
        try_parameters_crossbred(m, n, k_min, k_max, fn)

def parametres_crossbred_sbc(min_n, max_n):
    """
    Tries multiple valid parameters for SBC
    where m is prime and m = n+1
    """
    m = next_prime(min_n) 
    m_max = next_prime(max_n + 1)
    while m <= m_max:
        n = m - 1
        if n%2 != 0:
            m = next_prime(m + 1)
            continue
        print(f"n: {n} m: {m}")
        k_min = n - 90
        k_max = n - 30
        fn = f"parametres_admissibles_crossbred_sbc/{m}_{n}.csv"
        try_parameters_crossbred(m, n, k_min, k_max, fn)
        m = next_prime(m + 1)

#print(H_k_m_n(24, 160, 80))
#print(f"Série génératrice pour paramètres admissibles: \n{G_k_m_n(24, 160, 80).monomial_coefficients()}")
#print(f"Paramètres admissibles pour k=24, m=160, n=80:\n{parametres_admissibles(128, 257, 256)}")
#try_parameters(257, 256, 64, 128, "parametres_admissibles_256_257.csv")
#print(f"J(24, 160, 80): {J_k_m_n(24, 160, 80)}")
#print(f"Parametres admissibles CrossBred: {parametres_admissibles_crossbred(24, 160, 80)}")
#try_parameters_crossbred(257, 256, 64, 128, "parametres_admissibles_crossbred_256_257.csv")
#parametres_crossbred_large(10, 256)
parametres_crossbred_sbc(256, 257)