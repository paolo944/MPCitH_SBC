def H_k_m_n(k, m, n):
    R.<X,Y> = PowerSeriesRing(ZZ)
    termX = (1+X)^(n-k)
    termXY1 = (1+X*Y)^k
    termXY2 = (1+X^2 * Y^2)^m
    H = (termX - (termXY1 * termX) / termXY2)/Y
    return H

def G_k_m_n(k, m, n):
    R.<X,Y> = PowerSeriesRing(ZZ)
    termX = (1+X)^(n-k)
    termXY1 = (1+X)^k
    termXY2 = (1+X^2)^m
    H_1 = (termX - (termXY1 * termX) / termXY2)
    G = -(Y*H_k_m_n(k, m, n) - H_1)/((1-X)*(1-Y))
    return G

def parametres_admissibles(k, m, n):
    parametres = G_k_m_n(k, m, n).monomial_coefficients()
    parametres_admissibles = {}
    for (d1, d2), value in parametres.items():
        if value > 0:
            if d2 <= d1:
                parametres_admissibles[(d1, d2)] = value
    return parametres_admissibles

def J_k_m_n(k, m, n):
    R.<X,Y> = PowerSeriesRing(ZZ)
    term1 = 1 / ((1 - X)*(1 - Y))
    termX = (1 + X)^(n-k)
    termXY1 = (1 + X*Y)^k
    termXY2 = (1 + X^2 * Y^2)^m
    grosTerme = (termX * termXY1) / termXY2
    J = term1 * (grosTerme - ((1 + X)^n / (1 + X^2)^m) - ((1 + Y)^k / (1 + Y^2)^m)) + R.O(10)
    return J

def parametres_admissibles_crossbred(k, m, n):
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
    #file header
    with open(fn, "w") as f:
        f.write("d1,d2,k,m,n,coeff\n")

        for k in range(k_min, k_max+1):
            #print(f"k: {k}")
            p_admi = parametres_admissibles_crossbred(k, m, n)
            for (d1, d2), value in p_admi.items():
                f.write(f"{d1},{d2},{k},{m},{n},{value}\n")

def parametres_crossbred_large(min_n, max_n):
    for n in range(min_n, max_n + 1, 2):
        print(f"n: {n}")
        m = n + 1
        k_min = n // 4
        k_max = n // 2
        fn = f"parametres_admissibles_crossbred/{m}_{n}.csv"
        try_parameters_crossbred(m, n, k_min, k_max, fn)

def parametres_crossbred_sbc(min_n, max_n):
    m = next_prime(min_n) 
    m_max = next_prime(max_n + 1)
    while m <= m_max:
        n = m - 1
        if n%2 != 0:
            m = next_prime(m + 1)
            continue
        print(f"n: {n} m: {m}")
        k_min = n // 4
        k_max = n // 2
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
parametres_crossbred_sbc(10, 300)