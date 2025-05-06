import time

def format_bytes(size):
    size = float(size)
    for unit in ['o', 'Ko', 'Mo', 'Go', 'To']:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} Po"

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
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=30)
    term1 = 1 / ((1 - X)*(1 - Y))
    termX = (1 + X)^(n-k)
    termXY1 = (1 + X*Y)^k
    termXY2 = (1 + X^2 * Y^2)^m
    grosTerme = (termX * termXY1) / termXY2
    J = term1 * (grosTerme - ((1 + X)^n / (1 + X^2)^m) - ((1 + Y)^k / (1 + Y^2)^m))
    return J

def nb_monomials(d1, d2, k, n):
    term1 = binomial(k, d2 + 1)
    if d2 >= d1 - 1:
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
    and writes the result in fn, sorted by complexity_pre.
    """
    with open(fn, "w") as f:
        # File header
        f.write("d1,d2,value_J,nb_cols,complexity_pre,estimated_footprint\n")
        f.write(f"#n = {n} m = {m}\n")
        sizes = []

        for k in range(k_min, k_max + 1):
            if(k <= 0):
                continue
            print(f"k: {k}")
            complex_exhaustive = ceil(log(n - k)).bit_length() + (n - k)
            p_admi = parametres_admissibles_crossbred(k, m, n)

            f.write(f"\n\nexhaustive search over {n - k} bits costing 2^{complex_exhaustive}:\n")

            data_rows = []

            for (d1, d2), value in p_admi.items():
                if d2 >= d1 - 1:
                    continue
                nb_cols = nb_monomials(d1, d2, k, n)
                if nb_cols == 0:
                    continue
                if nb_cols > (value + m):
                    complexity_pre = nb_cols.bit_length() + 2
                else:
                    complexity_pre = (value + m).bit_length() + 2
                #if complexity_pre >= 128:
                #    continue

                footprint = ceil(nb_cols/8)*(value+m) // 10 // 5

                row = (footprint, d1, d2, value, nb_cols, complexity_pre)
                data_rows.append(row)

            # Sort rows by complexity_pre
            data_rows.sort()
            for i in range(len(data_rows)):
                if data_rows[i][0] > 0:
                    sizes.append((data_rows[i][0], data_rows[i][5], data_rows[i][1], data_rows[i][2], n-k))
                    break

            for footprint, d1, d2, value, nb_cols, complexity_pre in data_rows[:10]:
                f.write(f"{d1},{d2},{value},{nb_cols},{complexity_pre},{format_bytes(footprint)}\n")

        sizes.sort()
        print(f"{format_bytes(sizes[0][0])},{sizes[0][1:]}")
        print(f"{format_bytes(sizes[1][0])},{sizes[1][1:]}")
        print(f"{format_bytes(sizes[2][0])},{sizes[2][1:]}")
        print(f"{format_bytes(sizes[3][0])},{sizes[3][1:]}")
        print(f"{format_bytes(sizes[4][0])},{sizes[4][1:]}")
        print(f"{format_bytes(sizes[5][0])},{sizes[5][1:]}")

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
#parametres_crossbred_sbc(256, 257)

#start = time.time()
#try_parameters_crossbred(101, 100, 10, 90, "parametres_admissibles_crossbred_sbc/101_100.csv")
#end = time.time()
#print(f"{end - start} s")