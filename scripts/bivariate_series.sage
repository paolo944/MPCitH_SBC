import time
import sys
from sage.rings.polynomial.polydict import ETuple

prec = 35

def format_bytes(size):
    size = float(size)
    for unit in ['o', 'Ko', 'Mo', 'Go', 'To']:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} Po"

def print_progress_bar(current, total, bar_length=40):
    percent = float(current) / total
    arrow = '=' * int(round(percent * bar_length) - 1) + '>' if percent > 0 else ''
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f'\rProgress: [{arrow + spaces}] {int(percent * 100)}%')
    sys.stdout.flush()

def H_k_m_n(k, m, n):
    """
    Function to generate the series H in https://ia.cr/2024/992
    To change the precision, change the parameter default_prec
    """
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=prec)
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
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=prec)
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
    R.<X,Y> = PowerSeriesRing(ZZ, default_prec=prec)
    term1 = 1 / ((1 - X)*(1 - Y))
    termX = (1 + X)^(n-k)
    termXY1 = (1 + X*Y)^k
    termXY2 = (1 + X^2 * Y^2)^m
    grosTerme = (termX * termXY1) / termXY2
    J = term1 * (grosTerme - ((1 + X)^n / (1 + X^2)^m) - ((1 + Y)^k / (1 + Y^2)^m))
    return J

def nb_monomials_M_D_d(d1, d2, k, n):
    nb_monomials = 0
    for i in range(d2 + 1, d1 + 1):
        for j in range(0, d1 - i + 1):
            nb_monomials += binomial(k, i) * binomial(n - k, j)
    return nb_monomials

def nb_monomials_Mac_d(d, k, n):
    nb_monomials = 0
    for i in range(d + 1):
        nb_monomials += binomial(n-k, i)
    return nb_monomials

def nb_rows_M_D_d(n, m, k, D, d):
    nb_rows = 0
    for i in range(D+1):
        nb_rows += binomial(n, i)
    for i in range(d-1):
        nb_rows -= binomial(k, i)
    return nb_rows*m

def nb_rows_Mac_d(n, m, k, d):
    nb_rows = 0
    for i in range(d+1):
        nb_rows += binomial(n, i)
    return nb_rows*m

def block_lancszos_complexity(r, W, n, c):
    nnz_c = n*n // 4 + n / 2 + 1
    tmp = max(nnz_c*c^2 / W, c^2)
    return r/W * tmp

def block_wiendemann_complexity(n, mp, np, d, N):
    tmp1 = (1 + np / mp)
    tmp2 = (n*n // 4 + n // 2 + 1 + mp)
    return 2 * tmp1 * N * tmp2

def sparse_factor(n, d1, d2, k):
    return (n*n // 4 + n // 2 + 1) / nb_monomials(d1, d2, k, n)

def get_footprint(nb_rows, nnz):
    return 4 * (nnz + nb_rows + 1)

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
        f.write(f"#n = {n} m = {m}\nGenerating series precision {prec}")
        sizes = []
        fastest_data = []
        fastest_combined_data = []
        fastest_d2 = []

        total_k = k_max - k_min + 1
        progress_counter = 0


        for k in range(k_min, k_max + 1):
            if(k <= 0):
                continue
            progress_counter += 1
            print_progress_bar(progress_counter, total_k)

            #print(f"k: {k}")
            p_admi = parametres_admissibles_crossbred(k, m, n)

            data_rows = []

            g_serie = G_k_m_n(k, m, n)
            g_dict = g_serie.dict()

            for d, value in p_admi.items():
                (d1, d2) = d
                if d2 >= d1 - 1:
                    continue
                if d2 == 0:
                    continue
                nb_cols = nb_monomials_M_D_d(d1, d2, k, n)
                if nb_cols == 0:
                    continue

                h_coeff = g_dict[ETuple(list(d))]
                nb_rows = nb_rows_M_D_d(n, m, k, d1, d2)
                nnz = (n*n // 4 + n // 2 + 1) * (nb_rows)

                complexity_lancszos = float(log(block_lancszos_complexity(256, 64, n, nb_rows), 2))
                complexity_strassen = float(log(max(nb_cols, nb_rows)^2.81, 2))
                complexity_pre = min(complexity_lancszos, complexity_strassen)
                nb_cols_mac_d = nb_monomials_Mac_d(d2, k, n)
                bw1 = block_wiendemann_complexity(n, 512, 512, d1, nb_cols_mac_d)
                strassen = float(log(nb_cols_mac_d^2.81, 2))
                complex_exhaustive = float(log(2^(n-k) * min(bw1, strassen), 2))

                footprint = get_footprint(nb_rows, nnz)

                row = (footprint, d1, d2, nb_rows, nb_cols, complexity_pre)
                data_rows.append(row)

            data_rows.sort()
            for i in range(len(data_rows)):
                if data_rows[i][0] > 0:
                    sizes.append((data_rows[i][0], data_rows[i][5], data_rows[i][1], data_rows[i][2], n-k))
                    fastest_data.append((data_rows[i][5], complex_exhaustive, data_rows[i][1], data_rows[i][2], data_rows[i][0], n-k))
                    if(data_rows[i][2] == 1):
                        fastest_d2.append((data_rows[i][5], complex_exhaustive, data_rows[i][1], data_rows[i][2], data_rows[i][0], n-k))

        f.write("\n\nBest memory footprint:\n")
        f.write("Memory footprint | Memory footprint log_2 | Complexity pre | D | d | n-k\n")
        sizes.sort()
        for data in sizes[:6]:
            if len(data) >= 5:
                f.write(f"{format_bytes(data[0])} | {float(log(data[0], 2))} | {data[1]} | {data[2]} | {data[3]} | {data[4]}\n")
            else:
                print("Skipping data: not enough elements", data)

        f.write("\n\nBest complexity\ncomplexity_pre| complexity_ex| d1 | d2 | Memory footprint | Memory footprint log_2 | n-k:\n")
        fastest_data.sort()
        for data in fastest_data[:6]:
            if len(data) >= 5:
                f.write(f"{data[0]} | {data[1]} | {data[2]} | {data[3]} | {format_bytes(data[4])} | {float(log(data[4], 2))} | {data[5]}\n")
            else:
                print("Skipping data: not enough elements", data)

        f.write("\n\nBest max complexity\ncomplexity_pre | complexity_ex| d1 | d2 | Memory footprint | Memory footprint log_2 | n-k:\n")
        fastest_data.sort(key=lambda x: max(x[0], x[1]))
        for data in fastest_data[:6]:
            if len(data) >= 5:
                f.write(f"{data[0]} | {data[1]} | {data[2]} | {data[3]} | {format_bytes(data[4])} | {float(log(data[4], 2))}| {data[5]}\n")
            else:
                print("Skipping data: not enough elements", data)
        
        f.write("\n\nBest max complexity for d==1\ncomplexity_pre | complexity_ex| d1 | d2 | Memory footprint | Memory footprint log_2 | n-k:\n")
        fastest_d2.sort(key=lambda x: max(x[0], x[1]))
        for data in fastest_d2[:6]:
            if len(data) >= 5:
                f.write(f"{data[0]} | {data[1]} | {data[2]} | {data[3]} | {format_bytes(data[4])} | {float(log(data[4], 2))}| {data[5]}\n")
            else:
                print("Skipping data: not enough elements", data)
    print()

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
        if n - 140 > 0:
            k_min =  n - 140
        else:
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

start = time.time()
try_parameters_crossbred(257, 256, 128, 240, "parametres_admissibles_crossbred_sbc/257_256_256.csv")
#parametres_crossbred_sbc(100, 256)
end = time.time()
print(f"{end - start} s")