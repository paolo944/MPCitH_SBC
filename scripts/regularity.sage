def doit(n, m):
    """
    Generate random system of n variables and m
    polynomials on GF(2)
    Stolen from hpXbred :)
    """
    # planted solution
    V = GF(2)**n 
    x = V.random_element() 
    I = [] 
    
    def random_quad_poly(R):
        K = R.base_ring()
        v = vector(R.gens())
        n = len(v) 
        Mq = matrix.random(K, n, n)
        Ml = matrix.random(K, 1, n)
        f = v * Mq * v + (Ml*v)[0] + K.random_element()
        return f
    
    # m random polynomials
    for _ in range(m): 
        f = random_quad_poly(R) 
        f += f(*x) 
        I.append(f)

    return I

def hilbert_function(n, mpol, deg):
    """
    En supossant que la suite est semi-regulière
    """
    if(degree < 2 or mpol == 0):
        return binomial(nvar,deg)

def Nm(n, m, t1, t2):
    """
    internal function for the hilbert_biseries
    """
    sum_l = 0
    for l in range(1, m - n + 1):
        sum_k = 0
        for k in range(1, n + 1):
            sum_k += (t1^(n - k))*binomial(l + n - k - 1, n - k)
        bracket = 1 - (1 - t1)**l * sum_k
        term = ((1 - t1*t2)^(m - n - l)) * t1*t2*((1 - t2)^(n))
        sum_l += term * bracket
    return sum_l

def hilbert_biseries(nx, ny, m):
    """
    Returns the hilbert bi-series for a quadratic
    system of nx + ny variables and m polynomials
    result from https://arxiv.org/abs/1001.4004
    """
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) +2)
    denom = ((1 - tx)^(nx)) * ((1 - ty)^(ny))
    num = (1 - tx*ty)^m + Nm(ny, m, tx, ty) + Nm(nx, m, ty, tx)
    return num / denom

def convert_bi_series(Hs):
    """
    Convert the hilbert_biseries with tx=ty
    into a univariate power series
    """
    coefficients = Hs.monomial_coefficients()
    prec = Hs.prec()
    R.<z> = PowerSeriesRing(ZZ, default_prec=prec)
    H = R(0)
    for (t1, t2), coeff in coefficients.items():
        H += coeff*z^(t1 + t2)
    return H

def write_homogenized_system_magma(fn, system):
    """
    Writes the homogenized system in a magma type file
    """
    n = system[0].parent().ngens()
    magma_str = "F := GaloisField(2);\n"
    variables = ', '.join(f'x{i+1}' for i in range(n/2))
    variables += ', ' + ', '.join(f'y{i+1}' for i in range(n/2))
    magma_str += f"Field<{variables}> := PolynomialRing(F, {n}, \"grevlex\");\n"
    
    for idx, f in enumerate(system, 1):
        magma_str += f"f{idx} := {f};\n"
    
    f_list = ', '.join(f"f{i+1}" for i in range(len(system)))
    magma_str += f"PolynomialSystem := [{f_list}];\n"
    
    # Sauvegarde dans un fichier texte
    with open(fn, "w") as fd:
        fd.write(magma_str)

def homogenized_ideal(system):
    """
    Returns the homogenized system that is supposed
    quadratic
    """
    system2 = []
    for i in system:
        try:
            system2.append(i.homogeneous_components()[2])
        except KeyError:
            system2.append(i.homogeneous_components()[1])

    return system2

def generating_bardet_series(system):
    """
    Returns the generating series of bardet of a supposed 
    semi-regular system without taking into account the 
    field equations
    """
    n = system[0].parent().ngens()
    term1 = 1
    for i in system:
        term1 *= (1-z**i.degree())

    term2 = (1-z)**n
    return term1 / term2

def print_help():
    print("Two use examples of the script:")
    print("\t$sage regularity.sage random n m")
    print("\t$sage regularity.sage system_file.sobj")

if __name__ == "__main__":
    import sys

    if(len(sys.argv) == 4 and sys.argv[1] == "random"):
        try:
            R = PolynomialRing(GF(2), int(sys.argv[2]), 'x')
            system = doit(int(sys.argv[2]), int(sys.argv[3]), R)
        except Exception as error:
            print("Error while generating the random system: ", error)
            print_help()
            sys.exit()
    elif(len(sys.argv) == 2):
        try:
            system = load(sys.argv[1])
        except Exception as error:
            print("Erreur during the loading of the system: ", error)
            print_help()
            sys.exit()
    else:
        print_help()
        sys.exit()

    system_homo = homogenized_ideal(system)
    I = ideal(system_homo)
    series_ring.<z> = PowerSeriesRing(ZZ)
    
    gen_serie = generating_bardet_series(system_homo)

    #To find the first non_postivie coefficient of the serie
    #n_neg = next(i for i in range(gen_serie.prec()) if gen_serie[i] <= 0)
    
    n = system[0].parent().ngens()

    bi_reg = hilbert_biseries(n//2, n//2, len(system))
    bi_reg_uni = convert_bi_series(bi_reg)
    hilbert_series = series_ring(I.hilbert_series())

    print(f"Series for a polynomial system defined over {system[0].base_ring()} with {n} variables and {len(system)} polynomials")

    print(f"The hilbert series of I computed by sage (using the grobner basis) is: {hilbert_series}\n")
    print(f"Bardet generating series of the sequence: {gen_serie}\n")
    #print(f"First degree of non positif term : {n_neg}\n",)
    print(f"Hilbert bi-series: {bi_reg}\n")
    print(f"Converted bi-series with tx=ty: {bi_reg_uni}\n")
    print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")

    write_homogenized_system_magma("system_homo.magma", system_homo)

    #equal_up_to = all(gen_serie[i] == hilbert_series[i] for i in range(n_neg))
    #print(f"Is the sequence semi-regular ? {"yes" if equal_up_to else "no"}")


    #This section contains some code about the Hilbert bi-series of 128 128 257
    """
    neg_coeffs = load("scripts/neg_coeffs_128_128_257.sobj")
    for i, j in neg_coeffs.items():
        print(f'{i}: {j}')
    #
    min_deg_deg = 130
    min_deg_monomial = 0
    min_deg_coeff = 0

    min_coeff_coeff = -2895640507456856514805834469828912633846067050403218588917176002391601569024
    min_coeff_monomial = 0
    min_coeff_deg = 0

    for monomial, coeff in neg_coeffs.items():
        if monomial.degree() < min_deg_deg:
            min_deg_deg = monomial.degree()
            min_deg_monomial = monomial
            min_deg_coeff = coeff
        if coeff > min_coeff_coeff:
            min_coeff_coeff = coeff
            min_coeff_monomial = monomial
            min_coeff_deg = monomial.degree()

    print(min_deg_deg)
    print(min_deg_monomial)
    print(min_deg_coeff)
    print()
    print(min_coeff_deg)
    print(min_coeff_monomial)
    print(min_coeff_coeff)

    serie = load("scripts/bi_serie_128_128_257.sobj")
    print(serie)
    """