def doit(n, m, R):    
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

def Md(n):
    """
    Retourne la série qui représente le nombre de mônomes en degré D
    """
    poly_ring.<z> = PolynomialRing(QQ)
    return (1 + z)**n

def hilbert_function(n, mpol, deg):
    """
    En supossant que la suite est semi-regulière
    """
    if(degree < 2 or mpol == 0):
        return binomial(nvar,deg)

def homogenize(system, R):
    system_homo = []
    for poly in system:
        poly_homo = R(0)
        for monomial in poly.monomials():
            if monomial.deg() == 2:
                poly_homo += monomial
        system_homo.append(poly_homo)
    return system_homo

def hilbert_biseries(n, m):
    R.<t1,t2> = PowerSeriesRing(ZZ)
    denom = (1 - t1)^(n + 1) * (1 - t2)^(n+1)
    Nm = (1 - t1*t2)^m
    S1 = 0
    S2 = 0
    for l in range(1, m - (n + 1) + 1):
        term1 = (1 - t1*t2)^(m - (n + 1) - l) * t1*t2*(1-t2)^(n + 1)
        term2 = (1 - t1*t2)^(m - (n + 1) - l) * t1*t2*(1-t1)^(n + 1)

        Brackets1 = 1 - (1 - t1)^l
        Brackets2 = 1 - (1 - t1)^l
        S1_in = 0
        S2_in = 0

        for k in range(1, n + 1):
            S1_in += t1^(n + 1 - k)*binomial(l + n - k, n + 1 - k)
            S1_in += t2^(n + 1 - k)*binomial(l + n - k, n + 1 - k)
        
        Brackets1 *= S1_in
        Brackets2 *= S2_in

        S1 += term1 * S1_in
        S2 += term2 * S2_in
    Nm += S1 + S2
    return Nm / denom

def convert_bi_series(Hs):
    coefficients = Hs.monomial_coefficients()
    prec = Hs.prec()
    R.<z> = PowerSeriesRing(ZZ, default_prec=prec)
    H = R(0)
    for (t1, t2), coeff in coefficients.items():
        H += coeff*z^(t1 + t2)
    return H

if(sys.argv[1] == "random"):
    try:
        R = PolynomialRing(GF(2), int(sys.argv[2]), 'x')
        system = doit(int(sys.argv[2]), int(sys.argv[3]), R)
    except Exception as error:
        print("Erreur pendant la génération du système aléatoire: ", error)
        sys.exit(1)
else:
    try:
        system = load(sys.argv[1])
    except Exception as error:
        print("Erreur pendant le load: ", error)
        sys.exit(1)

#print(system)
#system2 = homogenize(system, R)
#print(system2)

system2 = []
for i in system:
    try:
        system2.append(i.homogeneous_components()[2])
    except KeyError:
        system2.append(i.homogeneous_components()[1])

I = ideal(system2)

field = ZZ
poly_ring.<t> = PolynomialRing(field)
series_ring.<z> = PowerSeriesRing(field)

n = system[0].parent().ngens()

#magma_str = "F := GaloisField(2);\n"
#variables = ', '.join(f'x{i+1}' for i in range(n/2))
#variables += ', ' + ', '.join(f'y{i+1}' for i in range(n/2))
#magma_str += f"Field<{variables}> := PolynomialRing(F, {n}, \"grevlex\");\n"
#
#for idx, f in enumerate(system2, 1):
#    magma_str += f"f{idx} := {f};\n"
#
#f_list = ', '.join(f"f{i+1}" for i in range(len(system2)))
#magma_str += f"PolynomialSystem := [{f_list}];\n"
#
#print(sys.argv[-1])
#
## Sauvegarde dans un fichier texte
#with open(sys.argv[-1], "w") as f:
#    f.write(magma_str)

denom = 1
for i in system2:
    denom *= (1-t**i.degree())

denom = series_ring(denom)
num = series_ring((1-t)**n)
res = num.inverse_of_unit() * denom
print(f"Generating series of the sequence: {res}")

bi_reg = hilbert_biseries(3, 7)
bi_reg_uni = convert_bi_series(bi_reg)
hilbert_series = I.hilbert_series()
hilbert_series_denom = hilbert_series.denominator()
hilbert_series_num = hilbert_series.numerator()
res2 = series_ring(hilbert_series_num) * series_ring(hilbert_series_denom).inverse_of_unit()
print(f"The hilbert series of I is: {res2}")
print(f"Série bi-regulière: {bi_reg}")
print(f"Série bi-regulière convertie: {bi_reg_uni}")
print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")
