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

def Nm(n, m, t1, t2):
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
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) +2)
    denom = ((1 - tx)^(nx)) * ((1 - ty)^(ny))
    num = (1 - tx*ty)^m + Nm(ny, m, tx, ty) + Nm(nx, m, ty, tx)
    return num / denom

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

series_ring.<z> = PowerSeriesRing(ZZ)

n = system[0].parent().ngens()
#
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
#
term1 = 1
for i in system2:
    term1 *= (1-z**i.degree())

term2 = (1-z)**n
gen_serie = term1 / term2
#n_neg = next(i for i in range(gen_serie.prec()) if gen_serie[i] <= 0)
#
bi_reg = hilbert_biseries(4, 4, 9)
bi_reg_uni = convert_bi_series(bi_reg)
hilbert_series = series_ring(I.hilbert_series())
print(f"The hilbert series of I is: {hilbert_series}\n")
print(f"Generating series of the sequence: {gen_serie}\n")
#print(f"First degree of non positif term : {n_neg}\n",)
print(f"Série bi-regulière: {bi_reg}\n")
print(f"Série bi-regulière convertie: {bi_reg_uni}\n")
print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")

#equal_up_to = all(gen_serie[i] == hilbert_series[i] for i in range(n_neg))
#print(f"Is the sequence semi-regular ? {"yes" if equal_up_to else "no"}")


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