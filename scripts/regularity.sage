def doit(n, m):
    R = PolynomialRing(GF(2), n, 'x')
    
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

if(sys.argv[1] == "random"):
    try:
        system = doit(int(sys.argv[2]), int(sys.argv[3]))
    except Exception as error:
        print("Erreur pendant la génération du système aléatoire: ", error)
        sys.exit(1)
else:
    try:
        system = load(sys.argv[1])
    except Exception as error:
        print("Erreur pendant le load: ", error)
        sys.exit(1)

system2 = []
for i in system:
    try:
        system2.append(i.homogeneous_components()[2])
    except KeyError:
        system2.append(i.homogeneous_components()[1])
#print(system2)
I = ideal(system2)

field = ZZ
poly_ring.<t> = PolynomialRing(field)
series_ring.<z> = PowerSeriesRing(field)

n = system[0].parent().ngens()

denom = 1
for i in system2:
    denom *= (1-t**i.degree())

denom = series_ring(denom)
num = series_ring((1-t)**n)
res = num.inverse_of_unit() * denom
print(f"Generating series of the sequence: {res}")

hilbert_series = I.hilbert_series()
hilbert_series_denom = hilbert_series.denominator()
hilbert_series_num = hilbert_series.numerator()
res2 = series_ring(hilbert_series_num) * series_ring(hilbert_series_denom).inverse_of_unit()
print(f"The hilbert series of I is: {res2}")
print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")
