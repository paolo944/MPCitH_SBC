import time

def read_privkey(f):
    vec = []
    B = GF(2)
    
    with open(f, 'r') as f:
        for ligne in f:
            elements = ligne.split()
            taille = int(elements[0])
            caracteristique = int(elements[1])
            if caracteristique != 2:
                raise ValueError(f"Caractéristique attendue 2, mais trouvée {caracteristique} dans la ligne : {ligne}")
            
            coefficients = [int(c) for c in elements[2:]]
            elem = sum([coefficients[i] for i in range(len(coefficients))])
            vec.append(B(elem))
    
    return vector(vec)

def read_pubkey(f):
    polynomes = []
    B.<t> = GF(2^257, modulus=x^257 + x^12 + 1)
    
    with open(f, 'r') as f:
        for ligne in f:
            elements = ligne.split()
            taille = int(elements[0])
            caracteristique = int(elements[1])
            if caracteristique != 2:
                raise ValueError(f"Caractéristique attendue 2, mais trouvée {caracteristique} dans la ligne : {ligne}")
            
            coefficients = [int(c) for c in elements[2:]]
            poly = sum([coefficients[i] * t^i for i in range(len(coefficients))])
            polynomes.append(B(poly))
    
    return (vector(polynomes), B)


def construct_g(R, monomials, u, v, n):
    ux = R(u[n-2])
    vx = R(v[n-2])
    uy = R(u[n-1])
    vy = R(v[n-1])

    for i in range(n-2):
        ux += R(u[i])*monomials[i]
        vx += R(v[i])*monomials[i]
        uy += R(u[i])*monomials[i+n-2]
        vy += R(v[i])*monomials[i+n-2]

    monomes = set()
    monomes.update((ux*vy).monomials())
    monomes.update((uy*vx).monomials())

    print(f"nb monomes {len(monomes)}") 

    return ux*vy - uy*vx

def decompose_g(R, g, basis, monomials_str, k):
    p = R.characteristic()
    Rprime = PolynomialRing(GF(p), monomials_str)
    system = [Rprime(0) for i in range(k)]
    monomials = g.monomials()
    g_coeffs = [g.monomial_coefficient(i) for i in monomials]
    for i in range(len(g_coeffs)):
        tmp_coeffs = g_coeffs[i].list()
        for j in range(len(tmp_coeffs)):
            if(tmp_coeffs[j] == 1):
                system[j] += Rprime(monomials[i])
    return system

def generate_system(n, field_eq_op, verbose):
    if n < 1:
        print("the vector's size must be positive")
        sys.exit(1)
    
    if verbose:
        print(f"generating system for n={n}")

    x = load("keys/x.sobj")
    y = load("keys/y.sobj")
    #x = read_privkey("keys/x")
    #y = read_privkey("keys/y")

    u = load("keys/u.sobj")
    v = load("keys/v.sobj")
    #u, F_qn = read_pubkey("keys/u.pub")
    #v, _ = read_pubkey("keys/v.pub")

    assert x[n-1] == 0 and x[n-2] == 1 and y[n-1] == 1 and y[n-2] == 0, "NSBC not right"
    assert len(u) == len(v) == len(x) == len(y), "Sizes do not match"
    expr = u.dot_product(x) * v.dot_product(y) - u.dot_product(y) * v.dot_product(x)  
    assert expr == 0, "Keys does not match"

    xy = vector(x[:-2].list() + y[:-2].list())

    q = 2  # Field characteristic
    k = 2*(n-2)+1 # Degree of field extension
    nvars = k-1

    F_q = GF(q)
    F_qn = GF(q^k, 't')

    if verbose:
        start_time = time.time()

    if verbose:
        print("assert ok")

    monomials_str = ['x'+str(i) for i in range(1, n-1)] + ['y'+str(i) for i in range(1, n-1)]

    # Créer l'anneau de polynômes en x1, x2,..., xn-2 et y1, y2,..., yn-2
    R = PolynomialRing(F_qn, monomials_str)
    Rp = PolynomialRing(F_q, monomials_str)

    # Convertir les monomiaux en polynômes dans R
    monomials = [R(monom) for monom in monomials_str]
    
    # Construire g
    g = construct_g(R, monomials, u, v, n)

    print(f"n_x = {nvars // 2} n_y = {nvars // 2} m = {k}")
    print(f"number of monomials in g: {len(g.monomials())}")

    if verbose:
        print("computed g")

    substitutions = {monomials[i]: xy[i] for i in range(k-1)}

    evaluation = g.subs(substitutions)
    if verbose:
        print(f"g(x, y) == {evaluation}")

    coefficients = []

    #quotient, reminder = g.quo_rem(t**2)

    basis = [F_qn.gen()**i for i in range(k)]

    system = decompose_g(R, g, basis, monomials_str, k)

    monomials = [Rp(monom) for monom in monomials_str]

    if(field_eq_op):
        if verbose:
            print("ajout des équations du corps")
        for i in range(nvars):
            system.append(monomials[i]**2 - monomials[i])

    t = k
    if(field_eq_op):
        t += nvars

    file = f"system/sage/system_bilin_{nvars}_{t}.sobj"
    save(system, file)
    if verbose:
        print("system computed")
        end_time = time.time()
        print(f"Temps d'exécution : {end_time - start_time} secondes")
        print(f"system written in {file}")

    return (system, file)


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("Usage: sage modelisation.py <n> <field_eq (0 or 1)>")
        sys.exit(1)

    n = int(sys.argv[1])
    field_eq = bool(int(sys.argv[2]))
    generate_system(n, field_eq, True)