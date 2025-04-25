import time

def construct_g(R, monomials, u, v, n):
    ux = R(u[n-2])
    vx = R(v[n-2])
    uy = R(u[n-1])
    vy = R(v[n-1])

    for i in range(n-2):
        ux += u[i]*monomials[i]
        vx += v[i]*monomials[i]
        uy += u[i]*monomials[i+n-2]
        vy += v[i]*monomials[i+n-2]

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
    if verbose:
        print(f"generating system for n={n}")

    x = load("keys/x.sobj")
    y = load("keys/y.sobj")

    xy = vector(x[:-2].list() + y[:-2].list())

    q = 2  # Field characteristic
    k = 2*(n-2)+1 # Degree of field extension
    nvars = k-1

    if n < 1:
        print("the vector's size must be positive")
        sys.exit(1)

    F_q = GF(q)
    F_qn.<t> = GF(q^k, 't')

    u = load("keys/u.sobj")
    v = load("keys/v.sobj")

    if verbose:
        start_time = time.time()

    assert len(u) == len(v) == len(x) == len(y), "Sizes do not match"

    assert u.dot_product(x) * v.dot_product(y) == u.dot_product(y) * v.dot_product(x), "NSBC not right"

    if verbose:
        print("assert ok")

    monomials_str = ['x'+str(i) for i in range(1, n-1)] + ['y'+str(i) for i in range(1, n-1)]

    # Créer l'anneau de polynômes en x1, x2,..., xn-2 et y1, y2,..., yn-2
    R = PolynomialRing(F_qn, monomials_str)

    # Convertir les monomiaux en polynômes dans R
    monomials = [R(monom) for monom in monomials_str]
    
    # Construire g
    g = construct_g(R, monomials, u, v, n)

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


#if __name__ == '__main__':
#    import sys
#    if len(sys.argv) != 3:
#        print("Usage: sage modelisation.py <n> <field_eq (0 or 1)>")
#        sys.exit(1)
#
#    n = int(sys.argv[1])
#    field_eq = bool(int(sys.argv[2]))
#    generate_system(n, field_eq)