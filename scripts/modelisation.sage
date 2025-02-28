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

#def eval_g(g, x, y):

x = load("keys/x.sobj")
y = load("keys/y.sobj")

xy = vector(x[:-2].list() + y[:-2].list())

q = 2  # Field characteristic
n = 130  # Vectors' size
k = 257  # Degree of field extension

F_q = GF(q)
F_qn.<t> = GF(q^k, 't')

u = load("keys/u.sobj")
v = load("keys/v.sobj")

assert len(u) == len(v) == len(x) == len(y), "Sizes do not match"

assert u.dot_product(x) * v.dot_product(y) == u.dot_product(y) * v.dot_product(x), "NSBC not right"

print("assert ok")

monomials = ['x'+str(i) for i in range(1, n-1)] + ['y'+str(i) for i in range(1, n-1)]

# Créer l'anneau de polynômes en x1, x2,..., xn-2 et y1, y2,..., yn-2
R = PolynomialRing(F_qn, monomials)

# Convertir les monomiaux en polynômes dans R
monomials = [R(monom) for monom in monomials]

# Construire g
g = construct_g(R, monomials, u, v, n)

print("computed g")

print(f"len(xy): {len(xy)} type: {type(xy)}")

substitutions = {monomials[i]: xy[i] for i in range(k-1)}

evaluation = g.subs(substitutions)
print(evaluation)

coefficients = []

#quotient, reminder = g.quo_rem(t^(k-1))

#print("Décomposition de g dans la base des puissances de t:", quotient)