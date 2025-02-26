def load_vector(file_name):
    with open(file_name, "r") as f:
        variables_line = f.readline().strip()
        variables = tuple(map(str, variables_line.split(',')))

        characteristic_line = f.readline().strip()
        characteristic = int(characteristic_line)
        
        vector = []
        for line in f:
            line = line.strip().strip(',\n').replace('^', '**')
            vector.append(eval(line))
        
        F = GF(characteristic, len(variables))
        loaded_vector = vector
    
    return loaded_vector, F

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

x, _ = load_vector("x")
y, _ = load_vector("y")

xy = x[:-2] + y[:-2]

q = 2  # Field characteristic
n = 130  # Vectors' size
k = 257  # Degree of field extension

F_q = GF(q)
F_qn.<t> = GF(q^k, 't')

generator = F_qn.gen()

basis = [generator^i for i in range(k)]

u, _ = load_vector("u.pub")
v, _ = load_vector("v.pub")

monomials = ['x'+str(i) for i in range(1, n-1)] + ['y'+str(i) for i in range(1, n-1)]

# Créer l'anneau de polynômes en x1, x2,..., xn-2 et y1, y2,..., yn-2
R = PolynomialRing(F_qn, monomials)

# Convertir les monomiaux en polynômes dans R
monomials = [R(monom) for monom in monomials]

# Construire g
g = construct_g(R, monomials, u, v, n)

print("computed g")

substitutions = {monomials[i]: xy[i] for i in range(2*(n-2))}

evaluation = g.subs(substitutions)
print(evaluation)

coefficients = []

#quotient, reminder = g.quo_rem(t^(k-1))

#print("Décomposition de g dans la base des puissances de t:", quotient)