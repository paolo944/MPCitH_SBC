system = load("system.sobj")

system2 = []

print(system)

for i in system:
    try:
        system2.append(i.homogeneous_components()[2])
    except KeyError:
        system2.append(i.homogeneous_components()[1])

print(system2)

I = ideal(system2)

K = FractionField(PolynomialRing(QQ, 't'))

t = K.gen()

n = system[0].parent().ngens()

denom = 1
for i in system2:
    print(i.degree())
    denom *= (1+t**i.degree())

hilbert = (1+t)**n/denom

print((hilbert - I.hilbert_series(algorithm="singular")).is_zero())

print(I.degree_of_semi_regularity())
