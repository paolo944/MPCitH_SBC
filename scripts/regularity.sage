try:
    system = load(sys.argv[1])
except:
    print("Erreur pendant le load")
    sys.exit(1)

system2 = []

print(system)

for i in system:
    try:
        system2.append(i.homogeneous_components()[2])
    except KeyError:
        system2.append(i.homogeneous_components()[1])

print(system2)

I = ideal(system2)

field = GF(2)
ring.<t> =  field[]

n = system[0].parent().ngens()

denom = 1
for i in system2:
    print(i.degree())
    denom *= (1+t**i.degree())

num = (1+t)**n
q, r = num.quo_rem(denom)

print(f"test de vérification: {q*denom + r == num}")
print(f"q: {q}")
print(f"r: {r}")
print(f"The hilbert series of I is: {I.hilbert_series()}")
print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")
print(f"Hilbert polynomial of I: {I.hilbert_polynomial()}")
