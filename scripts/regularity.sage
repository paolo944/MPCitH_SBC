try:
    system = load(sys.argv[1])
except:
    print("Erreur pendant le load")
    sys.exit(1)
system2 = []
#print(system)
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
    denom *= (1+t**i.degree())

denom = series_ring(denom)
num = series_ring((1+t)**n)
res = num * denom.inverse_of_unit()
print(f"Generating series of the sequence: {res}")

hilbert_series = I.hilbert_series()
hilbert_series_denom = hilbert_series.denominator()
hilbert_series_num = hilbert_series.numerator()
res2 = series_ring(hilbert_series_num) * series_ring(hilbert_series_denom).inverse_of_unit()
print(f"The hilbert series of I is: {res2}")
print(f"degree of semi-regularity of I: {I.degree_of_semi_regularity()}")
