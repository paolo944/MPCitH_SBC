import sympy as sym


x = sym.Symbol("x")
y = sym.Symbol('y')
yy = sym.Symbol('yy')

""" Paramètres """
k = 24
o2 = 31
n = 80
m = 160

d1 = 7
d2 = 6


Sx = sym.series(1/(1-x), x=x, n=d1)
Sy = sym.series(1/(1-x), x=x, n=d1).subs(x, y)
Syy = sym.series(1/(1-yy), x=yy, n=d1)
Sxy = sym.expand(Sx*Syy)

""" 
    Cette fonction prend en entré une série permet d'enlever les coefficients qui ne correspondent pas à des paramètres 
    réaliste.
    Exemple : X**3 * Y**5 ne correspond pas à un paramètre possible.
"""
def ser(tmp):
    s2 = 0
    flag = True
    i = 1
    while flag:
        if tmp.coeff(y, i) == 0:
            flag = False
        else:
            j = i+1
            while j < d1:
                s2 = sym.expand(s2 + tmp.coeff(y, i).coeff(x, j) * x ** j * y ** i)
                j += 1
            i += 1
    return s2


"""
    Cette fonction prend en entré une série et sort la forme dé-homogéniser de la série
    i.e. pour une série S en entrée, la fonction renvoit Y/[(1-X)(Y-1)]*S
"""
def susu(s):
    return ser(sym.expand(Sxy * ser(s)).subs(yy, 1 / y))


s = sym.series(1/(1+x**2)**m, x=x, n=d1).subs(x, x*y)

t = sym.series((1+x)**k, x=x, n=d1).subs(x, x*y)
tt = sym.series((1+x)**(n-k), x=x, n=d1)
M = t * tt

H = sym.expand((tt - s * t * tt)/y)
""" Print la série H(X,Y) """
print(H)
print(ser(H))
print("")

G = susu(H)

tmp = sym.series((1+y)**k/(1+y**2)**m, x=y, n=d1)
J = ser(sym.expand(Sx * Sy * tmp))

""" Print la série G(X,Y) """
print(ser(G))
print()

""" Print la série J(X,Y) """
print(f"J: {ser(G-J)}")
print(" ")

"""
    Print la série correspondant au degré de régularité d'un système homogène de k variables et m polynômes 
    i.e. la fonction renvoit la série : (1+Y)**k/(1+Y**2)**m
"""
print(tmp)
print("")
HilbertSer = sym.series((1+x)**n/(1+x**2)**m, x=x, n=d1)

""" 
    Print la série correspondant au degré de régularité d'un système homogène de n variables et m polynômes
    i.e. la fonction renvoit la série : (1+X)**n/(1+X**2)**m
"""
print(HilbertSer)

