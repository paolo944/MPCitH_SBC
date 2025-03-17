def mac_sparse_dict(list_poly, monom):
    r"""
    Build the Macaulay matrix for the polynomials of degree d in list_poly. It is represented as a couple of a matrix and the list of the signatures of its rows, assuming compatibility.
    The matrix is sparse.
 
    INPUT:
    - ``list_poly`` -- a list of degree d polynomials
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    OUTPUT:
    The Macaulay matrix defined by list_poly
    EXAMPLES::
        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.F5Trop_doctest import Make_Macaulay_Matrix_sparse_dict
        sage: list_poly = [[0,1,x**2+y**2], [1,1,x*y], [2,1,y*z], [3,1,z**2]]
        sage: Mac = Make_Macaulay_Matrix_sparse_dict(list_poly,[x**2, x*y, y**2, x*z, y*z, z**2])
        sage: Mac[0]
        [1 0 1 0 0 0]
        [0 1 0 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: Mac[1]
        [[0, 1], [1, 1], [2, 1], [3, 1]]
        
    """

    R = list_poly[0][2].parent()
    d = list_poly[0][2].degree()
    l = len(monom)

    dict_monom = {monom[aa].exponents()[0]:aa for aa in range(l)}
 
 
    nrows = len(list_poly)
 
    listsign = []
    Mac = MatrixSpace(R.base_ring(),nrows,l,sparse=True)(0)
    for u in range(nrows):
        fuple=list_poly[u]
        listsign.append([fuple[0],fuple[1]])
        f = fuple[2]
        list = f.exponents()
        for mon in list:
            if dict_monom.has_key(mon):
                j = dict_monom[mon]
                Mac[u,j] = f.monomial_coefficient(R({mon:1}))
    Mac = [Mac, listsign]
    return Mac

def echelon(M):
    n = A.nrows()

    for i in range(n):
        for j in range(i+1, n):
            if A[i, i] != 0:
                factor = A[j, i] / A[i, i]
                A.set_row(j, A.row(j) - factor * A.row(i))

def f5(F, D):
    R = F[0].parent()
    m = len(F)
    l = 
    for d in range(2, D+1):
        Mac = MatrixSpace(R,0,,sparse=True)(0)



system = load("system.sobj")
f5(system, 4)