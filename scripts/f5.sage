import time

class SignedMatrix:
    # Matrix together with hashmap associating singature (index) to each row - has special rref function that respects signatures.
    def __init__(self, mat, sgn, d, parent):
        self.mat = mat
        self.signature = sgn
        self.d = d
        self.parent = parent

    def LT(self):
        # returns the leading terms of the (polynomials represented by) rows of self.mat
        monomials = self.monomials()
        leading_terms = []
        for row in self.mat.rows():
            for i in range(len(row)):
                if row[i] != 0:
                    leading_terms.append(monomials[i]*row[i])
                    break
        return set(leading_terms)

    def monomials(self):
        # returns monomials of degree self.d in a list, sorted in decreasing order
        R = self.parent
        monomials = R.monomials_of_degree(self.d)
        monomials.sort(reverse=True)
        return monomials

    def rows(self):
        # return set of (polynomials represented by) rows of self.mat
        monomials = self.monomials()
        r = []
        for row in self.mat.rows():
            polynomial = 0
            for j in range(len(row)):
                polynomial += row[j]*monomials[j]
            r.append(polynomial)
        return r 

    def add_row(self, f, index):
        # returns a new matrix which is self with a row added corresponding to polynomial f with signature index
        row = [f.monomial_coefficient(mon) for mon in self.monomials()]
        self.mat = self.mat.augment(vector(row))
        self.signature.append(index)
        return self		

    # use position over term ordering
    def row_echelon_form_by_position(self):
        # returns a pair (M, n) where M is a new signed matrix which is the row-reduction of self via a sequence of
        # elementary row operations
        copy_mat = copy(self.mat)
        eliminated = True
        while eliminated:
            eliminated = False
            for i, row in enumerate(copy_mat.rows()):
                for j in range(len(row)):
                    if row[j] != 0:
                        # j is the leading term of this row, so use it to kill everything with higher signature
                        for new_i, new_row in enumerate(copy_mat.rows()):
                            if new_row[j] != 0 and self.signature[i] < self.signature[new_i]:
                                # we can reduce
                                copy_mat.add_multiple_of_row(new_i, i, 1)
                                eliminated = True
                        break
        return SignedMatrix(copy_mat, self.signature, self.d, self.parent)

def F5(F, D):
    # F=(f_1,...,f_m) is a set of polynomials with degere d_1 <= d_2 <= ... <= d_m
    # D is maximal degree
    # returns the set of elements of degree at most D of reduced Grobner bases of (f_1,...,f_i) for each i
    M_red = [SignedMatrix(matrix(GF(2)), dict(), 2, F[0].parent()), SignedMatrix(matrix(GF(2)), dict(), 3, F[0].parent())] # initialize reduced Macaulay matrices
    variables = list(F[0].parent().gens())
    variables.sort(reverse=True)
    for d in range(F[0].degree(),D+1):
        M = SignedMatrix(matrix(GF(2)), dict(), d, F[0].parent())
        for i in range(len(F)):
            if d < F[i].degree(): continue
            elif d == F[i].degree(): # Case 1: the degree of f_i is exactly d
                M = M.add_row(F[i], (i, 1))
            else: # Case 2: the degree of f_i is less than d
                Crit = M_red[2-F[i].degree()].LT() # build F_5 criterion list
                assert(M_red[1].mat.nrows() == len(M_red[1].signature))
                for j in range(M_red[1].mat.nrows()):
                    if M_red[1].signature[j][0] == i:
                        _,e = M_red[1].signature[j]
                        f = M_red[1].rows()[j]
                        if e == 1:
                            largest_var_in_e = 0
                        else:
                            largest_var_in_e = variables.index(e.variables()[-1]) # select which row to use to build new row
                        for k in range(largest_var_in_e,len(variables)):
                            if e*variables[k] not in Crit: # avoid signatures which F_5 criterion tells us are useless
                                M = M.add_row(variables[k]*f, (i,e*variables[k]))
        # reduce Macaulay matrix
        print("computing reducation of Mac")
        M_red[0] = copy(M_red[1])
        M_red[1] = M.row_echelon_form_by_position()
        print(f"computed M_tilde for d={d}")
    return (M.mat, M.signature, M_red[1])


system = load("system.sobj")

system = sorted(system, key=lambda p: p.degree())

start = time.time()

Mac, sig, Mac_red = F5(system, 4)

end = time.time()

lignes_a_zero = 0

for i in range(Mac_red.mat.nrows()):
    if all(Mac_red.mat[i, j] == 0 for j in range(Mac_red.mat.ncols())):  # Vérifie si toute la ligne est à 0
        lignes_a_zero += 1         

print(f"nombres de réductions à 0: {lignes_a_zero} / {Mac_red.mat.nrows()}")
print(f"temps: {end - start}")
