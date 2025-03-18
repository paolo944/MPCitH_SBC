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
        copy_mat = copy(self.mat)
        copy_signature = copy(self.signature)
        copy_mat = matrix(copy_mat.rows()+[row])
        copy_signature[copy_mat.nrows()-1] = index
        return SignedMatrix(copy_mat, copy_signature, self.d, self.parent)

    # use position over term ordering
    def row_echelon_form_by_position(self):
        # returns a pair (M, n) where M is a new signed matrix which is the row-reduction of self via a sequence of
        # elementary row operations
        # keep track of number of operations
        num_operations = 0
        copy_mat = copy(self.mat)
        eliminated = True
        first_reduction = True
        # keep track of reductions
        rdxn = dict()
        for i in range(len(copy_mat.rows())):
            rdxn[i] = []
        while eliminated:
            eliminated = False
            for i, row in enumerate(copy_mat.rows()):
                for j in range(len(row)):
                    if row[j] != 0:
                        # j is the leading term of this row, so use it to kill everything with higher signature
                        for new_i, new_row in enumerate(copy_mat.rows()):
                            if new_row[j] != 0 and self.signature[i] < self.signature[new_i]:
                                # we can reduce
                                lam = -(new_row[j]/row[j])
                                copy_mat.add_multiple_of_row(new_i, i, lam)
                                eliminated = True
                                if first_reduction: # only count top-reductions
                                    num_operations += len(new_row)
                                rdxn[new_i].append((i,lam))
                        break
            first_reduction = False # stop counting arithmetic operations

        for i, row in enumerate(copy_mat.rows()):
            for j in range(len(row)):
                if row[j] != 0:
                    # j is the coefficient of the leading term of this row, so divide this row by it
                    copy_mat.rescale_row(i,1/row[j])
                    break

        return (SignedMatrix(copy_mat, self.signature, self.d, self.parent), num_operations, rdxn)

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
        M_red[0] = copy(M_red[1])
        M_red[1], _, _ = M.row_echelon_form_by_position()
        print(f"computed M_tilde for d={d}")
    return (M.mat, M.signature, M_red[1])


system = load("system.sobj")

system = sorted(system, key=lambda p: p.degree())

Mac, sig, Mac_red = F5(system, 4)

lignes_a_zero = 0

for i in range(Mac_red.mat.rows()):
    if i.is_zero():
        lignes_a_zero += 1

print(f"nombres de réductions à 0: {lignes}")