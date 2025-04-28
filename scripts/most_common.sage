import subprocess
from sage.all import *
import random

#load('scripts/key_init.sage')
#load('scripts/modelisation.sage')

def analyze_system_with_subs(system, k):
    all_monomials = set()
    for poly in system:
        all_monomials.update(poly.monomials())

    #print(all_monomials)

    all_vars = system[0].parent().gens()

    target_vars = all_vars[:k]

    assignments = {v: random.randint(0, 1) for v in target_vars}

    stripped_monomials = set()
    for poly in system:
        poly_spe = poly.specialization(assignments)
        stripped_monomials.update(poly_spe.monomials())

    #print(f"Original number of monomials: {len(all_monomials)}")
    #print(f"Remaining monomials after setting top {k} vars to 1: {len(stripped_monomials)}")

    return len(all_monomials), len(stripped_monomials)


def analyze_system_with_subs_v2(system, k):
    all_monomials = set()
    all_vars = system[0].parent().gens()
    count_vars = {}

    for i in all_vars:
        count_vars[i] = 0

    for poly in system:
        tmp_monomials = poly.monomials()
        all_monomials.update(tmp_monomials)
        for monomial in tmp_monomials:
            for variables in monomial.variables():
                for _, i in variables:
                    count_vars[i] += 1

    count_vars = {k: v for k, v in sorted(count_vars.items(), key=lambda item: item[1], reverse=True)}
    vars_sorted = list(count_vars)

    most_common_vars = vars_sorted[-k:]
    assignments = {v: random.randint(0, 1) for v in most_common_vars}

    stripped_monomials = set()
    for poly in system:
        poly_spe = poly.specialization(assignments)
        stripped_monomials.update(poly_spe)

    #print(f"Original number of monomials: {len(all_monomials)}")
    #print(f"Remaining monomials after setting top {k} vars to 0 or 1: {len(stripped_monomials)}")

    return len(all_monomials), len(stripped_monomials)

def run_batch(n, output_file="most_common_results.txt"):
    # Open output file for writing results
    with open(output_file, 'w') as f:
        f.write(f"n = {2*(n-2)} m = {2*(n-2) + 1} n_x = n_y = {int(n-2)}\n")
        f.write("k,total_monomials,remaining_monomials,nb_poly\n")
        n_y = n-2
        for k in range(n_y // 2, n_y+1):
                system = load("system/sage/system_bilin_96_193.sobj")
                total_monomials, remaining_monomials = analyze_system_with_subs(system, k)
                f.write(f"{k},{total_monomials},{remaining_monomials},{4*(n-2) + 1 - k}\n")
                print(f"k = {k} / {n_y} complete for {2*(n-2)} variables")

if __name__ == '__main__':
    num_vars = 50
    output_filename = f"analysis_results_{num_vars}_most.txt"
    #system = load("system/sage/system_bilin_56_113.sobj")
    #analyze_system_with_subs(system,1)
    run_batch(num_vars, output_filename)