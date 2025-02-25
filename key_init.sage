# File to create a pair of public and private keys for
# NSBC problem by J. Huth and A. Joux
# article: https://eprint.iacr.org/2023/1685.pdf

def linear_dependance(u, v):
    M = matrix([u, v])
    return M.rank() < 2

def save_vector(vector, file_name):
    F = vector[0].parent()
    variables = F.gens()
    with open(file_name, "w") as f:
        if(variables != (1)):
            for i in variables[:-1]:
                f.write(i + ",")
            f.write(str(variables[-1]) + "\n")
        f.write(str(F.characteristic()) + "\n")
        for v in vector[:-1]:
            f.write(str(v) + ",\n")
        f.write(str(v))

q = 2 # Field characteristic
n = 130 # Vectors's size
k = 257 # Degree of field extension

F_q = GF(q)
F_qn = GF(q^k, 't')

while(True):
    x_prime = [F_q.random_element() for _ in range(n-2)]
    y_prime = [F_q.random_element() for _ in range(n-2)]

    x = x_prime + [1, 0]
    x = vector(x)

    y = y_prime + [0, 1]
    y = vector(y)

    u = [F_qn.random_element() for _ in range(n)]
    v = [F_qn.random_element() for _ in range(n-1)]

    u = vector(u)

    u_y = u.dot_product(y)
    u_x = u.dot_product(x)

    v_x = 0
    v_y = 0

    for i in range(0, n-1):
        v_x += v[i]*x[i]
        v_y += v[i]*y[i]

    denominator = y[n-1] * u_x

    if(denominator == 0):
        continue

    v_n = (u_y*v_x - u_x*v_y) / denominator

    v.append(v_n)
    v = vector(v)
    
    # Check if u and v are linearly dependant
    if(linear_dependance(u, v)):
        continue

    # If all is good
    break

# Test the instance
v_x = v.dot_product(x)
v_y = v.dot_product(y)

assert u_x * v_y == u_y * v_x, "NSBC not right"

save_vector(v, "v.pub")
save_vector(u, "u.pub")
save_vector(x, "x")
save_vector(y, "y")

print("generated and saved private and public keys")