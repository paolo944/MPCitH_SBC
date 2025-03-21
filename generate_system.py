import sys
import subprocess

formats = ["hpXbred", "msolve", "magma"]

if(len(sys.argv) != 4):
    print("Il manque des arguments")
    print("Il faut lancer comme ceci:")
    print("\tpython3 generate_system.py 4 hpXbred 1")

try:
    n = int(sys.argv[1])
except ValueError:
    print("Erreur, la taille des vecteur doit être un entier")
    sys.exit(1)

formatting = sys.argv[2]

try:
    field_eq = int(sys.argv[3])
except ValueError:
    print("Erreur, le booleen pour inclure ou non les équations du corps doit être un entier")
    sys.exit(1)

if(n < 1):
    print("La taille des vecteurs doit être supérieur ou égale à 1")
    sys.exit(1)

if(formatting not in formats):
    print("Le nom du format doit correspondre à un des formats pris en compte:")
    print(f"\t{formats}")
    sys.exit(1)

if(field_eq != 0 and field_eq != 1):
    print("Le booleen pour inclure les équations du corps doit être 0 ou 1")
    sys.exit(1)

launch = ["sage", "scripts/key_init.sage", str(n)]

result = subprocess.run(launch) 

if result.returncode != 0:
    sys.exit(1)

launch = ["./model", str(n), formatting, str(field_eq)]

result = subprocess.run(launch, text=True)

if result.returncode != 0:
    print(result.stderr)
    sys.exit(1)
