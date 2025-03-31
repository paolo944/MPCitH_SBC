def construct_mac(F):
    m = len(F)
    n = F[0].parent().ngens()

if(len(sys.argv) != 3):
    print("il manque des arguments")
    sys.exit(1)

n = int(sys.argv[1])
m = int(sys.argv[2])

F = load(f"system/sage/system_bilin_{n}_{m}.sobj")

construct_mac(F)
