# MPCItH_SBC

Library to generate instances of NSBC(normalized subfield bilinear collision
) Problem according to the paper: Huth, J., Joux, A.: MPC in the head using
subfield bilinear collision problem. Cryptology ePrint Archive, Paper 2023/1685
(2023), https://eprint.iacr.org/2023/1685

The system equations are affine because of the normalisation but they can easily
become homogeneous and bilinear by adding 4 homogenization variables.

## Dependencies:

- C compiler(gcc or clang)
- Make
- sagemath
- FLINT
- python3

## Quick Install:

make

## Usage:

You can generate a pair of public/private keys as the following:

```console
$ sage scripts/key_init.sage {n}
```

where n is the size of the vectors according to the paper.
The number of variables nvar is equal to 2 * (n-2) and the number of equations
m is equal to nvar + 1

This pair of public/private keys are stored in human friendly format and in 
.sobj format to use within sagemath.

The private key is in $(F_2)^n$. The public key is in $(F_{2^{m}})^n$.

To generate the system of equations on $F_2$ by Weil descent, you have 2
options, either use sagemath to have a system written in .sobj format but it's
very slow for big systems, or use the C compiled code using FLINT to have 1 of
the 3 supported formats.

For the first option, use:

```console
$ sage scripts/modelisation.sage {n} {0 or 1}
```

where {n} is the same vector size as for the key generation and the second 
option is a boolean for wether you want to include the field equations or not.

For the second option, use:
```console
$ ./model {n} {format} {0 or 1}
```

It has the same options except for the format option which indicates the output
format. You can choose between the following formats:

- msolve
- hpXbred
- magma

It's case sensitive.

The output system will be generated in the directory 
/system/{format}/system_bilin_n_m.{format}

Finally, if you only want the final system, you can juste type:

```console
python3 generate_system.py {n} {format} {0 or 1}
```

Where n is the vector size, format is one of the following: 
[sage, msolve, hpXbred, magma] and 0 or 1 for if you want to include the field's
equations.

## Bug Report:

For any bug report or any question:
paul.mekhail@etu.sorbonne-universite.fr
