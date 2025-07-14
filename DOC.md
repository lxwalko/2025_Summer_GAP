# LVC MathPhys GAP Code Documentation

### Preface

Each file will have its own header, underneath which the methods contained within are documented.

### Table of Contents

1. load.gap
2. utilities.gap
3. dits.gap
4. werner.gap
5. orbdim.gap
6. ncp_hyper.gap
7. lin_from_eigen.gap

### Common Syntax / Symbols

Below are some common names, syntax, etc. that are frequently used along with what they mean. I hope this saves you the
confusion that some of these terms caused me.

- `binlist`: A binary list
- `lsts`: A list of lists; often used to describe the chords / polygons of an NCC or NCP diagram. This may also refer to a matrix
- `dm`: Density matrix
- `cvec`: Sometimes `colvec` means column vector; represented by a list
- `res`: The result / return value of a function
- `rootUnity`: A root of unity
- `perm`: A permutation; written in cycle notation, as this is the format recognized by GAP
- `dim`: The local dimension of a qudit / qudit system
- `k` or `n`: A number, generally specifying a number of elements or iterations

## load.gap

`load.gap` is a file purely written for convienience. Instead of reading each file in individually, all the files
in this repo are collected in one place so only one file needs to be read. No methods are defined here.

## utilities.gap
