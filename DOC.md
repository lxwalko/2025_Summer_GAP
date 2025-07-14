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

`utilities.gap` begins with the definitions of some commonly used constants and operators, such as the Paulis.  
  
`BlanksString( n )`:  
- Input: `n` is a nonnegative integer  
- Output: A list with length `n`, each element is the space character, `" "`
  
`PrintMatrix( matrix )`:
- Input: `matrix` is a matrix, or a list of lists
- Output: A "neatly" formatted representation of `matrix` is printed to the terminal
- Note: On large matrices, `PrintMatrix` makes the output uglier than it would be normally
  
`NonZeroLocations( matrix )`:
- Input: `matrix` is a matrix, or a list of lists
- Output: A list of lists which logs the locations and values of each non-zero entry in `matrix`. The result list has
entries of the form `[ row, col, value ]` where `row` and `col` are nonnegative integers, and `value` is the entry
- Note: Outputs assume indexing starts at 0, as this code was made as a helper function for use outside of GAP. GAP
begins indexing at 1.
  
`IsIndependent( states )`:
- Input: `states` is a matrix, assumed to be made of either state vectors from NCC diagrams, or flattened density matrices
from NCP diagrams.
- Output: Either `true` or `false`, depending on if the rows of the matrix are linearly independent
  
