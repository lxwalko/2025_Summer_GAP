# LVC MathPhys GAP Code Documentation

### Preface

Each file will have its own section. Within each section the functions, methods, and other items are documented.

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
- `cvec`: A vector; represented by a list
- `res`: The result / return value of a function
- `rootUnity`: A root of unity
- `perm`: A permutation; written in cycle notation, as this is the format recognized by GAP
- `dim`: The local dimension of a qudit / qudit system
- `k` or `n`: A number, generally specifying a number of elements or iterations
- `mat`: A matrix
  
#### A note on terminology

Frequently list and vector are used interchangably. In similar fashion, so are matrix and list of lists.

## load.gap

`load.gap` is a file purely written for convienience. Instead of reading each file in individually, all the files
in this repo are collected in one place so only one file needs to be read. No methods are defined here.

## utilities.gap

`utilities.gap` begins with the definitions of some commonly used constants and operators, such as the Paulis.  
The functions found in this file are commonly used throughout the rest of the code in the repo.  
  
### Functions
  
`BlanksString( n )`:  
- Input: `n` is a nonnegative integer.
- Output: A list with length `n`, each element is the space character, `" "`.
  
`PrintMatrix( matrix )`:
- Input: `matrix` is a matrix, or a list of lists.
- Output: A "neatly" formatted representation of `matrix` is printed to the terminal.
- Note: On large matrices, `PrintMatrix` makes the output uglier than it would be normally.
  
`NonZeroLocations( matrix )`:
- Input: `matrix` is a matrix, or a list of lists.
- Output: A list of lists which logs the locations and values of each non-zero entry in `matrix`. The result list has
entries of the form `[ row, col, value ]` where `row` and `col` are nonnegative integers, and `value` is the entry.
- Note: Outputs assume indexing starts at 0, as this code was made as a helper function for use outside of GAP. GAP
begins indexing at 1.
  
`IsIndependent( states )`:
- Input: `states` is a matrix, assumed to be made of either state vectors from NCC diagrams, or flattened density matrices
from NCP diagrams.
- Output: Either `true` or `false`, depending on if the rows of `states` are linearly independent.
  
`Swap( lst, index1, index2 )`:
- Input: `lst` is a list. `index1` is the index of the first element. `index2` is the index of the second element.
- Output: Swaps the elements at `index1` and `index2` in place. Returns nothing.
  
`ConjugateTranspose( mat )`:
- Input: `mat` is a matrix.
- Output: The conjugate transpose (or dagger) of a matrix.
  
`InnerProductMats( a, b )`:
- Input: `a` and `b` are square matrices.
- Output: A scalar value, the inner product of `a` and `b`.
  
`virginiaCycleList( xs )`:
- Input: `xs` is a list.
- Output: A list where the first entry of `xs` is moved to the end, and all other elements are moved up one position.
  
`VirginiaCycleK( k, lst )`:
- Input: `k` is a positive integer, `lst` is a list.
- Output: A list where the first element of `lst` has been cycled to the back `k` times.
- Note: `lst` is changed in place, meaning `lst` will not be the same after being passed to `VirginiaCycleK`.
  
`MakeBasisVector( n, pos1, bit1, pos2, bit2 )`:
- Input: `n`, `pos1`, and `pos2` are positive integers, `bit1` and `bit2` are integers.
- Output: A vector of length `n` where `pos1` and `pos2` have values `bit1` and `bit2` respectively.
  
`FindMax( matrix )`:
- Input: `matrix` is a matrix.
- Output: The maximum value of `mat`.
  
`InnerProduct( cvec1, cvec2 )`:
- Input: `cvec1` and `cvec2` are vectors. 
- Output: A scalar value --  the inner product of the two vectors.
  
`OuterProduct( cvec1, cvec2 )`:
- Input: `cvec1` and `cvec2` are vectors.
- Output: A matrix -- the outer product of the two vectors.
  
`DM( cvec )`:
- Input: `cvec` is a vector.
- Output: The density matrix formed by the outer product of `cvec` and itself.
  
`norm( cvec )`:
- Input: `cvec` is a vector.
- Output: The normalization factor of `cvec`.
  
`normalize( cvec )`:
- Input: `cvec` is a vector.
- Output: The normalized form of `cvec`.
  
`normalizeDM( dm )`:
- Input: `dm` is a density matrix.
- Output: The normalized form of `dm`.

## dits.gap

`dits.gap` contains methods that pertain mostly to bits, bitstrings, dits, and ditstrings. Methods for binary or n-ary 
counting, ditstring generation, bitstring tensors, etc. are found here.  
  
### Functions
  
`BitstringTensor( bitstring )`:
- Input: `bitstring` is a list of zeros and ones.
- Output: The normalized state vector born out of ket-`bitstring`.
  
`BaseNCounting( base, numDigits )`:
- Input: `base` is an integer denoting the base of the counting (such as base 2, base 3, etc.). `numDigits` is a nonnegative
integer specifying the number of digits in each ditstring.
- Output: A list of ditstrings that count from `0` to `base^(numDigits - 1)`.
  
`BitSize( lst )`:
- Input: `lst` is a list.
- Output: log<sub>2</sub>(`Length( lst )`) rounded down to the nearest integer.
  
`DitSize( d, lst )`:
- Input: `d` is an integer specifying the base. `lst` is a list.
- Output: log<sub>d</sub>(`Length( lst )`) rounded down to the nearest integer.
