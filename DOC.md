# LVC MathPhys GAP Code Documentation

### Preface

Each file will have its own section. Within each section the functions, methods, and other items are documented. Assume
that any function from a file lower in the table of contents / load order depends on functions found in file(s) above it.

### Table of Contents

1. [load.gap](#loadgap)
2. [utilities.gap](#utilitiesgap)
3. [dits.gap](#ditsgap)
4. [werner.gap](#wernergap)
5. [orbdim.gap](#orbdimgap)
6. ncp_hyper.gap
7. lin\_from\_eigen.gap
8. coefficients.gap

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
- `dim`: The local dimension of a qubit / qudit system
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
  
`Not( b )`:
- Input: `b` is a boolean value.
- Output: The opposite boolean value of `b`.
  
`BlanksString( n )`:  
- Input: `n` is a nonnegative integer.
- Output: A list with length `n`, each element is the space character, `" "`.
  
`PrintMatrix( matrix )`:
- Input: `matrix` is a matrix, or a list of lists.
- Output: A "neatly" formatted representation of `matrix` is printed to the terminal.
- Note: On large matrices, `PrintMatrix` makes the output uglier than the unmodified output.
  
`NonZeroLocations( matrix )`:
- Input: `matrix` is a matrix, or a list of lists.
- Output: A list of lists which logs the locations and values of each non-zero entry in `matrix`. The result list has
entries of the form `[ row, col, value ]` where `row` and `col` are nonnegative integers, and `value` is the entry.
- Note: Outputs assume indexing starts at 0, as this code was made as a helper function for use outside of GAP. GAP
begins indexing at 1.
  
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

`kron( vec1, vec2 )`:
- Input: `vec1` and `vec2` are vectors.
- Output: The kronecker ( or tensor ) product of `vec1` and `vec2`.

`KronVec( list )`:
- Input: `list` is a list of vectors.
- Output: The kronecker product of every vector in `list`.

`Kron( list )`:
- Input: `list` is a list of matrices.
- Output: The kronecker product of every matrix in `list`.

`cdr( lst )`:
- Input: `lst` is a list.
- Output: A copy of `lst` with the first element removed.

`GramSchmidt( lst, inner )`:
- Input: `lst` is a list

## dits.gap

`dits.gap` contains methods that pertain mostly to bits, bitstrings, dits, and ditstrings. Methods for binary or n-ary 
counting, ditstring generation, bitstring tensors, etc. are found here.  
  
### Functions
  
`BitstringTensor( bitstring )`:
- Input: `bitstring` is a list of zeros and ones.
- Output: The normalized state vector born out of ket-`bitstring`.
  
`BitSize( lst )`:
- Input: `lst` is a list.
- Output: <code>log<sub>2</sub>(Length( lst ))</code> rounded down to the nearest integer.
  
`DitSize( d, lst )`:
- Input: `d` is an integer specifying the base. `lst` is a list.
- Output: <code>log<sub>d</sub>(Length( lst ))</code> rounded down to the nearest integer.
  
`BitComplement( binlist )`:
- Input: `binlist` is a binary list; a list of zeros and ones.
- Output: A binary list where each one in `binlist` becomes a zero in the return list, and vice versa. 
For example, `BitComplement( [1, 0, 0, 1] )` -> `[0, 1, 1, 0]`.
  
`dec2binlist( dec, numqubits )`:
- Input: `dec` and `numqubits` are nonnegative integers.
- Output: A binary list of length `numqubits`, which evaluates to `dec` when converted to base 10.
  
`dec2base4list( dec, numqubits )`:
- Input: `dec` and `numqubits` are nonnegative integers.
- Output: A base 4 list of length `numqubits`, which evaluates to `dec` when converted to base 10.
  
`dec2baseNlist( dec, N, numdigits )`:
- Input: `dec`, `N`, and `numdigits` are nonnegative integers.
- Output: A base `N` list of length `numdigits`, which evaluates to `dec` when converted to base 10.
  
`DecOfBinstring( binstring )`:
- Input: `binstring` is a binary string.
- Output: The base 10 value of `binstring`.
  
`CompositeBitlist( subsystem, list0, list1 )`:
- Input: `subsystem` is a binary list. `list0` and `list1` are lists. All lists have equal length.
- Output: A list with the same length as `subsystem`, created by inserting entries from `list0` when the corresponding bit
in `subsystem` is zero, and inserting from `list1` when the corresponding bit is one.
  
`OneBinList( k, n )`:
- Input: `k` and `n` are nonnegative integers.
- Output: A binary list of length `n` with a one in position `k`, all zeros elsewhere.
  
`e( num, bin )`:
- Input: `num` and `bin` are nonnegative integers.
- Output: A list of length <code>2<sup>num</sup></code> where the entry in position `bin + 1` is one, all zeros elsewhere.

`eDit( dim, pow, pos )`:
- Input: `dim`, `pos`, and `pos` are nonnegative integers.
- Output: A list of length <code>dim<sup>pow</sup></code> where the entry in position `pos + 1` is one, all zeros elsewhere.

`ee( lst )`:
- Input: `lst` is a list.
- Output: A list of length <code>2<sup>Length( lst )</sup></code> where the entry in position `Bin( lst ) + 1` is one, all zeros elsewhere.

`eeDit( dim, lst )`:
- Input: `dim` is a nonnegative integer. `lst` is a list.
- Output: A list of length <code>dim<sup>Length( lst )</sup></code> where the entry in position `BaseN( lst ) + 1` is one, all zeros elsewhere.

`bitflip( l, bitlist )`:
- Input: `l` is a positive integer. `bitlist` is a binary list.
- Output: `bitlist` where the bit in the `l`th position has been flipped iff `bitlist[ l ]` is one.

`Bin( vec )`:
- Input: `vec` is a vector of integers.
- Output: The binary sum of the entries in `vec`.
- Note: `Bin` does not require the entries of `vec` to be ones and zeros. `Bin` is the sum of the entries in `vec`
multiplied by increasing powers of two starting from the right.

`BaseN( base, vec )`:
- Input: `vec` is a vector of integers. `base` is an integer.
- Output: The base `base` sum of the entries in `vec`.
- Note: `BaseN` is the general form of `Bin`.
  
`BaseNCounting( base, numDigits )`:
- Input: `base` is an integer denoting the base of the counting (such as base 2, base 3, etc.). `numDigits` is a nonnegative
integer specifying the number of digits in each ditstring.
- Output: A list of ditstrings that count from `0` to <code>base<sup>numDigits - 1</sup></code>.

## werner.gap

`werner.gap` contains the code necessary to create and experiment upon Werner states. Unlike other sections in `DOC.md`, the
functions of `werner.gap` are documented in an order that attempts to reflect their importance to the user, rather than order of
appearance.

### Functions

`WernerDiagram( lsts )`:
- Input: `lsts` is a list of lists, which describes an NCP diagram.
- Output: The density matrix for the state described by `lsts` in local dimension two.

`WernerDiagramQudits( dim, lsts )`:
- Input: `dim` is a nonnegative integer, describing the local dimension of each constituent qudit. `lsts` is a list of lists
which describes an NCP digram.
- Output: The density matrix for the state described by `lsts` in local dimension `dim`.

`PermMatGen( n, perm )`:
- Input: `n` is a nonnegative integer which describes the number of qubits in the system. `perm` is a GAP permutation 
written in cycle notation.
- Output: A <code>2<sup>n</sup></code> by <code>2<sup>n</sup></code> permutation matrix.

`PermMatGenQudits( dim, n, perm )`:
- Input: `dim` is a positive integer which describes the local dimension of each qudit. `n` is a nonnegative integer which describes the number of qudits. 
`perm` is a GAP permutation written in cycle notation.
- Output: A <code>dim<sup>n</sup></code> by <code>dim<sup>n</sup></code> permutation matrix.

`WernerGonV( n )`:
- Input: `n` is a nonnegative integer that describes the number of qubits.
- Output: The density matrix of the `n`-gon in local dimension two.
- Note: There is a twin function, `WernerGon`, which has the same input but returns the complex conjugate of `WernerGonV`.

`WernerGonQudits( dim, numQudits )`:
- Input: `dim` is an nonnegative integer which describes the local dimension of each qudit. `numQudits` is a nonnegative
integer which describes the number of qudits.
- Output: The density matrix of the `numQudits`-gon in local dimension `dim`.

`NCPartitionsSet( xs )`:
- Input: `xs` is a list from one to *n*, where *n* is the number of qudits.
- Output: A list of lists describing all possible NCP diagrams with *n* qudits.

`DotlessNCPartitionsSet( xs )`:
- Input: `xs` is a list from one to *n*, where *n* is the number of qudits.
- Output: A list of lists describing all dotless NCP diagrams with *n* qudits.

## orbdim.gap

`orbdim.gap` holds code from summers past in one very long, mostly unorganized file. The infrastructure found in
`orbdim.gap` that was relevant to the endevours of summer 2025 was moved to respective files. The remaining code is
mostly uncommented and can be difficult to understand. If you have a need not met with the other files, it may be found within 
`orbdim.gap`. Parsing that search is an exercise left for the user.

## ncp_hyper.gap

`ncp_hyper.gap` handles hypergraph states, and their relationships to NCC diagrams.
