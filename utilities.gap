#############################################
### Utilities, Constants, and Helper Methods
### Written Summer 2025 by Lukas Walko
#############################################

Print( "Read utilities.gap" );

########## Basic operators, constants, and values ##########

# Hadamard gate
hadamard := 1/Sqrt( 2 ) * [ [1,1], [1,-1] ];

# Pauli X
pauliX := [ [0,1], [1,0] ];

# Pauli Y
pauliY := [ [0, -Sqrt( -1 )], [Sqrt( -1 ), 0] ];

# Pauli Z
pauliZ := [ [1, 0], [0, -1] ];

# 1 / Sqrt( 2 )
invroot2 := 1 / Sqrt( 2 );

# |0〉
ketzero := [ [1], [0] ];

# |1〉
ketone := [ [0], [1] ];

# |01〉- |10〉
singlet := invroot2 * [ 0, 1, -1, 0 ];

######
### Methods
######

# GAP's not is not a function
Not := b -> not b;

# BlanksString() is a helper function for PrintMatrix
# Returns a list of n spaces
BlanksString := function( n )
    return Concatenation( List( [1..n], i -> " " ) );
end;

###
# PrintMatrix() takes a list of lists and prints the
# formatted version to the console
# Output looks like:     [[...]
#                         [...]
#                         [...]]
# Function returns nothing
###
PrintMatrix := function( matrix )
    local row, elem, maxLen, rowStr, strElem;

    # Determine the maximum width of any element when converted to string
    maxLen := 0;

    for row in matrix do
        for elem in row do
            strElem := String( elem );

            if Length( strElem ) > maxLen then
                maxLen := Length( strElem );
            fi;
        od;
    od;

    # Print each row with elements left-aligned
    for row in matrix do

        rowStr := "[";

        for elem in row do
            strElem := String( elem );
            rowStr := Concatenation(rowStr,
                        strElem,
                        BlanksString( maxLen - Length( strElem ) ),
                        " ");
        od;

        Print( rowStr, "]\n" );
    od;

    Print( "\n" );

end;

###
# Function find all the non-zero entries of a matrix
# and logs their locations and values in the format
# [ row, col, value ] ***Assuming indexing starts at zero!
###
NonZeroLocations := function( matrix )
    local row, col, rowDim, colDim, res;

    rowDim := Length( matrix );
    colDim := Length( matrix[ 1 ] );

    res := [];

    for row in [1..rowDim] do
        for col in [1..colDim] do

            if matrix[ row ][ col ] <> 0 then
                Add( res, [ row - 1, col - 1, matrix[ row ][ col ] ] );
            fi;
        od;
    od;

    return res;

end;

# Helper function, swaps two objects in a list
Swap := function( lst, index1, index2 )
    local temp;

    temp := lst[ index1 ];
    lst[ index1 ] := lst[ index2 ];
    lst[ index2 ] := temp;
end;

# Returns the conjugate transpose of a matrix
ConjugateTranspose := mat -> ComplexConjugate( TransposedMat( mat ) );

###
# REQUIRES "orbdim.gap"
# Returns <a|b>
###
InnerProductMats := function( a, b )
    return Trace( ConjugateTranspose( a ) * b );
end;

###
# virginiaCycleList( [] ) takes a list, and moves the first entry to the end while advancing all other entries up one position
###
virginiaCycleList := xs -> Concatenation( xs{ [2..Length( xs )] }, xs{ [ 1 ] } );

###
# REQUIRES "orbdim.gap"
# Virginia cycles list k times
#
# Function changes the list given to it and
# returns the result
###
VirginiaCycleK := function( k, lst )
    local i;

    for i in [1..k] do
        lst := virginiaCycleList( lst );
    od;

    return lst;
end;

# Helper function to create a basis vector of length n with specified entries
MakeBasisVector := function( n, pos1, bit1, pos2, bit2 )
    local vec;

    vec := List( [1..n], i -> 0 );  # initialize to all 0s
    vec[ pos1 ] := bit1;
    vec[ pos2 ] := bit2;

    return vec;
end;

# Finds the maximum value in a matrix
FindMax := function( matrix )
    local row, col, big, current;

    # initialize value big to be the first entry in the matrix
    big := matrix[ 1 ][ 1 ];

    for row in [1..Length( matrix )] do
        for col in [1..Length( matrix[ row ] )] do
            current := matrix[ row ][ col ];

            # make current value big if it's larger than all previous values
            if big < current then
                big := current;
            fi;
        od;
    od;

    return big;
end;

# Returns the inner product of two vectors
InnerProduct := function( cvec1, cvec2 )
    return ComplexConjugate( cvec1 ) * cvec2;
end;

# Returns the outer product of two vectors
OuterProduct := function( cvec1, cvec2 )
    # Depends on ConjugateTranspose
    local colvec1, rowvec2;
    
    colvec1 := TransposedMat( [ cvec1 ] );
    rowvec2 := [ ComplexConjugate( cvec2 ) ];
    
    return colvec1 * rowvec2;
end;

###
# DM for density matrix
# constructs density matrix for vector cvec
###
DM := function( cvec )
    # depends on ConjugateTranspose
    local colvec;
    
    colvec := TransposedMat( [ cvec ] );
    
    return colvec * ConjugateTranspose( colvec );
end;

norm := function( cvec )
    local sum, elt;
    
    sum := 0;
    for elt in cvec do
        sum := sum + elt * ComplexConjugate( elt );
    od;
    
    return Sqrt( sum );
end;

normalize := function( cvec )
    # uses norm
    return cvec / norm( cvec );
end;

###
# Returns normalized density matrix
# Trace( normalizeDM( matrix ) ) = 1
###
normalizeDM := function( dm )
    Assert( 0, Trace( dm ) <> 0 );
    return dm / TraceMat( dm );
end;

# Kronecker product for vectors
kron := function( vec1, vec2 )
    return KroneckerProduct( [ vec1 ], [ vec2 ] )[ 1 ];
end;

# Konecker product for a list of vectors
KronVec := function( list )
    local vec, sumVec;
    
    sumVec := [ 1 ];
    for vec in list do
        sumVec := kron( sumVec, vec );
    od;
    
    return sumVec;
end;

# Kronecker product for a list of matrices
Kron := function( list )
    local mat, sumMat;
    
    sumMat := [ [ 1 ] ];
    for mat in list do
        sumMat := KroneckerProduct( sumMat, mat );
    od;
    
    return sumMat;
end;

# Returns lst with the first entry removed
cdr := function( lst )
  return lst{ [2..Length( lst )] };
end;

### Gram-Schmidt orthonormalization
GramSchmidt := function( lst )
    local helper;
    
    helper := function( ortho, rest )
        local new, term;
        
        if IsEmpty( rest ) then
            return ortho;
        else
            new := rest[ 1 ];
            
        for term in ortho do
            new := new - term * InnerProduct( term, rest[ 1 ] );
        od;
        
        new := new / Sqrt( InnerProduct( new, new ) );
        
        return helper( Concatenation( ortho, [ new ] ), cdr( rest ) );
        fi;
    end;
    
    return helper([],lst);
end;
