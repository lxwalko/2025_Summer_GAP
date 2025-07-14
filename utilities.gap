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

###
# IsIndependent() checks if a list of states is linearly independent
# states is a list of lists, like that produced by PairsToHypergraph()
#
# Returns true if linearly independent
# Returns false if not
###
IsIndependent := function( states )
    if RankMat( states ) = Length( states ) then
        return true;
    fi;

    return false;
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
