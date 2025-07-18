#######################################################
### Non-Crossing Polygons, Hypergraph and Werner States
### Written Summer 2025 by Lukas Walko
#######################################################

### MUST READ "orbdim.gap", etc. FOR SOME FUNCTIONS TO WORK --
### (Run this function)--> Read("orbdim.gap");
### Dependencies are marked

# In the process of refactoring / reorganizing. Dependencies will change.
# TODO: Update dependencies

### Some functions require "Float" package
### This must be downloaded and installed
### Visit https://gap-packages.github.io/float/
### for more information

Print( "Read ncp_hyper.gap" );

########## Hypergraph and Non-Crossing Chord Methods ##########

###
# NonCrossingPairs() takes a list of points and
# returns a list of all non-crossing pairs
# For example: lst := [1,2,3,4,5,6]
#
# Returns [[ [ 1, 2 ], [ 3, 4 ], [ 5, 6 ] ]
#          [ [ 1, 2 ], [ 3, 6 ], [ 4, 5 ] ]
#          [ [ 1, 4 ], [ 2, 3 ], [ 5, 6 ] ]
#          [ [ 1, 6 ], [ 2, 3 ], [ 4, 5 ] ]
#          [ [ 1, 6 ], [ 2, 5 ], [ 3, 4 ] ]]
###
NonCrossingPairs := function( lst )
    local res, i, left, right, leftPairs, rightPairs, lp, rp;

    if Length( lst ) = 0 then
        return [ [] ];
    fi;

    res := [];

    for i in [2,4..Length( lst )] do
        # Pair lst[ 1 ] with lst[ i ]
        left := lst{ [2..i-1] };
        right := lst{ [i+1..Length( lst )] };

        leftPairs := NonCrossingPairs( left );
        rightPairs := NonCrossingPairs( right );

        for lp in leftPairs do
            for rp in rightPairs do
                Add( res, Concatenation( [[ lst[ 1 ], lst[ i ] ]], lp, rp ));
            od;
        od;
    od;

    return res;
end;


###
# REQUIRES "orbdim.gap"
# PairsToHypergraph() turns lists of non-crossing pairs
# and turns them into hypergraph states
#
# numQubits is the number of qubits for the hypergraph state
# lsts is a list of lists, like that supplied by the function NonCrossingPairs()
#
# res is a list of states, where states are lists / vectors themselves
###
PairsToHypergraph := function( numQubits, lsts )
    local sublst, i, res, subres;

    res := [];

    for i in [1..Length( lsts )] do
        sublst := lsts[ i ];
        subres := HypergraphState( numQubits, sublst );

        Add( res, subres );
    od;

    return res;
end;

###
# StatesToMatrices() takes states produced by
# PairsToHypergraph() or HypergraphState() and displays
# the length n^2 state vectors as nxn matrices,
# with "1" in place of 1, and "-" in place of -1
#
# Function returns nothing
###
StatesToMatrices := function( states )
    local i, j, n, sublst, count ;

    n := Int( Sqrt( Length( states[ 1 ]))); # Sqrt(n^2), each state has n^2 elements

    # Write matrix
    for sublst in states do
        for i in [1..n] do
            Print( "[ " );

            for j in [1..n] do
                count := (i-1) * n + j;

                if sublst[ count ] = 1 then
                    Print( "1 " );
                else
                    Print( "- " );
                fi;
            od;
            Print( "]\n" );
        od;
        Print( "\n\n" );
    od;

end;

###
# REQUIRES "orbdim.gap"
#
# Bitstring is a list of 0s and 1s such as
# [ 1, 0, 0, 1, 1]
#
# Returns C(I) where I = bitstring
###
CofI := function( bitstring )
    local k, sum, rootUnity, n;

    n := Length( bitstring );

    rootUnity := E( n );

    sum := 0;

    for k in [0..(n - 1)] do
        prod := BitstringTensor( VirginiaCycleK( k, bitstring ) );

        sum := sum + ( rootUnity^k * prod );
    od;

    return normalize( sum );
end;

###
# SingletProduct() takes a list of pairs of the form:
# [ [a, b], [c, d],... ]
# where each pair specifies the qubit positions that form the singlet
#
# Function returns list of lists of the form:
# [ [1,0,...,1], [0,0,...1],... ]
# where each sublist is a term in the expansion of the tensor product
###
SingletProduct := function( pairs )
    local n, pair, result, vec1, vec2, newResult, s;

    # Determine number of qubits (largest index)
    n := 2 * Length( pairs );

    # Start with the empty product: one all-zero vector
    result := [ List([1..n], i -> 0) ];

    for pair in pairs do
        newResult := [];

        for s in result do
            # For each existing term, create new terms from the singlet

            vec1 := ShallowCopy( s );
            vec1[ pair[ 1 ] ] := 0;
            vec1[ pair[ 2 ] ] := 1;

            vec2 := ShallowCopy( s );
            vec2[ pair[ 1 ] ] := 1;
            vec2[ pair[ 2 ] ] := 0;

            Add( newResult, vec1 );
            Add( newResult, vec2 );
        od;

        result := newResult;
    od;

    return result;
end;

###
# REQUIRES "orbdim.gap"
#
# SingletTensor() returns the state vector of a
# non-crossing-chord diagram
# whose chords are specified by the list of lists called pairs
###
SingletTensor := function( pairs )
    local singlets, base, permutation, res, i, digits;

    singlets := List( [1..Length( pairs )], i -> singlet );

    # initial vector before permutation
    base := normalize( KronVec( singlets ) );

    # Creates permutation by comparing [1..numQubits] to Flat( [[polygon1], [polygon2],...] )
    permutation := PartialPerm( List( [1..Length( pairs ) * 2] ), Flat( pairs ) );

    res := PermuteQubitsPsi( base, permutation );

    digits := [];

    for i in [1..Length( res )] do
        if res[ i ] <> 0 then
            Add( digits, res[ i ] );
        fi;
    od;

    if digits[ 1 ] < 0 then
        res := (-1) * res;
    fi;

    return res;
end;

###
# REQUIRES "orbdim.gap"
#
# CreateNCCMatrix() takes the number of qubits
# in the CHORD diagrams as its argument
# Returns a matrix where each row is an NCC diagram's state vector
###
CreateNCCMatrix := function( numQubits )
    local allPairSets, pairs, res, prod, singletProducts, i, val, state;

    res := [];

    allPairSets := NonCrossingPairs( [1..numQubits] );

    singletProducts := [];

    # Creates a list of all the singlet product expansions for each of the pairs
    for pairs in allPairSets do

        prod := SingletProduct( pairs );
        Add( singletProducts, prod );
    od;

    # For each singlet expansion, replace the expansion
    # with the base 10 value of the minimum
    # value in its set
    for i in [1..Length( singletProducts )] do

        val := Bin( singletProducts[i][1] );
        singletProducts[i] := val;
    od;

    # Sort the values of the singlet expansions and their corresponding pairs in parallel
    SortParallel( singletProducts, allPairSets );

    # Fill in the result matrix with the state vectors from each set of pairs
    # Due to the previous sorting, the matrix now places the states of the
    # NCC diagrams in their natural ordering
    for pairs in allPairSets do

        state := SingletTensor( pairs );
        Add( res, state );
    od;

    return res;
end;

########## Non-Crossing Polygon State Methods ##########

###
# Generates all the NCP states for n qubits
#
# Returns a vector of the *flattened* rhos
###
GenerateNCPs := function( numQubits )
    local diagrams, rhos, d;

    rhos:= [];

    # generate all the polygons that describe the NCPs
    diagrams := NCPartitionsSet( [1..numQubits] );

    # create a flattened rho for each diagram and add it to the list
    for d in diagrams do
        Add( rhos, Flat( WernerDiagram( d ) ) );
    od;

    # return flattened rhos
    return rhos;
end;

###
# SplitNCP() takes an n-qubit NCP diagram and
# returns the chord pairs of the corresponding 2n-qubit NCC diagram
#
# lsts is a list of lists describing the connected vertices
#
# Function returns a list of chord pairs
# (which are represented by two element lists)
###
SplitNCP := function( lsts )
    local sublst, i, res, subres, j, chord;

    res := [];

    for sublst in lsts do
        subres := [];

        for i in [1..Length( sublst )] do

            # Add 2a_i - 1 and 2a_i to the list
            Add( subres, ( (2 * sublst[ i ]) - 1 ) );
            Add( subres, ( 2 * sublst[ i ] ) );
        od;

        for j in [1..( Length( subres ) - 1 ) ] do
            chord := [];

            if j = 1 then

                # Create the chord between the first and last elements
                Add( chord, subres[ 1 ] );
                Add( chord, subres[ Length(subres) ] );

                # Add the chord to the output list
                Add( res, chord );
            else
                if j mod 2 = 0 then

                    # Create the chord between 2a_i and 2a_(i+1) -1
                    Add( chord, subres[ j ] );
                    Add( chord, subres[ j + 1 ] );

                    # Add the chord to the output list
                    Add( res, chord );
                fi;
            fi;

        od;
    od;

    return res;
end;

###
# REQUIRES "orbdim.gap" AND "werner.gap"
#
# Returns the coefficient of the (numQubits - 1)-gon
# As an answer to the question what linear combination
# of (numQubits - 1) NCPs equals the
# N-gon with the nth qubit traced out
###
TraceOutNGon := function( numQubits )
    local tracedNGon, nGon, underGons, targetQubits, gon, res, place, count;

    # creates the n-gon
    nGon := WernerDiagram( [ [1..numQubits] ] );

    # creates n-entry list of 0s, where final entry is 1
    targetQubits := List( [1..numQubits], i -> 0 );
    targetQubits[ numQubits ] := 1;

    # traces out the nth qubit of the n-gon
    tracedNGon := PartialTrace( nGon, targetQubits );

    # list of NCP(n-1) states
    underGons := [];

    count := 0;
    place := 0;
    for gon in NCPartitionsSet( [1..(numQubits - 1)] ) do
        # add each NCP(n-1) to underGons
        Add( underGons, Flat( WernerDiagram( gon ) ) );

        # keeps track of the position of the (n-1)-gon in underGons
        count := count + 1;
        if Length( gon ) = 1 then
            place := count;
        fi;
    od;

    if res = fail then
        Display( "TraceOutNGon: No linear combination of NCP( n - 1 ) exists" );
        return -1;
    fi;

    # solve for the linear combination of NCP(n-1)s that equal
    # traced out n-gon
    res := SolutionMat( underGons, Flat( tracedNGon ) );

    # returns the coefficient of the (n-1)-gon
    return [res[ place ], res];
end;

###
# REQUIRES FLOAT PACKAGE
# Returns the decimal representations
#
# TraceOutNGonDecimal() functions the same as
# TraceOutNGon(), except it returns the decmial
# representation of results
###
TraceOutNGonDecimal := function( numQubits )
    local exactLst, decLst, pos, value, count;

    exactLst := TraceOutNGon( numQubits );

    # Find and store position of (n-1)-gon
    pos := 0;
    count := 1;
    for value in exactLst[ 2 ] do
        if exactLst[ 1 ] = value then
            pos := count;
        fi;

        count := count + 1;
    od;

    # Convert all cyclotomic values to decimal approximation
    decLst := List( exactLst[ 2 ], Float );

    # Return decimals in same pattern as TraceOutNGon() output
    return [ decLst[ pos ], decLst ];
end;

###
# TraceOutQudits() is a more general form of TraceOutNGon()
# which allows any local dimension dim, and tracing any number of qudits
# specified by a binary list, targetqudits. Ex: [0,0,0,1] -> traces out 4th qudit
# numQudits is the number of qudits of the n-gon being partial traced
#
# Function returns a list of the form:
# [ coefficient of top-gon, [ coefficients of lower gons including top-gon ] ]
###
TraceOutQudits := function( dim, numQudits, targetQudits )
    local tracedNGon, nGon, underGons, gon, res, place, count, quditDiff;

    # creates the n-gon
    nGon := WernerDiagramQudits( dim, [ [1..numQudits] ] );

    # traces out the target qudits of the n-gon
    tracedNGon := PartialTraceQudits( dim, nGon, targetQudits );

    # list of NCP(n-1) states
    underGons := [];

    quditDiff := Count( 1, targetQudits );

    count := 0;
    place := 0;
    for gon in NCPartitionsSet( [1..(numQudits - quditDiff)] ) do
        # add each NCP(n-1) to underGons
        Add( underGons, Flat( WernerDiagramQudits( dim, gon ) ) );

        # keeps track of the position of the (n-1)-gon in underGons
        count := count + 1;
        if Length( gon ) = 1 then
            place := count;
        fi;
    od;

    # solve for the linear combination of NCP(n-1)s that equal
    # traced out n-gon
    res := SolutionMat( underGons, Flat( tracedNGon ) );

    if res = fail then
        Display( "TraceOutQudits: No linear combination of NCP( n - quditDiff ) exists" );
        return -1;
    fi;

    # returns the coefficient of the (n-1)-gon
    return [res[ place ], res];
end;
