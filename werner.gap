### code mostly by SNW, summer 2013

###
# Edited Summer 2025 by Lukas Walko
###

Print("Read werner.gap");

mix := n -> normalizeHS( Kron( replicate( n )( I ) ) );

############################################
# Constructing non-crossing diagram states #
############################################

adjacentPairs := xs -> zip( xs, Concatenation( cdr( xs ), [ car( xs ) ] ) );

filterByPair := pr -> function( xs )
    if pr[ 1 ] < pr[ 2 ] then
        return filter( xs, x -> pr[ 1 ] <= x and x < pr[ 2 ] );
    else
        return filter( xs, x -> pr[ 1 ] <= x or  x < pr[ 2 ] );
    fi;
end;

partitionByList := ms -> xs -> List( adjacentPairs( SortedList( ms ) ), pr -> filterByPair( pr )( xs ) );

nonCrossing := ms -> xs -> Length( filter( partitionByList( ms )( xs ), compose( Not )( IsEmpty ) ) ) = 1;

nonCrossingPair := pr -> nonCrossing( pr[ 1 ] )( pr[ 2 ] );

NCPartition := part -> And( fmap( nonCrossingPair )( Combinations( part, 2 ) ) );

# Non-crossing partitions
# Give a list [1..numQubits] -> returns all of the possible NCP diagrams
NCPartitionsSet := xs -> filter( PartitionsSet( xs ), NCPartition );

################################
# Constructing Werner Diagrams #
################################

###
# REQUIRES "orbdim.gap"
# Generates 2^n by 2^n permutation matrices
# perm is a permutation written in cycle notation
###
PermMatGen := function( n, perm )
    local M, N, row;

    M := IdentityMat( 2^n );

    N := [];

    for row in M do
        Add( N, PermuteQubitsPsi( row, perm ) );
    od;

    # return TransposedMat( N );
    return N;

end;

###
# REQUIRES "orbdim.gap"
# Generates dim^n by dim^n permutation matrices
# perm is a permutation written in cycle notation
###
PermMatGenQudits := function( dim, n, perm )
    local M, N, row;

    M := IdentityMat( dim^n );

    N := [];

    for row in M do
        Add( N, PermuteQudits( dim, row, perm ) );
    od;

    return N;

end;

# use (-j) power for Virginia reel cycling
RootCycleV := function( psi )
    local n, gen;

    n := BitSize( psi );
    gen := GeneratorsOfGroup( CyclicGroup( IsPermGroup, n ) )[ 1 ];

    return Sum( List([0..n-1], j -> E(n)^j * PermuteQubitsPsi( psi, gen^( -j ))));
end;

###
# General version of RootCycleV() for local dimension dim
###
RootCycleQudits := function( dim, lst )
    local n, gen;

    n := DitSize( dim, lst );
    gen := GeneratorsOfGroup( CyclicGroup( IsPermGroup, n ) )[ 1 ];

    return Sum( List( [0..n-1], j -> E(n)^j * PermuteQudits( dim, lst, gen^( -j ) )) );
end;

###
# WernerGon(n) returns the conjugate of the Werner state of an n-vertex Non-Crossing-Polygon diagram
# For example, WernerGon(3) returns the state represented by a triangle connecting qubits [1, 2, 3]
###
WernerGon := n ->
  normalizeDM(Sum(List([0..2^n-1],composeList( [DM, RootCycle, ee, x->dec2binlist( x, n )] ))));

# This is the complex conjugate of WernerGon
# ( So the version of the Werner State you'd expect )
WernerGonV := n ->
  normalizeDM( Sum( List( [0..2^n-1], composeList( [ DM, RootCycleV, ee, x->dec2binlist( x, n ) ]))));
  
###
# REQUIRES "orbdim.gap"
# Returns density matrix of a numQudits-Gon diagram in local dimension dim
###
WernerGonQudits := function( dim, numQudits )
    local res, sz;

    sz := dim^numQudits - 1;

    res := Sum( List( [0..sz], x->DM( RootCycleQudits( dim, eeDit( dim, dec2baseNlist( x, dim, numQudits ))))));

    return normalizeDM( res );
end;

###
# WernerDiagram([ [ ] ]) takes a list of lists, for example: [ [1,2], [3] ]
# Each sublist in the argument "lsts" describes the vertices of a polygon / chord / point
# Using the above example, [1,2] is a chord between vertices 1 and 2, and [3] is vertex 3 by itself
#
# WernerDiagram() returns the density operator of the state
###
WernerDiagram := function( lsts )
    local flat, n, dmlst, res;

    flat := Concatenation( lsts );
    n := Length( flat );
    Assert( 0, SortedList( flat ) = [1..n] );
    dmlst := List( lsts, compose( WernerGonV )( Length ));

    res := OrderedKronList( lsts, dmlst );

    return res;
end;

###
# REQUIRES "orbdim.gap"
# WernerDiagramQudits is the general form of WernerDiagram()
# lsts is a list of lists, dim is the local dimension
# Function returns rho of the specified diagram state
###
WernerDiagramQudits := function( dim, lsts )
    local flatLsts, numQudits, densityMatrices, res, gon, permMatrices, permMatrix, perm, perms;

    # Flatten lsts and find numQudits in diagram
    flatLsts := Flat( lsts );
    numQudits := Length( flatLsts );

    # Checks that lsts has as many qudits as the maximum of lsts
    Assert( 0, SortedList( flatLsts ) = [1..numQudits] );

    densityMatrices := [];

    # Adds the n-gon state for each polygon in the diagram to a list
    for gon in lsts do
        Add( densityMatrices, WernerGonQudits( dim, Length( gon ) ) );
    od;

    res := Kron( densityMatrices );

    # Creates partial permutation mapping [1..numQudits] to flatLsts
    # Ex: PartialPerm([1..5], Flat( [[1,2,5], [3,4]] )) = (1)(2)(3,5,4)
    # Technically the [1..5] is superfluous, but I've kept it for clarity
    perm := PartialPerm( [1..numQudits], flatLsts );

    permMatrix := PermMatGenQudits( dim, numQudits, perm );

    return permMatrix^(-1) * res * permMatrix;
end;

########################################
# Checking the Werner-basis conjecture #
########################################

WernerBasisCheck := n -> RankMat( List( NCPartitionsSet( [1..n] ), 
                            compose( flatten )( WernerDiagram ) ) );

####################
# Dotless diagrams #
####################

dotless := nss -> Count( 1, List( nss, Length ) ) = 0;

# Dotless Non-crossing partitions
DotlessNCPartitionsSet := xs -> filter( PartitionsSet( xs ), nss -> NCPartition( nss ) and dotless( nss ) );

DotlessNCPartitionsSetImperative := function( xs )
    local ps, result;
    
    ps := PartitionsSet( xs );
    result := [];
    
    for part in ps do
        if dotless( part ) and NCPartition( part ) then
            Add( result, part );
        fi;
    od;
    
    return result;
end;

############################
# Dotless Diagram Checking #
############################

DotlessWernerBasisCheck := n -> RankMat( List( DotlessNCPartitionsSet( [1..n] ), 
                                                compose( flatten )( WernerDiagram ) ) );

DotlessWernerBasisCheckImperative := n -> RankMat( List( DotlessNCPartitionsSetImperative( [1..n] ), 
                                                        compose( flatten )( WernerDiagram ) ) );
                                                        

################################
# Linear Combinations Checking #
################################

# Tracing out qudits and checking to see what linear combinations of NCP( n-1 ) create them

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

# numQudits is the number of qudits of the n-gon to trace down from. targetQudits is a list of the qudits to trace out.
# dim is the local dimension
# targetQudits must have 3 zeros!
# This code finds the coefficients of linear combinations in Werner and Eggeling's basis that give traced n-gons
# NOT A CHECK OF THEIR PAPER, THIS IS SOMETHING ELSE
WEBasisLinCombs := function( dim, numQudits, targetQudits )
    local nGon, state, weBasis, res, Id, V12, V23, V31, V123, V321, rPlus, rMinus, rZero, rOne, rTwo, rThree;
    # Check that targetQudits leaves 3 qubits untouched
    Assert( 0, Count( 0, targetQudits ) = 3 );
    
    # Basis for 3-qudit Werner space defined by Werner and Eggeling
    # GAP doesn't like the permutation (1,2,3), so it's written in the equivalent form (2,3,1)
    # Below are the necessary permutation matrices
    Id := IdentityMat(dim^3);  V12 := PermMatGenQudits(dim, 3, (1,2));  V23 := PermMatGenQudits(dim, 3, (2,3)); 
    V31 := PermMatGenQudits(dim ,3, (3,1)); V123 := PermMatGenQudits(dim, 3, (2,3,1)); 
    V321 := PermMatGenQudits(dim, 3, (3,2,1));

    rPlus := 1/6 * ( Id + V12 + V23 + V31 + V123 + V321 );

    rMinus := 1/6 * ( Id - V12 - V23 - V31 + V123 + V321 );

    rZero := 1/3 * ( 2 * Id - V123 - V321 );

    rOne := 1/3 * ( 2 * V23 - V31 - V12 );

    rTwo := 1/Sqrt( 3 ) * ( V12 - V31 );

    rThree := E(4)/Sqrt(3) * ( V123 - V321 );
    
    # Create numQudits-gon and trace out all but 3 qudits
    nGon := WernerDiagramQudits( dim, [[1..numQudits]] );
    state := PartialTraceQudits( dim, nGon, targetQudits );
        
    # weBasis is the basis defined by Werner and Eggeling
    weBasis := List( [ rPlus, rMinus, rZero, rOne, rTwo, rThree ], Flat );
    
    # Find coefficients of W.E. basis that return traced n-gon
    res := SolutionMat( weBasis, Flat( state ) );
    
    return res;
end;

WEBasisCoeffs := function( dim, numQudits, targetQudits )
    local nGon, state, weBasis, R, res, Id, V12, V23, V31, V123, V321, rPlus, rMinus, rZero, rOne, rTwo, rThree;
    # Check that targetQudits leaves 3 qubits untouched
    Assert( 0, Count( 0, targetQudits ) = 3 );
    
    # Basis for 3-qudit Werner space defined by Werner and Eggeling
    Id := IdentityMat(dim^3);  V12 := PermMatGenQudits(dim, 3, (1,2));  V23 := PermMatGenQudits(dim, 3, (2,3)); 
    V31 := PermMatGenQudits(dim ,3, (3,1)); V123 := PermMatGenQudits(dim, 3, (2,3,1)); 
    V321 := PermMatGenQudits(dim, 3, (3,2,1));

    rPlus := 1/6 * ( Id + V12 + V23 + V31 + V123 + V321 );

    rMinus := 1/6 * ( Id - V12 - V23 - V31 + V123 + V321 );

    rZero := 1/3 * ( 2 * Id - V123 - V321 );

    rOne := 1/3 * ( 2 * V23 - V31 - V12 );

    rTwo := 1/Sqrt( 3 ) * ( V12 - V31 );

    rThree := E(4)/Sqrt(3) * ( V123 - V321 );
    
    # Create numQudits-gon and trace out all but 3 qudits
    nGon := WernerDiagramQudits( dim, [[1..numQudits]] );
    state := PartialTraceQudits( dim, nGon, targetQudits );
    
    # Collect Werner and Eggeling Basis for 'for' loop
    weBasis := [rPlus, rMinus, rZero, rOne, rTwo, rThree];
    
    # Calculate Tr( state * R ) for all R in weBasis, collect into res list; preserves order of weBasis
    res := [];
    for R in weBasis do
        Add( res, Trace( state * R ) );
    od;
    
    return res;
end;

# Checks biseparability and triseparability. Takes the output of WEBasisCoeffs as input.
# Output is a list of "true" and "false" corresponding to each condition
CheckSeperability := function( lst )
    local res, rPlus, rMinus, rZero, rOne, rTwo, rThree, trisep, bisep;
    # Separate W.E. Basis into parts
    rPlus := lst[1]; rMinus := lst[2]; rZero := lst[3]; rOne := lst[4]; rTwo := lst[5]; rThree := lst[6];
    
    res := [];
    
    # Check triseparability
    trisep := [];
    if 0 <= rMinus and rMinus <= 1/6 then 
        Add( trisep, true );
    else
        Add( trisep, false );
    fi;
    if 1/4 * (1 - 2*rMinus) <= rPlus and rPlus <= 1 - 5*rMinus then
        Add( trisep, true );
    else
        Add( trisep, false );
    fi;
    if (3*(rThree^2)+(1-3*rPlus)^2)*(1-6*rMinus) <= (rOne+rPlus-rMinus)*((rMinus-2*rPlus+2*rMinus)^2 - 3*(rTwo^2)) then
        Add( trisep true );
    else
        Add( trisep false );
    fi;
    
    # Check biseparability
    bisep := [];
end;
