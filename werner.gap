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
                                                        
