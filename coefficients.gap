################################
# Linear Combinations Checking #
################################

Print("Read coefficients.gap");

# Tracing out qudits and checking to see what linear combinations of NCP( n-1 ) create them
# Computing other sets of coefficients.
# Triseparability and Biseparability checks found in here

###                                    ###
### FLOAT PACKAGE NEEDED FOR THIS FILE ###
###                                    ###

###
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

# Biseparability condition (a)
BisepCondA := function( rPlus, rMinus, rOne, rTwo, rThree )
    local pt1, pt2;
    # Condition one
    pt1 := 3*rMinus-1 <= 1 + rOne - rMinus -2*rPlus and 1 + rOne - rMinus -2*rPlus <= 0;
    
    # Condition two
    pt2 := 3*(rTwo^2)+3*(rThree^2)+(1+2*rOne+rMinus-rPlus)^2 <= (2+rOne-4*rMinus-2*rPlus)^2;
    
    # Both must be true for condition (a) to be true
    return pt1 and pt2;
end;

# Biseparability condition (b)
BisepCondB := function( rPlus, rMinus, rOne, rTwo, rThree )
    local pt1, pt2;
    # Condition one
    pt1 := 0 <= 1 + rOne - rMinus - 2*rPlus and 1 + rOne - rMinus - 2*rPlus <= 1-3*rMinus;
    
    # Condition two
    pt2 := 3*(rTwo^2)+3*(rThree^2)+(1-3*rMinus-3*rPlus)^2 <= (rOne+2*rMinus-2*rPlus)^2;
    
    # Both must be true for (b) to be true
    return pt1 and pt2;
end;

# Checks biseparability and triseparability. Takes the output of WEBasisCoeffs as input.
# Output is a list of "true" and "false" corresponding to each condition
CheckSeparability := function( lst )
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
        Add( trisep, true );
    else
        Add( trisep, false );
    fi;
    
    # Check biseparability
    bisep := [];
    if 0 <= rMinus and rMinus <= 1/3 then
        Add( bisep, true );
    else
        Add( bisep, false );
    fi;
    if BisepCondA(rPlus, rMinus, rOne, rTwo, rThree) then
        Add( bisep, true );
    else
        Add( bisep, false );
    fi;
    if BisepCondB(rPlus, rMinus, rOne, rTwo, rThree) then
        Add( bisep, true );
    else
        Add( bisep, false );
    fi;
    # Collect all conditions into one output
    Add( res, trisep ); Add( res, bisep );
    
    return res;
end;
