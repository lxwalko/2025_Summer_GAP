#######################################
## Linear States from Eigenvectors Code
#######################################

###
# REQUIRES "orbdim.gap"
# Returns the 'dopey string' of a given polygon diagram
# The 'dopey string' is a bitstring made of minimal sub-bitstrings
# in the positions of the polygons in the diagram
# Ex: DopeyString( [ [1,4], [2, 3], [5, 6, 7] ] ); -> [ 0, 0, 1, 1, 0, 0, 1]
###
DopeyString := function( diagram )
    local polygon, qubit, bitstring, numQubits, gonLength, count;

    numQubits := Length( Flat( diagram ) );
    bitstring := List( [1..numQubits], x -> 0 );

    for polygon in diagram do
        gonLength := Length( polygon );
        count := 1;

        for qubit in polygon do
            if count = gonLength then
                bitstring[ qubit ] := 1;
            fi;

            count := count + 1;
        od;
    od;

    return bitstring;
end;

# Calculates the normalization factor of psi_D
CDNorm := function( diagram )
    local polygon, normFactor;

    normFactor := 1;
    for polygon in diagram do
        normFactor := normFactor * Length( polygon );
    od;

    return 1 / Sqrt( normFactor );
end;

# Cycles the target qubits in a bitstring k times
CycleQubitsK := function( k, bitstring, targetQubits )
    local i, substring, newBitstring;

    substring := [];
    for i in targetQubits do
        Add( substring, bitstring[ i ] );
    od;

    substring := VirginiaCycleK( k, substring );

    newBitstring := ShallowCopy( bitstring );
    for i in [1..Length( targetQubits )] do
        newBitstring[ targetQubits[ i ] ] := substring[ i ];
    od;

    return newBitstring;
end;

# Handles the recursion of looping over each polygon Length( polygon ) times
PolygonIteration := function( bitstring, diagram, dCount )
    local targetQubits, i, rootUnity, totalVec, recursiveVec, newBitstring;

    if dCount > Length( diagram ) then
        return BitstringTensor( bitstring );
    fi;

    targetQubits := diagram[ dCount ];
    totalVec := [];

    for i in [1..Length( targetQubits )] do
        rootUnity := E( Length( targetQubits ) )^i;
        newBitstring := CycleQubitsK( i, bitstring, targetQubits );

        recursiveVec := rootUnity * PolygonIteration( newBitstring, diagram, dCount + 1 );

        totalVec := totalVec + recursiveVec;
    od;

    return totalVec;
end;

# Calculates the state vector for a diagram using the dopey string
CDofI := function( diagram, bitstring )
    local normFactor, rootUnity, polygon, qubit, resVec;

    normFactor := CDNorm( diagram );
    resVec := PolygonIteration( bitstring, diagram, 1 );

    return normFactor * resVec;
end;

# Comparison function
# Sorts so that N-gon is last,  all dots is first
CompareDiagrams := function( diagramOne, diagramTwo )
    if Length( diagramOne ) < Length( diagramTwo ) then
        return false;
    fi;
    return true;
end;

###
# Creates Linear States from Eigenstates Matrix
#
# Each i,j entry where D = diagram is <psi_Di | rho_Dj | psi_Di>
###
CreateLSEMatrix := function( numQubits )
    local diagrams, polygonOperators, stateVecs, matrix, row, operator, state;
    # init matrix
    matrix := [];
    # Create NCP diagrams and sort so that all dots is first, and n-gon is last
    diagrams := NCPartitionsSet( [1..numQubits] );
    StableSort( diagrams, CompareDiagrams );
    # Create rhos and psis for each NCP diagram
    polygonOperators := List( diagrams, WernerDiagram );
    stateVecs := List( diagrams, x -> CDofI( x, DopeyString( x ) ) );

    for state in stateVecs do
        row := [];
        for operator in polygonOperators do
            # Add <psi_Di | rho_Ej | psi_Di> to matrix
            Add( row, InnerProduct( state, (operator * state) ) );
        od;

        Add( matrix, row );
    od;

    return matrix;
end;

###
# Creates a normalized LSE Matrix. Each entry is divided
# by <psi_Di | rho_Di | psi_Di>
# **Returned Matrix is transpose of LSEMatrix output**
###
CreateNLSEMatrix := function( numQubits )
    local diagrams, polygonOperators, stateVecs, matrix, row, i, j, normFactor, value;

    matrix := [];

    diagrams := NCPartitionsSet( [1..numQubits] );
    StableSort( diagrams, CompareDiagrams );
    polygonOperators := List( diagrams, WernerDiagram );
    stateVecs := List( diagrams, x -> CDofI( x, DopeyString( x ) ) );
    
    for i in [1..Length( stateVecs )] do
        row := [];
        
        for j in [1..Length( polygonOperators )] do
            normFactor := InnerProduct( stateVecs[ i ], ( polygonOperators[ i ] * stateVecs[ i ] ) );
            value := InnerProduct( stateVecs[ i ], ( polygonOperators[ j ] * stateVecs[ i ] ) );
            
            Add( row, value / normFactor );
        od;
        
        Add( matrix, row );
    od;

    return TransposedMat( matrix );
end;
