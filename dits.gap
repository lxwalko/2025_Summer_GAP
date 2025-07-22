#############################################
## Bits and Bitstrings / Dits and Ditstrings
#############################################

Print("Read dits.gap");

###
# REQUIRES "orbdim.gap"
# Takes a list of 0s and 1s such as
# [ 1, 0, 0, 1, 1]
# Returns the tensor product of
# ketzeros and ketones in given order
###
BitstringTensor := function( bitstring )
    local res, bit;

    res := [ [1] ];

    for bit in bitstring do

        if bit = 1 then
            res := KroneckerProduct( res, ketone );
        else
            res := KroneckerProduct( res, ketzero );
        fi;
    od;

    return normalize( res );
end;

###
# BitSize() returns log_2( Length( list ) )
# or log base 2 of the length of a list
###
BitSize := function( lst )
# works for state vectors also
  return LogInt( Length( lst ), 2 );
end;

###
# DitSize() returns log_d( Length( list ) )
# or log base d of the length of a list
###
DitSize := function( d, lst )
    return LogInt( Length( lst ), d );
end;

###
# Given a binary list, BitComplement() returns a binary list
# where 1s are 0s and 0s are 1s
###
BitComplement := function( binlist )
    local bit, complist;

    complist:=[];

    for bit in binlist do
        if bit = 0 then
            Add( complist, 1 );
        elif bit = 1 then
            Add( complist, 0 );
        else
            Error("BitComplement:  bad binlist");
        fi;
    od;

    return complist;
end;

###
# dec2binlist returns the binary expansion of a number in a given number of qubits
# For example:
# gap> dec2binlist( 5, 5 );
# [ 0, 0, 1, 0, 1 ]
#
# As 5 = 00101 in binary using 5 places
###
dec2binlist := function( dec, numqubits )
    local bit, binlist, remainder, place;

    binlist := [];
    remainder := dec;

    for bit in [1..numqubits] do
        place := 2^( numqubits-bit );

        if remainder >= place then
            Append( binlist, [ 1 ] );
            remainder := remainder - place;
        else
            Append( binlist, [0] );
        fi;
    od;

    return binlist;
end;

dec2base4list := function( dec, numqubits )
    local bit, b4list, remainder, place, digit;
    
    b4list := [];
    remainder := dec;
    
    for bit in [1..numqubits] do
        place := 4^( numqubits - bit );
        digit := Int( remainder / place );
        remainder := remainder - digit * place;
        Append( b4list, [ digit ] );
    od;
    
    return b4list;
end;

###
# dec2baseNlist() works in the same manner as dec2binlist()
# However, it returns answers in base N rather than base 2
# numdigits is length of the ditstring
# The call dec2baseNlist( 3, 4, 6 ) means:
#   Express 3 in a base 4 ditstring of length 6
#
# See description of dec2binlist() for an example
###
dec2baseNlist := function( dec, N, numdigits )
    local pos, list, remainder, place, digit;

    list := [];
    remainder := dec;

    for pos in [1..numdigits] do
        place := N^( numdigits - pos );
        digit := Int( remainder / place );
        remainder := remainder - (digit * place);

        Append( list, [ digit ] );
    od;

    return list;
end;

# Evaluates a bitstring in base 10
DecOfBinstring := function( binstring )
    local place, bit, dec;

    dec := 0;
    place := Length( binstring );

    for bit in binstring do
        place := place - 1;
        dec := dec + Int( [bit] ) * 2^place;
    od;

    return dec;
end;

CompositeBitlist := function( subsystem, list0, list1 )
    local i0, i1, bit, comp;

    i0:=1;
    i1:=1;
    comp:=[];

    for bit in subsystem do
        if bit = 0 then
            Append( comp, [ list0[ i0 ] ] );
            i0 := i0 + 1;
        else # if bit = 1
            Append( comp, [ list1[ i1 ] ] );
            i1 := i1 + 1;
        fi;
    od;

    return comp;
end;

OneBinList := function( k, n )
# return binlist of length n with 1 in position k, 0 elsewhere
    local bit,binlist;

    binlist:=[];

    for bit in [1..n] do
        if bit = k then
            Append( binlist, [ 1 ] );
        else
            Append( binlist, [ 0 ] );
        fi;
    od;

    return binlist;
end;

###
# Bin() returns the binary** sum of a vector
# Bin() does not require 1s and 0s. It works for any vector with numbers
# For example: Bin( [1, 2] ) = (2^(1) * 1) + (2^(0) * 2) = 4
###
Bin := function( vec )
    local place, dec;

    dec := 0;

    for place in [1..Length( vec )] do
        dec := dec + 2^(Length( vec ) - place) * vec[ place ];
    od;

    return dec;
end;

###
# BaseN() is a generalized form of Bin() for base N
###
BaseN := function( base, vec )
    local place, total;

    total := 0;

    for place in [1..Length( vec )] do
        total := total + base^(Length( vec ) - place) * vec[place];
    od;

    return total;
end;

###
# e() returns a list with length 2^num where the
# entry in the bin + 1 position is 1, all else 0
##
e := function( num, bin )
    local v,count;

    v := [ 0 ];
    for count in [1..num] do
        Append( v, v );
    od;

    v[ bin + 1 ] := 1;

    return v;
end;

eDit := function( dim, pow, pos )
    local res;

    res := List( [1..dim^pow], i -> 0 );

    res[ pos + 1 ] := 1;

    return res;
end;

ee := function( lst )
#    depends on e
    local n;

    n := Length( lst );

    return e( n, Bin( lst ) );
end;

eeDit := function( dim, lst )
    local numQudits;

    numQudits := Length( lst );

    return eDit( dim, numQudits, BaseN( dim, lst ) );
end;

bitflip := function( l, bitlist )
    local result, j;
    
    result := [];
    for j in [1..Length( bitlist )] do
        if j = l then
            Add( result, 1 - bitlist[ j ] );
        else
            Add( result, bitlist[ j ] );
        fi;
    od;
    
    return result;
end;

###
# Retuns a list of ditstrings that count from 0 to base^numDigits - 1
# Function is a binary counting algorithm expanded to all bases
###
BaseNCounting := function( base, numDigits )
    local i, value, ditStrings;

    ditStrings := [];

    for i in [0..(base^numDigits - 1)] do
        value := dec2baseNlist( i, base, numDigits );
        Add( ditStrings, value );
    od;

    return ditStrings;
end;
