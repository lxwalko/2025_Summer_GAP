# To start gap, type "gap" at shell prompt:
# $ gap

# To load orbit dimension functions:
# gap> Read("orbdim.gap");

# or just
# $ gap orbdim.gap

# Example of orbit dimension calculation:
# gap> ghz := ee([0,0,0]) + ee([1,1,1]);
# gap> OrbitDim(ghz);

# A five-qubit example:
# gap> fivedude := 3 * ee([0,0,1,1,1]) + 2 * ee([1,1,1,0,0]);
# gap> OrbitDim(fivedude);

######## Basic Utilities ########
# (ones that are used by other functions must come first)


###
# Comments made with this format were added in Summer 2025
###



Print("Read orbdim.gap");

i := E(4);

Real:=function(z)
  return (z+ComplexConjugate(z))/2;
end;
Imag:=function(z)
  return (z-ComplexConjugate(z))/(2*i);
end;

bracket:=function(a,b)
  return a*b - b*a;
end;

car:=function(lst)
  return lst[1];
end;

# Returns a list from 2 to Length( lst )
cdr := function( lst )
  return lst{ [2..Length( lst )] };
end;

cons:=function(a,lst)
  return Concatenation([a],lst);
end;

# Gap has a built-in map function called List, except arguments are reversed
# List([1,4,9],Sqrt) --> [ 1, 2, 3 ]

map:=function(f,lst)
  if IsEmpty(lst) then
    return [];
  else
    return cons(f(car(lst)),map(f,cdr(lst)));
  fi;
end;

fmap := f -> xs -> List(xs,f);

curry2 := f -> x -> y -> f(x,y);

flip := f -> x -> y -> f(y)(x);

compose := f -> g -> x -> f(g(x));

const := x -> y -> x;

replicate := num -> item -> map(const(item),[1..num]);

# recursive form
# zipWith := f -> xs -> function(ys)
#   if IsEmpty(xs) then
#     return [];
#   elif IsEmpty(ys) then
#     return [];
#   else
#     return Concatenation([f(xs[1])(ys[1])],zipWith(f)(cdr(xs))(cdr(ys)));
#   fi;
# end;

# non-recursive version for speed
# 2013 Aug 6
zipWith := f -> xs -> function(ys)
  local j,result;
  result := [];
  for j in [1..Minimum(Length(xs),Length(ys))] do
    Add(result,f(xs[j])(ys[j]));
  od;
  return result;
end;

# older recursive version (that I assume is slower and/or a memory hog)
# zip:=function(lst1,lst2)
#   local helper;
#   helper:=function(zipped,lst1,lst2)
#     if Length(lst1)=0 then
#       return zipped;
#     elif Length(lst2)=0 then
#       return zipped;
#     else
#       return helper(Concatenation(zipped,[[lst1[1],lst2[1]]]),
#                     cdr(lst1),cdr(lst2));
#     fi;
#   end;
#   return helper([],lst1,lst2);
# end;

# non-recursive version for speed
# 2013 Aug 5
zip := function(lst1,lst2)
  local j,result;
  result := [];
  for j in [1..Minimum(Length(lst1),Length(lst2))] do
    Add(result,[lst1[j],lst2[j]]);
  od;
  return result;
end;

SetMinus := xs -> ys -> Filtered(xs,x -> not (x in ys));

prod := x -> y -> x*y;

foldr := f -> v ->
  function(lst)
    if IsEmpty(lst) then
      return v;
    else
      return f(car(lst))(foldr(f)(v)(cdr(lst)));
    fi;
  end;

composeList := foldr(compose)(IdFunc);

if3 := function(pred,t,f)
  if pred then
    return t;
  else
    return f;
  fi;
end;

if12 := function(pred)
  local rest;
  rest := function(t,f)
    if pred then
      return t;
    else
      return f;
    fi;
  end;
  return rest;
end;

If := pred -> curry2(if12(pred));

And := function(bs)
  if IsEmpty(bs) then
    return true;
  else
    return (car(bs) and And(cdr(bs)));
  fi;
end;

ReplaceFirstOccurrence:=function(lst,target,replacement)
  if IsEmpty(lst) then
    return [];
  else
    if lst[1] = target then
      return cons(replacement,cdr(lst));
    else
      return cons(lst[1],ReplaceFirstOccurrence(cdr(lst),target,replacement));
    fi;
  fi;
end;

ReplaceAllOccurrences:=function(lst,target,replacement)
  if IsEmpty(lst) then
    return [];
  else
    if lst[1] = target then
      return cons(replacement,
                  ReplaceAllOccurrences(cdr(lst),target,replacement));
    else
      return cons(lst[1],ReplaceAllOccurrences(cdr(lst),target,replacement));
    fi;
  fi;
end;

IsSimpleList:=function(lst)
  return IsList(lst) and (IsEmpty(lst) or not IsList(lst[1]));
end;

MapRecursiveOntoSimpleList:=function(lst,func)
# func takes a simple list as argument
  local funcRecursive;
  funcRecursive:=function(lst)
    if not IsList(lst) then
      Display("Error:  MapRecursiveOntoSimpleList:  lst not a list");
      return 0;
    elif IsSimpleList(lst) then
      return func(lst);
    else
      return List(lst,funcRecursive);
    fi;
  end;
  return funcRecursive(lst);
end;

filter:=function(lst,pred)
  if IsEmpty(lst) then
    return lst;
  elif pred(lst[1]) then
    return cons(lst[1],filter(cdr(lst),pred));
  else
    return filter(cdr(lst),pred);
  fi;
end;

Filter := flip(curry2(filter));

# imperative filter
FilterImp := pred -> function(xs)
  local result,x;
  result := [];
  for x in xs do
    if pred(x) then Add(result,x);
    fi;
  od;
  return result;
end;

sort:=function(lst)
  local LessThanHead,GrEqHead;
  if IsEmpty(lst) then
    return lst;
  else
    LessThanHead:=function(num)
      if num < lst[1] then
        return true;
      else
        return false;
      fi;
    end;
    GrEqHead:=function(num)
      if num >= lst[1] then
        return true;
      else
        return false;
      fi;
    end;
    return Concatenation(sort(filter(cdr(lst),LessThanHead)),
                         [lst[1]],
                         sort(filter(cdr(lst),GrEqHead)));
  fi;
end;

SortingBySize := function(thing1,thing2)
    if IsList(thing1) and IsList(thing2) then
        if Length(thing1) < Length(thing2) then
            return true;
        elif Length(thing1) > Length(thing2) then
            return false;
        else
            return thing1 < thing2;
        fi;
    else
        return thing1 < thing2;
    fi;
end;

SortBySize := function(xss)
    local yss;
    yss := StructuralCopy(xss);
    Sort(yss,SortingBySize);
    return yss;
end;

Count:=function(num,bitlist)
  local sum,bit;
  sum:=0;
  for bit in bitlist do
    if bit=num then
      sum:=sum+1;
    fi;
  od;
  return sum;
end;

# Count2 is a functional version of Count
Count2:=function(num,bitlist)
  if IsEmpty(bitlist) then
    return 0;
  else
    if bitlist[1] = num then
      return 1 + Count2(num,cdr(bitlist));
    else
      return Count2(num,cdr(bitlist));
    fi;
  fi;
end;

CountNonZero := compose(Sum)(fmap(x -> If(x <> 0)(1)(0)));

RandList:=function(len)
  if len=0 then
    return [];
  else
    return cons(Random(Integers)+Random(Integers)*i,RandList(len-1));
  fi;
end;

RandomQubits:=function(n)
  return RandList(2^n);
end;

display:=function(data)
  local info,row,max,col;
  max:=[1,1,1,1,1,1,1,1];
  for row in data do
    for col in [1..Length(row)] do
      if Length(String(row[col])) > max[col] then
        max[col]:=Length(String(row[col]));
      fi;
    od;
  od;
  for row in data do
    for col in [1..Length(row)] do
      Print(String(row[col],max[col]+1));
    od;
    Print("\n");
  od;
end;

# functional version of OneBinList
OneBinListF:=function(k,n)
# return binlist of length n with 1 in position k, 0 elsewhere
  local bit,binlist,worker;
  worker:=function(k,n,lst)
    if n=0 then
      return lst;
    elif k=1 then
      return cons(1,worker(k-1,n-1,lst));
    else
      return cons(0,worker(k-1,n-1,lst));
    fi;
  end;
  return worker(k,n,[]);
end;

# A QubitSet is a subset of 1..n
# A BinarySubsystem is a list of n elements, each 0 or 1.

QubitSet2BinarySubsystem := function( n, bitset )
    if IsEmpty( bitset ) then
      return OneBinListF( n+1, n ); # list of n zeros
    else
      return OneBinListF( bitset[ 1 ], n ) + QubitSet2BinarySubsystem( n, cdr( bitset ) );
    fi;
end;

# gap> QubitSet2BinarySubsystem(3,[2,3]);
# [ 0, 1, 1 ]

# does FlattenOneLevel do anything different from Concatenation?
# FlattenOneLevel:=function(ll)
#   if IsEmpty(ll) then
#     return ll;
#   else
#     return Concatenation(ll[1],FlattenOneLevel(cdr(ll)));
#   fi;
# end;
# I think no, so:
FlattenOneLevel:=Concatenation;

Flatten2LL:=function(lst)
  if not IsList(lst) or IsEmpty(lst) then
    return lst;
  elif not IsList(lst[1]) or IsEmpty(lst[1]) then
    return lst;
  elif not IsList(lst[1][1]) then
    return lst;
  else
    return Flatten2LL(FlattenOneLevel(lst));
  fi;
end;

# Doesn't give quite the order I want.
# CombinationsBySize :=
#   set ->
#   FlattenOneLevel(List([0..Length(set)],size -> Combinations(set,size)));

CombinationsBySize := set -> SortBySize(Combinations(set));

SubsetList:=function(n)
  local combos;
  combos:=function(size)
    return Combinations([1..n],size);
  end;
  return FlattenOneLevel(List([0..n],combos));
end;
# Example:
# gap> SubsetList(3);
# [ [  ], [ 1 ], [ 2 ], [ 3 ], [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 1, 2, 3 ] ]

flatten:=function(mat)
  local vec,row;

  vec:=[];

  for row in mat do
    Add(vec,row);
  od;

  return vec;
end;

switchplace:=function(list,m,n)
  local newlist,holdthis,num;
  newlist:=[];
  for num in [1..Length(list)] do
    if num=m then
      Append(newlist,[list[n]]);
    elif num=n then
      Append(newlist,[list[m]]);
    else
      Append(newlist,[list[num]]);
    fi;
  od;
  return newlist;
end;

ComponentWithBitString := function(bitString,psi)
    if Length(bitString) = BitSize(psi) then
        return psi[Bin(bitString)+1];
    else
        Error("ComponentWithBitString: size mismatch");
    fi;
end;

AppendMat:=function(mat1,mat2)
# this function changes mat1
  local row,rownum;
  if Length(mat1) <> Length(mat2) then
    Display("Error:  AppendMat:  number of rows differs");
    return 0;
  fi;
  rownum:=0;
  for row in mat1 do
    rownum:=rownum+1;
    Append(row,mat2[rownum]);
  od;
end;

complement:=function(ns,subspace)
  local a,row;
  a:=[];
  for row in [1..Length(subspace)] do
    Add(a,SolutionMat(ns,subspace[row]));
  od;
  Display(a);
  return NullspaceMat(TransposedMat(a))*ns;
end;

square:=function(x)
  return x*x;
end;

FullSelfContraction:=function(tensor)
  if not IsList(tensor) then
    return square(tensor);
  else
    return Sum(map(FullSelfContraction,tensor));
  fi;
end;

######## end of Basic Utilities ########

######## Basic Objects ########

I:=[[1,0],[0,1]];
A:=[[i,0],[0,-i]];
B:=[[0,1],[-1,0]];
C:=[[0,i],[i,0]];
XX:=[A,B,C];

sx:=-i * C;
sy:=-i * B;
sz:=-i * A;

# here are some matrices useful for producing new state vectors
EA4:=[[1+i,0],[0,1-i]]/Sqrt(2);  # this is exp(pi A/4) = exp(i pi sz/4) = sqrt(A)
EB4:=[[1,1],[-1,1]]/Sqrt(2);     # this is exp(pi B/4) = exp(i pi sy/4) = sqrt(B)
EC4:=[[1,i],[i,1]]/Sqrt(2);      # this is exp(pi C/4) = exp(i pi sx/4) = sqrt(C)

# Hadamard gate
uhad := 1/Sqrt(2) * [[1,1],[1,-1]];

# gap> uhad * uhad = I;
# true
# gap> uhad * sx * uhad = sz;
# true

# Phase gate
uS := [[1,0],[0,i]];

# Swap gate (2-qubit)
swap := [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]];

# cos(2 pi m/n) = Real(E(n)^m)
cos:=function(deg)
  return Real(E(360*DenominatorRat(deg))^NumeratorRat(deg));
end;

# sin(2 pi m/n) = Imag(E(n)^m)
sin:=function(deg)
  return Imag(E(360*DenominatorRat(deg))^NumeratorRat(deg));
end;

x:=Indeterminate(Rationals,"x");
y:=Indeterminate(Rationals,"y");
z:=Indeterminate(Rationals,"z");

######## end of Basic Objects ########

########### Intermediate Utilities ####################

DMn:=function(cvec)
  return normalizeDM(DM(cvec));
end;

#DMn:=function(cvec)
#  return DM(normalize(cvec));
#end;

GrElt:=function(num,den,vec)
# this function returns exp(i 2 pi num/den normalize(vec) . sigma)
  local unitvec;
  unitvec:=normalize(vec);
  return Real(E(den)^num) * I + i * Imag(E(den)^num)
    * (unitvec[1] * sx + unitvec[2] * sy + unitvec[3] * sz);
end;

GeneralizedPauli := function( dim, a, b )
    local X, Z, i, j, Xmat, Zmat;

    # Shift matrix X
    Xmat := NullMat(dim, dim);
    for i in [1..dim] do
        Xmat[i][ (i mod dim) + 1 ] := 1;
    od;

    # Diagonal matrix Z
    Zmat := NullMat(dim, dim);
    for i in [1..dim] do
        Zmat[i][i] := E(dim)^(i - 1);
    od;

    return Xmat^a * Zmat^b;
end;

PauliQudits := function( dim, num )
    local a, b;

    a := Quotient(num, dim);   # First index
    b := num mod dim;          # Second index

    return GeneralizedPauli(dim, a, b);
end;

Pauli := function( num )
#    depends on Basic Objects
  if num = 0 then
    return I;
  elif num = 1 then
    return sx;
  elif num = 2 then
    return sy;
  elif num = 3 then
    return sz;
  else
    Error("Error:  Pauli:  bad index", num);
  fi;
end;

HalfPauli := function( num )
    return Pauli( num ) / 2;
end;

HalfPauliQudits := function( dim, num )
    return PauliQudits( dim, num ) / 2;
end;

realvec:=function(cvec)
#    depends on Real,Imag
  local rvec,c;
  rvec:=[];
  for c in cvec do
    Append(rvec,[Real(c),Imag(c)]);
  od;
  return rvec;
end;

PartialTrace := function( dm, subsystem )
#    depends on Count
# subsystem is a binlist
    local rdm, bit, tracesize, row, col, sumi, keepnbits, tracenbits;
    if Length( subsystem ) <> Count( 0, subsystem ) + Count( 1, subsystem ) then
        return "Error:  PartialTrace:  subsystem not a binlist";
    elif LogInt( Length( dm ), 2 ) <> Length( subsystem ) then
        return "Error:  PartialTrace:  sizes don't match";
    else
        rdm := [ [ 0 ] ];
        for bit in subsystem do
            if bit = 0 then
	            rdm := KroneckerProduct( rdm, [ [0,0], [0,0] ] );
            fi;
        od;
        tracenbits := Count( 1, subsystem );
        tracesize := 2^tracenbits;
        keepnbits := Length( subsystem ) - tracenbits;
        for row in [1..Length( rdm )] do
            for col in [1..Length( rdm )] do
	            for sumi in [1..tracesize] do
	                rdm[ row ][ col ] := rdm[ row ][ col ]
	                + dm[Bin(CompositeBitlist(subsystem,dec2binlist(row-1,keepnbits),dec2binlist(sumi-1,tracenbits)))+1]
		            [Bin(CompositeBitlist(subsystem,dec2binlist(col-1,keepnbits),dec2binlist(sumi-1,tracenbits)))+1];
	            od;
            od;
        od;

      return rdm;
    fi;
end;

# General form of PartialTrace()
PartialTraceQudits := function( dim, densityMatrix, targetQudits )
    local reducedMatrix, dit, traceSize, row, col, sumi, keepNDits, traceNDits;

    if Length( targetQudits ) <> Count( 0, targetQudits ) + Count( 1, targetQudits ) then
        Error( "Error: PartialTraceQudits: targetQudits not a binary list" );
    elif LogInt( Length( densityMatrix ), dim ) <> Length( targetQudits ) then
        Error( "Error: PartialTraceQudits: densityMatrix and targetQudits sizes don't match" );
    else

        traceNDits := Count( 1, targetQudits );
        traceSize := dim^traceNDits;
        keepNDits := Length( targetQudits ) - traceNDits;

        reducedMatrix := NullMat( dim^keepNDits, dim^keepNDits );

        for row in [1..Length( reducedMatrix )] do
            for col in [1..Length( reducedMatrix )] do
                for sumi in [1..traceSize] do
                    reducedMatrix[ row ][ col ] := reducedMatrix[ row ][ col ] +
                        densityMatrix[
                        BaseN( dim, CompositeBitlist( targetQudits, dec2baseNlist( row - 1, dim, keepNDits ),
                          dec2baseNlist( sumi - 1, dim, traceNDits ))) + 1 ]
                        [ BaseN( dim, CompositeBitlist( targetQudits, dec2baseNlist( col - 1, dim, keepNDits ),
                          dec2baseNlist( sumi - 1, dim, traceNDits ))) + 1 ];
                od;
            od;
        od;

        return reducedMatrix;
    fi;
end;

# gap> PartialTrace(Kron([DMn(e0+e1),DMn(e0)]),[0,1]);
# [ [ 1/2, 1/2 ], [ 1/2, 1/2 ] ]
# gap> PartialTrace(Kron([DMn(e0+e1),DMn(e0)]),[1,0]);
# [ [ 1, 0 ], [ 0, 0 ] ]
# gap> PartialTrace(Kron([DMn(e0+e1),DMn(e0)]),[1,1]);
# [ [ 1 ] ]
# gap> PartialTrace(Kron([DMn(e0+e1),DMn(e0)]),[0,0]);
# [ [ 1/2, 0, 1/2, 0 ], [ 0, 0, 0, 0 ], [ 1/2, 0, 1/2, 0 ], [ 0, 0, 0, 0 ] ]

rdm:=function(qubitset,dm)
  local n;
  n:=BitSize(dm);
  return PartialTrace(dm,QubitSet2BinarySubsystem(n,Difference([1..n],qubitset)));
end;

TensorCoeffList:=function(bitlist)
  local ReplaceExpand,MultReplaceExpand;
  MultReplaceExpand:=function(n,lst)
    if n=0 then
      return lst;
    else
      return MultReplaceExpand(n-1,ReplaceExpand(lst));
    fi;
  end;
  ReplaceExpand:=function(lst)
    if IsList(lst) then
      if not IsEmpty(lst) and IsList(lst[1]) then
        return List(lst,ReplaceExpand);
      else
        return [ReplaceFirstOccurrence(lst,-1,1),
                ReplaceFirstOccurrence(lst,-1,2),
                ReplaceFirstOccurrence(lst,-1,3)];
      fi;
    else
      Display("Error:  ReplaceExpand:  argument not a list");
      return 0;
    fi;
  end;
  return MultReplaceExpand(Count2(1,bitlist),
                           ReplaceAllOccurrences(bitlist,1,-1));
end;

FlatTensorCoeffList:=function(bitlist)
  local result;
  result:=TensorCoeffList(bitlist);
  if not IsList(result) or IsEmpty(result) then
    Display("Error:  FlatTensorCoeffList:  result not a list");
    return 0;
  elif not IsList(result[1]) then
    return [result];
  else
    return Flatten2LL(result);
  fi;
end;

# coefflist is a list of numbers
# dm is a density matrix
PauliCoeff := function( coefflist, dm )
#    depends on Kron,Pauli
    return TraceMat( Kron( List( coefflist, Pauli ) ) * dm );
end;

PauliCoeffQudits := function( dim, coefflist, densityMatrix )
    return TraceMat( Kron( List( coefflist, x->PauliQudits( dim, x ) ) ) * densityMatrix );
end;

PauliIndexList:=function(n)
# n = nbits
  local listlength,list,term,bit,listelt,coefflist,counter;
  counter:=0;
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for listlength in [0..n] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
        od;
        Print(counter,"  ",coefflist,"\n");
        counter:=counter+1;
      od;
    od;
  od;
end;

LexiPauliNTupleList:=function(n)
  local helper;
  helper:=function(j)
    return [0..3];
  end;
  return Cartesian(List([1..n],helper));
end;

# gap> LexiPauliNTupleList(2);
# [ [ 0, 0 ], [ 0, 1 ], [ 0, 2 ], [ 0, 3 ], [ 1, 0 ], [ 1, 1 ], [ 1, 2 ],
#   [ 1, 3 ], [ 2, 0 ], [ 2, 1 ], [ 2, 2 ], [ 2, 3 ], [ 3, 0 ], [ 3, 1 ],
#   [ 3, 2 ], [ 3, 3 ] ]

PauliNTupleList:=function(n)
  local qs2bs,helper,maphelper;
  qs2bs:=function(lst)
    return QubitSet2BinarySubsystem(n,lst);
  end;
  helper:=function(j)
    if j=0 then
      return [0];
    else
      return [1..3];
    fi;
  end;
  maphelper:=function(lst)
    return Cartesian(List(lst,helper));
  end;
  return Concatenation(List(List(SubsetList(n),qs2bs),maphelper));
end;

# gap> PauliNTupleList(2);
# [ [ 0, 0 ], [ 1, 0 ], [ 2, 0 ], [ 3, 0 ], [ 0, 1 ], [ 0, 2 ], [ 0, 3 ],
#   [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 1 ], [ 2, 2 ], [ 2, 3 ], [ 3, 1 ],
#   [ 3, 2 ], [ 3, 3 ] ]

PauliGroup:=function(n)
  local addphase;
  addphase:=function(pnt)
    return [cons(1,pnt),cons(i,pnt),cons(-1,pnt),cons(-i,pnt)];
  end;
  return FlattenOneLevel(List(PauliNTupleList(n),addphase));
end;

# PauliGroup(6) runs out of memory on T400.
# PauliGroup(5) is ok.
# But then Length(PauliGroup(6)) correctly reported 16384,
# and then PauliGroup(6) returned a huge list.  Hmmmm.

PauliGroupNoPhase:=PauliNTupleList;

DMTensor:=function(dm,subsystem)
  local pauliC;
  pauliC:=function(coefflist)
    return PauliCoeff(coefflist,dm);
  end;
  return MapRecursiveOntoSimpleList(TensorCoeffList(subsystem),pauliC);
end;

# gap> DMTensor(DMn(singlet),[1,1]);
# [ [ -1, 0, 0 ], [ 0, -1, 0 ], [ 0, 0, -1 ] ]

DMTensors:=function(dm)
  local func,n;
  n:=BitSize(dm);
  func:=function(subsetlst)
    return [subsetlst,DMTensor(dm,QubitSet2BinarySubsystem(n,subsetlst))];
  end;
  return List(SubsetList(n),func);
end;

# gap> DMTensors(DMn(singlet));
# [ [ [  ], 1 ], [ [ 1 ], [ 0, 0, 0 ] ], [ [ 2 ], [ 0, 0, 0 ] ],
#   [ [ 1, 2 ], [ [ -1, 0, 0 ], [ 0, -1, 0 ], [ 0, 0, -1 ] ] ] ]

DMContracts:=function(dm)
  local contract;
  contract:=function(pair)
    return [pair[1],FullSelfContraction(pair[2])];
  end;
  return map(contract,DMTensors(dm));
end;

DisplayListOfPairs:=function(lstofpairs)
  if not IsEmpty(lstofpairs) then
    Display(lstofpairs[1][1]);
    Display(lstofpairs[1][2]);
    Display(" ");
    DisplayListOfPairs(cdr(lstofpairs));
  fi;
end;

ShowTensors:=function(dm)
  DisplayListOfPairs(DMTensors(dm));
end;

# MakePIList is deprecated.  Use LieAlgList instead.
MakePIList:=function(n)
# n = nbits
  local listlength,list,term,bit,listelt,coefflist,result;
  result:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for listlength in [0..n] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
#          Print("listelt = ",listelt," list[listelt] = ",list[listelt],
#                " coefflist[list[listelt]] = ",coefflist[list[listelt]],"\n");
#          Print("coefflist = ",coefflist,"\n");
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
#          Print("coefflist = ",coefflist,"\n");
        od;
#        Print("final coefflist = ",coefflist,"\n");
        Add(result,StructuralCopy(coefflist));
#        Print("final coefflist = ",coefflist,"\n");
#        Print("result = ",result,"\n");
      od;
    od;
  od;
  return result;
end;

LieAlgList:=function(n,type)
# n = nbits
  local listlength,list,term,bit,listelt,coefflist,result,startll,endll;
  result:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  if type="u" then
    startll:=0;
    endll:=n;
  elif type="su" then
    startll:=1;
    endll:=n;
  elif type="lu" then
    startll:=1;
    endll:=1;
  else
    Display("Error:  LieAlgList:  bad type");
    return 0;
  fi;
  for listlength in [startll..endll] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
#          Print("listelt = ",listelt," list[listelt] = ",list[listelt],
#                " coefflist[list[listelt]] = ",coefflist[list[listelt]],"\n");
#          Print("coefflist = ",coefflist,"\n");
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
#          Print("coefflist = ",coefflist,"\n");
        od;
#        Print("final coefflist = ",coefflist,"\n");
        Add(result,StructuralCopy(coefflist));
#        Print("final coefflist = ",coefflist,"\n");
#        Print("result = ",result,"\n");
      od;
    od;
  od;
  return result;
end;

# functional version of LieAlgList
LieAlgListF:=function(n,startll,endll)
  local combos,qs2bs;
  combos:=function(size)
    return Combinations([1..n],size);
  end;
  qs2bs:=function(bitset)
    return QubitSet2BinarySubsystem(n,bitset);
  end;
#  return Flatten2LL(List(List(List([startll..endll],
#                                   combos),
#                              qs2bs),
#                         TensorCoeffList));
  return Flatten2LL(List(List(Flatten2LL(List([startll..endll],
                                         combos)),
                              qs2bs),
                         FlatTensorCoeffList));
end;

fillPauliForm:=function(b4list,pf,dm)
  local pi;
  if IsList(pf) then
#    Display(["Starting non-final loop",b4list,pf]);
    for pi in [0..3] do
      pf[pi+1]:=fillPauliForm(Concatenation(b4list,[pi]),pf[pi+1],dm);
    od;
#    Display(["Finishing non-final loop",b4list,pf]);
  else
    pf:=PauliCoeff(b4list,dm);
#    Display(["Final loop",b4list,pf]);
  fi;
  return pf;
end;

PauliForm:=function(dm)
  local pf,nbits,bit;
  pf:=0;
  nbits:=LogInt(Length(dm),2);
  for bit in [1..nbits] do
    pf:=[StructuralCopy(pf),StructuralCopy(pf),StructuralCopy(pf),StructuralCopy(pf)];
  od;
  fillPauliForm([],pf,dm);
  return pf;
end;

PauliListForm:=function(dm)
  local plf,nbits,bit,pil,coefflist;
  nbits:=LogInt(Length(dm),2);
#  pil:=MakePIList(nbits);
  pil:=LieAlgList(nbits,"u");
  plf:=[];
  for coefflist in pil do
    Add(plf,PauliCoeff(coefflist,dm));
  od;
  return plf;
end;

# The following function, PauliListForm2, is intended to do the same
# job as PauliListForm, but coded in a more functional way.

PauliListForm2:=function(dm)
  local pil2coeff;
  pil2coeff:=function(coefflist)
    return PauliCoeff(coefflist,dm);
  end;
  return List(LieAlgList(BitSize(dm),"u"),pil2coeff);
end;

PauliListForm3:=function(dm)
  local pil2coeff;
  pil2coeff:=function(coefflist)
    return [coefflist,PauliCoeff(coefflist,dm)];
  end;
  return List(LieAlgList(BitSize(dm),"u"),pil2coeff);
end;

DisplayPLF:=function(dm)
  display(PauliListForm3(dm));
end;

DMFromPLF:=function(plf)
  local n,dm,count,pil,b4list;
  n:=LogInt(Length(plf),4);
  dm:=0;
  count:=0;
  pil:=LieAlgList(n,"u");
  for b4list in pil do
    count:=count+1;
    dm:=dm+plf[count]*Kron(List(b4list,HalfPauli));
  od;
  return dm;
end;

helperDMFromPF:=function(list,pf)
  if IsList(pf) then
    return helperDMFromPF(Concatenation(list,[I/2]),pf[1])
         + helperDMFromPF(Concatenation(list,[sx/2]),pf[2])
         + helperDMFromPF(Concatenation(list,[sy/2]),pf[3])
         + helperDMFromPF(Concatenation(list,[sz/2]),pf[4]);
  else
    return pf * Kron(list);
  fi;
end;

DMFromPF:=function(pf)
  return helperDMFromPF([],pf);
end;

fillPTPF:=function(b4list,pf)
  local pi;
  if IsList(pf) then
#    Display(["Starting non-final loop",b4list,pf]);
    for pi in [1..3] do
      pf[pi+1]:=fillPTPF(Concatenation(b4list,[pi]),pf[pi+1]);
    od;
#    Display(["Finishing non-final loop",b4list,pf]);
  else
    pf:=0;
#    Display(["Final loop",b4list,pf]);
  fi;
  return pf;
end;

DMFromPLF3 := compose(Sum)(fmap(x -> x[2] * Kron(fmap(HalfPauli)(x[1]))));

# density matrix -> density matrix
# PauliWeightProject := w -> rho -> DMFromPLF3(Filter(x -> Count(0,x[1]) = BitSize(rho)-w)(PauliListForm3(rho)));
# This should be the same as PauliWeight.
# PauliWeight was written earlier and looks more cleanly written.

Pad := function( binlist, dm )
    # for every 1 in binlist, insert I/2 in the corresponding slot
    local n, term, b4list, coeff, sum, bigb4list, b4listcounter, bit;

    n := BitSize( dm );

    if Count( 0, binlist ) <> n then
        Display("Error:  Pad:  size mismatch");
        return 0;
    fi;

    sum := 0;
    for term in [0..4^n-1] do
        b4list := dec2baseNlist( term, 4, n );
        coeff := PauliCoeff( b4list, dm );
        b4listcounter := 0;
        bigb4list := [];

        for bit in binlist do
            if bit = 0 then
                b4listcounter := b4listcounter + 1;
                Add( bigb4list, b4list[ b4listcounter ] );
            elif bit = 1 then
                Add( bigb4list, 0 );
            else
                Display("Error:  Pad:  bad binlist");
            fi;
        od;
        sum := sum + coeff * Kron( List( bigb4list, HalfPauli ) );
    od;

    return sum;
end;

###
# Debugging tool, checks if matrix is square
###
IsSquareMat := function( matrix )
    if Length( matrix ) = Length( matrix[ 1 ] ) then
        return true;
    fi;

    return false;
end;

OrderedKron:=function(binlist,dm0,dm1)
  local n0,n1;
  n0:=BitSize(dm0);
  n1:=BitSize(dm1);
  if Count(0,binlist)=n0 and Count(1,binlist)=n1 and Length(binlist)=n0+n1 then
    return 2^n0 * 2^n1 * Pad(binlist,dm0) * Pad(BitComplement(binlist),dm1);
  else
    Error("OrderedKron:  wrong sizes.");
  fi;
end;

IsPartitionOfN := lsts -> SortedList(Concatenation(lsts)) = [1..Length(Concatenation(lsts))];
###
# dmlst is a list of density matrices
# partition is a list of lists that describes an NCP diagram
###
OrderedKronList := function( partition, dmlst )
    local n, i, result;

    # raises an error if partition and dmlst have different lengths
    Assert( 0, Length( partition ) = Length( dmlst ) );


    for i in [1..Length( partition )] do
      # checks that each gon in partition and dmlst are the same size
      Assert( 0, Length( partition[ i ] ) = BitSize( dmlst[ i ] ) );
    od;

    Assert( 0, IsPartitionOfN( partition ) );
    n := Length( Concatenation( partition ) );

    result := IdentityMat( 2^n );

    for i in [1..Length( partition )] do
      result := result * 2^(n-Length( partition[ i ] ))
                       * Pad( BitComplement( QubitSet2BinarySubsystem( n, partition[ i ] )), dmlst[ i ]);


    od;

    return result;
end;

# gap> OrderedKronList([[1,3],[2]],[DMn(singlet),DMn(e0)]);
# [ [ 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 1/2, 0, 0, -1/2, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0 ],
#   [ 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, -1/2, 0, 0, 1/2, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0 ],
#   [ 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0 ] ]

ordered_kron := function(binlist,psi0,psi1)
    local n0,n1,psi,j,l,II,J0,J1;
    n0 := BitSize(psi0);
    n1 := BitSize(psi1);
    if Count(0,binlist) = n0 and Count(1,binlist) = n1 and
       Length(binlist) = n0+n1 then
        psi:=[];
        for j in [1..2^(n0+n1)] do
            II := dec2binlist(j-1,n0+n1);
            J0 := [];
            J1 := [];
            for l in [1..n0+n1] do
                if binlist[l] = 0 then
                    Add(J0,II[l]);
                else
                    Add(J1,II[l]);
                fi;
            od;
            Add(psi,psi0[Bin(J0)+1] * psi1[Bin(J1)+1]);
        od;
        return psi;
    else
        Error("ordered_kron:  wrong sizes.");
    fi;
end;

OrderedConcatList:=function(partition,lst)
  local i;
  Assert(0,Length(partition)=Length(lst));
  for i in [1..Length(partition)] do
    Assert(0,Length(partition[i])=Length(lst[i]));
  od;
  Assert(0,IsPartitionOfN(partition));
  return List(SortedList(Concatenation(ListN(partition,lst,zip))),pr -> pr[2]);
end;

PartialTracePF:=function(pf)
#
#  Caution!:  PartialTracePF does not always produce (the Pauli form of) a legitimate
#             density matrix.
#
  local newpf;
  newpf:=StructuralCopy(pf);
  return fillPTPF([],newpf);
end;

SwitchBitsCvec:=function(cvec,m,n)
#    depends on switchplace,dec2binlist,Bin
  local nbits,item,newvec;
  nbits:=LogInt(Length(cvec),2);
  newvec:=[];
  for item in [1..Length(cvec)] do
    Append(newvec,[cvec[Bin(switchplace(dec2binlist(item-1,nbits),m,n))+1]]);
  od;
  return newvec;
end;

###
# perm is a GAP permutation
# GAP permutations are written in cycle notation
###
PermuteList := function( lst, perm )
    local newlst, i;
    newlst := [];

    for i in [1..Length( lst )] do
        Append( newlst, [lst[ i^( perm^( -1 ))]]);
    od;

    return newlst;
end;

###
# cvec is a vector such as [1,0,0,0,0,0,0,0]
# Returns a vector with qubits permuted according to perm
# Perm is a gap permutation (Written in cycle notation)
###
PermuteQubitsPsi := function( cvec, perm )
    local n, newvec, item;

    n := BitSize( cvec );
    newvec := [];

    for item in [1..Length( cvec )] do
        Append( newvec, [cvec [Bin( PermuteList( dec2binlist(item - 1, n), perm^(-1) ) ) + 1 ]]);
    od;

    return newvec;
end;

###
# General version of PermuteQubitsPsi for local dimension dim
###
PermuteQudits := function( dim, cvec, perm )
    local n, newvec, item;

    n := DitSize( dim, cvec );
    newvec := [];

    for item in [1..Length( cvec )] do
        Append( newvec, [cvec [BaseN( dim, PermuteList( dec2baseNlist( item - 1, dim, n ), perm^( -1 ) ) ) + 1 ]] );
    od;

    return newvec;
end;

# we should have
# SwitchBitsCvec(cvec,m,n) = PermuteQubitsPsi(cvec,(m,n))

# gap> DMTensors(DMn(kron(singlet,e0)));
# [ [ [ ], 1 ], [ [ 1 ], [ 0, 0, 0 ] ], [ [ 2 ], [ 0, 0, 0 ] ], [ [ 3 ], [ 0, 0, 1 ] ],
#   [ [ 1, 2 ], [ [ -1, 0, 0 ], [ 0, -1, 0 ], [ 0, 0, -1] ] ],
#   [ [ 1, 3 ], [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0] ] ],
#   [ [ 2, 3 ], [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0] ] ],
#   [ [ 1, 2, 3 ], [ [ [ 0, 0, -1 ], [ 0, 0,  0 ], [ 0, 0,  0 ] ],
#                    [ [ 0, 0,  0 ], [ 0, 0, -1 ], [ 0, 0,  0 ] ],
#                    [ [ 0, 0,  0 ], [ 0, 0,  0 ], [ 0, 0, -1 ] ] ] ] ]
# gap> DMTensors(DMn(PermuteQubitsPsi(kron(singlet,e0),(1,2,3))));
# [ [ [  ], 1 ], [ [ 1 ], [ 0, 0, 1 ] ], [ [ 2 ], [ 0, 0, 0 ] ], [ [ 3 ], [ 0, 0, 0 ] ],
#   [ [ 1, 2 ], [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0 ] ] ],
#   [ [ 1, 3 ], [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0 ] ] ],
#   [ [ 2, 3 ], [ [ -1, 0, 0 ], [ 0, -1, 0 ], [ 0, 0, -1 ] ] ],
#   [ [ 1, 2, 3 ], [ [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0 ] ],
#                    [ [  0, 0, 0 ], [ 0,  0, 0 ], [ 0, 0,  0 ] ],
#                    [ [ -1, 0, 0 ], [ 0, -1, 0 ], [ 0, 0, -1 ] ] ] ] ]

# gap> PermuteQubitsPsi(ee([0,0,1,1,1,1]),(2,3,4)) = ee([0,1,0,1,1,1]);
# true
# gap> (2,3,4) = (2,3)*(2,4);
# true
# gap> (2,3,4) = (2,3)*(4,2);
# true
# gap> (2,3,4) = (2,4)*(2,3);
# false
# How I think:  (2,3,4) means send 2 -> 3, 3 -> 4, 4 -> 2
#       (a,b,c,d) === (a,b,c,d)
# (2,3)     |             |
#           V             |
#       (a,c,b,d)         |     (2,3,4)
# (2,4)     |             |
#           V             V
#       (a,d,b,c) === (a,d,b,c)

PermStabPsi := cvec -> Stabilizer(SymmetricGroup(BitSize(cvec)),cvec,PermuteQubitsPsi);

# PermuteQubitsRho:=function(rho,perm)
#   local n,row,col,newrho,newrow,newcol;
#   n:=BitSize(rho);
# #  Print(n,"\n");
#   newrho:=StructuralCopy(rho);
#   for row in [1..Length(rho)] do
#     for col in [1..Length(rho)] do
#       newrow:=Bin(PermuteList(dec2binlist(row-1,n),perm))+1;
#       newcol:=Bin(PermuteList(dec2binlist(col-1,n),perm))+1;
# #      Print(newrow,newcol,"\n");
#       newrho[newrow][newcol]:=rho[row][col];
#     od;
#   od;
#   return newrho;
# end;

PermuteQubitsRho:=function(rho,perm)
  local n,row,col,newrho,newrow;
  n:=BitSize(rho);
#  Print(n,"\n");
  newrho:=[];
  for row in [1..Length(rho)] do
    newrow:=[];
    for col in [1..Length(rho)] do
      Append(newrow,[rho[Bin(PermuteList(dec2binlist(row-1,n),perm^(-1)))+1]
                        [Bin(PermuteList(dec2binlist(col-1,n),perm^(-1)))+1]]);
    od;
    Append(newrho,[newrow]);
  od;
  return newrho;
end;

PermStabRho := rho -> Stabilizer(SymmetricGroup(BitSize(rho)),rho,PermuteQubitsRho);

SymmetrizePsi:=function(psi)
  local result,perm;
  result:=0;
  for perm in SymmetricGroup(BitSize(psi)) do
    result:=result+PermuteQubitsPsi(psi,perm);
  od;
  return result;
end;

AntiSymmetrizePsi:=function(psi)
  local result,perm;
  result:=0;
  for perm in SymmetricGroup(BitSize(psi)) do
    result:=result+SignPerm(perm)*PermuteQubitsPsi(psi,perm);
  od;
  return result;
end;

unnormedSymmetrizeRho:=function(rho)
  local result,perm;
  result:=0;
  for perm in SymmetricGroup(BitSize(rho)) do
    result:=result+PermuteQubitsRho(rho,perm);
  od;
  return result;
end;

SymmetrizeRho:=function(rho)
  local result,perm,Sn;
  Sn:=SymmetricGroup(BitSize(rho));
  result:=0;
  for perm in Sn do
    result:=result+PermuteQubitsRho(rho,perm);
  od;
  return result/Order(Sn);
#  return normalizeDM(result);
end;

AntiSymmetrizeRho:=function(rho)
  local result,perm;
  result:=0;
  for perm in SymmetricGroup(BitSize(rho)) do
    result:=result+SignPerm(perm)*PermuteQubitsRho(rho,perm);
  od;
  return result;
end;

# not normalized
CycleSymmetrizeRho:=function(rho)
  local result,perm;
  result:=0;
  for perm in CyclicGroup(IsPermGroup,BitSize(rho)) do
    result:=result+PermuteQubitsRho(rho,perm);
  od;
  return result;
end;

# returns the first partition it finds, or false
QuickProductPsi:=function(psi)
  local n,rdm,part;
  n:=BitSize(psi);
  for part in PartitionsSet([1..n],2) do
    rdm:=PartialTrace(DMn(psi),QubitSet2BinarySubsystem(n,part[1]));
    if TraceMat(rdm*rdm)=1 then
      return part;
    fi;
  od;
  return false;
end;

# gap> QuickProductPsi(e000);
# [ [ 1 ], [ 2, 3 ] ]

ProductPsi:=function(psi)
  local helper;
  helper:=function(rho,bitLabelList)
    local n,rdm1,rdm2,part,relabel;
    n:=BitSize(rho);
    if n <> Length(bitLabelList) then
      Error("ProductPsi:  length mismatch; n,bitLabelList = ",n,bitLabelList);
    fi;
    if n=1 then
      return [[bitLabelList[1]]];
    fi;
    relabel:=function(index)
      return bitLabelList[index];
    end;
    for part in PartitionsSet([1..n],2) do
      # tracing OVER part[2]
      rdm1:=PartialTrace(rho,QubitSet2BinarySubsystem(n,part[2]));
      if TraceMat(rdm1*rdm1)=1 then
        rdm2:=PartialTrace(rho,QubitSet2BinarySubsystem(n,part[1]));
        # Print(rdm1,List(part[1],relabel));
        # Print(rdm2,List(part[2],relabel));
        return Concatenation(helper(rdm1,List(part[1],relabel)),
	                     helper(rdm2,List(part[2],relabel)));
      fi;
    od;
    return [List([1..n],relabel)];
  end;
  return helper(DMn(psi),[1..BitSize(psi)]);
end;

# pauli X on qubit l
fastX := function(l,cvec)
    local newvec,item,n,bl;
    n:=BitSize(cvec);
    newvec:=[];
    for item in [1..Length(cvec)] do
        bl:=dec2binlist(item-1,n);
        Append(newvec,[cvec[Bin(bitflip(l,bl))+1]]);
    od;
    return newvec;
end;

# pauli Y on qubit l
fastY := function(l,cvec)
    local newvec,item,n,bl;
    n:=BitSize(cvec);
    newvec:=[];
    for item in [1..Length(cvec)] do
        bl:=dec2binlist(item-1,n);
        if bl[l] = 0 then
            Append(newvec,[-i * cvec[Bin(bitflip(l,bl))+1]]);
        else
            Append(newvec,[ i * cvec[Bin(bitflip(l,bl))+1]]);
        fi;
    od;
    return newvec;
end;

# pauli Z on qubit l
fastZ := function(l,cvec)
  local newvec,item,n,bl;
  n:=BitSize(cvec);
  newvec:=[];
  for item in [1..Length(cvec)] do
    bl:=dec2binlist(item-1,n);
    if bl[l]=1 then
      Append(newvec,[-cvec[item]]);
    else
      Append(newvec,[ cvec[item]]);
    fi;
  od;
  return newvec;
end;

controlledZ := [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]];

cnot := [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]];

# controlled Z on qubits l and m of an n-qubit state
ControlledZ := function(n,l,m)
  return OrderedKron(QubitSet2BinarySubsystem(n,[l,m]),IdentityMat(2^(n-2)),controlledZ);
end;

# GenControlledZ is defined after HypergraphState below

Controlled := function(u)
    return Kron([[[1,0],[0,0]],IdentityMat(Length(u))])
           + Kron([[[0,0],[0,1]],u]);
end;

fastControlledZ := function(l,m,cvec)
  local newvec,item,n,bl;
  n:=BitSize(cvec);
  newvec:=[];
  for item in [1..Length(cvec)] do
    bl:=dec2binlist(item-1,n);
    if bl[l]=1 and bl[m]=1 then
      Append(newvec,[-cvec[item]]);
    else
      Append(newvec,[ cvec[item]]);
    fi;
  od;
  return newvec;
end;

fastControlledZs := function(pairList,cvec)
  if pairList = [] then
    return cvec;
  else
    return fastControlledZs(cdr(pairList),fastControlledZ(pairList[1][1],pairList[1][2],cvec));
  fi;
end;

fastControlledNot := function(l,m,cvec)
    local newvec,item,n,bl;
    n:=BitSize(cvec);
    newvec:=[];
    for item in [1..Length(cvec)] do
        bl:=dec2binlist(item-1,n);
        if bl[l]=1 then
            Append(newvec,[cvec[Bin(bitflip(m,bl))+1]]);
        else
            Append(newvec,[cvec[item]]);
        fi;
    od;
    return newvec;
end;

fastControlledNots := function(pairList,cvec)
  if pairList = [] then
    return cvec;
  else
    return fastControlledNots(cdr(pairList),fastControlledNot(pairList[1][1],pairList[1][2],cvec));
  fi;
end;

# Generalized controlled-Z (for hypergraph states)
fastGenControlledZ := function(qubitset,cvec)
    local newvec,item,n,bl;
    n:=BitSize(cvec);
    newvec:=[];
    for item in [1..Length(cvec)] do
        bl:=dec2binlist(item-1,n);
        if ForAll(qubitset,qb -> bl[qb]=1) then
            Append(newvec,[-cvec[item]]);
        else
            Append(newvec,[ cvec[item]]);
        fi;
    od;
    return newvec;
end;

fastGenControlledZs := function(qubitsetList,cvec)
    if qubitsetList = [] then
        return cvec;
    else
        return fastGenControlledZs(cdr(qubitsetList),fastGenControlledZ(qubitsetList[1],cvec));
    fi;
end;

Operator := function(n,op,qubitorder)
    local m;
    m := BitSize(op);
    return PermuteQubitsRho(Kron(cons(op,replicate(n-m)(I))),
                   MappingPermListList([1..m],qubitorder));
end;

# Figure 2.2 from Mermin's Quantum Computer Science
# gap> Operator(2,cnot,[1,2]) = Operator(2,uhad,[2]) * Operator(2,controlledZ,[1,2]) * Operator(2,uhad,[2]);
# true
# gap> Operator(2,controlledZ,[1,2]) = Operator(2,controlledZ,[2,1]);
# true

# Figure 4.8 of Nielsen and Chuang
# gap> GenControlledZ(3,[1,2,3]) = Operator(3,Controlled(uS),[2,3]) * Operator(3,cnot,[1,2]) * Operator(3,ConjugateTranspose(Controlled(uS)),[2,3]) * Operator(3,cnot,[1,2]) * Operator(3,Controlled(uS),[1,3]);
# true

# lonely column first
fastKPsi := function(psi)
    local n,mat,j;
    n := BitSize(psi);
    mat := [realvec(psi)];
    for j in [1..n] do
        Add(mat,realvec(fastX(j,psi)));
        Add(mat,realvec(fastY(j,psi)));
        Add(mat,realvec(fastZ(j,psi)));
    od;
    return NullspaceMat(mat);
end;

DMFromPolyCoeffs:=function(n,coefflist)
  local maketerm;
  maketerm:=function(term)
    local zeronum;
    zeronum := n - term[2] - term[3] - term[4];
    return term[1] * SymmetrizeRho(Kron(Concatenation([replicate(zeronum)(I)
                                                      ,replicate(term[2])(sx)
                                                      ,replicate(term[3])(sy)
                                                      ,replicate(term[4])(sz)])));
  end;
  return Sum(List(coefflist,maketerm))/2^n;
end;

# can handle constants
ePCOP:=function(poly,var)
  if poly in Rationals then
    return [poly];
  else
    return PolynomialCoefficientsOfPolynomial(poly,var);
  fi;
end;

# can handle constants
eCOUP:=function(poly)
  if poly in Rationals then
    return [poly];
  else
    return CoefficientsOfUnivariatePolynomial(poly);
  fi;
end;

SymRhoFromPoly:=function(n,poly)
  local xlist,xylist,xyzlist,result,coeff,xdeg,ydeg,zdeg;
  xlist:=ePCOP(poly,x);
  xylist:=List(xlist,term -> ePCOP(term,y));
  xyzlist:=List(xylist,termlist -> List(termlist,term -> eCOUP(term)));
  result:=DM(replicate(n)(0));
#  result:=Kron(replicate(n)([[0,0],[0,0]]));
  for xdeg in [0..Length(xyzlist)-1] do
    for ydeg in [0..Length(xyzlist[xdeg+1])-1] do
      for zdeg in [0..Length(xyzlist[xdeg+1][ydeg+1])-1] do
        coeff:=xyzlist[xdeg+1][ydeg+1][zdeg+1];
        result:=result+coeff*SymmetrizeRho(Kron(Concatenation(
          [replicate(n-xdeg-ydeg-zdeg)(I)
          ,replicate(xdeg)(sx)
          ,replicate(ydeg)(sy)
          ,replicate(zdeg)(sz)])));
      od;
    od;
  od;
  return result/2^n;
end;

PolyFromSymRho:=function(rho)
  local poly,n,pil,coefflist;
  if PermStabRho(rho) = SymmetricGroup(BitSize(rho)) then
    n:=BitSize(rho);
    pil:=LieAlgList(n,"u");
    poly:=0;
    for coefflist in pil do
      poly:=poly+PauliCoeff(coefflist,rho)*x^(Count(1,coefflist))
                                          *y^(Count(2,coefflist))
                                          *z^(Count(3,coefflist));
    od;
    return poly;
  else
    return "PolyFromSymRho:  rho is not symmetric";
  fi;
end;

### State makers

GraphState := function(n,edgeList)
  local cz;
  cz := function(edge)
    return ControlledZ(n,edge[1],edge[2]);
  end;
  return Product(List(edgeList,cz)) * ListWithIdenticalEntries(2^n,1);
end;

SymmetricProduct := lst -> Sum(List(PermutationsList(lst),KronVec));

# n = # of qubits, k = # of ones
Dicke:=function(n,k)
  return Sum(List(PermutationsList(
                    Concatenation(ListWithIdenticalEntries(n-k,0),
			          ListWithIdenticalEntries(k,1))
		  ),ee));
end;

UniformDickeMixture :=
  n -> Sum(List([0..n],compose(DMn)(curry2(Dicke)(n))))/(n+1);

Mickey:=function(n,k)
  return normalizeDM(Sum(List(List(PermutationsList(
		    Concatenation(ListWithIdenticalEntries(n-k,0),
			          ListWithIdenticalEntries(k,1))
		  ),ee),DM)));
end;

TotallyMixed:=function(n)
  return normalizeDM(IdentityMat(2^n));
end;

# One-qubit state vectors

e0:=ee([0]);
e1:=ee([1]);
xp:=e0+e1;
xm:=e0-e1;
yp:=e0+i*e1;
ym:=e0-i*e1;
nice:=[5/Sqrt(26),(3+4*i)/5 * 1/Sqrt(26)];
xp120:=[cos(60),sin(60)];
xm120:=[cos(60),-sin(60)];
d3:=SymmetricProduct([e0,xp120,xm120]);

# One-qubit density matrices

dme0:=DM(normalize(e0));
dme1:=DM(normalize(e1));
dmxp:=DM(normalize(xp));
dmxm:=DM(normalize(xm));
dmyp:=DM(normalize(yp));
dmym:=DM(normalize(ym));
dmzp:=DM(normalize(e0));
dmzm:=DM(normalize(e1));
dm1mixed:=normalizeDM(dmzp+dmzm);
dm1mid:=normalizeDM(3*dmzp+dmzm);
dm1nice:=DM(nice);   # nice is already normalized

# Two-qubit state vectors

e00:=e(2,Bin([0,0]));
e01:=e(2,Bin([0,1]));
e10:=e(2,Bin([1,0]));
e11:=e(2,Bin([1,1]));
part:=4/5*e00 + 3/5*e11;
twocat:=e00+e11;
singlet:=e01-e10;
tilted_singlet:=Kron([[[1,0],[0,E(16)]],[[cos(28),sin(28)],[-sin(28),cos(28)]]])*singlet;
tworandom:=RandomQubits(2);

# Two-qubit density matrices

dme00:=DM(e00);
dme11:=DM(e11);
dmpart:=DM(normalize(part));
dm2cat:=DM(normalize(twocat));
dmsinglet:=DM(normalize(singlet));
dm2mixed:=Kron([I/2,I/2]);
dm2nk:=normalizeDM(Kron([dmzp,dmxp])+Kron([dmxp,dmzp]));   # has trivial stabilizer subalgebra
dm2rank2asym:=3/4*dme00+1/4*dme11;
dm2random:=DMn(tworandom);
dm2discrete:=dm2mixed+1/64*Kron([sy,sx])+1/32*Kron([sy,sz])+1/64*Kron([sz,sx]);
  # discrete stabilizer but no continuous stabilizer
rr:=Kron([sx,sx])+Kron([sy,sy])+Kron([sz,sz]);
  # not a state, but Hermitian

# Three-qubit state vectors

e000:=e(3,Bin([0,0,0]));
e001:=e(3,Bin([0,0,1]));
e010:=e(3,Bin([0,1,0]));
e011:=e(3,Bin([0,1,1]));
e100:=e(3,Bin([1,0,0]));
e101:=e(3,Bin([1,0,1]));
e110:=e(3,Bin([1,1,0]));
e111:=e(3,Bin([1,1,1]));

fullab:=(e000+e110)/Sqrt(2);
fullac:=(e000+e101)/Sqrt(2);
fullbc:=(e000+e011)/Sqrt(2);
partab:=4/5*e000+3/5*e110;
partac:=4/5*e000+3/5*e101;
partbc:=4/5*e000+3/5*e011;
ghz:=e000+e111;
partghz:=4/5*e000+3/5*e111;
ghzt1:=e000+e100+e111;
ghzt2:=e000+e101+e111;
ghzt3:=e000+e110+e111;
ghzt12:=e000+e100+e101+e111;
ghzt13:=e000+e100+e110+e111;
ghzt23:=e000+e100+e101+e110+e111;
w:=(e100+e010+e001)/Sqrt(3);
# w is three choose two (same as three choose one)
genw1:=e100/Sqrt(2)+e010/2+e001/2;
genw2:=e100/Sqrt(2)+e010/Sqrt(3)+e001/Sqrt(6);
type4a:=e000/2+e100/2+e101/2+e110/2;
generic:=(2*e000+e001+e010+e011+e100+e101+e110+e111)/Sqrt(11);
#dan:=[-3*i,3+i,-3-3*i,4-2*i,-3*i,3+i,9,-3+9*i];
ghzx:=e000+e011+e101+e110;
ghzy:=e000-e001-e010-e011-e100-e101-e110+e111;
# ghzx = e^i(pi/4)(sigmay,sigmay,sigmay) ghz (named x because of K_rho)
# ghzy = e^i(pi/4)(sigmax,sigmax,sigmax) ghz (named y because of K_rho)
diosi1:=Sqrt(1/8)*e100+Sqrt(1/8)*e010+Sqrt(3/8)*e001+Sqrt(3/8)*e111;
# a state with two totally mixed 1-qubit RDMs and one not totally mixed
diosi2:=Sqrt(1/8)*e100+Sqrt(1/8)*e010+Sqrt(3/8)*e001-Sqrt(3/8)*e111;
hss3:=function(lambda3,phase1,phase2,phase3)
  return Sqrt(lambda3/2)*e100 + phase1*Sqrt(lambda3/2)*e010
    + phase2*Sqrt((1-lambda3)/2)*e001 + phase3*Sqrt((1-lambda3)/2)*e111;
end;
# lambda3 is lower of two eigenvalues
# length of Bloch vector is max-min = 1 - 2*lambda

# Three-qubit density matrices

dme000:=DM(normalize(e000));
dmghz:=DM(normalize(ghz));
dmfullab:=DM(normalize(fullab));
dmw:=DM(normalize(w));
dm3mixed:=Kron([I/2,I/2,I/2]);
#dm3dan:=DM(normalize(dan));

# Four-qubit state vectors

e0000:=ee([0,0,0,0]);
e0001:=ee([0,0,0,1]);
e0010:=ee([0,0,1,0]);
e0011:=ee([0,0,1,1]);
e0100:=ee([0,1,0,0]);
e0101:=ee([0,1,0,1]);
e0110:=ee([0,1,1,0]);
e0111:=ee([0,1,1,1]);
e1000:=ee([1,0,0,0]);
e1001:=ee([1,0,0,1]);
e1010:=ee([1,0,1,0]);
e1011:=ee([1,0,1,1]);
e1100:=ee([1,1,0,0]);
e1101:=ee([1,1,0,1]);
e1110:=ee([1,1,1,0]);
e1111:=ee([1,1,1,1]);

fourcat:=e0000+e1111;
fourchooseone:=e1000+e0100+e0010+e0001;
fourchoosetwo:=e1100+e1010+e1001+e0110+e0101+e0011;
# four choose two was called "toto"
singpair:=normalize(kron(singlet,singlet));
sing12sing34:=singpair;
sing13sing24:=SwitchBitsCvec(singpair,2,3);
sing14sing23:=PermuteQubitsPsi(sing12sing34,(2,4,3));
sing23sing41:=PermuteQubitsPsi(sing12sing34,(1,2,3,4));

# gap> PermuteQubitsPsi(sing23sing41,(1,2,3,4)) = sing12sing34;
# true

# David's crossing identity:
# gap> sing12sing34 + sing14sing23 = sing13sing24;
# true

# old definition:
# sing14sing23:=SwitchBitsCvec(singpair,2,4);
# gap> sing12sing34-sing13sing24=sing14sing23;
# true

quad:=sing12sing34+sing13sing24;
quad2:=sing12sing34+sing14sing23;
quadcheck:=e0011+e0101-2*e0110-2*e1001+e1010+e1100;
twotwo:=e0000-e0011-e1100-e1111;    # is this a product?
m4   :=e0011+e1100+E(3)  *(e0101+e1010)+E(3)^2*(e1001+e0110);
m4bar:=e0011+e1100+E(3)^2*(e0101+e1010)+E(3)  *(e1001+e0110);
psi4:=e0000+kron(normalize(xp),e011)+e1101+kron(normalize(xm),e110);
hss4:=function(lambda3,lambda4,phase1,phase2,phase3,phase4)
  local phi,psi;
  phi:=Sqrt(lambda3-lambda4)*e10
                   + phase1 * Sqrt(1-lambda3-lambda4)*e01;
  psi:=Sqrt(lambda3+lambda4)*e00
                   + phase2 * Sqrt(1-lambda3+lambda4)*e11;
  return kron(e0,kron(phi,e1))
        + phase3 * Sqrt(lambda4/(1/2+lambda4)) * kron(e0,kron(psi,e0))
        + phase4 * Sqrt((1/2)  /(1/2+lambda4)) * kron(e1,kron(psi,e1));
end;
twoof4totallymixed:=(Sqrt(3)*e0000 + Sqrt(25)*e0011
  + Sqrt(5)*e0101 + Sqrt(7)*e0110
  + Sqrt(12)*e1001 + Sqrt(28)*e1111)/Sqrt(80);
TotallyMixedNoStabilizer:=e1001+e0110+e1010+e0101+e0011+e1100+e0000-e1111;
TMNS2:=kron(e0,e000+e011+e101+e110)-kron(e1,e001+e010+e100-e111);
stab1d:=e0000+e0010+e0001+e1111;
roch4:=kron(e0,ghzt3);

# gap> m4 + E(6)*m4bar = -E(3)*i*2*Sqrt(3)*sing13sing24;
# true
# gap> m4 + E(2)*m4bar = i*2*Sqrt(3)*sing12sing34;
# true
# gap> m4 + E(6)^5*m4bar = E(3)^2*i*2*Sqrt(3)*sing14sing23;
# true
# gap> m4 + E(6)^5*m4bar = -E(3)^2*i*2*Sqrt(3)*sing23sing41;
# true

order3:=SymmetricProduct([e0,xp120,xm120,yp]);
g3:=[[-1/2,Sqrt(3)/2],[-Sqrt(3)/2,-1/2]];
order3g:=E(3)^(-1)*Kron([g3,g3,g3,g3]);

order3b:=SymmetricProduct([e0,xp,[1,E(3)],[1,E(3)^2]]);
g3b:=[[E(3),0],[0,E(3)^2]];
order3bg:=E(3)^(-1)*Kron([g3b,g3b,g3b,g3b]);

# 2x2 cluster state ("the square")
cluster22:=fastControlledZs([[1,2],[2,3],[3,4],[4,1]],ListWithIdenticalEntries(2^4,1));

tetrahedron:=SymmetricProduct([e0,Sqrt(1/3)*e0+Sqrt(2/3)*e1,
                                  Sqrt(1/3)*e0+E(3)*Sqrt(2/3)*e1,
				  Sqrt(1/3)*e0+E(3)^2*Sqrt(2/3)*e1]);

superDicke01:=3*Dicke(4,0)+4*Dicke(4,1);

# Four-qubit density matrices

dm4mixed:=Kron([I/2,I/2,I/2,I/2]);
dm4w0:=DMn(e0000);
dm4w1:=DMn(e1000+e0100+e0010+e0001);
dm4w2:=DMn(e1100+e1010+e1001+e0110+e0101+e0011);
dm4w3:=DMn(e1110+e1101+e1011+e0111);
dm4w4:=DMn(e1111);

dm4symmix:=dm4w0/5+dm4w1/5+dm4w2/5+dm4w3/5+dm4w4/5;

# Five-qubit state vectors

fivecat:=ee([0,0,0,0,0])+ee([1,1,1,1,1]);
w5:=ee([1,0,0,0,0])+ee([0,1,0,0,0])+ee([0,0,1,0,0])+ee([0,0,0,1,0])+ee([0,0,0,0,1]);
# w5 is five choose one
fivechooseone:=w5;
fivechoosetwo:=ee([1,1,0,0,0])+ee([1,0,1,0,0])+ee([1,0,0,1,0])+ee([1,0,0,0,1])+ee([0,1,1,0,0])
              +ee([0,1,0,1,0])+ee([0,1,0,0,1])+ee([0,0,1,1,0])+ee([0,0,1,0,1])+ee([0,0,0,1,1]);
roch5:=ee([0,0,0,0,0])+ee([1,0,1,0,0])+ee([0,1,0,0,1])+ee([1,1,1,0,1])+ee([0,1,0,1,1])+ee([1,1,1,1,1]);
pentagon:=GraphState(5,[[1,2],[2,3],[3,4],[4,5],[5,1]]);
fivesu2:=kron(m4,e0)+kron(m4bar,e1);

# Six-qubit state vectors

sixchooseone:=ee([1,0,0,0,0,0])+ee([0,1,0,0,0,0])+ee([0,0,1,0,0,0])
             +ee([0,0,0,1,0,0])+ee([0,0,0,0,1,0])+ee([0,0,0,0,0,1]);
sixchoosetwo:=ee([1,1,0,0,0,0])+ee([1,0,1,0,0,0])+ee([1,0,0,1,0,0])+ee([1,0,0,0,1,0])+ee([1,0,0,0,0,1])
             +ee([0,1,1,0,0,0])+ee([0,1,0,1,0,0])+ee([0,1,0,0,1,0])+ee([0,1,0,0,0,1])+ee([0,0,1,1,0,0])
             +ee([0,0,1,0,1,0])+ee([0,0,1,0,0,1])+ee([0,0,0,1,1,0])+ee([0,0,0,1,0,1])+ee([0,0,0,0,1,1]);
sixchoosethree:=ee([1,1,1,0,0,0])+ee([1,1,0,1,0,0])+ee([1,1,0,0,1,0])+ee([1,1,0,0,0,1])+ee([1,0,1,1,0,0])
               +ee([1,0,1,0,1,0])+ee([1,0,1,0,0,1])+ee([1,0,0,1,1,0])+ee([1,0,0,1,0,1])+ee([1,0,0,0,1,1])
               +ee([0,1,1,1,0,0])+ee([0,1,1,0,1,0])+ee([0,1,1,0,0,1])+ee([0,1,0,1,1,0])+ee([0,1,0,1,0,1])
               +ee([0,1,0,0,1,1])+ee([0,0,1,1,1,0])+ee([0,0,1,1,0,1])+ee([0,0,1,0,1,1])+ee([0,0,0,1,1,1]);
sixrandom:=     ee([1,1,1,0,0,0])+ee([1,1,0,1,0,0])+ee([1,1,0,0,1,0])+ee([1,1,0,0,0,1])+ee([1,0,1,1,0,0])
               +ee([1,0,1,0,1,0])+ee([1,0,1,0,0,1])+ee([1,0,0,1,1,0])+ee([1,0,0,1,0,1])+ee([1,0,0,0,1,1])
               +ee([0,1,1,1,0,0])+ee([0,1,1,0,1,0])+ee([0,1,1,0,0,1])+ee([0,1,0,1,1,0])+ee([0,1,0,1,0,1])
               +ee([0,1,0,0,1,1])+ee([0,0,1,1,1,0])+ee([0,0,1,1,0,1])+ee([0,0,1,0,1,1]);
singthree:=kron(singlet,kron(singlet,singlet));
# same as singthree1:=kron(kron(singlet,singlet),singlet);
# same as singthree2:=kron(singlet,kron(singlet,singlet));

sing12sing34sing56:=singthree;
sing12sing35sing46:=SwitchBitsCvec(singthree,4,5);
sing13sing25sing46:=SwitchBitsCvec(sing12sing35sing46,2,3);

ghz123ghz456:=kron(ghz,ghz);
ghz124ghz356:=SwitchBitsCvec(ghz123ghz456,3,4);

ghz123ghz546:=SwitchBitsCvec(ghz123ghz456,4,5);  # this doesn't do anything
ghz125ghz436:=SwitchBitsCvec(ghz123ghz456,3,5);
ghz125ghz346:=SwitchBitsCvec(ghz125ghz436,4,5);  # this does something
ghz126ghz345:=SwitchBitsCvec(ghz125ghz346,3,6);

# a better product notation

ghzghzaaabbb:=kron(ghz,ghz);
ghzghzaababb:=SwitchBitsCvec(ghzghzaaabbb,3,4);
ghzghzaabbab:=SwitchBitsCvec(ghzghzaaabbb,3,5);
ghzghzaabbba:=SwitchBitsCvec(ghzghzaaabbb,3,6);

fourchunktwo:=kron(quad,e00)+kron(quad2,e11);

octahedron:=SymmetricProduct([e0,e1,xp,xm,yp,ym]);
rot120111:=cos(120)*I+i*sin(120)*(sx+sy+sz)/Sqrt(3);
rot1206:=Kron([rot120111,rot120111,rot120111,rot120111,rot120111,rot120111]);
# gap> rot1206*octahedron=octahedron;
# true                                    ; this has an order 3 element
# gap> rot1206*rot1206*rot1206=IdentityMat(64);
# true
# gap> rot1206=IdentityMat(64);
# false

# Seven-qubit state vectors

# seven-qubit Steane code from Nielsen-Chuang, p. 453
steane0 := ee([0,0,0,0,0,0,0]) + ee([1,0,1,0,1,0,1]) + ee([0,1,1,0,0,1,1]) + ee([1,1,0,0,1,1,0])
         + ee([0,0,0,1,1,1,1]) + ee([1,0,1,1,0,1,0]) + ee([0,1,1,1,1,0,0]) + ee([1,1,0,1,0,0,1]);
steane1 := Kron([sx,sx,sx,sx,sx,sx,sx]) * steane0;
# steane0 is a stabilizer state! (it has a 128-element Pauli group)

# Nine-qubit state vectors

# 3x3 cluster state
cluster33:=fastControlledZs([[1,2],[2,3],[4,5],[5,6],[7,8],[8,9],
                             [1,4],[4,7],[2,5],[5,8],[3,6],[6,9]],ListWithIdenticalEntries(2^9,1));

# Some density matrices that rely on the partial trace

dm2fromghz:=PartialTrace(DM(normalize(ghz)),[1,0,0]);   # two qubit dm

# n-qubit state vectors

ngon:=function(n)
  local makeState;
  makeState:=function(k)
    return e0+E(n)^k*e1;
  end;
  return SymmetricProduct(List([0..(n-1)],makeState));
end;

shiftedngon:=function(n,shift)
  local makeState;
  makeState:=function(k)
    return e0+E(shift)*E(n)^k*e1;
  end;
  return SymmetricProduct(List([0..(n-1)],makeState));
end;

# gap> ngon(4);
# [ 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -24 ]
# gap> shiftedngon(4,8);
# [ 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24 ]

# General functions

T:=function(n,bit,abc)
  # n = number of qubits
  # bit = 1..n
  # abc = 1 for A, 2 for B, 3 for C
  local b,mat;
  if bit=1 then
    mat:=XX[abc];
  else
    mat:=I;
  fi;
  for b in [2..n] do
    if b=bit then
      mat:=KroneckerProduct(mat,XX[abc]);
    else
      mat:=KroneckerProduct(mat,I);
    fi;
  od;
  return mat;
end;

triprank:=function(cvec,list)
  local m,n,qubit,rankwo,rankw;
  n:=LogInt(Length(cvec),2);
  m:=[];
  for qubit in list do
    Append(m,[realvec(T(n,qubit,1)*cvec),
              realvec(T(n,qubit,2)*cvec),
              realvec(T(n,qubit,3)*cvec)]);
  od;
  if m=[] then
    rankwo:=0;
  else
    rankwo:=RankMat(m);
  fi;
  Add(m,realvec(-i*cvec));
  rankw:=RankMat(m);
  return [rankwo,rankw];
end;

M:=function(cvec)
  local m,n,qubit;
  n:=LogInt(Length(cvec),2);
  m:=[];
  for qubit in [1..n] do
    Append(m,[realvec(T(n,qubit,1)*cvec),
              realvec(T(n,qubit,2)*cvec),
              realvec(T(n,qubit,3)*cvec)]);
  od;
  Add(m,realvec(-i*cvec));
  return m;
end;

SmallM:=function(cvec)
  return M(cvec)*TransposedMat(M(cvec));
end;

M3n:=function(cvec)
  local m,n,qubit,s,t1,t2,t3;
  n:=LogInt(Length(cvec),2);
  m:=[];
  s:=realvec(-i*cvec);
  for qubit in [1..n] do
    t1:=realvec(T(n,qubit,1)*cvec);
    t2:=realvec(T(n,qubit,2)*cvec);
    t3:=realvec(T(n,qubit,3)*cvec);
    Append(m,[t1-(s*t1)*s,
              t2-(s*t2)*s,
              t3-(s*t3)*s]);
  od;
  return m;
end;

SmallM3n:=function(cvec)
  return M3n(cvec)*TransposedMat(M3n(cvec));
end;

triprank3n:=function(cvec,list)
  local m,n,qubit,rank,s,t1,t2,t3;
  n:=LogInt(Length(cvec),2);
  m:=[];
  s:=realvec(-i*cvec);
  for qubit in list do
    t1:=realvec(T(n,qubit,1)*cvec);
    t2:=realvec(T(n,qubit,2)*cvec);
    t3:=realvec(T(n,qubit,3)*cvec);
    Append(m,[t1-(s*t1)*s,
              t2-(s*t2)*s,
              t3-(s*t3)*s]);
  od;
  if m=[] then
    rank:=0;
  else
    rank:=RankMat(m);
  fi;
  return rank;
end;

OrbitDim:=function(cvec)
  local n;
  n:=LogInt(Length(cvec),2);
  return triprank(cvec,[1..n])[2]-1;
end;

LUOrbitDimDM:=function(dm)
  local n,flatbrackmat,coefflist,bit,l;
  n:=LogInt(Length(dm),2);
  flatbrackmat:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for l in [1..3*n] do
    bit:=Int((l-1)/3)+1;
    coefflist[bit]:=((l-1) mod 3) + 1;
    Append(flatbrackmat,[flatten(bracket(Kron(List(coefflist,Pauli)),dm))]);
    coefflist[bit]:=0;
  od;
  return RankMat(flatbrackmat);
end;

LUNullDM:=function(dm)
  local n,flatbrackmat,coefflist,bit,l;
  n:=LogInt(Length(dm),2);
  flatbrackmat:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for l in [1..3*n] do
    bit:=Int((l-1)/3)+1;
    coefflist[bit]:=((l-1) mod 3) + 1;
    Append(flatbrackmat,[flatten(bracket(Kron(List(coefflist,Pauli)),dm))]);
    coefflist[bit]:=0;
  od;
  return NullspaceMat(flatbrackmat);
end;

KRho:=LUNullDM;

LUNullityDM:=function(dm)
  return RankMat(LUNullDM(dm));
end;

Kcomp:=function(dm)
  Display("K_rho");
  Display(LUNullDM(dm));
  Display("K_Scooped(rho)");
  Display(LUNullDM(DMFromPF(PartialTracePF(PauliForm(dm)))));
end;

SUOrbitDimDM:=function(dm)
  local n,flatbrackmat,coefflist,bit,l;
  n:=LogInt(Length(dm),2);
  flatbrackmat:=[];
  for l in [1..4^n-1] do
    coefflist:=dec2baseNlist(l,4,n);
    Append(flatbrackmat,[flatten(bracket(Kron(List(coefflist,Pauli)),dm))]);
  od;
  return RankMat(flatbrackmat);
end;

SUNullDMold:=function(dm)
  local n,flatbrackmat,coefflist,bit,l;
  n:=LogInt(Length(dm),2);
  flatbrackmat:=[];
  for l in [1..4^n-1] do
    coefflist:=dec2baseNlist(l,4,n);
    Append(flatbrackmat,[flatten(bracket(Kron(List(coefflist,Pauli)),dm))]);
  od;
  return NullspaceMat(flatbrackmat);
end;

SUNullDM:=function(dm)
  local n,flatbrackmat,coefflist,bit,l,listlength,list,term,listelt;
  n:=LogInt(Length(dm),2);
  flatbrackmat:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for listlength in [1..n] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
        od;
        Append(flatbrackmat,[flatten(bracket(Kron(List(coefflist,Pauli)),dm))]);
      od;
    od;
  od;
  return NullspaceMat(flatbrackmat);
end;

SUNullityDM:=function(dm)
  return RankMat(SUNullDM(dm));
end;

LinearIrredVec:=function(cvec)
  local n,worker;
  n:=BitSize(cvec);
  worker:=function(coefflist)
    return realvec(Kron(List(coefflist,Pauli))*cvec);
  end;
  return NullspaceMat(TransposedMat(List(LieAlgListF(n,0,n-1),worker)));
end;

kerproj:=function(cvec,list)
  local ker,m,n,qubit,rankwo,rankw;
  n:=LogInt(Length(cvec),2);
  ker:=TransposedMat(NullspaceMat(M(cvec)));
  m:=[];
  if ker=[] then
    rankwo:=0;
    rankw:=0;
  else
    for qubit in list do
      Append(m,ker{[3*(qubit-1)+1..3*(qubit-1)+3]});
    od;
#    Display(m);
    if m=[] then
      rankwo:=0;
    else
      rankwo:=RankMat(m);
    fi;
    Add(m,ker[3*n+1]);
#    Display(m);
    rankw:=RankMat(m);
  fi;
  return [rankwo,rankw];
end;

freport:=function(cvec)
  local n,listlength,list,result,row,rank3n,kerinfo;
  n:=LogInt(Length(cvec),2);
#  result:=[];
  result:=[[""      ,"Rank","Rank","KerP","KerP"],
           [""      ,"w/o" ,"with","w/o" ,"with"],
           ["Qubits","lc"  ,"lc"  ,"lc"  ,"lc"  ]];
  for listlength in [0..n] do
    for list in Combinations([1..n],listlength) do
      row:=triprank(cvec,list);
      rank3n:=triprank3n(cvec,list);
      kerinfo:=kerproj(cvec,list);
      Add(result,[list,row[1],row[2],kerinfo[1],kerinfo[2]]);
#      Add(result,[list,row[1],row[2],row[2]-1,rank3n]);
    od;
  od;
  return result;
end;

freportDM:=function(dm)
  local n,listlength,list,result,row,rank3n,kerinfo,pt,bit,coefflist,listelt;
  n:=LogInt(Length(dm),2);
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  result:=[["Trace","Bits","LU"     ,"SU"     ],
           ["over","left","Nullity","Nullity"]];
  for listlength in [0..n-1] do
    for list in Combinations([1..n],listlength) do
      for bit in [1..n] do
        coefflist[bit]:=0;
      od;
      for listelt in [1..listlength] do
        coefflist[list[listelt]]:=1;
      od;
      pt:=PartialTrace(dm,coefflist);
      Add(result,[list,Difference([1..n],list),LUNullityDM(pt),SUNullityDM(pt)]);
    od;
  od;
  return result;
end;

info:=function(cvec)
  display(freport(cvec));
  return;
end;

infoDM:=function(dm)
  display(freportDM(dm));
  return;
end;

nullinfo:=function(cvec)
  Display(NullspaceMat(M(cvec)));
  return;
end;

KR:=function(dm)
# kernel report
  local n,brackmat,coefflist,bit,listlength,bigmat,list,listelt,l,term,bigmatrow;
  n:=LogInt(Length(dm),2);
  bigmat:=[];
  brackmat:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  for l in [1..3*n] do
    bit:=Int((l-1)/3)+1;
    coefflist[bit]:=((l-1) mod 3) + 1;
    Append(brackmat,[bracket(Kron(List(coefflist,Pauli)),dm)]);
    coefflist[bit]:=0;
  od;
  for listlength in [0..n] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
        od;
        bigmatrow:=[];
        for l in [1..3*n] do
          Append(bigmatrow,[PauliCoeff(coefflist,brackmat[l])]);
        od;
        Append(bigmat,[bigmatrow]);
      od;
    od;
    Display(listlength);
    Display(NullspaceMat(TransposedMat(bigmat)));
  od;
end;

KRdim:=function(dm)
# kernel report
  local n,brackmat,coefflist,bit,listlength,bigmat,list,listelt,l,term,bigmatrow,result;
  n:=LogInt(Length(dm),2);
  bigmat:=[];
  brackmat:=[];
  coefflist:=[];
  for bit in [1..n] do
    Append(coefflist,[0]);
  od;
  result:=[["Bits","LU"     ,"SU"     ],
           ["left","Nullity","Nullity"]];
  for l in [1..3*n] do
    bit:=Int((l-1)/3)+1;
    coefflist[bit]:=((l-1) mod 3) + 1;
    Append(brackmat,[bracket(Kron(List(coefflist,Pauli)),dm)]);
    coefflist[bit]:=0;
  od;
  for listlength in [0..n] do
    for list in Combinations([1..n],listlength) do
      for term in [0..3^listlength-1] do
        for bit in [1..n] do
          coefflist[bit]:=0;
        od;
        for listelt in [1..listlength] do
          coefflist[list[listelt]]:=dec2baseNlist(term,3,listlength)[listelt]+1;
        od;
        bigmatrow:=[];
        for l in [1..3*n] do
          Append(bigmatrow,[PauliCoeff(coefflist,brackmat[l])]);
        od;
        Append(bigmat,[bigmatrow]);
      od;
    od;
    Add(result,[listlength,RankMat(NullspaceMat(TransposedMat(bigmat)))]);
  od;
  display(result);
end;

MakeBlock:=function(listb4list,dm)
  local b4list,block;
  block:=[];
  for b4list in listb4list do
    Add(block,PauliListForm(bracket(i*Kron(List(b4list,Pauli)),dm)));
  od;
  return block;
end;

irred:=function(dm)
  local n,result,lun0,sun0,luns0,suns0,lun1,sun1,luns1,suns1,bit,matlu,matsu,dmrp;
  n:=BitSize(dm);

#  luns0:=NullspaceMat(MakeBlock(LieAlgList(n,"lu"),dm));
#  Display("LU Null Space spanned by");
#  Display(luns0);
#  lun0:=RankMat(luns0);
#  Print("LU Nullity is ",lun0,"\n");
#  Display("SU MakeBlock");
#  Display(MakeBlock(LieAlgList(n,"su"),dm));
#  suns0:=NullspaceMat(MakeBlock(LieAlgList(n,"su"),dm));
#  Display("SU Null Space spanned by");
#  Display(suns0);
#  sun0:=RankMat(suns0);
#  Print("SU Nullity is ",sun0,"\n");

  dmrp:=Pad(OneBinList(1,n),PartialTrace(dm,OneBinList(1,n)));
  matlu:=MakeBlock(LieAlgList(n,"lu"),dmrp);
  matsu:=MakeBlock(LieAlgList(n,"su"),dmrp);
  for bit in [2..n] do
    dmrp:=Pad(OneBinList(bit,n),PartialTrace(dm,OneBinList(bit,n)));
    AppendMat(matlu,MakeBlock(LieAlgList(n,"lu"),dmrp));
    AppendMat(matsu,MakeBlock(LieAlgList(n,"su"),dmrp));
  od;
  luns1:=NullspaceMat(matlu);
#  Display("RDM LU Null Space spanned by");
#  Display(luns1);
  lun1:=RankMat(luns1);
  suns1:=NullspaceMat(matsu);
  Display("RDM SU Null Space spanned by");
  Display(suns1);
  sun1:=RankMat(suns1);
  AppendMat(matlu,MakeBlock(LieAlgList(n,"lu"),dm));
  AppendMat(matsu,MakeBlock(LieAlgList(n,"su"),dm));
  luns0:=NullspaceMat(matlu);
  lun0:=RankMat(luns0);
  suns0:=NullspaceMat(matsu);
  Display("RDM+DM SU Null Space spanned by");
  Display(suns0);
  sun0:=RankMat(suns0);
  result:=[["","LU"     ,"SU"     ],
           ["","Nullity","Nullity"],
           ["RDMs",lun1,sun1],
           ["RDMs+DM",lun0,sun0]];
  display(result);
  return [luns1,luns0,suns1,suns0];
end;

displayfile:=function(data,filename)
  local info,row,max,col;
  max:=[7,10,6,5];
  for row in data do
    for col in [1..4] do
      if Length(String(row[col])) > max[col] then
        max[col]:=Length(String(row[col]));
      fi;
    od;
  od;
  AppendTo(filename,String("Triples",max[1]+1),String("Rank w/o s",max[2]+3),
        String("with s",max[3]+3),String("ws -1",max[4]+3),"\n");
  for row in data do
    AppendTo(filename,String(row[1],max[1]+1),String(row[2],max[2]+3),
          String(row[3],max[3]+3),String(row[4],max[4]+3),"\n");
  od;
end;

threeqbportfolio:=function(filename)
  AppendTo(filename,"\n\ne000\n\n");
  displayfile(freport(e000),filename);
  AppendTo(filename,"\n\nfullab:=(e000+e110)/Sqrt(2)\n\n");
  displayfile(freport(fullab),filename);
  AppendTo(filename,"\n\nfullac:=(e000+e101)/Sqrt(2)\n\n");
  displayfile(freport(fullac),filename);
  AppendTo(filename,"\n\nfullbc:=(e000+e011)/Sqrt(2)\n\n");
  displayfile(freport(fullbc),filename);
  AppendTo(filename,"\n\npartab:=4/5*e000+3/5*e110\n\n");
  displayfile(freport(partab),filename);
  AppendTo(filename,"\n\npartac:=4/5*e000+3/5*e101\n\n");
  displayfile(freport(partac),filename);
  AppendTo(filename,"\n\npartbc:=4/5*e000+3/5*e011\n\n");
  displayfile(freport(partbc),filename);
  AppendTo(filename,"\n\nghz:=(e000+e111)/Sqrt(2)\n\n");
  displayfile(freport(ghz),filename);
  AppendTo(filename,"\n\npartghz:=4/5*e000+3/5*e111\n\n");
  displayfile(freport(partghz),filename);
  AppendTo(filename,"\n\nghzt1:=e000+e100+e111\n\n");
  displayfile(freport(ghzt1),filename);
  AppendTo(filename,"\n\nghzt2:=e000+e101+e111\n\n");
  displayfile(freport(ghzt2),filename);
  AppendTo(filename,"\n\nghzt3:=e000+e110+e111\n\n");
  displayfile(freport(ghzt3),filename);
  AppendTo(filename,"\n\nghzt12:=e000+e100+e101+e111\n\n");
  displayfile(freport(ghzt12),filename);
  AppendTo(filename,"\n\nghzt13:=e000+e100+e110+e111\n\n");
  displayfile(freport(ghzt13),filename);
  AppendTo(filename,"\n\nghzt23:=e000+e100+e101+e110+e111\n\n");
  displayfile(freport(ghzt23),filename);
  AppendTo(filename,"\n\nw:=(e100+e010+e001)/Sqrt(3)\n\n");
  displayfile(freport(w),filename);
  AppendTo(filename,"\n\ngenw1:=e100/Sqrt(2)+e010/2+e001/2\n\n");
  displayfile(freport(genw1),filename);
  AppendTo(filename,"\n\ngenw2:=e100/Sqrt(2)+e010/Sqrt(3)+e001/Sqrt(6)\n\n");
  displayfile(freport(genw2),filename);
  AppendTo(filename,"\n\ntype4a:=e000/2+e100/2+e101/2+e110/2\n\n");
  displayfile(freport(type4a),filename);
  AppendTo(filename,"\n\ngeneric:=(2*e000+e001+e010+e011+e100+e101+e110+e111)/Sqrt(11)\n\n");
  displayfile(freport(generic),filename);
end;

# 16 Apr 2008


###
# NChooseK() returns a list of all possible combinations where k elements are chosen from a list of n elements
###
NChooseK:=function(n,k)
  local func;
  func:=function(bitset)
    return QubitSet2BinarySubsystem(n,bitset);
  end;
  return List(Combinations([1..n],k),func);
end;

###
# KPsi() finds a basis for K-Psi
# KPsi() returns a list with 7 entries, which mean
# [ First qubit Pauli Z, First qubit Pauli Y, First qubit Pauli X, Second qubit Pauli Z, Second qubit Pauli Y, Second qubit Pauli X, Phase ]
# For example, the output [ 1, 0, 0, 0, 1, 0, 0 ]
# would indicate Z^(1) + Y^(2)
###
KPsi:=function(cvec)
  return NullspaceMat(M(cvec));
end;

### GAP doesn't like this functional code.  It won't do more than 5 qubits.
# PauliStabPsi:=function(psi)
#   local stabpred;
#   stabpred:=function(pge)
#     local pnt,phase;
#     phase:=car(pge);
#     pnt:=cdr(pge);
#     return phase*Kron(List(pnt,Pauli))*psi = psi;
#   end;
#   return filter(PauliGroup(BitSize(psi)),stabpred);
# end;

# gap> PauliStabPsi(ghz);
# [ [ 1, 0, 0, 0 ], [ 1, 3, 3, 0 ], [ 1, 3, 0, 3 ], [ 1, 0, 3, 3 ],
#   [ 1, 1, 1, 1 ], [ -1, 1, 2, 2 ], [ -1, 2, 1, 2 ], [ -1, 2, 2, 1 ] ]
# first entry is a phase, remaining three are Pauli numbers 0..3

PauliStabPsi:=function(psi)
  local pge,results,phase,pnt;
  results:=[];
  for pge in PauliGroup(BitSize(psi)) do
    phase:=car(pge);
    pnt:=cdr(pge);
    if phase*Kron(List(pnt,Pauli))*psi = psi then
      Append(results,[pge]);
    fi;
  od;
  return results;
end;

PauliStabRho:=function(rho)
  local stabpred;
  stabpred:=function(pnt)
    local pgeMatNoPhase;
    pgeMatNoPhase:=Kron(List(pnt,Pauli));   # this is its own inverse
    return pgeMatNoPhase*rho*pgeMatNoPhase = rho;
  end;
  return filter(PauliGroupNoPhase(BitSize(rho)),stabpred);
end;

# gap> PauliStabRho(DMn(ghz));
# [ [ 0, 0, 0 ], [ 3, 3, 0 ], [ 3, 0, 3 ], [ 0, 3, 3 ],
#   [ 1, 1, 1 ], [ 1, 2, 2 ], [ 2, 1, 2 ], [ 2, 2, 1 ] ]
# no phases here, just Pauli numbers 0..3

SinglePauliProd:=function(a,b)
  local phasea,phaseb,paulia,paulib,extraphase;
  phasea:=a[1];
  paulia:=a[2];
  phaseb:=b[1];
  paulib:=b[2];
  if paulia=0 then
    return [phasea*phaseb,paulib];
  elif paulib=0 then
    return [phasea*phaseb,paulia];
  elif paulia=paulib then
    return [phasea*phaseb,0];
  else
    if (paulib-paulia) mod 3 = 1 then
      extraphase:=i;
    elif (paulib-paulia) mod 3 = 2 then
      extraphase:=-i;
    else
      Display("Error:  SinglePauliProd:  [paulia,paulib] = ");
      Display([paulia,paulib]);
      return 0;
    fi;
    return [phasea*phaseb*extraphase,6-paulia-paulib];
  fi;
end;

PauliProd:=function(a,b)
  local helper;
  helper:=function(result,NonPhaseList1,NonPhaseList2)
    local mult,newphase,newpauli,accumulatedPhase,accumulatedPauli;
    if Length(NonPhaseList1)=0 then
      return result;
    else
      mult:=SinglePauliProd([1,NonPhaseList1[1]],[1,NonPhaseList2[1]]);
      newphase:=mult[1];
      newpauli:=mult[2];
      accumulatedPhase:=newphase*result[1];
      accumulatedPauli:=Concatenation(cdr(result),[newpauli]);
      return helper(cons(accumulatedPhase,accumulatedPauli),
                    cdr(NonPhaseList1),cdr(NonPhaseList2));
    fi;
  end;
  return helper([a[1]*b[1]],cdr(a),cdr(b));
end;

PauliMat:=function(pge)
  local pnt,phase;
  phase:=car(pge);
  pnt:=cdr(pge);
  return phase*Kron(List(pnt,Pauli));
end;

# angle in degrees
SU2RotAxisAngle:=function(axis,angle)
  local unitv,sig;
  if Length(axis) = 3 then
    unitv:=normalize(axis);
    sig:=unitv[1]*sx + unitv[2]*sy + unitv[3]*sz;
    return cos(angle/2) * I - i * sin(angle/2) * sig;
  else
    return "SU2RotAxisAngle: axis must be a vector with 3 elements";
  fi;
end;

# examples
rotY45:=SU2RotAxisAngle([0,1,0],45);
rotX90:=SU2RotAxisAngle([1,0,0],90);

U2ToSO3:=function(op)
  local opd;
  opd:=ConjugateTranspose(op);
  return 1/2*[[Trace(sx*op*sx*opd),Trace(sx*op*sy*opd),Trace(sx*op*sz*opd)]
             ,[Trace(sy*op*sx*opd),Trace(sy*op*sy*opd),Trace(sy*op*sz*opd)]
             ,[Trace(sz*op*sx*opd),Trace(sz*op*sy*opd),Trace(sz*op*sz*opd)]];
end;

SO3Poly:=function(mat,poly)
  return Value(poly,[x,y,z],[x,y,z]*mat);
end;

U2Poly:=function(g,poly)
  return SO3Poly(U2ToSO3(g),poly);
end;

RotZ := deg -> [[cos(deg),-sin(deg),0]
               ,[sin(deg), cos(deg),0]
               ,[       0,        0,1]];

# The 48-element non-abelian subgroup of SU(2) that does rotations
# of the cube.
SU2Cube:=Group(EA4,EB4,EC4);

# gap> Order(SU2Cube);
# 48
# gap> IsAbelian(SU2Cube);
# false

### Table formatting stuff ###

PadRight:=function(str,len)
  local padlength;
  padlength:=len-Length(str);
  if padlength <= 0 then
    return str;
  else
    return PadRight(Concatenation(str," "),len);
  fi;
end;

ListString:=function(lst)
  return List(lst,String);
end;

ListListString:=function(lst)
  return List(lst,ListString);
end;

MaxLength:=function(lst)
  return Maximum(List(lst,Length));
end;

PrepareLine:=function(lst,format)
  local data,makestring;
  data:=zip(lst,format);
  makestring:=function(pair)
    return PadRight(pair[1],pair[2]);
  end;
  return Concatenation(Concatenation(List(data,makestring)),"\n");
end;

PrintLine:=function(lst,format)
  Print(PrepareLine(lst,format));
end;

PrepareTable:=function(tab)
  local colwidths,addwidth,printfunc;
  colwidths:=List(TransposedMat(tab),MaxLength);
  addwidth:=function(n)
    return n+2;
  end;
  printfunc:=function(lst)
    return PrepareLine(lst,List(colwidths,addwidth));
  end;
  return Concatenation(List(tab,printfunc));
end;

# print a table of strings, left justified
PrintTable:=function(tab)
  Print(PrepareTable(tab));
end;

### Pauli multiplication table

PrintPauliMultTable:=function()
  local pl,pt,multlist;
  pl:=[[1,0],[1,1],[1,2],[1,3]];
  multlist:=function(pe)
    local multbypeonLeft;
    multbypeonLeft:=function(pr)
      return SinglePauliProd(pe,pr);
    end;
    return List(pl,multbypeonLeft);
  end;
  pt:=List(pl,multlist);
  PrintTable(ListListString(pt));
end;

### Hilbert-Schmidt inner product

HilbertSchmidtInnerProduct:=function(a,b)
  return Trace(ConjugateTranspose(a)*b);
end;

HilbertSchmidtNorm := rho -> HilbertSchmidtInnerProduct(rho,rho);

### Vector spaces stabilized by stuff

VecSpaceStabByGroup:=function(glist)   # glist is a list of operators
  local n,subtractIdentity;
  n := BitSize(glist[1]);
  subtractIdentity:=function(op)
    return op - IdentityMat(2^n);
  end;
  return NullspaceMat(TransposedMat(Concatenation(List(glist,subtractIdentity))));
end;

OpFromABCLinAlg:=function(abcelem)
  local n,bigTs;
  n := (Length(abcelem) - 1)/3;
  bigTs:=function(bitnum,abclist)
    local len;
#    Display(["bitnum,abclist: ",bitnum,abclist]);
    len := (Length(abclist) - 1)/3;
    if len <= 0 then
      return abclist[1]*(-i)*IdentityMat(2^n);
    else
      return abclist[1]*T(n,bitnum,1) + abclist[2]*T(n,bitnum,2) + abclist[3]*T(n,bitnum,3)
              + bigTs(bitnum+1,cdr(cdr(cdr(abclist))));
    fi;
  end;
  return bigTs(1,abcelem);
end;

VecSpaceStabByAlg:=function(alist)
  return NullspaceMat(TransposedMat(Concatenation(List(alist,OpFromABCLinAlg))));
end;

# gap> VecSpaceStabByGroup(List(PauliStabPsi(ghz),PauliMat));
# [ [ 1, 0, 0, 0, 0, 0, 0, 1 ] ]

# gap> VecSpaceStabByAlg(KPsi(ghz));
# [ [ 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1 ] ]

# gap> Length(VecSpaceStabByAlg(KPsi(sixchoosethree)));
# 20

# We want the vector space of Hermitian operators stabilized by a list of group elements.
# First, we need a way (I think) to convert an operator on Hilbert space,
# which acts on density matrices by conjugation, to an operator on Pauli List Forms,
# which are vectors representing Hermitian matrices.

HilbertOpToPLFOp:=function(op)
  local n,actOnPLF;
  n := BitSize(op);
  actOnPLF:=function(plfvec)
    return PauliListForm(op*DMFromPLF(plfvec)*ConjugateTranspose(op));
  end;
  return TransposedMat(List(IdentityMat(4^n),actOnPLF));
end;

PLFStabByGroup:=function(glist)   # glist is a list of operators (on Hilbert space)
  return VecSpaceStabByGroup(List(glist,HilbertOpToPLFOp));
end;

# for use with density matrices (no lonely column)
OpFromXYZLinAlg:=function(xyzelem)
  local n,bigTs;
  n := Length(xyzelem)/3;
  return 2^n*DMFromPLF(Concatenation([0],xyzelem,replicate(4^n-3*n-1)(0)));
end;

HilbertOpToPLFOpLieAlg:=function(op)
  local n,actOnPLF;
  n := BitSize(op);
  actOnPLF:=function(plfvec)
    local dm;
    dm := DMFromPLF(plfvec);
    return PauliListForm(op*dm-dm*op);
  end;
  return TransposedMat(List(IdentityMat(4^n),actOnPLF));
end;

PLFStabByAlg:=function(alistxyz)
  return NullspaceMat(TransposedMat(Concatenation(
           List(alistxyz,compose(HilbertOpToPLFOpLieAlg)(OpFromXYZLinAlg)))));
end;

# I want to get results like (-2)*Y*X*Z
PauliAlgebra:=FreeAssociativeAlgebra(Rationals,"I","X","Y","Z");
GPA:=GeneratorsOfAlgebra(PauliAlgebra);
S0:=GPA[1];
S1:=GPA[2];
S2:=GPA[3];
S3:=GPA[4];
PauliAlg:=j->GPA[j+1];
PAFFromPLF:=function(plf)
  local n,lst,result,j;
  n:=LogInt(Length(plf),4);
  Assert(0,Length(plf) = 4^n,"Length of plf not a multiple of 4 ");
  lst:=LieAlgList(n,"u");
  result:=plf[1]*Product(List(lst[1],PauliAlg));
  for j in [2..Length(plf)] do
    result:=result + plf[j]*Product(List(lst[j],PauliAlg));
  od;
  return result;
end;
PAFStabByAlg:=function(alistxyz)
  return List(PLFStabByAlg(alistxyz),PAFFromPLF);
end;
RR:=S1*S1+S2*S2+S3*S3;
XYZ:=S1*S2*S3+S2*S3*S1+S3*S1*S2-S1*S3*S2-S2*S1*S3-S3*S2*S1;

PAFFromDM := compose(PAFFromPLF)(PauliListForm);
PAF := PAFFromDM;

sr := (sx + i*sy)/Sqrt(2);
sl := (sx - i*sy)/Sqrt(2);

RZLMat:=function(num)
#    depends on Basic Objects
  if num = 0 then
    return I;
  elif num = 1 then
    return sr;
  elif num = 2 then
    return sz;
  elif num = 3 then
    return sl;
  else
    Error("Pauli:  bad index");
  fi;
end;

RZLAlgebra:=FreeAssociativeAlgebra(Rationals,"I","R","Z","L");
GRZL:=GeneratorsOfAlgebra(RZLAlgebra);
RZLAlg:=j->GRZL[j+1];
RZLAF:=function(rho)
  local n,lst,result;
  n:=BitSize(rho);
  lst:=LieAlgList(n,"u");
  return Sum(fmap(t -> HilbertSchmidtInnerProduct(Kron(fmap(RZLMat)(t)),rho) * Product(fmap(RZLAlg)(t)))(lst));
end;

# n is the number of qubits
PermutationMatrix := n -> perm ->
  TransposedMat(List(IdentityMat(2^n),flip(curry2(PermuteQubitsPsi))(perm)));

PermutationSymmetricHilbertSpace := n ->
  VecSpaceStabByGroup(List(SymmetricGroup(n),PermutationMatrix(n)));

PermutationSymmetricHermitianSpace := n ->
  PLFStabByGroup(List(SymmetricGroup(n),PermutationMatrix(n)));

### Stabilizer code stuff

# Five Qubit Code
FiveQubitCodeGenerators:=[[1,1,3,3,1,0],
                          [1,0,1,3,3,1],
			  [1,1,0,1,3,3],
			  [1,3,1,0,1,3]];   # phase then five qubits

FiveQubitLogicalCodewords:=VecSpaceStabByGroup(List(FiveQubitCodeGenerators,PauliMat));
fiveQubitCodeword0:=FiveQubitLogicalCodewords[1];   # probably not a codeword
fiveQubitCodeword1:=FiveQubitLogicalCodewords[2];   # probably not a codeword
# pentagon lives in codespade of 5-qubit code:
# gap> pentagon+fiveQubitCodeword0+fiveQubitCodeword1;
# [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

# Nielsen-Chuang codewords (page 469)
NCFiveQubitCodeword0 := -fiveQubitCodeword0;
NCFiveQubitCodeword1 :=  fiveQubitCodeword1;

# My favorite codeword convention
WFiveQubitCodeword0 := -fiveQubitCodeword0;
WFiveQubitCodeword1 := -fiveQubitCodeword1;
WFiveQubitCodewordPlus  := WFiveQubitCodeword0 + WFiveQubitCodeword1;
WFiveQubitCodewordMinus := WFiveQubitCodeword0 - WFiveQubitCodeword1;
# gap> HypergraphState(5,[[1,2],[2,3],[3,4],[4,5],[5,1]]) = WFiveQubitCodewordPlus;
# true
# gap> HypergraphState(5,[[1],[2],[3],[4],[5],[1,2],[2,3],[3,4],[4,5],[5,1]]) = WFiveQubitCodewordMinus;
# true
# gap> pentagon = WFiveQubitCodewordPlus;
# true

# my convention (order of bits is original n, then new 4n)
EncodeWithFiveQubitCodeOrig := function(psi)
    local n,encodedPsi,b,j,l;
    n := BitSize(psi);
    encodedPsi := psi;
    for j in [1..n] do
        encodedPsi := kron(encodedPsi,KronVec(replicate(4)([1,1])));
    od;
    for j in [1..n] do
        for l in [n + 1 + 4*(j-1)..n + 4 + 4*(j-1)] do
            encodedPsi := fastControlledNot(l,j,encodedPsi);
        od;
    od;
    for j in [1..n] do
        b := Concatenation([j],[n + 1 + 4*(j-1)..n + 4 + 4*(j-1)]);
        encodedPsi := fastGenControlledZs(
                              [[b[1],b[2]],
                               [b[2],b[3]],
                               [b[3],b[4]],
                               [b[4],b[5]],
                               [b[5],b[1]]],encodedPsi);
    od;
    return encodedPsi;
end;
# gap> l0:=WFiveQubitCodeword0;
# gap> l1:=WFiveQubitCodeword1;
# gap> PermuteQubitsPsi(EncodeWithFiveQubitCode(e00+e01+e10-e11),(2,6,5,4,3)) = kron(l0,l0) + kron(l0,l1) + kron(l1,l0) - kron(l1,l1);
# true

# my convention (order of bits is 1 original then 4 new, 1 original then 4 new, etc.)
EncodeWithFiveQubitCode := function(psi)
    local n,encodedPsi,b,j,l;
    n := BitSize(psi);
    encodedPsi := ordered_kron(Concatenation(replicate(n)([0,1,1,1,1])),
                          psi,KronVec(replicate(4*n)([1,1])));
    for j in [1..n] do
        for l in [5*(j-1) + 2..5*(j-1) + 5] do
            encodedPsi := fastControlledNot(l,5*(j-1)+1,encodedPsi);
        od;
    od;
    for j in [1..n] do
        b := 5*(j-1);
        encodedPsi := fastGenControlledZs(
                              [[b+1,b+2],
                               [b+2,b+3],
                               [b+3,b+4],
                               [b+4,b+5],
                               [b+5,b+1]],encodedPsi);
    od;
    return encodedPsi;
end;
# gap> l0:=WFiveQubitCodeword0;
# gap> l1:=WFiveQubitCodeword1;
# gap> EncodeWithFiveQubitCode(e00+e01+e10-e11) = kron(l0,l0) + kron(l0,l1) + kron(l1,l0) - kron(l1,l1);
# true
# gap> lp:=WFiveQubitCodewordPlus;
# gap> EncodeWithFiveQubitCode(HypergraphState(3,[[1,2,3]])) = KronVec([lp,lp,lp]) - 2 * KronVec([l1,l1,l1]);
# true

# Can five-qubit code be extended into a stabilizer state, or is it "stuck"?
extend5QubitCode:=function()
  local gen5,vecspace,pauli;
  for pauli in PauliGroup(5) do
    gen5:=Concatenation(FiveQubitCodeGenerators,[pauli]);
    vecspace:=VecSpaceStabByGroup(List(gen5,PauliMat));
    Print(pauli,"  ",Length(vecspace),"\n");
    if Length(vecspace) = 1 then return "!";
    fi;
  od;
end;
# it looks like it can be extended into a stabilizer state
# It absolutely can be.  pentagon lives in five-qubit code codespace.

# Nine-qubit Shor Code

ShorCodeGenerators:=[[1,3,3,0,0,0,0,0,0,0],
                     [1,0,3,3,0,0,0,0,0,0],
		     [1,0,0,0,3,3,0,0,0,0],
		     [1,0,0,0,0,3,3,0,0,0],
		     [1,0,0,0,0,0,0,3,3,0],
		     [1,0,0,0,0,0,0,0,3,3],
		     [1,1,1,1,1,1,1,0,0,0],
		     [1,0,0,0,1,1,1,1,1,1]];  # phase then nine qubits

# The two logical code words are a basis for this space, but not the one given by GAP
NineQubitCodeSpace:=function()
  return VecSpaceStabByGroup(List(ShorCodeGenerators,PauliMat));
end;

# extendToStabState:=function(gens)
#   n:=Length(gens[1])-1;  # first element is phase
#   vecspace:=VecSpaceStabByGroup(List(gens,PauliMat));
#   if
#     Length(vecspace) <> 2^(n-Length(gens))
#   then
#     return "given gens not independent";
#   fi;
#   for pauli in PauliGroup(n) do
#     newvecspace:=VecSpaceStabByGroup(List(gens,PauliMat));

extendOneGen:=function(gens)
  local n,vecspace,extendPred;
  n:=Length(gens[1])-1;  # first element is phase
  vecspace:=VecSpaceStabByGroup(List(gens,PauliMat));
  if
    Length(vecspace) <> 2^(n-Length(gens))
  then
    return "given gens not independent";
  fi;
  extendPred:=function(pge)
    local newvecspace;
    newvecspace:=VecSpaceStabByGroup(List(Concatenation(gens,[pge]),PauliMat));
    return (2 * Length(newvecspace) = Length(vecspace));
  end;
  return filter(PauliGroup(n),extendPred);
end;

### Are Dicke States determined by their 2-qubit RDMs?

d52rdm3 := rdm([1,2,3],DMn(Dicke(5,2)));
d52pooper := 9/10*DMn(Dicke(3,1))+1/10*DMn(e111);

### General Vector Space stuff

OrthogonalComplement:=compose(NullspaceMat)(TransposedMat);

###############################
# Positive-Semidefinite stuff #
###############################

NegateAlternateSigns:=function(lst)
  if lst=[] then
    return [];
  elif Length(lst)=1 then
    return lst;
  else
    return Concatenation([lst[1],-lst[2]],NegateAlternateSigns(lst{[3..Length(lst)]}));
  fi;
end;

andList:=function(lst)
  if lst=[] then
    return true;
  else
    return lst[1] and andList(cdr(lst));
  fi;
end;

isHermitian := rho -> rho = ConjugateTranspose(rho);

# Proposition 8.2.7 on page 463
# of Matrix Mathematics, second edition, by Dennis S. Bernstein
PositiveSemidefinite:=function(mat)
  return isHermitian(mat) and
    andList(List(NegateAlternateSigns(Reversed(CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(mat)))),x->x>=0));
end;

ExtremalDM := rho -> Determinant(rho) = 0 and PositiveSemidefinite(rho);

### Gram-Schmidt orthonormalization

GramSchmidt:=function(lst,inner)
  local helper;
  helper:=function(ortho,rest)
    local new,term;
    if IsEmpty(rest) then
      return ortho;
    else
      new:=rest[1];
      for term in ortho do
        new:=new - term * inner(term,rest[1]);
      od;
      new:=new/Sqrt(inner(new,new));
      return helper(Concatenation(ortho,[new]),cdr(rest));
    fi;
  end;
  return helper([],lst);
end;

# No normalization version
GramSchmidtNN:=function(lst,inner)
  local helper;
  helper:=function(ortho,rest)
    local new,term;
    if IsEmpty(rest) then
      return ortho;
    else
      new:=rest[1];
      for term in ortho do
        new:=new - term * inner(term,rest[1])/inner(term,term);
      od;
      return helper(Concatenation(ortho,[new]),cdr(rest));
    fi;
  end;
  return helper([],lst);
end;

###
# DotProduct() takes two vectors
# and returns their dot product
###

DotProduct:=function(cvec1,cvec2)
  local colvec1,colvec2,mat;
  colvec1:=TransposedMat([cvec1]);
  colvec2:=TransposedMat([cvec2]);
  mat:=ConjugateTranspose(colvec1)*colvec2;
  return mat[1][1];
end;

Overlap:=function(lst,inner)
  local mat,i,j;
  mat := NullMat(Length(lst),Length(lst));
  for i in [1..Length(lst)] do
    for j in [1..Length(lst)] do
      mat[i][j] := inner(lst[i],lst[j]);
    od;
  od;
  return mat;
end;


#################
### Werner stuff
#################

# gap> virginiaCycleList([1,1,1,1,0,0,0,0]);
# [ 1, 1, 1, 0, 0, 0, 0, 1 ]
# gap> PermuteList([1,1,1,1,0,0,0,0],GeneratorsOfGroup(CyclicGroup(IsPermGroup,8))[1]^(-1));
# [ 1, 1, 1, 0, 0, 0, 0, 1 ]

# period of a list
period := function(xs)
  local j,vxs;
  j := 1;
  vxs := virginiaCycleList(xs);
  while xs <> vxs do
    j := j + 1;
    vxs := virginiaCycleList(vxs);
  od;
 return j;
end;

RootCycle:=function(psi)
  local n,gen;
  n:=BitSize(psi);
  gen:=GeneratorsOfGroup(CyclicGroup(IsPermGroup,n))[1];
  return Sum(List([0..n-1],j->E(n)^j*PermuteQubitsPsi(psi,gen^j)));
end;

PauliWeights := composeList([Set,fmap(x -> CountNonZero(x[1])),Filter(x -> x[2] <> 0),PauliListForm3]);

PauliWeight := w -> composeList([DMFromPLF3,Filter(x -> CountNonZero(x[1]) = w),PauliListForm3]);

PW := w -> compose(PAF)(PauliWeight(w));

ExpandInOrthogonalBasis := innerProduct -> basis -> vec
  -> fmap(x -> innerProduct(x,vec)/innerProduct(x,x))(basis);

ExpandInPauliBasis := ExpandInOrthogonalBasis(HilbertSchmidtInnerProduct)([I,sx,sy,sz]);

# This will return an expansion even if it is not unique.
ExpandInBasis := basis -> vec -> SolutionMat(fmap(flatten)(basis),flatten(vec));

LinearCombinationsGivingZero := vecs -> NullspaceMat(fmap(flatten)(vecs));

OperatorFromBracket := orthogonalBasis -> w
  -> TransposedMat(fmap(x -> ExpandInOrthogonalBasis(HilbertSchmidtInnerProduct)(orthogonalBasis)(bracket(w,x)))(orthogonalBasis));

OneQubitLOperator := OperatorFromBracket([sx,sy,sz])(sx-i*sy);
OneQubitROperator := OperatorFromBracket([sx,sy,sz])(sx+i*sy);
OneQubitZOperator := OperatorFromBracket([sx,sy,sz])(sz);

# gap> Overlap([sx,sy,sz],HilbertSchmidtInnerProduct);
# [ [ 2, 0, 0 ], [ 0, 2, 0 ], [ 0, 0, 2 ] ]
# gap> Overlap([sr,sz,sl],HilbertSchmidtInnerProduct);
# [ [ 2, 0, 0 ], [ 0, 2, 0 ], [ 0, 0, 2 ] ]

# h's are normalized in the Hilbert-Schmidt norm
# good for density matrix basis states
h0 :=  I/Sqrt(2);
h1 := sx/Sqrt(2);
h2 := sy/Sqrt(2);
h3 := sz/Sqrt(2);

hr := sr/Sqrt(2);
hz := sz/Sqrt(2);
hl := sl/Sqrt(2);

# b's are normalized to give nice bracket relations
# good for use as operators
b0 :=  I/2;
b1 := sx/2;
b2 := sy/2;
b3 := sz/2;

normalizeHS := rho -> rho / Sqrt(HilbertSchmidtNorm(rho));

# using half the lowering operator to make the coefficients smaller
LowerOneStep := function(rho)
  local n,lowerOp;
  n := BitSize(rho);
  lowerOp := Sum(fmap(compose(Kron)(fmap(RZLMat)))(3*IdentityMat(BitSize(rho))));
  return bracket(lowerOp/2,rho);
end;

LowerListDM := function(rho)
  local result,term;
  term := rho;
  result := [];
  while not IsZero(term) do
    result := Concatenation(result,[term]);
    term := LowerOneStep(term);
  od;
  return result;
end;

LowerList := function(rho)
  local result,term;
  term := rho;
  result := [];
  while not IsZero(term) do
    result := Concatenation(result,[term]);
    term := LowerOneStep(term);
  od;
  return fmap(RZLAF)(result);
end;

LowerAndNormalizeOneStep := function(rho)
  local n,lowerOp;
  n := BitSize(rho);
  lowerOp := Sum(fmap(compose(Kron)(fmap(RZLMat)))(3*IdentityMat(BitSize(rho))));
  return normalizeHS(bracket(lowerOp,rho));
end;

CheckLList := function(lst)
  local j;
  for j in [1..Length(lst)] do
    Assert(0,IsInt(lst[j]));
    Assert(0,lst[j] >= 0);
    if j = 1 then
      Assert(0,lst[1] = 1);
    else
      Assert(0,AbsoluteValue(lst[j-1]-lst[j]) <= 1);
      if lst[j-1] = 0 then
        Assert(0,lst[j] = 1);
      fi;
    fi;
  od;
end;

MatrixOfInnerProducts := function(avecs,bvecs,inner)
  local result,row,col;
  result:=NullMat(Length(avecs),Length(bvecs));
  for row in [1..Length(avecs)] do
    for col in [1..Length(bvecs)] do
      result[row][col] := inner(avecs[row],bvecs[col]);
    od;
  od;
  return result;
end;

# linear combo of lst1 orthogonal to lst2
# assumes there is only one
LinearComboOfOrthogonalTo := function(dmlst1,dmlst2)
  local matIP,results;
  matIP := MatrixOfInnerProducts(dmlst1,dmlst2,HilbertSchmidtInnerProduct);
#  Display(matIP);
  results := NullspaceMat(matIP);
  Assert(0,Length(results) = 1);
  return Sum(zipWith(prod)(results[1])(dmlst1));
end;

IrrepsOfVlKronV1 := function(dmlst)
  local bigrep,medrep,smlrep;
  bigrep := LowerListDM(Kron([dmlst[1],sr]));
  if Length(dmlst) = 1 then
    return [bigrep];
  else
    medrep := LowerListDM(LinearComboOfOrthogonalTo([Kron([dmlst[1],sz]),Kron([dmlst[2],sr])],[bigrep[2]]));
    smlrep := LowerListDM(LinearComboOfOrthogonalTo([Kron([dmlst[1],sl]),Kron([dmlst[2],sz]),Kron([dmlst[3],sr])],[bigrep[3],medrep[2]]));
    return [bigrep,medrep,smlrep];
  fi;
end;

# lst = irrepDiagramList
IrrepBasis := function(lst)
  local bigmedsml;
  CheckLList(lst);
  if IsEmpty(lst) then
    return [[1]];  # zero qubits
  elif Length(lst) = 1 then
    return LowerListDM(sr);
  else
    bigmedsml := IrrepsOfVlKronV1(IrrepBasis(lst{[1..Length(lst)-1]}));
    if lst[Length(lst)]-lst[Length(lst)-1] = 1 then
      return bigmedsml[1];
    elif lst[Length(lst)]-lst[Length(lst)-1] = 0 then
      return bigmedsml[2];
    else
      return bigmedsml[3];
    fi;
  fi;
end;

#############################
# Werner State Verification #
#############################

# Projector onto states with an even number of zeros and ones.
BalancedProjector01 := function(n)
  return Sum(List(List(NChooseK(n,n/2),ee),DMn));
end;

HadamardAll := function(n)
  return Kron(replicate(n)(uhad));
end;

BalancedProjectorPM := function(n)
  return HadamardAll(n) * BalancedProjector01(n) * HadamardAll(n);
end;

BothBalancedBasis := function(n)
  return NullspaceMat(Kron(replicate(n)(I)) - BalancedProjectorPM(n) * BalancedProjector01(n));
end;

#############################################
# Werner pure state orthonormal basis ideas #
#############################################

# gap> List(CyclicGroup(IsPermGroup,4));
# [ (), (1,4,3,2), (1,3)(2,4), (1,2,3,4) ]

# gap> List(Stabilizer(CyclicGroup(IsPermGroup,4),sing12sing34-sing14sing23,PermuteQubitsPsi));
# [ (), (1,3)(2,4), (1,4,3,2), (1,2,3,4) ]
# gap> List(Stabilizer(CyclicGroup(IsPermGroup,4),sing12sing34+sing14sing23,PermuteQubitsPsi));
# [ (), (1,3)(2,4) ]

# 6-qubit werner

sing12sing34sing56 := normalize(KronVec([singlet,singlet,singlet]));
sing23sing45sing61 := PermuteQubitsPsi(sing12sing34sing56,(1,2,3,4,5,6));

sing14sing23sing56 := PermuteQubitsPsi(sing12sing34sing56,(2,4,3));
sing36sing45sing12 := PermuteQubitsPsi(sing14sing23sing56,(1,2,3,4,5,6)^2);
sing52sing61sing34 := PermuteQubitsPsi(sing14sing23sing56,(1,2,3,4,5,6)^4);

werner620 := sing12sing34sing56 + sing23sing45sing61;
werner621 := sing12sing34sing56 - sing23sing45sing61;

werner630 := sing14sing23sing56 +          sing36sing45sing12 +          sing52sing61sing34;
werner631 := sing14sing23sing56 + E(3)   * sing36sing45sing12 + E(3)^2 * sing52sing61sing34;
werner632 := sing14sing23sing56 + E(3)^2 * sing36sing45sing12 + E(3)   * sing52sing61sing34;

werner6basis := [werner620,werner621,werner630,werner631,werner632];

# gap> PermuteQubitsPsi(werner620,(1,2,3,4,5,6)) = werner620;
# true
# gap> PermuteQubitsPsi(werner621,(1,2,3,4,5,6)) = -werner621;
# true
# gap> PermuteQubitsPsi(werner630,(1,2,3,4,5,6)) = -werner630;
# true
# gap> PermuteQubitsPsi(werner631,(1,2,3,4,5,6)) = E(6)^5*werner631;
# true
# gap> PermuteQubitsPsi(werner632,(1,2,3,4,5,6)) = E(6)*werner632;
# true

werner6a := normalize(werner620);
werner6b := normalize(werner621+2/3*werner630);
werner6c := normalize(werner630);
werner6d := normalize(werner631);
werner6e := normalize(werner632);

# gap> Overlap([werner6a,werner6b,werner6c,werner6d,werner6e],DotProduct);
# [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ]

# gap> PermuteQubitsPsi(werner630,(2,4,6)) = werner630;
# true

###################################
# Positive semidefinite searching #
###################################

PSDSearch := function(dm1,dm2,dm3,plist,qlist,basefilename)
  local p,q;
  for p in plist do
    for q in qlist do
      if PositiveSemidefinite(p*dm1+q*dm2+(1-p-q)*dm3) then
        AppendTo(Concatenation(basefilename,"-psd-true.dat"),p," ",q,"\n");
      else
        AppendTo(Concatenation(basefilename,"-psd-false.dat"),p," ",q,"\n");
      fi;
    od;
  od;
end;

#####################
# Hypergraph states #
#####################

###
# HypergraphState(n, [ [] ]) where n is the number of qubits and hyperedgeList is a list of lists. Each sublist is the qubits
# connected by a hyperedge
# For example, the input (4, [ [1,2], [3,4] ]) is a graph with 4 qubits, where [1,2] are connected by and edge, and [3,4] are connected by an edge
#
# HypergraphState() returns the row vector of 1s and -1s that describes the state
###

HypergraphState := function(n,hyperedgeList)
    return fastGenControlledZs(hyperedgeList,ListWithIdenticalEntries(2^n,1));
end;

###
# Given a state vector, Hypergraph() returns the list of lists that describe the edges of the hypergraph
# The state vector will be a list of 1s and -1s
###
Hypergraph := function(cvec)
    local cvec1,n,hyperedge,hyperedgeList;
    # might want some checks here
    hyperedgeList := [];
    n := BitSize(cvec);
    cvec1 := cvec;
    for hyperedge in SubsetList(n) do
        if cvec1[1 + Bin(QubitSet2BinarySubsystem(n,hyperedge))] < 0 then
            Add(hyperedgeList,hyperedge);
            cvec1 := fastGenControlledZ(hyperedge,cvec1);
        fi;
    od;
    return hyperedgeList;
end;

# HyperedgeSets := n -> Combinations(SetMinus(Combinations([1..n]))([[]]));
# HyperedgeSets := n -> CombinationsBySize(SetMinus(CombinationsBySize([1..n]))([[]]));
HyperedgeSets := n -> SortBySize(List(Combinations(SetMinus(CombinationsBySize([1..n]))([[]])),SortBySize));

SymmetricHypergraphState := function(n,completenessLevels)
    return HypergraphState(n,
                   Concatenation(
                           List(completenessLevels,
                                level -> Combinations([1..n],level))));
end;

GenControlledZ := function(n,qubitset)
    return DiagonalMat(HypergraphState(n,[qubitset]));
end;

# The product of C_e gates that forms a hypergraph state from |+>^n
HypergraphC := function(n,hyperedgeList)
    return DiagonalMat(HypergraphState(n,hyperedgeList));
end;

SymmetricHypergraphC := function(n,completenessLevels)
    return DiagonalMat(SymmetricHypergraphState(n,completenessLevels));
end;
