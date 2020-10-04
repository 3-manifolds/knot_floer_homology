The program uses Planar Diagram codes:
For example Trefoil = PD[X[4,2,5,1], X[2,6,3,5], X[6,4,1,3]]
It works for planar diagrams of knots.
A version for links will also be available later.

The X and the extra brackets are optional, the program will look for the PD start
and then read in the integers separated by other char, so 
PD[4,2,5,1,2,6,3,5,6,4,1,3], 
will also work.

You can give it several knots, just include PD in front of each:
PD[4,2,5,1,2,6,3,5,6,4,1,3], 
PD[12,2,13,1,7,3,8,2,3,15,4,14,8,14,9,13,15,11,16,10,4,10,5,9,11,7,12,6,16,6,1,5]

Compile with:
sh compile.sh

or equivalently:
g++ -std=c++11 -O3 Main.cpp Utility.cpp Min.cpp Max.cpp Crossing.cpp Simplify.cpp HomologyRank.cpp Report.cpp Diagrams.cpp Alternating.cpp KnotFloer.cpp 

It uses C++11,  and will produce an a.out file.
Run it  with ./a.out (depending on your setup a.out may also work).
It will ask for the name of an input file and a prime number p.
The answer will print on the screen, and to inputfile.modp, together
with additional info to inputfile.modp.Morse.
The second file contains a Morse presentation of the knot (that the program found and
used for the planar diagram). 

In the Morse presentation minimums are always added on the left. Max(i) indicates that the maximum was 
added to create the strands i and i+1. Crossing n means right-handed twist between strands n and n+1, and -n
the left handed twist.
Alternating projections are special, here we use the Alexander polynomial and signature.

For example Torus34.txt and 7 
will result in two output files:
Torus34.txt.mod7
Torus34.txt.mod7.Morse

