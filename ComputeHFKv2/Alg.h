#ifndef ALGARROW_H
#define ALGARROW_H
#include<vector>
#include<array>
#include<map>
#include<unordered_map>
#include<cstdint>
#include "Diagrams.h"

using namespace std;

//////////////////////
/// 1. KNOT INVARIANTS
//////////////////////

struct KnotFloerGen
{int Name; int Maslov; int Alexander;
};

struct ChainArrow
{int StartingGen; int EndingGen; int Coeff;
};

struct KnotFloerComplex
{vector<KnotFloerGen> Generators; vector<ChainArrow> Differential; int Prime;
}; 

  //Computes U\cdot V=0 version from a MorseList over F_{Prime}
KnotFloerComplex ComputingKnotFloer(MorseCode, int Prime, bool verbose=true);
  
int Genus(const KnotFloerComplex &);

bool Fibered(const KnotFloerComplex &); 

bool LSpaceKnot(const KnotFloerComplex &);

int Tau(const KnotFloerComplex &);

int Nu(const KnotFloerComplex &);

int Epsilon(const KnotFloerComplex &);

  //computes bigraded ranks:
map<pair<int, int>,int > KnotFloerRanks(const KnotFloerComplex &);

  //prints bigraded ranks:
void ReportKnotFloerRanks(const KnotFloerComplex &, ostream & os); 

  //An alternative  method for alternating projections:
void  KnotFloerForAlternatingKnots(PlanarDiagram, ostream & os);

struct ChainComplex
{vector<int> Generators; vector<ChainArrow> Differential; int Prime;
};

int HomologyRank(const ChainComplex & OldComplex);

/////////////////////////////
// 2. ACTIONS ON D-STRUCTURES
/////////////////////////////

  // Tensors with a Maximum bimodule.
  // Maximum is added to right of strand Position-1.
  // The new maximum is between the strands Position, Position+1.
  // Orientation=0, maximum is oriented left to right
  // Orientation=1, maximum is oriented right to left
void AfterMax(int Position, int Orientation);

  // Tensors with the minimum bimodule, where  minimum is added on the left
  // between strands 1 and 2: 
void AfterMin(); 

  // Tensors with a crossing bimodule
  // P^{C} if C>0  or N^{-C} if C<0:
void AfterCrossing(int Crossing); 

  // Repeatedly contracts short arrows
  // until there are no short differentials
  // in the D-structure:
void Simplify(); 

void Checking();

///////////////////////////////////////
// 3. EXTERNAL VARIABLES, STRUCTURES:
///////////////////////////////////////

  // Can be used for knots with  girth number <= 20: 
static const int MAXBRIDGE=10;     

  // Working over the field F_{Modulus}:
extern int Modulus;

  // Ex: idempotent {1,3,5} in algebra B(6,3) is represented as  2+8+32. 
  // Holds values between 0 and 2^(2*MAXBRIDGE)-1:
typedef uint32_t idem;       

  // Represents a monomial in F[U_1,...,U_{20}]. 
  // Ex: {1,0,2,0,...0} = U_1\times U_3^2 : 
typedef array<unsigned char, 2*MAXBRIDGE> monomial; 

  // Corresponds to 1 in F[U_1,...,U_{20}]:
extern const monomial MonomialOne; 

  // Each monomial that the program encounters has an integer identifier. 
  // Ex: 0 corresponds to MonomialOne. MonomialStore[0]=MonomialOne:
extern vector<monomial> MonomialStore; 

struct Hash
{  int operator() (const monomial& M) const
  {int ans=0; int i=2*MAXBRIDGE; while(i--) ans=(ans<<3)+M[i]; return ans; }
};                                     

  // Used for looking up the integer identifier:
extern unordered_map<monomial, int, Hash> MonomialMap; 

  // Finds the integer that represents X. 
  // If X is a new monomial updates MonomialStore and MonomialMap.
int MonomialLookUp(const monomial& X); 

  // The algebra elements are in B(2*Bridge,Bridge). 
  // Bridge is updated when we go through Minimums and Maximums:
extern int Bridge;               

  // In UpwardList 1 means upward, 0 downward. Ex: {1,0,0,1} means 
  // strands 1 and 4 are oriented upward, 
  // 2 and 3 downward.  Used for local Alexander and Maslov grading updates:
extern vector<int> UpwardList;   

  // Represents the current matching M in algebra A(Bridge, M). 
  // Ex: {4,3,2,1,6,5},  4 is matched with 1, 3 with 2, 5 with 6:
extern vector<int> MatchingList; 

  // One of the generators of a D-structure:
struct Gen { 
   int Name=-1; idem Idem=0; int Maslov=0; int Alexander=0; 
};

  // A compact representation of a term in delta ^1 in a D structure.
  //There could be multiple arrows between the same generators
  //Coeff is between 0 and Modulus-1:
struct Arrow {
    int StartingGen; int EndingGen; int MonomialIndex; int Coeff; 
};

  // The list of all the generators of a D structure
  // After updates (tensoring with bimodules or
  // simplifying) GeneratorList[i].Name should be i: 
extern vector<Gen> GeneratorList;

  // A temporary list of generators. 
  // Used while tensoring with bimodules.
  // After tensoring  cycles NewGeneratorList should be empty:
extern vector<Gen> NewGeneratorList;

  // The list of all the terms in delta ^1.
  // After tensoring with bimodules or simplifying: 
  // if there are multiple terms with 
  // a fixed StartingGen, EndingGen
  // and MonomialIndex, they are combined by adding up
  // their MonomialIndex mod Modulus.
  // Terms with MonomialIndex=0 are deleted:  
extern vector<Arrow> ArrowList;

  // A temporary list of arrows.
extern vector<Arrow> NewArrowList;

////////////////////////
// 4. UTILITY FUNCTIONS
//////////////////////// 

  // I1, I2 and MonIndex represents a pure algebra element
  // in B_0(2*Bridge, Bridge). NonZero checks if its image 
  // in the quotient algebra B(2*Bridge,Bridge) is non-trivial:  
bool NonZero(idem I1 , idem I2, int MonIndex);

  // At least one coordinate moves more than 1 unit, Ex: R2 \times R3. 
bool TooFar(idem I1, idem I2); 

  // A raw version of multiplication. 
  // It is used to compute the MonomialIndex of a\cdot b,
  // where "a" is represented by I1,I2,MonIndex1,
  // "b" is represented by I2,I3,MonIndex2, I1 and I2 are not TooFar,
  // I2 and I3 are not TooFar.
  // Note that I1 and I3 are allowed to be TooFar:
int Mult(idem I1, idem I2, idem I3, int MonIndex1, int MonIndex2); 

  // The number of coordinates left of the Wall:
int LeftCoord(idem I1, int Wall);

  // equals to 1 on R_Wall, -1 on L_Wall, 0 on U_Wall:
int LeftRight(idem I1, idem I2, int Wall);

  // Compares StartingGen, EndingGen, MonomialIndex (but not Coeff):
bool Strict(const Arrow & arrow1, const Arrow & arrow2); 

  // (StartingGen, EndingGen, MonomialIndex) same for arrow1 and arrow2:
bool Equal(const Arrow & arrow1, const Arrow & arrow2);  

  // Compares StartingGen: 
bool FirstLess(const Arrow & arrow1, const Arrow & arrow2);  

  // Compares EndingGen:
bool SecondLess(const Arrow & arrow1, const Arrow & arrow2); 

  // Combines multiple arrows between same generators and 
  // same algebra output. 
  // New coefficient is the sum mod Modulus. 
  // Terms with Coeff=0 are deleted. Used by all the tensor product operations
  // and  simplifying operations:
void RemoveMod(vector<Arrow> &); 

  //Renames the generators so that GeneratorList[i].Name=i, updates ArrowList:
void ReName();

template<typename T>
inline int sizeAsInt(const T& v) { return static_cast<int>(v.size()); }

#endif
