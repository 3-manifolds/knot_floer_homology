#include "Alg.h"
#include "Diagrams.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// Global variables used by functions called from this module.

vector<monomial> MonomialStore;
unordered_map <monomial,int,Hash> MonomialMap;
int Bridge;
int Modulus;
vector<int>  UpwardList;
vector<int>  MatchingList;
vector<Arrow> ArrowList;
vector<Arrow> NewArrowList;
vector<Gen>   GeneratorList;
vector<Gen>   NewGeneratorList;
monomial MonomialOne={0};

// Unexported declarations from Alternating.cpp 

struct Term{
  idem Idem;
  int Alexander;
  int Coeff;
};
extern vector<Term> AfterMaxAlt(vector<Term> Old, int Position);
extern vector<Term> AfterMinAlt(vector<Term> Old);
extern vector<Term> AfterCrossingAlt(vector<Term> Old, int Crossing);
extern int Signature (PlanarDiagram Diag);
  
// static helpers

#define ITEM(stream, key, value)					\
  stream << "  \"" << (key) << "\": " << (value) << "," << endl;

static bool isPrime(int n) {
  if (n < 2) {
    return false;
  }
  for (int i = 2; i < n; i++) {
    if (n % i == 0) {
      return false;
    }
  }
  return true;
}

static void createCString(char **p, ostringstream& s) {
  if (p) {
    *p = (char *) malloc(s.str().length() + 1);
    strcpy(*p, s.str().c_str());
  }
}

static void readPDCode(char *pd, string& s) {
  char *p;
  for (p = pd; *p != 0; p++) {
    if (*p == 'D') {
      break;
    }
  }
  for (; *p != 0; p++) {
    s.push_back(*p);
  }
}

static void KnotFloerRanksAsDict(const KnotFloerComplex& KFC, ostream& os)
{
  auto Map = KnotFloerRanks(KFC);

  os << "{" << endl;
  for(auto X: Map) {
    os << "  (" << X.first.first << "," << X.first.second << "): ";
    os << X.second << "," << endl;
  }
  os << "}" << endl;
}


static void MorseCodeAsEvents(MorseCode& code, ostream& os)
{
    int c;
    vector<int> Morse = code.GetMorseList();
    os << "{" << endl;
    ITEM(os, "girth", code.GetGirth());
    os << "  \"events\": [";
    for (size_t i = 0; i < Morse.size(); i++){
	if (Morse[i] >999){
	    c = Morse[++i];
	    os << "('cup', " << c - 1 << ", " << c << ")," << endl;
	}
	else if (Morse[i] >-1000){
	    c = Morse[i];
	    if (c > 0){
		os << "('cross', " << c - 1 << ", " << c << "), " << endl;
	    }
	    else{
                c = abs(c);
		os << "('cross', " << c << ", " << c - 1 << "), " << endl;
	    }
	}
	else{
	    os << "('cap', 0, 1)," << endl;
	}
    }
    os << "]}" << endl;
}

// Variant of KnotFloerForAlternatingKnots with a different output format.

static void KnotFloerForAlternatingKnotsAsDict(PlanarDiagram Diag, ostream& os) {
  vector<int> Morse = Diag.GetSmallGirthMorseCode(200).GetMorseList();
  Bridge=1;
  Term G1; G1.Alexander = 0; G1.Coeff = 1; G1.Idem = 2;
  vector<Term> Current; Current.push_back(G1);
  UpwardList.clear();
  if(Morse[0]==1000) {
    UpwardList.push_back(1);
    UpwardList.push_back(0);
  }
  if(Morse[0]==1001) {
    UpwardList.push_back(0);
    UpwardList.push_back(1);
  }
  int Steps=Morse.size();
  for(int i=2; i< Steps -1 ; i++) {
    if (Morse[i] == 1000) {
      int Position = Morse[i+1];
      Current = AfterMaxAlt(Current,Position);
      UpwardList[Position-1]=1;
      UpwardList[Position]=0;
      i++;
    } else if (Morse[i] == 1001) {
      int Position = Morse[i+1];
      Current = AfterMaxAlt(Current,Position);
      UpwardList[Position-1] = 0;
      UpwardList[Position] = 1;
      i++;
    } else if (Morse[i]==-1000) {
	Current=AfterMinAlt(Current);
    } else if (Morse[i] < 2*Bridge && Morse[i] > -(2*Bridge) && Morse[i] != 0) {
	Current=AfterCrossingAlt(Current, Morse[i]);
    }
  }

  int delta=-Signature(Diag)/2;
  map<int,int> Range;
  for(size_t i = 0; i < Current.size(); i++) {
    Term G = Current[i];
    Range[G.Alexander] = G.Coeff;
  }
  int TotalRank=0;
  bool LSpaceKnot=true;
  os << "{" << endl;
  os << "  \"ranks\": {" << endl;
  for(auto X: Range) {
    int a = X.first, b = abs(X.second);
    os << "    (" << a/2 << "," << a/2-delta << "): " << b << "," << endl;
    TotalRank += b;
    if (b != 1) {
      LSpaceKnot = false;
    }
  }
  os << "  }," << endl;
  int MaxAlex = -(*Range.begin()).first;
  int LeadingCoeff = (*Range.begin()).second;
  ITEM(os, "total_rank", TotalRank);
  ITEM(os, "seifert_genus", MaxAlex/2);
  ITEM(os, "fibered", (LeadingCoeff == 1 || LeadingCoeff == -1) ? "True" : "False"); 
  ITEM(os, "L_space_knot", LSpaceKnot ? "True" : "False");
  ITEM(os, "tau", delta);
  ITEM(os, "nu", delta);
  ITEM(os, "epsilon", (delta > 0) - (delta < 0)); 
  os << "}" << endl;
}

// The one and only function exported by this module.

extern "C" void
PDCodeToMorseAndHFK(
  char *pd,
  int prime,
  char **morse,
  char **hfk,
  char **error)
{
  ostringstream Error, Morse, HFK;
  string PDString;
  bool inputIsInvalid = false;

  if (error == NULL) {
    throw runtime_error("PDCodeToMorseAndHFK: error must not be null.");
  }
  MonomialStore.push_back(MonomialOne);
  MonomialMap.insert(make_pair(MonomialOne, 0));
  readPDCode(pd, PDString);
  PlanarDiagram Diag = PlanarDiagram(PDString);
  if (!isPrime(prime)) {
    Error << prime << " " << "is not a prime";
    inputIsInvalid = true;
  }
  if (!inputIsInvalid && Diag.NotValid()) {
    Error << "The PD code does not describe a knot projection.";
    inputIsInvalid = true;
  }
  if (!inputIsInvalid && Diag.R1Reducible()) {
    Error << "The PD code describes a knot projection with a reducing Reidemeister 1 move";
    inputIsInvalid = true;
  }
  if (!inputIsInvalid) {
    MorseCode LastCheckBeforeComputation = Diag.GetSmallGirthMorseCode(1);
    if (LastCheckBeforeComputation.GetMorseList().size() == 0) {
      Error << "Could not compute a small girth Morse code";
      inputIsInvalid = true;
    }
  }
  if (!inputIsInvalid) {
    if(Diag.Alternating()) {
      MorseCode M = Diag.GetSmallGirthMorseCode(10);
      MorseCodeAsEvents(M, Morse);
      if (hfk) {
	KnotFloerForAlternatingKnotsAsDict(Diag, HFK);
      }
    } else {
      MorseCode M = Diag.GetSmallGirthMorseCode();
      if (M.GetGirth() > 2*MAXBRIDGE) {
	Error << "Girth number exceeds " << 2*MAXBRIDGE;
      } else {
	MorseCodeAsEvents(M, Morse);
	if (hfk) {
	    KnotFloerComplex KFC=ComputingKnotFloer(M, prime);
	    HFK << "{" << endl;
	    ITEM(HFK, "modulus", prime);
	    HFK << "  \"ranks\": ";
	    KnotFloerRanksAsDict(KFC, HFK);
	    HFK << "," << endl;
	    ITEM(HFK, "total_rank", KFC.Generators.size());
	    ITEM(HFK, "seifert_genus", Genus(KFC));
	    ITEM(HFK, "fibered", (Fibered(KFC) ? "True" : "False"));
	    ITEM(HFK, "L_space_knot", (LSpaceKnot(KFC) ? "True" : "False"));
	    ITEM(HFK, "tau", Tau(KFC));
	    ITEM(HFK, "nu", Nu(KFC));
	    ITEM(HFK, "epsilon", Epsilon(KFC));
	    HFK << "}" << endl;
	}
      }
    }
  }
  createCString(error, Error);
  createCString(morse, Morse);
  createCString(hfk, HFK);
}
