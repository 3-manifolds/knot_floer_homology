#include "Alg.h"
#include "Diagrams.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

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

extern "C" void
PDCodeToMorseAndHFK(
  char *pd,
  int Prime,
  char **morse,
  char **hfk,
  char **error)
{
  ostringstream Error, Morse, HFK;
  bool isPrime = true;
  string S;
  char *p;

  if (Prime < 2) {
    isPrime = false;
  } else {
    for (int i = 2; i < Prime; i++) {
      if (Prime % i == 0) {
	isPrime = false;
	break;
      }
    }
  }
  if (!isPrime) {
    Error << Prime << " " << "is not a prime";
  }
  for (p = pd; *p != 0; p++) {
    if (*p == 'D') {
      break;
    }
  }
  for (; *p != 0; p++) {
    S.push_back(*p);
  }
  PlanarDiagram Diag = PlanarDiagram(S);
  if (Diag.NotValid()) {
    Error << "The PD code does not describe a knot projection.";
  }
  if (Diag.R1Reducible()) {
    Error << "The PD code describes a knot projection with a reducing Reidemeister 1 move";
  }
  MorseCode LastCheckBeforeComputation = Diag.GetSmallGirthMorseCode(1);
  if (LastCheckBeforeComputation.GetMorseList().size() == 0) {
    Error << "Could not compute a small girth Morse code";
  }
  if (Error.str().empty()) {
    if(Diag.Alternating()) {
      MorseCode M = Diag.GetSmallGirthMorseCode(10);
      M.Print(Morse);
      if (hfk) {
	KnotFloerForAlternatingKnots(Diag, HFK);
      }
    } else {
      MorseCode M = Diag.GetSmallGirthMorseCode();
      if (M.GetGirth() > 2*MAXBRIDGE) {
	Error << "Girth number exceeds " << 2*MAXBRIDGE;
      } else {
	M.Print(Morse);
	KnotFloerComplex KFC=ComputingKnotFloer(M, Prime);
	ReportKnotFloerRanks(KFC, HFK);
	HFK << "Total rank : " << KFC.Generators.size() << endl;
	int genus = Genus(KFC);
	HFK << "Seifert genus : " << genus << endl;
	HFK << "Fibered : " << (Fibered(KFC) ? "Yes" : "No") << endl;
	HFK << "L-space knot : " << (LSpaceKnot(KFC) ? "Yes" : "No") << endl;
	HFK << "Tau : " << Tau(KFC) << endl;
	HFK << "Nu  : " << Nu(KFC) << endl;
	int epsilon = Epsilon(KFC);
	HFK << "Epsilon : " << Epsilon(KFC) << endl;
      }
    }
  }
  string errorString = Error.str();
  *error = (char *) malloc(errorString.length() + 1);
  strcpy(*error, errorString.c_str());
  if (morse) {
    string morseString = Morse.str();
    *morse = (char *) malloc(morseString.length() + 1);
    strcpy(*morse, morseString.c_str());
  }
  if (hfk) {
    string hfkString = HFK.str();
    *hfk = (char *) malloc(hfkString.length() + 1);
    strcpy(*hfk, hfkString.c_str());
  }
}
