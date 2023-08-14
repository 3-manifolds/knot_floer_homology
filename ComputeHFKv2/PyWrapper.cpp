#include "Alg.h"
#include "Diagrams.h"
#include "PyWrapper.h"
#include "PyObjectLight.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <mutex>

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
const monomial MonomialOne={0};

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

static std::once_flag _monomialStoreAndMapInitialized;

static void _InitializeMonomialStoreAndMap()
{
    std::call_once(
        _monomialStoreAndMapInitialized,
        []() {
            // Shouldn't we clear this before?
            MonomialStore.push_back(MonomialOne);
            MonomialMap.insert(make_pair(MonomialOne, 0));
        });
}

// static helpers

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

static py::object MorseListAsEvents(const vector<int> &morseList)
{
    std::vector<py::object> events;

    for (int i = 0; i < sizeAsInt(morseList); i++){
        if (morseList[i] >999){
            const int c = morseList[++i];
            events.push_back(py::object("cup", c - 1, c));
        } else if (morseList[i] >-1000){
            const int c = morseList[i];
            if (c > 0){
                events.push_back(py::object("cross", c - 1, c));
            } else {
                events.push_back(py::object("cross", - c, - c - 1));
            }
        } else {
            events.push_back(py::object("cap", 0, 1));
        }
    }

    return events;
}

static py::object MorseCodeAsEvents(const MorseCode &code)
{
    return std::map<std::string, py::object>{
        { "events", MorseListAsEvents(code.GetMorseList()) },
        { "girth", code.GetGirth() } };
}

// The one and only function exported by this module.

PyObject *PDCodeToHFK(const char *pd, int prime, bool complex)
{
  const PlanarDiagram diag = PlanarDiagram(pd);

  if (diag.NotValid()) {
      py::RaiseValueError(
          "The PD code does not describe a knot projection.");
      return nullptr;
  }
  
  if (diag.R1Reducible()) {
      py::RaiseValueError(
          "The PD code describes a knot projection with a reducing "
          "Reidemeister 1 move");
      return nullptr;
  }

  const MorseCode lastCheck = diag.GetSmallGirthMorseCode(1);
  if (lastCheck.GetMorseList().size()==0) {
      py::RaiseValueError("This PD code cannot be handled, possibly a connected sum?");
      return nullptr;
  }

  if (prime > 32768) {
            py::RaiseValueError("The prime " + std::to_string(prime) +
                                " is greater than 2^15.");
      return nullptr;
  }

  
  if (!isPrime(prime)) {
      py::RaiseValueError(std::to_string(prime) + " is not prime");
      return nullptr;
  }

  _InitializeMonomialStoreAndMap();

  const MorseCode M = diag.GetSmallGirthMorseCode();
  if (M.GetMorseList().empty()) {
      py::RaiseValueError("Could not compute a small girth Morse code");
      return nullptr;
  }
  if (M.GetGirth() > 2*MAXBRIDGE) {
      py::RaiseValueError(
	  "Girth number exceeds " + std::to_string(2 * MAXBRIDGE));
      return nullptr;
  } else {
      KnotFloerComplex KFC = ComputingKnotFloer(M, prime, false);
      if (KFC.Prime == 0) {
	py::RaiseRuntimeError("Interrupted!");
	return nullptr;
      }
      if (complex){
	  py::object o(
	      std::map<std::string, py::object>{
		  { "modulus", prime },
		  { "ranks", KnotFloerRanks(KFC) },
		  { "total_rank", KFC.Generators.size() },
		  { "seifert_genus", Genus(KFC) },
		  { "fibered", Fibered(KFC) },
		  { "L_space_knot", LSpaceKnot(KFC) },
		  { "tau", Tau(KFC) },
		  { "nu", Nu(KFC) },
		  { "epsilon", Epsilon(KFC) },
		  { "generators", KnotFloerGenerators(KFC) },
		  { "differentials", KnotFloerDifferentials(KFC)}});
	  return o.StealObject();
      }
      else{
	  py::object o(
	      std::map<std::string, py::object>{
		  { "modulus", prime },
		  { "ranks", KnotFloerRanks(KFC) },
		  { "total_rank", KFC.Generators.size() },
		  { "seifert_genus", Genus(KFC) },
		  { "fibered", Fibered(KFC) },
		  { "L_space_knot", LSpaceKnot(KFC) },
		  { "tau", Tau(KFC) },
		  { "nu", Nu(KFC) },
		  { "epsilon", Epsilon(KFC) }});
	  return o.StealObject();
      }
  }
}

PyObject *PDCodeToMorse(const char *pd)
{
    const PlanarDiagram diag = PlanarDiagram(pd);

    if (diag.NotValid()) {
        py::RaiseValueError(
            "The PD code does not describe a knot projection.");
        return nullptr;
    }

    if (diag.R1Reducible()) {
        py::RaiseValueError(
            "The PD code describes a knot projection with a reducing "
            "Reidemeister 1 move");
        return nullptr;
    }

    const MorseCode lastCheck = diag.GetSmallGirthMorseCode(1);
    if (lastCheck.GetMorseList().size()==0) {
	py::RaiseValueError(
            "This PD code cannot be handled, possibly a connected sum?");
        return nullptr;
    }
    
    const bool alternating = diag.Alternating();

    // Trying to mimick the behavior PDCodeToHFK...
    const MorseCode morseCode = 
         alternating ? diag.GetSmallGirthMorseCode(10)
                     : diag.GetSmallGirthMorseCode();

    if (morseCode.GetMorseList().empty()) {
        py::RaiseValueError(
            "Could not compute a small girth Morse code");
        return nullptr;
    }

    if (alternating) {
        if (morseCode.GetGirth() > 2 * MAXBRIDGE) {
            py::RaiseValueError(
                "Girth number exceeds " + std::to_string(2 * MAXBRIDGE));
            return nullptr;
        }
    }

    py::object o = MorseCodeAsEvents(morseCode);
    return o.StealObject();
}
