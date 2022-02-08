#include <iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<iostream>
#include <utility>  // pair
#include "Alg.h"
#include "Diagrams.h"

using namespace std;

void MorseCode::Print(ostream & os) const {
  os<<"Morse Code: ";
  for (int i = 0; i < sizeAsInt(MorseList); i++)
    if (MorseList[i] >999)
      os<<"Max("<< MorseList[++i]<<"), ";
    else if (MorseList[i] >-1000)
      os<<MorseList[i]<<", ";
    else if (i < sizeAsInt(MorseList)-1)
      os<<"Min, ";   
    else 
      os<<"Min.\n";
  os<<"Girth: "<<Girth<<"\n\n";
}

void PlanarDiagram::Print(ostream & os) const {
  const vector<int> &PD = ListOfTuples;
  os<<"Planar Diagram: [";
  for (int i = 0; i < sizeAsInt(PD)/4; i++)
    os<<"["<<PD[4*i]<<", "<<PD[4*i+1]<<", "<<PD[4*i+2]<<", "<<PD[4*i+3]<<"], ";
  os<<"]\n";
} 

bool PlanarDiagram::NotValid() const { //partial check on data, also checks if it has more than 1 component
  const vector<int> &PD = ListOfTuples; 
  int y = sizeAsInt(PD);
  if (y == 0) 
    return true;
  if (y%4 != 0) 
    return true;
  vector<int> CrList = PD;
  sort(CrList.begin(), CrList.end());
  int x=y/4;
  for (int i = 0; i < 2*x; i++)
    if (CrList[i] != (i/2)+1) 
      return true;
  for (int i = 0; i < x; i++)
    if (PD[4*i]-PD[4*i+2] != -1 &&  PD[4*i]-PD[4*i+2] != 2*x-1) 
      return true;
    else if (abs(PD[4*i+1] -PD[4*i+3]) != 1 && abs(PD[4*i+1]-PD[4*i+3]) != 2*x-1)
      return true; 
  return false;
}
  
bool PlanarDiagram::Alternating() const {
  const vector<int> &PD = ListOfTuples;  
  int x = sizeAsInt(PD)/4; 
  int a = PD[0]%2;
  for (int i = 1; i < x; i++)
    if (a != (PD[4*i]%2)) 
      return false;
  return true;
}
  
bool PlanarDiagram::R1Reducible() const {
  const vector<int> &PD = ListOfTuples;
  int x = sizeAsInt(PD)/4;
  for (int i = 0; i < x; i++)
    if (PD[4*i] == PD[4*i+1] || PD[4*i+1] == PD[4*i+2] || PD[4*i+2]==PD[4*i+3] || PD[4*i+3]==PD[4*i]) 
      return true;
  return false;
}

PlanarDiagram::PlanarDiagram(const string &S){
  ListOfTuples = vector<int>();  
  int x = '0'; 
  int y = '9';
  int a = 0; 
  bool BuildingInt = false; 
  for (int i = 0; i < sizeAsInt(S); i++) 
      if (S[i] < x || S[i] > y) {
	if (BuildingInt == true) { //just finished reading an integer  
	  ListOfTuples.push_back(a);  a = 0; BuildingInt = false;
	}
      }
      else {  
          a = 10*a+(S[i]-x); BuildingInt = true;
      }
    
    if (ListOfTuples.size() > 0) { //normalizing the indexes for the edges
      int smallest=ListOfTuples[0];
      for (int i: ListOfTuples) 
        if (i < smallest) smallest = i;
      for (int & i: ListOfTuples) i = i-smallest+1;
    }
}

/* Auxiliary function for PlanarDiagram::GetSmallGirthMorseCode.
 * 
 * Find the crossings, not yet added to the Morse list, that share a maximal
 * number of edges with the bottom edges of the upper knot diagram.
 */
inline pair< int, vector< int > > GetMaxConnections(
  const vector< int >& PD,
  const vector< bool >& Added,
  const vector< int >& Edges
) {
  vector< int > MaxConnectedCrossings;
  int MaxConnections = 1;
  
  // Find maximally connected crossings
  for (int Crossing = 0; Crossing < (sizeAsInt(PD) / 4); ++Crossing) {
    if (Added[Crossing]) {
      continue;
    }
    
    /* Vectors of positions of edges that belong to the bottom edges of the
     * upper knot diagram and to the current crossing.
     * 
     * These vectors are static for memory optimization.
     */
    static vector< int > PosInEdges;
    static vector< int > PosInCrossing;
    PosInEdges.clear();
    PosInCrossing.clear();
    for (int i = 0; i < sizeAsInt(Edges); ++i) {
      int j = 0;
      while (j < 4 and PD[4 * Crossing + j] != Edges[i]) {
        ++j;
      }
      
      if (j != 4) {
        PosInEdges.push_back(i);
        PosInCrossing.push_back(j);
      }
    }
    
    int Connections = sizeAsInt(PosInEdges);
    if (PosInEdges.empty() or PosInEdges.back() - PosInEdges.front() > Connections) {
      continue;  // crossing attaches nowhere or in disjoint intervals
    }
    for (int i = 0; i + 1 < Connections; ++i) {
      if ((4 - PosInCrossing[i + 1] + PosInCrossing[i]) % 4 != 1) {
        goto CrossingLoopEnd;  // crossing does not attach counterclockwise
      }
    }
    
    if (Connections == MaxConnections) {
      MaxConnectedCrossings.push_back(Crossing);
    }
    else if (Connections > MaxConnections) {
      MaxConnections = Connections;
      MaxConnectedCrossings.clear();
      MaxConnectedCrossings.push_back(Crossing);
    }
    CrossingLoopEnd: ;
  }
  
  return make_pair(MaxConnections, MaxConnectedCrossings);
}

/* Auxiliary function for PlanarDiagram::GetSmallGirthMorseCode.
 * 
 * Extend the Morse list based on the new crossing and its connectivity.
 * Update the bottom edges of the upper knot diagram and the cost of the
 * current Morse list.
 * 
 * Depending on the connectivity, we may need to add local extrema. This is
 * detailed below.
 */
inline void ExtendMorseList(
  vector< int >& MorseList,
  vector< int >& Edges,
  long long& Cost,
  const vector< int >& PD,
  const int NextCrossing,
  const int Connectivity,
  const int nAdded
) {
  /* Find the position of the first bottom edge that connects to the next
   * crossing, as well as the position of the corresponding edge in the part
   * of the planar diagram with the next crossing's edges.
   */
  const vector< int > CrossingEdges(PD.begin() + 4 * NextCrossing, PD.begin() + 4 * (NextCrossing + 1));
  int FirstPos = 0;
  while (
    FirstPos < sizeAsInt(Edges)
    and Edges[FirstPos] != CrossingEdges[0]
    and Edges[FirstPos] != CrossingEdges[1]
    and Edges[FirstPos] != CrossingEdges[2]
    and Edges[FirstPos] != CrossingEdges[3]
  ) {
    ++FirstPos;
  }
  int CrossingFirstPos = 0;
  while (CrossingEdges[CrossingFirstPos] != Edges[FirstPos]) {
    ++CrossingFirstPos;
  }
  
  /* Compute cost of adding the crossing, in particular the cost of adding
   * local minima.
   */
  long long AddedCost = nAdded * sizeAsInt(Edges) * sizeAsInt(Edges);
  if (Connectivity == 3) {
    Cost += (FirstPos + 1) * AddedCost;
  }
  else if (Connectivity == 4) {
    Cost += (2 * FirstPos + 1) * AddedCost;
  }
  
  // Add Morse events according to the new crossing
  /* If Connectivity == 1, add a maximum and a crossing as follows:
   * 
   *     | |   _   |
   *     |  \ / \  |
   *     |   X   | |
   *     |  | |  | |
   * 
   *        L    R   Positions
   * 
   */
  if (Connectivity == 1) {
    int LeftPos = (CrossingFirstPos + 1) % 4;
    int RightPos = (CrossingFirstPos + 3) % 4;
    if ((sizeAsInt(PD) / 2 + CrossingEdges[RightPos] - CrossingEdges[LeftPos]) % (sizeAsInt(PD) / 2) == 1) {
      MorseList.push_back(1000);  // maximum oriented left to right
    }
    else {
      MorseList.push_back(1001);  // maximum oriented right to left
    }
    MorseList.push_back(FirstPos + 2);  // position of maximum
    if (CrossingFirstPos % 2 == 0) {
      MorseList.push_back(FirstPos + 1);  // positive crossing
    }
    else {
      MorseList.push_back(-FirstPos - 1);  // negative crossing
    }
  }
  
  // If Connectivity == 2, nothing complicated happens.
  else if (Connectivity == 2 and CrossingFirstPos % 2 == 0) {
    MorseList.push_back(FirstPos + 1);
  }
  else if (Connectivity == 2) {
    MorseList.push_back(-FirstPos - 1);
  }
  
  /* If Connectivity >= 3, then add a crossing and a minimum, and move
   * the minimum to the left:
   * 
   *     | |  | | | |
   *     | | /  _X  |
   *     |  X_ /  | |
   *     | /  X   | |
   *      X_ / |  | |
   *     /  X  |  | |
   *     \_/ | |  | |
   * 
   * If Connectivity == 4, then we need to do add and move two minima.
   */
  else if (Connectivity >= 3) {
    if (CrossingFirstPos % 2 == 0) {
      MorseList.push_back(-FirstPos - 2);  // negative crossing
    }
    else {
      MorseList.push_back(FirstPos + 2);  // positive crossing
    }
    for (int i = FirstPos; i > 0; --i) {
      MorseList.push_back(i);
      MorseList.push_back(i + 1);
    }
    MorseList.push_back(-1000);  // minimum
    
    // repeat last part if Connectivity == 4
    if (Connectivity == 4) {
      for (int i = FirstPos; i > 0; --i) {
        MorseList.push_back(i);
        MorseList.push_back(i + 1);
      }
      MorseList.push_back(-1000);
    }
  }
  
  // Update edge list
  for (int i = 0; i < Connectivity; ++i) {
    Edges.erase(Edges.begin() + FirstPos);  // erase CrossingEdges[(4 - i + CrossingFirstPos) % 4]
  }
  for (int i = Connectivity; i < 4; ++i) {
    Edges.insert(Edges.begin() + FirstPos, CrossingEdges[(4 - i + CrossingFirstPos) % 4]);
  }
}

/* Search for a small-girth Morse list.
 * 
 * We construct an upper knot diagram, adding crossings one by one. We choose a
 * crossing as the first crossing, and then we successively choose crossings
 * that are strongly connected to the upper knot diagram. In certain cases, we
 * will need to add local extrema.
 * 
 * The number of possible Morse lists may be exponentially large, so we limit
 * to a maximal number of tries and make some choices randomly.
 * 
 * To evaluate the size of the Morse list, we also use a cost function. This
 * cost increases when we need to add local minima, and is greater the more
 * crossings we have added.
 */
MorseCode PlanarDiagram::GetSmallGirthMorseCode(int MaxNumberOfTries) const {
  const vector< int > PD = ListOfTuples;  // local copy of ListOfTuples
  int nCrossings = sizeAsInt(PD) / 4;
  vector< int > SmallestMorseList;
  int MinGirth = 100;  // girth can't go over integer bit size, anyways
  long long MinCost = 1000000000;  // = 1e9
  MaxNumberOfTries = min(100 + nCrossings * nCrossings, MaxNumberOfTries);
  
  for (int Attempt = 0; Attempt < MaxNumberOfTries; ++Attempt) {
    /* INITIALIZATION STEP: select first crossing and orientation, and set
     * the initial values of variables: the edges at the bottom of the upper
     * knot diagram, the Morse list, the girth and cost of the Morse list.
     */
    auto FirstChoice = div(Attempt % sizeAsInt(PD), 4);
    int FirstCrossing = FirstChoice.quot;
    int Shift = FirstChoice.rem;
    vector< bool > Added(nCrossings, false);  // crossings added to the list
    Added[FirstCrossing] = true;
    vector< int > Edges = {
      PD[4 * FirstCrossing + Shift],
      PD[4 * FirstCrossing + (Shift + 1) % 4],
      PD[4 * FirstCrossing + (Shift + 2) % 4],
      PD[4 * FirstCrossing + (Shift + 3) % 4]
    };
    
    for (int x : Edges) { cout << (x - 1) << " "; } cout << endl;
    
    /* Initialize the Morse list with two maxima and the first crossing. We
     * also need the orientation of the maxima. The orientation of the knot
     * follows the edge numbering of the planar diagram, in increasing order:
     * 
     *             _
     *            / \   _
     *           |   \ / \
     *           |    X   |
     * 
     *     Edges 0   1 2  3
     *
     */
    vector< int > MorseList = {1001, 1, 1001, 3, -2};
    if ((2 * nCrossings + Edges[2] - Edges[0]) % (2 * nCrossings) == 1) {
      MorseList[0] = 1000;  // maximum oriented left to right
    }
    if ((2 * nCrossings + Edges[3] - Edges[1]) % (2 * nCrossings) == 1) {
      MorseList[2] = 1000;
    }
    if (Shift % 2 == 0) {
      MorseList[4] = 2;
    }
    int Girth = 4;
    long long Cost = 0;
    
    /* RECURSIVE STEP: iteratively add a maximally connected crossing */
    for (int nAdded = 1; nAdded < nCrossings; ++nAdded) {
      /* Get maximally connected crossings */
      auto Result = GetMaxConnections(PD, Added, Edges);
      auto MaxConnections = Result.first;
      auto MaxConnectedCrossings = Result.second;
      if (MaxConnectedCrossings.empty()) {
        goto AttemptEnd;  // give up on this attempt
      }
      
      /* Randomly select next crossing (almost uniformly) and update
       * MorseList, edges, cost, and girth accordingly.
       */
      int NextCrossing = MaxConnectedCrossings[rand() % sizeAsInt(MaxConnectedCrossings)];
      Added[NextCrossing] = true;
      ExtendMorseList(MorseList, Edges, Cost, PD, NextCrossing, MaxConnections, nAdded);
      if (sizeAsInt(Edges) > Girth) {
        Girth = sizeAsInt(Edges);
      }
      if (Girth > MinGirth) {
        goto AttemptEnd;  // give up on this attempt
      }
    }
    
    // Update smallest Morse list in the case of a successful attempt
    if (Girth < MinGirth or (Girth == MinGirth and Cost < MinCost)) {
      MinGirth = Girth;
      MinCost = Cost;
      SmallestMorseList = MorseList;
    }
    
    AttemptEnd: ;
  }
  
  // Catch the case where all attempts fail
  if (SmallestMorseList.empty()) {
    MinGirth = -1;
  }
  
  return MorseCode(SmallestMorseList, MinGirth);
}