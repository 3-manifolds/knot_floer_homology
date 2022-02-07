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

inline pair< int, vector< int > > get_max_connections(
  const vector< int >& pd,
  const vector< bool >& added,
  const vector< int >& edges
) {
  vector< int > max_connected_crossings;
  int max_connections = 1;
  
  //cout << "[Da] current edges " << edges << endl;
  // Find maximally connected crossings
  for (int crossing = 0; crossing < (sizeAsInt(pd) / 4); ++crossing) {
    if (added[crossing]) {
      continue;
    }
    
    // static variables for memory optimization
    static vector< int > pos_in_edges;
    static vector< int > pos_in_crossing;
    pos_in_edges.clear();
    pos_in_crossing.clear();
    
    for (int i = 0; i < sizeAsInt(edges); ++i) {
      int j = 0;
      while (j < 4 and pd[4 * crossing + j] != edges[i]) {
        ++j;
      }
      
      if (j != 4) {
        pos_in_edges.push_back(i);
        pos_in_crossing.push_back(j);
      }
    }
    int connections = sizeAsInt(pos_in_edges);
    
    if (connections == 0) continue;
    else if (pos_in_edges.back() - pos_in_edges.front() > connections) {
      continue;  // crossing attaches to temp in disjoint intervals
    }
    
    //std::cout << "pos in temp " << pos_in_edges << " pos in crossing " << pos_in_crossing << std::endl;
    
    for (int i = 0; i + 1 < connections; ++i) {
      if ((4 - pos_in_crossing[i + 1] + pos_in_crossing[i]) % 4 != 1) {
        // crossing does not attach counter-clockwise
        //std::cout << "[Da] data " << pos_in_crossing[i] << " " << pos_in_crossing[i + 1] << std::endl;
        //std::cout << "[Da] not counterclockwise " << crossing_edges << std::endl;
        goto crossing_loop_end;  // exit loop + jump using goto
      }
    }
    
    if (connections == max_connections) {
      max_connected_crossings.push_back(crossing);
    }
    else if (connections > max_connections) {
      max_connections = connections;
      max_connected_crossings.clear();
      max_connected_crossings.push_back(crossing);
    }
    crossing_loop_end: ;
  }
  
  return make_pair(max_connections, max_connected_crossings);
}


inline void extend_Morse_list(
  vector< int >& morse_list,
  vector< int >& edges,
  long long& cost,
  const vector< int >& pd,
  const int next_crossing,
  const int max_connections,
  const int n_added
) {
  const vector< int > crossing_edges(pd.begin() + 4 * next_crossing, pd.begin() + 4 * (next_crossing + 1));
  //cout << "[Da] next crossing " << next_crossing << ": " << crossing_edges << endl;
  int first_pos = 0;
  while (
    first_pos < sizeAsInt(edges)
    and edges[first_pos] != crossing_edges[0]
    and edges[first_pos] != crossing_edges[1]
    and edges[first_pos] != crossing_edges[2]
    and edges[first_pos] != crossing_edges[3]
  ) {
    ++first_pos;
  }
  int crossing_first_pos = 0;
  while (crossing_edges[crossing_first_pos] != edges[first_pos]) {
    ++crossing_first_pos;
  }
  
  // Compute cost of adding Morse event
  long long added_cost = n_added * sizeAsInt(edges) * sizeAsInt(edges);
  if (max_connections == 3) {
    cost += (first_pos + 1) * added_cost;
  }
  else if (max_connections == 4) {
    cost += (2 * first_pos + 1) * added_cost;
  }
  
  // Add Morse events according to the new crossing
  /* If max_connections == 1, add a maximum and a crossing as follows:
   * 
   *     | |   _   |
   *     |  \ / \  |
   *     |   X   | |
   *     |  | |  | |
   * 
   */
  if (max_connections == 1) {
    int left_pos = (crossing_first_pos + 1) % 4;
    int right_pos = (crossing_first_pos + 3) % 4;
    if ((crossing_edges[right_pos] - crossing_edges[left_pos]) % (sizeAsInt(pd) / 2) == 1) {
      morse_list.push_back(1000);  // maximum oriented left to right
    }
    else {
      morse_list.push_back(1001);  // maximum oriented right to left
    }
    morse_list.push_back(first_pos + 2);  // position of maximum
    if (crossing_first_pos % 2 == 0) {
      morse_list.push_back(first_pos + 1);  // positive crossing
    }
    else {
      morse_list.push_back(-first_pos - 1);  // negative crossing
    }
  }
  
  // If max_connections == 2, nothing complicated happens.
  else if (max_connections == 2 and crossing_first_pos % 2 == 0) {
    morse_list.push_back(first_pos + 1);
  }
  else if (max_connections == 2) {
    morse_list.push_back(-first_pos - 1);
  }
  
  /* If max_connections >= 3, then add a crossing and a minimum, and move
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
   * If max_connections == 4, then we need to do this process twice.
   */
  else if (max_connections >= 3) {
    if (crossing_first_pos % 2 == 0) {
      morse_list.push_back(-first_pos - 2);  // negative crossing
    }
    else {
      morse_list.push_back(first_pos + 2);  // positive crossing
    }
    for (int i = first_pos; i > 0; --i) {
      morse_list.push_back(i);
      morse_list.push_back(i + 1);
    }
    morse_list.push_back(-1000);  // minimum
    
    // repeat if max_connections == 4
    if (max_connections == 4) {
      for (int i = first_pos; i > 0; --i) {
        morse_list.push_back(i);
        morse_list.push_back(i + 1);
      }
      morse_list.push_back(-1000);
    }
  }
  
    
  // Update edge list
  for (int i = 0; i < max_connections; ++i) {
    edges.erase(edges.begin() + first_pos);  // erase crossing_edges[(4 - i crossing_first_pos) % 4]
  }
  for (int i = max_connections; i < 4; ++i) {
    edges.insert(edges.begin() + first_pos, crossing_edges[(4 - i + crossing_first_pos) % 4]);
  }
}

/* Search for a small-girth Morse list. Adapted from ComputeHFKv2/Diagrams.cpp
 * 
 * We construct an upper knot diagram, adding crossings one by one. We choose a
 * crossing as the first crossing, and then we successively choose crossings
 * that are strongly connected to the upper knot diagram. In certain cases, we
 * will need to add local extrema.
 * 
 * The number of possible Morse lists may be exponentially large, so we limit
 * to a maximal number of tries and test choices randomly.
 * 
 * To evaluate the size of the Morse list, we also use a cost function. This
 * cost increases when we need to add local minima, and is greater the more
 * crossings we have added.
 */
MorseCode PlanarDiagram::GetSmallGirthMorseCode(int max_number_of_tries) const {
  const vector< int > pd = ListOfTuples;  // local copy of ListOfTuples
  int n_crossings = sizeAsInt(pd) / 4;
  vector< int > smallest_morse_list;
  int min_girth = 100;  // girth can't go over integer bit size, anyways
  long long min_cost = 1000000000;  // = 1e9
  max_number_of_tries = min(100 + n_crossings * n_crossings, max_number_of_tries);
  
  for (int attempt = 0; attempt < max_number_of_tries; ++attempt) {
    /* STEP 1: Initialize the Morse list with the first crossing */
    
    int girth = 4;  // girth of current Morse list
    long long cost = 0;
    
    // Choose first crossing and orientation, in order
    auto first_choice = div(attempt % sizeAsInt(pd), 4);
    int first_crossing = first_choice.quot;
    int shift = first_choice.rem;
    
    // Edges at the bottom of the upper knot diagram
    vector< int > edges = {
      pd[4 * first_crossing + shift],
      pd[4 * first_crossing + (shift + 1) % 4],
      pd[4 * first_crossing + (shift + 2) % 4],
      pd[4 * first_crossing + (shift + 3) % 4]
    };
    
    /* Add two maxima and the first crossing to Morse list. We also calculate
     * the orientation of the maxima. The orientation of the knot follows the
     * edge numbering in increasing order.
     * There are twice as many edges as crossings.
     */
    vector< int > morse_list = {1001, 1, 1001, 3, -2};
    if ((edges[2] - edges[0]) % (2 * n_crossings) == 1) {
      morse_list[0] = 1000;
    }
    if ((edges[3] - edges[1]) % (2 * n_crossings) == 1) {
      morse_list[2] = 1000;
    }
    if (shift % 2 == 0) {
      morse_list[4] = 2;
    }
    
    vector< bool > added(n_crossings, false);  // added crossings
    added[first_crossing] = true;
    //cout << "\n[Da] first crossing " << first_crossing << ": " << edges << endl;
    
    /* STEP 2: Iteratively add a maximally connected crossing */
    for (int n_added = 1; n_added < n_crossings; ++n_added) {
      // Get maximally connected crossings
      auto result = get_max_connections(pd, added, edges);
      auto max_connections = result.first;
      auto max_connected_crossings = result.second;
      
      if (max_connected_crossings.empty()) {
        //cout << "[Da] no good crossings found" << endl;
        goto attempt_end;  // give up on this attempt
      }
      
      // Nearly uniform random choice of next crossing
      int next_crossing = max_connected_crossings[rand() % sizeAsInt(max_connected_crossings)];
      added[next_crossing] = true;
      
      // Update morse_list, edges, and cost
      extend_Morse_list(morse_list, edges, cost, pd, next_crossing, max_connections, n_added);
      
      if (sizeAsInt(edges) > girth) {
        girth = sizeAsInt(edges);
      }
      if (girth > min_girth) {
        //cout << "[Da] girth too big" << endl;
        goto attempt_end;  // give up on this attempt
      }
    }
    
    // Update smallest Morse list in the case of a successful attempt
    if (girth < min_girth or (girth == min_girth and cost < min_cost)) {
      min_girth = girth;
      min_cost = cost;
      smallest_morse_list = morse_list;
    }
    
    attempt_end: ;
  }
  
  // Catch the case where all attempts fail
  if (smallest_morse_list.empty()) {
    min_girth = -1;
  }
  
  return MorseCode(smallest_morse_list, min_girth);
}