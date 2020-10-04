#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include<vector>
#include<string>
#include<iostream>

using namespace std;


class MorseCode {
  private: 
  vector<int> MorseList;
  int Girth;

  public:  
  MorseCode(vector<int> X, int G) {
    MorseList=X;
    Girth=G;
  };
  
  vector<int> GetMorseList() {
    return MorseList;
  };
  
  int GetGirth() {
    return Girth;
  };

  void Print(ostream & os);
};


class PlanarDiagram {
  private:
  vector<int> ListOfTuples;

  public:

  vector<int> GetListOfTuples() {
    return ListOfTuples;
  };

  PlanarDiagram(string S);
  
  bool NotValid();
  
  bool Alternating();

  MorseCode GetSmallGirthMorseCode(int MaxNumberOfTries = 10000);

  bool R1Reducible();
  
  void Print(ostream & os);
};

#endif

