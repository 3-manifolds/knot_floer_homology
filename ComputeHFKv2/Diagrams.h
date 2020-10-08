#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include<vector>
#include<string>
#include<iostream>

class MorseCode {
  private: 
  std::vector<int> MorseList;
  int Girth;

  public:  
  MorseCode(const std::vector<int> X, int G)
      : MorseList(X), Girth(G)
  {
  };
  
  const std::vector<int> &GetMorseList() const {
    return MorseList;
  };
  
  int GetGirth() const {
    return Girth;
  };

  void Print(std::ostream & os) const;
};


class PlanarDiagram {
  private:
  std::vector<int> ListOfTuples;

  public:

  const std::vector<int> &GetListOfTuples() const {
    return ListOfTuples;
  };

  PlanarDiagram(const std::string &S);
  
  bool NotValid() const;
  
  bool Alternating() const;

  MorseCode GetSmallGirthMorseCode(int MaxNumberOfTries = 10000) const;

  bool R1Reducible() const;
  
  void Print(std::ostream & os) const;
};

#endif

