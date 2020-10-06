#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class MorseCode
{
public:  
    MorseCode(const std::vector<int> &morseList,
              const int girth)
      : _morseList(morseList), _girth(girth)
    {
    }
  
    const std::vector<int> &GetMorseList() const {
        return _morseList;
    };
    
    int GetGirth() const {
        return _girth;
    };
    
    void Print(std::ostream & os) const;
    
private: 
    std::vector<int> _morseList;
    int _girth;
};


class PlanarDiagram
{
public:
    PlanarDiagram(const std::string &);
  
    const std::vector<int> &GetListOfTuples() {
        return ListOfTuples;
    };

    bool NotValid();
    
    bool Alternating();
    
    MorseCode GetSmallGirthMorseCode(int MaxNumberOfTries = 10000);
    
    bool R1Reducible();
    
    void Print(ostream & os);
    
private:
    std::vector<int> ListOfTuples;
};

#endif

