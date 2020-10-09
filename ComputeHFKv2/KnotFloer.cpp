#include "Alg.h"
#include "Diagrams.h"
#include<vector>
#include<map>
#include<unordered_map>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

KnotFloerComplex  ComputingKnotFloer(MorseCode Morse, int Prime, bool verbose){    
    MatchingList.clear(); UpwardList.clear(); 
    GeneratorList.clear(); NewGeneratorList.clear(); 
    ArrowList.clear(); NewArrowList.clear();
    Bridge=1;
    MatchingList.push_back(2); MatchingList.push_back(1);
    Gen G; G.Name=0; 
    G.Idem=2; //idempotent={1}, 
    G.Maslov=0; G.Alexander=0;
    GeneratorList.push_back(G); //At this point ArrowList is empty. 
    Modulus=Prime; //finished creating the starting D-module
    vector<int> MorseList = Morse.GetMorseList(); 
    
    //orientation of the first maximum:
    if(MorseList[0]==1000) {UpwardList.push_back(1); UpwardList.push_back(0);} 
    if(MorseList[0]==1001) {UpwardList.push_back(0); UpwardList.push_back(1);}
     
    int Steps=sizeAsInt(MorseList);
    if (verbose) {
      cout<<"Computation is with mod "<<Modulus<<" coefficients"<<endl;
      cout<<"Steps to do:"<<endl;
      for(int i=2, B=1; i< Steps -1 ; i++)
	{if(MorseList[i]>999) {cout<<(++B);i++;}
	  else if(MorseList[i]<-999) cout<<(--B);
          else cout<<".";
	  cout.flush();}
      cout<<endl;
      cout<<"Steps in progress:"<<endl; 
    }
    for(int i=2; i< Steps -1 ; i++)
	{
	  if (verbose) {
	    if(MorseList[i]>999) {cout<<Bridge+1;}
	    else if(MorseList[i]<-999) cout<<Bridge-1;
	    else cout<<".";
	    cout.flush();
	  }

         if (MorseList[i]==1000) AfterMax(MorseList[++i], 0); 
	 else if (MorseList[i]==1001) AfterMax(MorseList[++i], 1);
         else if (MorseList[i]==-1000) AfterMin(); 
         else  AfterCrossing(MorseList[i]);
         
         Simplify(); 
	}
    if (verbose) {
      cout<<endl;//(GeneratorList,ArrowList) represents a D-Module over B(2,1);
    }
    ReName(); 
    KnotFloerComplex KFC; 
    KFC.Prime=Prime;
    for(auto G : GeneratorList)
	 {KnotFloerGen G2; G2.Name=G.Name; 
	  G2.Maslov=G.Maslov; G2.Alexander=G.Alexander/2; 
	  KFC.Generators.push_back(G2);}
    for(auto A: ArrowList)
	{ChainArrow A2; 
         A2.StartingGen=A.StartingGen;
         A2.EndingGen=A.EndingGen;
         A2.Coeff=A.Coeff;
         KFC.Differential.push_back(A2); }

    return KFC;
}           

 
