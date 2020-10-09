#include "Alg.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
using namespace std;

int Genus(const KnotFloerComplex & KFC)
{   int ans=KFC.Generators[0].Alexander;
    for(KnotFloerGen G : KFC.Generators) if(G.Alexander >ans) ans=G.Alexander;
    return ans;
}

map<pair<int,int>, int> KnotFloerRanks(const KnotFloerComplex & KFC) 
{   map<pair<int, int>, int> Map;
    for(KnotFloerGen G: KFC.Generators)
    Map[make_pair(G.Alexander, G.Maslov)]++; 
    return Map;
}

void ReportKnotFloerRanks(const KnotFloerComplex & KFC, ostream & os)
{   auto Map=KnotFloerRanks(KFC);        
    os<<"Ranks in Alexander, Maslov bigradings :"<<endl;
    for(auto X: Map)
    os<< X.second<<"   (" << X.first.first <<","<< X.first.second<<")" <<endl;
}    

bool Fibered(const KnotFloerComplex & KFC) 
{   int multiplicity=0; int genus=Genus(KFC); 
    for(KnotFloerGen G: KFC.Generators) if(G.Alexander==genus) multiplicity++;
    bool b= (multiplicity==1);
    return b;
}

bool LSpaceKnot(const KnotFloerComplex & KFC)
{   bool b=true;
    pair<int,int> A=make_pair(-10000,-10000);
    auto Map=KnotFloerRanks(KFC);
    for(auto X: Map)
       {pair<int,int> B=X.first;  
        if(X.second >1 || B.first <=A.first || B.second <=A.second) b=false;
        A=B;
       }
    return b;
}

KnotFloerComplex DualComplex(const KnotFloerComplex & KFC)
{   KnotFloerComplex C1; C1.Prime=KFC.Prime;
    for(KnotFloerGen G: KFC.Generators)
	{int a=G.Alexander; int m=G.Maslov;
	 G.Alexander=-a; G.Maslov=-m;
         C1.Generators.push_back(G);}
    for(ChainArrow A: KFC.Differential )
      {int x = A.StartingGen; int y= A.EndingGen; A.StartingGen=y; A.EndingGen=x; //reversing the arrow
         C1.Differential.push_back(A);}
    return C1;
}

int LowestFiltrationWithMaslovZero(const KnotFloerComplex & KFC)
{   int ans=Genus(KFC);
    for(KnotFloerGen G: KFC.Generators) 
         if(G.Alexander <ans && G.Maslov==0) ans=G.Alexander;
    return ans;
}

int Tau(const KnotFloerComplex & KFC)  
{   int x=(int)KFC.Generators.size();
    int Rank1=0; int Rank2=1; int genus=Genus(KFC); 
    int j= LowestFiltrationWithMaslovZero(KFC);       
    while(Rank2==Rank1+1 && j<=genus)
	 {ChainComplex C1;  C1.Prime=KFC.Prime; 
          //Constructing C_1, the Alexander <= j complex. 
	  for(KnotFloerGen G: KFC.Generators)
	    if(G.Alexander <= j) C1.Generators.push_back(G.Name); 
          for(ChainArrow A: KFC.Differential)
            {int a=KFC.Generators[A.StartingGen].Alexander;  
             int b=KFC.Generators[A.EndingGen].Alexander;  
             if(b<= a && a<= j) C1.Differential.push_back(A); }
      	             
          ChainComplex C2=C1; //C_2 will be a mapping cone from C1 to Filtered
          for(int i=0; i< x; i++) C2.Generators.push_back(x+i); //Adding the Filtered Generators
          for(ChainArrow A: KFC.Differential) //Adding the Filtered Differentials; 
	    {int a=KFC.Generators[A.StartingGen].Alexander; 
             int b=KFC.Generators[A.EndingGen].Alexander;   
             A.StartingGen+=x;
             A.EndingGen  +=x;
	     if(a >= b) C2.Differential.push_back(A);}
          for(KnotFloerGen G : KFC.Generators) //Adding the inclusion map  from C_1 to Filtered
	      {ChainArrow A; 
               A.StartingGen=G.Name; 
               A.EndingGen=G.Name+x; 
               if(G.Maslov %2==0) A.Coeff=1; else A.Coeff=KFC.Prime-1;        
	       if(G.Alexander<= j) C2.Differential.push_back(A);}
 	  Rank1=HomologyRank(C1);
          Rank2=HomologyRank(C2); // If Rank2=Rank1+1 f_\ast is trivial on homology
          j++;	
	 }
    return j-1;
}

int Nu(const KnotFloerComplex & KFC){ 
  int x=(int)KFC.Generators.size(); 
    int Rank1=0; int Rank2=1; 
    int j=Tau(KFC); int genus=Genus(KFC);
      
    while(Rank2==Rank1+1 && j<= genus)
      {ChainComplex C1; C1.Prime=KFC.Prime; 

       for(int t=0; t<x; t++) C1.Generators.push_back(t);
       for(ChainArrow  A: KFC.Differential)
	     {int a=KFC.Generators[A.StartingGen].Alexander; 
              int b=KFC.Generators[A.EndingGen].Alexander;
              if ((b<= a && a <=j) || (b>=a && a>=j ))
	         C1.Differential.push_back(A);}
      
       ChainComplex C2=C1;
       for(int i=0; i<x; i++) C2.Generators.push_back(i+x);//Adding the Filtered generators
       for(ChainArrow A: KFC.Differential) 
	     {int a=KFC.Generators[A.StartingGen].Alexander; 
              int b=KFC.Generators[A.EndingGen].Alexander;
              A.StartingGen+=x; 
              A.EndingGen+=x;
              if(a >= b) C2.Differential.push_back(A);}
       for(KnotFloerGen G: KFC.Generators) //Adding the projection map g from C1 to Filtered
	     {int a=G.Alexander;
              ChainArrow A; 
              A.StartingGen=G.Name; 
              A.EndingGen=G.Name+x; 
              if(G.Maslov %2==0) A.Coeff=1; else A.Coeff=KFC.Prime-1;          
	      if(a<= j) C2.Differential.push_back(A);}
	Rank1=HomologyRank(C1);  // Rank2=Rank1+1 if g_\ast is trivial
        Rank2=HomologyRank(C2);
        j++;  
	}
    return j-1;
}

int NuOfMirror(const KnotFloerComplex & KFC){
    return Nu(DualComplex(KFC));
}

int Epsilon(const KnotFloerComplex & KFC){ 
    int nu=Nu(KFC);
    int tau=Tau(KFC);
    int nu2=NuOfMirror(KFC); 
    int ans=0;
    if(nu==tau+1) ans=-1;
    else if(nu==tau && nu2== -tau) ans=0;
    else if(nu2== -tau+1) ans=1;
    return ans;
}
