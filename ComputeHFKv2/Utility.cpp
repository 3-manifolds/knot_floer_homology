#include "Alg.h"
#include<vector>
#include<array>
#include<unordered_map>
#include<algorithm>
#include<iostream>
using namespace std;

///////////////////////
// UTILITY FUNCTIONS
//////////////////////

  //Wall is from 1 to 2*Bridge.
  //computes the number of coordinates left of the wall:
int LeftCoord(idem I1, int Wall) 
{   int a=0; I1=I1&((1<<Wall)-1);
    while(I1) {a++, I1=I1&(I1-1);}
    return a;
}

  //evaluates 1 on R_{Wall}, -1 on L_{Wall}, and 0 on U_{Wall}:
int LeftRight(idem I1, idem I2, int Wall) 
{   return LeftCoord(I1, Wall) - LeftCoord(I2, Wall);
}                                         

int MonomialLookUp(const monomial& X)
{   auto iter=MonomialMap.find(X); 
    if(iter !=MonomialMap.end() ) 
          return (*iter).second;
    else {int z=sizeAsInt(MonomialStore); 
          MonomialMap.insert(make_pair(X,z) ); 
          MonomialStore.push_back(X);
          return z;} 
}

  //checks one of the relations in the ideal in B_0(2*Bridge, Bridge). 
  //True if going from I1 to I2 one of the coordinates moves more that 1 unit 
bool TooFar(idem I1, idem I2) 
{   int a=0; int i=2*Bridge; 
    while(i--)
      {if(a>0 && ((I2 &1)==0)) return true;
       if(a<0 && ((I1 &1)==0)) return true;
	    a=a+(I1 &1)-(I2& 1);I1=I1>>1; I2=I2>>1;}                    
    return false;   
}

  //Checks that the algebra element represented by (I1, I2, m) is non-zero 
  //in B(2*Bridge, Bridge). Precondition: I1 and I2 are not TooFar.
bool NonZero(idem I1, idem I2, int m)
{   monomial P=MonomialStore[m];      
    bool crossed=false;      //Measures if the  wall is crossed. 
    bool b, c, d;          
    bool r=true;             //If r is true we need to check the U_j power in P
    int j=0; int i=2*Bridge;          
    while(i--)
      {I1=I1>>1; I2=I2>>1; d=!!P[j]; j++;
       b=(I1 & 1); //b=true if I1 has a coordinate between walls j and j+1
       c=(I2 & 1);  
       if(r && (!b  || !c) && d) return false; //found the intervall
       else if((crossed  && (b!=c)) || (!crossed  && !b && !c)) r=true; 
       else if(r && !crossed && b && c && d) r=true; //grow the intervall 
       else r=false; //we look for a new intervall to be checked
       crossed= crossed ^ b ^ c;} 
    return true;
}

int Mult(idem I1, idem I2, idem I3, int m1, int m2){   
    int a1=0; int a2=0; int a3=0;                  
    monomial X1=MonomialStore[m1];                
    monomial X2=MonomialStore[m2];                 
    monomial X=MonomialOne; 
    int i=2*Bridge; int j=0;
    while(i--)
      {if(I1 & 1) a1++;
       if(I2 & 1) a2++;
       if(I3 & 1) a3++;
       if((a1>a2 && a2<a3) || (a1<a2 && a2>a3)) 
             X[j]=X1[j]+X2[j]+1;
       else  X[j]=X1[j]+X2[j];
       j++; I1=I1>>1; I2=I2>>1; I3=I3>>1;}
 
    return MonomialLookUp(X);    
}

bool Strict(const Arrow & A, const Arrow & B)
{   return(   (A.StartingGen < B.StartingGen) 
           || (A.StartingGen ==B.StartingGen &&  A.EndingGen < B.EndingGen) 
           || (A.StartingGen ==B.StartingGen  &&  A.EndingGen == B.EndingGen 
	       && A.MonomialIndex < B.MonomialIndex) ); 
} 

//Coeff could be different
bool Equal(const Arrow & A, const Arrow & B)
{   return ((A.StartingGen   == B.StartingGen)  &&  
            (A.EndingGen     == B.EndingGen) && 
	    (A.MonomialIndex == B.MonomialIndex) );
}

void  RemoveMod(vector<Arrow> & List)
{   if(List.size()==0) return;
    sort(List.begin(), List.end(), Strict);
    Arrow A=List[0]; A.Coeff=0; int Write=0;
    for(Arrow B: List) 
	{if(Equal(A,B)) A.Coeff=(A.Coeff+B.Coeff)%Modulus; 
         else if(A.Coeff==0) A=B;
         else {List[Write]=A; A=B; Write++;}
	}
    if(A.Coeff !=0) {List[Write]=A; Write++;} //adding the last element
    List.erase(List.begin()+Write, List.end());
}

void ReName() // This function arranges that for all i GeneratorList[i].Name=i
{   int Max=0;  int i=0;
    for(Gen G: GeneratorList) if(G.Name >Max) Max=G.Name;
    vector<int> Dictionary(Max+1,-1);
    for(Gen & G: GeneratorList) 
        {Dictionary[G.Name]=i; 
         G.Name=i;
         i++;}
    for(Arrow & A: ArrowList)
	 {A.StartingGen=Dictionary[A.StartingGen]; 
	  A.EndingGen  =Dictionary[A.EndingGen];}     
}





