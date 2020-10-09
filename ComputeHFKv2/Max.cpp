#include "Alg.h"
#include<vector>
#include<string>
#include<array>
#include<map>
#include<algorithm>
#include<iostream>
using namespace std;


/////////////////////////
//  MAXIMUM BIMODULE
////////////////////////

// X, Y, Z are the bimodule generators.
// Maximum is added between walls Position -1 and Position. 
// Position is between 1 and 2*Bridge+1, where Bridge is computed before the maximum is added
// After maximum is added, Position and Position +1 are matched.

idem ExtendX(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 % pow; 
    int y= (I1/(2*pow))*8*pow;
    return x + 3*pow + y;
}

idem ExtendY(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 %pow; 
    int y= (I1/(2*pow))*8*pow;
    return x + 6*pow + y;
}

idem ExtendZ(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 %pow ; 
    int y= (I1/(2*pow))*8*pow;
    return x + 2*pow + y;
}

void MaxM2Action(Arrow arrow1, int Position)//Bridge is increased before
{   idem I1=GeneratorList[arrow1.StartingGen].Idem;
    idem I2=GeneratorList[arrow1.EndingGen].Idem;
    bool XY= !!(I1 &(1<<(Position-1)));
    bool Z= !XY;

    monomial m=MonomialStore[arrow1.MonomialIndex];
    monomial NewM=MonomialOne; 
    for (int i=0; i<Position-1; i++) NewM[i]=m[i]; 
    for (int i=Position+1; i< 2*Bridge; i++) NewM[i]=m[i-2];
    NewM[Position-1]=0; NewM[Position]=0;
    int r=MonomialLookUp(NewM);

    int a=arrow1.StartingGen; int b=arrow1.EndingGen;
    int S =(LeftCoord(I1,Position)%2); //mod2 grading on bimodule generators: S on X; S+1 on Y, S on Z
    int Sign; if(S%2==0) Sign =1; else Sign =-1;
    int f1; int f2;          //Measures how \phi^(I1,I2)(1) crosses walls Position -1 and Position respectively
    if (Position ==1) f1=0; 
    else  f1=LeftRight(I1, I2, Position-1);
    f2=LeftRight(I1, I2, Position);
  
    if(f1==1 && f2==1) // Incoming algebra element's local weight=R2R1, Ex: R_{Position}R_{Position-1}
       {idem I3=ExtendY(I1, Position); idem I4=ExtendX(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen=3*a+1; Arr.EndingGen=3*b;
        Arr.Coeff=(arrow1.Coeff * (Modulus-Sign))%Modulus ;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1== -1 && f2== -1) //L1L2
       {idem I3=ExtendX(I1, Position); 
        idem I4=ExtendY(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a; Arr.EndingGen=3*b+1;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==1 && f2==0) //R1
       {idem I3=ExtendZ(I1, Position); 
        idem I4=ExtendX(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a+2; Arr.EndingGen=3*b;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==-1 && f2==0) //L1
       {idem I3=ExtendX(I1, Position); 
        idem I4=ExtendZ(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a; Arr.EndingGen=3*b+2;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==0 && f2==1) //R2
       {idem I3=ExtendY(I1, Position); 
        idem I4=ExtendZ(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a+1; Arr.EndingGen=3*b+2;
        Arr.Coeff=(arrow1.Coeff * (Modulus-Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==0 && f2== -1) //L2
       {idem I3=ExtendZ(I1, Position); 
        idem I4=ExtendY(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a+2; Arr.EndingGen=3*b+1;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==0 && f2==0 && Z) // Z to Z
       {idem I3=ExtendZ(I1, Position); idem I4=ExtendZ(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a+2; Arr.EndingGen=3*b+2;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==0 && f2==0 && XY) // X to X
       {idem I3=ExtendX(I1, Position); idem I4=ExtendX(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a; Arr.EndingGen=3*b;
        Arr.Coeff=(arrow1.Coeff * (Modulus+Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }

    if(f1==0 && f2==0 && XY) // Y to Y
       {idem I3=ExtendY(I1, Position); idem I4=ExtendY(I2,Position);
        Arrow Arr=arrow1; 
        Arr.MonomialIndex=r;
        Arr.StartingGen= 3*a+1; Arr.EndingGen=3*b+1;
        Arr.Coeff=(arrow1.Coeff * (Modulus-Sign))%Modulus;
        if(NonZero(I3,I4,r)) NewArrowList.push_back(Arr); }
}

//Tensoring with the maximum bimodule. 
//Precondition. GeneratorList[i].Name=i;
//Orientation=0 , means maximum is oriented left to right
void AfterMax(int Position, int Orientation) 
{   vector<int> NewMatch(2*Bridge +2); //Update MatchingList, UpwardList;
    vector<int> NewUp(2*Bridge+2);
  
    for(int i=0; i<Position -1;i++) 
      {int a=MatchingList[i]; 
       if(a < Position) NewMatch[i]=a; else NewMatch[i]=a+2;
       NewUp[i]=UpwardList[i];}
  
    for(int i=Position+1; i<2*Bridge+2; i++)
      {int a=MatchingList[i-2]; 
       if(a < Position) NewMatch[i]=a; else NewMatch[i]=a+2;
       NewUp[i]=UpwardList[i-2];}
   
    NewMatch[Position-1]=Position+1; NewMatch[Position]=Position;
    if(Orientation ==0) {NewUp[Position-1]=1; NewUp[Position]=0;}
    if(Orientation ==1) {NewUp[Position-1]=0; NewUp[Position]=1;}
    MatchingList=NewMatch;  UpwardList=NewUp;
  
    Bridge++;

    NewArrowList.clear();
    NewGeneratorList.clear();
     
    for(Gen G: GeneratorList)   //Names X=3*a, Y=3*a+1, Z=3*a+2
      {idem I1=G.Idem; int a=G.Name;
       bool XY= !!(I1 &(1<<(Position-1)));
       bool Z= !XY ;
       
       if(Z) {G.Idem=ExtendZ(I1, Position); G.Name=3*a+2; 
              NewGeneratorList.push_back(G); }
     
       if(XY) {G.Idem=ExtendX(I1,Position); G.Name=3*a;   NewGeneratorList.push_back(G);    
       	       G.Idem=ExtendY(I1,Position); G.Name=3*a+1; NewGeneratorList.push_back(G);
               
               Arrow Arr1; Arrow Arr2; // M1 actions
               Arr1.StartingGen=3*a;   Arr1.EndingGen=3*a+1; Arr1.MonomialIndex=0; Arr1.Coeff=1; //R_{i+1}\cdot R_i
               Arr2.StartingGen=3*a+1; Arr2.EndingGen=3*a;   Arr2.MonomialIndex=0; Arr2.Coeff=1; //L_i\cdot L_{i+1}
	       NewArrowList.push_back(Arr1); NewArrowList.push_back(Arr2);}}
    
    for(Arrow A : ArrowList)  MaxM2Action(A, Position);

    GeneratorList =NewGeneratorList;
    NewGeneratorList.clear();
    ArrowList=NewArrowList;
    NewArrowList.clear();
     
    RemoveMod(ArrowList);
    ReName();
}      








