#include "Alg.h"
#include<vector>
#include<string>
#include<array>
#include<map>
#include<algorithm>
#include<iostream>
using namespace std;


//////////////////////////
//  Minimum Bimodule
////////////////////////
 
  //Join is needed for the higher actions. 
  //It glues together two sets of arrows.
  //Used in the AfterMin function
vector<Arrow> Join(vector<Arrow> &  X, vector<Arrow> & Y)
{   vector<Arrow> Z; 
    sort(X.begin(),X.end(), [](Arrow a, Arrow b){return a.EndingGen < b.EndingGen;} ); 
    sort(Y.begin(),Y.end(), [](Arrow a, Arrow b){return a.StartingGen < b.StartingGen;} );
    int y1=0, y2=0, x1=0, x2=0;
    for(int Middle=0; Middle<sizeAsInt(GeneratorList); Middle++)
      {while(x2<sizeAsInt(X) && X[x2].EndingGen==Middle) x2++;
	while(y2<sizeAsInt(Y) && Y[y2].StartingGen==Middle) y2++;
        for(int i=x1; i<x2; i++)
          for(int j=y1; j<y2; j++)
		    {Arrow C; Arrow A=X[i]; Arrow B=Y[j];
                     int From=A.StartingGen; int To=B.EndingGen; 
		     idem I1=GeneratorList[From].Idem; 
                     idem I2=GeneratorList[Middle].Idem; 
                     idem I3=GeneratorList[To].Idem;
                     int m1=A.MonomialIndex; int m2=B.MonomialIndex;
	             C.StartingGen=From;
                     C.EndingGen=To;
                     C.MonomialIndex =Mult(I1,I2,I3,m1,m2);
                     C.Coeff      = (A.Coeff * B.Coeff)% Modulus;
		     if( ! TooFar(I1,I3) ) Z.push_back(C); }
	x1=x2; y1=y2;
       }
    RemoveMod(Z);
    return Z; 
}

void AfterMin()
{   int aa= MatchingList[0], bb= MatchingList[1];
    vector<Arrow> MoveU1;
    vector<Arrow> MoveU2;
    vector<Arrow> MoveL2;
    vector<Arrow> MoveR2;
    vector<Arrow> PartialSequence;
    vector<Arrow> FinishedSequence;
   
    NewArrowList.clear();
    for(Arrow A: ArrowList)
      {monomial M=MonomialStore[A.MonomialIndex];
       int a1=M[0], b1=M[1]; 
       int x1=GeneratorList[A.StartingGen].Idem, x2=GeneratorList[A.EndingGen].Idem; 
       if (x1 %8 ==4 && x2 %4==2)                      MoveL2.push_back(A); //local weight: L2\times U2^b1, b1 >= 0
       if (x1 %4 ==2 && x2 %8==4)                      MoveR2.push_back(A); //local weight: R2\times U2^b1  b1 >= 0
       if (x1 %4 ==2 && x2 %4==2 && a1>0  &&  b1==0) MoveU1.push_back(A); //local weight: U1^a1           a1 >  0
       if (x1 %4 ==2 && x2 %4==2 && a1==0 &&  b1>0 ) MoveU2.push_back(A); //local weight: U2^b1           b1 >  0
       if (x1 %8 ==4 && x2 %8==4 && a1==0)                                 //local weight: U2^b1,          b1 >= 0
	  {M[aa-1] =M[aa-1]+b1; // b1>0 is a "curved" action, (U_2,C_1), (U_2 ^2,C_1,C_1), (U_2^ b1, C_1,...,C_1)
	   idem I1 =(x1/4)-1; 
	   idem I2 =(x2/4)-1; 
           for(int i=0; i< 2*Bridge -2; i++) M[i]=M[i+2]; 
           M[2*Bridge -2]=0; M[2*Bridge -1]=0;
           int m=MonomialLookUp(M);
           A.MonomialIndex=m;
           Bridge--; //need to decrease Bridge for NonZero to work 
           if(!TooFar(I1,I2) && NonZero(I1,I2,m))   
	         NewArrowList.push_back(A);//Coeff is unchanged
           Bridge++;} //restoring Bridge
      }
    ArrowList.clear();//Removing excess
    vector<Arrow>().swap(ArrowList);  

    if(MoveL2.size() >0 && MoveU1.size() >0) PartialSequence=Join(MoveL2, MoveU1); 
    for(int Level=1; Level<200; Level++) //Level 1: (L2, U1, R2), Level 2, (L2, U1, U2, U1, R2). 
         {if (PartialSequence.size()==0) break;
          if (MoveR2.size()==0)          break;
          FinishedSequence=Join(PartialSequence,MoveR2); //we finish the sequence and add the resulting Arrows to 
          for(Arrow A: FinishedSequence)
             {monomial M=MonomialStore[A.MonomialIndex];
	      int a1=M[0]; int b1=M[1];
              M[aa-1] =M[aa-1]+b1-Level; M[bb-1]=M[bb-1]+a1-Level;  
              //b1>Level or a1>Level are curved actions. Ex: (L2, U1^2, C2, R2U2^3, C1, C1, C1) is Level 1
	      // we get this action from Join(Join(MoveL2, MoveU1), MoveR2)
              for(int i=0; i< 2*Bridge -2; i++) M[i]=M[i+2]; 
              M[2*Bridge -2]=0; M[2*Bridge -1]=0;
              int m=MonomialLookUp(M);
              A.MonomialIndex=m;
              int x1=GeneratorList[A.StartingGen].Idem; int x2=GeneratorList[A.EndingGen].Idem; 
              idem I1=(x1/4) -1; idem I2=(x2/4) -1;      
              if(Level %2 ==1) A.Coeff= (Modulus-A.Coeff)%Modulus; 
              //Ex: (L2, U1, R2) acts with a minus sign, (L2, U1, U2, U1, R2) acts with a plus sign
              Bridge--; 
              if(!TooFar(I1,I2) && NonZero(I1,I2,m) ) NewArrowList.push_back(A);
              Bridge++;
	     }
          if(MoveU2.size()==0) break;
          PartialSequence=Join(PartialSequence,MoveU2);
          if(PartialSequence.size() ==0) break;
          PartialSequence=Join(PartialSequence,MoveU1);
	 }
    MoveU1.clear(); MoveU2.clear(); MoveL2.clear(); MoveR2.clear(); PartialSequence.clear(); FinishedSequence.clear();
    vector<Arrow>().swap(MoveU1); vector<Arrow>().swap(MoveU2);  
    vector<Arrow>().swap(MoveL2);  vector<Arrow>().swap(MoveR2);    
    vector<Arrow>().swap(PartialSequence); vector<Arrow>().swap(FinishedSequence);//removing excess
    
    ArrowList.swap(NewArrowList);
    RemoveMod(ArrowList);
    int Write=0;
    for(int i=0; i<sizeAsInt(GeneratorList); i++ ) //Update GeneratorList
       {Gen G=GeneratorList[i]; 
        if (G.Idem %8 ==4) 
          {G.Idem=G.Idem/4 -1;  
	   GeneratorList[Write]=G;
           Write++;}
       }
    GeneratorList.erase(GeneratorList.begin()+Write, GeneratorList.end());
    
    vector<int> temp;  //Additional update     
    for(int i=0; i<sizeAsInt(MatchingList)-2; i++)
      {if (i==aa-3) temp.push_back(bb-2);
        if (i==bb-3) temp.push_back(aa-2);
       if (i !=aa-3 && i != bb-3) temp.push_back(MatchingList[i+2]-2); }
    MatchingList=temp;
    vector<int> temp2;
    for (int i=0; i<sizeAsInt(UpwardList)-2;i++) temp2.push_back(UpwardList[i+2]);
    UpwardList=temp2;
    Bridge=Bridge-1;

    ReName();
}










