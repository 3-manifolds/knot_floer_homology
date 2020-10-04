#include "Alg.h"
#include <iostream>
#include<vector>
#include<algorithm>
#include<map>
using namespace std;


//contracts short differentials in the D-module
//Precondition: GeneratorList[i].Name=i;
void Simplify()
{   int x=GeneratorList.size(); int y=ArrowList.size(); 
    if(y==0) return;
      
    //Preparation to contract.
    //Maps1[i] will hold indexes of arrows coming out of Generator i.
    //Maps2[i] will hold indexes of arrows to Generator i:
    vector<int>*  Maps1= new vector<int> [x];
    vector<int>*  Maps2= new vector<int> [x];
    for(int i=0;  i<y; i++)
         {int a=ArrowList[i].StartingGen;
          int b=ArrowList[i].EndingGen;
          (Maps1[a]).push_back(i);
          (Maps2[b]).push_back(i); }
      
    //Candidates holds the indexes of Generators that have short Arrow
    //coming out of them. It is sorted by the size of Maps1[i].
    vector<pair<int, int> > Candidates;
    for(int i=0; i<x; i++)
	{vector<int> X=Maps1[i]; if(X.size()==0) continue;
	 bool FoundShortArrow=false;  
	 for(int j=0; j<X.size(); j++)
	    {Arrow Y=ArrowList[X[j]]; int a=Y.StartingGen; int b=Y.EndingGen;
             if (Y.Coeff !=0 and Y.MonomialIndex==0 and
                 GeneratorList[a].Idem == GeneratorList[b].Idem)
	       FoundShortArrow=true;
	    }
	 if(FoundShortArrow) Candidates.push_back(make_pair(X.size(), i)); 
	}
     
    if(Candidates.size()>0)    
    sort(Candidates.begin(),Candidates.end());
      
    vector<int> DeletedNames(x,0);
    int deleted=0;
    int OldSize=y;
    int CurrentSize=y;//current size of ArrowList
    //we start the contracting algorithm:
    for(int i=0; i<Candidates.size(); i++)      
	{int From=(Candidates[i]).second; //Trying to eliminate "From"
	 vector<int> X=Maps1[From]; if( X.size()==0 ) continue;
         int index=-1; int To;  int Connectivity=1000000;
         //looking for a short differential out of "From"
         for(int j=0; j<X.size(); j++) 
	   {Arrow Arr=ArrowList[X[j]];  int c=Arr.EndingGen;
            if( Arr.Coeff !=0 
                and GeneratorList[From].Idem == GeneratorList[c].Idem 
	        and Arr.MonomialIndex==0  and (Maps2[c]).size() <Connectivity) 
	    {index=j; To=c; Connectivity=(Maps2[c]).size();}
	   }
         if(index==-1) continue; 
         Arrow Y=ArrowList[X[index]]; 
         //found a short differential Y between From and To
         //Contract Y and eliminate To and From:
         int Inverse=0; //used when contracting Y;
         for(int t=1; t<Modulus; t++) 
            if( (t* (Modulus-Y.Coeff))% Modulus==1) Inverse=t;
         
         vector<int> P=Maps2[To]; //indexes of Arrows that point to To
         //Update the arrows:
         for(int  i1 : P)//creating new arrows from zig-zags
	   for(int  i2: X) 
               {Arrow A=ArrowList[i1]; Arrow B=ArrowList[i2];
	        if(A.Coeff==0 or A.StartingGen==From 
                   or B.Coeff ==0 or B.EndingGen==To) continue;
	        int a=A.StartingGen; int b=B.EndingGen;
                idem I1=GeneratorList[a].Idem; 
                idem I2=GeneratorList[To].Idem; 
                idem I3=GeneratorList[b].Idem;
                if(TooFar(I1,I3) ) continue;
                int m=Mult(I1, I2, I3, A.MonomialIndex, B.MonomialIndex);
                if(! NonZero(I1, I3, m)) continue;
                Arrow Q; 
                //an extra term from a zig-zag that involve A, Y and B. 
                Q.StartingGen=a; 
                Q.EndingGen=b;
                Q.MonomialIndex=m;
                Q.Coeff=((A.Coeff)*(B.Coeff)*(Inverse)) % Modulus; 
                // Either  Q existed before (with perhaps different 
                //coefficient) or Q is new. We add or update:
	        int s=0; int f=Maps1[a].size();
                while(s<f and (ArrowList[Maps1[a][s]].EndingGen !=b 
                      or ! Equal(ArrowList[Maps1[a][s]], Q)) ) s++;
	        if(s==f) //Q is new
	          {ArrowList.push_back(Q); 
		   (Maps1[a]).push_back(CurrentSize); 
		   (Maps2[b]).push_back(CurrentSize);
	           CurrentSize++; }
                 else ArrowList[Maps1[a][s]].Coeff =
                       (ArrowList[Maps1[a][s]].Coeff+ Q.Coeff)% Modulus; 
	       }
		 
	 //at this point all the zig-zag arrows are added. 
         //Next step is to make the arrows that connect to To or From 
         //invisible in ArrowList:	
         vector<int> V1;  
         for(int r: Maps1[From]) V1.push_back(r);
         for(int r: Maps1[To] ) V1.push_back(r);
         for(int r: Maps2[From] ) V1.push_back(r);
         for(int r: Maps2[To] )   V1.push_back(r);
  
         for(int r: V1) ArrowList[r].Coeff=0; //now there are "invisible"
         DeletedNames[To]=1;
         DeletedNames[From]=1;
         deleted=deleted+2; 
         
         //This finished one cycle and eliminated 2 generators. 
         //Before starting a new cycle:  
         //If the size of the ArrowList increased too much a 
         //rearrange the data: 
         if(CurrentSize > 2*OldSize + 10000 and i<Candidates.size()-1) 
         // deleting those with Coeff=0:
	    {for(int j=0; j< x; j++) {Maps1[j].clear(); Maps2[j].clear();}
	     int Write=0;
	     for(Arrow A: ArrowList)
		{ if(A.Coeff !=0) 
                      {ArrowList[Write]=A; 
		       Maps1[A.StartingGen].push_back(Write); 
                       Maps2[A.EndingGen].push_back(Write);
                       Write++;} 
                }
             ArrowList.erase(ArrowList.begin()+Write, ArrowList.end() ); 
	     CurrentSize=Write; OldSize=CurrentSize; 
             
             //Candidates[0] to Candidates[i] are dealt with.
             //Reorder the remaining part of the Candidate list:
             for(int j=i+1; j<Candidates.size();j++)
	     Candidates[j].first=Maps1[Candidates[j].second].size();
             sort(Candidates.begin()+i+1, Candidates.end());
	    }
	}

    delete[] Maps1;
    delete[] Maps2;
      
    if(deleted>0) 
      {int Write=0;
       for(int i=0; i<x ; i++)
	   if(DeletedNames[i]==0) 
                {GeneratorList[Write]=GeneratorList[i]; Write++;}
       GeneratorList.erase(GeneratorList.begin()+Write, GeneratorList.end());        
       Write=0;
       for(int j=0; j<ArrowList.size(); j++)
	   if(ArrowList[j].Coeff !=0) {ArrowList[Write]=ArrowList[j]; Write++;}
        ArrowList.erase(ArrowList.begin()+Write, ArrowList.end());
      }
      
    ReName();//Arranging PostCondition: GeneratorList[i].Name=i.
 }




