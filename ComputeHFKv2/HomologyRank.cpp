#include "Alg.h"
#include <iostream>
#include<vector>
#include<algorithm>
#include<map>
using namespace std;

int  HomologyRank(const ChainComplex & OldComplex)//used for Tau, NuPlus, NuMinus
{   int x=OldComplex.Generators.size(); int y=OldComplex.Differential.size(); int deleted=0;
    if(y==0) return x;
    int OldSize=y;
    int CurrentSize=y;
 
    vector<int> GenList=OldComplex.Generators;
    vector<ChainArrow> ArrList= OldComplex.Differential;
    int prime=OldComplex.Prime;     
    int MaxName=0;
    for(int t=0; t<x; t++) 
         if(GenList[t] >MaxName) MaxName=GenList[t];
    vector<int> NewName(MaxName+1,-1);
    for(int t=0; t<x; t++)
	{ NewName[GenList[t]]=t; GenList[t]=t;}
    for(int t=0; t<y; t++)
	{ ArrList[t].StartingGen=NewName[ArrList[t].StartingGen];      
	  ArrList[t].EndingGen  =NewName[ArrList[t].EndingGen];}      
    
    vector<int>*  Maps1= new vector<int> [x];
    vector<int>*  Maps2= new vector<int> [x];
    vector<int> DeletedNames(x,0);

    for(int i=0;  i<y; i++)
         {int a=ArrList[i].StartingGen;
          int b=ArrList[i].EndingGen;
          (Maps1[a]).push_back(i);
          (Maps2[b]).push_back(i); }
      
    vector<int> X; 
    ChainArrow Y; 
    int From; int To; 
    vector<int> P; 
    vector<ChainArrow> NewDifferentials;
      
    vector<pair<int,int>> Candidates;
    for(int i=0; i<x; i++)
        if(Maps1[i].size()>0) Candidates.push_back( make_pair( Maps1[i].size() ,i)); 
    sort(Candidates.begin(),Candidates.end());
 
      //we start the contracting algorithm
    for(int i=0; i<Candidates.size(); i++)      
	{int r=(Candidates[i]).second;
         X=Maps1[r]; 
         if( X.size()==0 ) continue;
         int n=-1; int Connectivity=10000000;//looking for a short differential out of i  
         for(int j=0; j<X.size(); j++) 
	   {ChainArrow Arr=ArrList[X[j]];  //Arr connects i to ArrList[X[j]].EndingGen 
            if( Arr.Coeff !=0  and (Maps2[Arr.EndingGen]).size() <Connectivity) 
	      {n=j; Connectivity=(Maps2[Arr.EndingGen]).size();}}
         if(n==-1) continue; 
         else Y=ArrList[X[n]]; //found a  differential Y between From and To
         int Inverse=0; //used when contracting Y;
         for(int t=1; t<prime; t++) if( (t* (prime-Y.Coeff))% prime==1) Inverse=t;
         From =r; To= Y.EndingGen;
         P=Maps2[To]; //indexes of Arrows that point to To
         
         //Contract the short differential: Update the arrows and mark From and To as deleted:
         NewDifferentials.clear(); 
         for(int i1=0; i1< P.size(); i1++)//creating new arrows from zig-zags
               {ChainArrow temp1=ArrList[P[i1]]; 
		if (temp1.Coeff ==0 or temp1.StartingGen==From) continue;//we exclude Y and already deleted arrows  
	        for(int i2=0; i2< X.size(); i2++)
	              {ChainArrow temp2=ArrList[X[i2]];
		       if(temp2.Coeff ==0 or temp2.EndingGen==To) continue;//we exclude Y and already deleted arrows
	               int a=temp1.StartingGen; int b=temp2.EndingGen;
                       
                       ChainArrow Q; //an extra term from a zig-zag that involve temp1, Y and temp2. 
                       Q.StartingGen=a; 
                       Q.EndingGen=b;
                       Q.Coeff=((temp1.Coeff)*(temp2.Coeff)*(Inverse)) % prime; 
		       NewDifferentials.push_back(Q);
                      }
	       }
          
	 for(int m=0; m<NewDifferentials.size();m++)
	   {ChainArrow Q=NewDifferentials[m]; int a=Q.StartingGen; int b=Q.EndingGen;                   	               
                       // Either  Q existed before (with perhaps different coefficient) or Q is new. We add or update.
	               int s=0; int f=Maps1[a].size();
                       while(s<f and ArrList[Maps1[a][s]].EndingGen !=b ) s++;
		       if(s==f) //Q is new
	                     {ArrList.push_back(Q); 
			      (Maps1[a]).push_back(CurrentSize);// the arrow Q from a to b is in ArrList[CurrentSize] 
                              (Maps2[b]).push_back(CurrentSize);
	                      CurrentSize++; }
                       else ArrList[Maps1[a][s]].Coeff =(ArrList[Maps1[a][s]].Coeff+ Q.Coeff)% prime; //Q is old
	   }
	   
           //at this point all the zig-zag arrows are added. 
           //Next step is to make the arrows that connect to To or From invisible in ArrList:	
           vector<int> V1;  
           for(int i3=0; i3< Maps1[From].size(); i3++) V1.push_back(Maps1[From][i3]);
           for(int i3=0; i3< Maps1[To].size();   i3++) V1.push_back(Maps1[To][i3]);
           for(int i3=0; i3< Maps2[From].size(); i3++) V1.push_back(Maps2[From][i3]);
           for(int i3=0; i3<Maps2[To].size();    i3++) V1.push_back(Maps2[To][i3]);
  
           for(int i4=0; i4<V1.size();           i4++) ArrList[V1[i4]].Coeff=0; //now there are "invisible"
           DeletedNames[To]=1;
           DeletedNames[From]=1;
           deleted=deleted+2; 
           
           

           
           //If the size of the ArrowList increased too much and there are still elements to process in Candidates,
           //rearrange the data: 
           if(CurrentSize > OldSize + OldSize/3 + 10000 and i<Candidates.size()-1) // deleting those with Coeff=0.
	    {for(int j=0; j< x; j++) {Maps1[j].clear(); Maps2[j].clear();}
	     int v=0;
	     for(int j=0; j< CurrentSize; j++)
		{ChainArrow arrow1=ArrList[j]; 
		 if(arrow1.Coeff !=0) 
                      {ArrList[v]=arrow1; 
		       Maps1[arrow1.StartingGen].push_back(v); 
                       Maps2[arrow1.EndingGen].push_back(v);
                       v++;} 
                }
	      CurrentSize=v; OldSize=CurrentSize; ArrList.erase(ArrList.begin()+v, ArrList.end() ); 
	      //Then reorder the remaining part of the Candidate list:
              for(int j=i+1; j<Candidates.size();j++)
		Candidates[j].first=Maps1[Candidates[j].second].size();
              sort(Candidates.begin()+i+1, Candidates.end());
	    }

	}

    delete[] Maps1;
    delete[] Maps2;
    return x-deleted;
}

