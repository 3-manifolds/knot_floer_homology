#include <iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<iostream>
#include "Alg.h"
#include "Diagrams.h"

using namespace std;

void MorseCode::Print(ostream & os) const {
  os<<"Morse Code: ";
  for (int i = 0; i < sizeAsInt(MorseList); i++)
    if (MorseList[i] >999)
      os<<"Max("<< MorseList[++i]<<"), ";
    else if (MorseList[i] >-1000)
      os<<MorseList[i]<<", ";
    else if (i < sizeAsInt(MorseList)-1)
      os<<"Min, ";   
    else 
      os<<"Min.\n";
  os<<"Girth: "<<Girth<<"\n\n";
}

void PlanarDiagram::Print(ostream & os) const {
  const vector<int> &PD = ListOfTuples;
  os<<"Planar Diagram: [";
  for (int i = 0; i < sizeAsInt(PD)/4; i++)
    os<<"["<<PD[4*i]<<", "<<PD[4*i+1]<<", "<<PD[4*i+2]<<", "<<PD[4*i+3]<<"], ";
  os<<"]\n";
} 

bool PlanarDiagram::NotValid() const { //partial check on data, also checks if it has more than 1 component
  const vector<int> &PD = ListOfTuples; 
  int y = sizeAsInt(PD);
  if (y == 0) 
    return true;
  if (y%4 != 0) 
    return true;
  vector<int> CrList = PD;
  sort(CrList.begin(), CrList.end());
  int x=y/4;
  for (int i = 0; i < 2*x; i++)
    if (CrList[i] != (i/2)+1) 
      return true;
  for (int i = 0; i < x; i++)
    if (PD[4*i]-PD[4*i+2] != -1 &&  PD[4*i]-PD[4*i+2] != 2*x-1) 
      return true;
    else if (abs(PD[4*i+1] -PD[4*i+3]) != 1 && abs(PD[4*i+1]-PD[4*i+3]) != 2*x-1)
      return true; 
  return false;
}
  
bool PlanarDiagram::Alternating() const {
  const vector<int> &PD = ListOfTuples;  
  int x = sizeAsInt(PD)/4; 
  int a = PD[0]%2;
  for (int i = 1; i < x; i++)
    if (a != (PD[4*i]%2)) 
      return false;
  return true;
}
  
bool PlanarDiagram::R1Reducible() const {
  const vector<int> &PD = ListOfTuples;
  int x = sizeAsInt(PD)/4;
  for (int i = 0; i < x; i++)
    if (PD[4*i] == PD[4*i+1] || PD[4*i+1] == PD[4*i+2] || PD[4*i+2]==PD[4*i+3] || PD[4*i+3]==PD[4*i]) 
      return true;
  return false;
}

PlanarDiagram::PlanarDiagram(const string &S){
  ListOfTuples = vector<int>();  
  int x = '0'; 
  int y = '9';
  int a = 0; 
  bool BuildingInt = false; 
  for (int i = 0; i < sizeAsInt(S); i++) 
      if (S[i] < x || S[i] > y) {
	if (BuildingInt == true) { //just finished reading an integer  
	  ListOfTuples.push_back(a);  a = 0; BuildingInt = false;
	}
      }
      else {  
          a = 10*a+(S[i]-x); BuildingInt = true;
      }
    
    if (ListOfTuples.size() > 0) { //normalizing the indexes for the edges
      int smallest=ListOfTuples[0];
      for (int i: ListOfTuples) 
        if (i < smallest) smallest = i;
      for (int & i: ListOfTuples) i = i-smallest+1;
    }
}

MorseCode PlanarDiagram::GetSmallGirthMorseCode(int MaxNumberOfTries) const {
  const vector<int> &PD=ListOfTuples;
  int x=sizeAsInt(PD)/4;
  int SmallestGirth=10000; 
  long long Complexity=1000000000;
  vector<int> MorseList;
  int R=min(100+x*x, MaxNumberOfTries); //randomly looking for a reasonable MorseList presentation 
    for(int T=0; T< R; T++)
        {int B=4; //computes the maximal intersection number with y=t in this cycle
         int FirstC=rand()% x;  //the first crossing used the Morse presentation
         vector<int> TempMorseList(5);
         long long TempComplexity=0;
      
         vector<int> temp(4); //this will list the name of the strands that intersect a y=t line
         int Shift=rand()%4;
         temp[0]=PD[4*FirstC +(Shift%4) ]; //example, planar crossing had X[1,8,2,9] and Shift =0, Then temp={1,8,2,9}. 
         temp[1]=PD[4*FirstC+ ((Shift+1)%4)]; 
         temp[2]=PD[4*FirstC+ ((Shift+2)%4)];
         temp[3]=PD[4*FirstC+ ((Shift+3)%4)];
         int lastStrand=temp[3]; //if temp={1,8,2,9} then 9 will remain to be the last coordinate while temp grows and shrinks
         
         if(temp[2]% (2*x) == (temp[0]+1)%(2*x)) TempMorseList[0]=1000; else TempMorseList[0]=1001;
         TempMorseList[1]=1;
         if(temp[3]% (2*x) == (temp[1]+1)%(2*x)) TempMorseList[2]=1000; else TempMorseList[2]=1001; 
         TempMorseList[3]=3;
         if(Shift%2 ==0) TempMorseList[4]=2; else TempMorseList[4]=-2; //Ex. X(1,8,2,9) and shift 0,  gives {1000,1,1000,3,2}
 
         vector<int> CrossingUsed(x,0);//value 1 means it is used already
         CrossingUsed[FirstC]=1;
         for(int A=0;A<x-1;A++)  
           {vector<int> MostConnected; //we choose NextCrossing  
            int MaxCon=1; 
            for(int j=0;j<x; j++)
	       {if(CrossingUsed[j]==1) continue; 
		 if(A<x-2 &&  (PD[4*j]==lastStrand || PD[4*j+1]==lastStrand 
		 	       || PD[4*j+2]==lastStrand || PD[4*j+3]==lastStrand ) )
		             continue; //not gluing to the last strand
		  
		int Con=0; vector<int> Where;
                for(int k=0;k<4;k++)
		  {int q=(int)(temp.end()-find(temp.begin(), temp.end(), PD[4*j+k]));
                    if(q>0) {Con++; Where.push_back(q);}}
		if(Con==0) continue;                
                sort(Where.begin(), Where.end());
                if(Where[Con-1] - Where[0]> Con-1) continue; //crossing j attaches to temp in disjoint intervalls 
		if(Con==MaxCon)  MostConnected.push_back((int)j);
                else if(Con>MaxCon) {MaxCon=Con; MostConnected.clear(); MostConnected.push_back((int)j);}
	       }  
	    int aa=sizeAsInt(MostConnected);
             if (aa == 0) { //problem with planar diagram
               vector<int> Empty =vector<int>();
	       return MorseCode(Empty, -1);
	     }     
             int a1=rand()% aa;
             int NextC =MostConnected[a1]; 
             CrossingUsed[NextC]=1;
             vector<int> V(4); for(int k=0; k<4; k++) V[k]=PD[4*NextC +k];
             
             int t=0;
             while(t<sizeAsInt(temp) && temp[t] != V[0] && temp[t] != V[1] 
		   && temp[t] != V[2] && temp[t] != V[3])
	            t++;
             int FirstPosition =t; 
             int k=0; while(V[k] != temp[FirstPosition]) k++;//k is between 0 and 3;
	     
	     long long w=A*temp.size()*temp.size();
	     if(MaxCon==3) w=(FirstPosition+1)*w;
	     if(MaxCon==4) w=(2*FirstPosition+1)*w;
             TempComplexity+=w; //adding the "cost" of tensoring with  this bimodule

             if(MaxCon==2 && k%2==0) TempMorseList.push_back((int)(FirstPosition+1)); //adding a positive crossing 
             if(MaxCon==2 && k%2==1) TempMorseList.push_back(-(int)(FirstPosition+1));//or a negative crossing
             
             if(MaxCon==1)   //this adds a maximum and a crossing 
                {int a=(k+1)%4; int b=(k+3)%4;
                 if( (V[a]+1)%(2*x) == V[b]%(2*x) ) TempMorseList.push_back(1000); //maximum oriented to the right 
                 else   TempMorseList.push_back(1001); //or oriented to the left
	         TempMorseList.push_back((int)(FirstPosition+2)); //where to put the maximum
                 if(k%2==0) TempMorseList.push_back((int)(FirstPosition+1)); //adding a positive crossing 
                 if(k%2==1) TempMorseList.push_back(-(int)(FirstPosition+1));//or negative crossing 
		}

             if(MaxCon >2 && k%2==0)
	       {TempMorseList.push_back(-(int)(FirstPosition+2)); //adding a crossing
		 for(int a=(int)FirstPosition; a>0;a--)  //repeated Reidemeister 2 moves pushing the minimum to the left 
	          {TempMorseList.push_back(-a); TempMorseList.push_back(-a-1);}
	           TempMorseList.push_back(-1000); //adding the minimum at left
		}
              
	     if(MaxCon >2 && k%2==1) //similar as before, but crossing is positive
               {TempMorseList.push_back((int)(FirstPosition+2)); 
		 for(int a=(int)FirstPosition; a>0;a--) 
	         {TempMorseList.push_back(a); TempMorseList.push_back(a+1);}
	          TempMorseList.push_back(-1000);
	       }
	
            if(MaxCon==4) //adding an extra minimum
	      {for(int a=(int)FirstPosition; a>0;a--)
	          {TempMorseList.push_back(a);TempMorseList.push_back(a+1);}
	           TempMorseList.push_back(-1000);
	       }
               
            for(int i=0;i<MaxCon;i++) //updating temp
                temp.erase(temp.begin()+FirstPosition);
            for(int i=0;i<4-MaxCon ;i++)
            temp.insert(temp.begin()+FirstPosition+i, V[(k+i+1)%4 ]); 
                
            if (sizeAsInt(temp)>B) B=sizeAsInt(temp); 
            if(B> SmallestGirth) break; //if this partial sequence is already bad, we start the next cycle
	   }
            	
	 if(B < SmallestGirth || (B==SmallestGirth  &&  TempComplexity< Complexity ))//saving the impoved MorseList:
           {SmallestGirth=B;  Complexity=TempComplexity; MorseList=TempMorseList;}
	}
        
    return MorseCode(MorseList, SmallestGirth);
}

