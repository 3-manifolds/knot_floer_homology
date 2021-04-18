#include "Alg.h"
#include "Diagrams.h"
#include<vector>
#include<string>
#include<set>
#include<map>
#include<algorithm>
#include<iostream>

using namespace std;

  // The contribution of all upper Kauffman states in a given idempotent and Alexander grading
struct Term{idem Idem; int Alexander; int Coeff;
};

void Update(vector<Term> & Old)
{   if(Old.size()==0) return;
    sort(Old.begin(), Old.end(), [](Term a, Term b) 
      {return (a.Idem<b.Idem || (a.Idem==b.Idem && a.Alexander<b.Alexander));
      } ); //sorting by idempotent and Alexander grading.

    Term a=Old[0]; a.Coeff=0; int Write=0;
    for(Term b: Old)
       if(a.Idem==b.Idem && a.Alexander==b.Alexander) a.Coeff+=b.Coeff;
       else if(a.Coeff==0) a=b;
       else {Old[Write]=a; Write++; a=b;}
    if(a.Coeff!=0)
	Old[Write]=a;
    Write++;
    Old.erase(Old.begin()+Write, Old.end()); 
}

//enum{North, East, West, South};

bool ExtendableA(idem x, int n, int Cor) 
{   if(Cor==0)  return  !!(x & (1<< n));
    if(Cor==3)  return !(x & (1<< n));
    if(Cor==2)  return ( !(x & (1<< n)) && (x & (1<<(n-1))) );
    if(Cor==1)  return ( !(x & (1<< n)) && (x & (1<<(n+1))) );
    return false;
}

idem ExtendA(idem x, int n, int Cor)
{   if(Cor==1)  return x-(1<<n);
    if(Cor==2)  return x+(1<<(n-1)); 
    if(Cor==0 || Cor==3)  return x;
    return false;
}

vector<Term> AfterCrossingAlt(vector<Term> Old, int Crossing)
{   vector<Term> New; int n=abs(Crossing);
    int PM; if(Crossing >0) PM=1; else PM=-1;
    int Sign1=UpwardList[n-1]; int Sign2=UpwardList[n];
    UpwardList[n-1]=Sign2; UpwardList[n]=Sign1;
    for(int Cor=0; Cor<4; Cor++)
      for(Term G: Old)
         {if( !ExtendableA(G.Idem, n, Cor)) continue; 
          G.Idem=ExtendA(G.Idem, n, Cor);
          int A=G.Alexander; int M=G.Coeff;
          if((Sign1==1 && Sign2==1 && Cor==0) || (Sign1==0 && Sign2==0 && Cor==3)) {A=A-PM; M=-M;}
          if((Sign1==1 && Sign2==0 && Cor==2)  || (Sign1==0 && Sign2==1 && Cor==1 )) {A=A+PM; M=-M;}
          if((Sign1==1 && Sign2==0 && Cor==1)  || (Sign1==0 && Sign2==1 && Cor==2 )) {A=A-PM;}                               
          if((Sign1==1 && Sign2==1 && Cor==3) || (Sign1==0 && Sign2==0 && Cor==0)) {A=A+PM;}                    
          G.Alexander=A; G.Coeff=M;  
          New.push_back(G);
	 }
    Update(New);   
    return New; 
}  

vector<Term> AfterMinAlt(vector<Term> Old)
{   vector<Term> New;  
    for(Term G: Old )
       {if(G.Idem %8 !=4) continue; 
        else  {G.Idem= G.Idem/4 -1; New.push_back(G);} 
       }
    vector<int> temp;
    for (int i=0; i<sizeAsInt(UpwardList)-2;i++) 
       temp.push_back(UpwardList[i+2]);
    UpwardList=temp;
    Bridge=Bridge-1;
    Update(New);
    return New;
}

idem ExtendXA(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 % pow; 
    int y= (I1/(2*pow))*8*pow;
    return x + 3*pow + y;
}

idem ExtendYA(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 %pow; 
    int y= (I1/(2*pow))*8*pow;
    return x + 6*pow + y;
}

idem ExtendZA(idem I1, int Position)
{   int pow=1<<(Position-1); 
    int x= I1 %pow ; 
    int y= (I1/(2*pow))*8*pow;
    return x + 2*pow + y;
}

 vector<Term> AfterMaxAlt(vector<Term> Old, int Position)
{   vector<Term> New;
    vector<int> NewUp(2*Bridge+2);
    for(int i=0; i<Position -1;i++) NewUp[i]=UpwardList[i];
    for(int i=Position+1; i<2*Bridge+2; i++)  NewUp[i]=UpwardList[i-2];
    NewUp[Position-1]=0; NewUp[Position]=0;  
    UpwardList=NewUp;
    Bridge++;
    for(Term G: Old)   
       {idem I1=G.Idem; 
        bool XY = !!(I1 &(1<<(Position-1))); bool Z = !XY ;
        if(Z) {G.Idem=ExtendZA(I1, Position);  New.push_back(G);}
        if(XY) {G.Idem=ExtendXA(I1,Position);  New.push_back(G); 
	       G.Idem=ExtendYA(I1,Position);  New.push_back(G); }
       }    
    Update(New);  
    return New;;
}

int Signature (PlanarDiagram Diag)//for an alternating projection
{   vector<int> PD = Diag.GetListOfTuples(); 
  int x=sizeAsInt(PD)/4;
    int Positive=0; //number of positive crossings
    for(int i=0; i<x; i++) if((PD[4*i+1]+2*x-PD[4*i+3])%(2*x)==1) Positive++; 

    int Black=0; //number of black regions
    set<int> Free; for(int i=1; i<=2*x; i++) Free.insert(i); //Free contains all the edges, ( PD code uses edges from 1 to 2*x)
    while(Free.size()>0)
      {Black++; int E=*(Free.begin()); Free.erase(E); int Old=E; int New=0;
       while(E !=New ) //finding the contours of the the new black region, and erasing them from Free
	  {int j=0; while(PD[2*j] != Old) j++; 
	  New=PD[2*j+1]; Free.erase(New); Old=New;
	  }
      }
    return Black-Positive-1;
}
  
void KnotFloerForAlternatingKnots(PlanarDiagram Diag, ostream & os)
{   vector<int> Morse = Diag.GetSmallGirthMorseCode(200).GetMorseList();
    Bridge=1; Term G1;  G1.Alexander=0; G1.Coeff=1; G1.Idem=2;
    vector<Term> Current; Current.push_back(G1);       
    UpwardList.clear(); 
    if(Morse[0]==1000) {UpwardList.push_back(1); UpwardList.push_back(0);}
    if(Morse[0]==1001) {UpwardList.push_back(0); UpwardList.push_back(1);}

    int Steps=sizeAsInt(Morse);
    for(int i=2; i< Steps -1 ; i++)
	{if (Morse[i]==1000)                                           
	    {int Position= Morse[i+1]; Current=AfterMaxAlt(Current,Position); 
	     UpwardList[Position-1]=1; UpwardList[Position]=0; i++;}
          
         else if (Morse[i]==1001)                                           
	    {int Position=Morse[i+1]; Current=AfterMaxAlt(Current,Position);
	     UpwardList[Position-1]=0; UpwardList[Position]=1; i++; }

         else if (Morse[i]==-1000)                                          
             Current=AfterMinAlt(Current);  
		
         else if  (Morse[i]< 2*Bridge && Morse[i] >-(2*Bridge) && Morse[i] != 0) 
	     Current=AfterCrossingAlt(Current, Morse[i]);   
	}

    int delta=-Signature(Diag)/2;
    map<int,int> Range;
    for(int i=0; i<sizeAsInt(Current); i++) 
	{Term G=Current[i]; Range[G.Alexander]=G.Coeff;}
      
    int TotalRank=0;  
    bool LSpaceKnot=true;
    os<<"Ranks in Alexander, Maslov bigradings :"<<endl;
    
    for(auto X: Range) 
	{int a=X.first; int b=X.second; 
         if(b<0) b=-b; 
 	 os<<b<<"   ("<<a/2<<","<<a/2-delta<<")"<<endl; 
         TotalRank+=b;    
         if(b !=1) LSpaceKnot=false;
        }
    os<<"Total rank : "<<TotalRank<<endl;
    int MaxAlex= -(*Range.begin()).first;
    int LeadingCoeff= (*Range.begin()).second;
    os<<"Seifert genus : "<<MaxAlex/2<<endl;
    if(LeadingCoeff==1 || LeadingCoeff==-1) os<<"Fibered : Yes\n";
    else os<<"Fibered : No\n";       
    if(LSpaceKnot) os<<"L-space knot : Yes\n";
    else os<<"L-space knot : No\n";
    int epsilon, tau, nu;
    tau = delta;
    if(delta > 0){
	epsilon = 1;
	nu = delta;
    }
    if(delta == 0){
	epsilon = 0;
	nu = 0;
    }
    if(delta < 0){
	epsilon = -1;
	nu = delta + 1;
    }
    os<<"Tau : "<<tau<<endl;
    os<<"Nu : "<<nu<<endl;
    os<<"Epsilon : "<<epsilon<<endl;
}











