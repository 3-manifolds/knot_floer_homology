#include "Alg.h"
#include<vector>
#include<string>
#include<array>
#include<map>
#include<algorithm>
#include<iostream>
using namespace std;


/////////////////////////
//  CROSSING BIMODULES
////////////////////////

// n is a crossing between 1 and 2*Bridge-1
// Bimodule Generators are North, East, West, South

enum{North, East, West, South};

inline  bool Extendable(idem x, int n, int Cor) 
{   if(Cor==North) return  !!(x & (1<< n));
    if(Cor==South) return !(x & (1<< n));
    if(Cor==West)  return ( !(x & (1<< n)) && (x & (1<<(n-1))) );
    if(Cor==East)  return ( !(x & (1<< n)) && (x & (1<<(n+1))) );
    return false; /* Avoid compiler warning */
}

inline  idem Extend(idem x, int n, int Cor)
{   if(Cor==North || Cor==South) return x;
    if(Cor==East) return x-(1<<n);
    if(Cor==West) return x+(1<<(n-1));
    return false; /* Avoid compiler warning */
}

////////////
//M1 Actions
////////////
    
void PosM1(Gen G, int n, int Cor) //Cor=East or West
{   idem I1=G.Idem;
    if(Extendable(I1, n, Cor))
      {Arrow ans;  
       ans.MonomialIndex=0; 
       ans.StartingGen=4*(G.Name)+Cor;
       ans.EndingGen  =4*(G.Name)+South;
       if(Cor==East) ans.Coeff=1; 
       if(Cor==West) ans.Coeff=Modulus-1; 
       NewArrowList.push_back(ans);}
}

void NegM1( Gen G, int n, int Cor) //Cor=East or West
{   idem I1=G.Idem;
    if(Extendable(I1, n, Cor))
      {Arrow ans;  
       ans.MonomialIndex=0; 
       ans.StartingGen=4*(G.Name)+South;
       ans.EndingGen  =4*(G.Name)+Cor;
       if(Cor==East) ans.Coeff=1; 
       if(Cor==West) ans.Coeff=Modulus-1; 
       NewArrowList.push_back(ans);}
}

void CurvedPosM1(Gen G, int n,  int k, int Cor) //Cor=East or West
{   idem I1=G.Idem; monomial Y=MonomialOne;
    if( ! Extendable(I1, n, Cor)) return; 
    Arrow ans;  
    idem I2=Extend(I1, n, Cor);
    ans.StartingGen  =4*(G.Name)+South;
    ans.EndingGen    =4*(G.Name)+Cor;
    if(Cor==East) ans.Coeff=1; 
    if(Cor==West) ans.Coeff=Modulus-1;
    Y[k-1]=1;
    int m=MonomialLookUp(Y);
    ans.MonomialIndex=m; 
    if(NonZero(I1,I2,m) ) NewArrowList.push_back(ans);
}

void CurvedNegM1(Gen G, int n,  int k, int Cor) //Cor=East or West
{   idem I1=G.Idem; monomial Y=MonomialOne;
    if(! Extendable(I1, n, Cor)) return;  
    Arrow ans; 
    idem I2=Extend(I1, n, Cor);
    ans.StartingGen  =4*(G.Name)+Cor;
    ans.EndingGen    =4*(G.Name)+South;
    if(Cor==East) ans.Coeff=1; 
    if(Cor==West) ans.Coeff=Modulus-1;
    Y[k-1]=1;
    int m=MonomialLookUp(Y); 
    ans.MonomialIndex=m;
    if(NonZero(I2,I1,m) ) NewArrowList.push_back(ans);
}

void M1Actions(int Crossing, int k1, int k2)//k1 matches with n, k2  with n+1
{   int n=abs(Crossing); bool Pos=(Crossing >0); bool Neg (Crossing <0);
    for (Gen G: GeneratorList)
      {if(Pos) {PosM1(G,n,East); PosM1(G,n,West);}
       if(Pos && k1 !=n+1) {CurvedPosM1(G,n,k1,East); CurvedPosM1(G,n,k2,West);}
       if(Neg) {NegM1(G,n,East); NegM1(G,n,West);}
       if(Neg && k1 !=n+1) {CurvedNegM1(G,n,k1,East); CurvedNegM1(G,n,k2,West);}
      }
}

////////////
//M2 Actions
////////////

  //works for both positive and negative actions:
void M2S(Arrow arrow1, int n)//Local weights are (U_n\cdot U_{n+1})^t 
{   idem I1=GeneratorList[arrow1.StartingGen].Idem; 
    idem I2=GeneratorList[arrow1.EndingGen].Idem; 
    if(Extendable(I1, n, South) && Extendable(I2, n, South))
       {monomial X=MonomialStore[arrow1.MonomialIndex]; 
        if (X[n]==X[n-1]) 
	   {arrow1.StartingGen=4*(arrow1.StartingGen)+South; 
            arrow1.EndingGen  =4*(arrow1.EndingGen)+South; 
	    arrow1.Coeff      = (Modulus- arrow1.Coeff)%Modulus; 
            NewArrowList.push_back(arrow1);}
       }
}     

  //For Positive Actions: 
  //local weights and corners
  //f1=1 means R_n, f1=-1 means L_n. 
  //f2 corresponds to wall n+1.
  //U1 and U2 gives powers of U_n and U_{n+1}
  //X=0 North, X=1 East, X=2 West
int LookBack(int f1, int f2, int U1, int U2, int X)
{   switch(X) {
    case 0: { 
      if  (            f1==f2            ) return 0;
      if  (f1==1  &&  f2==0   && U1< U2) return 1;
      if  (f1==1  &&  f2==0   && U1>=U2) return 2;
      if  (f1==0  &&  f2==-1  && U1<=U2) return 1;
      if  (f1==0  &&  f2==-1  && U1> U2) return 2;}
 
    case 1: {   
      if  (f1==0  &&  f2==0   && U1<=U2) return 1;
      if  (f1==0  &&  f2==0   && U1> U2) return 2;
      if  (f1==-1 &&  f2==0             ) return 0;
      if  (f1==0  &&  f2==1             ) return 0;}
  
    case 2: {
      if  (f1==0  &&  f2==0   && U1>=U2) return 2;
      if  (f1==0  &&  f2==0   && U1< U2) return 1;
      if  (f1==-1 &&  f2==0             ) return 0;
      if  (f1==0  &&  f2==1             ) return 0;}}
  
    return -1;
}

void  PosM2(Arrow arrow1, int n, int EndCor)//EndCor =North, East or West
{   idem I1=GeneratorList[arrow1.StartingGen].Idem; 
    idem I2=GeneratorList[arrow1.EndingGen].Idem;  
    if(! Extendable(I2, n, EndCor)) return; 
    int f1=LeftRight(I1,I2,n); 
    int f2=LeftRight(I1,I2,n+1);
    
    monomial X=MonomialStore[arrow1.MonomialIndex];
    I2=Extend(I2, n, EndCor); 
    int U1=X[n-1]; int U2=X[n]; 
    int StartCor=LookBack(f1,f2,U1,U2,EndCor);
    if  (StartCor==-1) return;
    if(! Extendable(I1, n, StartCor)) return;  else I1=Extend(I1, n, StartCor);
    if(TooFar(I1,I2)) return;
    int H1=2*U1+abs(f1); int H2=2*U2+abs(f2);
    if(StartCor ==East) H2--; 
    if(EndCor   ==East) H2++; 
    if(StartCor ==West) H1--; 
    if(EndCor   ==West) H1++;
    X[n-1]=H2/2;  X[n]=H1/2; 
    arrow1.StartingGen=4*(arrow1.StartingGen)+StartCor;
    arrow1.EndingGen  =4*(arrow1.EndingGen)  +EndCor;
    arrow1.MonomialIndex =MonomialLookUp(X);
    if(NonZero(I1,I2,arrow1.MonomialIndex) ) NewArrowList.push_back(arrow1);
}

//For Negative Actions:
int LookForward(int f1, int f2, int U1, int U2, int X)
{   switch(X){
    case 0: {
      if  (            f1==f2            )   return 0;
      if  (f1==-1 &&  f2==0   && U1< U2)   return 1;
      if  (f1==-1 &&  f2==0   && U1>=U2)   return 2;
      if  (f1==0  &&  f2==1   && U1<=U2)   return 1;
      if  (f1==0  &&  f2==1   && U1> U2)   return 2;}
    
    case 1: {
      if  (f1==0  &&   f2==0  && U1<=U2)   return 1;
      if  (f1==0  &&   f2==0  && U1 >U2)   return 2;
      if  (f1==1  &&   f2==0            )   return 0;
      if  (f1==0  &&   f2==-1           )   return 0;}
  
    case 2: {
      if  (f1==0  &&   f2==0  && U1>=U2)   return 2;
      if  (f1==0  &&   f2==0  && U1< U2)   return 1;
      if  (f1==1  &&   f2==0            )   return 0;
      if  (f1==0  &&   f2==-1           )   return 0;}}
  
    return -1;
}

void  NegM2(Arrow arrow1, int n, int StartCor)//StartCor =North, East or West
{   idem I1=GeneratorList[arrow1.StartingGen].Idem; 
    idem I2=GeneratorList[arrow1.EndingGen].Idem; 
    if(! Extendable(I1, n, StartCor)) return; 
    int f1=LeftRight(I1,I2,n);
    int f2=LeftRight(I1,I2,n+1);

    monomial X=MonomialStore[arrow1.MonomialIndex];
    I1=Extend(I1, n, StartCor);
    int U1=X[n-1]; int U2=X[n];
    int EndCor=LookForward(f1,f2,U1,U2,StartCor);
    if  (EndCor==-1) return;
    if(! Extendable(I2, n, EndCor)) return; else I2=Extend(I2, n, EndCor);
    if(TooFar(I1,I2)) return;   
    int H1=2*U1+abs(f1); int H2=2*U2+abs(f2);    
    if(EndCor   ==East) H2--;
    if(StartCor ==East) H2++;
    if(EndCor   ==West) H1--;
    if(StartCor ==West) H1++;
    X[n-1]=H2/2;
    X[n]=H1/2;
    arrow1.StartingGen=4*(arrow1.StartingGen)+StartCor;
    arrow1.EndingGen  =4*(arrow1.EndingGen)  +EndCor;
    arrow1.MonomialIndex =MonomialLookUp(X);
    if(NonZero(I1,I2,arrow1.MonomialIndex) ) NewArrowList.push_back(arrow1);
}

void M2Actions(int Crossing)
{   int n=abs(Crossing); bool Pos=(Crossing >0); bool Neg=(Crossing <0);
    for(Arrow A : ArrowList){
       if(Pos) {M2S(A,n); PosM2(A,n,North); PosM2(A,n,West); PosM2(A,n,East);} 
       if(Neg) {M2S(A,n); NegM2(A,n,North); NegM2(A,n,West); NegM2(A,n,East);}}
}

////PosM3 Actions ////

void  PosM3(Arrow arrow1, Arrow arrow2, int n) //StartCor=South
{   idem I1=GeneratorList[arrow1.StartingGen].Idem;  
    idem I2=GeneratorList[arrow1.EndingGen].Idem;
    idem I3=GeneratorList[arrow2.EndingGen].Idem;
    if ( ! Extendable(I1, n, South)) return;
    monomial m1=MonomialStore[arrow1.MonomialIndex];
    monomial m2=MonomialStore[arrow2.MonomialIndex];     
     
    int f1=LeftRight(I1,I2,n); 
    int f2=LeftRight(I1,I2,n+1); 
    int U1=m1[n-1]; int U2=m1[n];
    int g1=LeftRight(I2,I3,n); 
    int g2=LeftRight(I2,I3,n+1); 
    int V1=m2[n-1]; int V2=m2[n];
      
    for(int EndCor=0; EndCor <3; EndCor++) 
       {if( ! Extendable(I3, n, EndCor)) continue;
        idem I4=Extend(I3, n, EndCor);
        if(TooFar(I1, I4))  continue; 
	int x1=LookBack(g1,g2,V1,V2,EndCor);   if(x1==-1) continue;
        int x2=LookBack(f1,f2,U1,U2,x1);       if(x2==-1 || x2==0) continue;
        int x3=LookBack(g1+f1, g2+f2, U1+V1+(abs(g1)+abs(f1))/2, U2+V2+(abs(g2)+abs(f2))/2,EndCor); 
        if(x3==0 || x2==x3) continue;
        if(x3==-1 && EndCor==North) continue;
        if(x3==-1 && EndCor==East && !(f1==1 && U1==0 && f2==0 && U2==0 && g1==0 && V1==0 && g2==1)) continue;
        if(x3==-1 && EndCor==West && !(f1==0 && U1==0 && f2==-1 && U2==0 && g1==-1 && g2==0 && V2==0)) continue;  
        monomial X=MonomialOne;
        Arrow ans;
        ans.StartingGen=4*arrow1.StartingGen+South;
        ans.EndingGen  =4*arrow2.EndingGen+  EndCor;
        int p1=0; int p2=0; int p3=0;
        for(int i=0; i<2* Bridge; i++) 
	   {if(I1 &(1<<i)) p1++;
	    if(I2 &(1<<i)) p2++;
            if(I3 &(1<<i)) p3++;
            if((p1>p2 && p2<p3) || (p1<p2 && p2>p3))
 	       X[i]=m1[i]+m2[i]+1;
	    else X[i]=m1[i]+m2[i];}

        int H1=2*U1+ 2*V1+ abs(f1)+abs(g1)-1;
        int H2=2*U2+ 2*V2+ abs(f2)+abs(g2)-1;
        if(EndCor==East) H2++;
        if(EndCor==West) H1++;
        X[n-1]=H2/2; X[n]=H1/2;
        ans.MonomialIndex =MonomialLookUp(X);
        if( NonZero(I1,I4,ans.MonomialIndex)==false) continue;
        if (x2==2) ans.Coeff=      ((arrow1.Coeff)*(arrow2.Coeff))%Modulus;                  
        if (x2==1) ans.Coeff= ( (Modulus-(arrow1.Coeff)) * (arrow2.Coeff) )%Modulus;                  
        NewArrowList.push_back(ans); 
       }
}

///Negative M3 actions

void  NegM3(Arrow arrow1, Arrow arrow2, int n)//EndCor=South
{   idem I1=GeneratorList[arrow1.StartingGen].Idem;  idem I2=GeneratorList[arrow1.EndingGen].Idem;
    idem I3=GeneratorList[arrow2.EndingGen].Idem;
    if ( ! Extendable(I3, n, South)) return;
      
    monomial m1=MonomialStore[arrow1.MonomialIndex];
    monomial m2=MonomialStore[arrow2.MonomialIndex];     
      
    int f1=LeftRight(I1,I2,n); 
    int f2=LeftRight(I1,I2,n+1); 
    int U1=m1[n-1]; int U2=m1[n];
    int g1=LeftRight(I2,I3,n);
    int g2=LeftRight(I2,I3,n+1);
    int V1=m2[n-1]; int V2=m2[n]; 

    for(int StartCor=0; StartCor <3; StartCor++) 
         {if( ! Extendable(I1, n, StartCor)) continue;
          idem I0=Extend(I1, n, StartCor);
          if( TooFar(I0, I3)) continue; 
          int x1=LookForward(f1,f2,U1,U2,StartCor);             if(x1==-1) continue;
          int x2=LookForward(g1,g2,V1,V2,x1);                   if(x2==-1 || x2==0) continue;
          int x3=LookForward(f1+g1, f2+g2, U1+V1+(abs(f1)+abs(g1))/2, U2+V2+(abs(f2)+abs(g2))/2,StartCor); 
          if((x3==0) || (x2==x3)) continue;
          if(x3==-1 && StartCor==North) continue;
          if(x3==-1 && StartCor==East && !(f1==0 && U1==0 && f2==-1 && g1==-1 && V1==0 && g2==0 && V2==0)) continue;
          if(x3==-1 && StartCor==West && !(f1==1 && U2==0 && f2==0  && g1==0  && V1==0 && g2==1 && V2==0)) continue;  
          monomial X=MonomialOne;
          Arrow ans;
          ans.StartingGen=4*arrow1.StartingGen+StartCor;
          ans.EndingGen  =4*arrow2.EndingGen+  South;
      
          int p1=0; int p2=0; int p3=0;
          for(int i=0; i<2* Bridge; i++) 
	    {if(I1 &(1<<i)) p1++;
	     if(I2 &(1<<i)) p2++;
             if(I3 &(1<<i)) p3++;
             if((p1>p2 && p2<p3) || (p1<p2 && p2>p3))
 	       X[i]=m1[i]+m2[i]+1;
	     else X[i]=m1[i]+m2[i];}
	  int H1=2*U1+ 2*V1+ abs(f1)+ abs(g1)-1;
          int H2=2*U2+ 2*V2+ abs(f2)+ abs(g2)-1;
          if(StartCor==East) H2++;
          if(StartCor==West) H1++;
          X[n-1]=H2/2;
          X[n]=H1/2;
          ans.MonomialIndex =MonomialLookUp(X);

          if(NonZero(I0,I3, ans.MonomialIndex)==false) continue;
          if(x2==2) ans.Coeff= ((arrow1.Coeff)*(arrow2.Coeff))%Modulus;                  
          if(x2==1) ans.Coeff= ((Modulus-(arrow1.Coeff))*(arrow2.Coeff))%Modulus;                  
          NewArrowList.push_back(ans); 
	 }
}


void  M3Actions(int Crossing)
{   int x=sizeAsInt(ArrowList); 
    if(x==0) return;
    int y=sizeAsInt(GeneratorList);
    vector<int> X(y+1,0);
    for(int i=0; i<x; i++) X[ArrowList[i].StartingGen]++;
    X[y]=x;
    for(int i=y-1; i>=0 ;i--) X[i]= X[i+1]-X[i]; 
    vector<int> Index(x,0);
    for(int i=0; i<x; i++) {int a=ArrowList[i].StartingGen; Index[X[a]]=i; X[a]++;}
    for(int i=0; i<x; i++) 
	  {Arrow arrow1=ArrowList[i]; int a=arrow1.EndingGen; 
	    int j=0; if(a>0) j=X[a-1];	    
	    while(j<X[a])
	      {Arrow arrow2=ArrowList[Index[j]]; j++;
             if(Crossing >0) PosM3(arrow1,arrow2, Crossing);
             if(Crossing <0) NegM3(arrow1,arrow2,-Crossing);} 
	  }
}






void AfterCrossing(int Crossing)
{   int n; if (Crossing >0) n=Crossing; else  n=-Crossing;
    int PM; if(Crossing >0) PM=1; else PM=-1;
    int aa=MatchingList[n-1]; int bb=MatchingList[n];
    if(aa != n+1) {MatchingList[n-1]=bb; MatchingList[n]=aa; MatchingList[aa-1]=n+1; MatchingList[bb-1]=n;}
    int Sign1=UpwardList[n-1]; int Sign2=UpwardList[n];
    UpwardList[n-1]=Sign2; UpwardList[n]=Sign1;
    
    int x=sizeAsInt(GeneratorList);
    for(int i=0; i<x; i++)
      for(int Cor=0; Cor<4; Cor++)
	{Gen G=GeneratorList[i]; 
	  if( !Extendable(G.Idem, n, Cor) ) continue; 
	  G.Idem=Extend(G.Idem, n, Cor);
          int A=G.Alexander; int M=G.Maslov;
          if((Sign1==1 && Sign2==1 && Cor==North) || (Sign1==0 && Sign2==0 && Cor==South)) {A=A-PM; M=M-PM;}
          if((Sign1==1 && Sign2==0 && Cor==West)  || (Sign1==0 && Sign2==1 && Cor==East )) {A=A+PM; M=M+PM;}
          if((Sign1==1 && Sign2==0 && Cor==East)  || (Sign1==0 && Sign2==1 && Cor==West )) {A=A-PM;}                    
          if((Sign1==1 && Sign2==1 && Cor==South) || (Sign1==0 && Sign2==0 && Cor==North)) {A=A+PM;}                    
          G.Alexander=A; G.Maslov=M; G.Name=4*G.Name+Cor;  
          NewGeneratorList.push_back(G);
	}
        
    NewArrowList.reserve(2*ArrowList.size());
    M1Actions(Crossing,aa,bb); 
    M2Actions(Crossing);  
    M3Actions(Crossing); 
    
    vector<Gen>().swap(GeneratorList); 
    GeneratorList.swap(NewGeneratorList);
    vector<Gen>().swap(NewGeneratorList);   
    
    vector<Arrow>().swap(ArrowList);   
    ArrowList.swap(NewArrowList);
    vector<Arrow>().swap(NewArrowList);
  
    ReName();
    RemoveMod(ArrowList);
}

