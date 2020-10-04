#include "Alg.h"
#include "Diagrams.h"
#include<vector>
#include<map>
#include<unordered_map>
#include<iostream>
#include<fstream>
using namespace std;

vector<monomial> MonomialStore;
unordered_map <monomial,int,Hash> MonomialMap;
int Bridge;
int Modulus;
vector<int>  UpwardList;
vector<int>  MatchingList;
vector<Arrow> ArrowList;
vector<Arrow> NewArrowList;
vector<Gen>   GeneratorList;
vector<Gen>   NewGeneratorList;
monomial MonomialOne={0};

int main() {
    string FileName;
    int Prime;

    cout<<"Enter a filename that contains one or more knots in Planar Diagram format:\n";
    cout.flush();
    cin>>FileName;
    ifstream inFile;
    inFile.open(FileName);
    if (inFile.is_open() == false) { 
      cout<<"File is not found\n"; 
      exit(1); 
    }

    cout<<"and a prime number:  "<<endl;
    cout.flush();
    cin>>Prime;
    bool isPrime = true;
    if (Prime < 2) 
      isPrime = false;
    for (int i = 2; i < Prime; i++)
      if (Prime % i == 0) 
        isPrime = false;
    if(isPrime == false) {
      cout<<"Not a prime!\n";
      inFile.close();
      exit(1);
    }

    string number= to_string(Prime);
    string OutPutFileName=FileName+".mod"+number;
    string OutPutFileName2=FileName+".mod"+number+".Morse";
  
    clock_t beginT=clock();
    MonomialStore.push_back(MonomialOne);
    MonomialMap.insert(make_pair(MonomialOne, 0) );
    //srand(time(NULL));
  

     //reading the PDCode from the file:
     
    char ch; 
    string S; 
    vector<string> ListOfKnots;
    char ch2;
    while(inFile>>ch) 
        if(ch=='P') {
	  inFile>>ch; 
          if(ch=='D') break;//discards info before the first PD
	}          
    while(inFile>>ch) 
       if(ch !=' ' and ch !='P') S.push_back(ch);
       else if(ch =='P') {ListOfKnots.push_back(S); S.clear();}
      
    ListOfKnots.push_back(S); S.clear(); 
    inFile.close();
    
    ofstream outFile;
    outFile.open(OutPutFileName);
    ofstream outFile2;
    outFile2.open(OutPutFileName2);
    
    outFile<<"Computations from file :  "<<FileName<<"\n";
    outFile<<"with Coefficient = "<<Prime<<"\n\n\n";
    
        
    for(int t=0; t<ListOfKnots.size(); t++) {
      S=ListOfKnots[t]; 
      PlanarDiagram Diag = PlanarDiagram(S);    
      cout<<"Knot # "<<t+1<<endl;
      outFile<<"Knot # "<<t+1<<endl;
       
      if(Diag.NotValid()) {
	cout<<"Knot # "<<t+1<< " was not a valid Planar Diagram for a knot\n";
        outFile<<"Knot # "<<t+1<< " was not a valid Planar Diagram for a knot\n";
        cout.flush();
        outFile.flush(); 
        outFile.close(); 
        exit(1);
      }
      if(Diag.R1Reducible()) {
	cout<<"Knot # "<<t+1<< " Can be simplified by a Reidemeister 1 move \n";
        outFile<<"Knot # "<<t+1<< " Can be simplified by a Reidemeister 1 move \n";
        cout.flush();
        outFile.flush(); 
        outFile.close(); 
        exit(1);
      }

      MorseCode LastCheckBeforeComputation = Diag.GetSmallGirthMorseCode(1);
      if (LastCheckBeforeComputation.GetMorseList().size() == 0) {
        cout<<"Knot # "<<t+1<< " Incorrect Planar Diagram\n";
        outFile<<"Knot # "<<t+1<< " Incorrect Planar Diagram\n";
        cout.flush();
        outFile.flush(); 
        outFile.close(); 
        exit(1);
      }


      if(Diag.Alternating()) {
        MorseCode M = Diag.GetSmallGirthMorseCode(10);
        outFile2<< "Knot # "<<t+1<< "\n\n";
        M.Print(outFile2);     
        cout<<"Projection is alternating.\n";
        KnotFloerForAlternatingKnots(Diag, cout);
        KnotFloerForAlternatingKnots(Diag, outFile); 
      }
      else { 
        MorseCode M = Diag.GetSmallGirthMorseCode();
	outFile2<< "Knot # "<<t+1<< "\n\n";
        M.Print(outFile2);     
        if (M.GetGirth() > 2*MAXBRIDGE) {
	  cout<<" Knot # "<<t+1<< "seems to have Girth number > "<< 2*MAXBRIDGE<< "\n";
	  outFile<<" Knot # "<<t+1<< "seems to have Girth number > "<< 2*MAXBRIDGE<< "\n"; 
          cout.flush();
          outFile.flush(); outFile.close(); exit(1);
	}
	  
        KnotFloerComplex KFC=ComputingKnotFloer(M, Prime);
         
	ReportKnotFloerRanks(KFC, cout);
        ReportKnotFloerRanks(KFC, outFile);
          
        cout<<"Total rank : "<<KFC.Generators.size()<<endl;
	outFile<<"Total rank : "<<KFC.Generators.size()<<endl;
	  
        int genus=Genus(KFC);
        cout<<"Seifert genus : "<<genus<<endl;
	outFile<<"Seifert genus : "<<genus<<endl;
          
        cout<<"Fibered : ";
        outFile<<"Fibered : ";
        if(Fibered(KFC)) {cout<<"Yes\n"; outFile<<"Yes\n";}
        else {cout<<"No\n"; outFile<<"No\n";}                    
         
        cout<<"L-space knot : ";
        outFile<<"L-space knot : ";
        if(LSpaceKnot(KFC)) {cout<<"Yes\n"; outFile<<"Yes\n";}
        else {cout<<"No\n"; outFile<<"No\n";}
	  
        int tau=Tau(KFC);
        cout<<"Tau : "<<tau<<endl;
        outFile<<"Tau : "<<tau<<endl;
          
        int nu=Nu(KFC);
        cout<<"Nu  : "<<nu<<endl;
        outFile<<"Nu  : "<<nu<<endl;
          
	int epsilon=Epsilon(KFC);
        cout<<"Epsilon : "<<epsilon<<"\n";
        outFile<<"Epsilon : "<<epsilon<<"\n\n\n"; 
      } 
    }
      
    outFile.close();
    outFile2.close();
    
    clock_t endT=clock();
    double ElapsedT=double(endT-beginT)/CLOCKS_PER_SEC ;
    cout<<"\nFinished in "<<ElapsedT<<" seconds"<< endl;
    cout<<"Saved to "<<OutPutFileName<<endl;
    cout.flush();
    return 0;
}









 

 


 
 



