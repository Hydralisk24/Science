#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <new>
#include <vector>
#include <algorithm>

#include <fstream>
#include <string>

#include <sstream> // double -> string-hez kell
#include <iomanip> // setprecision

#define Msol 1.989e33		/*Solar mass (g)*/

#include "Fejlec.cpp"


using namespace std;


//space-hez felbontas a gridhez, 200 jo kb
// set in MCMC parameter file
int nn=150;
//int whole=1;
//int burnin=0;







				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				/* Allithato parameterek / setable oaratemeters */
				int n=250;	 	// hany legjobb legyen kiirva? / how many best write out? 
				int m=50;		// scriptbe hany keruljon, ami abrazolva is lesz? / how many of the be plotted?








//Ha nem egyenlo a 2 regio, akkor az abrazlasnal nem egymas mellett van, ami baj mert arra van allitva (Eth skalaja nincs feliratkozva az a baj)
// from MCMC
/*
double tkezd=0;
double tmin=0;
double tmax=300;
double R0min=2e12;   	  // Acceptable minimum
double Mmin=3*Msol;   	 
double Ekmin=0.1e51;   	 
double Emin=0.1e51; 	
double R0max=8e13;   	  // Acceptable maximum
double Mmax=60*Msol;   	 
double Ekmax=30e51;   	 
double Emax=30e51; 	

double ER0min=3e65;
double ER0max=2e12;
double Epmin=1e48;
double Epmax=3e51;
double Mnimin=5e-4*Msol;  
double Mnimax=0.5*Msol;
double tpmin=864;
double tpmax=8640000;
*/

double norma=1;  // chi**2 normalas (amivel leosztunk, szum 1/hiba**2) (jelenleg unused)
//double delta=0;



/*
//stringbol kimasolja az a es b kozott szoveget
string copy(string str, int a, int b){
	   
	   if(a>b) return "";
	   
	   string c;
	   int i;
	   
	   for(i=a;i<=b;i++)
	   c=c+str[i];
	   
	   return c;
	   }
*/
	   
	   
//konvertal szamba	   
double conv(string a){
	   double temp;
	     
	   temp=atof(a.c_str());

	   return temp;
	   }
	   
	   
	   
// Misc +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   
	   
// square	  
double h(double x){
	  return x*x;
	  }
	  
	  
	  
// nem 0
double nen(double x){
	  if(x==0) return 1;
	  return x;
	  }
	  
	   
	   
// linear interpolation
double interpol(double x, double xm, double ym, double xp, double yp){
	  double y;
	  
	  y=ym+ (x-xm) * (yp-ym) / (xp-xm);
	  
	  return y;
	  }
	  
	  
	  
// gives back the array number
int inverztar(double x, double a, double b){
	double y;
	
	//tar[i][0]=a + i *(b-a)/double(nn);
	y = (x - a) * double(nn)/(b-a);
	
	return int(y);
    }
    
    
    
int randsor(int i, int n){
	//if(i==0) return 0;
	//if(i==1) return n-1;
	
	return int( 1.*n*rand()/RAND_MAX );
}
    
    
    
    
    
    
    
    
    
    
    
    
/*    
    
// c++ data reading
std::ifstream f("parametersMC.inp"); 
double read(){

	double y=-1000;
	string str;
	int i;
	int j;
	//c++ fajl beolvasas, a sort olvassa be 
	do{  std::getline(f, str);  } while(str[0]=='#');
	

        
		j=0;	
		while(str[j]==' ') j++;
		i=j;
		while(str[i]!=' ' && i<=str.length()) i++;
		
		y = atof(   copy(str,j,i-1)    .c_str());
		
 
    
    return y;
}



void data()			 //Reading input file
 {


 if( !f.is_open() ){printf("No MCMC parameter file! Using default values!\n"); return;}

   
   tkezd=read();
   tmin=read();
   read();
   read();
   read();
   read();
   read();
   read();
   tmax=read();
   read();
   read();
   read();
   
   R0min=read();
   R0max=read();
   Mmin=read();  
   Mmax=read(); 	 
   Ekmin=read();  
   Ekmax=read(); 	 
   Emin=read(); 
   Emax=read(); 
   
   //kappamin=read(); 
   //kappamax=read();	
   //smin=read(); 
   //smax=read(); 	  
   //Tionmin=read();
   //Tionmax=read();
   //Mnimin=read();
   //Mnimax=read();
   //Agmin=read();
   //Agmax=read();
   
   read();
   read();
   read();
   read();
   read();	
   read();
   read();
   read();
   read();
   read();
   read();
   read();
   read();
   read();
   read();
   delta=read();
   read();	
   read();
   
   whole=int(read());
   nn=int(read());
   burnin=int(read());
   //Na=int(read());
   
   
   
 if(R0min<=0  || R0max<=R0min)    { R0min=2e12;   R0max=80e12; }
 if(Mmin<=0   || Mmax<=Mmin)      { Mmin=3;   Mmax=60; }
 if(Ekmin<=0  || Ekmax<=Ekmin)    { Ekmin=0.1;   Ekmax=30; }
 if(Emin<=0   || Emax<=Emin)      { Emin=0.1;   Emax=30; }
 if(nn<=0) nn=80;

   
   
 Mmin=Msol*Mmin;
 Mmax=Msol*Mmax;
 Emin=1e51*Emin;
 Emax=1e51*Emax;
 Ekmin=1e51*Ekmin;
 Ekmax=1e51*Ekmax;
 
	   
   }
	   
*/	   
	   
	   
	   
// likelihood function	   
double jof(double x){
	if(x<0) return 0;
	
	double e;
	//x=x-108;
	
	//e=1;
	//e=exp(-x*x);
	//e=1./(x*x);
	//e=exp(-0.5/2.5*1e6*x*x);
	//e=exp(-0.5*norma*x*x);
	
	//e=exp(-0.5*x);
	e=1;       // autocorr!!!!!!
	
	
	return e;
}



	   





















int main (){
	
	srand(time(NULL));

int i;  // cikluas valt
int j;
int k;
int q; 
int qq;
int qqq;


int baj=0;   // kapcsolok
int mod=0;
    


double chimin=0;
double chimax=0;


double N=0;
double N2=0;
double N3[2];
//double Nke[5];
double jo;
double szumkonf=0;

vector <double> Adatk;
//vector <double> Adat0;
//vector <double> Adat1;
//vector <double> Adat2;
//vector <double> Adat3;
//vector <double> Adat4;

//vector <vector <double> > Adat(1000,vector<double>(1000,3.141));
//vector <vector <double> > Adat;
vector <double> Adat[6];



   				
   				
   				
   				
   				
   				
//whole=0;    // teljes szoras?   				

//mod=0;		// unused


// misc
chimin=0;
chimax=1e100;


if(whole>1) whole=0;
if(whole<0) whole=0;
if(chimax<=0 || chimax<chimin) chimax=1e100;
data();
nn=grid;
std::ifstream f("rand.out");  
std::ifstream f2("rand2.out");  


std::ifstream f3("rand.out");  










//space ++++++
FILE *fkik;
FILE *fkik2;
FILE *fkik3;
FILE *fkika;
ofstream gnus;
ofstream gnusa;

if(whole==1){
fkik=fopen("SpaceBestW.txt","w");
fkik2=fopen("SpaceBestW2.txt","w");
fkik3=fopen("SpaceMeanW.txt","w");
fkika=fopen("SpaceALLW.txt","w");
gnus.open ("spaceW.gp");
//gnusa.open ("spaceAllW.gp");
}
else{
fkik=fopen("SpaceBestT.txt","w");
fkik2=fopen("SpaceBestT2.txt","w");
fkik3=fopen("SpaceMeanT.txt","w");
fkika=fopen("SpaceALLT.txt","w");
gnus.open ("spaceT.gp");
//gnusa.open ("spaceAllT.gp");
}



// beolvasok
double be1;
double be2;
double be3;
double be4;
double be5;
double be6;
double be7;

// tomb intek
int tari;
int tarik;
int tari1;
int tarik1;


double tarill[nn][2];
double tarka[7][nn];
double tarkak[8][nn];







double X1;
double X2;
double X0s[8];
double X1s[8][2];
double X2s[8][2];

double a1;
double b1;
double a2;
double b2;
double a3;
double b3;
double a4;
double b4;
double a5;
double b5;
double a6;
double b6;
double a7;
double b7;

double aa[8];
double bb[8];



// hatarok
a1=log10(Mmin/Msol);  b1=log10(Mmax/Msol);
a2=log10(Ekmin/1e51); b2=log10(Ekmax/1e51);

a3=log10(Mnimin/Msol);  b3=log10(Mnimax/Msol);
a4=log10(Epmin/1e51);   b4=log10(Epmax/1e51);

//a3=log10(R0min/1e12); b3=log10(R0max/1e12);
//a4=log10(Emin/1e51);  b4=log10(Emax/1e51);

a5=8; b5=9;  //10000 km/s -> cm/s 
a6=log10(ER0min/1); b6=log10(ER0max/1);

a7=log10(tpmin/86400);   b7=log10(tpmax/86400);



aa[1]=a1;
aa[2]=a2;
aa[3]=a3;
aa[4]=a4;
aa[5]=a5;
aa[6]=a6;
aa[7]=a7;

bb[1]=b1;
bb[2]=b2;
bb[3]=b3;
bb[4]=b4;
bb[5]=b5;
bb[6]=b6;
bb[7]=b7;







//double tar[nn][nn];
//double tar2[nn][nn];

//double * tar  = new double[nn];
//double * tar2 = new double[nn];

// ez mind foglalas
double** tar = new double*[nn];
double** tar2 = new double*[nn];
double** tara = new double*[nn];
double** tara2 = new double*[nn];
double** tarb = new double*[nn];
double** tarb2 = new double*[nn];
    double *sorb;

    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tar[j]=sorb;
    }
    


    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tar2[j]=sorb;
    }
    
    
    
    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tara[j]=sorb;
    }
    


    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tara2[j]=sorb;
    }
    
    
    
    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tarb[j]=sorb;
    }
    


    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    tarb2[j]=sorb;
    }
    


// mindent tarolo, =1+2+3+... (szumma n-1 ig)  -2
double tarm[4][nn][nn];
double tarmn[4][nn][nn];

/*
double*** tarm = new double**[nn];
double*** tarmn = new double**[nn];
double **sorbb;

  for(int k=0;k<5;k++){  
    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    sorbb[j]=sorb;
    }
  
  tarm[k]=sorbb;
    
  }
  
  
  
  for(int k=0;k<5;k++){  
    for(int j=0;j<nn;j++)
    {
    sorb= new double[nn];
    for(int i=0;i<nn;i++)
        sorb[i]=(double)j+i;

    sorbb[j]=sorb;
    }
  
  tarmn[k]=sorbb;
    
  }
*/





// NULLAZAS
// 2d tar nullazas
for(i=0;i<nn;i++) for(j=0;j<nn;j++){

	tar[i][j]=1e100;
	tar2[i][j]=1e100;
	
	tara[i][j]=0;
	tara2[i][j]=0;
	
	tarb[i][j]=0;
	tarb2[i][j]=0;
	
	
	tarm[0][i][j]=0;
	tarmn[0][i][j]=0;
	
	tarm[1][i][j]=0;
	tarmn[1][i][j]=0;
	
	tarm[2][i][j]=0;
	tarmn[2][i][j]=0;
	
	tarm[3][i][j]=0;
	tarmn[3][i][j]=0;
		
		
}



// 1 d
for(i=0;i<nn;i++){
	
	
	tarill[i][1]=1e100;
	
	tarka[1][i]=0;
	tarka[2][i]=0;
	tarka[3][i]=0;
	tarka[4][i]=0;
	tarka[5][i]=0;
	tarka[6][i]=0;
	tarka[7][i]=0;
	
	tarkak[1][i]=0;
	tarkak[2][i]=0;
	tarkak[3][i]=0;
	tarkak[4][i]=0;
	tarkak[5][i]=0;
	tarkak[6][i]=0;
	tarkak[7][i]=0;
	
}










  
//best ++++++ 
ofstream g;
ofstream g2;
ofstream gbat;
FILE *gsh;
FILE *gkonf;
ofstream gnub;
if(whole==1){
g.open ("bestW.txt"); 
g2.open ("brandW.txt"); 
gbat.open ("bestbatW.bat"); 
gsh=fopen("bestscriptW.sh","w");
//gkonf=fopen("KonfW.txt","w");
gnub.open ("bestW.gp");
}
else{
g.open ("bestT.txt"); 
g2.open ("brandT.txt"); 
gbat.open ("bestbatT.bat"); 
gsh=fopen("bestscriptT.sh","w");
//gkonf=fopen("KonfT.txt","w");
gnub.open ("bestT.gp");
}

string str;
string savestr[n];
string savestr2[n];

double has; 
double has2;
double has2best;
double savehas[n];

double conf95[2];
double conf67[2];
double conf95a;
double conf67a;
conf95[0]=1e100; conf95[1]=1e100;
conf67[0]=1e100; conf67[1]=1e100;

for(i=0;i<n;i++) savehas[i]=1e100+i*1e95;






// mean ++++++++++
double szum[5];
double szumk[5][5];
double be[5];;

double atl[5];
double szoras[5];

double korr;
    
for(k=0;k<5;k++)  szum[k]=0;
    
    
    
    
    
    
    
    
    
	// BEOLVASO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    

/* 
Hatulrol indul:
k=0 a legutolso
k=utolso szama addig  el kell vinni


16 ID 	 
15 ER0 [cm]    
14 M [M_sol]    
13 Ek [erg]    
12 Ep [erg]     
11 tp []     
10 v [km/s]       
9 T0 [K]      
8 M_Ni [M_sol]  
7 Ag [day^2]     
6 kappa [g/cm^2]    
5 s    
4 Tion [K]      
scatter [mag], chi^2; whole and without tail
3 scat w
2 scat t
1 chi w
0 chi t

   
 */  

    
qqq=0;    
cout << "Beolvasas, readin:\n" ;
q=0;    N=0;	N2=0;  qq=0;
//c++ fajl beolvasas, a sort olvassa be 
for(int cikl=0;cikl<    3 /*set to 1 for faster run*/   ;cikl++){
	q=0;
//while (std::getline(cikl==1?f2:f, str))
while (std::getline(      cikl==1?f2:( cikl==0?f:f3 )    , str))
    if(str[0]!='#') { if(q>burnin){
    	
     baj=0;
     
     
	 //read   
     i=str.length();   
	 for(k=0;k<=     16      ;k++){   //k !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 
	   while(str[i]==' '){i--;}
	   j=i;
	   while(str[i]!=' ' && i>0){i--;}
	   //if(i<1 && k!=14){ baj=1; break; cout << q << "  baj van!\n"; }
	
	
	   if(k==whole){
	     has=conv( copy(str, i, j) );
	     if(has*0!=0) has=1e100;
	     if(has<=0) has=1e100;
       }
       
       
       if(k==whole+2){
	     has2=conv( copy(str, i, j) );
	     if(has2*0!=0) has=1e100;
	     if(has2<=0) has=1e100;
       }
    
    
    
     if(cikl<2){
	
       //k !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       //Ep
       if(k==12){
	   be[4]=be4=conv( copy(str, i, j) );
       }
    
       //Ekin
       if(k==13){
	   be[2]=be2=conv( copy(str, i, j) );
       }
    
       //Mej
       if(k==14){
	   be[1]=be1=conv( copy(str, i, j) );
       } 
    
       //Mni
       if(k==8){
	   be[3]=be3=conv( copy(str, i, j) );
       }
       
       //v
       if(k==10){
	   be5=conv( copy(str, i, j) );
	   }
	   
	   //ER0
       if(k==15){
	   be6=conv( copy(str, i, j) );
	   }
	   
	   //tp
       if(k==11){
	   be7=conv( copy(str, i, j) );
       }
       
     }
    
	
    }
    
    if(baj==1) has=1e100;
    
    if(has<chimin) has=1e100;
    if(has>chimax) has=1e100;

	
	jo=jof(has); 
	if(cikl==0){ N=N+jo;  N2=N2+1; }
	
	
	if(cikl==0){
	// save
	//Adat.push_back(jo); // ki kell kommentelni
	Adat[0].push_back(jo );
	Adat[1].push_back(log10(be1) );
	Adat[2].push_back(log10(be2) );
	Adat[3].push_back(log10(be3) );
	Adat[4].push_back(log10(be4) );
	if(1) Adat[5].push_back(has );
	}
	
	
	//cout << has <<endl; getchar();
	
	
	
	
	
													
	//space
	if(cikl<2){
	// M & Ek
	tari =inverztar( log10(be1) , a1, b1);    
    tarik=inverztar( log10(be2/1e51) , a2, b2);              
    if(tari>nn-1)  tari=nn-1;   if(tari<0)  tari=0;
    if(tarik>nn-1) tarik=nn-1;  if(tarik<0) tarik=0;
	tari1 =tari;
	tarik1=tarik;
		
	
	  // nn*nn		   
	  if(has2<tar[tari][tarik]){ tar[tari][tarik]=has2;   }  // legjobb khi
	}
	  		
	
	  if(cikl==0){
	  tara[tari][tarik]+=jo;  // szum likelihood
	  tarb[tari][tarik]+=1;   // szum N
	
	  // marginalizacio
	  tarka[1][tari]+=jo;   tarkak[1][tari]+=1;
	  tarka[2][tarik]+=jo;  tarkak[2][tarik]+=1;
	  }
	
	
	if(cikl<2){
	// Mni & Ep
	tari =inverztar( log10(be3/1) , a3, b3);    
    tarik=inverztar( log10(be4/1e51) , a4, b4);               
    if(tari>nn-1) tari=nn-1;  if(tari<0) tari=0;
    if(tarik>nn-1) tarik=nn-1;  if(tarik<0) tarik=0;	
	
	  // nn*nn		   
	  if(has2<tar2[tari][tarik]){ tar2[tari][tarik]=has2;   }
	}
	  
	  
	  if(cikl==0){
	  tara2[tari][tarik]+=jo;
	  tarb2[tari][tarik]+=1;
	
	  // marginalizacio
	  tarka[3][tari]+=jo;   tarkak[3][tari]+=1;
	  tarka[4][tarik]+=jo;  tarkak[4][tarik]+=1;
	  
	  if(has<tarill[tarik][1]){ tarill[tarik][1]=has;  tarill[tarik][0]=be3;  } // spaceEgyenes-hez: corr illesztes
	  }
	  
	  
	  
	  
	

	
	// Tobbi tag, minden mindennel, korrelacios matrix
	if(cikl==0){
	 //1:tari1:M 
	 //2:tarik1:Ekin 
	 //3:tari:R0 
	 //4:tarik:Eth 
	
	 //1:2  tar !!!
	 tarm[0][tari1][tari  ]+=jo; tarmn[0][tari1][tari  ]+=1; //1:3
	 tarm[1][tari1][tarik ]+=jo; tarmn[1][tari1][tarik ]+=1; //1:4
	
	 tarm[2][tarik1][tari ]+=jo; tarmn[2][tarik1][tari ]+=1; //2:3
	 tarm[3][tarik1][tarik]+=jo; tarmn[3][tarik1][tarik]+=1; //2:4
	
	 //3:4  tar2 !!!
    
    
    
    
	
	




	
	 //best
	 
	 if(has<savehas[n-1])
	 for(i=0;i<n;i++) if(has<savehas[i]){ 
	 if(i==0){ X0s[1]=be1; X0s[2]=be2; X0s[3]=be3; X0s[4]=be4; X0s[5]=be5*1e5; X0s[6]=be6; X0s[7]=be7; has2best=has2; }
	 for(j=n-2;j>=i;j--){	savehas[j+1]=savehas[j]; savestr[j+1]=savestr[j]; }    
	 savehas[i]=has; savestr[i]=str;
	 break;
	
	 }
	 
	 
	
	
	 // mean
	 for(k=1;k<5;k++)  szum[k]=szum[k]+log10(be[k])   *jo;
	 // + save
	 
	 
	}
	
	
	if(cikl==2)
	 if(has<conf95a)
	 if( 1.*rand()/RAND_MAX < 1.*n/N){ 
	 /*savehas[qqq]=has;*/ if(qqq<n){ savestr2[qqq]=str; qqq++;}
	 //cout << qqq << "\n";
	 }
	
	
	// v & Eth*R0
	if(cikl==0){
	tari =inverztar( log10(be5*1e5) , a5, b5);    
    tarik=inverztar( log10(be6) , a6, b6);              
    if(tari>nn-1)  tari=nn-1;   if(tari<0)  tari=0;
    if(tarik>nn-1) tarik=nn-1;  if(tarik<0) tarik=0;
	
	  // marginalizacio
	  tarka[5][tari]+=jo;   tarkak[5][tari]+=1;
	  tarka[6][tarik]+=jo;  tarkak[6][tarik]+=1;
	  }
	  
	  
	
	// tp
	if(cikl==0){
	tari =inverztar( log10(be7) , a7, b7);    
    //tarik=inverztar( log10(be8/1) , a8, b8);              
    if(tari>nn-1)  tari=nn-1;   if(tari<0)  tari=0;
    //if(tarik>nn-1) tarik=nn-1;  if(tarik<0) tarik=0;
	
	  // marginalizacio
	  tarka[7][tari]+=jo;   tarkak[7][tari]+=1;
	  //tarka[8][tarik]+=jo;  tarkak[8][tarik]+=1;
	  }
	
	
	
	
	
	
        
    //cout << str << "  " << str.length() << endl;
	}
    q++; if(q%50000==0) printf("Line:%d \t %d. time\n",q,cikl+1) ; 
	//if(q>190584){ cout << q << "  " << str <<  endl; getchar(); }
	//q++; 
	
        
          
  } 
   
   
	
	
	
	
	
	
	  
    else {    // norma-t olvassa: szoras es khinegyzet kozotti
    	//csak a legelejen mukodik!!!
    	if(qq==1+whole){
    		
     i=str.length();    
	 while(str[i]==' '){i--;}
	 j=i;
	 while(str[i]!=' ' && i>0){i--;}
	
	
	 has=conv( copy(str, i, j) );
	 if(has*0!=0) has=1;
	 if(has<=0) has=1;
	 
	 norma=has;
	 
     }
    		
    	qq++;
	}
	
	
	
if(cikl==0){
sort(Adat[5].begin(), Adat[5].end());	
conf95a=Adat[5][ int(0.95*Adat[5].size()) ];
conf67a=Adat[5][ int(0.667*Adat[5].size()) ];
}

printf("REread: \n") ;
/*if(cikl<0) q=0;*/  


//r.rewind();


}
	
cout << "readin done\n" ;
	
	
	
	



	
	
	
	
	
	
	
// mean calculate +++++++++++++++++++++++++++++++++	
	
// atlag
for(k=1;k<5;k++)  atl[k]=szum[k]/N;    

//nullazas
for(k=0;k<5;k++)  szum[k]=0;
for(i=0;i<5;i++) for(k=0;k<5;k++) szumk[i][k]=0;

// std
for(j=0;j<N2;j++) for(k=1;k<5;k++)  szum[k]=szum[k]+h( Adat[k][j]-atl[k] )*Adat[0][j]; 
// correlation
for(j=0;j<N2;j++) for(i=1;i<5;i++) for(k=1;k<5;k++) szumk[i][k]=szumk[i][k] + ( Adat[i][j]-atl[i] ) * ( Adat[k][j]-atl[k] ) *Adat[0][j];
	
// norma
for(k=1;k<5;k++)  szoras[k]=sqrt(szum[k]/N);
for(i=1;i<5;i++) for(k=1;k<5;k++) szumk[i][k]=szumk[i][k]/N/szoras[i]/szoras[k];









    



// confidencia 95% and 67% calculate ++++++++++  
cout << endl << "sum likelihood: " << N << endl << endl;
for(int q=0;q<2;q++){
N3[q]=0;
for(i=0;i<nn;i++) for(j=0;j<nn;j++){
   //Adatk.push_back(tara[i][j]/  nen(tarb[i][j])  );
   //N3=N3+tara[i][j]/nen(tarb[i][j]);
   if(q==0) Adatk.push_back(tara[i][j] );
   if(q==0) N3[q]=N3[q]+tara[i][j];
   if(q==1) Adatk.push_back(tara2[i][j] );
   if(q==1) N3[q]=N3[q]+tara2[i][j];
   }
    


sort(Adatk.begin(), Adatk.end());

	
	//for(i=0;i<1e2;i++) cout << Adatk[nn*nn-1-i]/N << endl; //getchar();
	//for(i=0;i<1e2;i++) cout << Adatk[i] << endl; getchar();
	
	has=0; k=0; j=0; i=0;
	//while(has<0.995 && i<q-1 ){ has=has+Adatk[q-1-i]/N; i++;   if(has>0.67 && k==0){ conf67=Adatk[q-1-i]/N; k=1;}  if(has>0.95 && j==0){ conf95=Adatk[q-1-i]/N; j=1;} }
	while(has<0.995 && i<nn*nn-1 ){ 
									has=has+Adatk[nn*nn-1-i]/N3[q];    
									if(has>0.67 && k==0){ conf67[q]=Adatk[nn*nn-1-i]/N3[q]; k=1;}  
									if(has>0.95 && j==0){ conf95[q]=Adatk[nn*nn-1-i]/N3[q]; j=1;} 
									i++; }
	//conf995=Adatk[q-1-i]/N;
	//cout << has << endl;
    


/*
printf("   normalt ,   valodi\n");
printf("67\% = %e ,   %e\n",conf67[q], conf67[q]*N3[q]); //getchar();
printf("95\% = %e ,   %e\n",conf95[q], conf95[q]*N3[q]); //getchar();
*/

Adatk.clear();
    
}
    
    // set contour base 
	// set cntrparam levels discrete 0.95, 1
	// unset surface
















	//Konfidencia and best Calculate
	
// gkonf=fopen("KonfW.txt","w");


//X0s[2]=X0s[2]/1e51;
//X0s[3]=X0s[3]/1e12;
//X0s[4]=X0s[4]/1e51;


cout << endl;

cout << "Best values:\n";
cout << "1:Mej  2:Ekin  3:Mni  4:Ep  5:v  6:R0*Eth  7:tp\n";
for(j=1;j<8;j++) cout << j << "  " << X0s[j] << "  " << endl;
cout << endl;



/*
// 2 szigma
//be1 es be2 tempek!!!!
for(j=1;j<7;j++){

  be2=1e100;
  for(be1=0.002;be1<0.05;be1=be1+0.001){
    szumkonf=0;
    for(i=0;i<nn;i++){ szumkonf+=tarka[j][i]/N; if(szumkonf>=be1){ X1 = pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]) ); break; }}

    szumkonf=0;
    for(i=nn-1;i>=0;i--){ szumkonf+=tarka[j][i]/N; if(szumkonf>=0.05-be1){ X2 = pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]) ); break; }}

    if(X2-X1<be2){ be2=X2-X1; X1s[j][0]=X2; X2s[j][0]=X1; } 
  }
}


cout << "2 sigma\n";
for(j=1;j<7;j++) cout << j << "  " << X1s[j][0] << "  " << X2s[j][0] << endl;
cout << endl;



// 1 szigma
//be1 es be2 tempek!!!!
for(j=1;j<7;j++){

  be2=1e100;
  for(be1=0.002;be1<0.166;be1=be1+0.001){
    szumkonf=0;
    for(i=0;i<nn;i++){ szumkonf+=tarka[j][i]/N; if(szumkonf>=be1){ X1 = pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]) ); break; }}

    szumkonf=0;
    for(i=nn-1;i>=0;i--){ szumkonf+=tarka[j][i]/N; if(szumkonf>=0.166-be1){ X2 = pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]) ); break; }}

    if(X2-X1<be2){ be2=X2-X1; X1s[j][1]=X2; X2s[j][1]=X1; } 
  }
}


cout << "1 sigma\n";
for(j=1;j<7;j++) cout << j << "  " << X1s[j][1] << "  " << X2s[j][1] << endl;
cout << endl;
*/


/*for(i=0;i<nn;i++){
for(j=1;j<5;j++) fprintf(gkonf,"%e  %f    ",pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]) )  ,   tarka[j][i]/N   );
fprintf(gkonf,"\n");	
}


fclose(gkonf);*/



for(j=1;j<8;j++){

  vector <double> adattemp;
  double conf67;
  double conf95;
  double me[8];
  me[1]=1;
  me[2]=1e-51;
  me[3]=1;
  me[4]=1e-51;
  me[5]=1;
  me[6]=1;
  me[7]=1;
  
  szumkonf=0;
  for(i=0;i<nn;i++) adattemp.push_back(tarka[j][i]/N);
  sort(adattemp.begin(), adattemp.end());
  
  
    k=0; i=0;
  	while(szumkonf<0.95 && i<nn ){ 
									szumkonf=szumkonf+adattemp[nn-1-i];   
									if(szumkonf>0.67 && k==0){ conf67=adattemp[nn-1-i]; k=1;}  
									i++; }
									
    conf95=adattemp[nn-1-i];
    
    
    
    
    
    
    
    k=0;
    for(i=0;i<nn;i++){
    	
    	if(tarka[j][i]/N>conf95 && k==0){ X1s[j][0]=pow(10, aa[j]+ (i-1)/double(nn) *(bb[j]-aa[j])); k=1; }
    	if(tarka[j][i]/N>conf67        ){ X1s[j][1]=pow(10, aa[j]+ (i-1)/double(nn) *(bb[j]-aa[j])); break; }
    	if( pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]))  >X0s[j]*me[j]){ X2s[j][1]=pow(10, aa[j]+ (i-1)/double(nn) *(bb[j]-aa[j])); break; }
	} 
	
	
	k=0;
    for(i=nn-1;i>=0;i--){
    	
    	if(tarka[j][i]/N>conf95 && k==0){ X2s[j][0]=pow(10, aa[j]+ (i+1)/double(nn) *(bb[j]-aa[j])); k=1; }
    	if(tarka[j][i]/N>conf67        ){ X2s[j][1]=pow(10, aa[j]+ (i+1)/double(nn) *(bb[j]-aa[j])); break; }
    	if( pow(10, aa[j]+ i/double(nn) *(bb[j]-aa[j]))  <X0s[j]*me[j]){ X2s[j][1]=pow(10, aa[j]+ (i+1)/double(nn) *(bb[j]-aa[j])); break; }
	} 
  
  
  
  
  
  adattemp.clear();
}










// valtas
X1s[2][0]=X1s[2][0]*1e51;
X1s[3][0]=X1s[3][0]*1;
X1s[4][0]=X1s[4][0]*1e51;
X1s[6][0]=X1s[6][0]*1;

X2s[2][0]=X2s[2][0]*1e51;
X2s[3][0]=X2s[3][0]*1;
X2s[4][0]=X2s[4][0]*1e51;
X2s[6][0]=X2s[6][0]*1;

X1s[2][1]=X1s[2][1]*1e51;
X1s[3][1]=X1s[3][1]*1;
X1s[4][1]=X1s[4][1]*1e51;
X1s[6][1]=X1s[6][1]*1;

X2s[2][1]=X2s[2][1]*1e51;
X2s[3][1]=X2s[3][1]*1;
X2s[4][1]=X2s[4][1]*1e51;
X2s[6][1]=X2s[6][1]*1;





cout << "err: 2 sigma (also, felso ertek)\n";
for(j=1;j<8;j++) cout << j << "  " << X1s[j][0] << "  " << X2s[j][0] << endl;
cout << endl;

cout << "err: 1 sigma (also, felso ertek)\n";
for(j=1;j<8;j++) cout << j << "  " << X1s[j][1] << "  " << X2s[j][1] << endl;
cout << endl;
















// kiiras +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// Mean +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// kepernyore
printf("\n     : Mej   Ekin   Mni   Ep\n");
printf("Mean:  %f   %e   %e   %e\n",pow(10,atl[1]),pow(10,atl[2]),pow(10,atl[3]),pow(10,atl[4]));
printf("std.dev: %f   %f   %f   %f   (multiplier)\n",pow(10,szoras[1]),pow(10,szoras[2]),pow(10,szoras[3]),pow(10,szoras[4]));
printf("Lower: %f   %e   %e   %e\n",pow(10,atl[1])/pow(10,szoras[1]),pow(10,atl[2])/pow(10,szoras[2]),pow(10,atl[3])/pow(10,szoras[3]),pow(10,atl[4])/pow(10,szoras[4]));
printf("Upper: %f   %e   %e   %e\n",pow(10,atl[1])*pow(10,szoras[1]),pow(10,atl[2])*pow(10,szoras[2]),pow(10,atl[3])*pow(10,szoras[3]),pow(10,atl[4])*pow(10,szoras[4]));

printf("\nKorrelacios matrix:\n");
printf(" Mej   Ekin   Mni   Ep\n");
for(i=1;i<5;i++){ for(k=1;k<5;k++)
				printf("%f ",szumk[i][k]);
				printf("\n");
				}








// mean OUTPUT +++++++

FILE *gmean;
gmean=fopen("mean.log","a");



if(whole==1) fprintf(gmean,"\nWhole with tail\n"); else fprintf(gmean,"\nNo tail\n");
fprintf(gmean,"loop number: %d\n\n",q);

fprintf(gmean,"Confidence interval limit (in the plot) using %d grid (technical)\n",nn);
fprintf(gmean,"   normalised ,   actual\n");
fprintf(gmean,"Ek-M parameter space\n");
fprintf(gmean,"67\% = %e ,   %e\n",conf67[0], conf67[0]*N3[0]);
fprintf(gmean,"95\% = %e ,   %e\n\n",conf95[0], conf95[0]*N3[0]);
fprintf(gmean,"Ep-Mni  parameter space\n");
fprintf(gmean,"67\% = %e ,   %e\n",conf67[1], conf67[1]*N3[1]);
fprintf(gmean,"95\% = %e ,   %e\n\n\n",conf95[1], conf95[1]*N3[1]);

fprintf(gmean,"Best fit and Confidence intervals: (s: sigma)\n");
fprintf(gmean,"Mej [Msol]: \t%2.3f  \t\t -%2.3f  +%2.3f (2s)  \t\t\t -%2.3f  +%2.3f (1s)\n", X0s[1],X0s[1]-X1s[1][0],X2s[1][0]-X0s[1],X0s[1]-X1s[1][1],X2s[1][1]-X0s[1]);
fprintf(gmean,"Ekin [erg]: \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n",X0s[2],X0s[2]-X1s[2][0],X2s[2][0]-X0s[2],X0s[2]-X1s[2][1],X2s[2][1]-X0s[2]);
fprintf(gmean,"Mni [Msol]: \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n",  X0s[3],X0s[3]-X1s[3][0],X2s[3][0]-X0s[3],X0s[3]-X1s[3][1],X2s[3][1]-X0s[3]);
fprintf(gmean,"Ep [erg]: \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n", X0s[4],X0s[4]-X1s[4][0],X2s[4][0]-X0s[4],X0s[4]-X1s[4][1],X2s[4][1]-X0s[4]);
fprintf(gmean,"v [cm/s]: \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n",   X0s[5],X0s[5]-X1s[5][0],X2s[5][0]-X0s[5],X0s[5]-X1s[5][1],X2s[5][1]-X0s[5]);
fprintf(gmean,"R0*Eth : \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n",X0s[6],X0s[6]-X1s[6][0],X2s[6][0]-X0s[6],X0s[6]-X1s[6][1],X2s[6][1]-X0s[6]);
fprintf(gmean,"[cm*erg]\n");
fprintf(gmean,"tp [day] : \t%2.3e  \t -%2.3e  +%2.3e (2s)  \t -%2.3e  +%2.3e (1s)\n",X0s[7],X0s[7]-X1s[7][0],X2s[7][0]-X0s[7],X0s[7]-X1s[7][1],X2s[7][1]-X0s[7]);



fprintf(gmean,"\n     \t\t Mej      Ekin         Mni           Ep\n");
fprintf(gmean,"Mean:    \t %2.3f   %2.3e   %2.3e   %2.3e\n",pow(10,atl[1]),pow(10,atl[2]),pow(10,atl[3]),pow(10,atl[4]));
fprintf(gmean,"Std dev: \t %2.3f   %2.3f   %2.3f   %2.3f   (multiplier!)\n",pow(10,szoras[1]),pow(10,szoras[2]),pow(10,szoras[3]),pow(10,szoras[4]));
fprintf(gmean,"Lower limit: \t %2.3f   %2.3e   %2.3e   %2.3e\n",pow(10,atl[1])/pow(10,szoras[1]),pow(10,atl[2])/pow(10,szoras[2]),pow(10,atl[3])/pow(10,szoras[3]),pow(10,atl[4])/pow(10,szoras[4]));
fprintf(gmean,"Upper limit: \t %2.3f   %2.3e   %2.3e   %2.3e\n",pow(10,atl[1])*pow(10,szoras[1]),pow(10,atl[2])*pow(10,szoras[2]),pow(10,atl[3])*pow(10,szoras[3]),pow(10,atl[4])*pow(10,szoras[4]));

fprintf(gmean,"\n\nCorrelation matrix:\n");
fprintf(gmean,"Mej     Ekin     Mni     Ep\n");
for(i=1;i<5;i++){ for(k=1;k<5;k++)
				fprintf(gmean,"%1.3f   ",szumk[i][k]);
				fprintf(gmean,"\n");
				}
fprintf(gmean,"\n");





fclose(gmean);













// autocorr ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(0){
cout << "\n\nautocorr\n";

int Na=0;
FILE *autoc;
autoc=fopen("auto.txt","w");
Na=6000;
k=1;


double autocorr[Na];
for(j=0;j<Na;j++) autocorr[j]=0;
double Norma[2*Na];

//cout << "cica"; getchar();

//double Norma2[2*Na];  for(j=0;j<2*Na;j++) Norma2[j]=0;



//N2 all
for(j=0;j<Na;j++){

  for(i=0;i<N2-j;i++){
	
 	 autocorr[j]=autocorr[j]+ (Adat[k][i]-atl[k]) *  (Adat[k][i+j]-atl[k])   *Adat[0][i]*Adat[0][i+j];   
	 //autocorr[j]=autocorr[j]+1;   *Adat[0][i]
	
	
	
  }
 
 if(j%1000==0) cout << j << endl;
}


Norma[0]=autocorr[0];
Norma[2*Na-1]=autocorr[0];

for(j=1;j<Na;j++){   Norma[j]        = Norma[j-1]        - h((Adat[k][j]-atl[k])*Adat[0][j]) ;
                     Norma[2*Na-1-j] = Norma[2*Na-1-j+1] - h((Adat[k][N2-1-j]-atl[k])*Adat[0][N2-1-j]);            
					 //Norma[j]        = 1 ;
					 //Norma[2*Na-1-j] = 1 ;
					 }
                     

/* // Keplet alapjan, tesztelve, a masik ugyan ezt adja
for(j=0;j<Na;j++){
				for(i=0;i<N2-j;i++)     Norma2[j]        = Norma2[j] + h(Adat[k][i]-atl[k]) ;
				for(i=N2-1;i>=j;i--)  Norma2[2*Na-1-j] = Norma2[2*Na-1-j] + h(Adat[k][i]-atl[k]) ;
				if(j%1000==0) cout << j << endl;
				
				}
*/


//for(j=0;j<Na;j++) fprintf(autoc,"%d  %f \n",j,autocorr[j]/autocorr[0]);

	for(j=0;j<Na;j++) fprintf(autoc,"%d  %f \n",j,autocorr[j]/sqrt(Norma[j]*Norma[2*Na-1-j]) );



//test
//for(j=0;j<Na;j++) fprintf(autoc,"%d  %f \n",j,Norma2[j]*Norma2[2*Na-1-j]/(Norma[j]*Norma[2*Na-1-j]) );
	


//cout << "\n\nautocorrelacio:\n\n";
//for(j=0;j<Na;j++) cout << j << "  " << autocorr[j]/autocorr[0] << endl;



fclose(autoc);
cout << "autocorr done\n\n";
}
// autocorr ------------------------------------------------------------




















    

















	//space OUTPUT ++++++ !!!!!!!!!!!!!
	
//fkik=fopen("SpaceBestW.txt","w");
//fkik2=fopen("SpaceBestW2.txt","w");
//fkik3=fopen("SpaceMeanW.txt","w");
//fkika=fopen("SpaceALLW.txt","w");
if(whole==1){
fprintf(fkik ,"#Whole with tail:\n");
fprintf(fkik2,"#Whole with tail:\n");
fprintf(fkik3,"#Whole with tail:\n");
}
else{
fprintf(fkik ,"#Without tail:\n");
fprintf(fkik2,"#Without tail:\n");
fprintf(fkik3,"#Without tail:\n");
}

fprintf(fkik ,"#2D margalisation:\n");
fprintf(fkik2,"#2D margalisation:\n");
fprintf(fkik3,"#1D margalisation:\n");

fprintf(fkik ,"#grid:%d*%d\n",nn,nn);
fprintf(fkik2,"#grid:%d*%d\n",nn,nn);
fprintf(fkik3,"#grid:%d\n",nn);

fprintf(fkik ,"#Number of all points:%d\n",q-burnin);
fprintf(fkik2,"#Number of all points:%d\n",q-burnin);
fprintf(fkik3,"#Number of all points:%d\n",q-burnin);

fprintf(fkik ,"#Burnin:%d\n",burnin);
fprintf(fkik2,"#Burnin:%d\n",burnin);
fprintf(fkik3,"#Burnin:%d\n",burnin);


fprintf(fkik ,"#Mej [Msun],   Ekin [foe], best khi^2, posterior, N \n");
fprintf(fkik2,"#Mni [1e12 cm],   Ep [foe], best khi^2, posterior, N \n");
fprintf(fkik3,"#Mej [Msun] posterior N,   Ekin [foe] posterior N,   Mni [Msun] posterior N,   Ep [foe] posterior N,   v [cm/s] posterior N,   Eth*R0 [erg*cm] posterior N,   tp [day] posterior N\n");

for(i=0;i<nn;i++){

  //2D
  for(j=0;j<nn;j++){
   fprintf(fkik ,"%e %e %e %e %d\n",  pow(10, a1+ i/double(nn) *(b1-a1) ), pow(10, a2+ j/double(nn) *(b2-a2) ), tar[i][j] , tara[i][j] , int(tarb[i][j]) );
   fprintf(fkik2,"%e %e %e %e %d\n",  pow(10, a3+ i/double(nn) *(b3-a3) ), pow(10, a4+ j/double(nn) *(b4-a4) ), tar2[i][j] , tara2[i][j] , int(tarb2[i][j]) );
  }
 
  fprintf(fkik ,"\n");
  fprintf(fkik2,"\n");
  
  // 1D marginalized
  fprintf(fkik3,"%2.3e %2.6f %d    %2.3e %2.6f %d    ",  pow(10, a1+ i/double(nn) *(b1-a1) ), tarka[1][i]/N , int(tarkak[1][i]) , pow(10, a2+ i/double(nn) *(b2-a2) ), tarka[2][i]/N , int(tarkak[2][i]) );
  fprintf(fkik3,"%2.3e %2.6f %d    %2.3e %2.6f %d    ",  pow(10, a3+ i/double(nn) *(b3-a3) ), tarka[3][i]/N , int(tarkak[3][i]) , pow(10, a4+ i/double(nn) *(b4-a4) ), tarka[4][i]/N , int(tarkak[4][i]) );
  fprintf(fkik3,"%2.3e %2.6f %d    %2.3e %2.6f %d    ",  pow(10, a5+ i/double(nn) *(b5-a5) ), tarka[5][i]/N , int(tarkak[5][i]) , pow(10, a6+ i/double(nn) *(b6-a6) ), tarka[6][i]/N , int(tarkak[6][i]) );
  fprintf(fkik3,"%2.3e %2.6f %d\n",                      pow(10, a7+ i/double(nn) *(b7-a7) ), tarka[7][i]/N , int(tarkak[7][i]) );//, pow(10, a6+ i/double(nn) *(b6-a6) ), tarka[6][i]/N , int(tarkak[6][i]) );
}






/*
	//1:2  tar !!!
	tarm[0][tari1][tari  ]+=jo; //1:3
	tarm[1][tari1][tarik ]+=jo; //1:4
	
	tarm[2][tarik1][tari ]+=jo; //2:3
	tarm[3][tarik1][tarik]+=jo; //2:4
	
	//3:4  tar2 !!!
*/


// Minden mindennel, 2D                       ++++++++ BROKEN ++++++
if(whole==1){
fprintf(fkika ,"#Whole with tail:\n");
}
else{
fprintf(fkika ,"#Without tail:\n");

}
fprintf(fkika ,"#Minden mindennel, Everything with everything, grid: %d*%d\n",nn,nn);
fprintf(fkika,"#Number of all points:%d\n",q-burnin);
fprintf(fkika,"#Burnin:%d\n",burnin);
fprintf(fkika ,"#X value, Y value, posterior, N\n");
fprintf(fkika ,"#Mej:Ekin (1-4), Mej:R0 (5-8), Mej:Eth (9-12), Ekin:R0 (13-16), Ekin:Eth (17-20), R0:Eth (21-24)\n");
for(i=0;i<nn;i++){

  for(j=0;j<nn;j++){
    fprintf(fkika ,"%e %e %e %d\t",  pow(10, aa[1]+ i/double(nn) *(bb[1]-aa[1]) ), pow(10, aa[2]+ j/double(nn) *(bb[2]-aa[2]) ), tara[i][j] , int(tarb[i][j]) );
    fprintf(fkika ,"%e %e %e %d\t",  pow(10, aa[1]+ i/double(nn) *(bb[1]-aa[1]) ), pow(10, aa[3]+ j/double(nn) *(bb[3]-aa[3]) ), tarm[0][i][j] , int(tarmn[0][i][j]) );
    fprintf(fkika ,"%e %e %e %d\t",  pow(10, aa[1]+ i/double(nn) *(bb[1]-aa[1]) ), pow(10, aa[4]+ j/double(nn) *(bb[4]-aa[4]) ), tarm[1][i][j] , int(tarmn[1][i][j]) );

    fprintf(fkika ,"%e %e %e %d\t",  pow(10, aa[2]+ i/double(nn) *(bb[2]-aa[2]) ), pow(10, aa[3]+ j/double(nn) *(bb[3]-aa[3]) ), tarm[2][i][j] , int(tarmn[2][i][j]) );
    fprintf(fkika ,"%e %e %e %d\t",  pow(10, aa[2]+ i/double(nn) *(bb[2]-aa[2]) ), pow(10, aa[4]+ j/double(nn) *(bb[4]-aa[4]) ), tarm[3][i][j] , int(tarmn[3][i][j]) );

    fprintf(fkika ,"%e %e %e %d\n",  pow(10, aa[3]+ i/double(nn) *(bb[3]-aa[3]) ), pow(10, aa[4]+ j/double(nn) *(bb[4]-aa[4]) ), tara2[i][j] , int(tarb2[i][j]) );
  }
 
 
  fprintf(fkika ,"\n");

}



fclose(fkik);
fclose(fkik2);
fclose(fkik3);
fclose(fkika);





















	//best ++++++   

// OUTPUT !!!!!!!
//g.open ("bestW.txt"); 
//g2.open ("brandW.txt")
g << "#The best " << n << " fit parameters (lot of them repeats)." << endl;
g2 << "#Randomly chosen " << qqq << " accepted 0.95 level fit parameters." << endl;
//best: legjobb n kiirasa
//brand (best rand) 95% intervalumba tartozóból kiirni veletlenszeruen valamennyit
for(i=0;i<n;i++) g << savestr[i] << endl;
for(i=0;i<qqq;i++) g2 << savestr2[i] << endl;




// SCRIPT
//gbat.open ("bestbatW.bat"); 
//gsh=fopen("bestscriptW.sh","w");
//gnub.open ("bestW.gnu");

/*
Hiba: a veletlen kivalasztaskor n/N (kivalasztas/osszes) esellyel valaszt ki.
Ez azt jelenti, hogy nem fog n darabot kivalasztani, Kevesebbet. Tobbet nem mert maximalizalva van
Am n>>m, tehat jo esellyel kiir mindent. De neha ismetlodhet ugyan az.
Ures sorkor nem ir echo-t es a modellezo kod nem fut le
Igy ahol ureset valasztott ki ott semmi nem lesz, es hianyozni is fog az a kimenet
*/

// le kell forditattni elotte!
//gbat << "gcc LC3p4.cpp -o" << endl;
string strtemp;
for(i=0;i<m;i++){
			strtemp.clear();
			gbat << "del par.inp" << endl;
 			if(i==0) gbat << "echo " << savestr[0] << " > par.inp" << endl;
 			if(i!=0) strtemp=savestr2[randsor(i,qqq)];
			if(i!=0 && strtemp.length()>1) gbat << "echo " << strtemp << " > par.inp" << endl;
 			gbat << "LC3p4.exe 1" << endl;
 			if(whole==1){
			gbat << "del kimenetW" << i << ".out" << endl;
 			gbat << "ren kimenet.out kimenetW" << i << ".out" << endl;
 			} else{
 			gbat << "del kimenetT" << i << ".out" << endl;
 			gbat << "ren kimenet.out kimenetT" << i << ".out" << endl;
 			}
 			
 			
}
;



fprintf(gsh,"#rm LC3p4.run\n");
fprintf(gsh,"g++ LC3p4.cpp -o LC3p4.run\n");
for(i=0;i<m;i++){
 			if(i==0) fprintf(gsh,"echo %s > par.inp\n",savestr[0].c_str());
 			if(i!=0) fprintf(gsh,"echo %s > par.inp\n",savestr2[randsor(i,qqq)].c_str());
 			fprintf(gsh,"./LC3p4.run 1\n");
 			if(whole==1) fprintf(gsh,"mv kimenet.out kimenetW%d.out\n",i);
 			else fprintf(gsh,"mv kimenet.out kimenetT%d.out\n",i);
 			
 			
}
//fprintf(gsh,"rm a.out\n");
// dos2unix bestscript.sh
// leveszi aszemt karraktereket! de telepeni kell





gnub << "#GNUPLOT script for plotting the best fit, and " << m << " random 0.95 level fit" << endl << endl;

for(j=0;j<2;j++){

if(j==0) gnub << "set terminal jpeg  size 1024,768" << endl;
if(j==1) gnub << "set key noautotitles samplen -1" << endl;
if(j==1) gnub << "set term postscript eps enhanced color 16" << endl;

if(whole==1) gnub << "set output \"bestW." << (j==0?"jpeg":"eps") << "\"" << endl;
else         gnub << "set output \"bestT." << (j==0?"jpeg":"eps") << "\"" << endl;
gnub << "reset" << endl;
if(whole==1) gnub << "set title \"Best fits (whole)\" font \",20\"" << endl;
else         gnub << "set title \"Best fits (no tail)\" font \",20\"" << endl;
gnub << "set xlabel \"Phase [days]\"" << endl;
gnub << "set ylabel \"Log of luminosity [erg]\"" << endl;
gnub << "unset key" << endl;
gnub << "a=" << tkezd << endl << endl;

if(whole==1){
gnub << "p [0:" << tmax << "]" ;
for(i=1;i<m;i++) gnub << "   \"kimenetW" << i << ".out\"  every ::"<<tmin<<" u 1:3 w l lw 2 lt " << (j==0?"8":"5") << " title \"" << i+1 << "\",\\" << endl;
gnub  << " \"kimenetW0.out\"  every ::"<<tmin<<" u 1:3 w l lw 3 lt " << (j==0?"1":"1") << " title \"1\",\\" << endl;
}
else{
//gnub << "p [0:" << tmax << "]" << " \"kimenetT0.out\"  every ::"<<tmin<<" u 1:3 w l lw 4 title \"1\",\\" << endl;
gnub << "p [0:" << tmax << "]" ;
for(i=1;i<m;i++) gnub << "   \"kimenetT" << i << ".out\"  every ::"<<tmin<<" u 1:3 w l lw 2 lt " << (j==0?"8":"5") << " title \"" << i+1 << "\",\\" << endl;
gnub  << " \"kimenetT0.out\"  every ::"<<tmin<<" u 1:3 w l lw 3 lt " << (j==0?"1":"1") << " title \"1\",\\" << endl;
}


//gnub <<   "   \"bol-gorbe.txt\" u ($1-a):(-0.4*$2+17.5+18):(0.4*$3)  ps 1.5 pt 7 lt -1  w yerr title \"bol\""  << endl << endl << "set out" << endl << endl;
gnub <<   "   \"bol-gorbe.txt\" u ($1-a):(-0.4*$4+17.5+18):(0.4*$5)  ps 1.5 pt 7 lt -1  w yerr title \"bol\""  << endl << endl << endl << endl;
 
}

gnub << endl << endl << "set out" << endl << endl;
    
    
g.close();  
g2.close();  
gbat.close();
gnub.close();
fclose(gsh);  


























	//space SCRIPT ++++++ !!!!!!!!!!
	
// gnus.open ("spaceW.gnu");
// gnusa.open ("spaceAllW.gnu");

//Ha nem egyenlo a 2 regio, akkor az abrazlasnal nem egymas mellett van, ami baj mert arra van allitva (Eth skalaja nincs feliratkozva az a baj)
if(Ekmin!=Emin || Ekmax!=Emax) cout << "Warning Energy region not the same, plotting problems!" << endl;
has=(log(has2best)-0.1);

gnus << "#GNUPLOT script for plotting the parameter space for Mej-Ekin and R0-Eth." << endl;
if(whole==1) gnus << "name=\"Parameter space (whole)\"" << endl;
else         gnus << "name=\"Parameter space (no tail)\"" << endl;
gnus << "ref=\"\"" << endl;
gnus << "#ref: file name tag" << endl;
gnus << "#Number of all points: " << q-burnin << endl;
gnus << "#Burnin: " << burnin << endl;
gnus << "#grid: " << nn << "*" << nn << endl << endl;

// Figyelem! a ket contour elter, az abrazolas miatt a 2. (R0 Eth) nem lett normalva, 
// hogy a tengelyen 1 alatti szamok legyenek
gnus << "#Calculated values:" << endl;
gnus << "a=" << has << endl;
gnus << "#cont67a=" << conf67[0] << endl;
gnus << "#cont95a=" << conf95[0] << endl;
gnus << "cont67b=" << conf67[1] << endl;
gnus << "cont95b=" << conf95[1] << endl;
gnus << "cont67a=" << conf67[0]*N3[0] << endl;
gnus << "cont95a=" << conf95[0]*N3[0] << endl;
gnus << "#cont67b=" << conf67[1]*N3[1] << endl;
gnus << "#cont95b=" << conf95[1]*N3[1] << endl;
gnus << endl << endl;



for(j=0;j<2;j++)
for(i=0;i<3;i++){
	


//set output "best.eps"

if(j==0) gnus << "set terminal jpeg  size 1024,768" << endl;
if(j==1) gnus << "set key noautotitles samplen -1" << endl;
if(j==1) gnus << "set term postscript eps enhanced color 16" << endl;

if(i==0)
if(whole==1) gnus << "set output ref.\"spaceW." << (j==0?"jpeg":"eps") << "\"" << endl;
else         gnus << "set output ref.\"spaceT." << (j==0?"jpeg":"eps") << "\"" << endl;
if(i==1)
if(whole==1) gnus << "set output ref.\"spaceWl." << (j==0?"jpeg":"eps") << "\"" << endl;
else         gnus << "set output ref.\"spaceTl." << (j==0?"jpeg":"eps") << "\"" << endl;
if(i==2)
if(whole==1) gnus << "set output ref.\"spaceWn." << (j==0?"jpeg":"eps") << "\"" << endl;
else         gnus << "set output ref.\"spaceTn." << (j==0?"jpeg":"eps") << "\"" << endl;

gnus << "reset" << endl;

//if(whole==1) gnus << "set multiplot layout 1, 2  title \"Parameter space (whole)\"" << endl;
//else         gnus << "set multiplot layout 1, 2  title \"Parameter space (no tail)\"" << endl;
gnus << "set multiplot layout 1, 2  title name  font \", 20\"" << endl;
gnus << "set tmargin at screen 0.95" << endl;
gnus << "set bmargin at screen 0.11" << endl;
gnus << "set lmargin at screen 0.08" << endl;
gnus << "set rmargin at screen 0.44" << endl << endl;


gnus << "set view map" << endl;
gnus << "set logscale xy" << endl;
gnus << "unset colorbox" << endl;
gnus << "set xlabel \"Ejected mass [M_{sol}]\"  font \", 20\"" << endl;
gnus << "set ylabel \"Kinetic energy [10^{51} erg]\"  font \", 20\" offset 1.5" << endl;
//gnus << "set ylabel \"Kinetic energy [1e51 erg]\"  font \", 20\" offset 46.5" << endl;
gnus << "unset key" << endl;
gnus << "set xtics 5, 3" << endl;
gnus << "set ytics 0.4, 5" << endl;
gnus << "set tics font \", 15\"" << endl;
gnus << "set ytics  offset 0.5" << endl;
gnus << "set cbtics font \", 1\"" << endl;
if(i==0) gnus << "set palette model HSV rgbformulae 3,2,2" << endl;
if(i!=0) gnus << "set palette model CMY rgbformulae 7,5,15" << endl;

if(i==1){
gnus << "set contour base" << endl;
gnus << "set cntrparam levels discrete cont67a, cont95a" /*<< conf67[0] << ", " << conf95[0]*/ << endl;
//gnus << "set cntrparam levels discrete " << conf67[0]*N3[0] << ", " << conf95[0]*N3[0] << endl;
gnus << "set cntrparam bspline" << endl;
gnus << "set cntrparam points 7" << endl;
gnus << "set cntrparam order 10" << endl;
gnus << "#unset surface" << endl;
}
// surface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(i==0){
gnus << "set cbrange [ 0.0 : a ]" << endl;
if(whole==1)
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][a:0]  \"SpaceBestW.txt\"  u 1:2:(log($3))  with pm3d  lw 2" << endl  << endl;
else
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][a:0]  \"SpaceBestT.txt\"  u 1:2:(log($3))  with pm3d  lw 2" << endl  << endl;
}

if(i==1){
gnus << "set cbrange [  ]" << endl;
if(whole==1)
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][]  \"SpaceBestW.txt\"  u 1:2:4  with pm3d  lw 2" << endl  << endl;
else
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][]  \"SpaceBestT.txt\"  u 1:2:4  with pm3d  lw 2" << endl  << endl;
}


if(i==2){
gnus << "set cbrange [  ]" << endl;
if(whole==1)
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][]  \"SpaceBestW.txt\"  u 1:2:5  with pm3d  lw 2" << endl  << endl;
else
gnus << "splot [" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "][]  \"SpaceBestT.txt\"  u 1:2:5  with pm3d  lw 2" << endl  << endl;
}



gnus << "set lmargin at screen 0.52" << endl;
gnus << "set rmargin at screen 0.88" << endl;
gnus << "set view map" << endl;
gnus << "set colorbox" << endl;
gnus << "set logscale xy" << endl;
gnus << "set xlabel \"Nickel mass [M_{sol}]\"  font \", 20\"" << endl;
gnus << "set ylabel \"Magnetar energy [10^{51} erg]\"  font \", 20\" offset 1.5" << endl;
if(i==0) gnus << "set cblabel \"best log of chi_N^2\"  font \", 20\" offset -1.5" << endl;
if(i==1) gnus << "set cblabel \"sum likelihood\"  font \", 20\" offset -1.5" << endl;
if(i==2) gnus << "set cblabel \"sum points number\"  font \", 20\" offset -1.5" << endl;
gnus << "unset key" << endl;
gnus << "unset y2label " << endl;
gnus << "set xtics .001, 5" << endl;
gnus << "set tics font \", 15\"" << endl;
gnus << "set ytics 0.004, 5  offset 0.5" << endl;
if(i==0) gnus << "set palette model HSV rgbformulae 3,2,2" << endl;
if(i!=0) gnus << "set palette model CMY rgbformulae 7,5,15" << endl;

if(i==1){
//gnus << "set cntrparam levels discrete " << conf67[1] << ", " << conf95[1] << endl;
//gnus << "#set cntrparam levels discrete " << conf67[1]*N3[1] << ", " << conf95[1]*N3[1] << endl;
gnus << "set cntrparam levels discrete cont67b, cont95b" << endl;
}


if(i==0){
if(whole==1)
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][a:0]  \"SpaceBestW2.txt\"  u 1:2:(log($3))  with pm3d  lw 2" << endl;
else
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][a:0]  \"SpaceBestT2.txt\"  u 1:2:(log($3))  with pm3d  lw 2" << endl;
}

if(i==1){
if(whole==1)
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][]  \"SpaceBestW2.txt\"  u 1:2:($4/" <<N3[1]<< ")  with pm3d  lw 2" << endl;
else
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][]  \"SpaceBestT2.txt\"  u 1:2:($4/" <<N3[1]<< ")  with pm3d  lw 2" << endl;
}

if(i==2){
if(whole==1)
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][]  \"SpaceBestW2.txt\"  u 1:2:5  with pm3d  lw 2" << endl;
else
gnus << "splot [" << Mnimin/Msol << ":" << Mnimax/Msol << "][" << Epmin/1e51 << ":" << Epmax/1e51 << "][]  \"SpaceBestT2.txt\"  u 1:2:5  with pm3d  lw 2" << endl;
}


gnus << "unset multiplot"  << endl << "set out" << endl << endl << endl << endl << endl << endl;


} // i

gnus.close();



/*  UNUSED */
/*
FILE * fkik4=fopen("Spaceegyenes.txt","w");

fprintf(fkik4,"fit m*x+b \"Spaceegyenes.txt\"  u (log10($1)):(log10($2)):3 via m,b\n");
fprintf(fkik4,"plot   \"Spaceegyenes.txt\"  u  (log10($1)):(log10($2)) w l,\\ \n");
fprintf(fkik4,"  m*x+b\n\n");

for(i=0;i<nn;i++) fprintf(fkik4,"%e %e  %e\n",  pow(10, a4+ i/double(nn) *(b4-a4) ), tarill[i][0]/1e12 , tarill[i][1] );
fclose(fkik4);
*/

// m*log(Eth)+b = log(R0)  ezt illesztjuk
// m=-1 lesz
// 10**b = Eth * R0
// 10**b Eth es R0 mertek egysegeben lesz megadva, 
// vagyis pl b=1 : 10 foe*(12cm) = 10**(1+12+51) erg*cm



















// [a:b][c:d]
// "[%f:%f][%f:f%]",a,b,c,d
// "[" << Mmin/Msol << ":" << Mmax/Msol << "][" << Ekmin/1e51 << ":" << Ekmax/1e51 << "]"


/*
	//1:2  tar !!!
	tarm[0][tari1][tari  ]+=jo; //1:3
	tarm[1][tari1][tarik ]+=jo; //1:4
	
	tarm[2][tarik1][tari ]+=jo; //2:3
	tarm[3][tarik1][tarik]+=jo; //2:4
	
	//3:4  tar2 !!!
*/
/*
   R0min
   R0max
   Mmin 
   Mmax	 
   Ekmin 
   Ekmax 	 
   Emin
   Emax 
*/


// hatarok a plotnal +++
/*
double minmax[6][4];
double minmax2[4][2];


minmax[0][0]=Mmin/Msol; minmax[0][1]=Mmax/Msol;
minmax[1][0]=Mmin/Msol; minmax[1][1]=Mmax/Msol;
minmax[2][0]=Mmin/Msol; minmax[2][1]=Mmax/Msol;

minmax[3][0]=Ekmin/1e51; minmax[3][1]=Ekmax/1e51;
minmax[4][0]=Ekmin/1e51; minmax[4][1]=Ekmax/1e51;

minmax[5][0]=R0min/1e12; minmax[5][1]=R0max/1e12;




minmax[0][2]=Ekmin/1e51; minmax[0][3]=Ekmax/1e51;
minmax[1][2]=R0min/1e12; minmax[1][3]=R0max/1e12;
minmax[2][2]=Emin/1e51;  minmax[2][3]=Emax/1e51;

minmax[3][2]=R0min/1e12; minmax[3][3]=R0max/1e12;
minmax[4][2]=Emin/1e51;  minmax[4][3]=Emax/1e51;

minmax[5][2]=Emin/1e51;  minmax[5][3]=Emax/1e51;





minmax2[0][0]=Mmin/Msol;  minmax2[0][1]=Mmax/Msol;
minmax2[1][0]=Ekmin/1e51; minmax2[1][1]=Ekmax/1e51;
minmax2[2][0]=R0min/1e12; minmax2[2][1]=R0max/1e12;
minmax2[3][0]=Emin/1e51;  minmax2[3][1]=Emax/1e51;
*/





/* UNUSED & BROKEN    */
/*

string minmax[6][4];
string minmax2[4][2];


minmax[0][0]="Mmin"; minmax[0][1]="Mmax";
minmax[1][0]="Mmin"; minmax[1][1]="Mmax";
minmax[2][0]="Mmin"; minmax[2][1]="Mmax";

minmax[3][0]="Ekmin"; minmax[3][1]="Ekmax";
minmax[4][0]="Ekmin"; minmax[4][1]="Ekmax";

minmax[5][0]="R0min"; minmax[5][1]="R0max";




minmax[0][2]="Ekmin"; minmax[0][3]="Ekmax";
minmax[1][2]="R0min"; minmax[1][3]="R0max";
minmax[2][2]="Emin";  minmax[2][3]="Emax";

minmax[3][2]="R0min"; minmax[3][3]="R0max";
minmax[4][2]="Emin";  minmax[4][3]="Emax";

minmax[5][2]="Emin";  minmax[5][3]="Emax";





minmax2[0][0]="Mmin";  minmax2[0][1]="Mmax";
minmax2[1][0]="Ekmin"; minmax2[1][1]="Ekmax";
minmax2[2][0]="R0min"; minmax2[2][1]="R0max";
minmax2[3][0]="Emin";  minmax2[3][1]="Emax";












gnusa << "#GNUPLOT script for plotting the parameter space between all 4 main parameter" << endl;
gnusa << "#Number of all points: " << q-burnin << endl;
gnusa << "#Burnin: " << burnin << endl;
gnusa << "#grid: " << nn << "*" << nn << endl << endl;

gnusa << minmax[0][0] << "=" << Mmin/Msol << endl;
gnusa << minmax[0][1] << "=" << Mmax/Msol << endl;
gnusa << minmax[3][0] << "=" << Ekmin/1e51 << endl;
gnusa << minmax[3][1] << "=" << Ekmax/1e51 << endl;
gnusa << minmax[5][0] << "=" << R0min/1e12 << endl;
gnusa << minmax[5][1] << "=" << R0max/1e12 << endl;
gnusa << minmax[5][2] << "=" << Emin/1e51 << endl;
gnusa << minmax[5][3] << "=" << Emax/1e51 << endl;

gnusa << endl << endl;




	

gnusa << "set terminal jpeg  size 1024,768" << endl;
gnusa << "#set key noautotitles samplen -1" << endl;
gnusa << "#set term postscript eps enhanced color 16" << endl;


if(whole==1) gnusa << "set output \"spaceAllW.jpeg\"" << endl;
else         gnusa << "set output \"spaceAllT.jpeg\"" << endl;
if(whole==1) gnusa << "#set output \"spaceAllW.eps\"" << endl;
else         gnusa << "#set output \"spaceAllT.eps\"" << endl;

gnusa << "reset" << endl;
gnusa << "#posterior plot limit:" << endl;
gnusa << "b=0.5" << endl;                  // posterior hatar!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
gnusa << "c=0.2" << endl;  
if(whole==1) gnusa << "set multiplot layout 2, 2  title \"Correlation space (whole)\"" << endl;
else         gnusa << "set multiplot layout 2, 2  title \"Correlation space (no tail)\"" << endl;



gnusa << "set view map" << endl;
gnusa << "set logscale x" << endl;
gnusa << "unset colorbox" << endl;
gnusa << "unset key" << endl;
gnusa << "set palette model CMY rgbformulae 7,5,15" << endl;

gnusa << "set format xy \"%.1f\"" << endl;
gnusa << "set format y \"%.2f\"" << endl;
gnusa << "#set contour base" << endl;
gnusa << "#set cntrparam levels discrete " << conf67[1] << ", " << conf95[1] << endl;
gnusa << "set cntrparam levels discrete " << conf67[1]*N3[1] << ", " << conf95[1]*N3[1] << endl;
gnusa << "set cntrparam bspline" << endl;
gnusa << "set cntrparam points 7" << endl;
gnusa << "set cntrparam order 10" << endl;
gnusa << "#unset surface" << endl;
gnusa << endl << endl << endl << endl << endl;



k=0;
//for(i=0;i<4;i++) for(j=0;j<=i;j++){
for(j=0;j<4;j++) for(i=j;i<4;i++){
	
//x:j y:i



gnusa << "set tmargin at screen " << 0.990-0.200*i << endl;
gnusa << "set bmargin at screen " << 0.790-0.200*i << endl;
gnusa << "set lmargin at screen " << 0.150+0.200*j << endl;
gnusa << "set rmargin at screen " << 0.350+0.200*j << endl << endl;


//gnusa << "set ylabel \"Kinetic energy [1e51 erg]\"  font \", 20\" offset 46.5" << endl;

//gnusa << "set xtics 5, 3" << endl;
//gnusa << "set ytics 0.4, 5" << endl;
//gnusa << "set tics font \", 20\"" << endl;
//gnusa << "set cbtics font \", 1\"" << endl;
//gnusa << "unset y2label " << endl;


	if(i!=3){
	    gnusa <<  "unset xlabel" << endl;
	    gnusa << "set xtics font \"time, 0.01\"" << endl;
		}
	else{
		if(j==0) gnusa <<  "set xlabel \"Ejected mass \\n[M_{sol}]\"  font \", 11\"   offset 0, 0.5" << endl;
		if(j==1) gnusa <<  "set xlabel \"Kinetic energy \\n[10^{51} erg]\"  font \", 11\"   offset 0, 0.5" << endl;
		if(j==2) gnusa <<  "set xlabel \"Radius \\n[10^{12} cm]\"  font \", 11\"   offset 0, 0.5" << endl;
		if(j==3) gnusa <<  "set xlabel \"Thermal energy \\n[10^{51} erg]\"  font \", 11\"   offset 0, 0.5" << endl;
		gnusa << "set xtics font \", 10\"" << endl;
	}
	
	
	if(j!=0){
		gnusa <<  "unset ylabel" << endl;
		gnusa << "set ytics font \"time, 0.01\"" << endl;
		}
	else{
		if(i==0) gnusa <<  "set ylabel \"Posterior\"  font \", 11\"" << endl;
		//if(i==1) gnusa <<  "set ylabel \"Ejected mass [M_{sol}]\"  font \", 18\"" << endl;
		if(i==1) gnusa <<  "set ylabel \"Kinetic energy \\n[10^{51} erg]\"  font \", 11\"" << endl;
		if(i==2) gnusa <<  "set ylabel \"Radius \\n[10^{12} cm]\"  font \", 11\"" << endl;
		if(i==3) gnusa <<  "set ylabel \"Thermal energy \\n[10^{51} erg]\"  font \", 11\"" << endl;
		gnusa << "set ytics font \", 10\"" << endl;
	}
	
	







    //general
	if(i==j){
		gnusa << "unset ytics" << endl;
		gnusa << "set y2tics font \", 12\"" << endl;
		gnusa << "unset logscale y" << endl;
		}
	
	
	else{
		gnusa << "set ytics" << endl;
		gnusa << "unset y2tics" << endl;
		gnusa << "set logscale y" << endl;
		}




if(whole==1){

	if(i==j){
		gnusa << "set xtics "<<minmax2[i][0]<<"**0.8*"<<minmax2[i][1]<<"**0.2, ("<<minmax2[i][1]<<"/"<<minmax2[i][0]<<")**0.3" << endl;
		gnusa << "set y2tics 0.2*"<<(i<2?"b":"c")<<", 0.3*"<<(i<2?"b":"c") << endl;
		gnusa << "plot [" << minmax2[i][0] << ":" << minmax2[i][1] << "][0:"<<(i<2?"b":"c")<<"]  \"SpaceMeanW.txt\"  u " << 3*i+1 << ":" << 3*i+2 << "  w l " << endl  << endl;
		}
	
	
	else{
		gnusa << "set xtics "<<minmax[k][0]<<"**0.8*"<<minmax[k][1]<<"**0.2, ("<<minmax[k][1]<<"/"<<minmax[k][0]<<")**0.3" << endl;
		gnusa << "set ytics "<<minmax[k][2]<<"**0.8*"<<minmax[k][3]<<"**0.2, ("<<minmax[k][3]<<"/"<<minmax[k][2]<<")**0.3" << endl;
		gnusa << "splot [" << minmax[k][0] << ":" << minmax[k][1] << "][" << minmax[k][2] << ":" << minmax[k][3] << "][]  \"SpaceALLW.txt\"  u " << 4*k+1 << ":" << 4*k+2 << ":($" << 4*k+3 << ")  with pm3d  lw 2" << endl  << endl;
		}


}

else{
//<<(i<2?"b":"c")<<
	if(i==j){
		gnusa << "set xtics "<<minmax2[i][0]<<"**0.8*"<<minmax2[i][1]<<"**0.2, ("<<minmax2[i][1]<<"/"<<minmax2[i][0]<<")**0.3" << endl;
		gnusa << "set y2tics 0.2*"<<(i<2?"b":"c")<<", 0.3*"<<(i<2?"b":"c") << endl;
		gnusa << "plot [" << minmax2[i][0] << ":" << minmax2[i][1] << "][0:"<<(i<2?"b":"c")<<"]  \"SpaceMeanT.txt\"  u " << 3*i+1 << ":" << 3*i+2 << "  w l " << endl  << endl;
		}
	
	
	else{
		gnusa << "set xtics "<<minmax[k][0]<<"**0.8*"<<minmax[k][1]<<"**0.2, ("<<minmax[k][1]<<"/"<<minmax[k][0]<<")**0.3" << endl;
		gnusa << "set ytics "<<minmax[k][2]<<"**0.8*"<<minmax[k][3]<<"**0.2, ("<<minmax[k][3]<<"/"<<minmax[k][2]<<")**0.3" << endl;
		gnusa << "splot [" << minmax[k][0] << ":" << minmax[k][1] << "][" << minmax[k][2] << ":" << minmax[k][3] << "][]  \"SpaceALLT.txt\"  u " << 4*k+1 << ":" << 4*k+2 << ":($" << 4*k+3 << ")  with pm3d  lw 2" << endl  << endl;
		}

		
	
}



//gnusa << "set xlabel \"Radius [1e12 cm]\"  font \", 20\"" << endl;
//gnusa << "set ylabel \"Kinetic energy [1e51 erg]\\nThermal energy [1e51 erg]\"  font \", 20\" offset 3.5" << endl;
//if(i==0) gnusa << "set cblabel \"log of chi**2\"  font \", 20\" offset -1.5" << endl;



gnusa << endl << endl << endl << endl;
if(i!=j) k++;

} // i



gnusa << "unset multiplot" << endl << "set out" << endl << endl << endl << endl;


gnusa.close();



*/











printf("\n\nEverything writen to outputs. (mean.log)\n");
printf("done\n"); 
//getchar();
return 0;   
    
}


