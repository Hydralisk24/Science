#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream> // double -> string-hez kell
#include <iomanip> // set precision

#define day 86400.0		/*1 day in seconds*/
//#define xmin 0.2		/*Minimum ionization radius*/	
#define Eni 3.89e10		/*Ni energy generation rate (erg/g/s)*/
#define Eco 6.8e9		/*Co energy generation rate (erg/g/s)*/
#define Ekin 2.18e8
#define Ean 3.63e8
#define c 2.99e10		/*Speed of light (cm/s)*/
#define Msol 1.989e33		/*Solar mass (g)*/
#define kgamma 0.028		/*Opacity of gamma-rays (cm^2/g)*/
#define kpoz 7.0		/*Opacity of pozitrons (cm^2/g)*/

#include "Fejlec.cpp"


using namespace std;


/*
Based on LC2.2 by Nagy Andrea
Upgrade, Marcow Chain Monte Carlo solutions by Jäger Zoltán

MCMC:
parameter file:
parametersMC.inp
comments after the values are enabled, # also

tkezd     Measure time shift 
tmin      Time where the scatter count starts
Tion 	  Ionization/recombination temperature (K)
			Mni 	  Initial nickel mass (M_sun), no longer used!
a         Exponential density profile exponent
s		  Power-low density profile exponent
kappa     Thomson scattering opacity (cm^2/g)
Ep		  Initial magnetar rotational energy (erg)
tp		  Characteristic time scale of magnetar spin-down (d)
			Ag		  Gamma-leak (d^2), no longer used!
tmax	  Final epoch (day)
NMC	      Number of MC loops
tfark     Time where the tail starts

tkezd: initial shift: JD of shock breakout (in the program, the shock breakout is at 0)
tmin: the fitting starts at this time. Necessery beacuse the first days may cannot be fitted (i.e. because of two component)
tfark: the start time of the tail
Tion, s, kappa: if negativ, then this is fitted, else fixed
a: constant, if a!=0 and s=0 the a used instead of s. No significant difference from s
Ep, tp: currently consant, and cannot be fitted

Mni, Ag fitting: these are independent, and subscribed well by the late light curve (after the plateu).
The of of these are also done by this program (with mcmc).
First time; if there is no NiAgBest.txt file, the niag section starts
fits the tail with Ni and Ag
Then gives the best Ni values for every Ag, or the Ni(Ag) function: NiAgBest.txt
if there is NiAgBest.txt the normal mod starts
M and Ekin gives Ag, and Ag gives Ni (from niag)

Light Curve file:
bol-gorbe.txt
no comment possible currently
format: time[day]  bolometric_absolute_mag  mag_err
error is necessery (if not exist, then be 1)

the program fits:
R0: progenitor radius
M: ejected mass
Ek: kinetic energy
Eth: thermal energy
(Tion, s, kappa is set to)


Marcow Chain Monte Carlo (mcmc), Metropolis–Hastings method
The fitted parameters are scattered randomly in the accaptable region
The program calculates the scatter/khi from the calculated light curve and the bolometric light curve
The algorithm gives more sample at the spots where the scatter is lower, and fewer at where it is higher
Then writes out the parameters and the scatter
The program calculates 2 scatter: one between tmin and tmax (all)
And one between tmin and tfark, and this one used at the mcmc algorithm!!!!!
Because this region described by the fitted parameters, while the tail is fitted independently 

Methods:
Lsn=L+Lion+Lpoz
L is full numeritic at first: 4th grade Runge-Kutta method
Lion is semi analitic: dr=dxi*R component is numerical (calculated at L)
Lpoz is full analitic
When xi=xmin, it becomes constant, dr=dxi=0 -> Lion=0, and L becomes full analitic 
(however if xmin>0.2, L has a transient phase after xi=xmin, and that needs to be calculated numerically)
In full numeric phase dt=40000s
In the numeric phase dt=20000s * fxi(xi)
fxi(xi)=1 if xi>0.4, at lower xi, better accuracy necessery, hence fxi(xi)<1
xmin should be 0.1, but 0.2 used because it don't have significant effect to the light curve, but gives numerical errors
If the scatter between tmin and tmin+50 is bigger then 0.5 mag, then xmin=0.3, this speeds up the program

mcmc:
x1 = x0 + rand, which is gauss(sigma , mu=0)
R0=tR0*pow(10.,RandomDoublegauss( Beta*sqrt(2.5*tszoras) , 0) );  // log10()
log(R0)=log(tR0)+RandomDoublegauss( Beta*sqrt(2.5*tszoras) , 0);  // log10 smooth
sigma = Beta*sqrt(2.5*tszoras), function of the scatter


standard output: prints how many done
+ if the scatter is lower then a threshold (used for monitoring, if there are changing +, you doing it well)

output:
R0 [cm]    M [M_sol]    Ek [erg]    Eth [erg]         v [km/s]       T0 [K]    Mni [M_sol]     kappa [g/cm^2]    s    Tion [K]      scatter:all [mag]      scatter:till tail [mag]

v: velocity

The output then can be plotted, to see where are the solutions

not independent:
Ek and M  (and kappa)
Eth and R0
*/



double xmin=0.2;		/*Minimum ionization radius*/	
//double txmin=0.2;		/*saved Minimum ionization radius*/

double dt=1.0;	    	/*Time step (s)*/
double segedvalt1=0;    /*Assistant vairables*/
double segedvalt2=0;
double segedvalt3=0;
double E,Mni,Ek,F,L,T0,T,Tc,Lion,Ith,IM,INi,Eth,Lsn;
double t,td,th,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,Z,A,Ag;
double M,v,R0,R,rho,ri,dr,Q,tni,tco,Em,Em1;
double Lpoz,Ag1;
/*
double E,Mni,Ek,F,L,T0,T,Tc,Tion,Lion,Ith,IM,INi,Eth,Lsn;
double t,td,th,tmax,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,a,Z,A,Ag,Ep,s;
double kappa,M,v,R0,R,rho,ri,dr,Q,tni,tco,tp,Em,Em1;
double Lpoz,Ag1;
*/
//double tmin, tkezd, tfark;
double ER0;

// MC komponent
//double NMC;   // Number of MC loops

double logLh=0;      // prev logL
double tprev=0;		 // prev t
double Tn;           // measurement number
double Tnf;			 // measurement number till tail
double Tnx;			 // measurement number till tmin + 50
double Adat[500][3]; // meauserement data

double szoras;    // scatter
double szorash;   // scatter prev
double szorast;   // scatter temp
double tszoras;   // scatter save
double szoraskuszob=0.2; // scatter threshold
double szorasf;  // scatter tail

// scatter = sqrt( chisquare/norma ), norma = Tn = 1/error^2
double khinegyzet;    // chisquare
double khinegyzetf;    // chisquare tail
double khinegyzeth;    // chisquare prev


// step paramters
double dR0;   	  /*Initial radius (cm)*/
double dM;   	  /*Ejecta mass (M_sun)*/
double dEk;   	  /*Initial kinetic energy (1e51 erg)*/ 
double dE; 		  /*Intital thermal energy (1e51 erg)*/ 
double dkappa; 	  /*Thomson scattering opacity (cm^2/g)*/
double ds; 	      /*Power-low density profile exponent*/	  
double dTion;	  /*Ionization/recombination temperature (K)*/
double dMni;   	  /*Initial nickel mass (M_sun) */
double dAg;   	  /*Gamma-leak (d^2)*/

double dEp=0;
double dtp=0;

double tR0=0;   	  // temp parameter (accepted mcmc chain holder)
double tM=0;   	 
double tEk=0;   	 
double tE=0; 	
double tkappa=0; 	
double ts=0; 	  
double tTion=0;
double tMni=0;
double tAg=0;
double tER0=0;
double tEp=0;
double ttp=0;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// only for store
double tszorasf; 
double tkhinegyzetf; 
double tkhinegyzet; 
double tv; 
double tT0;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*
double R0min=0;   	  // Acceptable minimum
double Mmin=0;   	 
double Ekmin=0;   	 
double Emin=0; 	
double kappamin=0; 	
double smin=0; 	  
double Tionmin=0;
double Mnimin=0;
double Agmin=0;
double ER0min=0;
double Epmin=0;
double tpmin=0;

double R0max=0;   	  // Acceptable maximum
double Mmax=0;   	 
double Ekmax=0;   	 
double Emax=0; 	
double kappamax=0; 	
double smax=0; 	  
double Tionmax=0;
double Mnimax=0;
double Agmax=0;
double ER0max=0;
double Epmax=0;
double tpmax=0;
*/

//double para=0;	  // parameter a for manual MNi(Ag)
//double parb=0;	  // parameter b


// Velocity priori (Gaussian) !!!!!!!!!!!!!!!!!!! (if 0, then unset) Hardcoded
//double vref=0;  // mean, in [cm/s] e.g: 5e8
//double vrefstd=0; // error, in [cm/s] e.g: 1e8




					//double Beta=0.1;	  // Metropolis–Hastings step parameter, default 0.1, set in parameter file
					//double Delta=30000;



int niag=0; // niag mode?
//int nifit=1; // fit ni; niag mode?



// tar[] always bigger than nn! May cause problem if not!!!!! (default 1000 and 1050)   !!!!!!!!!!!!!!!!!
// Do not change! This gives the resolution of the numerical Ni(Ag) function, and used there as a holder
int nn=1000;
double tar[1050][3];






// Random generator 64 bit
unsigned 
rand256()
{
    static unsigned const limit = RAND_MAX - RAND_MAX % 256;
    unsigned result = rand();
    while ( result >= limit ) {
        result = rand();
    }
    return result % 256;
}

unsigned long long
rand64bits()
{
    unsigned long long results = 0ULL;
    for ( int count = 8; count > 0; -- count ) {
        results = 256U * results + rand256();
    }
    return results;
}


double RandomDouble(double min, double max){
    if(min == max) return min;
    //What if the user reversed min and max?
    if(min > max) {double tmp = min; min = max; max = tmp;}
    
	double v;

    // 5.4210108624275221700372640043497e-20 = pow(2,-64)
	v=( 1.0*rand64bits()* 5.42101086242752217e-20 )*(max-min)+min;
	
	return v;
	}
	
	
double RandomDoublelog(double min, double max){
    if(min == max) return min;
    //What if the user reversed min and max?
    if(min > max) {double tmp = min; min = max; max = tmp;}
    
	double v;

	v=( exp( RandomDouble(0.,log(20)) )-1 )/19 *(max-min)+min;
	
	return v;
	}
	
	
		
double RandomDoubleGauss(double szigma, double mu){
	double v;

	v= sqrt( -2.0*log( RandomDouble(0.0, 1.) ) ) * cos(2.*M_PI * RandomDouble(0., 1.) ); 
	if(v>2)  v=RandomDouble(2., 10.);
	if(v<-2) v=-RandomDouble(2., 10.);
	
	if(v*0!=0) return -1000;
	return v*szigma+mu;
	}
	
	
	
	





/*


// stringbol kimasolja az a es b kozott szoveget
// copy the string between a and b
string copy(string str, int a, int b){
	   
	   if(a>b) return "";
	   
	   string y;
	   int i;
	   
	   for(i=a;i<=b;i++)
	   y=y+str[i];
	   
	   return y;
	   }



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




double read2(string str, int n){

	double y=-1000;
	//string str;
	int i;
	int j;
	int k;
	//c++ fajl beolvasas, a sort olvassa be 
	//do{  std::getline(f, str);  } while(str[0]=='#');
	

	j=0;
    for(k=0;k<n;k++){	
        if(k!=0) j=i+2;
		while(str[j]==' ') j++;
		i=j;
		while(str[i]!=' ' && i<=str.length()) i++;
	}
	
	y = atof(   copy(str,j,i-1)    .c_str());
	
    
    return y;
}






 void defaultdata()			 //Writeing default input file
 {
 
 FILE *f;

 f=fopen("parametersMC.inp","wt");	
 
   fprintf(f,"%d\t[Data time shift (day)]\n",0);   	  //Measure time shift
   fprintf(f,"%d\t[Time where the fitting starts (day)]\n",0);    	  //Time where the scatter count starts
   fprintf(f,"%d\t[Ionization/recombination temperature (K)(-: fit)]\n",5500); 	  //Ionization/recombination temperature (K)
   //fprintf(f,"%lf\t[Initial nickel mass]\n",Mni);  	  //Initial nickel mass (M_sun)   
   fprintf(f,"%d\t[Exponential density profile exponent]\n",0);            //Exponential density profile exponent
   fprintf(f,"%d\t[Power-low density profile exponent (-: fit)]\n",0);		  //Power-low density profile exponent
   fprintf(f,"%1.1f\t[Thomson scattering opacity (cm^2/g)(-: fit)]\n",0.3);     	  //Thomson scattering opacity (cm^2/g)
   fprintf(f,"%d\t[Initial magnetar rotational energy (erg)]\n",0);		  //Initial magnetar rotational energy (erg)
   fprintf(f,"%d\t[Characteristic time scale of magnetar spin-down (day)]\n",0);		  //Characteristic time scale of magnetar spin-down (d)
   //fprintf(f,"%lf\t[Gamma-leak]\n",Ag);		  //Gamma-leak (d^2)
   fprintf(f,"%d\t[Final epoch (day)]\n",300);	  //Final epoch (day)
   fprintf(f,"%d\t[Number of MC loops]\n",300000);	  //Number of MC loops
   fprintf(f,"%d\t[Time where the plateu ends (day)]\n\n",100);    	  //Time where the plateu ends
   
   fprintf(f,"%1.0e\t[MIN Initial radius (cm)]\n",2e12);   	  //R0min
   fprintf(f,"%1.0e\t[MAX Initial radius (cm)]\n",80e12);   	  //R0max
   fprintf(f,"%d\t[MIN Ejected mass (Msol)]\n",3);   	  //Mmin
   fprintf(f,"%d\t[MAX Ejected mass (Msol)]\n",60);   	  //Mmax
   fprintf(f,"%1.1f\t[MIN Kinetic energy (foe)]\n",0.1);   	  //Ekmin
   fprintf(f,"%d\t[MAX Kinetic energy (foe)]\n",30);   	      //Ekmax
   fprintf(f,"%1.1f\t[MIN Thermal energy (foe)]\n",0.1);   	  //Emin
   fprintf(f,"%d\t[MAX Thermal energy (foe)]\n",30);   	      //Emax
   fprintf(f,"%1.2f\t[MIN Thomson scattering opacity (cm^2/g)]\n",0.05);   	  //kappamin
   fprintf(f,"%1.1f\t[MAX Thomson scattering opacity (cm^2/g)]\n",0.4);   	  //kappamax
   fprintf(f,"%d\t[MIN Power-low density profile exponent]\n",0);   	  //smin
   fprintf(f,"%d\t[MAX Power-low density profile exponent]\n",3);   	  //smax
   fprintf(f,"%d\t[MIN Ionization/recombination temperature (K)]\n",0);   	      //Tionmin
   fprintf(f,"%d\t[MAX Ionization/recombination temperature (K)]\n",20000);   	  //Tionmax
   fprintf(f,"%1.0e\t[MIN Initial nickel mass (Msol)]\n",0.0005);   	  //Mnimin
   fprintf(f,"%1.0e\t[MAX Initial nickel mass (Msol)]\n",0.5);   	      //Mnimax
   fprintf(f,"%1.0e\t[MIN Gamma-leak (day^2)]\n",1e3);   	  //Agmin
   fprintf(f,"%1.0e\t[MAX Gamma-leak (day^2)]\n",1e7);   	  //Agmax
   fprintf(f,"%d\t[Expansion velocity mean priori (km/s)(0: not fit)]\n",0);   	  //vref
   fprintf(f,"%d\t[Expansion velocity sigma priori (km/s)(gauss)]\n",0);   	  //vrefstd
   fprintf(f,"%d\t[Determinate the nickel mass Ag function? (binary)]\n",1);	  //nifit
   fprintf(f,"%1.1f\t[Metropolis–Hastings step parameter, default 0.1]\n",0.1);	  //Beta
   fprintf(f,"%d\t[Step parameter (exponental) decay (loop number), default 30000]\n",30000);	  //Delta
   fprintf(f,"%1.1f\t[Minimum ionisation zone (between 0.1 and 1)]\n",0.2);	  //txmin
   fprintf(f,"%1.1f\t[Threshold for szoras, technical parameter, default 0.2]\n\n",0.2);	  //szoraskuszob
   
   fprintf(f,"%d\t[Whole (W) =1 or No Tail (T) =0 used?]\n",0);	  //whole
   fprintf(f,"%d\t[Plotting grid number]\n",80);	  //grid
   fprintf(f,"%d\t[burn in (steps)(=Step decay suggested)]\n",30000);	  //burnin
   


 fclose(f);
 	
 }
 
 
 
void data()			 //Reading input file
 {
 //FILE *f;
 //f=fopen("parametersMC.inp","rt");

 if( !f.is_open() ){printf("No parameter file! Generating a default! Exiting!\n"); defaultdata(); exit(1);}

   
   //fscanf(f,"%lf",&tkezd);   	  //Measure time shift
   //fscanf(f,"%lf",&tmin);    	  //Time where the scatter count starts
   //fscanf(f,"%lf",&Tion); 	  //Ionization/recombination temperature (K)
   //fscanf(f,"%lf",&Mni);  	  //Initial nickel mass (M_sun) 
   //fscanf(f,"%lf",&a);            //Exponential density profile exponent
   //fscanf(f,"%lf",&s);		  //Power-low density profile exponent
   //fscanf(f,"%lf",&kappa);     	  //Thomson scattering opacity (cm^2/g)
   //fscanf(f,"%lf",&Ep);		  //Initial magnetar rotational energy (erg)
   //fscanf(f,"%lf",&tp);		  //Characteristic time scale of magnetar spin-down (d)
   //fscanf(f,"%lf",&Ag);		  //Gamma-leak (d^2)
   //fscanf(f,"%lf",&tmax);	  //Final epoch (day)
   //fscanf(f,"%lf",&NMC);	  //Number of MC loops
   //fscanf(f,"%lf",&tfark);    	  //Time where the plateu ends
   
   tkezd=read();
   tmin=read();
   Tion=read();
   //Mni=read();
   a=read();
   s=read();
   kappa=read();
   Ep=read();
   tp=read();
   //Ag=read();
   tmax=read();
   NMC=read();
   tfark=read();
   read();
   
   R0min=read();
   R0max=read();
   Mmin=read();  
   Mmax=read(); 	 
   Ekmin=read();  
   Ekmax=read(); 	 
   Emin=read(); 
   Emax=read(); 	
   kappamin=read(); 
   kappamax=read();	
   smin=read(); 
   smax=read(); 	  
   Tionmin=read();
   Tionmax=read();
   Mnimin=read();
   Mnimax=read();
   Agmin=read();
   Agmax=read();
   vref=read();
   vrefstd=read();
   nifit=int(read());
   Beta=read();
   Delta=read();
   txmin=read();
   szoraskuszob=read();


 //++++++++++++++++++++++++++++++++++++++++++++
 
 if(s==3.0 || s==5.0){printf("Parameter error! The n = 3.0 and n = 5.0 are forbidden!\n"); exit(1);}
 
 // Default acceptable region set if the input is wrong (between physically possible range) 
 if(R0min<=0  || R0max<=R0min)    { R0min=2e12;   R0max=80e12; }
 if(Mmin<=0   || Mmax<=Mmin)      { Mmin=3;   Mmax=60; }
 if(Ekmin<=0  || Ekmax<=Ekmin)    { Ekmin=0.1;   Ekmax=30; }
 if(Emin<=0   || Emax<=Emin)      { Emin=0.1;   Emax=30; }
 if(kappamin<=0 || kappamax<=kappamin){ kappamin=0.05;   kappamax=0.4; }
 if(smin<0    || smax<=smin)      { smin=0;   smax=3; }
 if(Tionmin<0 || Tionmax<=Tionmin){ Tionmin=0;   Tionmax=20000; }
 if(Mnimin<=0 || Mnimax<=Mnimin)  { Mnimin=0.0005;   Mnimax=0.5; }
 if(Agmin<=0  || Agmax<=Agmin)    { Agmin=1e3;   Agmax=1e7; }
 if(Beta<=0){ Beta=0.1; }
 if(Delta<=0){ Delta=0.01; }
 if(txmin<0.1){ txmin=0.2; }
 if(szoraskuszob<=0){ szoraskuszob=0.2; }
 
 if(ER0min<=0  || ER0max<=ER0min)    { ER0min=0.1e63;   ER0max=300e63; }
 if(Epmin<=0   || Epmax<=Epmin)      { Epmin=0.01e51;   Epmax=30e51; }
 if(tpmin<=0   || tpmax<=tpmin)      { tpmin=864;   tpmax=100*86400; }
 
 //para=Mni;
 //parb=Ag;
 
 
 // initialization
 R0=40e12;
 M=30;
 Ek=5;
 E=5;
 
 xmin=txmin;
 
 // Positive: fixed. Negative: MC
 if(kappa<0){ dkappa=Beta*0.2; kappa=0.5*(kappamin+kappamax); }
 if(s<0){ ds=Beta*1.5; s=0.5*(smin+smax); }
 if(Tion<0){ dTion=Beta*10000; Tion=0.5*(Tionmin+Tionmax); }
 
 if(Ep<0) dEp=1; else dEp=0;
 if(tp<0) dtp=1; else dtp=0;
 
 //szoraskuszob=0.2;

 M=Msol*M;	 		
 Mni=Msol*Mni;
 E=1e51*E;
 Ek=1e51*Ek;
 Ep=1e51*Ep;
 tp=tp*86400.0;
 Ag=Ag*pow(86400.0,2);		
 T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 
 
 Mmin=Msol*Mmin;
 Mmax=Msol*Mmax;
 Emin=1e51*Emin;
 Emax=1e51*Emax;
 Ekmin=1e51*Ekmin;
 Ekmax=1e51*Ekmax;
 Mnimin=Msol*Mnimin;
 Mnimax=Msol*Mnimax;
 Agmin=Agmin*pow(86400.0,2);
 Agmax=Agmax*pow(86400.0,2);
 vref=vref*1e5;
 vrefstd=vrefstd*1e5;
 
 //Epmin=1e51*Epmin;
 //Epmax=1e51*Epmax;
 //tpmin=86400*tpmin;
 //tpmax=86400*tpmax;
 
 

 dR0=Beta*R0;   
 dM=Beta*M;   
 dEk=Beta*Ek;   	
 dE=Beta*E; 
 

 //fclose(f);
 f.close();
 }

*/

void meres()			 /*Reading input measurement file*/
 {
 int i;
 int j;
 int k;
 int l;
 FILE *f;

 f=fopen("bol-gorbe.txt","rt");
 if(f==NULL){printf("No measurement file!\n"); exit(1);}
 
 double temp;
 double csere[3];

 i=0; Tn=0; Tnf=0; Tnx=0;
 if(f!=NULL) while (!feof(f)){    
 								  // read in formats
 								  fscanf(f,"%lf %lf %lf %lf %lf\n",&Adat[i][0],&temp,&temp,&Adat[i][1],&Adat[i][2]);
 								  //fscanf(f,"%lf %lf %lf\n",&Adat[i][0],&Adat[i][1],&Adat[i][2]);
 
 								  //converting to log erg
 			 	   				  Adat[i][0]=Adat[i][0]-tkezd;
 			 	   				  Adat[i][1]=-0.4*Adat[i][1]+(71.21+17.5)*0.4;
 			 	   				  Adat[i][2]=0.4*Adat[i][2];
								  //printf("%f %f %f\n",Adat[i][0],Adat[i][1],Adat[i][2]); 
								  
								  // calculating the measurement numbers for scatter (all, befor tail, and between tmin and tmin+50 (technical))
								  if(Adat[i][0]>=tmin && Adat[i][0]<=tfark) Tn=Tn+1/Adat[i][2]/Adat[i][2];
								  if(Adat[i][0]>=tmin && Adat[i][0]<=tmax) Tnf=Tnf+1/Adat[i][2]/Adat[i][2];
								  if(Adat[i][0]>=tmin && Adat[i][0]<=tmin+50) Tnx=Tnx+1/Adat[i][2]/Adat[i][2];
								  //if(Adat[i][0]>=tmin && Adat[i][0]<=tmax) Tn=Tn+1;
								  //tmax here is in day, later will be in s
                              i++;}
                              

 // arrenge by time
 for(j=0;j<i;j++){
					 
					 for(k=0;k<i-1;k++)					  
                    	if(Adat[k][0]>Adat[k+1][0]) 
						    for(l=0;l<3;l++){ csere[l]=Adat[k][l]; Adat[k][l]=Adat[k+1][l]; Adat[k+1][l]=csere[l]; }
					
					}
					
 if(Adat[0][0]<0){printf("First date is negative!\n"); exit(1);}
 if(Adat[0][0]>tmax){printf("First date is out of limit!\n"); exit(1);}

 fclose(f);
 }
 
 
 
 void niagread()			 /*Reading input Niag fitting file if exist, else fits one*/
 {
 int i;
 int j;
 int k;
 int l;
 double temp;
 FILE *f;

 f=fopen("NiAgBest.txt","r");
 if(f==NULL){ niag=1; return;}
 else niag=0;

 
 i=0;
 while (!feof(f)){    
                   
				   fscanf(f,"%lf %lf %lf",&tar[i][0],&tar[i][1],&tar[i][2]);
				   i++;
                              				  				  
				   }


 nn=i;
 fclose(f);
 }







// Physx  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double psi(double x)			/*Temperature profile assuming uniform density*/
	   {
	   double z;
	   if(x==1.0) z=0.0; else if(x>0) z=sin(M_PI*x)/(M_PI*x); else z=1.0;
	   return z;
	   }



double eta(double x, double a1)		/*Exponential density profile*/
	   {
	   double z;
	   if(x<xmin) z=1; else z=exp(-a1*(x-xmin));
	   return z;
	   }



double theta(double x, double a2)	/*Power-low density profile*/
	   {
	   double z;
	   if(x<xmin) z=1; else z=pow(x/xmin,-a2);
	   return z;
	   }


	
	
double temp(double T1,double y0)	/*Find ionization radius*/
	   {
	   if(Tion==0) return y0;
	   double y,dy,zone;
	   double a,b;
	   a=xmin; b=y0;
	   dy=4e-6;

	   					
	   while(b-a>=dy)
	   	{
		 	   y=(a+b)/2;
			   T=T1*pow(psi(y),0.25);
			   if(T==Tion) break;
			   if(T>Tion){ a=y; }
			   else{ b=y; }
		} 
		
		y=(a+b)/2;
		

	if(y<y0) zone=y+0.5*dy; else zone=y0;
	return zone;	
	} 
	
	
	
/*
double temp(double T1,double y0)
{
double y,dy,zone;
y=y0;
dy=1e-9;
T=T1*pow(psi(1.0),0.25);

  while(y>=xmin && T<Tion)
  {
     T=T1*pow(psi(y),0.25);
     y=y-dy;
  }

if(y<y0) zone=y+0.5*dy; else zone=y0;
return zone;
}
*/


double IM_int(double b, double a1, double a2)
	   {
	   double sum=0.0,x=0.0,dx=1e-3;
	   while(x<b)
	   {
	   //sum=sum+(eta(x,a1)*pow(x,a2)+eta(x+dx,a1)*pow(x+dx,a2))*dx*0.5; 	/*Trapesoid rule*/
	   sum=sum+	   (eta(x,a1)*pow(x,a2)+4.*eta(x+0.5*dx,a1)*pow(x+0.5*dx,a2)+eta(x+dx,a1)*pow(x+dx,a2))*dx*1./6.; 	/*Simpson's rule*/
	   x=x+dx;
	   }
	   return sum;
	   }



double Ith_int(double b)
	   {
	   double sum=0.0,x=0.0,dx=1e-3;

	   while(x<b)
	   {
	   //sum=sum+(psi(x)*x*x+psi(x+dx)*pow(x+dx,2.0))*dx*0.5;
	   sum=sum+	   (psi(x)*x*x+4*psi(x+0.5*dx)*pow(x+0.5*dx,2.0)+psi(x+dx)*pow(x+dx,2.0))*dx*1./6.;
	   x=x+dx;
	   }
	   return sum;
	   }
	   
	   



double series(double x)			/*Taylor-series of sin(pi*x)/(pi*x)*/
	   {
	   double z;
	   z=1-pow(M_PI*x,2)/6.0+pow(M_PI*x,4)/120.0-pow(M_PI*x,6)/5040.0+pow(M_PI*x,8)/362880.0;
	   //z=sin(M_PI*x)/(M_PI*x);
	   return z;
	   }





	   
// F Runge-Kutta
double Frk(double t, double F, double i){
	  if(F<0) F=0;
	  double dF;
	  
	  double sigma;
	  double xil=xi;
	  double g;
	  double dxi;
	  double Em;
	  double R;
	  double Tc; 

  	  if(tp==0) Em=0;
  	  else Em=Ep/(tp*pow((1.0+t/tp),2.0));
  	  
	  R=R0+v*t;	   
	  sigma=R/R0;			
	  Tc=T0*pow(F,0.25)/sigma;
	 
	  if(xil>=xmin && Tc>Tion){ xil=temp(Tc,xh); }
	  else{ xil=xmin; }
	  segedvalt1=xil;
	 
	 
	  dxi=(series(xh)-series(xil))*xil/i;
	  segedvalt2=dxi;
	 
	  	
	  //Xni es Xco analitical
	  Xni=  exp( -t/tni );   
	  Xco= (-exp( -t/tni ) + exp( -t/tco )) / (-tni*(1./tco-1./tni));     	  	  
  	  g=Xni+p5*Xco;	  	  
	  
	  /*Runge-Kutta method*/	     
	  dF=sigma/pow(xil,3.0)*(p1*g-p2*F*xil-2.0*F*pow(xil,2.0)*dxi*tni/(sigma*dt)+p3*Em);
	  
	  
	  return dF;
	  }
	  
	  
	  
// dinamical dt from xi	 (empirical)
double fxi(double x){
	   double y=1;
	   
	   if(x==xmin) return y;
	   if(x<0.4) y=0.01+(x-0.1)*(x-0.1)*11.;
	   
	   return y;
	   }
	  






// Misc +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   
	   
// square	  
double h(double x){
	  return x*x;
	  }
	  
	   
	   
// linear interpolation
double interpol(double x, double xm, double ym, double xp, double yp){
	  if(xp==xm) return 0;
	  double y;
	  
	  y=ym+ (x-xm) * (yp-ym) / (xp-xm);
	  
	  return y;
	  }
	  
	  
	  
// gives back the array number, for Ni Ag, numerical
int inverztar(double x, double a, double b, double nn){
	double y;
	
	//tar[i][0]=a + i *(b-a)/double(nn);
	y = (x - a) * double(nn)/(b-a);
	
	return int(y);
    }

	  
	    
// funcion MNi(Ag), for best fit, empirical	  
/*
double fMni(double x){
	double y=0;
	
	if(x>parb) return para;
	if(x<parb/4.) y=pow( 10 , (  log10(para) + 0.75*log10(parb/2.2) - 0.75*log10(x))  );
	else          y=pow( 10 , (  log10(para) + 0.35*log10(parb)   - 0.35*log10(x))  );
	
	return y;
}
*/


// fitted values
double fMni(double x){
	double y=0;
	int i;
	
    i=inverztar( log10(x) ,log10(tar[0][0]),log10(tar[nn-2][0]),nn);
    if(i<0) i=0;  if(i>nn-1) i=nn-1;
    //y=tar[i][1];
    y=interpol(x, tar[i][0], tar[i][1], tar[i+1][0], tar[i+1][1]);
	
	return y;
}


	   
	   

	   
	   
	   
	   



int main()
{
	
srand(time(NULL));
	
int j=0;
int i=0; 
int lep=0; // for dinamical dt
int itt=0; // measure index
int ciklus;    // MC ciklus
int elso=0;      // first run?
int tilt=0;		 // After xi is constant, this paramter is 1, and the program becomes simplier
int accept=0;      // Is the mcmc step accepted?
//int skip=0;
double dF1,dF2,dF3,dF4;  // for RK
double opac,opac1,f,g1;

double acc=0;    // acceptance ratio helper
double accall=0;

tni=8.8*day;				/*Ni decay time (s)*/
tco=111.3*day;				/*Co decay time (s)*/
Z=1.0;					/*Average atomic charge in ejecta*/
A=1.0;					/*Average atomic mass in ejecta*/


data(); //parameter file read

 if(kappa<0){ dkappa=Beta*0.2; kappa=0.5*(kappamin+kappamax); }
 if(s<0){ ds=Beta*1.5; s=0.5*(smin+smax); }
 if(Tion<0){ dTion=Beta*10000; Tion=0.5*(Tionmin+Tionmax); }
 
 if(Ep<0) dEp=1; else dEp=0;
 if(tp<0) dtp=1; else dtp=0;
 
 T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 
 xmin=txmin;
 
 szoraskuszob=0.2;


// Chech if there are fitted Ni(Ag) function, if not, then the program fits IF not MCMC sampled
if(nifit!=0) niagread(); else{ niag=0; printf("No Ni(Ag) function seek! Ni MCMC sampled instead!\n"); }
/*
if nifit=0 then Ni mass is sampled in the MCMC method like everything else, and no further procedure happens
if nifit=1 then Ni mass fitted seperatly, determining the best Ni Ag pairs (function)
first run niag=0 and the program make this fitting in an MCMC fit but only in the tail
This make seperate block in the code
After this niag=1 and the normal MCMC fitting takes place
*/


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
string sor;  // seged
//FILE *fki;   // rand.out and rand2.out
//FILE *fki2;
//FILE *fkiS;  // seged
std::ifstream IN; // seged, used only here
std::ofstream fki;
std::ofstream fki2;
//std::ofstream fkiS;

if(niag!=1){
//fkiS=fopen("rand.out","r");
IN.open("rand.out");
//fki.open("rand.out");

/* Is this the first run? check */
if(/*fkiS==NULL*/ !IN.is_open()) elso=1; else{ 

/* Reading the last line of the previous runs. */
IN.seekg(-2,IN.end); // Az eof es soremeles ele lepunk vissza

     int i=IN.tellg(); //Innen lepunk vissza mindig egy pozicioval
     char ch = ' ';
     while(ch!='\n')
         {
            IN.seekg(--i,IN.beg); //A file elejetol szamoljuk a poziciot
            ch=IN.peek();         //Megnezzuk mi van ott, de nem mozditjuk a pointer
         }
         IN.seekg(++i,IN.beg);    //Mivel pont a sor emelesnel alt meg moge lepunk.
      //string sor;
      getline(IN,sor);
      cout << "Continuing previous work. Is this line correct? (Enter)\n\n";
      cout << sor << endl;
      getchar();

 		elso=0;
 }
 
 
//fki=fopen("rand.out","a");
fki.open("rand.out",std::ofstream::app);
//fki2=fopen("rand2.out","a"); 
fki2.open("rand2.out",std::ofstream::app);
}

//buffer for output
stringstream outtemp;
stringstream outtemp2;



if(niag==1){
//fki=fopen("rand-niag.out","a");
fki.open("rand-niag.out",std::ofstream::app);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(niag==1){ 
			//tkezd=0; tmin=150.0; tmax=400.0; 
			printf("\nSampling the nebular phase to determinate the Ni mass Ag function:\n300000 steps, please do not interrupt.\n");
			//printf("Initial shift:  "); scanf("%lf",&tkezd);
			//printf("Start time:  ");    scanf("%lf",&tmin);
			//printf("Max time:  ");      scanf("%lf",&tmax);
			//printf("\n\n");
			//printf("%f %f %f",tkezd,tmin,tmax);
			//getchar(); getchar();
			
			tmin=tfark+5;
			// shift it with 5 day!!!!!
			
}
if(niag==1){ tilt=1; NMC=300000; tfark=tmax; }

meres();  // measurement read


//getchar();


tmax=tmax*day; 			/*Final epoch (s)*/
tfark=tfark*day; 			/*Tail epoch (s)*/

szorash=1000;
khinegyzeth=1e100;
tszoras=1;

acc=0;
accall=0;


   //Mni  	  /*Initial nickel mass (M_sun) */
   //Ep		  /*Initial magnetar rotational energy (erg)*/
   //tp		  /*Characteristic time scale of magnetar spin-down (d)*/
   //Ag		  /*Gamma-leak (d^2)*/
   //R0   	  /*Initial radius (cm)*/
   //M    	  /*Ejecta mass (M_sun)*/
   //Ek   	  /*Initial kinetic energy (1e51 erg)*/ 
   //E 		  /*Intital thermal energy (1e51 erg)*/  
   
   
                 // First random value, between acceptable region (below)
    			 ER0=0.5*ER0max*RandomDoublelog(2.*ER0min/ER0max,2.); 
				 //R0=0.5*R0max*RandomDoublelog(2.*R0min/R0max,2.); 
   	  			 M=0.5*Mmax*RandomDoublelog(2.*Mmin/Mmax,2.);
				 Ek=0.5*Ekmax*RandomDoublelog(2.*Ekmin/Ekmax,2.);	  
				 //E=0.5*Emax*RandomDoublelog(2.*Emin/Emax,2.);
				 if(nifit==0) Mni=tMni=sqrt(Mnimin*Mnimax);
				 
				 if(dEp!=0) Ep=0.5*Epmax*RandomDoublelog(2.*Epmin/Epmax,2.);
				 if(dtp!=0) tp=0.5*tpmax*RandomDoublelog(2.*tpmin/tpmax,2.);
				 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
if(elso==1) fprintf(fki,"# Marcow Chain Monte Carlo solutions for Supernova light curve model based on Arnett & Fu ApJ 340, 396 (1989)\n");   

// multiply (and square!) the output chi with this to get the unnormalised chi square

if(elso==1) fprintf(fki,"# chi**2 norma (no tail) = %f\n",Tn/2.5/2.5);
if(elso==1) fprintf(fki,"# chi**2 norma (whole) = %f\n",Tnf/2.5/2.5);
  
if(elso==1) fprintf(fki,"# Start JD = %f\n",tkezd);
if(elso==1) fprintf(fki,"# First counted day = %f\n",tmin); 
if(elso==1) fprintf(fki,"# Exponential density profile exponent = %f\n",a); 
if(elso==1) fprintf(fki,"# Initial magnetar rotational energy (erg) = %e\n",Ep);
if(elso==1) fprintf(fki,"# Characteristic time scale of magnetar spin-down (d) = %e\n",tp);
if(elso==1) fprintf(fki,"# \n");
if(elso==1) fprintf(fki,"# ID 	 R0 [cm]     M [M_sol]    Ek [erg]      Eth [erg]          v [km/s]       T0 [K]        M_Ni [M_sol]    Ag [d^2]       kappa [g/cm^2]    s    Tion [K]       STD [mag], chi^2; whole and without tail\n");


if(elso==1) fprintf(fki2,"# Marcow Chain Monte Carlo solutions for Supernova light curve model which are not accapted\n");
if(elso==1) fprintf(fki2,"# Start JD = %f\n",tkezd);
if(elso==1) fprintf(fki2,"# ID 	 R0 [cm]     M [M_sol]    Ek [erg]      Eth [erg]          v [km/s]       T0 [K]        M_Ni [M_sol]    Ag [d^2]       kappa [g/cm^2]    s    Tion [K]       STD [mag], chi^2; whole and without tail\n");
*/

if(niag!=1){
if(elso==1)	cout << "Normal MCMC fitting\n\n";
if(elso==1) fki << "# Marcow Chain Monte Carlo solutions for Supernova light curve model based on Arnett & Fu ApJ 340, 396 (1989)\n";   

if(elso==1) fki << "# chi**2 norma (no tail) = "<< Tn/2.5/2.5 << endl;
if(elso==1) fki << "# chi**2 norma (whole) = "  << Tnf/2.5/2.5 << endl;
  
if(elso==1) fki << "# Start JD = "<<tkezd <<endl;  
if(elso==1) fki << "# First counted day = "<<tmin <<endl; 
if(elso==1) fki << "# Time where the plateu ends (day) = "<<tfark/86400 <<endl; 
if(elso==1) fki << "# Final epoch (day) = "<<tfark/86400 <<endl; 
if(elso==1) fki << "# Exponential density profile exponent = "<<a <<endl; 
//if(elso==1) fki << "# Initial magnetar rotational energy (erg) = "<<Ep <<endl; 
//if(elso==1) fki << "# Characteristic time scale of magnetar spin-down (d) = "<<tp <<endl; 
if(elso==1) fki << "# Expansion velocity mean priori (km/s) = "<< vref/1e5 <<endl;
if(elso==1) fki << "# Expansion velocity sigma priori (km/s) = "<< vrefstd/1e5 <<endl;
if(elso==1) fki << "# Nickel mass Ag function fit (0:MCMC fit): "<< nifit <<endl;
if(elso==1) fki << "# Metropolis-Hastings step parameter = "<< Beta <<endl;
if(elso==1) fki << "# Step parameter (exponental) decay = "<< Delta <<endl;
if(elso==1) fki << "# Minimum ionisation zone = "<< txmin <<endl;
if(elso==1) fki << "# \n";
if(elso==1) fki << "# ID 	 E*R0 [erg*cm]  M [M_sol]    Ek [erg]      Ep [erg]      tp [day]          v [km/s]       T0 [K]        M_Ni [M_sol]    Ag [d^2]       kappa [g/cm^2]    s    Tion [K]       STD [mag], chi^2; whole and without tail\n";


if(elso==1) fki2 << "# Marcow Chain Monte Carlo solutions for Supernova light curve model which are not accapted\n";
if(elso==1) fki2 << "# Start JD = "<<tkezd <<endl; 
if(elso==1) fki2 << "# ID 	 E*R0 [erg*cm]  M [M_sol]    Ek [erg]      Ep [erg]      tp [day]          v [km/s]       T0 [K]        M_Ni [M_sol]    Ag [d^2]       kappa [g/cm^2]    s    Tion [K]       STD [mag], chi^2; whole and without tail\n";
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


if(niag!=1) if(elso==0){
	tER0=read2(sor,2);
	//tR0=read2(sor,2);
	tM=read2(sor,3)*Msol;
	tEk=read2(sor,4);
	//tE=read2(sor,5);
	tEp=read2(sor,5);
	ttp=read2(sor,6)*86400;
	tv=read2(sor,7)*1e5;
	tT0=read2(sor,8);
	tMni=read2(sor,9)*Msol;
	tAg=read2(sor,10)*(day*day);
	tkappa=read2(sor,11);
	ts=read2(sor,12);
	tTion=read2(sor,13);
	tszorasf=read2(sor,14)/2.5;
	szorash=tszoras=read2(sor,15)/2.5;
	tkhinegyzetf=read2(sor,16);
	khinegyzeth=tkhinegyzet=read2(sor,17);
	
	if(tEp<=0 && dEp!=0) dEp=0;
	if(ttp<=0 && dtp!=0) dtp=0;
}




for( ( elso || niag!=0 ) ?ciklus=0 :ciklus=read2(sor,1)+1  ;ciklus<=int(NMC);ciklus++){

   t=0; tprev=0; itt=0; lep=0; 
   if(niag!=1) tilt=0;
   xmin=txmin;

	
   // Ni Ag search mode	
   if(niag==1){
   	
     // first value, fixed	
     if(ciklus==0){
     	Mni=tMni=sqrt(Mnimin*Mnimax);
   	    Ag=tAg=sqrt(Agmin*Agmax);
     }
   
   
     // This generates the random values
     if(ciklus>0){
				 
				 do{ Mni=tMni*pow(10.,RandomDoubleGauss(2*Beta*sqrt(2.5*tszoras),0) );   }         while(Mni<Mnimin || Mni>Mnimax);
				 do{ Ag=tAg*pow(10.,RandomDoubleGauss(2*Beta*sqrt(2.5*tszoras),0) );   }           while(Ag<Agmin   || Ag>Agmax );
				     
				 }
   
	}
   	
	
			 
				                                              							  
   if(niag!=1) if(ciklus!=0){									/* Acceptable region ! */
   				 // Smaller region, smooth in lin
   			   	 //do{ R0=tR0+dR0*RandomDoubleGauss(1.,0);   }  while(R0<2e12 || R0>20e12);   
   	  			 //do{ M=tM+dM*RandomDoubleGauss(1.,0);   }     while(M<6e33 || M>60e33);     
				 //do{ Ek=tEk+dEk*RandomDoubleGauss(1.,0); }    while(Ek<2e50 || Ek>20e50);    
				 //do{ E=tE+dE*RandomDoubleGauss(1.,0);   }     while(E<2e50 || E>20e50);
				 
				 // Bigger region, smooth in log
				 //   *sqrt(2.5*tszoras)
				 //do{ R0=tR0*pow(10.,RandomDoubleGauss(Beta,0) );   }       while(R0<R0min || R0>R0max);
   	  			 //do{ M=tM*pow(10.,RandomDoubleGauss(Beta,0) );   }         while(M<Mmin || M>Mmax);
				 //do{ Ek=tEk*pow(10.,RandomDoubleGauss(Beta,0) ); }         while(Ek<Ekmin || Ek>Ekmax);
				 //do{ E=tE*pow(10.,RandomDoubleGauss(Beta,0) );   }         while(E<Emin || E>Emax);
				 
				 //do{ R0=tR0*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras) ,0) );   }       while(R0<R0min || R0>R0max);
   	  			 //do{ M=tM*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras) ,0) );   }         while(M<Mmin || M>Mmax);
				 //do{ Ek=tEk*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras) ,0) ); }         while(Ek<Ekmin || Ek>Ekmax);
				 //do{ E=tE*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras) ,0) );   }         while(E<Emin || E>Emax);
				 
				 do{ ER0=tER0*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) );   }     while(ER0<ER0min || ER0>ER0max); 
				 //do{ R0=tR0*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) );   }       while(R0<R0min || R0>R0max);
   	  			 do{ M=tM*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) );   }         while(M<Mmin || M>Mmax);
				 do{ Ek=tEk*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) ); }         while(Ek<Ekmin || Ek>Ekmax);
				 //do{ E=tE*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) );   }         while(E<Emin || E>Emax);
				 
				 if(dEp!=0) do{ Ep=tEp*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) ); }         while(Ep<Epmin || Ep>Epmax);
				 if(dtp!=0) do{ tp=ttp*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) ); }         while(tp<tpmin || tp>tpmax);
				 
				 if(dkappa!=0) do{ kappa=tkappa+dkappa*RandomDoubleGauss(1.,0);   }     while(kappa<kappamin || kappa>kappamax); 
				 if(ds!=0)     do{ s=ts+ds*RandomDoubleGauss(1.,0);   }    			  	while(s<smin || s>=smax ); 
				 if(dTion!=0)  do{ Tion=tTion+dTion*RandomDoubleGauss(1.,0);   }        while(Tion<Tionmin || Tion>Tionmax ); 
				 
				 //if(nifit==0) do{ Mni=tMni*pow(10.,RandomDoubleGauss(Beta,0) );   }         while(Mni<Mnimin || Mni>Mnimax);
				 //if(nifit==0) do{ Mni=tMni*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras) ,0) );   }         while(Mni<Mnimin || Mni>Mnimax);
				 if(nifit==0) do{ Mni=tMni*pow(10.,RandomDoubleGauss(Beta*sqrt(2.5*tszoras)*(exp((-1.*ciklus/Delta)+1.)+1) ,0) );   }         while(Mni<Mnimin || Mni>Mnimax);
				     
				 }
				 
				 
				 E=Ek;
				 R0=ER0/E;
				 
		
	/*skip=0;
	
				 if(R0<R0min || R0>R0max) skip=1;
   	  			 if(M<Mmin || M>Mmax) skip=1;
				 if(Ek<Ekmin || Ek>Ekmax) skip=1;
				 if(E<Emin || E>Emax) skip=1;
				 
				 if(dkappa!=0) if(kappa<kappamin || kappa>kappamax) skip=1; 
				 if(ds!=0)     if(s<smin || s>=smax ) skip=1; 
				 if(dTion!=0)  if(Tion<Tionmin || Tion>Tionmax ) skip=1; 
				 
				 if(nifit==0) if(Mni<Mnimin || Mni>Mnimax) skip=1;
				 */
				 


   szoras=0; szorasf=0;
   
   
  if(niag!=1){
  
  T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25);
   

  t=0.0; F=1.0;
  xh=1.0; xi=1.0; 
  Xni=1.0; Xco=0.0;
  //Q=1.6e13*Z/A*pow(Z,1.333333);		/*Ionization energy release*/
  Q = 3.22e13;
  
  Ith=Ith_int(1.0);
  IM=IM_int(1.0,a,2.0);
  
  
  // xmin 0.1 here, because here needs the real value !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xmin=0.1;
  if(s==0.0) f=IM_int(1.0,a,2.0);
  else f=(3.0*pow(xmin,s)-s*pow(xmin,3.0))/(3.0*(3.0-s));
  
  if(s==0.0) g1=IM_int(1.0,a,4.0);
  else g1=(5.0*pow(xmin,s)-s*pow(xmin,5.0))/(5.0*(5.0-s));
  xmin=txmin;
  
  v=sqrt(2.0*Ek*f/(g1*M));	
  
  rho=M/(4.0*M_PI*pow(R0,3)*f);
  
  td=3.0*kappa*rho*R0*R0/(pow(M_PI,2.0)*c);
  th=R0/v;
    
  Eth=4.0*M_PI*pow(R0,3)*7.57e-15*pow(T0,4.0)*Ith;
  //t0=pow(2.0*th*td,0.5);	    // obsolete
	

  // Ag and Ni compution (depends on the other parameters)      
  Ag=kgamma*M/(4.0*M_PI*f*v*v);   
  if(nifit!=0) Mni=Msol*fMni(Ag/86400./86400.);

	   
  p1=Eni*Mni*tni/Eth;
  p2=tni/td;
  p3=tni/Eth;
  p4=tni/tco;
  p5=Eco/Eni; 
  }
 
	while(t<=tmax)
	{ 
	  opac= 1.-exp(-Ag/(t*t));
	  opac1=1.-exp(- kpoz*Ag/kgamma /(t*t));
  	  
	  R=R0+v*t;	   
	  sigma=R/R0;				  
	  
	  
	  Tc=T0*pow(F,0.25)/sigma;
	  //if(xi>=xmin && Tc>Tion){ xi=temp(Tc,xh); }
	  //else xi=xmin; 
	  
	  
	  // 4th grade Runge-Kutta
	  if(tilt==0){
	  dF1=Frk(t, F, 1.);
	  segedvalt3=segedvalt1; 
	  dxi=segedvalt2*1./6.;
	  
	  dF2=Frk(t+0.5*dt, F+0.5*dF1*dt/tni, 1.5);
	  dxi=dxi+segedvalt2*2./6.;
	  
	  dF3=Frk(t+0.5*dt, F+0.5*dF2*dt/tni, 1.5);
	  dxi=dxi+segedvalt2*2./6.;
	  
	  dF4=Frk(t+1.0*dt, F+1.0*dF3*dt/tni, 2.);
	  dxi=dxi+segedvalt2*1./6.;
	  
	  F=F+(dF1+2.*dF2+2.*dF3+dF4)*dt/(6.*tni);
	  }
	  
	  
	  xi=segedvalt3; 
	  
	  Xni=   exp( -t/tni );   
	  Xco= (-exp( -t/tni ) + exp( -t/tco )) / (-tni*(1./tco-1./tni));   	  	  
  	  g=Xni+p5*Xco;
  	  
  	  if(tp==0) Em=0;
  	  else Em=Ep/(tp*pow((1.0+t/tp),2.0));
  	  
	 	  	  	  	  
	  ri=xi*R;
	  dr=dxi*R; if(dr>0.0) dr=0.0;
	  
	  
	  Lion=0.;
	  if(tilt==0){  
	  if (s==0.0)Lion=-4.0*M_PI*rho*eta(xi,a)*pow(sigma,-3.0)*Q*ri*ri*dr/dt;
	  else       Lion=-4.0*M_PI*rho*theta(xi,s)*pow(sigma,-3.0)*Q*ri*ri*dr/dt; 
	  if (Lion<0) Lion=0.0;	 
	  L=xi*F*Eth*opac/td; 
	  }
	  if(tilt==1) L= Mni*opac*( Eni *Xni + Eco *Xco)  +  opac*Em;
	  
	  Lpoz=Mni*(Ekin+Ean*opac)*opac1 * exp(-t/tco)-exp(-t/tni);
	   	  	     	  
	  Lsn=L+Lion+Lpoz;
	  logL=log10(Lsn);  
	  

	  // scatter (here in log erg!)
	  while(Adat[itt][0]*day>tprev && Adat[itt][0]*day<=t){
      double y=interpol(Adat[itt][0]*day, tprev, logLh, t, logL);
	  if(t>=tmin*day && t<=tfark) szoras=szoras+ h( (y-Adat[itt][1])/Adat[itt][2] )   ;
	  if(t>=tmin*day && t<=tmax) szorasf=szorasf+ h( (y-Adat[itt][1])/Adat[itt][2] )   ;
	  //if(t>=tmin*day && t<=tmax) szoras=szoras+ h( (y-Adat[itt][1]) )   ;
	  itt++;
	  }

	  
	  	  
	  xh=xi;
	  logLh=logL;
	  tprev=t;
	  t=t+dt;
	  
      
      // dinamical dt (empirical)
      if(niag!=1) if(xi==xmin) lep++;
	  else lep=0;
	  if(niag!=1) if(lep> (Ek>M*M*M/4e50?1:40) ) tilt=1;  //2e49/Msol**3 = 2e49/8e99   above this line the numerical error is very big after plateu, so it is analitical instantly, else it analitical after 10 day
	  //if(lep> 400 ) tilt=1;   if xmin>0.2   higher numerical step required, but more stabil!
	  else tilt=0;
	  if(tilt==1) dt=40000.;
	  else dt=20000. * fxi(xi);
	  
	  // simplification, if the scatter is high before the end of the plateu, xmin becomes larger to fasten the run time
	  if(niag!=1) if(t/day>=tmin+50 && t/day<tmin+51 ) if(2.5*sqrt(szoras/Tnx) > szoraskuszob*2.5 && xi>txmin+0.1 ) xmin=txmin+0.1;
	  
	  
	} 	
  
  khinegyzet=szoras;
  khinegyzetf=szorasf;
  szoras=sqrt(szoras/Tn);      //std.dev between tmin and tfark; before tail
  szorasf=sqrt(szorasf/Tnf);   //whole std.dev
  
  
  
  
  accept=0;
  
  
  /* Metropolis–Hastings method */
  if(niag!=1){
  szorast=exp(-0.5* khinegyzet   +0.5* khinegyzeth ); 
  if(vrefstd>0 && vref>0) szorast=szorast*exp(-0.5*pow((v-vref)/vrefstd,2) +0.5*pow((tv-vref)/vrefstd,2));
  
  }
  else /*szorast=khinegyzeth/khinegyzet;*/ szorast=szorash/szoras;
  if(szorast>1) szorast=1;
  
   if(khinegyzet*0!=0 && khinegyzeth*0==0) szorast=-1;
   if(khinegyzeth*0!=0) szorast=1e100;
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(niag!=1){
    if(szorast > RandomDouble(0.,1.) || ciklus==0){ tR0=R0; tM=M; tEk=Ek; tE=E; tkappa=kappa; ts=s; tTion=Tion;   tszoras=szoras;   tMni=Mni;   tAg=Ag;   accept=1; 
    			 tszorasf=szorasf; tkhinegyzetf=khinegyzetf; tkhinegyzet=khinegyzet; tv=v; tT0=T0; tER0=ER0; tEp=Ep; ttp=tp; }  // Ths is only used for store
    /*if(ciklus==0){ tR0=R0; tM=M; tEk=Ek; tE=E; tkappa=kappa; ts=s; tTion=Tion;   tszoras=szoras;   tMni=Mni;   tAg=Ag;    accept=1;
	             tszorasf=szorasf; tkhinegyzetf=khinegyzetf; tkhinegyzet=khinegyzet; tv=v; tT0=T0; }*/  // Ths is only used for store
    } else {
      if(szorast > RandomDouble(0.,1.) ){ tszoras=szoras; tMni=Mni; tAg=Ag; }
      if(ciklus==0){ tszoras=szoras; tMni=Mni; tAg=Ag; }
      }
  
  szorash=tszoras;
  khinegyzeth=tkhinegyzet;
  if(tszoras*0!=0) tszoras=1;
  
  if(accept==1) acc++;
  accall++;
  
  
  //if(accept==1) fkiS=fki; else fkiS=fki2;
  
  
  
  
  
    if(niag!=1){
     //if(szorasf*0==0 && szoras*0==0) fprintf(fkiS,"%d \t %2.3e   %4.3f   %2.3e   %2.3e     %2.3e   %2.3e   %2.3e   %2.3e     %4.4f   %4.3f   %d      ",ciklus,R0,M/Msol,Ek,E,v,T0, Mni/Msol,Ag/(day*day) ,kappa,s,int(Tion));
     //if(szorasf*0==0 && szoras*0==0) fprintf(fkiS,"%f   %f      %f   %f\n",szorasf*2.5,szoras*2.5,khinegyzetf,khinegyzet);
	 if(szorasf*0==0 && szoras*0==0) ((accept==1)?outtemp:outtemp2) << setprecision(3) << ciklus 
	 << " \t " << scientific << ER0 << "   " << fixed << M/Msol << scientific << "   " << Ek << "   " << Ep   << "   " << fixed << tp/86400 
	 << "      " << int(v/1e5) << "   " << scientific << T0 << "   " <<  Mni/Msol << "   " << Ag/(day*day)
	 << "      " << fixed << kappa << "   " << s << "   " << int(Tion) 
	 << "      " << szorasf*2.5 << "   " << szoras*2.5 << "     " << scientific << khinegyzetf << "   " << khinegyzet << endl;
     
     if(accept==0){
     	//if(szorasf*0==0 && szoras*0==0) fprintf(fki,"%d \t %2.3e   %4.3f   %2.3e   %2.3e     %2.3e   %2.3e   %2.3e   %2.3e     %4.4f   %4.3f   %d      ",ciklus,tR0,tM/Msol,tEk,tE,tv,tT0, tMni/Msol,tAg/(day*day) ,tkappa,ts,int(tTion));
        //if(szorasf*0==0 && szoras*0==0) fprintf(fki,"%f   %f      %f   %f\n",tszorasf*2.5,tszoras*2.5,tkhinegyzetf,tkhinegyzet);
        if(szorasf*0==0 && szoras*0==0) outtemp << setprecision(3) << ciklus 
	 	<< " \t " << scientific << tER0 << "   " << fixed << tM/Msol << scientific << "   " << tEk << "   " << tEp   << "   " << fixed << ttp/86400 
	 	<< "      " << int(tv/1e5) << "   " << scientific << tT0 << "   " <<  tMni/Msol << "   " << tAg/(day*day)
	 	<< "      " << fixed << tkappa << "   " << ts << "   " << int(tTion) 
	 	<< "      " << tszorasf*2.5 << "   " << tszoras*2.5 << "     " << scientific << tkhinegyzetf << "   " << tkhinegyzet << endl;
	 }
     
     
   } else{
     //if(szoras*0==0) fprintf(fki,"%f    ",szoras*2.5);
     //if(szoras*0==0) fprintf(fki,"%4.5f   %2.3e\n",Mni/Msol,Ag/(86400.*86400.));
     
	 if(szoras*0==0) fki << setprecision(4) << scientific << szoras*2.5 << "    "
     << Mni/Msol << "   " << Ag/(86400.*86400.) << endl;
     }
  
  
  
 	 if(ciklus%(niag!=1?100:1000)==0){
  	  if(szoras<szoraskuszob && szoras>=0){printf("+ "); }
  	  printf("Ciklus: %d   ratio: %2.2f%c\n",ciklus,100.*acc/accall,'%');
  	  //if(acc/accall<0.2) Beta=Beta*0.99;   // numerical method to keep the accep ratio around 20%
  	  //if(acc/accall>0.25) Beta=Beta*1.01;
  	  acc=accall=0;
  	  if(niag!=1){
  	    fki  << outtemp.str();
  	    fki2 << outtemp2.str();
   	    outtemp.clear();
	    outtemp.str(std::string());
	    outtemp2.clear();
	    outtemp2.str(std::string());
	    }
     }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  
  
  
  }

  
  
    fki  << outtemp.rdbuf();
    fki2 << outtemp2.rdbuf();
    fki.close();
    fki2.close();
  //fclose(fki);
  //fclose(fki2);
  
  
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
// If Ni(Ag) fitted, this section makes the best function file, which the program reads
  
  
if(niag==1){
	
	//int nn=1000;
	double be;
	double be1;
	double bek;
	
	double a;
	double b;
	int tari;

	double tar[nn][3];
	
	FILE *fkik;
	FILE *fki;
	fki=fopen("rand-niag.out","r");
	rewind(fki);
	fkik=fopen("NiAgBest.txt","w");
	
	
	
a=log10(Agmin/86400./86400.); b=log10(Agmax/86400./86400.);
for(i=0;i<nn;i++){
	
	tar[i][0]=a+ i/double(nn) *(b-a);
	tar[i][1]=0;
	tar[i][2]=1000;
	
	
}

for(i=0;i<nn;i++) tar[i][0]=pow(10,tar[i][0]);






while (!feof(fki)){    
                   
				   fscanf(fki,"%lf %lf %lf",&bek,&be,&be1);
                              
                   tari=inverztar( log10(be1) ,a,b,nn);    
				   
				   if(bek<tar[tari][2]){ tar[tari][2]=bek;  tar[tari][1]=be;  }
                              				   
				  
				   }
			



printf("Writing best fits\n");
for(i=0;i<nn;i++){

//printf("%e %f %f\n",tar[i][0],tar[i][1],tar[i][2]);
fprintf(fkik,"%e %e %f\n",tar[i][0],tar[i][1],tar[i][2]);
 
}
	

fclose(fki);
fclose(fkik);

	
}


  
  printf("done!\n"); 
  //getchar();
  return 0;

}




