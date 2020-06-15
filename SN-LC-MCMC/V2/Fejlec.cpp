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

#define Msol 1.989e33		/*Solar mass (g)*/



using namespace std;




//double xmin=0.2;		/*Minimum ionization radius*/	
double txmin=0.2;		/*saved Minimum ionization radius*/


double Tion;
//double E,Mni,Ek,F,L,T0,T,Tc,Tion,Lion,Ith,IM,INi,Eth,Lsn;
//double t,td,th,tmax,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
double a;
double s;
double Ep;
//double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,a,Z,A,Ag,Ep,s;
double kappa;
double tp;
//double kappa,M,v,R0,R,rho,ri,dr,Q,tni,tco,tp,Em,Em1;
//double Lpoz,Ag1;
double tmin, tkezd, tfark;
double tmax;

double NMC;   // Number of MC loops

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


//double para=0;	  // parameter a for manual MNi(Ag)
//double parb=0;	  // parameter b


// Velocity priori (Gaussian) !!!!!!!!!!!!!!!!!!! (if 0, then unset) Hardcoded
double vref=0;  // mean, in [cm/s] e.g: 5e8
double vrefstd=0; // error, in [cm/s] e.g: 1e8




					double Beta=0.1;	  // Metropolis–Hastings step parameter, default 0.1, set in parameter file
					double Delta=30000;



//int niag=0; // niag mode?
int nifit=1; // fit ni; niag mode?





//space-hez felbontas a gridhez, 200 jo kb
// set in MCMC parameter file
int grid=150;
int whole=1;
int burnin=0;




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






 void defaultdata()			 /*Writeing default input file*/
 {
 
 FILE *f;

 f=fopen("parametersMC.inp","wt");	
 
   fprintf(f,"%d\t[Data time shift (day)]\n",0);   	  /*Measure time shift*/
   fprintf(f,"%d\t[Time where the fitting starts (day)]\n",0);    	  /*Time where the scatter count starts*/
   fprintf(f,"%d\t[Ionization/recombination temperature (K)(-: fit)]\n",5500); 	  /*Ionization/recombination temperature (K)*/
   //fprintf(f,"%lf\t[Initial nickel mass]\n",Mni);  	  /*Initial nickel mass (M_sun) */  
   fprintf(f,"%d\t[Exponential density profile exponent]\n",0);            /*Exponential density profile exponent*/
   fprintf(f,"%d\t[Power-low density profile exponent (-: fit)]\n",0);		  /*Power-low density profile exponent*/
   fprintf(f,"%1.1f\t[Thomson scattering opacity (cm^2/g)(-: fit)]\n",0.3);     	  /*Thomson scattering opacity (cm^2/g)*/
   fprintf(f,"%d\t[Initial magnetar rotational energy (erg)(-: fit)]\n",0);		  /*Initial magnetar rotational energy (erg)*/
   fprintf(f,"%d\t[Characteristic time scale of magnetar spin-down (day)(-: fit)]\n",0);		  /*Characteristic time scale of magnetar spin-down (d)*/
   //fprintf(f,"%lf\t[Gamma-leak]\n",Ag);		  /*Gamma-leak (d^2)*/
   fprintf(f,"%d\t[Final epoch (day)]\n",300);	  /*Final epoch (day)*/
   fprintf(f,"%d\t[Number of MC loops]\n",300000);	  /*Number of MC loops*/
   fprintf(f,"%d\t[Time where the plateu ends (day)]\n\n",100);    	  /*Time where the plateu ends*/
   
   fprintf(f,"%1.0e\t[MIN Initial Thermal energy*radius (erg*cm)]\n",0.1e63);   	  /*ER0min*/
   fprintf(f,"%1.0e\t[MAX Initial Thermal energy*radius (erg*cm)]\n",300e63);   	  /*ER0max*/
   fprintf(f,"%d\t[MIN Ejected mass (Msol)]\n",3);   	  /*Mmin*/
   fprintf(f,"%d\t[MAX Ejected mass (Msol)]\n",60);   	  /*Mmax*/
   fprintf(f,"%1.1f\t[MIN Kinetic energy (foe)]\n",0.1);   	  /*Ekmin*/
   fprintf(f,"%d\t[MAX Kinetic energy (foe)]\n",30);   	      /*Ekmax*/
   fprintf(f,"%1.3f\t[MIN Magnetar energy (foe)]\n",0.001);   	  /*Epmin*/
   fprintf(f,"%d\t[MAX Magnetar energy (foe)]\n",30);   	      /*Epmax*/
   fprintf(f,"%1.2f\t[MIN Magnetar spin down (day)]\n",0.01);   	  /*tpmin*/
   fprintf(f,"%d\t[MAX Magnetar spin down (day)]\n",100);   	      /*tpmax*/
   fprintf(f,"%1.2f\t[MIN Thomson scattering opacity (cm^2/g)]\n",0.05);   	  /*kappamin*/
   fprintf(f,"%1.1f\t[MAX Thomson scattering opacity (cm^2/g)]\n",0.4);   	  /*kappamax*/
   fprintf(f,"%d\t[MIN Power-low density profile exponent]\n",0);   	  /*smin*/
   fprintf(f,"%d\t[MAX Power-low density profile exponent]\n",3);   	  /*smax*/
   fprintf(f,"%d\t[MIN Ionization/recombination temperature (K)]\n",0);   	      /*Tionmin*/
   fprintf(f,"%d\t[MAX Ionization/recombination temperature (K)]\n",20000);   	  /*Tionmax*/
   fprintf(f,"%1.0e\t[MIN Initial nickel mass (Msol)]\n",0.0005);   	  /*Mnimin*/
   fprintf(f,"%1.0e\t[MAX Initial nickel mass (Msol)]\n",0.5);   	      /*Mnimax*/
   fprintf(f,"%1.0e\t[MIN Gamma-leak (day^2) in NiAg]\n",1e3);   	  /*Agmin*/
   fprintf(f,"%1.0e\t[MAX Gamma-leak (day^2) in NiAg]\n",1e7);   	  /*Agmax*/
   fprintf(f,"%d\t[Expansion velocity mean priori (km/s)(0: not fit)]\n",0);   	  /*vref*/
   fprintf(f,"%d\t[Expansion velocity sigma priori (km/s)(gauss)]\n",0);   	  /*vrefstd*/
   fprintf(f,"%d\t[Determinate the nickel mass Ag function? (binary)]\n",0);	  /*nifit*/
   fprintf(f,"%1.1f\t[Metropolis–Hastings step parameter, default 0.1]\n",0.1);	  /*Beta*/
   fprintf(f,"%d\t[Step parameter (exponental) decay (loop number), default 30000]\n",30000);	  /*Delta*/
   fprintf(f,"%1.1f\t[Minimum ionisation zone (between 0.1 and 1)]\n\n",0.2);	  /*txmin*/
   //fprintf(f,"%1.1f\t[Threshold for szoras, technical parameter, default 0.2]\n\n",0.2);	  /*szoraskuszob*/
   
   fprintf(f,"%d\t[Whole (W) =1 or No Tail (T) =0 used?]\n",1);	  /*whole*/
   fprintf(f,"%d\t[Plotting grid number]\n",150);	  /*grid*/
   fprintf(f,"%d\t[burn in (steps)(=Step decay suggested)]\n",50000);	  /*burnin*/
   


 fclose(f);
 	
 }
 
 
 
void data()			 /*Reading input file*/
 {
 //FILE *f;
 //f=fopen("parametersMC.inp","rt");

 if( !f.is_open() ){printf("No parameter file! Generating a default! Exiting!\n"); defaultdata(); exit(1);}


   
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
   
   ER0min=read();
   ER0max=read();
   Mmin=read();  
   Mmax=read(); 	 
   Ekmin=read();  
   Ekmax=read(); 	 
   Epmin=read(); 
   Epmax=read(); 
   tpmin=read(); 
   tpmax=read(); 	
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
   //szoraskuszob=read();
   read();
   
   whole=int(read());
   grid=int(read()); //nn
   burnin=int(read());


 //++++++++++++++++++++++++++++++++++++++++++++
 
 if(s==3.0 || s==5.0){printf("Parameter error! The n = 3.0 and n = 5.0 are forbidden!\n"); exit(1);}
 
 /* Default acceptable region set if the input is wrong (between physically possible range) */
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
 //if(szoraskuszob<=0){ szoraskuszob=0.2; }
 
 if(ER0min<=0  || ER0max<=ER0min)    { ER0min=0.1e63;   ER0max=300e63; }
 if(Epmin<=0   || Epmax<=Epmin)      { Epmin=0.001;   Epmax=30; }
 if(tpmin<=0   || tpmax<=tpmin)      { tpmin=0.01;   tpmax=100; }
 
 if(grid<=0) grid=80;
 
 //para=Mni;
 //parb=Ag;
 
 
 
 //xmin=txmin;
 
 // Positive: fixed. Negative: MC
 /*
 if(kappa<0){ dkappa=Beta*0.2; kappa=0.5*(kappamin+kappamax); }
 if(s<0){ ds=Beta*1.5; s=0.5*(smin+smax); }
 if(Tion<0){ dTion=Beta*10000; Tion=0.5*(Tionmin+Tionmax); }
 
 if(Ep<0) dEp=1; else dEp=0;
 if(tp<0) dtp=1; else dtp=0;
 */
 
 //szoraskuszob=0.2;

 //M=Msol*M;	 		
 //Mni=Msol*Mni;
 //E=1e51*E;
 //Ek=1e51*Ek;
 Ep=1e51*Ep;
 tp=tp*86400.0;
 //Ag=Ag*pow(86400.0,2);		
 //T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 
 
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
 
 
 Epmin=1e51*Epmin;
 Epmax=1e51*Epmax;
 tpmin=86400*tpmin;
 tpmax=86400*tpmax;
 
 
 if(nifit!=0 && Ep!=0 && tp!=0){ printf("Magnetar is not compatible with NiAg fitting. Setting NiAg fit to 0."); nifit=0;
 }
 
 

 //fclose(f);
 f.close();
 }



//int main(){ return 0; }


    

