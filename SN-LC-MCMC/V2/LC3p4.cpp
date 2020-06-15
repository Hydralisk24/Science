#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream> // double -> string-hez kell

#define day 86400.0		/*1 day in seconds*/
//#define xmin 0.1		/*Minimum ionization radius*/	
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


double dt=200.0;		/*Time step (s)*/
double segedvalt1=0;    /*Assistant vairables*/
double segedvalt2=0;
double segedvalt3=0;
//double E,Mni,Ek,F,L,T0,T,Tc,Tion,Lion,Ith,IM,INi,Eth,Lsn;
//double t,td,th,tmax,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
//double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,a,Z,A,Ag,Ep,s;
//double kappa,M,v,R0,R,rho,ri,dr,Q,tni,tco,tp,Em,Em1;
//double Lpoz,Ag1;

double E,Mni,Ek,F,L,T0,T,Tc,Lion,Ith,IM,INi,Eth,Lsn;
double t,td,th,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,Z,A,Ag;
double M,v,R0,R,rho,ri,dr,Q,tni,tco,Em,Em1;
double Lpoz,Ag1;

double ER0;

double xmin=0.1;		/*Minimum ionization radius*/	







			// fast running (but with bigger numerical errors) (this used in the mcmc)
			int fastrun=1;
	
	
	
	
	
	
	
			
// Data checking?
int csekk=0;  // defined by arguments!!!!
			
			
						
			
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
//std::ifstream f("parameters.inp"); 
std::ifstream f; 
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
			
			
			
			
			
*/			



 void defaultdataA()			 /*Writeing default input file*/
 {
 
 FILE *f;

 f=fopen("parameters.inp","wt");	
 
   
   fprintf(f,"%1.0e\t[Initial radius (cm)]\n",1e13);   	  /*Initial radius (cm)*/
   fprintf(f,"%1.1f\t[Ejecta mass (Msol)]\n",7.5);    	  /*Ejecta mass (M_sun)*/
   fprintf(f,"%d\t[Ionization/recombination temperature (K)]\n",5500); 	  /*Ionization/recombination temperature (K)*/
   fprintf(f,"%1.1f\t[Initial nickel mass (Msol)]\n",0.1);  	  /*Initial nickel mass (M_sun) */
   fprintf(f,"%1.1f\t[Initial kinetic energy (1e51 erg)]\n",1.2); 	  /*Initial kinetic energy (1e51 erg)*/ 
   fprintf(f,"%1.1f\t[Intital thermal energy (1e51 erg)]\n",1.2); 		  /*Intital thermal energy (1e51 erg)*/    
   fprintf(f,"%d\t[Exponential density profile exponent]\n",0);            /*Exponential density profile exponent*/
   fprintf(f,"%d\t[Power-low density profile exponent]\n",0);		  /*Power-low density profile exponent*/
   fprintf(f,"%1.1f\t[Thomson scattering opacity (cm^2/g)]\n",0.3);     	  /*Thomson scattering opacity (cm^2/g)*/
   fprintf(f,"%d\t[Initial magnetar rotational energy (1e51 erg)]\n",0);		  /*Initial magnetar rotational energy (erg)*/
   fprintf(f,"%d\t[Characteristic time scale of magnetar spin-down (day)]\n",0);		  /*Characteristic time scale of magnetar spin-down (d)*/
   fprintf(f,"%1.0e\t[Unused Gamma-leak (day^2)]\n",1e6);		  /*Gamma-leak (d^2)*/
   fprintf(f,"%d\t[Final epoch (day)]\n",300);	  /*Final epoch (day)*/
   
   


 fclose(f);
 	
 }
 

void dataA()			 /*Reading input file*/
 {
 //FILE *f;
 
 //std::ifstream f; 
 //f.open("parameters.inp", std::ifstream::in);
 f.close();
 f.clear();
 f.open("parameters.inp");
 //f=fopen("parameters.inp","rt");
 if( !f.is_open() ){printf("No standard parameter file! Generating a default! Exiting!\n"); defaultdataA(); exit(1);}


   //fscanf(f,"%lf",&R0);   	  /*Initial radius (cm)*/
   //fscanf(f,"%lf",&M);    	  /*Ejecta mass (M_sun)*/
   //fscanf(f,"%lf",&Tion); 	  /*Ionization/recombination temperature (K)*/
   //fscanf(f,"%lf",&Mni);  	  /*Initial nickel mass (M_sun) */
   //fscanf(f,"%lf",&Ek); 	  /*Initial kinetic energy (1e51 erg)*/ 
   //fscanf(f,"%lf",&E); 		  /*Intital thermal energy (1e51 erg)*/    
   //fscanf(f,"%lf",&a);            /*Exponential density profile exponent*/
   //fscanf(f,"%lf",&s);		  /*Power-low density profile exponent*/
   //fscanf(f,"%lf",&kappa);     	  /*Thomson scattering opacity (cm^2/g)*/
   //fscanf(f,"%lf",&Ep);		  /*Initial magnetar rotational energy (erg)*/
   //fscanf(f,"%lf",&tp);		  /*Characteristic time scale of magnetar spin-down (d)*/
   //fscanf(f,"%lf",&Ag1);		  /*Gamma-leak (d^2)*/
   //fscanf(f,"%lf",&tmax);	  /*Final epoch (day)*/
   
   R0=read();
   M=read();
   Tion=read();
   Mni=read();
   Ek=read();
   E=read();
   a=read();
   s=read();
   kappa=read();
   Ep=read();
   tp=read();
   Ag1=read();
   tmax=read(); 


 M=Msol*M;	 		
 Mni=Msol*Mni;
 E=1e51*E;
 Ek=1e51*Ek;
 Ep=1e51*Ep;
 tp=tp*86400.0;
 Ag1=Ag1*pow(86400.0,2);		
 T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 

 if(s==3.0 || s==5.0){printf("Parameter error! The n = 3.0 and n = 5.0 are forbidden!\n"); exit(1);}

 //fclose(f);
  f.close();
 }
 
 
 
 void data2()			 /*Reading input file for check*/
 {
 	
 //FILE *f;
 FILE *fg;
 
 //f.open("parametersMC.inp", std::ifstream::in);
 //f=fopen("parametersMC.inp","rt");
 if( !f.is_open() ){printf("No MCMC parameter file! Exiting!\n"); exit(1);}
 fg=fopen("par.inp","rt");
 if(fg==NULL){printf("No checking parameter file (par.inp)! Exiting!\n"); exit(1);}
 
 double temp;
 
 data();

   //fscanf(fg,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&temp,&R0,&M,&Ek,&E,&temp,&temp,&Mni,&temp,&kappa,&s,&Tion);
   fscanf(fg,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&temp,&ER0,&M,&Ek,&Ep,&tp,&temp,&temp,&Mni,&temp,&kappa,&s,&Tion);
   
   //fscanf(f,"%lf",&temp);   	  /*Measure time shift*/
   //fscanf(f,"%lf",&temp);    	  /*Time where the scatter count starts*/
   //fscanf(f,"%lf",&temp); 	  /*Ionization/recombination temperature (K)*/
   //fscanf(f,"%lf",&temp);  	  /*Initial nickel mass (M_sun) */  
   //fscanf(f,"%lf",&a);            /*Exponential density profile exponent*/
   //fscanf(f,"%lf",&temp);		  /*Power-low density profile exponent*/
   //fscanf(f,"%lf",&temp);     	  /*Thomson scattering opacity (cm^2/g)*/
   //fscanf(f,"%lf",&Ep);		  /*Initial magnetar rotational energy (erg)*/
   //fscanf(f,"%lf",&tp);		  /*Characteristic time scale of magnetar spin-down (d)*/
   //fscanf(f,"%lf",&Ag1);		  /*Gamma-leak (d^2)*/
   //fscanf(f,"%lf",&tmax);	  /*Final epoch (day)*/
   //fscanf(f,"%lf",&temp);	  /*Number of MC loops*/
   
   
   /*
   temp=read();
   temp=read();
   temp=read();
   //temp=read();
   a=read();
   temp=read();
   temp=read();
   temp=read();
   temp=read();
   //Ep=read();
   //tp=read();
   //Ag1=read();
   tmax=read();
   //temp=read();
   //temp=read(); 
   
   */  
   
   
   E=Ek;
   R0=ER0/E;


   

 M=Msol*M;	 		
 Mni=Msol*Mni;
 //E=1e51*E;
 //Ek=1e51*Ek;
 //Ep=1e51*Ep;
 //tp=tp*86400.0;
 Ag1=Ag1*pow(86400.0,2);		
 T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 
 tp=tp*86400;

 if(s==3.0 || s==5.0){printf("Parameter error! The n = 3.0 and n = 5.0 are forbidden!\n"); exit(1);}
 


 //fclose(f);
 fclose(fg);
 }



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


	
	
	
double temp(double T1,double y0)	/*Find ionization radius short*/
	   {
	   double y,dy,zone;
	   double a,b;
	   a=xmin; b=y0;
	   if(fastrun!=1) dy=1e-10;      // +++
	   else dy=4e-6;

	   					
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



double IM_int(double b, double a1, double a2)
	   {
	   double sum=0.0,x=0.0, dx=1e-6;    // +++
	   if(fastrun==1) dx=1e-3;
	   
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
	   double sum=0.0,x=0.0,dx=1e-6;     // +++
	   if(fastrun==1) dx=1e-3;

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
	 
	  	
	  //Xni and Xco analitical
	  Xni=  exp( -t/tni );   
	  Xco= (-exp( -t/tni ) + exp( -t/tco )) / (-tni*(1./tco-1./tni));     	  	  
  	  g=Xni+p5*Xco;	  	  
	  
	  /*Runge-Kutta method*/	     
	  dF=sigma/pow(xil,3.0)*(p1*g-p2*F*xil-2.0*F*pow(xil,2.0)*dxi*tni/(sigma*dt)+p3*Em);
	  
	  return dF;
	  }
	  
	  
	  
	 
double fxi(double x){
	   double y=1;
	   
	   if(x==xmin) return y;
	   if(x<0.4) y=0.01+(x-0.1)*(x-0.1)*11;
	   
	   return y;
	   } 
	   
	   







int main(int argc, char **argv)
{
	csekk=argc-1; if(csekk>1) csekk=1; if(csekk<0) csekk=0;
	//csekk=1;
	//file: no csekk, file 1: csekk
	//argc: number of arguments (default 1), file 1 1 -> argc=3

int j=0,i=0, lep=0;
double dF1,dF2,dF3,dF4;
double opac,opac1,f,g1;
int tilt=0;

tni=8.8*day;				/*Ni decay time (s)*/
tco=111.3*day;				/*Co decay time (s)*/
Z=1.0;					/*Average atomic charge in ejecta*/
A=1.0;					/*Average atomic mass in ejecta*/

if(fastrun!=1) printf("Slower precise run.\n"); else printf("Fast run.\n");
if(csekk!=1) printf("Standard mode\n\n"); else printf("Checking: MCMC parameter file used.\n\n");

FILE *fki;
if(csekk!=1) dataA(); else data2();
fki=fopen("kimenet.out","w");

if(fastrun!=1) xmin=0.1; else xmin=0.2;   // +++                !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tmax=tmax*day; 			/*Final epoch (s)*/
  t=0.0; F=1.0;
  xh=1.0; xi=1.0; 
  Xni=1.0; Xco=0.0;
  //Q=1.6e13*Z/A*pow(Z,1.333333);		/*Ionization energy release*/
  Q = 3.22e13;
  
  Ith=Ith_int(1.0);
  IM=IM_int(1.0,a,2.0);
  
  //printf("%1.9f\n%1.9f\n",Ith,IM);
  
  
  if(s==0.0) f=IM_int(1.0,a,2.0);
  else f=(3.0*pow(xmin,s)-s*pow(xmin,3.0))/(3.0*(3.0-s));
  
  if(s==0.0) g1=IM_int(1.0,a,4.0);
  else g1=(5.0*pow(xmin,s)-s*pow(xmin,5.0))/(5.0*(5.0-s));
  
  v=sqrt(2.0*Ek*f/(g1*M));	
  
  rho=M/(4.0*M_PI*pow(R0,3)*f);
  
  td=3.0*kappa*rho*R0*R0/(pow(M_PI,2.0)*c);
  th=R0/v;
    
  Eth=4.0*M_PI*pow(R0,3)*7.57e-15*pow(T0,4.0)*Ith;
  //t0=pow(2.0*th*td,0.5);	   // obsolete
  
  Ag=kgamma*M/(4.0*M_PI*f*v*v);     // +++
  
  fprintf(fki,"#dt=200s\n");
  fprintf(fki,"#Ag difference: %e (calculated/given)\n",Ag/Ag1);   
   fprintf(fki,"#Supernova light curve model based on Arnett & Fu ApJ 340, 396 (1989)\n");
   fprintf(fki,"#\n");
   fprintf(fki,"#Inital model parameters\n");
   fprintf(fki,"#R0 = %lg cm\n",R0);
   fprintf(fki,"#Mej = %lg M_sol\n",M/Msol);
   fprintf(fki,"#MNi = %lg M_sol\n",Mni/Msol);
   fprintf(fki,"#Eth = %lg erg\n",Eth);
   fprintf(fki,"#Ekin = %lg erg\n",Ek);
   fprintf(fki,"#Trec = %lg K\n",Tion);
   fprintf(fki,"#kappa = %lg cm^2/g\n",kappa);
   fprintf(fki,"#a = %lg \n",a);
   fprintf(fki,"#s = %lg \n",s);   
   fprintf(fki,"#Ep = %lg erg\n",Ep);
   fprintf(fki,"#tp= %lg d\n",tp/day);
   fprintf(fki,"#Ag = %lg d^2\n",Ag/(day*day)); 
   fprintf(fki,"#Ag+ = %lg d\n",kpoz*Ag/kgamma/(day*day));
   fprintf(fki,"#\n");
   fprintf(fki,"#Calaculated physical properties \n");
   fprintf(fki,"#v = %lg km/s\n", v/1e5);
   fprintf(fki,"#rho0 = %lg g/cm^3\n", rho);
   fprintf(fki,"#td = %lg d\n",td/day);
   fprintf(fki,"#th = %lg d\n",th/day);
   fprintf(fki,"#T0 = %lg K (init temp, not gamma leak)\n",T0);
   fprintf(fki,"#\n");
   fprintf(fki,"#t/day    Lsn    logLsn    MAG     Ldiff      Lrec 	  xi    F\n");
       
  p1=Eni*Mni*tni/Eth;
  p2=tni/td;
  p3=tni/Eth;
  p4=tni/tco;
  p5=Eco/Eni; 
 
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
	  //F=F+dF1*dt/(tni);
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
	  
	  
	  if(j==0){
	  if(logL*0==0) fprintf(fki,"%lf %lg %lg %lg %lg %lg %lg %f\n",t/day, Lsn, logL, -2.5*(logL-17.5-18), L, Lion, xi, F);
	  printf("t=%4.1f    dt=%4.1f   \n",t/day,dt);
	  }
	  j++;
	  if( fastrun!=1?j==500:j==5 ) j=0;     // +++
	  
	  	  
	  xh=xi;
	  t=t+dt;
	  
      
      // dinamical dt	 
	  if(fastrun!=1){
      if(xi==xmin) lep++;
	  else lep=0;
	  if(lep>10000) tilt=1;   // +++
	  else tilt=0;
	  if(tilt==1) dt=400.;    // +++
	  else dt=200.;           // +++
      }
      
      else{
      if(xi==xmin) lep++;
	  else lep=0;
	  if(lep>5) tilt=1;
	  else tilt=0;
	  if(tilt==1) dt=40000.;
	  else dt=20000. * fxi(xi);
	  }

	  
	} 	
  fclose(fki);

return 0;
}

