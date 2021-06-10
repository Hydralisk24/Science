#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define G 6.674e-11
#define sig 5.670367e-8




// sdss griz :0 (AB, F_f) vagy Johson BVRI :1 (Vega, F_l) bolometrikus integralas
int mod=1;


	  
	  
	  
// Fekete Test Sugarzas, F_f es F_l is   W/m^2/m --> erg/s/cm^2/A   es  erg/s/cm^2/Hz
// bemenet: Amstrong, Kelvin
float FTS(float l, float T){
	  float B;
	  l=l/1e10; //Amstron konverzio
	  
	  if(mod==0) B= 1000*      PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) * l*l/3e8;
	  if(mod==1) B= 1000/1e10* PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) ;
	  
	  return B;
	  }
	  
// Fekete Test Sugarzas, F_f es F_l is   W/m^2/m --> erg/s/cm^2/A   es  erg/s/cm^2/Hz    
// DOUBLE-val!!! Integralas-ra hasznalt!
// bemenet: meter, Kelvin
double FTSd(double l, float T){
	  double B;
	  
	  if(mod==0) B= 1000*      PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) * l*l/3e8;
	  if(mod==1) B= 1000/1e10* PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) ;
	  
	  return B;
	  }
	  

// FTS integralja, bemenet: Amstron es Kelvin
double FTSint(float s, float tig, float T){	  
	  if(T<=0) return 0;
	  
	  double v=0;
	  double a;
	  double ce;
	  double ck;
	  double ch;
	  double cn;
	  double i;
	
	//double dt=0.01e-9; //nm
	double dt=1e-9; //nm
	
	tig=tig/dt/1e10;
	s=s/dt/1e10;

	//RK ciklus
	for(i=s;i<tig;i++){

	   a=FTSd(i*dt,T);
	   ce=dt*a;

	   a=FTSd((i+0.5)*dt,T);
	   ck=dt*a;

   	   a=FTSd((i+0.5)*dt,T);
   	   ch=dt*a;

  	   a=FTSd((i+1.0)*dt,T);
  	   cn=dt*a;

  	   v=v+(1.0/6.0)*ce+(1.0/3.0)*ck+(1.0/3.0)*ch+(1.0/6.0)*cn;  // RungeKutta modszer
  	   
       }

return 1e10*v;  // mertek egyseg!
}



// erg/s/cm^2 ben jon ki a valtoszam 1000




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
	
	
	
	




// visszaadja azon tomb szamot, ami x erteknel van
// kell hozza az elso ertek (a) utolso ertek (b) es a lepes koz (dx)
// csak egyenletes dx lepes kozokre jo!!!!!!
int inverztar(double x, double n, double a, double b){
	double y;
	
	//if(dx==0) return 0;
	if(x>b) x=b;
	
	//tar[i][0]=a + i *(b-a)/double(nn);
	y = (x - a) * double(n)/(b-a);
	//y = (x - a) /dt;
	if(y<0) y=0;
	if(y>n-1) y=n-1;
	
	
	return int(y);
    }
	

int n=10000;	
double Adatcs[10000];
	
double intcosin(double a,double b){
	double y=0;
	
	y=Adatcs[inverztar(b,n,0,180/57.3)]-Adatcs[inverztar(a,n,0,180/57.3)];
	
	return y;
	
}














double fint(double x, double T){
	
	double y=1*cos(x)*sin(x);
	if(y<0) return 0;
	else return y;
}



// FTS integralja, bemenet: Amstron es Kelvin
void integral(float s, float tig, float T){	  
	  double v=0;
	  double a;
	  double ce;
	  double ck;
	  double ch;
	  double cn;
	  double i;

	//double dt=0.01e-9; //nm
	double dt=1e-5; //nm
	
	tig=tig/dt;
	s=s/dt;

	//RK ciklus
	for(i=s;i<tig;i++){

	   a=fint(i*dt,T);
	   ce=dt*a;

	   a=fint((i+0.5)*dt,T);
	   ck=dt*a;

   	   a=fint((i+0.5)*dt,T);
   	   ch=dt*a;

  	   a=fint((i+1.0)*dt,T);
  	   cn=dt*a;

  	   v=v+(1.0/6.0)*ce+(1.0/3.0)*ck+(1.0/3.0)*ch+(1.0/6.0)*cn;  // RungeKutta modszer
  	   
  	   Adatcs[int((i-s)*n/tig)]=v;
  	   
       }

Adatcs[n-1]=0.5;
//return v;  // mertek egyseg!
}
	
	
	
	




double h(double x){ return x*x; }


double eg(double x){ if(x>1) return x; else return 1; }


int main (){
	
	

	
srand(time(NULL));
printf("Integrating...\n");
integral(0, 180/57.3, 0);
printf("Integrating DONE\n\n");
    
//for(int i=0;i<=180;i++) printf("%d fok    ertek: %f     pos: %d\n",int(i),intcosin(0,i/57.3),  inverztar(i/57.3,n,0,180/57.3) );
//getchar();


double T;

double R;
double q;
double e;
double uQ;
double d;
double fi;

double M; //hold tomeg

double qB; // bolygo surusege
double MB; // bolygo tomege
double RB; // bolygo sugara

//double a;
//double M;
double P;
double dis=1;
double albedo=0.5;
double albedoB=0.5;
double Rsun,Tsun,Lsun,Msun;

double B=1;
double Tn;
double flux;
double fluxpart;
double fluxpartB;
double fluxpart1, fluxpart2;
double Mag;
double MagB;

double lmin;
double lmax;

double tlmax,slmax;

double Tlim=3000; // limit for temp



double Tbb=0;
double Tbbs=0;
double Ts=0;

double cc=1;

double Mc=20;
double Lkh=0; // Kelvin Heimholtz luminosity
double LB=0;

double TB=0; // temperature bolygo



FILE *f;
f=fopen("flux.txt","w");





d=1;
    



	//ELT METIS spektro
	//https://www.eso.org/public/teles-instr/elt/elt-instr/metis/
	lmin=139000; // um -> A
	lmax=140000; // um -> A
	
	//lmin=5000; lmax=6000; //johnson V
	//lmin=1000; lmax=1000000;
	
	dis=30; //pc
	dis=dis*3.26*1e13*1e5; // pc -> cm
	
	Rsun=696e6;
	//Tsun=5700;
	Msun=1;
	
	qB=1000; //Kg/m3  //samplingolva van, most ez a kezdo
	MB=500; //Mfold
	
	Mc=20; // maximum tomegarany
	
	albedo=0.5;  //hold albedo
	albedoB=0.5; //bolygo albedo
	
	Tlim=3000; // max limit for temperature
	//3000K
	//https://iopscience.iop.org/article/10.3847/1538-4357/ab9cba/pdf
	//https://ui.adsabs.harvard.edu/abs/2020ApJ...898..160E/abstract
	
	


tlmax=lmax;

MB=MB*6e24; //Mfold -> Kg
RB=pow(MB/qB/(PI*4./3),1./3);

Lsun=3.84e26*pow(Msun,3.5);
Msun=Msun*2e30;
//Lsun=4*PI*h(Rsun)*sig*pow(Tsun,4);
Tsun=pow(Lsun/sig/4/PI,0.25)/sqrt(Rsun);

//printf("%e  %e  %f  %f\n",Msun,Lsun,Tsun,Rsun/1000); exit(1);

    
/*    
B=a*pow(3*M/2/PI/q,-1./3) ;
a^3/P^2 = G*M / 4*PI

P^2/3 = a/(G*M / 4*PI)^1/3
B=a*(2/3*PI*q / M)^1/3

a/M^1/3 = P^2/3*(G/4PI)^1/3
B = P^2/3*(G/4PI)^1/3 * (2/3*PI*q)^1/3
*/



double sun=( 4*PI*h(Rsun)*FTSint(lmin, lmax, Tsun) );
//*1e4 meg kell szorozni, hogy erg/s legyen!!!  *1e3-al pedig hogy W
//sun/1000=3.84e26
 
//                             albedo helye
Tbbs=pow(Lsun/h(1*150e9)/4/PI*(1)  /4/sig,0.25);

 
double sR;
double sq;
double se;
double suQ;
double sP; 
double sd;
double sfi;
double sfluxpart=0;
double sfluxpartB=0;
double sMag=0;
double sMagB=0;
double sM=0;
double sMB=0;
double sqB=0;
 
double tR;
double tq;
double te;
double tuQ;
double tP; 
double td;
double tfi;
double tfluxpart=0;
double tfluxpartB=0;
double tM;

double tT;
double tTs;
double tMB;
double tqB;
 
double qq=0;
int accept=1;
int all=1;

double szorast=1;

double Plim=0;
double tPlim=1;
double Phill=25;
double tPhill=25;

double temp1=0;
double temp2=0;
double temp3=0;


	//allitandoak:
	double mcmctau=10;
	double mcmcmax=3;
	int lim=100000; //lim=1;
	double Rspotlim=13000e3;
 
 

tR=R=2*6378e3; //m
tq=q=5515; //KG/m3
te=e=0.1;
tuQ=uQ=pow(10,14.0);
tP=P=1*86400; //day, sec
td=d=1; //CSE/AU
tfi=fi=90;
tfluxpart=1e-300;
tfluxpartB=1e-300;

tMB=MB;
tqB=qB;


fprintf(f,"#ID  qq   flux (moon/star)[ppm]  (moon/planet)[1]  L(abs)[SI]  AB mag  Bol mag  \t RM[km]  qM[SI]  e  P[day]  log uQ  fi    Tbb  T  Ts  TB 	\t  d[AU]  MB[M_earth]  qB[SI]  LB(in band)[SI]    lmax\n");


for(int i=0;i<lim;i++){



  do{
  
  do{ lmax=tlmax*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(lmax<5000 || lmax>300000);
  lmin=0.99*lmax;
  sun=( 4*PI*h(Rsun)*FTSint(lmin, lmax, Tsun) );
  
  //do{ R=tR+10000e3*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(R<=100e3 || R>25000e3);
  //do{ q=tq+1000  *exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(q<=0 || q>30000);
  //do{ e=te+0.05  *exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(e<0 || e>1);
  //do{ uQ=tuQ*pow(10, 0.5*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(uQ<10e10 || uQ>10e22);
  //do{ P=tP+0.2*86400*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(P<0.25*86400 || P>50*86400);
  //do{ d=td+0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(d<=0.01 || d>50);
  
  do{ fi=tfi+10*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); } while(fi<0.5 || fi>180);
  //fi=180;
  
  
  //do{ R=tR*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(R<=500e3 || R>25000e3);
  //do{ q=tq*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(q<=100 || q>30000);
  //do{ e=te*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(e<0.0001 || e>1);

  //do{ P=tP*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(P<Plim*86400 || P>15*86400);
  //do{ d=td*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); } while(d<=0.01 || d>100);
  
  
  do{ MB=tMB*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) );
  	  qB=tqB+1000  *exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);
  	  RB=pow(MB/qB/(PI*4./3),1./3);
  	  d=td*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) );
  
  	  R=tR*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) );
  	  q=tq+1000  *exp(-qq/mcmctau)*RandomDoubleGauss(1.,0); 
      P=tP*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) ); 
      e=te*pow(10, 0.2*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0) );
	  temp1=0.0625/sqrt( qB/5500)/pow(1-e,1.5) *pow(1+R/RB,1.5);
	  temp2=0.0625/sqrt(3*q/5500)/pow(1-e,1.5);
	  M=4./3*PI*pow(R,3)*q;
	  (temp1>temp2?Plim=temp1:Plim=temp2); 
	  Phill=0.35*pow(d*150e9,1.5) * sqrt(4*PI*PI/G/3/Msun)/pow(1+e,1.5)/86400; 
	  } 
	  while(P<Plim*86400 || P>Phill*86400 || q<=2500 || q>10000 || e<0.0001 || e>1 || R<=100e3 || R>15000e3  || M>MB/Mc || MB>2*300*6e24 || MB<50*6e24 || qB<=500 || qB>10000 || d<=0.01 || d>100);

  //M=4./3*PI*pow(R,3)*q;
  
  //Felszinen a keringesi ido: P=gyok( 3PI/(G*q) ), avagy 1.5ora/gyok(q) [fold surusegben]
  //negyed akkora suru gazoriasnal a keringesi ido a felszinen 3 ora (1/8 nap), legyen a minimum ennek 2x
  //
  //a_roche = RB*kobgyok( 3*qB/q ) / (1-e)
  //a > RB*kobgyok( 3*qB/q ) / (1-e)
  //kobgyok( P^2 *G*MB/4/PI^2 ) > RB*kobgyok( 3*qB/q ) / (1-e)
  // P^2 *G*MB/4/PI^2  > RB^3*( 3*qB/q ) / (1-e)^3
  // P^2 *G/PI/3  > 3/q / (1-e)^3
  //sqrt( 3* 3*PI/(G*q) ) /(1-e)^1.5 < P
  
  
  //a^3/P^2 = G*MB/4/PI^2
  //a = kobgyok( P^2 *G*MB/4/PI^2 )
  //RB+R < kobgyok( P^2 *G*MB/4/PI^2 ) *(1-e)
  
  //1+R/RB < kobgyok( P^2 *G*MB/4/PI^2/RB^3 )  *(1-e)
  //1+R/RB < kobgyok( P^2 *G*MB/3/PI/(4/3*PI*RB^3) )  *(1-e)
  //1+R/RB < kobgyok( P^2 *G*qB/3/PI)  *(1-e)
  //(1+R/RB)^1.5 *sqrt(3*PI/(G*qB) ) /(1-e)^1.5 < P
  
  //RB = kobgyok( P^2 *G*MB/4/PI^2 )
  //RB^3 = P^2 *G*MB/4/PI^2 
  //RB^3*4/3*PI/MB = P^2 *G/3/PI 
  //1/qB = P^2 *G/3/PI 
  //1 = P^2 *G/3/PI *qB
  
  //r_hill = d * kobgyok(mB/3Mstar)  /(1+e)
  // a^3/P^2 = G/4/PI/PI*M   ->     P^2 = a^3/M * 4*PI*PI/G    ->    a^3/M  =  P^2 * G/4/PI/PI
  //r_hill/mB^1/3 = d * kobgyok(1/3Mstar)  /(1+e)
  //(P^2 * G/4/PI/PI)^1/3 = d * kobgyok(1/3Mstar)  /(1+e)
  //P^2/3 = d * kobgyok(4*PI*PI/3Mstar/G)  /(1+e)
  //P < P_hill = d^1.5 * gyok(4*PI*PI/3Mstar/G)   /(1+e)^1.5

  //kommenteld ki az uQ= reszt es be a break-et ogy visszakapd a fixedQ-t.
  //aproxximation for viscoelastic
  for(double x=250;;x=x-100*log(x/T)){
 
    if(x<1){ T=1; break; };
    uQ=pow(10,  -4.94e-14*R*R+9.61e-7*R+5.59 +0.017*x);
	//tidal T! 
    B = pow(P,2./3) *  pow(G/6*q,1./3);   
    Tn=  17.183*  sqrt(392.*pow(PI*G,5)/9747.)/sig * pow(R,5)*pow(q,4.5)/uQ  * e*e*pow(B,-7.5)    ;
    T=pow(    Tn    ,0.25);  
    
    if(x/T<1.01 && T/x<1.01) break;
    //break;
  
    }
  
  Tbb=Tbbs/sqrt(d)*pow(1-albedo,0.25);
  
  cc=h(sin(0.5*fi/57.3));
  Ts=pow(pow(Tbb,4)+pow(T,4)/cc,0.25);
  
  }while(Ts>Tlim || T<=0);
  
  //if(R>Rspotlim) fi=180;
  
  //printf("%1.2f\n",Phill); getchar();
  
  if(i==0){
   tT=T; tTs=Ts;
   tPlim=Plim;
   tPhill=Phill;
   }


  //relative!
  //fluxus arany a naphoz kepest!!!!
  //fluxpart=        4*PI*h(R)*FTSint(lmin, lmax, T)/sun;
  //fluxpart=( 0.9  *4*PI*h(R)*FTSint(lmin, lmax, Tbb/sqrt(d)+ca*T)  +0.1  *4*PI*h(R)*FTSint(lmin, lmax, Tbb/sqrt(d)+cb*T)  )/sun;
  //fluxpart=( 0.9  *4*PI*h(R)*FTSint(lmin, lmax, Tbb/sqrt(d)+ca*T)  +0.1  *4*PI*h(R)*FTSint(lmin, lmax, Tbb/sqrt(d)+cc*T)  )/sun;
  
  //fluxpart = 4*PI*h(R)*sig*Tn/3.84e26;
  temp1=FTSint(lmin, lmax, Ts);
  temp2=FTSint(lmin, lmax, Tbb);
  
  fluxpart=(intcosin(0,fi/57.3)*temp1 + intcosin(fi/57.3,180/57.3)*temp2 )*2*4*PI*h(R)/sun  +  albedo*h(R/d/150e9)/4;
  flux=(intcosin(0,fi/57.3)*pow(Ts,4) + intcosin(fi/57.3,180/57.3)*pow(Tbb,4) )*sig*4*PI*h(R)  +  Lsun*albedo*h(R/d/150e9)/4;
  
  Lkh=0.5*G*h(MB)/RB / (3.1e7 *3*1e10); //t Kelvin Heimholtz = 30 milliard ev
  // https://astronomy.stackexchange.com/questions/7970/what-is-the-long-term-fate-of-the-gas-giants
  // https://arxiv.org/abs/1405.3752
  // https://ui.adsabs.harvard.edu/abs/2014arXiv1405.3752G/abstract
  LB=Lkh+4*PI*h(RB)*sig*pow( Tbbs ,4) *(1-albedoB)/h(d);
  
  TB=pow( LB/sig/4/PI/h(RB) , 0.25) ;
  
  LB=4*PI*h(RB)*FTSint(lmin, lmax, TB);
  
  //sun*0.5*h(R/d)/4
  //                         albedo
  fluxpartB=fluxpart*sun/( sun*albedoB*h(RB/d/150e9)/4 + LB );   // 4*PI*h(RB)*temp2 + Lkh*1000
  //!!!!!!!!!!!! temp2-nel Tbb akkor igaz ha a hold es a bolygo albedoja megegyezik, kulonben mas a Tbb
  
  
  
  
  //MagB=4.75-2.5*log10(fluxpart*sun/1000/h(dis/3.26e19)/3.84e26);   //bolometric mag
  MagB=4.75-2.5*log10(flux/h(dis/3.26e19)/3.84e26);   //bolometric mag
  
  //Mag=-48.60-2.5*log10(   fluxpart*sun/(lmax-lmin) /3e18*h(lmax+lmin)/4   /h(dis)/4/PI   );  // AB mag
                         //erg/s / A          /  c/l/l               / distance(cm2)                 
  Mag=-48.60-2.5*log10(   fluxpart*  sun*1e4/(3e18/lmin-3e18/lmax)  /h(dis)/4/PI   );  // AB mag
  
  //AB forrasok:
  //http://mips.as.arizona.edu/~cnaw/sun.html
  //https://en.wikipedia.org/wiki/AB_magnitude
  




  // metropolis ratio
  szorast=fluxpart/tfluxpart;
  
  //szorast=szorast*exp(-0.5*pow((R-6000e3)/6000e3,2) +0.5*pow((tR-6000e3)/6000e3,2));
  szorast=szorast*(1-exp(-0.5*pow((R/15000e3-1)/0.125, 2 ) ) ) / ( 1-exp(-0.5*pow((tR/15000e3-1)/0.125, 2 ) ) );
  //szorast=szorast*exp(-0.5*pow((q-6000)/6000,2) +0.5*pow((tq-6000)/6000,2));
  szorast=szorast*(1-exp(-0.5*pow((Mc*M/MB-1)/0.25, 2 ) ) ) / ( 1-exp(-0.5*pow((Mc*M/MB-1)/0.25, 2 ) ) );
  
  
  //szorast=szorast*exp(-0.5*pow((e-0.0)/0.05,2) +0.5*pow((te-0.0)/0.05,2));
  //szorast=szorast*(1-exp(-0.5*pow((e-1)/0.4, 2 ) ) ) / ( 1-exp(-0.5*pow((te-1)/0.4, 2 ) ) );
  szorast=szorast*(1-exp(-0.5*pow((log10(e)-0)/.5, 2 ) ) ) / ( 1-exp(-0.5*pow((log10(te)-0)/.5, 2 ) ) );
  
  szorast=szorast*exp(-0.5*pow((log10(uQ)-15)/2,2) +0.5*pow((log10(tuQ)-15)/2,2));
  //szorast=szorast*exp(-0.5*pow((P-5*86400)/86400,2) +0.5*pow((tP-5*86400)/86400,2));
  //szorast=szorast*(1-exp(-(P-0)/86400 ) ) / ( 1-exp(-(tP-0)/86400 ) );
  szorast=szorast*(1-exp(-0.5*pow((P-Plim*86400)/86400, 2 ) ) ) / ( 1-exp(-0.5*pow((tP-tPlim*86400)/86400, 2 ) ) );
  
  szorast=szorast*(1-exp(-0.5*pow((d-0)/0.1, 2 ) ) ) / ( 1-exp(-0.5*pow((td-0)/0.1, 2 ) ) );
  
  szorast=szorast*exp(-0.5*pow((Ts-100)/(0.5*Tlim),2) +0.5*pow((tTs-100)/(0.5*Tlim),2));
  //szorast=szorast*(1-exp(-0.5*pow((T/Tlim-1)/0.35, 2 ) ) ) / ( 1-exp(-0.5*pow((tT/Tlim-1)/0.35, 2 ) ) );  if(T>Tlim) szorast=0;
  
  
  
  if(szorast*0!=0) szorast=-1;
  
  if( szorast >= RandomDouble(0.,1.) ){

	 				tlmax=lmax;
					 
					tfluxpart=fluxpart;    
	 				tfluxpartB=fluxpartB;
					accept++;
		 			tR=R; tq=q; te=e; tuQ=uQ; tP=P;
		 			td=d;
		 			//if(R<=Rspotlim) 
					tfi=fi;
		 			tT=T; tTs=Ts;
		 			tM=M;
		 			tPlim=Plim; tPhill=Phill;
		 			
		 			tMB=MB;
		 			tqB=qB;
		 			
		 			
	    			}
	    
  all++; // counting MCMC step

  if(i%100==0){
					//cout << "Ratio= " << 100.*accept/all <<"    qq= "<<qq<< endl;
					//printf("i:%d   ratio:%d   qq:%d     flux:%1.2f [ppm] \t R:%1.2f  q:%1.2f  e:%1.2f  t:%1.2f  log(uQ):%1.2f   d:%1.2f   fi:%1.2f \n",i,int(100.*accept/all),int(qq),tfluxpart*1e6,tR/1000,tq,te,tP/86400,log10(tuQ),td,tfi);
					printf("i:%d   ratio:%d   qq:%d \n",i,int(100.*accept/all),int(qq));
					
					if(1.*accept/all<0.1) qq++;
					if(1.*accept/all>0.5) qq--;
					if(qq<0) qq=0;
					if(qq>mcmctau*mcmcmax) qq=mcmctau*mcmcmax;
					if(1.*accept/all<0.04 && qq==mcmctau*mcmcmax){ //cout << "break at: " << i << endl; 
															lim=i;
															//break; 
															}									
					accept=0; all=0;
	    			}

  
  //legjobb ertekek
  if(sfluxpart<fluxpart){
  	slmax=lmax;
	sR=R; sq=q; se=e; suQ=uQ; sP=P;
  	sfluxpart=fluxpart; sfluxpartB=fluxpartB;
	sd=d; 
	sfi=fi;
	sMag=Mag; sMagB=MagB;
	
	sMB=MB;
	sqB=qB;
    }
  //          1   2     3       4     5      6      7           8      9     10      11    12      13       14    15      16    17            18    19      20     21      22
  fprintf(f,"%d  %d   %2.2f  %1.3f  %1.2e  %1.2f  %1.2f   \t  %1.1f  %1.1f  %1.5f  %1.2f  %1.2f  %1.1f    %1.2f  %1.2f  %1.2f  %1.2f 	\t  %1.3f  %1.2f  %1.1f  %1.2e   %1.2e\n",
             i,int(qq),fluxpart*1e6,fluxpartB,flux,Mag,MagB,R/1000,q,e,P/86400,log10(uQ),fi,Tbb,T,Ts,TB,d,MB/6e24,qB,LB/1000,lmax);
             //1   2       3            4      5    6   7    8     9 10   11     12      13 14 15 16 17 18   19   20   21     22
			//if(fluxpartB>1e5){	printf("wtf\n"); return 0;	}
  
  }
  
printf("best:\n");
printf("flux:%1.2f [ppm]  fluxB:%1.3f  Mag: %1.2f  %1.2f   R:%1.2f  q:%1.2f  e:%1.2f  t:%1.2f  log(uQ):%1.2f   fi:%1.2f \n",sfluxpart*1e6,fluxpartB,sMag,sMagB,sR/1000,sq,se,sP/86400,log10(suQ),sfi);
printf("d:%1.2f   MB:%1.2f   qB:%1.2f   lmax:%1.2e\n",sd,sMB/6e24,sqB,slmax);

fprintf(f,"#flux:%1.2f [ppm]  fluxB:%1.3f  Mag: %1.2f  %1.2f   R:%1.2f  q:%1.2f  e:%1.2f  t:%1.2f  log(uQ):%1.2f   fi:%1.2f \n",sfluxpart*1e6,fluxpartB,sMag,sMagB,sR/1000,sq,se,sP/86400,log10(suQ),sfi);
fprintf(f,"#d:%1.2f   MB:%1.2f   qB:%1.2f   lmax:%1.2e\n",sd,sMB/6e24,sqB,slmax);


fclose(f);








/*
{
FILE *f;
f=fopen("ampter2.txt","w");

printf("\nIntegraling:\n");

double amp=0;
double Tegyen=0;

double temp1, temp2;

Tbb=0;
	
for(double T=100;T<10600;T=T*1.01){
for(double fi=0.5;fi<=180;fi=fi+0.1){
	
  //amp=2*integral(0, fi/57.3, 0)/h(sin(0.5*fi/57.3));
  
  Tbb=T;
  
  //temp1=integral(0, fi/57.3, 0);
  //temp2=integral(fi/57.3, 180/57.3, 0);
  
  Tegyen=pow( pow(Tbb,4)+pow(T,4) ,0.25);
  cc=h(sin(0.5*fi/57.3));
  Ts=pow(pow(Tbb,4)+pow(T,4)/cc,0.25);
  
  //amp=(temp1*FTSint(lmin, lmax, Ts) + temp2*FTSint(lmin, lmax, Tbb) ) / (0.5*FTSint(lmin, lmax, Tegyen));
  amp=(intcosin(0,fi/57.3)*FTSint(lmin, lmax, Ts) + intcosin(fi/57.3,180/57.3)*FTSint(lmin, lmax, Tbb) )/  (0.5*FTSint(lmin, lmax, Tegyen));
	
  fprintf(f,"%f  %f  %f\n",T,fi,amp);
  }
  fprintf(f,"\n");
  printf("T=%d\n",int(T));
  }
	
fclose(f);	
}
*/





//printf("\n\ndone\n");  
//getchar();
return 0;   
    
}


//2*pi*R integral 0 to R/2   cos(pi/2*x/R)*sin(pi/2*x/R) dx
//2*pi*R integral 0 to R sin(pi/2*x/R) dx


