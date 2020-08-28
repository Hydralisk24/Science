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




//#include <Dense>


using namespace std;

//using namespace Eigen;





/*
Teendo:
Gx/ref arany ne legyen nagy? (2nel nagyobb???)
telluric vonalakat is belerakni?
Extinkciora/vorosedre korrigalni a bemenetet? (alul a kodban vannak megoldasok ra)
a hiba: err ne allando legyen hamem valami fuggvenyt illeszteni arra is (de az fix) (pl ext/vorosedes korrigalo fuggvenyt) (ez viszont a kimenetet is valtoztatja, mert most le van osztva err-rel)?
c0 szamiasnak ne 100, 200, ... elem hanem lmax 1/5, 2/5, ... kent szamolni?

Likelihood-t atgndolni/ellenorozni (most khinegyzet/gyok(N) van ami jonak tunik, de tenyleg jo?)?
*/





// l : hullamhossz, wavelenght

// meres, data: spektrum (spec)
// csak egyszer van beolvasva
double specdl;  // lepeskoz, dl of data
double specl0;  // kezdo ertek, start l of data
double speclmax; // utolso ertek, end l of data

double extdl;  // lepeskoz, dl of extintion model
double extl0;  // kezdo ertek, start l of extintion model
double extlmax; // utolso ertek, end l of extintion model



// modelhez kello dolgok, used by the models
// minden loopban bevan olvasva!
double dt;  // lepes koz idoben, difference in time of the model
double t0;  // start t of model
double tmax; // end  t of model

double dl; // dl of model
double l0; // start l of model
double lmax; // end l of model

int g; // galaxy ID


vector <double> ext;  // extinction model




// galaxy modelhez kello dolgok, used by the galaxy models
// csak egyszer van beolvasva
vector <double> lgxmax;   // dl of GX
vector <double> lgx0;     // l0 of GX
vector <double> dlgx;     // lmax of GX



// parameter (a, b, c, z, gx ...) tarolok az alapok mellett, globalisan jol jon, grad = gradiantre ugyan ez
vector <double> parameters;
vector <double> grad;



// normalisation variables
int norman=0;
double norma=0;




				/* declaration */


					//egyenlore hardcoded, de majd belesz vive config fajlba !!!!
					// hiba erteke, error
						//double err=0.5e-16;  //PSN NGC6412
					double err=0.5e-14;
					int uniformerr=1;   // uniform error, or use a model?
					int correarth=0;   // apply the earth extionction model?
					double ebv=0;     // interstellar extiontion E(B-V)
					//fuggetlen voroseltolodas, red shift
						//double z0=0.0044; //PSN NGC6412  
					double z0=0.00013;   // min value, host galaxy z
					double zmax=0.001;   // maximum alloved
					int fixzg=0;  // fit z of galax too? (0 to yes)
					//idok, priori of time
					double priort0=1;   // minimum alloved epoch for the SN (0 is OK)
					double priortmax=90; // maximum alloved epoch for the SN (high values have no effects)
					// limits, parabola parameter limits
					double amax=6;
					double amin=-6;
					double bmax=6;
					double bmin=-3.5;

					
					//technical:
					//step number, decay numbers, step when to change to grad from mcmc
					double dstep=0.01;   // grad modszernel mekkora lepesekkel mozduljon odebb
					double gradtau=50.;  // grad modszernel a lepeskoz csokkenese, dstep exp(-3.5*q/gradtau)-val csokken, 
										 // q: akkor no ha nem csokken a khi negyzet, -3.5 egy atlag, random szam 2 es 7 kozott (hardcoded)
					double gradmax=150;  // ha q>gradmax akkor abba hagyja a gradiens modszert
					
					double Beta=0.1;     // mcmc modszernel a monte carlo parameter (linearisban es logaritmusban is)
					double mcmctau=5.;   // mcmc modszernel a MC parameter csokkenese, Beta exp(-qq/mcmctau)-val csokken, 
										 // qq: akkor no ha az elfogadasi arany kisebb mint 10% (100bol), hardcoded
					// hany MCMC lepes legyen? elso: GX nelkuli eset (legelso), masodik GX-val eset)
					double mcmclimit1=2000;
					double mcmclimit2=1000;
					double mcmcmax=3;    // ha qq>gradmax akkor abba hagyja a mcmc modszert
					
					double lx=6300; // desired wavelength for err autosearch
  				    int nl=300;  // err autosearch sample number
					



int samplederr=0; //seged valtozo






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
std::ifstream f("paramSP.inp"); 
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






 void defaultdata()			 /*Writeing default input file*/
 {
 
 FILE *f;

 f=fopen("paramSP.inp","wt");	
 

   fprintf(f,"%1.0e\t[The uniform error of the data, use - to fit]\n",-1e-15);   	  /*err*/
   fprintf(f,"%d\t[Use uniform error (1) or a model (0)?]\n",0);   	  /*uniformerr*/
   fprintf(f,"%d\t[Apply the extinction model (1)?]\n",0);   	  /*correarth*/
   fprintf(f,"%d\t[Interstellar E(B-V)]\n",0);   	      /*ebv*/
   fprintf(f,"%d\t[Red shift (z) of the host]\n",0);   	  /*z0*/
   fprintf(f,"%1.3f\t[Maximum allowed red shift (z)]\n",0.001);   	  /*zmax*/
   fprintf(f,"%d\t[z of galaxy fix (=z0) (1) or fited as well (0)?]\n",1);   	  /*fixzg*/
   fprintf(f,"%d\t[MIN SN epoch (day)]\n",0);   	  /*priort0*/
   fprintf(f,"%d\t[MAX SN epoch (day)]\n\n",90);   	  /*priortmax*/
   
   fprintf(f,"%d\t[MIN value for a (a*x*x)]\n",-6);   	  /*amin*/
   fprintf(f,"%d\t[MAX value for a (a*x*x)]\n",6);   	  /*amax*/
   fprintf(f,"%d\t[MIN value for b (b*x)]\n",-3);   	  /*bmin*/
   fprintf(f,"%d\t[MAX value for b (b*x)]\n\n",6);   	  /*bmax*/
   
   fprintf(f,"%1.2f\t[Grad method step]\n",0.04);   	  /*dstep*/
   fprintf(f,"%d\t[Grad MAX precision, exp(- _) model]\n",3);      /*gradmax*/
   fprintf(f,"%d\t[Grad increment number for exp model]\n",10);   	  /*gradtau*/
   fprintf(f,"%1.1f\t[Initial Metropolis-Hastings step parameter]\n",0.1);   	  /*Beta*/
   fprintf(f,"%d\t[Maximal M-H step precision, exp(- _) model]\n",2);      /*mcmcmax*/
   fprintf(f,"%1.1f\t[MCMC increment number for exp model]\n",1.5);   	  /*mcmctau*/
   fprintf(f,"%d\t[Number of max MCMC loops (first case)]\n",4000);   	  /*mcmclimit1*/
   fprintf(f,"%d\t[Number of max MCMC loops (other case)]\n",2000);   	  /*mcmclimit2*/
   fprintf(f,"%d\t[Desired wavelength for err autosearch]\n",6300);   	  /*lx*/
   fprintf(f,"%d\t[Sample number for err autosearch]\n\n",300);   	  /*nl*/
   
   
   
   

 fclose(f);
 	
 }
 
 
 
 
 
 
void data()			 /*Reading input file*/
 {
 //FILE *f;
 //f=fopen("parametersMC.inp","rt");

 if(f==NULL){printf("No parameter file! Generating a default! Exiting!\n"); defaultdata(); exit(1);}

   

   
   
    err=read();   // error value
    uniformerr=read(); // uniform error?
	correarth=read(); // apply the earth extionction model?
	ebv=read();   // interstellar extiontion E(B-V)
	z0=read();   // min value, host galaxy z
	zmax=read();   // maximum alloved
	fixzg=read();  // galaxy z being fitted? (0 to yes)
	priort0=read();   // minimum alloved epoch for the SN (0 is OK)
	priortmax=read(); // maximum alloved epoch for the SN (high values have no effects)
	
	read();
	amin=read();  // limits, parabola parameter limits
	amax=read();
	bmin=read();
	bmax=read();
	read();

					
	//technical:
	//step number, decay numbers, step when to change to grad from mcmc
	dstep=read();   // grad modszernel mekkora lepesekkel mozduljon odebb
	gradmax=read(); // ha q>gradmax*gradtau akkor abba hagyja a gradiens modszert
	gradtau=read();  // grad modszernel a lepeskoz csokkenese, dstep exp(-3.5*q/gradtau)-val csokken, 
	// q: akkor no ha nem csokken a kgi negyzet, -3.5 egy atlag, random szam 2 es 7 kozott (hardcoded)
	
	Beta=read();     // mcmc modszernel a monte carlo parameter (linearisban es logaritmusban is)
	mcmcmax=read();   // ha qq>mcmcmax*mcmctau akkor abba hagyja a mcmc modszert
	mcmctau=read();   // mcmc modszernel a MC parameter csokkenese, Beta exp(-qq/mcmctau)-val csokken, 
	// qq: akkor no ha az elfogadasi arany kisebb mint 10% (100bol), hardcoded
	// hany MCMC lepes legyen? elso: GX nelkuli eset (legelso), masodik GX-val eset)
	mcmclimit1=read();
	mcmclimit2=read();
	
	lx=read(); // desired wavelength for err autosearch
  	nl=read();  // err autosearch sample number
	
	
	
	
	if(err<=0) cout << "Note: err <= 0 given, err will be calculated." << endl;
	if(z0<0) cout << "WARNING! z0 < 0 !" << endl;
	if(zmax<z0){ cout << "ERROR! z0 > zmax ! Exiting!" << endl; exit(0); }
	if(zmax==z0){ cout << "WARNING! z0 = zmax !" << endl; }
	if(priortmax<priort0){ cout << "ERROR! t0 > tmax ! Exiting!" << endl; exit(0); }
	if(priortmax==priort0){ cout << "WARNING! t0 = tmax !" << endl; }
	
	if(amax<amin){ cout << "ERROR! amin > amax ! Exiting!" << endl; exit(0); }
	if(amax==amin){ cout << "WARNING! amin = amax !" << endl; }
	if(bmax<bmin){ cout << "ERROR! bmin > bmax ! Exiting!" << endl; exit(0); }
	if(bmax==bmin){ cout << "WARNING! bmin = bmax !" << endl; }
	
	if(amax<0){ cout << "ERROR! amax < 0 ! Exiting!" << endl; exit(0); }
	if(amin>0){ cout << "ERROR! amin > 0 ! Exiting!" << endl; exit(0); }
	if(bmax<0){ cout << "ERROR! bmax < 0 ! Exiting!" << endl; exit(0); }
	if(bmin>0){ cout << "ERROR! bmin > 0 ! Exiting!" << endl; exit(0); }
	
	if(dstep<=0){ cout << "ERROR! dstep <= 0 ! Exiting!" << endl; exit(0); }
	if(Beta<=0){ cout << "ERROR! Beta <= 0 ! Exiting!" << endl; exit(0); }
	if(gradtau<=0){ cout << "ERROR! gradtau <= 0 ! Exiting!" << endl; exit(0); }
	if(mcmctau<=0){ cout << "ERROR! mcmctau <= 0 ! Exiting!" << endl; exit(0); }
	
	if(gradmax<0 && ( mcmclimit1<0 || mcmclimit2<0 || mcmcmax<0 )  ){ cout << "ERROR! Both method turned off ! Exiting!" << endl; exit(0); }
	
	if(uniformerr<0) uniformerr=0;
	if(uniformerr>1) uniformerr=1;
	
	//default if wrong value
	if(lx<=0) lx=6300;
	if(nl<=0) nl=300;
   
   
 

 //fclose(f);
 f.close();
 }


	
	
	
	
	
	
	
	
	
	
	
	
	
	










// visszaadja azon tomb szamot, ami x erteknel van
// kell hozza az elso ertek (a) utolso ertek (b) es a lepes koz (dx)
// csak egyenletes dx lepes kozokre jo!!!!!!
int inverztar(double x, double dx, double a, double b){
	double y;
	
	if(dx==0) return 0;
	
	//tar[i][0]=a + i *(b-a)/double(nn);
	//y = (x - a) * double(nn)/(b-a);
	if(x>b) x=b;
	y = (x - a) / dx;
	if(y<0) y=0;
	
	
	return int(y);
    }
    
    
    
// tanulsag: csonkol: 2.99 => 2
/*
cout << "1.5: " << inverztar(1.5,dt,t0,tmax) << endl;
cout << "1.75: " << inverztar(1.75,dt,t0,tmax) << endl;
cout << "1.99: " << inverztar(1.99,dt,t0,tmax) << endl;
cout << "2: " << inverztar(2,dt,t0,tmax) << endl;
cout << "2.25: " << inverztar(2.25,dt,t0,tmax) << endl;
cout << "2.5: " << inverztar(2.5,dt,t0,tmax) << endl;
cout << "2.75: " << inverztar(2.75,dt,t0,tmax) << endl;
cout << "2.95: " << inverztar(2.95,dt,t0,tmax) << endl;
cout << "3: " << inverztar(3,dt,t0,tmax) << endl << endl;
*/
    
    
    
       




//melyik nagyobb, akkor az
double mn(double a, double b){
	return (a>b?a:b);
}

//melyik kisebb, akkor az
double mk(double a, double b){
	return (a<b?a:b);
}







// MISC
double h(double x){ if(x*0==0) return x*x;    return 0; }




// !!!!!!!!!!!!!!!!!
//double sgn(double x){ return x; }
//double sgn(double x){ if(x==0) return 0; if(x>0 && x<0.00001) return 0.00001; if(x<0 && x>-0.00001) return -0.00001; return x; }
//double sgn(double x){ if(x>1) return 1; if(x<-1) return -1; return x; }

//double sgn(double x){ if(x>0) return 1; if(x<0) return -1; if(x==0) return 0; }
double sgn(double x){ if(x==0) return 0; if(x>0 && x<1) return 1; if(x<0 && x>-1) return -1; return x; }
//double sgn(double x){ if(x==0) return 0; if(x>0 && x<0.1) return 0.1; if(x<0 && x>-0.1) return -0.1; return x; }



// kiszuri a nem szamat (filters nan)
double number(double x){ if(x*0==0) return x;    return 0; }



 // illesztesi fuggveny   
// a*l*l+b*l+c      
	// c*(a*l*l+b*l+1)   = c*a*x*x + a*b*x + c
	// c*(a*h(l-b)+1)    =  c*a*h(x-b)+c  = c*a*x*x-c*a*b*2x  +c*a*b*b +c
	//  osztani lmax-xal
	//  a+b>0
	//  min: x=-b/2a   -> -b*b/4a  +c (v 1)  >0  ha x e 0:1 
double fl(double x, double a, double b, double c){
	return ( c*(a*x*x+b*x+1) );
}



// ext[]  1-t ad vissza, ha nincs beolvasva semmi, vagy uniform error van


// chi square with gx
double chi(vector <double> ref, vector <double> spec, vector <double> gx){
//double chi(double ref[], double spec[], double gx[]){
	
	double khi;
	double l;
	
	double z=parameters[0];
	double a=parameters[1];
	double b=parameters[2];
	double c=parameters[3];
	double ag=parameters[4];
	double bg=parameters[5];
	double cg=parameters[6];
	double zg=parameters[7];
	
	double lalso;
	double lfelso;
	lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
	lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	
	khi=0; norma=0;  norman=0; 
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	// khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		//khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*ref[0][inverztar(t,dt,t0,tmax)][inverztar(l*(1.+z),dl,l0,lmax)]
		khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*     ref[inverztar(l*(1.+z),dl,l0,lmax)]
		+            ( fl(1-l/lfelso,ag,bg,cg) )*  gx[inverztar(l*(1.+zg),dlgx[g],lgx0[g],lgxmax[g])]
		-spec[inverztar(l,specdl,specl0,speclmax)]    )/h(err * ext[inverztar(l,extdl,extl0,extlmax)] );
		norma=norma+1./h(err * ext[inverztar(l,extdl,extl0,extlmax)] ); norman++;
	}
	// gx kaphat az intervalluman kivuli szamot. ekkor 0-t ad vissza!!!!
	
	//cout << h(err * ext[inverztar(l,extdl,extl0,extlmax)] ) << endl;
	return khi;
	
	
}




// chi square without gx
double chiNG(vector <double> ref, vector <double> spec){
//double chi(double ref[], double spec[], double gx[]){
	
	double khi;
	double l;
	
	double z=parameters[0];
	double a=parameters[1];
	double b=parameters[2];
	double c=parameters[3];
	//double ag=parameters[4];
	//double bg=parameters[5];
	//double cg=parameters[6];
	//double zg=parameters[7];
	
	double lalso;
	double lfelso;
	lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
	lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	
	khi=0; norma=0;  norman=0; 
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	// khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		//khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*ref[0][inverztar(t,dt,t0,tmax)][inverztar(l*(1.+z),dl,l0,lmax)]
		khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*     ref[inverztar(l*(1.+z),dl,l0,lmax)]
		//+            ( fl(1-l/lfelso,ag,bg,cg) )*  gx[inverztar(l*(1.+zg),dlgx[g],lgx0[g],lgxmax[g])]
		-spec[inverztar(l,specdl,specl0,speclmax)]    )/h(err * ext[inverztar(l,extdl,extl0,extlmax)] );
		norma=norma+1./h(err * ext[inverztar(l,extdl,extl0,extlmax)] ); norman++;
		
	}
	// gx kaphat az intervalluman kivuli szamot. ekkor 0-t ad vissza!!!!
	
	return khi;
	
	
}







// chi square without gx
double chiNSN(vector <double> spec, vector <double> gx){
//double chi(double ref[], double spec[], double gx[]){
	
	double khi;
	double l;
	
	//double z=parameters[0];
	//double a=parameters[1];
	//double b=parameters[2];
	//double c=parameters[3];
	double ag=parameters[4];
	double bg=parameters[5];
	double cg=parameters[6];
	double zg=parameters[7];
	
	double lalso;
	double lfelso;
	//lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
	//lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	lalso=specl0;
	lfelso=speclmax;
	
	khi=0; norma=0;  norman=0; 
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	// khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		//khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*ref[0][inverztar(t,dt,t0,tmax)][inverztar(l*(1.+z),dl,l0,lmax)]
		khi=khi+h(   //( fl(1-l/lfelso,a,b,c) )*     ref[inverztar(l*(1.+z),dl,l0,lmax)]
		            ( fl(1-l/lfelso,ag,bg,cg) )*  gx[inverztar(l*(1.+zg),dlgx[g],lgx0[g],lgxmax[g])]
		-spec[inverztar(l,specdl,specl0,speclmax)]    )/h(err * ext[inverztar(l,extdl,extl0,extlmax)] );
		norma=norma+1./h(err * ext[inverztar(l,extdl,extl0,extlmax)] ); norman++;
	}
	// gx kaphat az intervalluman kivuli szamot. ekkor 0-t ad vissza!!!!
	
	return khi;
	
	
}












// Galaxy/SN model flux ratio
double SNR(vector <double> ref, vector <double> gx){
//double chi(double ref[], double spec[], double gx[]){
	
	double khi;
	double l;
	int norman=0; 
	
	double z=parameters[0];
	double a=parameters[1];
	double b=parameters[2];
	double c=parameters[3];
	double ag=parameters[4];
	double bg=parameters[5];
	double cg=parameters[6];
	double zg=parameters[7];
	
	double lalso;
	double lfelso;
	lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
	lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	
	khi=0;
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	// khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		//khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*ref[0][inverztar(t,dt,t0,tmax)][inverztar(l*(1.+z),dl,l0,lmax)]
		khi=khi+ number(    //( fl(1-l/lfelso,a,b,c) )*     ref[inverztar(l*(1.+z),dl,l0,lmax)]
		             ( fl(1-l/lfelso,ag,bg,cg) )*  gx[inverztar(l*(1.+zg),dlgx[g],lgx0[g],lgxmax[g])]
		/			 ( fl(1-l/lfelso,a,b,c) ) /    ref[inverztar(l*(1.+z),dl,l0,lmax)]  );
		norman++;
		//cout << gx[inverztar(l*(1.+z0),dlgx[g],lgx0[g],lgxmax[g])] << endl; getchar();
	}
	// gx kaphat az intervalluman kivuli szamot. ekkor 0-t ad vissza!!!!
	
	return khi/norman;
	
	
}








	// RESET:
	//cout.unsetf(ios::fixed);
	//cout.unsetf(ios::showpoint);
	//cout << setprecision(-1);







/*                                                           MAIN               */



int main(int argc, char **argv){
	
	
	//cout << setprecision(6) << showpoint << 1380.50 << endl;
	//cout << setprecision(6) << showpoint << 46.8385 << endl;
	//cout << setprecision(6) << showpoint << 7.43354e-013 << endl;
	//cout << setprecision(6) << showpoint << 0.270551 << endl;
	//return 0;
	
	/*	
MatrixXd m(2,2);
m(0,0) = 3;
m(1,0) = 2.5;
m(0,1) = -1;
m(1,1) = m(1,0) + m(0,1);
cout << m << std::endl << endl;
   
MatrixXd M(2,3);
M << 1,2,3,1,2,3;
MatrixXd MT=M.transpose();
cout << M << endl << endl;
cout << MT << endl<< endl;
cout << M*MT << endl<< endl;
cout << m.inverse()*m << endl<< endl;
cout << m.determinant() << endl << endl;

Vector2d v(1,2);
cout << v << endl << endl;
cout << m*v << endl << endl;


MatrixXd m2(3,3);
Vector3d v2(7,3,6);
m2 << 2,3,2,
      1,1,1,
      2,2,3;
cout << m2.inverse()*v2 << endl << endl;
//megolddas: 2 1 0
   

MatrixXcd mc(1,1);
mc << 1;
cout << mc << endl;

   
   exit(1);
	
*/	
	
	
	
	
	
	
	srand(time(NULL));
	data();

	
	
	
// global variable initialisatin	
parameters.push_back(0); // for z
parameters.push_back(0); // for a
parameters.push_back(0); // for b
parameters.push_back(1); // for c
parameters.push_back(0); // for ag
parameters.push_back(0); // for bg
parameters.push_back(0); // for cg
parameters.push_back(z0); // for zg


//gradient
grad.push_back(0); // for z
grad.push_back(0); // for a
grad.push_back(0); // for b
grad.push_back(0); // for c
grad.push_back(0); // for ag
grad.push_back(0); // for bg
grad.push_back(0); // for cg
grad.push_back(0); // for zg
	
	
	
	
	
	
    

// illesztesi valtozok, lasd fl() fuggveny fennt    
double khi=0; // chi square
double a,b,c; // for polinom
double ag,bg,cg;  // for GX polinom
double z;  // red shift
double c0;  // initial guess for c
double c0g;
double zg=z0;
					
					

// MCMC cuccai, azon ertekek ahonnan tovabb lepunk, accepted values for MCC
double tkhi=1e300;
double ta,tb,tc;
double tz;  
double tt;
double tag,tbg,tcg;
double tzg;


int tj;   // save for j


double t; // ido, time
double l; // hullamhossz, wavelenght




// legjobb ertekeket tarolja!!!!    saves the best values
vector <double> khiS;
vector <double> khiNS;
vector <double> aS;
vector <double> bS;
vector <double> cS;
vector <double> zS;
vector <double> tS;  
vector <double> targx; 

vector <double> agS;
vector <double> bgS;
vector <double> cgS;
vector <double> zgS;

vector <double> tSNR;

    
    
    
    
string sor;  // seged string beolvasashoz
double mezo1,mezo2,mezo3;  // beolvashoz segedek

int q=0;


// beolvasashoz seged tarolok
double tsave=-1000;
double lsave=0;



// list !!!!
// meres/data es lista vektorok
vector <double> spec;
vector <string> list;
vector <string> list2;
//vector <double> ext;   // global






			// meres es lista fileok
			//int argc, char **argv
			//std::ofstream fkispec;
			//fkispec.open("specki.dat",std::ofstream::out);
			//std::ifstream specIN("has.dat");  // data file
			string dict="Nugent_template/";
			string dict2="GX/";
			std::ifstream specIN;  // data file
			std::ifstream listIN("list");   // for SN model
			std::ifstream list2IN("list2"); // for GX model
			std::ifstream extIN("Ext.dat");  // data file
			
			if(argc>1){ specIN.open(argv[1],std::ifstream::in); cout << "Input: " << argv[1] << endl << endl; }
			else{        specIN.open("has.dat",std::ifstream::in); cout << "Input: " << "has.dat" << endl << endl; }
			
			//getchar();







// exit if these are not existing
if(!specIN.is_open()){ cout << "No data file! Exiting!\n"; exit(0); }
if(!listIN.is_open()){ cout << "No SN list file! (list) Exiting!\n"; exit(0); }
if(!list2IN.is_open()){ cout << "No Galaxy list file! (list2) Exiting!\n"; exit(0); }

if(!extIN.is_open() && uniformerr==0 ){ cout << "WARNING! No Extinction model file! (Ext.dat) Uniform used. \n"; extdl=0; extl0=0; extlmax=0; uniformerr=-1; }



// kimenet, output
// append!
std::ofstream fki;
//fki.open("ki.out",std::ofstream::app);

if(argc>1){ sor=argv[1]; sor=sor+".out"; fki.open((sor.c_str()),std::ofstream::app); }
else{        fki.open("ki.out",std::ofstream::app); }













// beolvasasok, read in


q=0;
//   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  DATA, meres beolvasas !!!!!!!!


while(getline(specIN,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    stringstream strs(sor);

    //strs >> mezo1;
    strs >> mezo2;
    strs >> mezo3;

  //Feldolgozzuk a mezoket.
  //cout  << mezo1 << "\t" << mezo2 << "\t" << mezo3 << endl;
  //cout  << mezo2 << "\t" << mezo3 << endl;
  
  if(q==0){ lsave=specl0=mezo2;}   // elso hullamhossz
  if(q==1){ specdl=mezo2-lsave; }  // masodik hullamhossz, dl kiszamitas a kulonbsegbol
  
  //spec.push_back( mezo3 );
  //double temp=(mezo3+0.5*err*RandomDoubleGauss(1, 0) )*  fl(1-mezo2/25000.,.1,.3,1);
  double temp=mezo3;
  //if(temp<0.25*err) temp=0.25*err;
  if(temp<0) temp=0;  // negativ ertek nem lehet
  spec.push_back(  temp   );  // vegso tarolas
  
  
  
  q++;


}

speclmax=mezo2;  // utolso hullamhossz


specIN.close();








//                                                    lista beolvasasok !!!!!!!!!!!!!
while(getline(listIN,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    //stringstream strs(sor);

    //strs >> mezo1;
    //strs >> mezo2;
    //strs >> mezo3;

  //Feldolgozzuk a mezoket.
  //cout  << mezo1 << "\t" << mezo2 << "\t" << mezo3 << endl;
  
  //list.push_back( "Nugent_template/" + sor );
  list.push_back( sor );
  

}


listIN.close();


cout << "read in files:" << endl;
for(int k=0;k<list.size();k++){
	
	cout << (dict+list[k]).c_str() << endl;
	
	}
	
	




	
	
	

// lista 2 beolvasas
while(getline(list2IN,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    //stringstream strs(sor);

    //strs >> mezo1;
    //strs >> mezo2;
    //strs >> mezo3;

  //Feldolgozzuk a mezoket.
  //cout  << mezo1 << "\t" << mezo2 << "\t" << mezo3 << endl;
  
  //list2.push_back( "GX/" + sor + ".sp" ); 
  list2.push_back( sor + ".sp" );
  

}


list2IN.close();
	
	
cout << endl;	
cout << "read in files 2:" << endl;
for(int k=0;k<list2.size();k++){
	
	cout << (dict2+list2[k]).c_str() << endl;
	
	}
	
	
	
	
	
	
	
	
if(extIN.is_open() && uniformerr==0 ){

q=0;	
while(getline(extIN,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    stringstream strs(sor);

    //strs >> mezo1;
    strs >> mezo2;
    strs >> mezo3;

  //Feldolgozzuk a mezoket.
  //cout  << mezo1 << "\t" << mezo2 << "\t" << mezo3 << endl;
  
  if(q==0){ lsave=extl0=mezo2;}   // elso hullamhossz
  if(q==1){ extdl=mezo2-lsave; }  // masodik hullamhossz, dl kiszamitas a kulonbsegbol
  
  
  double temp=mezo3;
  if(temp<=0) temp=1;  // negativ ertek nem lehet
  ext.push_back(  temp   );  // vegso tarolas
  
  
  
  q++;


}

extlmax=mezo2;  // utolso hullamhossz


extIN.close();
if(ext.size()<1){ cout << "ERROR! Extinction model wrong readin! Exiting!" << endl; exit(0); }


 }
 else{ ext.push_back( 1 ); extdl=0; extl0=0; extlmax=0; }
 // ha nincs beolvasas, ez az egyetlen 1es lesz a tombben es ez lest mindig meghivva
 
 
 
 
 


// Applies Earth Extinction correction
if(correarth==1 && ext.size()>1){
	
	
	double lalso;
	double lfelso;
	lalso=specl0;
	lfelso=speclmax;
	
	
	for( l=lalso  ;l<=lfelso  ; l=l+specdl ){
		
		 spec[inverztar(l,specdl,specl0,speclmax)]=spec[inverztar(l,specdl,specl0,speclmax)]*ext[inverztar(l,extdl,extl0,extlmax)];
		
	}
	
	
	
}





//computes reddened fluxes from Fitzpatrick & Massa astro-ph 0705.0154
// https://arxiv.org/pdf/0705.0154.pdf
// Applies Interstellar Extinction correction
if(ebv!=0){
	
	
	double RV=3.1; double o1=2.055; double o2=1.322; double o3=0.00; double kir=1.057;
	double x0=4.592; double g=0.922; double c1=-0.175; double c2=0.807; double c3=2.991; double c4=0.319; double c5=6.097;
	double DD=0; double A=0; double x=1; double k;
	double corr=1;
	
	
	double lalso;
	double lfelso;
	lalso=specl0;
	lfelso=speclmax;
	
	
	for( l=lalso  ;l<=lfelso  ; l=l+specdl ){
		
		x=1e4/l;
		if(l<3300) { DD = x*x/( h(x*x-x0*x0) + h(g*x));
                  if(x<=c5){ k = c1+c2*x+c3*DD; } 
				  else {k = c1+c2*x+c3*DD+c4*h(x-c5); }
		          }
    	if(l>=3300 && l<4000) { k = o1+(o2-o1)/(4000.0-3300.0)*(l-3300); }
    	if(l>=4000 && l<=5500) { k = o2+(o3-o2)/(5430.0-4000.0)*(l-4000); }
    	if(l>5500) { k = kir*pow(1e4/l,1.84) - RV; } 
    	A = ebv * (RV + k);
    	corr=pow(10,-0.4*A);
    	// m0=-2.5*log(f)/log(10)+A; f0=10^(-0.4*m0);
    	// f0=f*10^(-0.4*A);
		
		
		spec[inverztar(l,specdl,specl0,speclmax)]=spec[inverztar(l,specdl,specl0,speclmax)]/corr;
		
	}	
	
	
	
	lalso=extl0;
	lfelso=extlmax;
	// also done to the ext[] as it is used as sigma!
	//not applies twice, beacause this is AFTER it
	if(ext.size()>1) for( l=lalso  ;l<=lfelso  ; l=l+extdl ){
		
		x=1e4/l;
		if(l<3300) { DD = x*x/( h(x*x-x0*x0) + h(g*x));
                  if(x<=c5){ k = c1+c2*x+c3*DD; } 
				  else {k = c1+c2*x+c3*DD+c4*h(x-c5); }
		          }
    	if(l>=3300 && l<4000) { k = o1+(o2-o1)/(4000.0-3300.0)*(l-3300); }
    	if(l>=4000 && l<=5500) { k = o2+(o3-o2)/(5430.0-4000.0)*(l-4000); }
    	if(l>5500) { k = kir*pow(1e4/l,1.84) - RV; } 
    	A = ebv * (RV + k);
    	corr=pow(10,-0.4*A);
		
		
		ext[inverztar(l,extdl,extl0,extlmax)]=ext[inverztar(l,extdl,extl0,extlmax)]/corr;
		
	}
	
	
}

 
 
 
 
 
 
   // auto search for err
   //spec[inverztar(l,specdl,specl0,speclmax)]
  if(err<=0){
  	//double lx=6300;
  	//int nl=300;
  	samplederr=1;


  	
  	if(lx>speclmax) lx=speclmax-300;
  	if(lx-nl*specdl<specl0) lx=specl0+300+nl*specdl;
  	
  	err=0;
  	int j=inverztar(lx,specdl,specl0,speclmax);
  	
  	if(nl>j) if(nl>spec.size()-1){ nl=spec.size()-1; nl=j; }
  			else{ nl=j; }
  	
  	//failsafe
  	if( spec.size()<j ){ cout << "ERROR in the auto search for err! (big) Exiting!" << endl; exit(1); }
  	if(  0>j-nl        ){ cout << "ERROR in the auto search for err! (low) Exiting!" << endl; exit(1); }
  	//cout << "OK" << endl; return 0;
  	
  	double atlag=0;
  	for(int i=0;i<nl;i++) atlag=atlag+   spec[j-i];
  	atlag=atlag/nl;
  	
  	for(int i=0;i<nl;i++) err=err+h(spec[j-i]-atlag);
  	err=sqrt(err/nl); 	
  	
  	//err=0.5e-14;
  	
  }

	

//for(int i=0;i<ext.size();i++)	cout << ext[i] << endl;
//cout << (ext[inverztar(8000,extdl,extl0,extlmax)] ) << endl;
//cout << h(err * ext[inverztar(4000,extdl,extl0,extlmax)] ) << endl;
//return 0;
	
cout << "\nerror value: " << err << endl;
cout << "\nContinuing.\n" << endl;
//cout << "\nPress enter to continue.\n" << endl;
//getchar();













//                                         Galaxys beolvasas         !!!!!!!!

vector<vector <double> > gx[2];   // galaxy hasznalt taroloja, galaxy GX used variable
//vector <double> lgxmax;

for(int k=0;k<list2.size();k++){
	
	std::ifstream IN( (dict2+list2[k]).c_str());
	
	if(!IN.is_open()){ cout << "Error in galaxy list file! Invalid file name:\n"<< ( (dict2+list2[k]).c_str() ) <<"\nPlease check list2 file! Exiting!\n"; exit(0); }
	
	double mezo1,mezo2,mezo3;
	
	
	gx[0].push_back( vector <double> () );
	gx[1].push_back( vector <double> () );
	
	
	
q=0;	
while(getline(IN,sor))
     {

      if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec

      stringstream strs(sor);

      //strs >> mezo1;
      strs >> mezo2;
      strs >> mezo3;
    

	  gx[0][k].push_back(mezo3);   // saveing the flux
	  gx[1][k].push_back(mezo2);   // saveing the wavelength too just in case

	  q++;
	
      }

	
	// elso es utolso erteket 0-nak adjuk, ugyhogy tartjuk a rendet a hullamhosszban!
	// OK: igy ha inverz tar visszahivaskor olyan erteket adunk ami nincs benne, akkor 0 lesz az erteke
	// egy plusz sort hozza adunk +dl hullamhosszal, de 0 ertekkel
	// FIGYELEM!!!!    gx[1][k][11] -nek leteznie kell!!!
	
	if(gx[1][k].size()<12) cout << (dict2+list2[k]).c_str() << "  WARNING!  GX model lenght is smaller than 11, this may cause a crash! Checking is recommended! " << endl;
	
	dlgx.push_back(gx[1][k][11]-gx[1][k][10]);
	
	mezo2=mezo2+gx[1][k][11]-gx[1][k][10];
	gx[0][k].push_back(0);
	gx[1][k].push_back(mezo2);
	lgxmax.push_back(mezo2);
	
	// egy plusz sort elore tolunk -dl hullamhosszal, de 0 ertekkel
	gx[0][k].insert(gx[0][k].begin(),0);
	gx[1][k].insert(gx[1][k].begin(),   2*gx[1][k][0]-gx[1][k][1]   );
	
	lgx0.push_back(gx[1][k][0]);

	
	}
		
	
	//int g=0; cout <<  gx[0][g][inverztar(0,gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] << endl; exit(1);
	























/*                       REFERENCE READ IN AND THE MINIMALISATION METHOD !!!!!!!!!!!!!!!!!!                           */

    

//list.push_back("Nugent_template/hyper_flux.v1.2.dat");
//list.push_back("Nugent_template/sn1a_flux.v1.2.dat");

// search for SN model list, done one by one (grid method)
// +1 part is for the galaxy only fit !!!!!
// k<list.size()  when the are actual SN fit
for(int k=0;k<list.size()+1;k++){
	


  int x=0;
  
  std::ifstream IN( (k<list.size()?(dict+list[k]).c_str():"") );

  //   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
  vector<vector <double> > ref[1];  // SN model array
  vector <double> reft;    // SN model to store time

  double mezo1,mezo2,mezo3; // seged beolvaskor
  int i=-1;  // beolvasasi seged index, -1 rol indul mert elobb leptetjuk


  vector <double> parasave;  // adott idoben a legjobb fit.

  parasave.push_back(0); // z
  parasave.push_back(0); // a
  parasave.push_back(0); // b
  parasave.push_back(1); // c
  parasave.push_back(0); // ag
  parasave.push_back(0); // bg
  parasave.push_back(0); // cg
  parasave.push_back(1e300); // khi
  parasave.push_back(1e300); // sqrt(khi/norman)
  parasave.push_back(z0); // zg







  tsave=-1000;  // -1000 hogy mindig nagyobb legyen
  lsave=0;
  q=0;
  t0=0; tmax=0; dt=0;


  if(k<list.size()){
  
    // beolvassa a filet! Mindig egyet olvass csak be, es azzal dolgozik, majd torli!!!!!!!!!!!!!!!!
    cout << endl << list[k].c_str() << endl;

    while(getline(IN,sor))
       {

      if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec

      stringstream strs(sor);

      strs >> mezo1;
      strs >> mezo2;
      strs >> mezo3;
    
      if(mezo1>=(priort0>0?priort0:1) && mezo1<=(priortmax>priort0?priortmax:1000)){

    	//Feldolgozzuk a mezoket.
      	//cout  << mezo1 << "\t" << mezo2 << "\t" << mezo3 << endl;
      	if(q==0){ t0=mezo1; lsave=l0=mezo2;}
      	if(q==1){ dl=mezo2-lsave; }
  
  
  	  	if(mezo1!=tsave){
	     ref[0].push_back( vector <double> () );
	     //ref[1].push_back( vector <double> () );
	     reft.push_back(mezo1);
	     i++;
	     dt=mezo1-tsave;
	     tsave=mezo1;
	     }

	  	ref[0][i].push_back(mezo3);
	  	//ref[1][i].push_back(mezo1);
	  	

	  	q++;
	  	tmax=mezo1;
	  	lmax=mezo2;
	
      	}
      
     }


    //   UNUSED
    // kiserlet arra, hogy otvozzunk ketto meglevo adatot, de linearis kombinacio kellene inkabb, es ez megneheziti, nem biztos hogy fizikailag megalapozott
    /*
    if(dt!=1){
    ref[0].push_back( vector <double> () );
    ref[1].push_back( vector <double> () );

    for(int j=0;j<q;j++){
    ref[0][2].push_back( 0.5*(800*ref[0][0][j]+1*ref[0][1][j]) );
    ref[1][2].push_back( int(0.5*(ref[1][0][j]+ref[1][1][j])+0.5) );
    }
    }*/
  
  
  




    // readin infos
    cout << "t0:  " << t0 << "  tmax:  " << tmax  << "  last dt:  "  << dt  << endl;
    cout << "t samples:  " << ref[0].size() << endl;
    //cout << "mas tmax:   " << ref[1][ref[0].size()-1][1] << endl;
    cout <<"read in done" << endl << endl;

    fki << endl <<list[k].c_str() << endl;
    fki << "t0:  " << t0 << "  tmax:  " << tmax  << "  last dt:  "  << dt  << endl;
    fki << "t samples:  " << ref[0].size() << endl;
    //exit(1);
  
  
  }
  else{
  	cout << "Galaxy only" << endl;
  	fki  << "Galaxy only" << endl;
  }










  // kezdo ertekek!!
  tt=0.5*(t0+tmax);
  a=ta=1e-6;
  b=tb=1e-6;
  z=tz=z0+1e-6;
  tkhi=1e300;
  zg=tzg=z0;

  ag=tag=1e-6;
  bg=tbg=1e-6;

  c=tc=c0=1;
  cg=tcg=c0g=1;


  // best fit save initialization
  // every SN model storen in a new array
  // szukseges itt initializalni (pushbackelni) oket, kulonben elcsuszik a dolog
  khiS.push_back(tkhi);
  khiNS.push_back(tkhi);
  aS.push_back(ta);
  bS.push_back(tb);
  cS.push_back(tc);
  zS.push_back(tz);
  tS.push_back(tt); 
  targx.push_back(-1);  // galaxy ID
  agS.push_back(tag);
  bgS.push_back(tbg);
  cgS.push_back(tcg);
  zgS.push_back(tzg);

  tSNR.push_back(0);  // SN/galaxy model flux ratio
  
  
  
  
  
  parameters[0]=tz;
  parameters[1]=ta;
  parameters[2]=tb;
  parameters[3]=tc;
  parameters[4]=tag;
  parameters[5]=tbg;
  parameters[6]=tcg;
  parameters[7]=tzg;

  grad[0]=0;
  grad[1]=0;
  grad[2]=0;
  grad[3]=0;
  grad[4]=0;
  grad[5]=0;
  grad[6]=0;
  grad[7]=0;


 


  // errors
  if(!IN.is_open() && k<list.size()){ 
	  cout << "ERROR! Invalid file name. Please check list file! Skipping..." << endl << endl; 
	  fki  << "ERROR! Invalid file name. Please check list file! Skipping..." << endl;
	  continue; 
	  }

  
  if(ref[0].size()<=0 && k<list.size()){
	  cout << "ERROR! No model data in this given time period. Skipping..." << endl << endl;   
	  fki  << "ERROR! No model data in this given time period. Skipping..." << endl;
	  continue;
	  }
	  
  
  if(k<list.size()) if(ref[0][0].size()<501) cout << list[k].c_str() << "    SN model lenght is smaller than 500, this may cause a crash! Checking is recommended! " << endl;
  if(spec.size()<501) cout << "    Data lenght is smaller than 500, this may cause a crash! Checking is recommended! " << endl;



  //std::ofstream ftemp;
  //ftemp.open(k==0?"temp.dat":"",std::ofstream::out);








  /*
  Ha dt>1 akkor minden lehetoseg vegig probal mert nincs sok
  Ha dt=1 akkor 5-osevel lekphed vegig, igy gyorsabb
  de a legjobb ertek korul egy +-9 ertekben azert megegyszer vegig nezi 1es lepes kozzel
  (-9 fix, de mindig a legjobb+9 a max)
  */

  /*
  Minimalizacio:
  https://neutrium.net/mathematics/least-squares-fitting-of-a-polynomial/

  Megjegyzes:
  ket valtozo T plank fuggveny aranya kozelitoleg exp fuggveny
  exp((h*c/l/k)*(dT/(T1*T2)))  = exp(C*dT/l)
  T=5000 korul olyan 2 mikronnal kezd serulni a kozelites
  300nm-nel kb 5, 2000nm-nel kb 1.5, ha dT=1000 (5000 es 6000)
  = [1:6] exp(1.6*(dT/1000)/x)  kb ugyan ezt adja
  */





  // search for time in grid method!
  
  //if dt=1 then do it twice
  for(int o=0;o<=(dt==1?1:0);o++)
   // j stand for the index of the time, tj saves the best fit j
   for(int j=(o==0?0:mn(tj-9,0));     j<=(o==0?ref[0].size()-1+(k>=list.size()):mk(tj+9,ref[0].size()-1));  ){



 	
	 if(k<list.size()){
	
	  //t=ref[1][j][1];  // stores the time (not the index), get it from the array
	  t=reft[j];  // stores the time (not the index), get it from the array
 
 	 // initial guess for the c
	 c=tc=c0=0.2*(
 	 spec[100]/ref[0][j][100]+
 	 spec[200]/ref[0][j][200]+
 	 spec[300]/ref[0][j][300]+
 	 spec[400]/ref[0][j][400]+
 	 spec[500]/ref[0][j][500]
 	 );
 	 cout << "t= " << int(t) << "  start c:   " << tc << endl ;//<< endl; //getchar();
 	 fki << "t= " << int(t) << "  start c:   " << tc << endl ;
    }
    else{
    	t=0; c=tc=c0=1;
	}


 	 cg=tcg=c0g=1;

 
 
 
     // kezdeti ertek adas looponkent is
	 parasave[0]=z0+1e-6;
	 parasave[1]=1e-6;
	 parasave[2]=1e-6;
	 parasave[3]=c0;
	 parasave[4]=1e-6;
	 parasave[5]=1e-6;
	 parasave[6]=c0g;
	 parasave[7]=1e300;
	 parasave[8]=1e300;
	 parasave[9]=z0;
 
	 tkhi=1e300;



	 /* itt kell a int(list2.size() )  kulonben az allitas nem igaz!!  ha az osszehasonlito negativ akkor nem jo, pozitivra jo*/
	 // search for galaxy model (grid method), note: g=-1 tands for no galaxy model used
	 for(g=-1;g<int(list2.size() );g++){
	 //for(g=-1;g<0;g++){
	 
	 
	 	 if(k>=list.size() && g<0) continue;  // skipping the no SN and no GX part
	 
 		 //cout << "g: " << g << endl;
 
 		 if(g>=0) cg=tcg=c0g=0.2*0.1*(
 		 spec[100]/gx[0][g][100]+
 		 spec[200]/gx[0][g][200]+
 		 spec[300]/gx[0][g][300]+
 		 spec[400]/gx[0][g][400]+
 		 spec[500]/gx[0][g][500]
 		 );

 		 //cout <<"Galaxy ID: "<<g<< "  start c_g:   " << tcg << endl; //getchar();
 		 //fki <<"Galaxy ID: "<<g<< "  start c_g:   " << tcg << endl;
 
 
 		 double q=0;
 		 int qq=0;
 		 int accept=0; int all=0;
 		 
 		 int lim1=mcmclimit1;
 		 int lim2=mcmclimit2;



		 /* The minimalisation part! 
		  Search for the best khi 
		  k,t,g are searched in grid method, but the other parameters (z,a,b,c,ag,bg,cg) is searched here
		  two (!) methods used: first MCMC, then deepest step
		 */
		 // numbers here for i are max values after the search stops, a guarantee not to become an infinite loop
		 for(int i=0;i<(g<0?20000:10000);i++){
	

	
	
			 // In the second iteration (second g) and later on, the first values are the previous iterations best fits
			 if(g>=0 && i==0){
				if(k<list.size()){ ta=parasave[1]; tb=parasave[2]; tc=parasave[3]; tz=parasave[0]; }         //tkhi=parasave[7];
				if(g>0){  tag=parasave[4]; tbg=parasave[5]; tcg=parasave[6]; if(fixzg==0) tzg=parasave[9]; }
				parameters[0]=tz;
	 			parameters[1]=ta;
	 			parameters[2]=tb;
	 			parameters[3]=tc;
	 			parameters[4]=tag;
	 			parameters[5]=tbg;
	 			parameters[6]=tcg;
	 			if(fixzg==0) parameters[7]=tzg;
	 			if(k<list.size()) tkhi=chi(ref[0][j],spec,gx[0][g]);
	 			else tkhi=chiNSN(spec,gx[0][g]);
	 			
	 			if(fixzg==0 && g>=0) if(tzg>tz){ tzg=tz; /*cout<<"tzg>tz "<<tzg<<" "<<tz<<endl; getchar();*/ }
	    		}
	    		
	    	if(k>=list.size() && g==0 && i==0){ tcg=c0g; parameters[6]=tcg; tkhi=1e300; }  // galaxy only elso esetben kell egy kezdeti ertek
	
	
			
	
	
			 //  MCMC
			 //  for the first time it will do it until mcmclimit1, else mcmclimit2
			 if(i< (g<0?lim1:lim2) ){

				 if(k<list.size()){
				 	// osszelettek kotve, hogy kudarc eseten mindkettot ujra generalja, ez lehet lassabb de a szabalyzask miatt nem okoz vegtelen loopot
				 	x=0; do{ b=tb+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   
	    		 	a=ta+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   
					x++; if(x>1000) fki << "LOOP: a b ta tb: " << a <<" " << b <<" " << ta <<" " << tb << endl; 
					}    			  	
				 	while((b>bmax || b<bmin || a>amax || a<amin || a+b<-1 || (b<-2 && a<0.25*b*b) ) && x<1010 );
	 			
	 			 	// linear
				 	//do{ b=tb+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    			  	while(b>bmax || b<bmin);
				 	//do{ a=ta+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    			  	while(a>amax || a<amin || a+b<-1 || (b<-2 && a<0.25*b*b) );
				 	//do{ c=tc+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    				while(c<=0 );
				 	x=0; do{              z=tz+Beta*(zmax-z0)*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);  
					if(fixzg==0 && g>=0) zg=tzg+Beta*(zmax-z0)*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);
					x++; if(x>1000) fki << "LOOP: z tz zg tzg: " << z <<" " << tz <<" " << zg <<" " << tzg << endl; 
					}    			  	while((z<z0 || z>zmax || ((fixzg==0 && g>=0) && (zg<z0 || zg>z)) ) && x<1010 );

				 	// log
				 	//do{ a=ta*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));   }    			  			while(a<1e-6 || a>4);
				 	//do{ b=tb*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));   }    			  			while(b<1e-6 || b>4);
				 	x=0; do{ c=tc*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));  
					   x++; if(x>1000) fki << "LOOP: c tc: " << c <<" " << tc << endl; }		    			  	while(0); //while(c<c0/10000 || c>c0*10000 );
				 	//do{ z=tz*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));   }    			 		 	while(z<z0+1e-6 || z>0.2 );
				 	}
	 
				 if(g>=0){
	   
	   				x=0; do{ bg=tbg+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   
	       			ag=tag+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   
					x++; if(x>1000) fki << "LOOP: ag bg tag tbg: " << ag <<" " << bg <<" " << tag <<" " << tbg << endl;
					}    			  					
	       			while((bg>bmax || bg<bmin || ag>amax || ag<amin || ag+bg<-1 || (bg<-2 && ag<0.25*bg*bg) ) && x<1010 );
	   
	   				//do{ bg=tbg+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    			  					while(bg>bmax || bg<bmin);
	   				//do{ ag=tag+Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    			  					while(ag>amax || ag<amin || ag+bg<-1 || (bg<-2 && ag<0.25*bg*bg) );
	 	
	   				//do{ ag=tag*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));   }    			  		while(ag<1e-6 || ag>4);
	   				//do{ bg=tbg*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));   }    			  		while(bg<1e-6 || bg>4);
	   				x=0; do{ cg=tcg*pow(10.,Beta*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0));  
					   x++; if(x>1000) fki << "LOOP: cg tcg: " << cg <<" " << tcg << endl; }		    			  	while(0); //while(cg<c0g/10000 || cg>c0g*10000 );
	   				//if(fixzg==0) do{ zg=tzg+Beta*(zmax-z0)*exp(-qq/mcmctau)*RandomDoubleGauss(1.,0);   }    			  	while(zg<z0 || zg>z );   // vegtelen loopot dobhat???????
	   				//cout << zg << endl;
	   				
	   				}

	
	
	
				 // khi calculation, unused but use the same method as the chi() function
				 // gx kaphat az intervalluman kivuli szamot. ekkor 0-t ad vissza!!!!
				 /*
				 double lalso;
				 double lfelso;
				 lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
				 lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	
				 khi=0; norma=0;  norman=0; 
				 //   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
				 // khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
				 for(l=lalso  ;l<=lfelso  ; l=l+specdl){
					//khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*ref[0][inverztar(t,dt,t0,tmax)][inverztar(l*(1.+z),dl,l0,lmax)]
					khi=khi+h(   ( fl(1-l/lfelso,a,b,c) )*     ref[0][j][inverztar(l*(1.+z),dl,l0,lmax)]
					+            ( fl(1-l/lfelso,ag,bg,cg) )*  gx[0][g][inverztar(l*(1.+z0),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])]
					-spec[inverztar(l,specdl,specl0,speclmax)]    )/h(err);
					norma=norma+1./h(err); norman++;
					}*/
				 
	 
				 // chi, khi calculation (parameters a global parameters)
				 parameters[0]=z;
				 parameters[1]=a;
				 parameters[2]=b;
				 parameters[3]=c;
				 parameters[4]=ag;
				 parameters[5]=bg;
				 parameters[6]=cg;
				 if(fixzg==0) parameters[7]=tzg;
				 if(k<list.size()){
				 	if(g<0)  khi=chiNG(ref[0][j],spec); 
				 	if(g>=0) khi=chi(ref[0][j],spec,gx[0][g]);  
				 }
				 else khi=chiNSN(spec,gx[0][g]);
				 // azonos a kikommentteltel
	 
	
	
				 // metropolis ratio
				 if( exp(-0.5* khi   +0.5* tkhi ) >= RandomDouble(0.,1.) ){
				 //if( tkhi/khi >= RandomDouble(0.,1.) ){
	 				if(k<list.size()){tt=t; ta=a; tb=b; tc=c; tz=z;}   tkhi=khi;    accept++;
		 			if(g>=0){tag=ag; tbg=bg; tcg=cg;   if(fixzg==0) tzg=zg; }    //cout << "accept\n";
	    			}
	    
				 all++; // counting MCMC step
	
				 if(i%100==0){
					//cout << "Ratio= " << 100.*accept/all <<"    qq= "<<qq<< endl;
					if(1.*accept/all<0.1) qq++;
					if(1.*accept/all>0.5) qq--;
					if(qq<0) qq=0;
					if(qq>mcmctau*mcmcmax) qq=mcmctau*mcmcmax;
					if(1.*accept/all<0.04 && qq==mcmctau*mcmcmax){ //cout << "break at: " << i << endl; 
															(g<0?lim1=i:lim2=i);
															//break; 
															}									
					accept=0; all=0;
	    			}
	    

			 } // MCMC end
	
	
			 // seting the parameters for the best values or the last accepted values, used for first values for the grad method
			 // azert last accepted, mert GX loopra nem orizzuk a legjobb erteket, igy nem tudjuk mi a legjobb ertek aktualisan
			 if(i== (g<0?lim1:lim2) ){
				//a=parasave[1]; 
				//b=parasave[2]; 
				//c=parasave[3]; 
				//z=parasave[0];    
				//if(g>0){  ag=parasave[4]; bg=parasave[5]; cg=parasave[6];   }  
				//tkhi=parasave[7];
				
				parameters[0]=tz;
	 			parameters[1]=ta;
	 			parameters[2]=tb;
	 			parameters[3]=tc;
	 			parameters[4]=tag;
	 			parameters[5]=tbg;
	 			parameters[6]=tcg;
	 			if(fixzg==0) parameters[7]=tzg;
	 			
	 			
	 			if(k<list.size()){
				 	if(g<0)  khi=chiNG(ref[0][j],spec); 
				 	if(g>=0) khi=chi(ref[0][j],spec,gx[0][g]);  
				}
				else khi=chiNSN(spec,gx[0][g]);
	 			
	 			//cout << "MCMC legjobb khi: " << parasave[8] << endl;
			 }
	
	
	
			 // Steepest step method (grad)
			 if(i> (g<0?lim1:lim2) ){
	

				 /*parameters[0]=z;
	 			 parameters[1]=a;
	 			 parameters[2]=b;
	 			 parameters[3]=c;
	 			 parameters[4]=ag;
	 			 parameters[5]=bg;
	 			 parameters[6]=cg;*/
	

	   
	   
	   			 //deriv: f(10**x)
				 //f'(10**x)*10**x*ln(10)
				 /*parameters.size()*/
				 
				 // calculating the gradient numerical!
				 for(int e=(k<list.size()?0:4);  e<(g<0?4: (fixzg==0?8:7) );  e++){
					double temp1=parameters[e];
					double temp2;
					double temp3;
					
					// step forward
					if(e!=3 && e!=6)  parameters[e]+=dstep*exp(-q/gradtau);  //lin
					else              parameters[e]*=pow(10,dstep*exp(-q/gradtau) );  //log (c es cg)
					
					 if(k<list.size()){
					 	if(g<0)  temp2=chiNG(ref[0][j],spec); 
	    			 	if(g>=0) temp2=chi(ref[0][j],spec,gx[0][g]);  
	    			 }
	    			 else temp2=chiNSN(spec,gx[0][g]);
	     
	    			// step backward
	    			//temp3=tkhi;
					parameters[e]=temp1;
					if(e!=3 && e!=6)  parameters[e]-=dstep*exp(-q/gradtau);
					else              parameters[e]*=pow(10,-1.*dstep*exp(-q/gradtau) );
					 
					 if(k<list.size()){
					 	if(g<0)  temp3=chiNG(ref[0][j],spec); 
	    			 	if(g>=0) temp3=chi(ref[0][j],spec,gx[0][g]);
	    			 }
	    			 else temp3=chiNSN(spec,gx[0][g]);
	    			 
	    
	    			if(tkhi<=temp2 && tkhi<=temp3) grad[e]=0;  // if the actual point is a minimum, the the derivate is 0
	    			else{
						if(e!=3 && e!=6)  grad[e]=( sqrt(temp2)-sqrt(temp3) ) / (2*dstep*exp(-q/gradtau) ) / norman;
						else 			  grad[e]=( sqrt(temp2)-sqrt(temp3) ) / (pow(10,2*dstep*exp(-q/gradtau) )) / norman;
						}
					parameters[e]=temp1; // reset parameter back to original value
	  
    			 }
    			 
    			 
	 			 // executeing the steps
	 			 if(k<list.size()){
	    		 	parameters[0]=tz=z=z-dstep*exp(-q/gradtau)*sgn(grad[0]);
	    		 	parameters[1]=ta=a=a-dstep*exp(-q/gradtau)*sgn(grad[1]);
	    		 	parameters[2]=tb=b=b-dstep*exp(-q/gradtau)*sgn(grad[2]);
	    		 	parameters[3]=tc=c=c*pow(10,-dstep*exp(-q/gradtau) *sgn(grad[3]) );
	    		 	}
	    		 if(g>=0){
	    		   parameters[4]=tag=ag=ag-dstep*exp(-q/gradtau)*sgn(grad[4]);
	    		   parameters[5]=tbg=bg=bg-dstep*exp(-q/gradtau)*sgn(grad[5]);
	    		   parameters[6]=tcg=cg=cg*pow(10,-dstep*exp(-q/gradtau) *sgn(grad[6]) );
	    		   if(fixzg==0) parameters[7]=tzg=zg=zg-dstep*exp(-q/gradtau)*sgn(grad[7]);
	    		   }
	   
	   
	 
	 
	 
	 
	   
	   
	 
	 
	 


				 // limits (priori)
				 if(k<list.size()){
				 if(b<bmin) parameters[2]=tb=b=bmin;
				 if(b>bmax) parameters[2]=tb=b=bmax;
				 if(a<amin) parameters[1]=ta=a=amin;
				 if(a+b<-1) parameters[1]=ta=a=-1-b;
				 if(b<-2 && a<0.25*b*b) parameters[1]=ta=a=0.25*b*b;
				 if(a>amax) parameters[1]=ta=a=amax;
				 if(c<c0/10000) parameters[3]=tc=c=c0/10000;
				 if(c>c0*10000 ) parameters[3]=tc=c=c0*10000;
				 if(z<z0) parameters[0]=tz=z=z0;
				 if(z>zmax) parameters[0]=tz=z=zmax;
				 }
	 
				 // GX
				 if(g>=0){
	   				if(bg<bmin) parameters[5]=tbg=bg=bmin;
	   				if(bg>bmax) parameters[5]=tbg=bg=bmax;
	   				if(ag<amin) parameters[4]=tag=ag=amin;
	   				if(ag+bg<-1) parameters[4]=tag=ag=-1-bg;
	   				if(bg<-2 && ag<0.25*bg*bg) parameters[1]=tag=ag=0.25*bg*bg;
	   				if(ag>amax) parameters[4]=tag=ag=amax;
	   				if(cg<c0g/10000) parameters[6]=tcg=cg=c0g/10000;
	   				if(cg>c0g*10000 ) parameters[6]=tcg=cg=c0g*10000;
	   				
	   				if(fixzg==0) if(zg<z0) parameters[7]=tzg=zg=z0;
				 	if(fixzg==0) if(zg>z) parameters[7]=tzg=zg=z;
	   				}
	   
	   
				 // gradient absolute length
				 //double tempt=0;
				 //for(int e=0;e<(g<0?4:7);e++) tempt+=h(grad[e]); 
	   
				 if(k<list.size()){
				 	if(g<0)  khi=chiNG(ref[0][j],spec); 
				 	if(g>=0) khi=chi(ref[0][j],spec,gx[0][g]);  
				}
				else khi=chiNSN(spec,gx[0][g]);
				 //if(tempt==0) khi=tkhi;
	 
	
	 
				 // dynamical step change: if the new step leads to a bigger (or equal) khi, then decrease the step length
				 if( int(khi*10+0.5)>=int(tkhi*10+0.5) ){ q=q+2+RandomDouble(0.,3.); /*q=q+2+int(RandomDouble(0.,3.)+0.5);*/ }
				 if( int(khi*10+0.5)<int(tkhi*10+0.5)  ){ q=q-1; }
				 if(q<0) q=0;
				 tkhi=khi;
				 if(q>gradmax*gradtau){ 
				      //cout <<"Converge at: " << (g<0?lim1:lim2) <<"+"<<i-(g<0?lim1:lim2) << endl; 
					  break; }  // criteria to exit the searcing
				 
	 
	 
	 
				 //cout << a <<"   "<< b <<"   "<< log10(c) <<endl;
				 //cout << setprecision(6) << fixed << khi <<"  q= " << q << "   grad= " << sqrt(tempt) << endl;
				 //cout << khi <<"  q= " << q << "   grad= " << sqrt(tempt) << endl;
				 //for(int e=0;e<(g<0?4:7);e++) cout <<e << "  grad= " << grad[e] << "  ";
				 //cout << endl;
	 
	 
	 
			 } // end of grad
	
	
	
	
	
	
	
	
	
			 // legjobb ertekek memorizalasa, store the best values
			 if(khi<khiS[k]){
				khiS[k]=khi; khiNS[k]=sqrt(khi/norman); tS[k]=t; aS[k]=a; bS[k]=b; cS[k]=c; zS[k]=z; targx[k]=g;      tj=j;
				if(g>=0){ agS[k]=ag; bgS[k]=bg; cgS[k]=cg;  if(fixzg==0) zgS[k]=zg; }
				//cout << khiNS[k] << "   " << c << endl;
				}
	
	
	
			 // aktualis (t-re) legjobb ertek memorizalasa, store the best values at current t 
			 if(khi<parasave[7]){
				parasave[7]=khi;  parasave[8]=sqrt(khi/norman);  parasave[1]=a; parasave[2]=b; parasave[3]=c; parasave[0]=z;   //targx[k]=g;      
				if(g>=0){ parasave[4]=ag;   parasave[5]=bg; parasave[6]=cg;   if(fixzg==0) parasave[9]=zg; }
				//cout << khiNS[k] << "   " << c << endl;
				}
	
	

	
	
	
			 //test kiiras
			 /*
			 if(k==0) ftemp << a <<" \t  "<< ta <<" \t  "<< b <<"  \t "<< tb <<" \t  "<< log10(c) <<"  \t  "<< log10(tc) <<"  \t  " 
			 << ag <<" \t  "<< tag <<" \t  "<< bg <<"  \t "<< tbg <<" \t  "<< log10(cg) <<"  \t  "<< log10(tcg) <<"  \t  " 
			 <<khi<<"  \t  " <<tkhi   <<"  \t  " << parasave[7] <<"  \t  " << q <<endl;
			 */
			 //cout << q << endl;
			

	
 		 } // end of minimalisation
 		 
 		 
 		 

	 } // end of galaxy gridding





   cout << "t= " << int(t) << "      Khi_N= " << parasave[8] << endl << endl;
   fki << "t= " << int(t) << "      Khi_N= " << parasave[8] << endl;

   // steps in time
   if(dt!=1 || o>0) j++;
   else j=j+5;

   }  // end of time gridding



  cout << endl << endl;






  //kepernyore iras, writing to standard output
  if(k<list.size()) cout << list[k].c_str() << endl;
  else cout << "Galaxy only" << endl;
  if(targx[k]<0) cout << "without added galaxy values" << endl;
  else cout << "with:  " << list2[targx[k]].c_str() << endl;
  cout << "best: " << endl;
  cout << "khi: " << khiS[k] << "  \t  khi_norm: " << khiNS[k] <<endl;//<< "   L: " << exp(-0.5*khiS[k]) << endl; 
  if(k<list.size()){
  	cout << "t: " << int(tS[k]) << endl; 
  	//cout << "t: " << inverztar(tS[0],dt,t0,tmax) << endl;
  	cout << "a: " << aS[k] << endl; 
  	cout << "b: " << bS[k] << endl; 
  	cout << "c: " << cS[k] << endl; 
  	cout << "z: " << zS[k] << endl << endl;
  }

  if(targx[k]>=0){
	  parameters[0]=zS[k];
	  parameters[1]=aS[k];
	  parameters[2]=bS[k];
	  parameters[3]=cS[k];
	  parameters[4]=agS[k];
	  parameters[5]=bgS[k];
	  parameters[6]=cgS[k];
	  if(fixzg==0) parameters[7]=zgS[k];
	  g=targx[k];
	  cout << "ag: " << agS[k] << endl; 
	  cout << "bg: " << bgS[k] << endl; 
	  cout << "cg: " << cgS[k] << endl; 
	  cout << "zg: " << zgS[k]; if(fixzg!=0) cout << "   (fixed) "; cout << endl; 
	  cout << "Galaxy/Model flux ratio:  "; 
	  if(k<list.size()) tSNR[k] = SNR(ref[0][tj],gx[0][targx[k]]) ; 
	  else tSNR[k] = -1;
	  cout << tSNR[k] << endl;
	  cout << endl; 
      }
    
    cout << endl;




  /*Ki irja a legjobb fitet, a spektrumot*/
  /* Writing out the best fit spectra  */
  // REF konyvtarat nem hozza letre magatol, se a nugent-et!
  // REF and nugent folder not created autmatically!
  if(k<list.size()){
	
	  std::ofstream fkispec;
	  //fkispec.open(("REF/"+list[k]+".dat2").c_str(),std::ofstream::out);
	  
	  if(argc>1){ sor=argv[1]; sor="REF/"+sor+"."+list[k]+".dat2"; fkispec.open((sor.c_str()),std::ofstream::out); }
      else{        fkispec.open(("REF/has.dat."+list[k]+".dat2").c_str(),std::ofstream::out); } 
	  
	  std::ofstream fkispec2;
	  std::ofstream fkispec3;
	  //fkispec2.open(("REF/"+list[k]+".dat2s").c_str(),std::ofstream::out);
	  //fkispec3.open(("REF/"+list[k]+".dat2g").c_str(),std::ofstream::out);
	  
	  if(argc>1){ sor=argv[1]; sor="REF/"+sor+"."+list[k]+".s.dat2"; fkispec2.open((sor.c_str()),std::ofstream::out); }
      else{        fkispec2.open(("REF/has.dat."+list[k]+".s.dat2").c_str(),std::ofstream::out); } 
      
      if(argc>1){ sor=argv[1]; sor="REF/"+sor+"."+list[k]+".g.dat2"; fkispec3.open((sor.c_str()),std::ofstream::out); }
      else{        fkispec3.open(("REF/has.dat."+list[k]+".g.dat2").c_str(),std::ofstream::out); } 
	  
      
      

	
	  double l;
	  g=targx[k];
	  z=zS[k];
	  zg=zgS[k];
	
	
	  /*
	  double z=parameters[0];
	  double a=parameters[1];
	  double b=parameters[2];
	  double c=parameters[3];
	  double ag=parameters[4];
	  double bg=parameters[5];
	  double cg=parameters[6];
	  */
	
	  double lalso;
	  double lfelso;

	  lalso=(l0/(1.+z)>specl0?l0/(1.+z):specl0);
	  lfelso=(lmax/(1.+z)<speclmax?lmax/(1.+z):speclmax);
	
	  //   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	  // khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	  // with galaxy
	  if(g>=0) for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		  fkispec << l << " \t " << (
		     ( fl(1-l/lfelso,aS[k],bS[k],cS[k]) )*     ref[0][tj][inverztar(l*(1.+z),dl,l0,lmax)]
		  +  ( fl(1-l/lfelso,agS[k],bgS[k],cgS[k]) )*  gx[0][g][inverztar(l*(1.+zg),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] )/err
		  << endl;	
	      }
	  // without galaxy
	  else for(l=lalso  ;l<=lfelso  ; l=l+specdl){
	  	  fkispec << l << " \t " <<
		     ( fl(1-l/lfelso,aS[k],bS[k],cS[k]) )*     ref[0][tj][inverztar(l*(1.+z),dl,l0,lmax)]/err
		  //+  ( fl(1-l/lfelso,agS[k],bgS[k],cgS[k]) )*  gx[0][g][inverztar(l*(1.+z0),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])]
		  << endl;
		  }
		  
		  
		  
		  
		  
	  // spectra
	  for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		  fkispec2 << l << " \t " << (
		     ( fl(1-l/lfelso,aS[k],bS[k],cS[k]) )*     ref[0][tj][inverztar(l*(1.+z),dl,l0,lmax)]
		  //+  ( fl(1-l/lfelso,agS[k],bgS[k],cgS[k]) )*  gx[0][g][inverztar(l*(1.+zg),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] 
		  )/err
		  << endl;	
	      }
	      
	      
	  // galaxy
	  if(g>=0) for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		  fkispec3 << l << " \t " << (
		  //   ( fl(1-l/lfelso,aS[k],bS[k],cS[k]) )*     ref[0][tj][inverztar(l*(1.+z),dl,l0,lmax)]
		    ( fl(1-l/lfelso,agS[k],bgS[k],cgS[k]) )*  gx[0][g][inverztar(l*(1.+zg),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] 
		  )/err
		  << endl;	
	      }
	      
	

	  fkispec.close();
	  
	  fkispec2.close();
	  fkispec3.close();
	
  }
  else{
  	
  	  std::ofstream fkispec;
	  //fkispec.open(("REF/Nugent_template/GalaxyOnly.dat.dat2"),std::ofstream::out);
	  
	  if(argc>1){ sor=argv[1]; sor="REF/"+sor+".GalaxyOnly.dat2"; fkispec.open((sor.c_str()),std::ofstream::out); }
      else{        fkispec.open("REF/has.dat.GalaxyOnly.dat2",std::ofstream::out); } 
      

	
	  double l;
	  g=targx[k];
	  //z=zS[k];
	  zg=zgS[k];
	
	
	  /*
	  double z=parameters[0];
	  double a=parameters[1];
	  double b=parameters[2];
	  double c=parameters[3];
	  double ag=parameters[4];
	  double bg=parameters[5];
	  double cg=parameters[6];
	  */
	
	  double lalso;
	  double lfelso;

	  lalso=(specl0);
	  lfelso=(speclmax);
	
	  //   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas
	  // khi osszeadas: referencia+galaxis-meres, majd negyzetre emel es leoszt a hiba negyzettel
	  // with galaxy
	  for(l=lalso  ;l<=lfelso  ; l=l+specdl){
		  fkispec << l << " \t " << (
		  //   ( fl(1-l/lfelso,aS[k],bS[k],cS[k]) )*     ref[0][tj][inverztar(l*(1.+z),dl,l0,lmax)]
		    ( fl(1-l/lfelso,agS[k],bgS[k],cgS[k]) )*  gx[0][g][inverztar(l*(1.+zg),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] )/err
		  << endl;	
	      }

	

	  fkispec.close();
  	
  }



}
// end of SN modell and minimalisation loop







































//sqrt( khiS[k]) * khiNS[k]  =  khi / sqrt(N)


// OSSZEGZO KIMENET IRASA
// WIRITING OUT THE SUMMARISING OUTPUT
double szum=0;
double best=1e300;
double temp;
int tempe;
int rendezo[list.size()+1];
double rendezoS[list.size()+1];
for(int k=0;k<list.size()+1;k++) rendezoS[k]=khiS[k];
for(int k=0;k<list.size()+1;k++) rendezo[k]=k;

//search for best khi value
//for(int k=0;k<list.size();k++) if(khiS[k]<best) best=khiS[k];
//for(int k=0;k<list.size();k++) if(h(khiNS[k])<best) best=h(khiNS[k]);
for(int k=0;k<list.size()+1;k++) if( sqrt( khiS[k]) * khiNS[k] <best) best= sqrt( khiS[k]) * khiNS[k] ;

//summing the likelihood (to normalise)
//for(int k=0;k<list.size();k++) szum=szum+ exp(-0.5*khiS[k]+0.5*best) ;
//for(int k=0;k<list.size();k++) szum=szum+ exp(-0.5*h(khiNS[k])+0.5*best) ;
for(int k=0;k<list.size()+1;k++) szum=szum+ exp(-0.5*( sqrt( khiS[k]) * khiNS[k] )+0.5*best) ;
//cout << "norman: " <<norman << endl;





 // arrange in order by khi
 for(int j=0;j<list.size()+1;j++){	 
			for(int k=0;k<list.size()-1+1;k++)					  
                    	if(rendezoS[k]>rendezoS[k+1]) 
						    { 
							  temp=rendezoS[k]; rendezoS[k]=rendezoS[k+1]; rendezoS[k+1]=temp; 
						      tempe=rendezo[k]; rendezo[k]=rendezo[k+1];  rendezo[k+1]=tempe;
							}
				
					
			}





// append!
fki << endl;
fki << endl <<endl << endl << "final:" << endl;
fki << "error value: " << err << endl;
if(samplederr==0) fki << "(fixed)" << endl;
else fki << "Sampled between: " << lx-nl*specdl << " and " << lx << endl;
if(uniformerr==1) fki << "Uniform error model used." << endl;
if(uniformerr==0) fki << "NON uniform error model used." << endl;
if(uniformerr==-1) fki << "WARNING! Error model not found. Uniform error model used." << endl;
if(uniformerr!=0) fki << "Applied error model to data: " << correarth << endl;
fki << "E(B-V): " << ebv << endl;
fki << "z0 (red shift null): " << z0 << endl;
fki << "z max: " << zmax << endl;
fki << "t min: " << priort0 << endl;
fki << "t max: " << priortmax << endl;
fki << "a min: " << amin << endl;
fki << "a max: " << amax << endl;
fki << "b min: " << bmin << endl;
fki << "b max: " << bmax << endl;
fki << "dstep: " << dstep << endl;
fki << "gradmax: " << gradmax << endl;
fki << "gradtau: " << gradtau << endl;
fki << "Beta: " << Beta << endl;
fki << "mcmcmax: " << mcmcmax << endl;
fki << "mcmctau: " << mcmctau << endl;
fki << "mcmclimit1: " << mcmclimit1 << endl;
fki << "mcmclimit2: " << mcmclimit2 << endl;
fki << endl;

//for(int j)
//for(int k=0;k<list.size();k++){
for(int i,k=0;i<list.size()+1;i++){
	k=rendezo[i];

	if(k<list.size()) fki << list[k].c_str() << endl;
	else fki << "Galaxy Only" << endl;
	if(khiS[k]==1e300) fki << "WARNING! No fit!" << endl;
	if(targx[k]<0) fki << "without added galaxy values" << endl;
	else fki << "with:  " << list2[targx[k]].c_str() << endl;
	fki << "khi= " << khiS[k] << endl;
	fki << "khi_norm= " << khiNS[k] << endl;
	//fki << "L= " << 100.*exp(-0.5*khiS[k]+0.5*best)/szum <<"% "<< endl; 
	//fki << "L= " << 100.*exp(-0.5*h(khiNS[k])+0.5*best)/szum <<"% "<< endl; 
	fki << setprecision(2) << fixed << "L= " << 100.*exp(-0.5*( sqrt( khiS[k]) * khiNS[k]  )+0.5*best)/szum <<"% "<< endl;  // << setprecision(5) << fixed
	fki.unsetf(ios::fixed);
	fki << setprecision(-1);
	if(k<list.size()){
		fki << "t= " << int(tS[k]) << endl; 
		fki << "a= " << aS[k] << endl; 
		fki << "b= " << bS[k] << endl; 
		fki << "c= " << cS[k] << endl; 
		fki << "z= " << zS[k] << endl;
		}
	
	if(targx[k]>=0){
		fki << "ag= " << agS[k] << endl; 
		fki << "bg= " << bgS[k] << endl; 
		fki << "cg= " << cgS[k] << endl;
		fki << "zg= " << zgS[k]; if(fixzg!=0) fki << "   (fixed) "; fki << endl;
		fki << "Galaxy/Model flux ratio=  "; 
		fki << tSNR[k] << endl; 
		fki << endl << endl; 
		}
	
	}
	
fki << endl << endl << endl;

fki.close();


// tSNR[k] : galaxy/SN model flux ratio
// SN fluxus relativ hozzajarulas  ( SNflux/(SNflux+GXflux) ):
//  1/(tSNR[k]+1) 

















/*TEMP?: megvaltoztatjuk a bemenetet, es msot ki is irjuk hogy konzisztens legyen*/
// write off the omorised data too 
{
	
	std::ofstream fkispec;
	//fkispec.open("specki.dat",std::ofstream::out);
	//fkispec.open(("REF/Nugent_template/specki.dat.dat2"),std::ofstream::out);
	
	if(argc>1){ sor=argv[1]; sor="REF/"+sor+".ref.dat2"; fkispec.open((sor.c_str()),std::ofstream::out); }
    else{        fkispec.open("REF/has.dat.ref.dat2",std::ofstream::out); }

	
	double l;
	
	/*
	double z=parameters[0];
	double a=parameters[1];
	double b=parameters[2];
	double c=parameters[3];
	double ag=parameters[4];
	double bg=parameters[5];
	double cg=parameters[6];
	*/
	
	double lalso;
	double lfelso;
	lalso=(specl0);
	lfelso=(speclmax);
	
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas

	for(l=lalso  ;l<=lfelso  ; l=l+specdl){

		fkispec << l << " \t " << spec[inverztar(l,specdl,specl0,speclmax)]/err    << endl;
		
		//if(l==1000){  cout << l << " \t " << spec[inverztar(l,specdl,specl0,speclmax)]  << endl; getchar(); }
	}

	
	

	fkispec.close();
	
}








//http://faculty.cs.niu.edu/~mcmahon/CS241/c241man/node83.html
{

int elso=0;

std::ofstream fkibulk;

sor="bulk.outall";


std::ifstream temp(sor.c_str());
//temp.open(sor,std::ifstream::in);
if(!temp.is_open()){ elso=1; }
else{ elso=0; }
temp.close();



fkibulk.open(sor.c_str(),std::ofstream::app);

if(elso==1){
	fkibulk << "# chi^2 \t chi^2_N \t GX/SN flux ratio \t GX type" << endl;
	fkibulk << "# Name \t ";
	for(int k=0;k<list.size();k++) fkibulk << list[k].c_str() << " \t ";
	fkibulk << "Galaxy Only " << endl;
}


if(argc>1){ sor=argv[1]; }
       else{ sor="has.dat"; } 

fkibulk << sor << " \t ";


for(int k=0;k<list.size()+1;k++){
	fkibulk << showpoint << setprecision(6) << khiS[k] << " \t ";
}

fkibulk << " \t ";

for(int k=0;k<list.size()+1;k++){
	fkibulk << showpoint << setprecision(6) << khiNS[k] << " \t ";
}

fkibulk << " \t ";

for(int k=0;k<list.size()+1;k++){
	if(tSNR[k]!=-1) fkibulk << showpoint << setprecision(6) << tSNR[k] << " \t ";
	else            fkibulk << fixed << setprecision(0) << tSNR[k] << " \t ";
}

fkibulk << " \t ";

for(int k=0;k<list.size()+1;k++){
	if(targx[k]<0) fkibulk << "no galaxy" << " \t ";
	else fkibulk << list2[targx[k]].c_str() << " \t ";
}


fkibulk << endl;

fkibulk.close();

}





/*
{
	
	std::ofstream fkispec;
	fkispec.open("alt.dat",std::ofstream::out);

	
	double l;
	
	
	//double z=parameters[0];
	//double a=parameters[1];
	//double b=parameters[2];
	//double c=parameters[3];
	//double ag=parameters[4];
	//double bg=parameters[5];
	//double cg=parameters[6];
	
	
	double lalso;
	double lfelso;
	lalso=(specl0);
	lfelso=(speclmax);
	g=1;
	
	//   tomb:    elso: 0/1 ertek vagy ido?    masodik: t szerinti tarolas     harmadik: hullamhossz szerinti tarolas

	for(l=lalso  ;l<=lfelso  ; l=l+specdl){

		fkispec << l << " \t " << ( 0.01*spec[inverztar(l,specdl,specl0,speclmax)]    
		+  ( fl(1-l/lfelso,0,0,3) )*  gx[0][g][inverztar(l*(1.+z0),gx[1][g][11]-gx[1][g][10],gx[1][g][0],lgxmax[g])] )
		<< endl;
		
		//if(l==1000){  cout << l << " \t " << spec[inverztar(l,specdl,specl0,speclmax)]  << endl; getchar(); }
	}

	
	

	fkispec.close();
	
}
*/




//delete ref;    
cout << endl;
printf("\n\ndone\n");  
//getchar();
return 0;   
    
}




/*
https://www.eso.org/sci/facilities/paranal/decommissioned/isaac/tools/spectroscopic_standards.html
telluric IR



#/bin/bash
# computes reddened fluxes from Fitzpatrick & Massa astro-ph 0705.0154
# usage: sedex_der.a file EBV
#   where file has two columns: lambda  flux
#
grep -v '#' $1 >sedex.in
awk 'BEGIN{ EBV='$2'; RV=3.1; o1=2.055; o2=1.322; o3=0.00; kir=1.057;
     x0=4.592;g=0.922;c1=-0.175;c2=0.807;c3=2.991;c4=0.319;c5=6.097;}
     {x=1e4/$1; 
     if($1<3300) { DD = x*x/( (x*x-x0*x0)^2 + (g*x)^2);
                  if(x<=c5){ k = c1+c2*x+c3*DD; } else {k = c1+c2*x+c3*DD+c4*(x-c5)^2; }
		 }
     if($1>=3300 && $1<4000) { k = o1+(o2-o1)/(4000.0-3300.0)*($1-3300); }
     if($1>=4000 && $1<=5500) { k = o2+(o3-o2)/(5430.0-4000.0)*($1-4000); }
     if($1>5500) { k = kir*(1e4/$1)^1.84 - RV; } 
     A = EBV * (RV + k);
     f=$2; m0=-2.5*log(f)/log(10)-A; f0=10^(-0.4*m0);
     print $1, $2, m0, f0;
}' sedex.in 
rm sedex.in






az IRAF-ban van légtömeg-korrekció a onedspec csomagban:
https://iraf.net/irafhelp.php?val=calibrate&help=Help+Page
ezzel próbáld meg kiszámolni a talált levegotömegre.

Így kell eljárnod (t.f.h. a spektrum az sp.dat file-ban van):
rspect sp.dat sp.fits
calibrate sp.fits sp_cal.fits extinct+ flux- extinct=/iraf/iraf/noao/lib/onedstds/kpnoextinct.dat airmass=1.182 exptime=<exptime a fits headerbol>
wspectext sp_cal.fits sp_cal.dat
így az ererdmény az sp_cal.dat file-ba kerül.
*/

