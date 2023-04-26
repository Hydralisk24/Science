#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <fstream>
#include <string>
#include <sstream> 

#define PI 3.14159265358979

using namespace std;


/*
Ez a program a beviteli magnitudokat integralja ki, hogy bolometrikus fenygorbe legyen belole.
Supported:
sdss: g r i z u
johnson: V R I B U
Swift: w2 m2 w1 u b v
johnson: J H K
sdss es johnson kozott a mod-dal kell valtani. johnson JHK-ra nem vonatkozik
U u kesobb lett belerakva, ezert kezeli a program utolsokent

- johnson/sdss l-gorbe.txt a bevitel, JD mag magerr magerr formaban! l helyere a szuro neve (g-gorbe.txt pl)
A 2. mag err a valodi hiba!   kulso==0!
Ha kulso==1 akkor all-gorbe.txt es ebben legyenek egymas mellett a johnson/sdss szurok
JD u g r i z , JD U B V R I !   es ertelem szeruen szuro mellett a hiba
- Swift swiftbe.txt a bevitel, JD w2 w2err m2 m2err w1 w1err u uerr b berr v verr
formatumban kell szerepelnie. Ha valamelyik nincs, akkor oda irj 0/-1 eket!
Ha nincs Swift meres/adat, akkor egyszeruen kihagyja
- JHK JHK-gorbe.txt a bevitel, JD J Jerr H Herr K Kerr 
formatumban kell szerepelnie. Ha valamelyik nincs, akkor oda irj 0/-1 eket!
- 1 sdss/johnson-ra is lefut, a tobbit kihagyja. 
0-val bar lefut, de HIBAS LESZ!!! mert a JD-t innen szamitja

JD-ket sorba rendezi mindegyiknel.
Elobb a johnson/sdss -el dolgozik, rendezi, egy adat tombe rakja
(TILTAS-sal kiveheto az ami nem kell)
Ennek a JD-je lesz az etalon, ezt fogja a tobbinel hasznalni JD[] !!!
Interpolal azokra az adatok, ahol van szomszedos adat, de nincs adat (interpol() )
Interpol: y=ym+ (x-xm) * (yp-ym) / (xp-xm)
Megadhato tiltott sav erre

Fluxust szamol: 
johnson: flux_l=10.0^(-0.4*(mag+48.60-nedl))*300.0/( (l*1e-8)*(l*1e-8) )   (JHK is)
sdss   : flux_l=10.0^(-0.4*(mag+K_l  -nedl))
l hullamhossz, nedl ehhez a hullamhosszhoz tartozo hullamhossz, K l fuggo konstants (Vega mag)
Berakja adatflux adat tipusba

Swift JD rendezes
Etalon JD[]-re interpolalja az adatokat!
+2 napot meg extrapolal, a tobbi napot ervenyteleniti
Fluxust szamol:
flux=texp( -0.4* (mag+ZP-(C*EBV) ) );   ZP pozitiv, C pozitiv, texp 10^()
Berakja adatflux adat tipusba

JHK is van rendezes
interpolal==1-nal etalon JD[]-re interpolal
interpolal==0-nal etalon JD[]-vel megegegyezo napokat rak bele
flux szamol (fenti keplet), es berakja az adatfluxba

Integralas, bolint(), bolhib()
Elotte l szerint rendezi az adatokat
Csak akkor csinalja ha van legalabb 3 szuro (ez most mar allithato)
A fuggveny trapezoid modszerrel integral a szurok kozott.
IR:
- Rayleigh-Jeans kozelitessel
- Fekete test integral
UV: pedig empirikus kozelitessel veszi figyelembe:
- Legrovidebb l es lecsap kozott derkszogu 3szog terulete
- lecsapig extrapolalja a fluxust (2 legrovidebbol), aztan 0, ennek a terulete
- Fekete test integral
Hiba hibaterjedessel, FTS numerikusan
Ez utan bolometrikus magniba konvertal
EZT A SZEKCIONAL LEHET ALLITANI, KICSIT LEJEBB!!!!
! parameter fileban mar allithato kivulrol

Mindig a letezo szuroket integralja ki, ha 3 van akkor csak annyit, ha 10 akkor annyit!
Masik adat: kiirja a BC egy szint es azok hibajat is, hogy a szinbol valo szamitassal osszevetheto legyen
sdss:     BC=bol-g    szin=g-r
johnson:  BC=bol-B    szin=B-I
Az ebbol szamitja a bolt:
BC_g = 0.053 - 0.089 × (g - r) - 0.736 × (g - r)^2
BC_B = 0.004 - 0.297 × (B - I) - 0.149 × (B - I)^2
Lymann et al. 2014


Extinkcio: (nedx) R=
U:4.9
B:4.1
V:3.1
R:2.4
I:1.7
J:0.8
H:0.5
K:0.3
*E(B-V) = EBV


V3p2
Integracio fejlesztese


V3p3, UI upgrade
Adat olvasas cserelve c++ hogy a #-t is birja
Kiirashoz lett rakva felirat, angolul is
Parameter file hasznalat/beolvasas, default letrehozasa ha nincs
Igy nem a forras kodba kell irni, hanem parameter filebol lehet allitani azt amit kell
Interpolalas kicsit valtoztatva lett: kriteriumot kapott, hogy mikor ne inteprolaljon, ez allithato
Integracional szepitettem kicsit, foleg hogy allithato legyen melyik modell legyen hasznalva
Mason nem valtoztattam, a lenyeg ugyanaz
*/



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// swift JD-hez hozza adja ezt a szamot (mert elterhet a masikhoz kepest)
//maradjon 0, az az intuitiv
//float datenull=-2400000;
float datenull=0;


//maunal
int tiltas(float JD){
int e=1;	

			// tiltás megadása: (nincs tiltás esetén 1)
	      	// if(JD=="") e=0;    "" <- ide jon a tiltott JD
			
			//if(JD==56701.500000) e=0;
			//if(JD==56784.000000) e=0;

return e;			
}




	  		// sdss griz :0 vagy Johson BVRI :1 bolometrikus integralas
	    	int mod=1;
	    	// kulso forras fajl legyen?
	    	int kulso=1; //default 1 lett, mindig all-gorbe-bol olvas igy
	    	// legyen interpolacio?
	    	int interpolal=1;
	    	/*
	    	0 eseten semmilyen adatkozi interpolacio nincs
	    	1 eseten interpolal az adatok kozott ott ahol hiany van, pl ha van V R meres akkor probal oda csinalni I merest is
	    	2 eseten force interpolacio van, lent kifejtem
	    	ha a hezeg 30 napnal tobb akkor nem interpolal kivege ha 2-es
	    	tiltott resznel (tiltom es tiltop) semmikor nem interpolal, akkor sem ha a pont amivel interpolalra atnyulik ezen a szakaszon
	    	SN plato fazis ilyen tiltott resz: nem akarjuk hogy egy plato utani ponthoz fel hasznaljon egy plato elotti pontot, mert az biztos nem lesz jo
	    	*/
	    	
	    	//interpolalas szempontbol tiltott resz
            float tiltom=-1;
            float tiltop=-1;
	  
	  
	        // NED korrekciók megadása: g/V r/R i/I z/B
	        
	        // NGC 3448
	        // g r i z
	        //float nedu=0.;
            //float nedg=0.113;
            //float nedr=0.092;
            //float nedi=0.068;
            //float nedz=0.051;
            
            // NGC 6412
	        // V R I B
	        //float nedu=0.;
            //float nedg=0.111;
            //float nedr=0.087;
            //float nedi=0.061;
            //float nedz=0.146;
            
            // M51
	        // V R I B
	        float nedu=0.152;
            float nedg=0.096;
            float nedr=0.076;
            float nedi=0.053;
            float nedz=0.127;
            
            // JHK extinkciok, azonos modon
            // M51
            float nedj=0.025;
            float nedh=0.016;
            float nedk=0.011;
            
            // swift-hez az extinkcio, B-V extinkcioja az EBV, johnson modban tulajdonkeppen nedz-nedg
            //float EBV=0.011;    // NGC 3448
            float EBV=nedz-nedg;
            
            // tavolsag megadasa:
			// NGC 6412
            //float tav=23.5  * 3.086*pow(10,22);   // Mpc --> m
            
            // M51
            float tav=7.1  * 3.086*pow(10,22);   // Mpc --> m
            
            
            //integracional hasznalt modellek
            int minn=3; // mennyinel pont kell hogy illesszen?
            int UVmod=1;
            int IRmod=1;
            /*
            UVmod=-1: semmi
            UVmod=0: levagas, 2000 A-nel nulla, e pont es elsopont kozott trapez
            UVmod=1: extrapolacio, elso ket pont alapjan linearisan extrapolal
            UVmod=2: Fekete test illesztes
            
            IRmod=-1: semmi
            IRmod=0: Rayleigh-Jeans (RJ) kozelites
            IRmod=1: Fekete test illesztes
            */
            
            string szabadls="345";
            /*
			Engedelyezett hullamhosszok indexei az FTS illeszteshez, novekvo sorrendben
			eleg csak benne szerepelni a stringben, nem kell sorrendben, es megy
			1: u/U
			2: g/B
			3: r/V
			4: i/R
			5: z/I
			6: J
			7: H
			8: K
			*/
   
            
            

            
int szabadlt[9];  // Engedelyezett hullamhosszok tomb-ben         
			
    
// felhigult fekete test (FTS) paramterei, diluted Black Body (BB) parameters        
//Dessart, Hillier 2005
//BVI
float ad=0.63241;
float bd=-0.38375;
float cd=0.28425;            








FILE *Trki;

float lsu=0.;
float lsg=0.; // g r i z hullamhossz valtozo
float lsr=0.;
float lsi=0.;
float lsz=0.;

float lju=0.;
float ljv=0.; // V R I B hullamhossz valtozo
float ljr=0.;
float lji=0.;
float ljb=0.;

float Thib=0; // T FTS hibaja
float Thibk=0; // T FTS hibaja 2



// abszolet ertek
float absz(float x){
	  if(x<0) x=-x;
	  
	  return x;	  
	  }

// negyzet	  
float h(float x){
	  return x*x;
	  }
	
// negyzet osszeg gyoke (a hibahoz)  
float h(float x, float y){
	  return sqrt(x*x+y*y);
	  }



// linearis interpolacio, kriteriummal
float interpol(float x, float xm, float ym, float xp, float yp){
	
	if(interpolal!=2) if(xp-xm>30) return -1;
	//interpolalas szempontbol tiltott idoszakok
	if(x>tiltom && x<tiltop)  return -1;
	if(x<tiltom && xp>tiltom)  return -1;
	if(x>tiltop && xm<tiltop)  return -1;

	
	  float y;
	  
	  y=ym+ (x-xm) * (yp-ym) / (xp-xm);
	  
	  return y;
	  }
	  
	  
// linearis interpolacio kritierium nelkul
float interpolf(float x, float xm, float ym, float xp, float yp){
	  float y;
	  
	  y=ym+ (x-xm) * (yp-ym) / (xp-xm);
	  
	  return y;
	  }
	  

// reverz interpolacio
float revinterpol(float y, float xm, float ym, float xp, float yp){
	  float x;
	  
	  x=xm+ (y-ym) * (xp-xm) / (yp-ym);
	  
	  return x;
	  }
	  
	  
	  
// linearis interpolacio, kriteriummal
float interpolhib(float x, float xm, float yme, float xp, float ype){
	
	if(interpolal!=2) if(xp-xm>30) return -1;
	//interpolalas szempontbol tiltott idoszakok
	if(x>tiltom && x<tiltop)  return -1;
	if(x<tiltom && xp>tiltom)  return -1;
	if(x>tiltop && xm<tiltop)  return -1;
	
	
	  float yh;
	  
	  yh=h( (1 -(x-xm)  / (xp-xm))*yme  )
	  +  h( (   (x-xm)  / (xp-xm))*ype  );
	  
	  return sqrt(yh);
	  }
	  
	  
	  
// linearis interpolacio kritierium nelkul
float interpolhibf(float x, float xm, float yme, float xp, float ype){
	  float yh;
	  
	  yh=h( (1 -(x-xm)  / (xp-xm))*yme  )
	  +  h( (   (x-xm)  / (xp-xm))*ype  );
	  
	  return sqrt(yh);
	  }
	  
	  
void interpolszol(float x, float xm, float ym, float xp, float yp, string c){
	 
	 if(xp-xm>30) cout << "I.pol.: 30nap+ lyuk!  " << c << ":  @: " << x << "      " << xm << "  " << xp << " kozt\n";
	 
	 }
	 
	  
	  
	  
	  
// Fekete Test Sugarzas, F_f es F_l is   W/m^2/m --> erg/s/cm^2/A   es  erg/s/cm^2/Hz
float FTS(float l, float T){
	  float B;
	  l=l/1e10;
	  
	  if(mod==0) B= 1000*      PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) * l*l/3e8;
	  if(mod==1) B= 1000/1e10* PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) ;
	  
	  return B;
	  }
	  
// Fekete Test Sugarzas, F_f es F_l is   W/m^2/m --> erg/s/cm^2/A   es  erg/s/cm^2/Hz    
// DOUBLE-val!!! Integralas-ra hasznalt!
double FTSd(double l, float T){
	  double B;
	  
	  if(mod==0) B= 1000*      PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) * l*l/3e8;
	  if(mod==1) B= 1000/1e10* PI* 2*3e8*3e8*6.626e-34* pow(l,-5) /(exp(3e8*6.626e-34/l/T/1.38e-23)-1) ;
	  
	  return B;
	  }
	  

// FTS integralja 
float FTSint(float s, float tig, float T){	  
	  double v=0;
	  double a;
	  double ce;
	  double ck;
	  double ch;
	  double cn;
	  double i;

	double dt=0.01e-9;
	
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

return float(1e10*v);  // mertek egyseg!
}













// flux adat
struct adatflux{
	float l;
	float flux;
	float err;
	
};


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// mely hullamhosszakra engedelyezze a FTS illesztest? (e=1 mind)
int szabadl(float l){
	//return 1;
	int e=0;
	
	if( szabadlt[1]==1 ) if(l==lsu) e=1;
	if( szabadlt[2]==1 ) if(l==lsg) e=1;
	if( szabadlt[3]==1 ) if(l==lsr) e=1;
	if( szabadlt[4]==1 ) if(l==lsi) e=1;
	if( szabadlt[5]==1 ) if(l==lsz) e=1;
	
	if( szabadlt[1]==1 ) if(l==lju) e=1;
	if( szabadlt[2]==1 ) if(l==ljb) e=1;
	if( szabadlt[3]==1 ) if(l==ljv) e=1;
	if( szabadlt[4]==1 ) if(l==ljr) e=1;
	if( szabadlt[5]==1 ) if(l==lji) e=1;

	if( szabadlt[6]==1 ) if(l==12600) e=1;
	if( szabadlt[7]==1 ) if(l==16000) e=1;
	if( szabadlt[8]==1 ) if(l==22200) e=1;
	//e=1;

	
	return e;
	}

/*
MAIN INTEGRATION OF THE FLUXES!!!

Fluxus integralas, ha nincs valahol eleg adat kilep es 0 (3 kell legalabb, allithato de 3 az alap)
adatflux-t hasznal, ez tartalmazza a hullamhosszat, fluxust, es hibat is
a program ugy csinalja, hogy a bemenetben csak a jo adatok legyenek es azok is sorban
kikommentelhetoek egyes sorok, hogy azokat ne vegye szamitsba

trapezoid modszer: ismert adatokra. 2 szuro kozti reszt integralja
nem ismert tartomanyok a trukkosek. ezekbol lehet ki kommentelgetni
Illeszt egy fekete test-t is a vegen (ez kikapcsolhato)
IR:
Rayleigh-Jeans kozelites     gyors mert analitikus, de kevesbe pontos, de elfogadhato lehet
feteke test valodi integralja    lassu, mert illeszt, de pontos, es jo felteves
UV:
levagas: legrovidebb szuro, es a lecsap kozotti terulet fele      gyors, egyszeru
extrapolacio: extrapolalja a fluxust (2 legrovidebbol) lecsapig, majd lecsaptol 0, es ennek az integralja
  gyors, valamivel eszszerubb model
fekete test integralja UV-ra      lassu, csak nagyon korai fazisban igaz, utana nem az (sok fem vonal)

fekete test illesztes: racs modszer, vegig nezi az osszeset, adaptiv lepeskozokkel
valtokozo lepeskozokkel szorastol fuggoen
T es a szorzot is illeszti (Tallgyarto.cpp analog modon)
stefan boltzmann torveny (szigma*T^4)-el szamolja ki a teljes integralt
utana Runge-Kutta modszerrel integralja l1-tol l2-ig
l2 a leghosszabb hullamhosszu szuro
l1 vagy a legrovidebb, vagy 1A, attol fuggoen UV fekete testet is akarunk, vagy valami mas modszert alkalmazunk
az integralt utana kivonjuk a teljes fekete testbol (amit a stefan boltzmannbol kaptunk)
majd utana hozza adjuk a trapezoiddal kapott erteket

szabadl megadja mely l-ekre illesszen
*/
// hol vagja le a spektrumot, ahol mar 0, at this point the SED will be rendered 0
float lecsap=2000.; // by Lyman et al. 2014
float UVhib=0;
float Tsave=12000;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
float bolint(adatflux flux[], int n){
	  if(n<minn) return 0;
	  
	  float bol=0;
	  float extra=0;
	  float lecsapc=0;
	  int i;
	  
	  extra=interpolf(lecsap , flux[0].l, flux[0].flux, flux[1].l, flux[1].flux);
	  lecsapc=revinterpol(0., flux[0].l, flux[0].flux, flux[1].l, flux[1].flux);
	  if(extra>0) lecsapc=lecsap; else extra=0;
	  
	  // 2 szuro kozotti trapezoid integral
	  for(i=0;i<n-1;i++) bol=bol+( flux[i].flux+0.5*(flux[i+1].flux-flux[i].flux) )*( flux[i+1].l-flux[i].l );
	  
	  // Rayleigh-Jeans kozelito farok int
	  if(IRmod==0) bol=bol+(flux[n-1].l*flux[n-1].flux)/3.0;    // RJ
	  
	  // Empirikus kozeleites az UV-ra, ha nagy a lecsap, akkor nem szamolja
	  UVhib=0;
      if(UVmod==0) if(lecsap<flux[0].l) UVhib=0.5*( flux[0].l-lecsap )*flux[0].flux;                   // levagas
	  if(UVmod==1) if(lecsapc<flux[0].l) UVhib=( extra+0.5*(flux[0].flux-extra) )*( flux[0].l-lecsapc );   // extrapolacio
	  bol=bol+UVhib;


	  if(IRmod>0 || UVmod>1){ // kikapcsolo
	  int j=0;
	  int k=0;
	  float e;     // C(R)
	  float s;     // szoras
	  float smag;     // szoras magnitudo
	  float sz;    // szoras tarolo
	  float szz;   // seged
	  float szzmag;   // seged magnitudo
	  float Ttemp=0; // T futo
	  float T;     // T
	  float Terr;  // T hiba
	  float Terrk;  // T hiba kisebb
	  float Tsz;   // T szoras
	  float Tszmag;   // T szoras mag
	  float Te;    // T hiba
	  float en;    // e-nek a kezdetije
	  float R;	   // sugar	
	  float Re;	   // sugar	hiba
	  float Tesz;  // Te seged
	  
	  double N;
	  double tar[2048][2];
	  for(k=0; k<2048; k++){ tar[k][0]=20100; tar[k][1]=0; }


	  // engedelyezett l, ha keves, kilep
	  for(i=0;i<n;i++) if( szabadl(flux[i].l) ) j++; ;
	  if(j<3){ Thib=0; Thibk=0; return bol; }
	  
	  szz=0;
	  //for(i=0;i<n;i++) if( szabadl(flux[i].l) ) szz=szz+1/h(2.5/log(10)*flux[i].err/flux[i].flux);
	  //for(i=0;i<n;i++) if( szabadl(flux[i].l) ) szz=szz+1/h(flux[i].err);
	  for(i=0;i<n;i++) if( szabadl(flux[i].l) ) szz=szz+1;
	  
	  szzmag=0;
	  for(i=0;i<n;i++) if( szabadl(flux[i].l) ) szzmag=szzmag+1/h(2.5/log(10)*flux[i].err/flux[i].flux);
	  
	  
	  en=-2.5*log10(flux[1].flux)+2.5*log10(FTS(flux[1].l,Tsave));

	  Tsz=10000; Tszmag=1000;
	  Ttemp=1500;
	  j=0;  
	  N=0; k=0;
	  Terr=Terrk=-1;
				 
	  while(Ttemp<20100){
				  sz=10000;					
				  e=en-9;	
				  tar[k][0]=Ttemp; tar[k][1]=0;			 
				  while(e<en+6){
 			 	 			 
							 s=0; smag=0;
							 // hiba negyzettel sulyozott szoras, legkisebb negyzetek
							 // s: flux, smag: mangitudoban, s-re illeszt, de smag-t is hasznal pl a dinamikus lepesnel (es ezt irja ki)
							 for(i=0;i<n;i++)	if( szabadl(flux[i].l) ){
 			 	 			 smag=smag+h( 2.5*log10(flux[i].flux)-2.5*log10(FTS(flux[i].l,Ttemp))+e )  /h(2.5/log(10)*flux[i].err/flux[i].flux); 
 			 	 			 s=s+h( flux[i].flux-FTS(flux[i].l,Ttemp)*pow(10,-0.4*e) )  /h(flux[i].err);
 			 	 		     }
							 
							 N=N+exp(-0.5*(s-30));
							 tar[k][1]=tar[k][1]+exp(-0.5*(s-30));
							 
							 
							 //s=sqrt(s/n);
						     s=sqrt( s/szz ); 
						     smag=sqrt( smag/szzmag );
							 //s=-2.5*log10(s)/300;
 			 	 			 
 			 	 			 
				 			 if( s<Tsz ){ T=Ttemp; Tsz=s;  Te=e; Tszmag=smag; }
				 			 if( s<sz ){ sz=s; } // adott T-n a legjobb erteket tarolja
				 			 
				 			 
				 			 e=e+smag/50;
				 			 //e=e+0.01;
				 			 
	
				 			 //if(s<1)  e=e+0.01; 
				 			 //if(s>=1) e=e+0.1;
				 			 //if(s>=8) e=e+0.9;
				 			 
				 			 
							 }
							 
				 //if(sz>1.1*Tsz)  j++;
				 //if(sz<=1.1*Tsz){ j=0; Terr=Ttemp-T; }
				 
				 //Ttemp=Ttemp+500*sz;
				 Ttemp=Ttemp+50;	 
				 
				 //if(j<11)  Ttemp=Ttemp+50; 
				 //if(j>=11) Ttemp=Ttemp+300; 
				 //printf("%f %f   \n",s,Ttemp);
				 k++;
				 }
	
				 
	  s=0;
	  for(j=k;j>=0;j--){
	  	s=s+tar[j][1]/N; 
		  //printf("%e\n",N);
	  	if(s>0.025){ Terr=tar[j][0]-T; break; }
	  }
	  
	  s=0;
	  for(j=0;j<=k;j++){
	  	s=s+tar[j][1]/N; 
		  //printf("%e\n",N);
	  	if(s>0.025){ Terrk=T-tar[j][0]; break; }
	  }
	
	  
	  // 100nal ne legyen kisebb a hiba
	  if(Terr<100) Terr=100;
	  if(Terrk<100) Terrk=100;
	  
	  if(T*0==0) Tsave=T;
	  
	  
	  //printf("%f\n",T); getchar();
	  
				 
	
	  // sugar szamolas			 
	  R=pow(10,-0.2*Te) / (ad+bd*1e4/T+cd*1e8/T/T) *tav;
      Re=h( -0.2*log(10)*pow(10,-0.2*Te)*Tszmag / (ad+bd*1e4/T+cd*1e8/T/T) ) ;
      Re=sqrt(Re)*tav;
	  			 
      printf("T:%d  \terr:%d   \tC(R):%f       szoras:%f    R:%d\n",int(T+0.5),int(((Terr>Terrk)?Terr:Terrk)+0.5),Te,Tszmag,int(R/150e9+0.5));
      fprintf(Trki,"%d  %d  %f  %f  %f  %f\n",int(T+0.5),int(Terr+0.5),Te,Tszmag,R/150e9,Re/150e9);
      Tesz=pow(10,-0.4*(Te-Tszmag));  Te=pow(10,-0.4*Te);
      
      
      
	
	  // FTS-ek hozzadasa	
	  if(IRmod<1 && UVmod>1) Ttemp=Te*FTSint(1,flux[0].l,T);                                          // FTS csak az elso szuro elottre
	  if(IRmod>0 && UVmod>1) Ttemp=Te*5.67e-8*T*T*T*T*1000-Te*FTSint(flux[0].l,flux[n-1].l,T);     // FTS mindenre ami nincs merve
	  if(IRmod>0 && UVmod<2) Ttemp=Te*5.67e-8*T*T*T*T*1000-Te*FTSint(1,flux[n-1].l,T);               // FTS csak az utolso szuro utan
	  bol=bol+Ttemp;
	  
	  // FTS hiba szamitasa
	  //Thibk=(Tesz/Te -1)*Ttemp;      Thibk=0;        
	  
	  T=T+Terr;
	  if(IRmod<1 && UVmod>1) Thib=Te*FTSint(1,flux[0].l,T)-Ttemp;                                        // FTS csak az elso szuro elottre
	  if(IRmod>0 && UVmod>1) Thib=Te*5.67e-8*T*T*T*T*1000-Te*FTSint(flux[0].l,flux[n-1].l,T)-Ttemp;     // FTS mindenre ami nincs merve
	  if(IRmod>0 && UVmod<2) Thib=Te*5.67e-8*T*T*T*T*1000-Te*FTSint(1,flux[n-1].l,T)-Ttemp;             // FTS csak az utolso szuro utan
	  
	  T=T-Terr-Terrk;
	  if(IRmod<1 && UVmod>1) Thibk=Ttemp-Te*FTSint(1,flux[0].l,T);                                        // FTS csak az elso szuro elottre
	  if(IRmod>0 && UVmod>1) Thibk=Ttemp-Te*5.67e-8*T*T*T*T*1000-Te*FTSint(flux[0].l,flux[n-1].l,T);     // FTS mindenre ami nincs merve
	  if(IRmod>0 && UVmod<2) Thibk=Ttemp-Te*5.67e-8*T*T*T*T*1000+Te*FTSint(1,flux[n-1].l,T);             // FTS csak az utolso szuro utan
	  
	  if(Thibk>Thib) Thib=Thibk;    // Nagyobbat veszi csak
	  
	  // Modell hiba !!!!  Illesztes szorasa (magban) ad egy arany elterest, ekkora +% hiba lehet meg pluszban a fluxusban
	  Thib=Thib*Thib+h(    pow(10,0.4*(Tszmag))-1    )*h(Ttemp);   
	  
	  
	  }
	  
	  
	  return bol;
	  }


// hiba: hibaterjedes: ha nincs valahol eleg adat kilep es 0 (itt is 3 kell)
// ahol nincs adat, ott 0 a hiba, igy az nincs figyelembe veve: nem okoz problemat
// adatflux-t hasznal, ez tartalmazza a hullamhosszat, fluxust, es hibat is
// a program ugy csinalja, hogy a bemenetben csak a jo adatok legyenek es azok is sorban
// kikommentelhetoek egyes sorok, hogy azokat ne vegye szamitsba
// FTS T illesztest is figyelembe veszi, de azt maga az integral szamolja, es adja at ennek
// Vigyazat, csak akkor mukodik ha a main-ban egy for cikluson bellul van a ketto, sorrendben!!!
// trapezoid integrlasasok hibaterjedesel, magyarul a hiba
float bolhib(adatflux flux[], int n){
	  if(n<minn) return 0;
	  
	  float hib=0;
	  float extrah=0;
	  float lecsapc;
	  int i;
	  
	  lecsapc=revinterpol(0., flux[0].l, flux[0].flux, flux[1].l, flux[1].flux);
	  if( interpolf(lecsap , flux[0].l, flux[0].flux, flux[1].l, flux[1].flux) >0) lecsapc=lecsap;
	  extrah=interpolhibf(lecsapc, flux[0].l, flux[0].err, flux[1].l, flux[1].err);
	  
	  
	  //2. es utolso elottit rakja ossze
	  for(i=1;i<n-1;i++) hib=hib+ h( 0.5*flux[i].err *( flux[i+1].l-flux[i-1].l ) );
	  
	  //Utolso kulon + Rayleigh-Jeans (2. sort usd ki ha nem kell)
	  hib=hib+ h(  0.5*flux[n-1].err *( flux[n-1].l-flux[n-2].l ) +  
	  (IRmod==0) * flux[n-1].l*flux[n-1].err/3.0  );   // RJ
	  
	  
	  //Elsot kulon + empirikus (2. sort usd ki ha nem kell)
	  hib=hib+ h(  0.5*flux[0].err *( flux[1].l-flux[0].l ) +
	  (UVmod==0)* 0.5*flux[0].err *(flux[0].l-lecsap)   *(lecsap<flux[0].l) +     // levagas
	  (UVmod==1)* 0.5*flux[0].err *(flux[0].l-lecsapc)   *(lecsapc<flux[0].l)    // extrapolacio
	  );
	  
	  
	  //UV flux extrapolalas hibaja
	  if(UVmod==1) hib=hib+ h(  0.5*extrah*( flux[0].l-lecsapc )  )     *(lecsapc<flux[0].l);    // extrapolacio
	  
	  
	  //printf("T_hiba:%e alap_hiba:%e\n",h(Thib),hib);
	  // T FTS hiba
	  
	  //printf("%e  %e  %e\n",h(Thib),h(Thibk),hib); getchar();
	  hib=hib+Thib+h(UVhib*0.5);
	  
	   
	   
	  return sqrt(hib);
	  }











//swift adat
struct adatswift{
	float JD;
	float w2;
	float m2;
	float w1;
	float u;
	float b;
	float v;
	
};

//+
adatswift operator+(adatswift a,adatswift b){
	adatswift c;
	
	 c.JD=a.JD;//+b.JD;
	 c.w2=a.w2+b.w2;
	 c.m2=a.m2+b.m2;
	 c.w1=a.w1+b.w1;
	 c.u=a.u+b.u;
	 c.b=a.b+b.b;
	 c.v=a.v+b.v;
	
	return c;
}

//*
adatswift operator*(adatswift a,adatswift b){
	adatswift c;
	
	 c.JD=a.JD;
	 c.w2=a.w2*b.w2;
	 c.m2=a.m2*b.m2;
	 c.w1=a.w1*b.w1;
	 c.u=a.u*b.u;
	 c.b=a.b*b.b;
	 c.v=a.v*b.v;
	
	return c;
}

//*
adatswift operator*(adatswift a,float b){
	adatswift c;
	
	 c.JD=a.JD;
	 c.w2=a.w2*b;
	 c.m2=a.m2*b;
	 c.w1=a.w1*b;
	 c.u=a.u*b;
	 c.b=a.b*b;
	 c.v=a.v*b;
	
	return c;
}

//*
adatswift operator*(float b, adatswift a){
	adatswift c;
	
    c=a*b;

    return c;	
}

//-
adatswift operator-(adatswift a,adatswift b){
	adatswift c;
	
	c=a+(-1*b);
	
	return c;
}

//10^()
// -0.4 * -1 = 0.4 esetben az ertek nem jo, tehat 0
// ez a fajta megoldas azert kell, mert kulonben nem igaz, a 9. szamjegyben elteres van
adatswift texp(adatswift a){
	adatswift c;
	
	 c.JD=a.JD;
	 
	 if(int(a.w2*10000+0.5)!=4000) c.w2=pow(10,a.w2);
	 if(int(a.m2*10000+0.5)!=4000) c.m2=pow(10,a.m2);
	 if(int(a.w1*10000+0.5)!=4000) c.w1=pow(10,a.w1);
	 if(int(a.u *10000+0.5)!=4000) c.u=pow(10,a.u);
	 if(int(a.b *10000+0.5)!=4000) c.b=pow(10,a.b);
	 if(int(a.v *10000+0.5)!=4000) c.v=pow(10,a.v);
	 
	 if(int(a.w2*10000+0.5)==4000) c.w2=0;
	 if(int(a.m2*10000+0.5)==4000) c.m2=0;
	 if(int(a.w1*10000+0.5)==4000) c.w1=0;
	 if(int(a.u *10000+0.5)==4000) c.u=0;
	 if(int(a.b *10000+0.5)==4000) c.b=0;
	 if(int(a.v *10000+0.5)==4000) c.v=0;
	
	return c;
}

//csak a legnagyobb marad meg, a tobbi 0 lesz
adatswift errormax(adatswift b){
		  adatswift a;
		  a=b;
		  float max=a.w2;
		  int s=1;
		  
		  if(max<a.m2){ max=a.m2; s=2; }
		  if(max<a.w1){ max=a.w1; s=3; }
		  if(max<a.u ){ max=a.u ; s=4; }
		  if(max<a.b ){ max=a.b ; s=5; }
		  if(max<a.v ){ max=a.v ; s=6; }
		  
		  if(s!=1) a.m2=0;
		  if(s!=2) a.m2=0;
		  if(s!=3) a.w1=0;
		  if(s!=4) a.u =0;
		  if(s!=5) a.b =0;
		  if(s!=6) a.v =0;
		  
		  
		  return a;
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
std::ifstream adat("parametersBOL.inp"); 
double read(){

	double y=-1000;
	string str;
	int i;
	int j;
	//c++ fajl beolvasas, a sort olvassa be 
	do{  std::getline(adat, str);  } while(str[0]=='#');
	

        
		j=0;	
		while(str[j]==' ') j++;
		i=j;
		while(str[i]!=' ' && i<=str.length()) i++;
		
		y = atof(   copy(str,j,i-1)    .c_str());
 
    
    return y;
}


string reads(){

	string str;
	int i;
	int j;
	//c++ fajl beolvasas, a sort olvassa be 
	do{  std::getline(adat, str);  } while(str[0]=='#');
	
        
		j=0;	
		while(str[j]==' ') j++;
		i=j;
		while(str[i]!=' ' && i<=str.length()) i++;
		
 
    
    return copy(str,j,i-1);
}


void defaultdata()			 /*Writeing default input file*/
 {
 
 FILE *f;

 f=fopen("parametersBOL.inp","wt");	
 
    fprintf(f,"%d \t[Data is? sdss griz =0 or Johson BVRI =1]\n",0);   //mod // sdss griz :0 vagy Johson BVRI :1 bolometrikus integralas	    	
	//kulso=read(); // kulso forras fajl legyen? Ha kulso akkor all-gorbe
	fprintf(f,"%d \t[Interpolate for data gaps? (2: force interpolate)]\n",1);  //interpolal // legyen interpolacio?
	//interpolalas szempontbol tiltott idoszakok
    fprintf(f,"%1.2f \t[Forbidden interpolation time start]\n",-1.);  //tiltom
    fprintf(f,"%1.2f \t[Forbidden interpolation time end]\n",-1.);  //tiltop
	// NED korrekciók megadása: u/U g/V r/R i/I z/B
	fprintf(f,"%1.2f \t[Extinction value for u/U]\n",0);  //nedu
    fprintf(f,"%1.2f \t[Extinction value for g/V]\n",0);  //nedg
    fprintf(f,"%1.2f \t[Extinction value for r/R]\n",0);  //nedr
    fprintf(f,"%1.2f \t[Extinction value for i/I]\n",0);  //nedi
    fprintf(f,"%1.2f \t[Extinction value for z/B]\n",0);  //nedz
    // JHK extinkciok, azonos modon
    fprintf(f,"%1.2f \t[Extinction value for J]\n",0);  //nedj
    fprintf(f,"%1.2f \t[Extinction value for H]\n",0);  //nedh
    fprintf(f,"%1.2f \t[Extinction value for K]\n",0);  //fprintf(f,"%f\t[Extinction value for u/U]\n",0);  nedk 
    // swift-hez az extinkcio, B-V extinkcioja az EBV, johnson modban tulajdonkeppen nedz-nedg
    fprintf(f,"%1.2f \t[E(B-V) Extinction value for SDSS]\n",0);  //EBV // only for sdss
    //tavolasg:
    fprintf(f,"%1.2f \t[Distance of the object [Mpc]]\n",10.);  //tav   // Mpc
    fprintf(f,"%d \t[MIN valid wavelength points for BB fit]\n",3);  //minn
    //illesztesi modellek
    fprintf(f,"%d \t[UV fit model: -1 ignore, 0 descent to zero at 2000 A, 1 extrapolation (from first 2 point), 2 BB fit]\n",1);  //UVmod;
    fprintf(f,"%d \t[IR fit model: -1 ignore, 0 Rayleigh-Jeans approx., 1 BB fit]\n",1);  //IRmod
    //engedelyezett hullamhosszok amire az FTS illesztes megy
    fprintf(f,"%d \t[Allowed wavelength for BB fit, 1:u/U, 2:g/B, 3:r/V, ...]\n",345);  //szabadls  

   


 fclose(f);
 	
 }
 

void data()			 /*Reading input file*/
 {
 //FILE *f;
 //f=fopen("parametersMC.inp","rt");
 
 float temp;

 if( !adat.is_open() ){printf("No parameter file! Generating a default! Exiting!\n"); defaultdata(); exit(1);}

   
   
	mod=read(); // sdss griz :0 vagy Johson BVRI :1 bolometrikus integralas	    	
	//kulso=read(); // kulso forras fajl legyen? Ha kulso akkor all-gorbe
	interpolal=read(); // legyen interpolacio?
	//interpolalas szempontbol tiltott idoszakok
    tiltom=read();
    tiltop=read();
	// NED korrekciók megadása: u/U g/V r/R i/I z/B
	nedu=read();
    nedg=read();
    nedr=read();
    nedi=read();
    nedz=read();
    // JHK extinkciok, azonos modon
    nedj=read();
    nedh=read();
    nedk=read(); 
    // swift-hez az extinkcio, B-V extinkcioja az EBV, johnson modban tulajdonkeppen nedz-nedg
    EBV=read(); // only for sdss
    //tavolasg:
    tav=read();   // Mpc
    minn=read();
    //illesztesi modellek
    UVmod=read();
    IRmod=read();
    //engedelyezett hullamhosszok amire az FTS illesztes megy
    szabadls=reads();  
    
    		/*
            UVmod=-1: semmi
            UVmod=0: levagas, 2000 A-nel nulla, e pont es elsopont kozott trapez
            UVmod=1: extrapolacio, elso ket pont alapjan linearisan extrapolal
            UVmod=2: Fekete test illesztes
            
            IRmod=-1: semmi
            IRmod=0: Rayleigh-Jeans (RJ) kozelites
            IRmod=1: Fekete test illesztes
            */
    



 //++++++++++++++++++++++++++++++++++++++++++++
 
 
 if(mod=1) if(nedz>0 && nedg>0) EBV=nedz-nedg;
 tav=tav * 3.086*pow(10,22); // Mpc --> m
 
 if(interpolal>2) interpolal=1;
 if(interpolal<0) interpolal=0;
 
 if(tiltom>tiltop){ temp=tiltom; tiltom=tiltop; tiltop=temp; } 
 
 if(minn<2) minn=2;
 
 for(int i=0;i<9;i++) szabadlt[i]=0;
 if(szabadls.find('0')!=string::npos) szabadlt[0]=1;
 if(szabadls.find('1')!=string::npos) szabadlt[1]=1;
 if(szabadls.find('2')!=string::npos) szabadlt[2]=1;
 if(szabadls.find('3')!=string::npos) szabadlt[3]=1;
 if(szabadls.find('4')!=string::npos) szabadlt[4]=1;
 if(szabadls.find('5')!=string::npos) szabadlt[5]=1;
 if(szabadls.find('6')!=string::npos) szabadlt[6]=1;
 if(szabadls.find('7')!=string::npos) szabadlt[7]=1;
 if(szabadls.find('8')!=string::npos) szabadlt[8]=1;
 
 temp=0;
 for(int i=0;i<9;i++) temp=temp+szabadlt[i];
 if(temp==0) cout << "Nincs ervenyes hullamhossz az illeszteshez / No valid wavelength for BB fit" << endl;
 
 
 //cout << szabadls << endl;
 //for(int i=0;i<9;i++)  cout << i <<"  " << szabadlt[i] << endl;   getchar();

 

 //fclose(f);
 adat.close();
 }
 
 





			        







main ()
{

// ciklus vált
int i;
int j;
int k;

// paraméterek
int n;  // nem használt
int m;

// segéd számok
int szum;  // JD táblázat adatainak számát adja meg
int szumk; // fluxus táblázat adatainak száma
int csekk;
int csekkb;

// adatok száma. EOFnál megszámolja
int uszam=0;
int gszam=0;
int rszam=0;
int iszam=0;
int zszam=0;

float csere; // csere seged

string sor;  // seged string beolvasashoz


//elsõként közelitõleg az éjszakák száma, biztonsági (nagy) szám, hogy biztosan benne legyen
m=200;



data();
if(mod!=0 && mod!=1){ printf("Rossz/Bad mod parameter!\n");  exit(1); }
if(tav<=0){ printf("Tavolsag/Distance nem megfelelo/inadequate!\n");  exit(1); }

if(mod==0) printf("sdss griz mod\n\n");
if(mod==1) printf("Johnson BVRI mod\n\n");



// ide olvas be
float Tu[m][4];
float Tg[m][4];
float Tr[m][4];
float Ti[m][4];
float Tz[m][4];
for(i=0;i<m;i++) for(j=0;j<4;j++) {Tu[i][j]=0; Tg[i][j]=0; Tr[i][j]=0; Ti[i][j]=0; Tz[i][j]=0;}

// összes JD
float JD[m];
// adatflux elem taroloja
int Tn[m];
for(i=0;i<m;i++) Tn[i]=0;

// osszeset tartalmazo adat tömb
float Adat[m][11];	
float cserea[11];






// NED korrekciók, es tav helyben, de globalisan vannak megadva!!!
//float nedg;
//float nedr;
//float nedi;
//float nedz;
//float tav;


//lejebb vannak deklaralva:
//float flux[szum][9];
//adatflux fluxusok[szumk][20];
//adatflux cseref;



FILE *usz;
FILE *gsz;
FILE *rsz;
FILE *isz;
FILE *zsz;
FILE *f;
FILE *bol;
FILE *bol2;
FILE *bc;


//bemenet: 4 szürõ fényesség adatai
if(mod==0){
// ssds griz szurok hullamhosszai (A)
lju=lsu=3543.;
ljv=lsg=4770.;
ljr=lsr=6230.;
lji=lsi=7630.;
ljb=lsz=9130.;
usz=fopen("u-gorbe.txt","rb");
gsz=fopen("g-gorbe.txt","rb");
rsz=fopen("r-gorbe.txt","rb");
isz=fopen("i-gorbe.txt","rb");
zsz=fopen("z-gorbe.txt","rb");
}

if(mod==1){
// johnson VRIB szurok hullamhosszai (A)
lju=lsu=3600.;
ljv=lsg=5450.;
ljr=lsr=6410.;
lji=lsi=7980.;
ljb=lsz=4380.;
usz=fopen("U-gorbe.txt","rb");
gsz=fopen("V-gorbe.txt","rb");
rsz=fopen("R-gorbe.txt","rb");
isz=fopen("I-gorbe.txt","rb");
zsz=fopen("B-gorbe.txt","rb");
}





if(usz == NULL) printf("u/U file nincs\n");
if(gsz == NULL) printf("g/V file nincs\n");
if(rsz == NULL) printf("r/R file nincs\n");
if(isz == NULL) printf("i/I file nincs\n");
if(zsz == NULL) printf("z/B file nincs\n");


//kimenet
f=fopen("AdatSor.txt","wt");
bol=fopen("bolF-gorbe.txt","wt");
bol2=fopen("bol-gorbe.txt","wt");
bc=fopen("BC-gorbe.txt","wt");
Trki=fopen("Tr-gorbe.txt","wt"); // Trki globalisan deklaralva, kiirja a T R ertekeket is kulon ( bolint() irja ki)

fprintf(Trki,"#JD \t Temp [K]   err \t C(R) \t fit err \t R [AU] \t err \n");


if(kulso==1){ usz = NULL; gsz = NULL; rsz = NULL; isz = NULL; zsz = NULL; }


// beolvasás

i=0;
if(usz!=NULL) while (!feof(usz)){    fscanf(usz,"%f %f %f %f\n",&Tu[i][0],&Tu[i][1],&Tu[i][2],&Tu[i][3]); 
                              uszam++; i++;}
i=0;
if(gsz!=NULL) while (!feof(gsz)){    fscanf(gsz,"%f %f %f %f\n",&Tg[i][0],&Tg[i][1],&Tg[i][2],&Tg[i][3]); 
                              gszam++; i++;}
i=0;
if(rsz!=NULL) while (!feof(rsz)){	  fscanf(rsz,"%f %f %f %f\n",&Tr[i][0],&Tr[i][1],&Tr[i][2],&Tr[i][3]);  
                              rszam++; i++;}
i=0;
if(isz!=NULL) while (!feof(isz)){	  fscanf(isz,"%f %f %f %f\n",&Ti[i][0],&Ti[i][1],&Ti[i][2],&Ti[i][3]);  
                              iszam++; i++;}
i=0;
if(zsz!=NULL) while (!feof(zsz)){	  fscanf(zsz,"%f %f %f %f\n",&Tz[i][0],&Tz[i][1],&Tz[i][2],&Tz[i][3]);  
                              zszam++; i++;}

printf("\nBeolvasott r:\n");                              
for(i=0;i<rszam;i++) printf("%f %f %f %f\n",Tr[i][0],Tr[i][1],Tr[i][2],Tr[i][3]);
printf("\n"); 

//gszam++; rszam++; iszam++; zszam++; 
// Nem kell mert a beolvasás még egy nem létezõ sort is kiolvas. (eof file failorrel zárul)
// A SZÁMÁT adja meg! Nem az utolsó elem tömb számát (1 a különbség)                           

printf("beolvasott u/U=%d\n",uszam);
printf("beolvasott g/V=%d\n",gszam);
printf("beolvasott r/R=%d\n",rszam);
printf("beolvasott i/I=%d\n",iszam);
printf("beolvasott z/Z=%d\n",zszam);


if(usz!=NULL) fclose(usz);	  
if(gsz!=NULL) fclose(gsz);	
if(rsz!=NULL) fclose(rsz);
if(isz!=NULL) fclose(isz);
if(zsz!=NULL) fclose(zsz);




// itt fordul az U, es a vegere kerul, ott kezeli!


   
// JD-be íratás. Megkeresi milyen JD-k vannak 
// r/R, kezdõ
for(i=0;i<rszam;i++){
    
    JD[i]=Tr[i][0];
    
    }
szum=rszam;
  
   
printf("szum r/R-nel=%d\n",szum);



// g/V
for(i=0;i<gszam;i++){
				 
     csekk=0;
     for(j=0;j<szum;j++){
	 					 if( JD[j]-Tg[i][0]<0.26 && JD[j]-Tg[i][0]>-0.26 ){ csekk=1; break; }
	 					 }
						 
	if(csekk==0 && Tg[i][0]!=0){
				 szum=szum+1;
				 JD[szum-1]=Tg[i][0];
				 }					 				 
    }
    
printf("szum g/V-nel=%d\n",szum);
    

// i/I    
for(i=0;i<iszam;i++){
				 
     csekk=0;
     for(j=0;j<szum;j++){
	 					 if( JD[j]-Ti[i][0]<0.26 && JD[j]-Ti[i][0]>-0.26 ){ csekk=1; break; }
						 }
						 
	if(csekk==0 && Ti[i][0]!=0){
				 szum=szum+1;
				 JD[szum-1]=Ti[i][0];
				 }					 						 
    }
    
printf("szum i/I-nel=%d\n",szum);


// z/B    
for(i=0;i<zszam;i++){
				 
     csekk=0;
     for(j=0;j<szum;j++){
	 					 if( JD[j]-Tz[i][0]<0.26 && JD[j]-Tz[i][0]>-0.26 ){ csekk=1; break; }
						 }
						 
	if(csekk==0 && Tz[i][0]!=0){
				 szum=szum+1;
				 JD[szum-1]=Tz[i][0];
				 }					 						 
    }

printf("szum z/B-nel=%d\n\n",szum);


// u/U   
for(i=0;i<uszam;i++){
				 
     csekk=0;
     for(j=0;j<szum;j++){
	 					 if( JD[j]-Tu[i][0]<0.26 && JD[j]-Tu[i][0]>-0.26 ){ csekk=1; break; }
						 }
						 
	if(csekk==0 && Tu[i][0]!=0){
				 szum=szum+1;
				 JD[szum-1]=Tu[i][0];
				 }					 						 
    }

printf("szum u/U-nel=%d\n\n",szum);



// sorba rendezi a JD-t
for(i=0;i<szum;i++){
					 
					 for(j=0;j<szum-1;j++)					  
                    	if(JD[j]>JD[j+1]){ csere=JD[j]; JD[j]=JD[j+1]; JD[j+1]=csere; }
					
					}


// kiiras
/*
printf("rendezett JD:\n");
for(i=0;i<szum;i++){
					printf("%f\n",JD[i]);
					}
*/					
					

		
					
// Egységesités, Adatba vitel					
for(i=0;i<szum;i++){
					Adat[i][0]=JD[i];
					
					csekk=0;
					for(j=0;j<gszam;j++){
										 if(JD[i]==Tg[j][0]){ csekk++; Adat[i][1]=Tg[j][1];
										 Adat[i][2]=Tg[j][3]; }
										 }
				    if(csekk==0){ Adat[i][1]=-1.0;  Adat[i][2]=0.0;  }
				    if(csekk>1) printf("Hiba: tobb adat (%d) g-ben JD=%f-nel",csekk,JD[i]);
					
					
					csekk=0;
					for(j=0;j<rszam;j++){
										 if(JD[i]==Tr[j][0]){ csekk++; Adat[i][3]=Tr[j][1];
										 Adat[i][4]=Tr[j][3]; }
										 }
				    if(csekk==0){ Adat[i][3]=-1.0;  Adat[i][4]=0.0;  }
				    if(csekk>1) printf("Hiba: tobb adat (%d) r-ben JD=%f-nel",csekk,JD[i]);
					
					
					csekk=0;
					for(j=0;j<iszam;j++){
										 if(JD[i]==Ti[j][0]){ csekk++; Adat[i][5]=Ti[j][1];
										 Adat[i][6]=Ti[j][3]; }
										 }
				    if(csekk==0){ Adat[i][5]=-1.0;  Adat[i][6]=0.0;  }
				    if(csekk>1) printf("Hiba: tobb adat (%d) i-ben JD=%f-nel",csekk,JD[i]);
				    
				    
   					csekk=0;
					for(j=0;j<zszam;j++){
										 if(JD[i]==Tz[j][0]){ csekk++; Adat[i][7]=Tz[j][1];
										 Adat[i][8]=Tz[j][3]; }
										 }
				    if(csekk==0){ Adat[i][7]=-1.0;  Adat[i][8]=0.0;  }
				    if(csekk>1) printf("Hiba: tobb adat (%d) z-ben JD=%f-nel",csekk,JD[i]);
				    
				    
				    csekk=0;
					for(j=0;j<uszam;j++){
										 if(JD[i]==Tu[j][0]){ csekk++; Adat[i][9]=Tu[j][1];
										 Adat[i][10]=Tu[j][3]; }
										 }
				    if(csekk==0){ Adat[i][9]=-1.0;  Adat[i][10]=0.0;  }
				    if(csekk>1) printf("Hiba: tobb adat (%d) u-ben JD=%f-nel",csekk,JD[i]);
				    
				    
					}

				
// Alternativ bevitel +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//FILE *fall;
//fall=fopen("all-gorbe.txt","rb");
std::ifstream fall("all-gorbe.txt"); 

if(/*fall!=NULL*/ fall.is_open() && kulso==1){
i=0; 
szum=0;

/*
while (!feof(fall)){  
			   //if(mod==0) fscanf(fall,"%f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n",&Adat[i][0],&Adat[i][9],&Adat[i][10],&Adat[i][1],&Adat[i][2],&Adat[i][3],&Adat[i][4],&Adat[i][5],&Adat[i][6],&Adat[i][7],&Adat[i][8]);
			   //if(mod==1) fscanf(fall,"%f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n",&Adat[i][0],&Adat[i][9],&Adat[i][10],&Adat[i][7],&Adat[i][8],&Adat[i][1],&Adat[i][2],&Adat[i][3],&Adat[i][4],&Adat[i][5],&Adat[i][6]);
			   
			   fscanf(fall,"%f   %f  %f  %f   %f   %f   %f   %f   %f   %f   %f\n",&Adat[i][0],&Adat[i][9],&Adat[i][10],&Adat[i][7],&Adat[i][8],&Adat[i][1],&Adat[i][2],&Adat[i][3],&Adat[i][4],&Adat[i][5],&Adat[i][6]);
			   
			   szum++; i++;}*/
			   
	
	
			   
while(getline(fall,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    stringstream strs(sor);

    //strs >> mezo1;
    strs >> Adat[i][0];
    strs >> Adat[i][9];
    strs >> Adat[i][10];
    strs >> Adat[i][7];
    strs >> Adat[i][8];
    strs >> Adat[i][1];
    strs >> Adat[i][2];
    strs >> Adat[i][3];
    strs >> Adat[i][4];
    strs >> Adat[i][5];
    strs >> Adat[i][6];



   szum++; i++;

}



fall.close();

	

for(i=0;i<szum;i++){
					if(Adat[i][1]==0) Adat[i][1]=-1; 
					if(Adat[i][3]==0) Adat[i][3]=-1;
					if(Adat[i][5]==0) Adat[i][5]=-1;
					if(Adat[i][7]==0) Adat[i][7]=-1;
					if(Adat[i][9]==0) Adat[i][9]=-1;
					}

printf("All szam:%d\n\n",szum);


//rendezes
//sorba rendezi a JD-t			 
for(i=0;i<szum;i++){
					 for(j=0;j<szum-1;j++)					  
                    	if(Adat[j][0]>Adat[j+1][0]) for(k=0;k<11;k++){ cserea[k]=Adat[j][k]; Adat[j][k]=Adat[j+1][k]; Adat[j+1][k]=cserea[k]; }
					
					 }
					 
					 


}
// alternativ bevitel ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++				
				
			
					
// fájlba íratás
if(mod==0) fprintf(f,"JD             g mag       g err       r mag       r err       i mag       i err       z mag       z err       u mag      u err\n");				
if(mod==1) fprintf(f,"JD             V mag       V err       R mag       R err       I mag       I err       B mag       B err       U mag      U err\n");				

for(i=0;i<szum;i++){
					fprintf(f,"%f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n",Adat[i][0],Adat[i][1],Adat[i][2],Adat[i][3],Adat[i][4],Adat[i][5],Adat[i][6],Adat[i][7],Adat[i][8],Adat[i][9],Adat[i][10]);
										
					
					}				


// felszabaditja a beolvasott adatokat, mert ezek már a Adat-ban vannak					
//free(Tg);
//free(Tr);
//free(Ti);
//free(Tz);



		   		   	
				
// interpol(float x, float xm, float ym, float xp, float yp)					
// lináris interpoláció, NÉZD MEG!!!  Elso es utolsot nem csinalja
// Ha elso vagy utolsoknal nincs adat, azt nem interpolalja, az ugy marad
if(interpolal>=1) for(i=1;i<szum-1;i++){
					
            //tiltás:
			if(  tiltas(Adat[i][0])  )
			
					
					//g/V
					csekk=0; csekkb=0;
					if(Adat[i][1]==-1.0){  
										   while(Adat[i-1-csekk][1]==-1.0 && i-1-csekk>0) csekk++;
					                       while(Adat[i+1+csekkb][1]==-1.0 && i+1+csekkb<szum-1) csekkb++;
					                       
										   if( Adat[i-1-csekk][1]!=-1.0 && Adat[i+1+csekkb][1]!=-1.0 ){
						   	   			   interpolszol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][1],Adat[i+1+csekkb][0],Adat[i+1+csekkb][1],"g/V");
										   Adat[i][1]=interpol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][1],Adat[i+1+csekkb][0],Adat[i+1+csekkb][1]);
										   Adat[i][2]=interpolhib(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][2],Adat[i+1+csekkb][0],Adat[i+1+csekkb][2]);                                          
										   //Adat[i][2]=h(Adat[i-1-csekk][2],Adat[i+1+csekkb][2]);	
										   															}
										   }
					
					//r/R
					csekk=0; csekkb=0;
					if(Adat[i][3]==-1.0){  
                                           while(Adat[i-1-csekk][3]==-1.0 && i-1-csekk>0) csekk++;
					                       while(Adat[i+1+csekkb][3]==-1.0 && i+1+csekkb<szum-1) csekkb++;
					                       
					                       if( Adat[i-1-csekk][3]!=-1.0 && Adat[i+1+csekkb][3]!=-1.0 ){
										   interpolszol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][3],Adat[i+1+csekkb][0],Adat[i+1+csekkb][3],"r/R");
										   Adat[i][3]=interpol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][3],Adat[i+1+csekkb][0],Adat[i+1+csekkb][3]);
										   Adat[i][4]=interpolhib(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][4],Adat[i+1+csekkb][0],Adat[i+1+csekkb][4]);					                       
										   //Adat[i][4]=h(Adat[i-1-csekk][4],Adat[i+1+csekkb][4]);	
										   															}
										   }

                    //i/I
                    csekk=0; csekkb=0;
					if(Adat[i][5]==-1.0){  
										   while(Adat[i-1-csekk][5]==-1.0 && i-1-csekk>0) csekk++;
					                       while(Adat[i+1+csekkb][5]==-1.0 && i+1+csekkb<szum-1) csekkb++;
					                       
					                       if( Adat[i-1-csekk][5]!=-1.0 && Adat[i+1+csekkb][5]!=-1.0 ){
										   interpolszol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][5],Adat[i+1+csekkb][0],Adat[i+1+csekkb][5],"i/I");
										   Adat[i][5]=interpol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][5],Adat[i+1+csekkb][0],Adat[i+1+csekkb][5]);
										   Adat[i][6]=interpolhib(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][6],Adat[i+1+csekkb][0],Adat[i+1+csekkb][6]);
										   //Adat[i][6]=h(Adat[i-1-csekk][6],Adat[i+1+csekkb][6]);    
										   															}
										   }

					//z/B
					csekk=0; csekkb=0;
					if(Adat[i][7]==-1.0){  
										   while(Adat[i-1-csekk][7]==-1.0 && i-1-csekk>0) csekk++;
					                       while(Adat[i+1+csekkb][7]==-1.0 && i+1+csekkb<szum-1) csekkb++;
					                       
					                       if( Adat[i-1-csekk][7]!=-1.0 && Adat[i+1+csekkb][7]!=-1.0 ){
										   interpolszol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][7],Adat[i+1+csekkb][0],Adat[i+1+csekkb][7],"z/B");
										   Adat[i][7]=interpol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][7],Adat[i+1+csekkb][0],Adat[i+1+csekkb][7]);
										   Adat[i][8]=interpolhib(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][8],Adat[i+1+csekkb][0],Adat[i+1+csekkb][8]);
										   //Adat[i][8]=h(Adat[i-1-csekk][8],Adat[i+1+csekkb][8]);    
										   															}
										   }
										   
		            //u/U
					csekk=0; csekkb=0;
					if(Adat[i][9]==-1.0){  
										   while(Adat[i-1-csekk][9]==-1.0 && i-1-csekk>0) csekk++;
					                       while(Adat[i+1+csekkb][9]==-1.0 && i+1+csekkb<szum-1) csekkb++;
					                       
					                       if( Adat[i-1-csekk][9]!=-1.0 && Adat[i+1+csekkb][9]!=-1.0 ){
										   interpolszol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][9],Adat[i+1+csekkb][0],Adat[i+1+csekkb][9],"u/U");
										   Adat[i][9]=interpol(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][9],Adat[i+1+csekkb][0],Adat[i+1+csekkb][9]);
										   Adat[i][10]=interpolhib(Adat[i][0],Adat[i-1-csekk][0],Adat[i-1-csekk][10],Adat[i+1+csekkb][0],Adat[i+1+csekkb][10]);
										   //Adat[i][8]=h(Adat[i-1-csekk][8],Adat[i+1+csekkb][8]);    
										   															}
										   }

					
	
					
		}

	
// lineáris extrapoláció, ide csinald majd meg ha kell megint!
//(extra hiba: meredekség hibájával számitott érték 56700tõl: ah*(JD-56700) )
//if( Adat[i][0]<56784.000000 )     // meddig inter !!!				

/*		lineáris extrapoláció elemei:		
a=0.037349 +- 0.000829 (2.22%)
b=15.1459  +- 0.12     (0.08%)
-56700-el eltolva
eltolás nélkül:
b=-2102.5424
*/



// fájlba íratás
if(interpolal==0) fprintf(f,"\nNincs Interpolalas!");
if(interpolal==0) fprintf(f,"\nNo Interpolation!");
if(mod==0) fprintf(f,"\nAfter interpolation:\n");				
if(mod==1) fprintf(f,"\nAfter interpolation:\n");	
if(mod==0) fprintf(f,"Interpolalas utan:\nJD             g mag       g err       r mag       r err       i mag       i err       z mag       z err       u mag      u err\n");				
if(mod==1) fprintf(f,"Interpolalas utan:\nJD             V mag       V err       R mag       R err       I mag       I err       B mag       B err       U mag      U err\n");	
		

for(i=0;i<szum;i++){
					if(  tiltas(Adat[i][0])  )
					fprintf(f,"%f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n",Adat[i][0],Adat[i][1],Adat[i][2],Adat[i][3],Adat[i][4],Adat[i][5],Adat[i][6],Adat[i][7],Adat[i][8],Adat[i][9],Adat[i][10]);
										
					
					}
					
					

// fluxusok
float flux[szum][11];
for(i=0;i<szum;i++) for(j=0;j<11;j++) flux[i][j]=0; // nullazza



					
            



// F [w/m*m] = 10^( -0.4*(m+19.26) )
// 1 w/m*m = 10000 erg/s*cm*cm



// fluxus trafo (griz:AB magni --> erg/cm/cm/s/A), hiba: hibaterjedessel
// a *c/l/l (fenysebesseg/hullahossz^2) szorzo az atvaltas: erg/cm/cm/s/Hz -> erg/cm/cm/s/A
if(mod==0){ // szekcio sdss
csekk=0;
for(i=0;i<szum;i++){
										
					if(  tiltas(Adat[i][0])  ){
					
					flux[i-csekk][0]=Adat[i][0];
					
					
					if(Adat[i][1]!=-1.0){
					Adat[i][1]=Adat[i][1]-nedg;
					flux[i-csekk][1]= pow(10.0,-0.4*(Adat[i][1]+48.60))*300.0/( (lsg*pow(10,-8))*(lsg*pow(10,-8)) );
					flux[i-csekk][2]= pow(10.0,-0.4*(Adat[i][1]+48.60))*300.0/( (lsg*pow(10,-8))*(lsg*pow(10,-8)) ) * -0.4*log(10)*Adat[i][2];	 
					}
					
					
					if(Adat[i][3]!=-1.0){
					Adat[i][3]=Adat[i][3]-nedr;
					flux[i-csekk][3]= pow(10.0,-0.4*(Adat[i][3]+48.60))*300.0/( (lsr*pow(10,-8))*(lsr*pow(10,-8)) );	 
					flux[i-csekk][4]= pow(10.0,-0.4*(Adat[i][3]+48.60))*300.0/( (lsr*pow(10,-8))*(lsr*pow(10,-8)) ) * -0.4*log(10)*Adat[i][4];
					}
					
					
					if(Adat[i][5]!=-1.0){
					Adat[i][5]=Adat[i][5]-nedi;
					flux[i-csekk][5]= pow(10.0,-0.4*(Adat[i][5]+48.60))*300.0/( (lsi*pow(10,-8))*(lsi*pow(10,-8)) );	 
	    			flux[i-csekk][6]= pow(10.0,-0.4*(Adat[i][5]+48.60))*300.0/( (lsi*pow(10,-8))*(lsi*pow(10,-8)) ) * -0.4*log(10)*Adat[i][6];
					}
					
					
					if(Adat[i][7]!=-1.0){
					Adat[i][7]=Adat[i][7]-nedz+0.02;
					flux[i-csekk][7]= pow(10.0,-0.4*(Adat[i][7]+48.60))*300.0/( (lsz*pow(10,-8))*(lsz*pow(10,-8)) );	 
				    flux[i-csekk][8]= pow(10.0,-0.4*(Adat[i][7]+48.60))*300.0/( (lsz*pow(10,-8))*(lsz*pow(10,-8)) ) * -0.4*log(10)*Adat[i][8];
					}
					
					if(Adat[i][9]!=-1.0){
					Adat[i][9]=Adat[i][9]-nedu;
					flux[i-csekk][9]= pow(10.0,-0.4*(Adat[i][9]+48.60))*300.0/( (lsu*pow(10,-8))*(lsu*pow(10,-8)) );	 
				    flux[i-csekk][10]= pow(10.0,-0.4*(Adat[i][9]+48.60))*300.0/( (lsu*pow(10,-8))*(lsu*pow(10,-8)) ) * -0.4*log(10)*Adat[i][10];
					}
					
										
					}
					else csekk++;
             
			 }
			 
szumk=szum-csekk;
} // szekcio veg


// fluxus trafo (Johnson:Vega magni --> erg/cm/cm/s/A), hiba: hibaterjedessel
// m=0 a flux= 	U:417.5 	B:632 	V:363.1 	R:217.7 	I:112.6 	J:31.47 	H:11.38 	K:3.961   x 1e-11 erg/cm/cm/s/A
//http://www.astronomy.ohio-state.edu/~martini/usefuldata.html & https://en.wikipedia.org/wiki/Apparent_magnitude
if(mod==1){ // szekcio Johnson
csekk=0;
for(i=0;i<szum;i++){
										
					if(  tiltas(Adat[i][0])  ){
					
					flux[i-csekk][0]=Adat[i][0];
					
					
					if(Adat[i][1]!=-1.0){
					Adat[i][1]=Adat[i][1]-nedg;
					flux[i-csekk][1]= pow(10.0,-0.4*(Adat[i][1]+21.10));	 
					flux[i-csekk][2]= pow(10.0,-0.4*(Adat[i][1]+21.10)) * 0.4*log(10)*Adat[i][2];
					}
					
					
					if(Adat[i][3]!=-1.0){
					Adat[i][3]=Adat[i][3]-nedr;
					flux[i-csekk][3]= pow(10.0,-0.4*(Adat[i][3]+21.65));	
					flux[i-csekk][4]= pow(10.0,-0.4*(Adat[i][3]+21.65)) * 0.4*log(10)*Adat[i][4];
					}
					
					
					if(Adat[i][5]!=-1.0){
					Adat[i][5]=Adat[i][5]-nedi;
					flux[i-csekk][5]= pow(10.0,-0.4*(Adat[i][5]+22.37));	 
					flux[i-csekk][6]= pow(10.0,-0.4*(Adat[i][5]+22.37)) * 0.4*log(10)*Adat[i][6];
					}
					
					
					if(Adat[i][7]!=-1.0){
					Adat[i][7]=Adat[i][7]-nedz;
					flux[i-csekk][7]= pow(10.0,-0.4*(Adat[i][7]+20.50));	 
					flux[i-csekk][8]= pow(10.0,-0.4*(Adat[i][7]+20.50)) * 0.4*log(10)*Adat[i][8];
					}
					
					
					if(Adat[i][9]!=-1.0){
					Adat[i][9]=Adat[i][9]-nedu;
					flux[i-csekk][9]= pow(10.0,-0.4*(Adat[i][9]+20.95));	 
					flux[i-csekk][10]= pow(10.0,-0.4*(Adat[i][9]+20.95)) * 0.4*log(10)*Adat[i][10];
					}
					
										
					}
					else csekk++;
             
			 }
			 
szumk=szum-csekk;
} // szekcio veg

// ahol nincs adat ott 0 a fluxus es a hibaja


// fájlba íratás
if(mod==0) fprintf(f,"\nFluxusok/Fluxes: JD g r i z u es hibaik (errors) [erg/cm/cm/s/A]\n");				
if(mod==1) fprintf(f,"\nFluxusok/Fluxes: JD V R I B U es hibaik (errors) [erg/cm/cm/s/A]\n");	


for(i=0;i<szumk;i++){
					fprintf(f,"%4.3f   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e\n",flux[i][0],flux[i][1],flux[i][2],flux[i][3],flux[i][4],flux[i][5],flux[i][6],flux[i][7],flux[i][8],flux[i][9],flux[i][10]);
										
					
					}
					

//bolomnak atadott fluxusok
adatflux fluxusok[szumk][25];
adatflux cseref;
			

// adatflux tipusba taroljuk tovabb, + JD-ben az idot
// JD es Tn is itt kezdodik
//lsx formatum van, de lsx=ljx, igy van megadva
n=0;
for(i=0;i<szumk;i++){
					 JD[i]=flux[i][0];
					 
					 if(flux[i][1]!=0){ 
					 fluxusok[i][Tn[i]].l=lsg;
					 fluxusok[i][Tn[i]].flux=flux[i][1];
					 fluxusok[i][Tn[i]].err=flux[i][2];
					 Tn[i]++;}
					 
					 if(flux[i][3]!=0){ 
					 fluxusok[i][Tn[i]].l=lsr;
					 fluxusok[i][Tn[i]].flux=flux[i][3];
					 fluxusok[i][Tn[i]].err=flux[i][4];
					 Tn[i]++;}
					 
					 if(flux[i][5]!=0){ 
					 fluxusok[i][Tn[i]].l=lsi;
					 fluxusok[i][Tn[i]].flux=flux[i][5];
					 fluxusok[i][Tn[i]].err=flux[i][6];
					 Tn[i]++;}
					 
					 if(flux[i][7]!=0){ 
					 fluxusok[i][Tn[i]].l=lsz;
					 fluxusok[i][Tn[i]].flux=flux[i][7];
					 fluxusok[i][Tn[i]].err=flux[i][8];
					 Tn[i]++;}
					 
					 if(flux[i][9]!=0){ 
					 fluxusok[i][Tn[i]].l=lsu;
					 fluxusok[i][Tn[i]].flux=flux[i][9];
					 fluxusok[i][Tn[i]].err=flux[i][10];
					 Tn[i]++;}
					 
					 }
					
					

// Swift+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//FILE *fsw;
//fsw=fopen("swiftbe.txt","rb");
std::ifstream fsw("swiftbe.txt"); 
//if(fsw==NULL) printf("\nNincs Swift file!\n\n");
if(!fsw.is_open()) printf("\nNincs Swift file!\n\n");

//if(fsw!=NULL){
if(fsw.is_open()){
   
   
int swszam=0;

   //swift adat tomb
adatswift swbe[100];    // beolvasasi tomb, MAG!!!
adatswift swbeerr[100]; // hibat beolvasva, mag
adatswift sw[100];      // konvertalo, FLUXUS!!!
adatswift swerr[100];  //e felettinek a hibaja, fluxus
adatswift sws[100];  // fluxus
adatswift swse[100]; // fluxus hiba
adatswift swcsere; //cserek
adatswift swcsereb;
adatswift ZP;   // swift zeruspontok
adatswift hh;   // swift hullamhosszak
adatswift null; // null "vektor"
adatswift nullm; // -1 "vektor"
adatswift Sned; // extinkcio

   nullm.JD=0;
   nullm.w2=-1.0;
   nullm.m2=-1.0;
   nullm.w1=-1.0;
   nullm.u =-1.0;
   nullm.b =-1.0;
   nullm.v =-1.0;
   null.JD=0;
   null.w2=0;
   null.m2=0;
   null.w1=0;
   null.u =0;
   null.b =0;
   null.v =0;

for(i=0;i<100;i++){
	swbe[i]=nullm;
	swbeerr[i]=null;
	sw[i]=nullm;
	swerr[i]=null;
	sws[i]=null;
	swse[i]=null;
}
      
   
   
   //zeruspontok (2. az extinkcio, NED)
   ZP.JD=0;
   ZP.w2=20.66;
   ZP.m2=20.85;
   ZP.w1=21.01;
   ZP.u =21.13;
   ZP.b =20.47;
   ZP.v =21.06;
   
   //zeruspontok (A)
   hh.JD=0;
   hh.w2=2030.0;
   hh.m2=2231.0;
   hh.w1=2634.0;
   hh.u =3465.0;
   hh.b =4392.0;
   hh.v =5468.0;
   
   //extinkcio
   Sned.JD=0;
   Sned.w2=7.63*EBV;
   Sned.m2=8.37*EBV;
   Sned.w1=6.18*EBV;
   Sned.u =5.00*EBV;
   Sned.b =4.16*EBV;
   Sned.v =3.16*EBV;  
   
   ZP=ZP-Sned;
   
   //swift adat beolvasas
   //printf("Beolvasott Swift adatok: (v hiba nelkul)\n");
   printf("\nBeolvasott Swift adatok: \n");
   i=0;
   /*
   while (!feof(fsw)){	  
   		 				fscanf(fsw,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n",&swbe[i].JD,&swbe[i].w2,&swbeerr[i].w2,&swbe[i].m2,&swbeerr[i].m2,&swbe[i].w1,&swbeerr[i].w1,&swbe[i].u,&swbeerr[i].u,&swbe[i].b,&swbeerr[i].b,&swbe[i].v,&swbeerr[i].v);  
                        //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n",swbe[i].JD,swbe[i].w2,swbeerr[i].w2,swbe[i].m2,swbeerr[i].m2,swbe[i].w1,swbeerr[i].w1,swbe[i].u,swbeerr[i].u,swbe[i].b,swbeerr[i].b,swbe[i].v,swbeerr[i].v);  
                        cout << swbe[i].JD << " " << swbe[i].w2 << " " << swbeerr[i].w2 << " " << swbe[i].m2 << " " << swbeerr[i].m2 << " " << swbe[i].w1 << " " << swbeerr[i].w1 << " " << swbe[i].u << " " << swbeerr[i].u << " " << swbe[i].b  << " " << swbeerr[i].b << " " << swbe[i].v << " "  << "\n";
                        // << swbeerr[i].v
                        
						swbeerr[i].JD=swbe[i].JD; swszam++; i++; 
						}*/
                              
   //swift-re ugyan az mint az alapra
   
   
   
   		   
while(getline(fsw,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    stringstream strs(sor);

    //strs >> mezo1;
    strs >> swbe[i].JD;
    strs >> swbe[i].w2;
    strs >> swbeerr[i].w2;
    strs >> swbe[i].m2;
    strs >> swbeerr[i].m2;
    strs >> swbe[i].w1;
    strs >> swbeerr[i].w1;
    strs >> swbe[i].u;
    strs >> swbeerr[i].u;
    strs >> swbe[i].b;
    strs >> swbeerr[i].b;
    strs >> swbe[i].v;
    strs >> swbeerr[i].v;
    

    
    
    
    cout << swbe[i].JD << " " << swbe[i].w2 << " " << swbeerr[i].w2 << " " << swbe[i].m2 << " " << swbeerr[i].m2 << " " << swbe[i].w1 << " " << swbeerr[i].w1 << " " << swbe[i].u << " " << swbeerr[i].u << " " << swbe[i].b  << " " << swbeerr[i].b << " " << swbe[i].v << swbeerr[i].v  << " "  << "\n";

    swbeerr[i].JD=swbe[i].JD; swszam++; i++; 
   
}



   
   
   
   printf("\nsw szam=%d\n\n",swszam);
                         
   //fclose(fsw);
   fsw.close();
   
   
   
   
   //JD-hez hozza ad egy szamot
   for(i=0;i<swszam;i++) swbe[i].JD=datenull+swbe[i].JD;
   
   
   // 0 se szamit jonak!!!
   for(i=0;i<swszam;i++){
   		 if(swbe[i].w2==0) swbe[i].w2=-1;
   		 if(swbe[i].m2==0) swbe[i].m2=-1;
   		 if(swbe[i].w1==0) swbe[i].w1=-1;
   		 if(swbe[i].u ==0) swbe[i].u =-1;
   		 if(swbe[i].b ==0) swbe[i].b =-1;
   		 if(swbe[i].v ==0) swbe[i].v =-1;
   		 }


   // sorba rendezi a JD-t a swiftre
   for(i=0;i<swszam;i++){
					 
					 for(j=0;j<swszam-1;j++)					  
                    	if(swbe[j].JD>swbe[j+1].JD){ swcsere=swbe[j]; swbe[j]=swbe[j+1]; swbe[j+1]=swcsere; 
						swcsereb=swbeerr[j]; swbeerr[j]=swbeerr[j+1]; swbeerr[j+1]=swcsereb; }
					
					}



																					
// interpol(float x, float xm, float ym, float xp, float yp)		  
//swift interpolacio azonos JD-re															
for(j=0;j<szumk;j++){
   sw[j].JD=swerr[j].JD=JD[j];
   if(swbe[swszam-1].JD>=JD[j])
                    for(i=0;i<swszam-1;i++){

	
	                //w2
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].w2==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].w2==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       //printf("%d %d\n",csekk,csekkb);
					                       
					                       if( swbe[i-csekk].w2!=-1.0 && swbe[i+1+csekkb].w2!=-1.0 ){
										   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].w2,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].w2,"w2");
					                       sw[j].w2=interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].w2,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].w2);
					                       swerr[j].w2=interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].w2,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].w2);
										   //swerr[j].w2=h(swbeerr[i-csekk].w2,swbeerr[i+1+csekkb].w2); 	
										   																}   
										   }

	                //m2
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].m2==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].m2==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       
					                       if( swbe[i-csekk].m2!=-1.0 && swbe[i+1+csekkb].m2!=-1.0 ){
										   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].m2,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].m2,"m2");
										   sw[j].m2=interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].m2,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].m2);
										   swerr[j].m2=interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].m2,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].m2);					                       
					                       //swerr[j].m2=h(swbeerr[i-csekk].m2,swbeerr[i+1+csekkb].m2); 	
										   																}
										   }
										   
					                    
	                //w1
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].w1==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].w1==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       
					                       if( swbe[i-csekk].w1!=-1.0 && swbe[i+1+csekkb].w1!=-1.0 ){
 	   									   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].w1,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].w1,"w1");
										   sw[j].w1=interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].w1,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].w1);
										   swerr[j].w1=interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].w1,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].w1);					                       
					                       //swerr[j].w1=h(swbeerr[i-csekk].w1,swbeerr[i+1+csekkb].w1); 	
										   																}
										   }
										   
					                    
	                //u
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].u==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].u==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       
					                       if( swbe[i-csekk].u!=-1.0 && swbe[i+1+csekkb].u!=-1.0 ){
 	   									   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].u ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].u ,"S u");
										   sw[j].u =interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].u ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].u );
										   swerr[j].u =interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].u ,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].u );					                       
										   //swerr[j].u =h(swbeerr[i-csekk].u,swbeerr[i+1+csekkb].u);	
										   			  												}    
										   }

	                //b
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].b==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].b==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       
					                       if( swbe[i-csekk].b!=-1.0 && swbe[i+1+csekkb].b!=-1.0 ){
 	   									   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].b ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].b ,"S b");	
										   sw[j].b =interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].b ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].b );	
										   swerr[j].b =interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].b ,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].b );				                       
										   //swerr[j].b =h(swbeerr[i-csekk].b,swbeerr[i+1+csekkb].b);	
										   			  												}    
										   }
					                    
	                //v
					csekk=0; csekkb=0;
					if(swbe[i].JD<=JD[j] && swbe[i+1].JD>JD[j]){  
                                           while(swbe[i-csekk].v==-1.0 && i-csekk>0) csekk++;
					                       while(swbe[i+1+csekkb].v==-1.0 && i+1+csekkb<swszam-1) csekkb++;
					                       
					                       if( swbe[i-csekk].v!=-1.0 && swbe[i+1+csekkb].v!=-1.0 ){
 	   									   interpolszol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].v ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].v ,"S v");	
										   sw[j].v =interpol(JD[j],swbe[i-csekk].JD,swbe[i-csekk].v ,swbe[i+1+csekkb].JD,swbe[i+1+csekkb].v );	
										   swerr[j].v =interpolhib(JD[j],swbe[i-csekk].JD,swbeerr[i-csekk].v ,swbe[i+1+csekkb].JD,swbeerr[i+1+csekkb].v );				                       
										   //swerr[j].v =h(swbeerr[i-csekk].v,swbeerr[i+1+csekkb].v);	
										   			  												}    
										   }
						
						
						
										   
					
					//extra -1 nap, elore
					if(swbe[0].JD- 1 <=JD[j] && swbe[0].JD>JD[j]){  
					                       
					                       if( swbe[0].w2!=-1.0 && swbe[1].w2!=-1.0 ){
					                       sw[j].w2=interpol(JD[j],swbe[0].JD,swbe[0].w2,swbe[1].JD,swbe[1].w2);
					                       swerr[j].w2=interpolhib(JD[j],swbe[0].JD,swbeerr[0].w2,swbe[1].JD,swbeerr[1].w2);
										   //swerr[j].w2=h(swbeerr[0].w2,swbeerr[1].w2); 	
										   												}   
										     
					                       
					                       if( swbe[0].m2!=-1.0 && swbe[1].m2!=-1.0 ){
					                       sw[j].m2=interpol(JD[j],swbe[0].JD,swbe[0].m2,swbe[1].JD,swbe[1].m2);
					                       swerr[j].m2=interpolhib(JD[j],swbe[0].JD,swbeerr[0].m2,swbe[1].JD,swbeerr[1].m2);
										   //swerr[j].m2=h(swbeerr[0].m2,swbeerr[1].m2); 	
										   												}   
										   
					                       
					                       if( swbe[0].w1!=-1.0 && swbe[1].w1!=-1.0 ){
					                       sw[j].w1=interpol(JD[j],swbe[0].JD,swbe[0].w1,swbe[1].JD,swbe[1].w1);
					                       swerr[j].w1=interpolhib(JD[j],swbe[0].JD,swbeerr[0].w1,swbe[1].JD,swbeerr[1].w1);
										   //swerr[j].w1=h(swbeerr[0].w1,swbeerr[1].w1); 	
										   												}   
										    
					                       
					                       if( swbe[0].u!=-1.0 && swbe[1].u!=-1.0 ){
					                       sw[j].u=interpol(JD[j],swbe[0].JD,swbe[0].u,swbe[1].JD,swbe[1].u);
					                       swerr[j].u=interpolhib(JD[j],swbe[0].JD,swbeerr[0].u,swbe[1].JD,swbeerr[1].u);
										   //swerr[j].u=h(swbeerr[0].u,swbeerr[1].u); 	
										   												}   
										   
					                       
					                       if( swbe[0].b!=-1.0 && swbe[1].b!=-1.0 ){
					                       sw[j].b=interpol(JD[j],swbe[0].JD,swbe[0].b,swbe[1].JD,swbe[1].b);
					                       swerr[j].b=interpolhib(JD[j],swbe[0].JD,swbeerr[0].b,swbe[1].JD,swbeerr[1].b);
										   //swerr[j].b=h(swbeerr[0].b,swbeerr[1].b); 	
										   												}   
										   
					                       
					                       if( swbe[0].v!=-1.0 && swbe[1].v!=-1.0 ){
					                       sw[j].v=interpol(JD[j],swbe[0].JD,swbe[0].v,swbe[1].JD,swbe[1].v);
					                       swerr[j].v=interpolhib(JD[j],swbe[0].JD,swbeerr[0].v,swbe[1].JD,swbeerr[1].v);
										   //swerr[j].v=h(swbeerr[0].v,swbeerr[1].v); 	
										   												}   
										   }
					
										   
		            //extra 2 nap
					if(swbe[swszam-1].JD<=JD[j] && swbe[swszam-1].JD+ 2 >JD[j]){  
					                       
					                       if( swbe[swszam-1].w2!=-1.0 && swbe[swszam-2].w2!=-1.0 ){
					                       sw[j].w2=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].w2,swbe[swszam-1].JD,swbe[swszam-1].w2);
					                       swerr[j].w2=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].w2,swbe[swszam-1].JD,swbeerr[swszam-1].w2);
										   //swerr[j].w2=h(swbeerr[swszam-1].w2,swbeerr[swszam-2].w2); 	
										   																}   
										     
					                       
					                       if( swbe[swszam-1].m2!=-1.0 && swbe[swszam-2].m2!=-1.0 ){
					                       sw[j].m2=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].m2,swbe[swszam-1].JD,swbe[swszam-1].m2);
					                       swerr[j].m2=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].m2,swbe[swszam-1].JD,swbeerr[swszam-1].m2);
										   //swerr[j].m2=h(swbeerr[swszam-1].m2,swbeerr[swszam-2].m2); 	
										   																}   
										   
					                       
					                       if( swbe[swszam-1].w1!=-1.0 && swbe[swszam-2].w1!=-1.0 ){
					                       sw[j].w1=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].w1,swbe[swszam-1].JD,swbe[swszam-1].w1);
					                       swerr[j].w1=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].w1,swbe[swszam-1].JD,swbeerr[swszam-1].w1);
										   //swerr[j].w1=h(swbeerr[swszam-1].w1,swbeerr[swszam-2].w1); 	
										   																}   
										    
					                       
					                       if( swbe[swszam-1].u!=-1.0 && swbe[swszam-2].u!=-1.0 ){
					                       sw[j].u=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].u,swbe[swszam-1].JD,swbe[swszam-1].u);
					                       swerr[j].u=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].u,swbe[swszam-1].JD,swbeerr[swszam-1].u);
										   //swerr[j].u=h(swbeerr[swszam-1].u,swbeerr[swszam-2].u); 	
										   															}   
										   
					                       
					                       if( swbe[swszam-1].b!=-1.0 && swbe[swszam-2].b!=-1.0 ){
					                       sw[j].b=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].b,swbe[swszam-1].JD,swbe[swszam-1].b);
					                       swerr[j].b=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].b,swbe[swszam-1].JD,swbeerr[swszam-1].b);
										   //swerr[j].b=h(swbeerr[swszam-1].b,swbeerr[swszam-2].b); 	
										   															}   
										   
					                       
					                       if( swbe[swszam-1].v!=-1.0 && swbe[swszam-2].v!=-1.0 ){
					                       sw[j].v=interpol(JD[j],swbe[swszam-2].JD,swbe[swszam-2].v,swbe[swszam-1].JD,swbe[swszam-1].v);
					                       swerr[j].v=interpolhib(JD[j],swbe[swszam-2].JD,swbeerr[swszam-2].v,swbe[swszam-1].JD,swbeerr[swszam-1].v);
										   //swerr[j].v=h(swbeerr[swszam-1].v,swbeerr[swszam-2].v); 	
										   															}   
										   
										   }
		            
										   

                     }
}   

   
swszam=szumk;

   // fájlba íratás
   //printf("\nswift JD:\n");
   fprintf(f,"\nSwift, After interpolation:\n");
   fprintf(f,"Swift Interpolalas utan:\nJD       w2 mag    w2 err    m2 mag    m2 err    w1 mag    w1 err    u mag    u err    b mag    b err    v mag    v err\n");	
				
   for(i=0;i<swszam;i++){
					fprintf(f,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n",sw[i].JD,sw[i].w2,swerr[i].w2,sw[i].m2,swerr[i].m2,sw[i].w1,swerr[i].w1,sw[i].u,swerr[i].u,sw[i].b,swerr[i].b,sw[i].v,swerr[i].v);
					//printf("%f\n",sw[i].JD);
										
					}
   printf("\n");


   // fluxust szamol [erg/s/cm2/A] (2. a hiba)
   // -1 es erteket nullazza
   for(i=0;i<swszam;i++){
   sws[i] =texp(-0.4*sw[i])*texp(-0.4*ZP);
   swse[i]=texp(-0.4*sw[i])*texp(-0.4*ZP) * 0.4*log(10)*swerr[i];
   }
   
      // fájlba íratás, fluxus
   fprintf(f,"\nSwift fluxus/flux: [erg/cm/cm/s/A]\nJD       w2 mag    w2 err    m2 mag    m2 err    w1 mag    w1 err    u mag    u err    b mag    b err    v mag    v err   \n");				
   for(i=0;i<swszam;i++){
					fprintf(f,"%4.1f  %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e\n",sws[i].JD,sws[i].w2,swse[i].w2,sws[i].m2,swse[i].m2,sws[i].w1,swse[i].w1,sws[i].u,swse[i].u,sws[i].b,swse[i].b,sws[i].v,swse[i].v);
										
					}
					
					
					
					
				// beleiras az adatflux tipusba
				for(i=0;i<szumk;i++){
												 
					if(sws[i].w2!=0){ 
					fluxusok[i][Tn[i]].l=hh.w2;
					fluxusok[i][Tn[i]].flux=sws[i].w2;
					fluxusok[i][Tn[i]].err=swse[i].w2;
					Tn[i]++;}
					 
					if(sws[i].m2!=0){  
					fluxusok[i][Tn[i]].l=hh.m2;
					fluxusok[i][Tn[i]].flux=sws[i].m2;
					fluxusok[i][Tn[i]].err=swse[i].m2;
					Tn[i]++;}
					
					if(sws[i].w1!=0){  
					fluxusok[i][Tn[i]].l=hh.w1;
					fluxusok[i][Tn[i]].flux=sws[i].w1;
					fluxusok[i][Tn[i]].err=swse[i].w1;
					Tn[i]++;}
					
					if(sws[i].u!=0){  
					fluxusok[i][Tn[i]].l=hh.u;
					fluxusok[i][Tn[i]].flux=sws[i].u;
					fluxusok[i][Tn[i]].err=swse[i].u;
					Tn[i]++;}
					
					if(sws[i].b!=0){  
					fluxusok[i][Tn[i]].l=hh.b;
					fluxusok[i][Tn[i]].flux=sws[i].b;
					fluxusok[i][Tn[i]].err=swse[i].b;
					Tn[i]++;}
					
					if(sws[i].v!=0){  
					fluxusok[i][Tn[i]].l=hh.v;
					fluxusok[i][Tn[i]].flux=sws[i].v;
					fluxusok[i][Tn[i]].err=swse[i].v;
					Tn[i]++;}
	
					}
   
}
// Swift+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// JHK ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//FILE *fjhk;
//fjhk=fopen("JHK-gorbe.txt","rb");
std::ifstream fjhk("JHK-gorbe.txt"); 
//if(fjhk==NULL) printf("\nNincs JHK file!\n\n");
if(!fjhk.is_open()) printf("\nNincs JHK file!\n\n");

//if(fjhk!=NULL){
if(fjhk.is_open()){
			   
int szumj;  // JHK szam
float Adatj[m][7];	// JHK adattomb
float Adatjk[m][7]; // JHK adattomb 2
float fluxjhk[m][7]; // JHK fluxusok
//nullazza
for(i=0;i<m;i++){ Adatjk[i][1]=-1; Adatjk[i][3]=-1; Adatjk[i][5]=-1;
				  Adatjk[i][2]=0; Adatjk[i][4]=0; Adatjk[i][6]=0; Adatjk[i][0]=0;
				  }
for(i=0;i<m;i++) for(k=0;k<7;k++) fluxjhk[i][k]=0;
				  


i=0; 
szumj=0;

/*while (!feof(fjhk)){  
			   //fscanf(fall,"%f   %f   %f   %f   %f   %f   %f \n",&Adatj[i][0],&Adatj[i][1],&Adatj[i][2],&Adatj[i][3],&Adatj[i][4],&Adatj[i][5],&Adatj[i][6]);
			   
			   fscanf(fjhk,"%f/%f/%f  %f   %f   %f   %f   %f   %f    %f    %f  \n",&csere,&csere,&csere,&csere,&Adatj[i][0],&Adatj[i][1],&Adatj[i][2],&Adatj[i][3],&Adatj[i][4],&Adatj[i][5],&Adatj[i][6]);
			   
			   szumj++; i++;}*/
			   
			   
			   		   
while(getline(fjhk,sor))
     {

    if(sor.find('#')!=string::npos) continue; //Ha egy sor # kezdodik kihagyjuk, mivel comment vagy fejlec
    
    stringstream strs(sor);

    //strs >> mezo1;
    //strs >> csere;
    //strs >> csere;
    //strs >> csere;
    //strs >> csere;
    strs >> Adatj[i][0];   
    strs >> Adatj[i][1];   
    strs >> Adatj[i][2];   
    strs >> Adatj[i][3];   
    strs >> Adatj[i][4];   
    strs >> Adatj[i][5];
    strs >> Adatj[i][6];


    szumj++; i++;

}



fjhk.close();
			   
	

for(i=0;i<szumj;i++){
					if(Adatj[i][1]==0) Adatj[i][1]=-1; 
					if(Adatj[i][3]==0) Adatj[i][3]=-1;
					if(Adatj[i][5]==0) Adatj[i][5]=-1;
					}

printf("JHK szam:%d\n\n",szumj);


//rendezes
//sorba rendezi a JD-t			 
for(i=0;i<szumj;i++){
					 for(j=0;j<szumj-1;j++)					  
                    	if(Adatj[j][0]>Adatj[j+1][0]) for(k=0;k<7;k++){ cserea[k]=Adatj[j][k]; Adatj[j][k]=Adatj[j+1][k]; Adatj[j+1][k]=cserea[k]; }
					
					 }


// azonos JD-juket beteszi, interpolalas nelkul
if(interpolal==0)
for(i=0;i<szumk;i++){ 
					  Adatjk[i][0]=JD[i];
					  for(j=0;j<szumj;j++) if( int(Adatj[j][0])==int(JD[i]) ) for(k=0;k<7;k++) Adatjk[i][k]=Adatj[j][k];
					  }



																	
// interpol(float x, float xm, float ym, float xp, float yp)		  
//JHK interpolacio azonos JD-re	
if(interpolal>=1)														
for(j=0;j<szumk;j++){
   Adatjk[j][0]=JD[j];
   if(Adatj[szumj-1][0]>=JD[j] && !(53642.5<JD[j] && 53677.40>JD[j]) )
                    for(i=0;i<szumj-1;i++){

	                //J
					csekk=0; csekkb=0;
					if(Adatj[i][0]<=JD[j] && Adatj[i+1][0]>JD[j]){  
                                           while(Adatj[i-csekk][1]==-1.0 && i-csekk>0) csekk++;
					                       while(Adatj[i+1+csekkb][1]==-1.0 && i+1+csekkb<szumj-1) csekkb++;
					                       
					                       if( Adatj[i-csekk][1]!=-1.0 && Adatj[i+1+csekkb][1]!=-1.0 ){
										   interpolszol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][1] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][1] ,"J");
										   Adatjk[j][1] =interpol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][1] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][1] );
										   Adatjk[j][2] =interpolhib(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][2] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][2] );					                       
										   //Adatjk[j][2] =h(Adatj[i-csekk][2],Adatj[i+1+csekkb][2]);	
										   															}    
										   }

	                //H
					csekk=0; csekkb=0;
					if(Adatj[i][0]<=JD[j] && Adatj[i+1][0]>JD[j]){  
                                           while(Adatj[i-csekk][3]==-1.0 && i-csekk>0) csekk++;
					                       while(Adatj[i+1+csekkb][3]==-1.0 && i+1+csekkb<szumj-1) csekkb++;
					                       
					                       if( Adatj[i-csekk][3]!=-1.0 && Adatj[i+1+csekkb][3]!=-1.0 ){
										   interpolszol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][3] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][3] ,"H");
										   Adatjk[j][3] =interpol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][3] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][3] );
										   Adatjk[j][4] =interpolhib(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][4] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][4] );					                       
										   //Adatjk[j][4] =h(Adatj[i-csekk][4],Adatj[i+1+csekkb][4]);	
										   															}    
										   }
					                    
	                //K
					csekk=0; csekkb=0;
					if(Adatj[i][0]<=JD[j] && Adatj[i+1][0]>JD[j]){  
                                           while(Adatj[i-csekk][5]==-1.0 && i-csekk>0) csekk++;
					                       while(Adatj[i+1+csekkb][5]==-1.0 && i+1+csekkb<szumj-1) csekkb++;
					                       
					                       if( Adatj[i-csekk][5]!=-1.0 && Adatj[i+1+csekkb][5]!=-1.0 ){
										   interpolszol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][5] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][5] ,"K");
										   Adatjk[j][5] =interpol(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][5] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][5] );					                       
										   Adatjk[j][6] =interpolhib(JD[j],Adatj[i-csekk][0],Adatj[i-csekk][6] ,Adatj[i+1+csekkb][0],Adatj[i+1+csekkb][6] );
										   //Adatjk[j][6] =h(Adatj[i-csekk][6],Adatj[i+1+csekkb][6]);	
										   															}    
										   }
										   
		   }
}




fprintf(f,"\nInterp.:JD       J                  H                    K\n");
for(i=0;i<szumk;i++) fprintf(f,"%f  %f  %f  %f  %f  %f  %f\n",Adatjk[i][0],Adatjk[i][1],Adatjk[i][2],Adatjk[i][3],Adatjk[i][4],Adatjk[i][5],Adatjk[i][6]);



// J:31.47 	 H:11.38 	K:3.961   x 1e-11 erg/cm/cm/s/A
for(i=0;i<szumk;i++){

					fluxjhk[i][0]=Adatjk[i][0];
					
					
					if(Adatjk[i][1]!=-1.0){
					Adatjk[i][1]=Adatjk[i][1]-nedj;
					fluxjhk[i][1]= pow(10.0,-0.4*(Adatjk[i][1]+23.76));	 
					fluxjhk[i][2]= pow(10.0,-0.4*(Adatjk[i][1]+23.76)) * 0.4*log(10)*Adatjk[i][2];
					}
					
					
					if(Adatjk[i][3]!=-1.0){
					Adatjk[i][3]=Adatjk[i][3]-nedh;
					fluxjhk[i][3]= pow(10.0,-0.4*(Adatjk[i][3]+24.86));	
					fluxjhk[i][4]= pow(10.0,-0.4*(Adatjk[i][3]+24.86)) * 0.4*log(10)*Adatjk[i][4];
					}
					
					
					if(Adatjk[i][5]!=-1.0){
					Adatjk[i][5]=Adatjk[i][5]-nedk;
					fluxjhk[i][5]= pow(10.0,-0.4*(Adatjk[i][5]+26.00));	 
					fluxjhk[i][6]= pow(10.0,-0.4*(Adatjk[i][5]+26.00)) * 0.4*log(10)*Adatjk[i][6];
					}
															
             
			 }



//flux kiiras
fprintf(f,"\nFluxusok/Fluxes: JD  J H K es hibaik (errors) [erg/cm/cm/s/A]\n");

for(i=0;i<szumk;i++){
					fprintf(f,"%4.3f   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e   %4.3e \n",fluxjhk[i][0],fluxjhk[i][1],fluxjhk[i][2],fluxjhk[i][3],fluxjhk[i][4],fluxjhk[i][5],fluxjhk[i][6]);
										
					
					}
					

			

// adatflux tipusba taroljuk tovabb, + JD-ben az idot
// JD es Tn is itt kezdodik
//lsx formatum van, de lsx=ljx, igy van megadva
for(i=0;i<szumk;i++){
					 
					 if(fluxjhk[i][1]!=0){ 
					 fluxusok[i][Tn[i]].l=12600;
					 fluxusok[i][Tn[i]].flux=fluxjhk[i][1];
					 fluxusok[i][Tn[i]].err=fluxjhk[i][2];
					 Tn[i]++;}
					 
					 if(fluxjhk[i][3]!=0){ 
					 fluxusok[i][Tn[i]].l=16000;
					 fluxusok[i][Tn[i]].flux=fluxjhk[i][3];
					 fluxusok[i][Tn[i]].err=fluxjhk[i][4];
					 Tn[i]++;}
					 
					 if(fluxjhk[i][5]!=0){ 
					 fluxusok[i][Tn[i]].l=22200;
					 fluxusok[i][Tn[i]].flux=fluxjhk[i][5];
					 fluxusok[i][Tn[i]].err=fluxjhk[i][6];
					 Tn[i]++;}
					 
					 }





}
// JHK +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					



//rendezes
//fluxusok[szumk][20];
// sorba rendezi a hullamhosszat-t
for(k=0;k<szumk;k++)
					 
            for(i=0;i<Tn[k];i++){
					 for(j=0;j<Tn[k]-1;j++)					  
                    	if(fluxusok[k][j].l>fluxusok[k][j+1].l){ cseref=fluxusok[k][j]; fluxusok[k][j]=fluxusok[k][j+1]; fluxusok[k][j+1]=cseref; }
					
					 }

printf("\nJD, letezo szurok szama, es hullamhosszaik (A):\n");
printf("JD, number of valid bands, and their wavelength (A):\n");
for(i=0;i<szumk;i++){
					 printf("%4.2f Jo:%d   ",JD[i],Tn[i]);
					 for(j=0;j<Tn[i];j++) cout << fluxusok[i][j].l << " ";
					 printf("\n");
					 }
printf("\n");


fprintf(f,"\nJD, letezo szurok szama, es hullamhosszaik (A):\n");
fprintf(f,"JD, number of valid bands, and their wavelength (A):\n");
for(i=0;i<szumk;i++){
					 fprintf(f,"%4.2f Jo:%d    ",JD[i],Tn[i]);
					 for(j=0;j<Tn[i];j++) fprintf(f,"%d ",int(fluxusok[i][j].l));
					 fprintf(f,"\n");
					 }
fprintf(f,"\n");




csere=0;
for(int i=0;i<9;i++) csere=csere+szabadlt[i];
if(csere==0) cout << "!!! Nincs ervenyes hullamhossz az illeszteshez / No valid wavelength for BB fit !!!" << endl;



// fluxus integralas +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// bolometrikus fényeséggek
double bolom[szumk][5];


// bolometrikus fényesség elkészitése, integrálása, PROBLEMA: ahol nincs adat, ott 0
for(i=0;i<szumk;i++){					 
					 
					bolom[i][0]=flux[i][0];
					 					 
					//if(mod==0) bolom[i][1]= bolints(flux[i][1],flux[i][3],flux[i][5],flux[i][7]);
					//if(mod==1) bolom[i][1]= bolintj(flux[i][1],flux[i][3],flux[i][5],flux[i][7]);
					if(Tn[i]>2) fprintf(Trki,"%f ",bolom[i][0]);
					bolom[i][1]= bolint(fluxusok[i],Tn[i]);
					
			        //if(mod==0) bolom[i][2]= bolhibs(flux[i][1],flux[i][3],flux[i][5],flux[i][7],flux[i][2],flux[i][4],flux[i][6],flux[i][8]);
					//if(mod==1) bolom[i][2]= bolhibj(flux[i][1],flux[i][3],flux[i][5],flux[i][7],flux[i][2],flux[i][4],flux[i][6],flux[i][8]);
					bolom[i][2]= bolhib(fluxusok[i],Tn[i]);
								        	        
			        bolom[i][3]=-2.5*log10(4*PI*tav*tav*bolom[i][1]/1000.0)+71.21;
			        			        
			        bolom[i][4]=2.5*(bolom[i][2]/bolom[i][1])/(log(10));
					 
					}
					 


// bolometrikus fényességek kiiratása
fprintf(bol,"#JD \t bol flux \t err \t bol mag \t err \n");
for(i=0;i<szumk;i++){
					 if(bolom[i][1]!=0) fprintf(bol,"%f %e %e %f %f\n",bolom[i][0],bolom[i][1],bolom[i][2],bolom[i][3],bolom[i][4]);
					 
					 }
					 
					 
// CSAK a bolometrikus magnik -> kompatibilis az illesztovel
fprintf(bol2,"#JD \t bol mag \t err \n");
for(i=0;i<szumk;i++){
					 if(bolom[i][1]!=0) fprintf(bol2,"%f %f %f\n",bolom[i][0],bolom[i][3],bolom[i][4]);
					 
					 }
					 

// BC=bol-g    szin=g-r
// BC=bol-B    szin=B-I
// BC es szint is kiirja pluszba meg				 
csekk=0;
if(mod==0) fprintf(bc,"#JD \t bol mag \t err \t BC=bol-g \t err \t g-r \t err \n");
if(mod==1) fprintf(bc,"#JD \t bol mag \t err \t BC=bol-B \t err \t B-I \t err \n"); 
for(i=0;i<szumk;i++){
					 while( 1-tiltas(Adat[i+csekk][0]) ) csekk++;
					 //if ( tiltas(Adat[i+csekk][0]) ) ;  else csekk++;
					 if(mod==0) if(Adat[i][1]!=-1 && Adat[i][3]!=-1) fprintf(bc,"%f %f %f %f %f %f %f\n",bolom[i][0],bolom[i][3],bolom[i][4],bolom[i][3]-Adat[i+csekk][1]+5*log10(tav /3.086/pow(10,16))-5,h(bolom[i][4],Adat[i+csekk][2]),Adat[i+csekk][1]-Adat[i+csekk][3],h(Adat[i+csekk][2],Adat[i+csekk][4]));
					 if(mod==1) if(Adat[i][7]!=-1 && Adat[i][5]!=-1) fprintf(bc,"%f %f %f %f %f %f %f\n",bolom[i][0],bolom[i][3],bolom[i][4],bolom[i][3]-Adat[i+csekk][7]+5*log10(tav /3.086/pow(10,16))-5,h(bolom[i][4],Adat[i+csekk][8]),Adat[i+csekk][7]-Adat[i+csekk][5],h(Adat[i+csekk][8],Adat[i+csekk][6]));
					 
					 }


// Kimenet: MJD, fluxus [erg/s/cm2], fluxus hiba [erg/s/cm2], abszolút bolometrikus fényesség [watt, magni], hiba [-II-]

fclose(f);
fclose(bol);
fclose(bol2);
fclose(bc);
fclose(Trki);

printf("\nLefutott\n");
		
//getchar();
	
	
	
	
}




