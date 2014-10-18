//Nygren-Keener_model 0d
//Carlos-Lugo UPC 2010
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
//Domain number of nodes and variables// 
#define Lx 300    //Nodes in x
#define Ly 150    //Nodes in y
#define lx 2   //Initial conf Lx
#define ly 1   //Initial conf Ly
//Initial Conditions
const double vo=-74.2525; //mV
const double cnaco=130.0; //130.0110; //mmol/L
const double ckco=5.4;    //5.3581; //mmol/L
const double ccaco=1.8;   //1.8147; //mmol/L
const double cnaio=8.5547; //mmol/L
const double ckio=129.4350; //mmol/L
const double ccaso=6.72495e-5; //mmol/L FORMER ccaio
const double ccado=7.2494e-5; //mmol/L
const double ccaupo=0.664564;//0.6646; //mmol/L
const double ccarelo=0.646458;//0.6465; //mmol/L
const double mo=3.2017e-3; 
const double h1o=0.881421; //0.8814;
const double h2o=0.874233; //0.8742;
const double dlo=1.3005e-5;
const double fl1o=0.998639; //0.9986;
const double fl2o=0.998604; //0.9986;
const double ro=1.0678e-3;
const double so=0.948975; //0.949;
const double rsuso=1.5949e-4;
const double ssuso=0.991169;//0.9912;
const double no=4.8357e-3;
const double pao=5.16114e-5; //0.0001;
const double F1o=0.428396;//0.4284;
const double F2o=0.0028005;  //0.0028;
const double OCo=0.0274971;  //0.0275;
const double OTCo=0.0132801; //0.0133;
const double OTMgCo= 0.196085; //0.1061;
const double OTMgMgo=0.709417; //0.7094;
const double Ocsqo=0.43686; //0.4369;
const double x00o=0.998;
const double x10o=0.0;
const double x11o=0.0;
const double x01o=0.002;
const double ccaio=6.72495e-5; //New Caio
//////////////////////////
//PARAMETERS
/////////////////////////
const double cnab=130.0; //mmol/L
const double ckb=5.4;//mmol/L
const double ccab=1.8;//mmol/L
const double cmgi=2.5;//mmol/L
const double ecaapp=60.0;//mV
const double kca=0.025; //mmol/L
const double R=8314.0; //mJ/molK
const double T=306.15; //K
const double F=96487; //C/mol
const double Cm=1.0;  //¿pF?,¿nF?
const double voli=0.11768; //nL
const double volrel=0.000882;//nL
//const double volup=0.007938;//nL
const double tauna=14300;//msec
const double tauk=10000;//msec
const double tauca=24700;//msec
const double taudi=1.0; //msec
const double inakb=1.416506;//¿pA?
const double knakk=1.0; //mmol/L
const double knakna=11.0;//mmol/L
const double icapb=0.08; //¿pA?
const double kcap=0.0002;//mmol/L
const double knaca=0.000749684; //0.00079684;//pA/(mmol/L)^4?
const double gama=0.45;
const double dnaca=0.0003;//(mmol/L)^-4
const double phinaen=-0.0336;//pA
const double iupb=56;//pA
const double kcyca=0.0003; //mmol/L
const double ksrca=0.5;//mmol/L
const double kxcs=0.4;//
const double tautr=10.0; //msec
//%const double alpharel=4000;//pA*L/mmol
const double arel=0.3;
const double arel2=0.3;
const double kreli=0.0003;//mmol/L
const double kreld=0.003;//mmol/L
const double rrecov=0.815e-3;
/******************************************/
const double km2=0.06; //ms^-1
const double k2=45.0; //mmol^-2msec^-1
const double km6=1.0/200; //msec^-1
const double k6=0.3; //mmol^-1 msec^-1
/******************************************/
const double k7=0.5; //mM^-1 msec^-1
const double Kd3=10.0;//mmol
const double Kd8=0.6383;//mmol
const double pna=3.2e-5;
const double gcal1=0.135; //Nygren
const double gcal2=0.135;  //Nygren-Keener Alternans
const double gt=0.15;
const double gsus=0.055;
const double gks=0.02;
const double gkr=0.01;
const double gk1=0.06;
const double gbna=0.00121198;
const double gbca=0.00157362;
const double CSQtot=1.5; //Keener total CSQ mmol
////////////////////////////////////////
//NEW DIFUSSIVE TIME SARCOLEMMA-CYTOSOL/
const double tausi=0.1;
////////////////////////
const double dt=0.01;
const double dx=0.025;
const double diff=0.0005;
//const double diff=0.00;
////////////////////////
const double tf=10000;//50e3; //Final Time
const double Ts=400; //Stimulus Period msec
const double tso=0.0; //Stimulus Initial time msec
const double tsf=51e3;//Stimulus Final time
const double dts=3.0;//6.0;//stimulus duration msec
const double Amp=20.0; //Amplitude

////////////////////////
const double ld1=1.0; //
const double ld2=0.50; //
const double ld3=1.3; //
const double ld4=70.0; //
const double ld5=1.0/100; 
const double ld6=500.0;
////////////////////////
gsl_rng * rnum;
const gsl_rng_type * Trnum;
////////////////////////
double volup,volc,vold,Rtf,t,vols,volrel2;
double k2d,k2i,km6d,km6i,k6d,k6i,km2i,km2d;
double v[Lx][Ly],vt[Lx][Ly],Nac[Lx][Ly],Kc[Lx][Ly],Cac[Lx][Ly],Nai[Lx][Ly],Ki[Lx][Ly],Cai[Lx][Ly];
double Cad[Lx][Ly],Caup[Lx][Ly],Carel[Lx][Ly],m[Lx][Ly],h1[Lx][Ly],h2[Lx][Ly],DL[Lx][Ly],fL1[Lx][Ly];
double fL2[Lx][Ly],r[Lx][Ly],s[Lx][Ly],rsus[Lx][Ly],ssus[Lx][Ly],n[Lx][Ly],Pa[Lx][Ly],F1[Lx][Ly],F2[Lx][Ly];
double OC[Lx][Ly],OTC[Lx][Ly],OTMgC[Lx][Ly],OTMgMg[Lx][Ly],OCSQ[Lx][Ly],X00[Lx][Ly],X10[Lx][Ly],X11[Lx][Ly],X01[Lx][Ly];
double CT[Lx][Ly];
double Y00[Lx][Ly],Y10[Lx][Ly],Y11[Lx][Ly],Y01[Lx][Ly],Carel2[Lx][Ly],OCSQ2[Lx][Ly];
///NEW Cas THE OLD Cai is now Cas
double Cas[Lx][Ly]; 
////////////////////////
double istim,ena,ek,eca;
double mbar,hbar,auxna1,auxna2,auxna3,taum,tauh1,tauh2,ina;
double  DLbar,fLbar,tauDL,taufL1,taufL2,auxca1,auxca2,auxca3,fCa,ical;
double rbar,sbar,auxr1,auxs1,taur,taus,rsusbar,ssusbar,taursus,taussus,it,isus; 
double nbar,auxks1,taun,iks,pabar,pi,auxkr1,taupa,ikr;
double ik1,ibna,ibca,inak,icap,inacaux1,inacaux2,inacaux3,inaca;
double  DOc,DOTc,DOTMgC,DOTMgMg;
double Iup,Itr,DOcsq,rax1,rax2,ract,rinact,dF1,dF2,rax3,irel;
double q,cadi,act1,act2,inact1,inact2,dx00,dx10,dx01,dx11;
double q2,cadi2,act12,act22,inact12,inact22,dy00,dy10,dy01,dy11,irel2,Itr2,DOcsq2;
double idi,DO,Dv;
double isi; //NEW isi CURRENT
///////////////////////////////////
void inicond(void);//Initial 0d State Setup
void inicond2(void);//Initial 2d State Setup
void prtstuff(double X[Lx][Ly],double tp, int k);//Printout things
void prtstuff2(double X[Lx][Ly],double tp, int k);//Print system snapshots
void Eulerupdate(int eudt, int rnc);
//////////////////////////////////
void Eulerupdate(int eudt, int rnc){
  int i,j,nx,ny;
  double gnx,gnx1,gnx2; 
  if(eudt==0){nx=lx;ny=ly;}
  if(eudt==1){nx=Lx;ny=Ly;}
  if(rnc==1){gnx1=0.6;
             gnx2=0.6;}
  else{gnx1=0.6;
       gnx2=0.6;}  
  
  //x and y loops start
  for(i=0;i<nx;i++){
   for(j=0;j<ny;j++){
     /*printf("************************\n"); 
     printf("Node= %d, %d\n",i,j);
     printf("Vini=%20.19e\n",v[i][j]);
     printf("Time:%lf\n",t);
     printf("Nac:%lf, Nai:%lf, Kc:%lf, Ki:%lf, Cac:%lf, Cai:%lf\n",Nac[i][j],Nai[i][j],Kc[i][j],Ki[i][j],Cac[i][j],Cai[i][j]);*/
     //getchar(); 

     /***********/
     /*Calc-Stim*/
     /***********/
       if((fmod(t,Ts)<dts)&&((i<10)&&(j<Ly))){istim=Amp;}
      // if((fmod(t,Ts)<dts)&&(i<5)){istim=Amp;}
       else{istim=0.0;}
       //printf("t:%10.9e istim:%10.9e\n",t,istim);
     /***************/
     /*E_Na,E_Ca,E_K*/
     /***************/
       ena=Rtf*gsl_sf_log(Nac[i][j]/Nai[i][j]);
       eca=(0.5*Rtf)*gsl_sf_log(Cac[i][j]/Cas[i][j]);
       ek=Rtf*gsl_sf_log(Kc[i][j]/Ki[i][j]); 
       //printf("Ena=%20.19e\nEka=%20.19e\nEca=%20.19e\n",ena,ek,eca);
       //printf("Istim=%lf\n",istim);
     /*************/
     /*Na+ Current*/
     /*************/
       mbar=1.0/(1.0+gsl_sf_exp((v[i][j]+27.12)/(-8.21)));
       //printf("mbar:%10.9e\n",mbar);
       hbar=1.0/(1.0+gsl_sf_exp((v[i][j]+63.6)/(5.3)));
       //printf("hbar:%10.9e\n",hbar);
       ///getchar();
       auxna1=-gsl_sf_pow_int(((v[i][j]+25.57)/28.8),2);
       taum=0.024+gsl_sf_exp_mult(auxna1,0.042);
       tauh1=0.3+30.0/(1.0+gsl_sf_exp((v[i][j]+35.1)/3.2));
       tauh2=3.0+120/(1.0+gsl_sf_exp((v[i][j]+35.1)/3.2));
       //getchar();
       m[i][j]+=dt*(mbar-m[i][j])/taum;
       h1[i][j]+=dt*(hbar-h1[i][j])/tauh1;
       h2[i][j]+=dt*(hbar-h2[i][j])/tauh2;
       auxna2=(v[i][j]-ena)/Rtf;
       auxna3=v[i][j]/Rtf;
       ina=(pna*gsl_sf_pow_int(m[i][j],3))*(F/Rtf)*(v[i][j]*Nac[i][j]);
       ina*=(0.9*h1[i][j]+0.1*h2[i][j])*((gsl_sf_expm1(auxna2))/(gsl_sf_expm1(auxna3)));
       /*printf("m=%10.9e\n",m[i][j]);
       printf("h1=%10.9e\n",h1[i][j]);
       printf("h2=%10.9e\n",h2[i][j]);
       printf("ina=%10.9e\n",ina);*/
       //getchar();
     /*******************/
     /*L-type Ca current*/
     /*******************/
       DLbar=1.0/(1.0+gsl_sf_exp((v[i][j]+9.0)/(-5.8)));
       fLbar=1.0/(1.0+gsl_sf_exp((v[i][j]+27.4)/7.1));
       auxca1=-gsl_sf_pow_int(((v[i][j]+35.0)/30.0),2);
       auxca2=-gsl_sf_pow_int(((v[i][j]+40.0)/14.4),2);
       auxca3=-gsl_sf_pow_int(((v[i][j]+40.0)/14.2),2);
       tauDL=2.0+gsl_sf_exp_mult(auxca1,2.7);
       taufL1=10.0+gsl_sf_exp_mult(auxca2,161.0);
       taufL2=62.6+gsl_sf_exp_mult(auxca3,1332.3);
       fCa=Cad[i][j]/(kca+Cad[i][j]);

       DL[i][j]+=dt*(DLbar-DL[i][j])/tauDL;
       fL1[i][j]+=dt*(fLbar-fL1[i][j])/taufL1;    
       fL2[i][j]+=dt*(fLbar-fL2[i][j])/taufL2;
       ical=DL[i][j]*(v[i][j]-ecaapp)*(fCa*fL1[i][j]+(1.0-fCa)*fL2[i][j]);
       if(CT[i][j]==0.0){ical*=gcal1;}
       else{
       if(CT[i][j]==1.0){gnx=gnx1;}
       if(CT[i][j]==2.0){gnx=gnx2;}
       ical*=(gcal2*gnx);}
       //printf("dl=%10.9E\n",DL[i][j]);
       //printf("fl1=%10.9e\n",fL1[i][j]);
       //printf("fl2=%10.9e\n",fL2[i][j]);
       //printf("ical=%10.9e\n",ical);
     /*****************************************************/
     /*Transient and sustained outward K+ currents It,Isus*/
     /*****************************************************/
       rbar=1.0/(1.0+gsl_sf_exp((v[i][j]-1.0)/(-11.0)));
       sbar=1.0/(1.0+gsl_sf_exp((v[i][j]+40.5)/11.5));
       auxr1=-gsl_sf_pow_int((v[i][j]/30.0),2);
       taur=1.5+gsl_sf_exp_mult(auxr1,3.5);
       auxs1=-gsl_sf_pow_int(((v[i][j]+52.45)/14.97),2);
       taus=14.14+gsl_sf_exp_mult(auxs1,481.2);
       r[i][j]+=dt*(rbar-r[i][j])/taur;
       s[i][j]+=dt*(sbar-s[i][j])/taus;  
       it=gt*r[i][j]*s[i][j]*(v[i][j]-ek);
       
       rsusbar=1.0/(1.0+gsl_sf_exp((v[i][j]+4.3)/(-8.0)));
       ssusbar=0.6+(0.4/(1.0+gsl_sf_exp((v[i][j]+20.0)/10.0)));
       taursus=0.5+(9.0/(1.0+gsl_sf_exp((v[i][j]+5.0)/12.0)));
       taussus=300.0+(47.0/(1.0+gsl_sf_exp((v[i][j]+60.0)/10.0)));
       rsus[i][j]+=dt*(rsusbar-rsus[i][j])/taursus;
       ssus[i][j]+=dt*(ssusbar-ssus[i][j])/taussus;
       isus=gsus*rsus[i][j]*ssus[i][j]*(v[i][j]-ek);
       /*printf("r=%10.9E\n",r[i][j]);
       printf("s=%10.9E\n",s[i][j]);
       printf("rsus=%10.9E\n",rsus[i][j]);
       printf("ssus=%10.9E\n",ssus[i][j]);
       printf("it=%10.9e\n",it);
       printf("isus=%10.9e\n",isus);*/
     /*****************************************/
     /*Delayed Rectifier K+ currents: Iks, Ikr*/
     /*****************************************/
       nbar=1.0/(1.0+gsl_sf_exp((v[i][j]-19.9)/(-12.7)));
       auxks1=-gsl_sf_pow_int(((v[i][j]-20.0)/20.0),2);
       taun=700+gsl_sf_exp_mult(auxks1,400); 
       n[i][j]+=dt*(nbar-n[i][j])/taun;
       iks=gks*n[i][j]*(v[i][j]-ek);
       
       pabar=1.0/(1.0+gsl_sf_exp((v[i][j]+15.0)/(-6.0)));
       pi=1.0/(1.0+gsl_sf_exp((v[i][j]+55.0)/24.0));
       auxkr1=-gsl_sf_pow_int(((v[i][j]+20.1376)/22.1996),2);
       taupa=31.18+gsl_sf_exp_mult(auxkr1,217.18);
       Pa[i][j]+=dt*(pabar-Pa[i][j])/taupa;
       ikr=gkr*Pa[i][j]*pi*(v[i][j]-ek);
       //printf("pa=%10.9E\n",Pa[i][j]);
       //printf("n=%10.9E\n",n[i][j]);
       //printf("iks=%10.9e\n",iks);
       //printf("ikr=%10.9e\n",ikr);     
     /**********************************/
     /*Inward Rectifier K+ current:I_K1*/   
     /**********************************/
       ik1=gk1*(v[i][j]-ek)/(1.0+gsl_sf_exp((1.5/Rtf)*(v[i][j]-ek+3.6)));
       ik1*=pow(Kc[i][j],0.4457);
       //printf("ik1=%10.9e\n",ik1);
     /********************************************/
     /*Background Inward Currents: I_B,Na, I_B,Ca*/
     /********************************************/
       ibna=gbna*(v[i][j]-ena);
       ibca=(gbca*(v[i][j]-eca));
       //printf("ibna=%10.9e\n",ibna);
       //printf("ibca=%10.9e\n",ibca);
     /************************************************/
     /*Pump and exchanger Currents I_Nak,I_CaP,I_NaCa*/
     /************************************************/
       inak=inakb*(Kc[i][j]/(Kc[i][j]+knakk));
       inak*=(sqrt(gsl_sf_pow_int(Nai[i][j],3))/(sqrt(gsl_sf_pow_int(knakna,3))+sqrt(gsl_sf_pow_int(Nai[i][j],3))));
       inak*=((v[i][j]+150.0)/(v[i][j]+200.0));
       //icap=icapb*(Cai[i][j]/(Cai[i][j]+kcap));
       icap=(icapb)*(Cas[i][j]/(Cas[i][j]+kcap));
       inacaux1=knaca*gsl_sf_pow_int(Nai[i][j],3)*Cac[i][j]*gsl_sf_exp((gama*v[i][j]/Rtf));
       //inacaux2=knaca*gsl_sf_pow_int(Nac[i][j],3)*Cai[i][j]*gsl_sf_exp(((gama-1.0)*v[i][j])/Rtf);
       inacaux2=knaca*gsl_sf_pow_int(Nac[i][j],3)*Cas[i][j]*gsl_sf_exp(((gama-1.0)*v[i][j])/Rtf);
       //inacaux3=1.0+dnaca*(gsl_sf_pow_int(Nac[i][j],3)*Cai[i][j]+gsl_sf_pow_int(Nai[i][j],3)*Cac[i][j]);
       inacaux3=1.0+dnaca*(gsl_sf_pow_int(Nac[i][j],3)*Cas[i][j]+gsl_sf_pow_int(Nai[i][j],3)*Cac[i][j]);
       //inaca=(inacaux1-inacaux2)/inacaux3;
       inaca=(inacaux1-inacaux2)/(1.0*inacaux3);

       //inaca=2.0*(inacaux1-inacaux2)/inacaux3;
       //printf("icap=%10.9e\n",icap);
       //printf("inak=%10.9e\n",inak);
       //printf("inaca=%10.9e\n",inaca);
     /****************************/
     /*Intracellular Ca buffering*/
     /****************************/
       DOc=200*Cai[i][j]*(1.0-OC[i][j])-0.476*OC[i][j];
       DOTc=78.4*Cai[i][j]*(1.0-OTC[i][j])-0.392*OTC[i][j];
       DOTMgC=200*Cai[i][j]*(1.0-OTMgC[i][j]-OTMgMg[i][j])-0.0066*OTMgC[i][j];
       DOTMgMg=2.0*cmgi*(1.0-OTMgC[i][j]-OTMgMg[i][j])-0.666*OTMgMg[i][j];

       OC[i][j]+=dt*DOc;
       OTC[i][j]+=dt*DOTc;
       OTMgC[i][j]+=dt*DOTMgC;
       OTMgMg[i][j]+=dt*DOTMgMg;

       //printf("OC:%10.9E\n",OC[i][j]); 
       //printf("OTC:%10.9E\n",OTC[i][j]); 
       //printf("OTMgC:%10.9E\n",OTMgC[i][j]); 
       //printf("OTMgMg:%10.9E\n",OTMgMg[i][j]); 
     /*******************************************/
     /*Cleft Space Ion concentrations:Nac,Kc,Cac*/
     /*******************************************/
       Nac[i][j]+=dt*(((cnab-Nac[i][j])/tauna)+(ina+ibna+3.0*inak+3.0*inaca+(0.0*phinaen)-istim)/(volc*F));
       Kc[i][j]+=dt*(((ckb-Kc[i][j])/tauk)+(it+isus+ik1+iks+ikr-2.0*inak)/(volc*F));
       Cac[i][j]+=dt*(((ccab-Cac[i][j])/tauca)+(ical+ibca+icap-2.0*inaca)/(2.0*volc*F));
       //printf("Nac=%10.9e\nKc=%10.9e\nCac=%10.9e\n",Nac[i][j],Kc[i][j],Cac[i][j]);
     /***********************/
     /*Ca Handling by the SR*/
     /***********************/
      Iup=iupb*(Cai[i][j]/kcyca-(kxcs*kxcs)*(Caup[i][j]/ksrca));
      //Iup=iupb*(Cas[i][j]/kcyca-(kxcs*kxcs)*(Caup[i][j]/ksrca));
      Iup/=(((Cai[i][j]+kcyca)/kcyca)+((kxcs/ksrca)*(Caup[i][j]+ksrca)));
      //Iup/=(((Cas[i][j]+kcyca)/kcyca)+((kxcs/ksrca)*(Caup[i][j]+ksrca)));
      Itr=(Caup[i][j]-Carel[i][j])*(2.0*F*volrel/(tautr));
      Itr2=(Caup[i][j]-Carel2[i][j])*(2.0*F*volrel2/tautr);
      DOcsq=0.48*(Carel[i][j])*(1.0-OCSQ[i][j])-0.4*OCSQ[i][j];
      DOcsq2=0.48*(Carel2[i][j])*(1.0-OCSQ2[i][j])-0.4*OCSQ2[i][j];
      OCSQ[i][j]+=dt*(DOcsq);
      OCSQ2[i][j]+=dt*(DOcsq2);
      //printf("%15.10e %15.10e\n",DOcsq,DOcsq2);
      //getchar();
     ///STANDARD NYGREN
     /*  if(CT[i][j]==0.0){
        rax1=Cai[i][j]/(Cai[i][j]+kreli);
        rax2=Cad[i][j]/(Cad[i][j]+kreld);
        ract=0.2038*(gsl_sf_pow_int(rax1,4)+gsl_sf_pow_int(rax2,4));
        rinact=0.03396+0.3396*gsl_sf_pow_int(rax1,4);
        dF1=rrecov*(1-F1[i][j]-F2[i][j])-ract*F1[i][j];
        dF2=ract*F1[i][j]-rinact*F2[i][j];
        F1[i][j]+=dt*(dF1); 
        F2[i][j]+=dt*(dF2);
        rax3=F2[i][j]/(F2[i][j]+0.25);
        irel=alpharel*gsl_sf_pow_int(rax3,2)*(Carel[i][j]-Cai[i][j]); 
       }*/
      //////////////////
      //////////////////////////////
      ///NEW COMPARTMENTS Cai DYNAMICS   
        if(CT[i][j]>=1.0){

         /**************************************/
         if(CT[i][j]==1.0){
		           /*Dyadic*/
		           km6d=(1.0/200); //1/200
		           k2d=k2;  //45.0
			   k6d=k6;    //0.3
			   km2d=2.0*km2;  //0.06
                          /*Cytosol*/
                           km6i=(1.0/200); //1/200
 			   k6i=k6;   //144-0
			   km2i=2.0*km2;  //0.06          
			   k2i=k2; //9000.0
			   
                         }
         else{
	         /*Dyadic*/
       		 km6d=(1.0/(200+ld6));		          		
       		 //km6d=(1.0/150);		          		
		 k2d=k2;  
		 k6d=k6;
		 km2d=2.0*km2;
		 /*Cytosol*/
		 km6i=(1.0/(200+ld6));                                     
		 //km6i=(1.0/150);                                     
                 k6i=k6;
                 km2i=2.0*km2;           
		 k2i=k2; 	 
             }
        //q=CSQtot*(1.0-OCSQ[i][j]);
        q=0.0;
        cadi=gsl_sf_pow_int(Cad[i][j],2);
        act1=k2d*(ld4*ld4);//(k2d*Kd1+k11*q)/(Kd1+q);        
        act2=k2d*(ld4*ld4);// (k2d*Kd8+k11*q)/(Kd8+q);
        inact1=k6d*ld4;//(k6d*Kd3+k7*q)/(Kd3+q);
        inact2=k6d*ld4;// (k6d*Kd1+k7*q)/(Kd1+q);
        dx00=km2d*X10[i][j]-act1*cadi*X00[i][j]+km6d*X01[i][j]-inact2*Cad[i][j]*X00[i][j];
        dx10=act1*cadi*X00[i][j]-km2d*X10[i][j]+km6d*X11[i][j]-inact1*Cad[i][j]*X10[i][j];
        dx11=inact1*Cad[i][j]*X10[i][j]-km6d*X11[i][j] + act2*cadi*X01[i][j]-km2d*X11[i][j];
        dx01=km2d*X11[i][j]-act2*cadi*X01[i][j]+inact2*Cad[i][j]*X00[i][j]-km6d*X01[i][j];
        X00[i][j]+=dt*(dx00);
        X10[i][j]+=dt*(dx10);
        X11[i][j]+=dt*(dx11);
        X01[i][j]+=dt*(dx01);
        
       /*printf("x00=%10.9E\n",X00[i][j]); 
       printf("x10=%10.9E\n",X10[i][j]); 
       printf("x11=%10.9E\n",X11[i][j]); 
       printf("x01=%10.9E\n",X01[i][j]);*/ 
        
       irel=(ld3*arel)*(2.0*volrel*F)*X10[i][j]*(Carel[i][j]-Cad[i][j]);
       /***************************************/
        //q=CSQtot*(1.0-OCSQ[i][j]);
        q2=0.0;
        cadi2=gsl_sf_pow_int(Cai[i][j],2);
        act12=k2i*(ld4*ld4); //(k2ai*Kd1+k11*q)/(Kd1+q);        
        act22=k2i*(ld4*ld4);// (k2ai*Kd8+k11*q)/(Kd8+q);
        inact12=k6i*ld4;//(k62*Kd3+k7*q)/(Kd3+q);
        inact22=k6i*ld4; // (k62*Kd1+k7*q)/(Kd1+q);
        dy00=km2i*Y10[i][j]-act12*cadi2*Y00[i][j]+km6i*Y01[i][j]-inact22*Cai[i][j]*Y00[i][j];
        dy10=act12*cadi2*Y00[i][j]-km2i*Y10[i][j]+km6i*Y11[i][j]-inact12*Cai[i][j]*Y10[i][j];
        dy11=inact12*Cai[i][j]*Y10[i][j]-km6i*Y11[i][j] + act2*cadi2*Y01[i][j]-km2i*Y11[i][j];
        dy01=km2i*Y11[i][j]-act2*cadi2*Y01[i][j]+inact22*Cai[i][j]*Y00[i][j]-km6i*Y01[i][j];
        Y00[i][j]+=dt*(dy00);
        Y10[i][j]+=dt*(dy10);
        Y11[i][j]+=dt*(dy11);
        Y01[i][j]+=dt*(dy01);
       /*printf("x00=%10.9E\n",X00[i][j]); 
       printf("x10=%10.9E\n",X10[i][j]); 
       printf("x11=%10.9E\n",X11[i][j]); 
       printf("x01=%10.9E\n",X01[i][j]);*/ 
       //irel2=(10.0*arel2)*Y10[i][j]*(Carel2[i][j]-Cai[i][j]);
       irel2=(ld3*arel2)*(2.0*volrel2*F)*Y10[i][j]*(Carel2[i][j]-Cai[i][j]);

      }
           
	
       DO=0.045*DOc+0.08*DOTc+0.16*DOTMgC;

       //Carel[i][j]+=dt*((Itr-irel)/(2.0*volrel*F));
       Carel[i][j]+=dt*((Itr-irel)/(2.0*volrel*F)-ld1*DOcsq);
       //Carel[i][j]+=dt*((Itr-irel)/(2.0*volrel*F));//-31.0*DOcsq);
       Carel2[i][j]+=dt*((Itr2-irel2)/(2.0*volrel2*F)-ld1*DOcsq2);
       //Carel2[i][j]+=dt*((Itr2-irel2)/(2.0*volrel2*F)); //-31.0*DOcsq2);
       Caup[i][j]+=dt*((Iup-Itr-Itr2)/(2.0*volup*F));
       //Caup[i][j]+=dt*((Iup-Itr)/(2.0*volup*F));
       
       //printf("irel=%10.9e\n",irel);
       //printf("%10.9e %lf\n",t,fmod(t,0.1));
       //getchar();
       //if((i==0)&&(fmod(t,0.1)<dt)){printf("%10.9e %10.9e %10.9e %10.9e %10.9e\n",t,Itr,Itr2,irel,irel2);
       //                                getchar();}
       //printf("Carel=%10.9e\nCaup:%10.9e\n",Carel[i][j],Caup[i][j]);
    /**********************************************/
    /*NEW EQUATIONS FOR Cai                       */
    /**********************************************/
      isi=(Cas[i][j]-Cai[i][j])*((2.0*vols*F)/tausi);
      //Cai[i][j]+=dt*(((isi+irel2-Iup)/(2.0*voli*F))-((voli/(vols+voli))*DO));
      if(CT[i][j]==1.0){
      Cai[i][j]+=dt*(((isi+irel2-Iup)/(2.0*voli*F))-((1.0*ld2*voli)/(voli+vols))*DO);}
      else{
	   Cai[i][j]+=dt*(((isi+irel2-Iup)/(2.0*voli*F))-((1.0*ld2*voli)/(voli+vols))*DO);}
      //Cai[i][j]+=dt*(((isi+irel2-Iup)/(2.0*voli*F))-DO);
       
      //if((i==0)&&(fmod(t,0.1)<dt)){printf("%10.9e %10.9e %10.9e %10.9e %10.9e %10.9e\n",t,Itr,Itr2,irel,irel2,isi);}
      //printf("isi=10.9e\n",isi);
      //printf("Cai=10.9e\n",Cai[i][j]);
    /**********************************************/
    /*Intracellular Ion Concentrations: Nai,Ki,Cai*/
    /**********************************************/
    //  DO=0.045*DOc+0.08*DOTc+0.16*DOTMgC;

      ///STANDARD NYGREN///////////////////////////
      /*if(CT[i][j]==0.0){
         
        //idi=(Cad[i][j]-Cai[i][j])*(2.0*vold*F/(10*taudi)); //taudi*10
        idi=(Cad[i][j]-Cas[i][j])*(2.0*vold*F/(10*taudi)); //taudi*10
        Cad[i][j]+=(-1.0)*dt*(ical+idi)/(2.0*F*vold);
        //Cai[i][j]+=(-1.0)*dt*((-idi+ibca+icap-2.0*inaca+Iup-irel)/(2*voli*F)+DO);
        Cas[i][j]+=(-1.0)*dt*((-idi+ibca+icap-2.0*inaca+isi)/(2*vols*F)+DO);
      }*/
      ////////////////////////////////////////////// 
      ///MODIFIED NYGREN-KEENER
      if(CT[i][j]>=1.0){
        //idi=(Cad[i][j]-Cai[i][j])*(2.0*vold*F/taudi);
        //idi=(Cad[i][j]-Cas[i][j])*(2.0*vold*F/(10.*taudi));
        idi=(Cad[i][j]-Cas[i][j])*(2.0*vold*F/(ld5*taudi));

        Cad[i][j]+=(-1.0)*dt*(ical+idi-irel)/(2.0*F*vold);
        //Cai[i][j]+=(-1.0)*dt*((-idi+ibca+icap-2.0*inaca+Iup)/(2*voli*F)+DO);
        //Cas[i][j]+=(-1.0)*dt*( ((-idi+ibca+icap-2.0*inaca+isi)/(2*vols*F))+((vols/(vols+voli))*DO) );
        //Cas[i][j]+=(-1.0)*dt*( ((-idi+ibca+icap-2.0*inaca+isi)/(2*vols*F))+DO);
        Cas[i][j]+=(-1.0)*dt*( ((-idi+ibca+icap-2.0*inaca+isi)/(2*vols*F)));
        //Cas[i][j]+=(-1.0)*dt*( ((-idi+ibca+icap-2.0*inaca+isi)/(2*vols*F)));
      }
      ////////////////////////////////////////////// 
        Nai[i][j]+=(-1.0)*dt*((ina+ibna+3.0*inak+3.0*inaca+(0.0*phinaen)-istim)/(voli*F));
        Ki[i][j]+=(-1.0)*dt*((it+isus+ik1+iks+ikr-2.0*inak)/(voli*F));
        /*printf("idi=%10.9e\n",idi);
        printf("Cad=%10.9e\nCas:%10.9e\n",Cad[i][j],Cas[i][j]);
        printf("Nai=%10.9e\nKi:%10.9e\n",Nai[i][j],Ki[i][j]);*/
   /********************/
   /*LAPLACIAN COUPLING*/
   /********************/
      Dv=0.0;
      if(eudt==1){
      if(i>0){Dv+=(v[i-1][j]-v[i][j]);}
      if(i<(Lx-1)){Dv+=(v[i+1][j]-v[i][j]);}
      if(j>0){Dv+=(v[i][j-1]-v[i][j]);}
      if(j<(Ly-1)){Dv+=(v[i][j+1]-v[i][j]);}
      Dv*=(diff*dt/(dx*dx));
      }
      vt[i][j]=v[i][j]-((dt/Cm)*(ina+ical+it+isus+ik1+ikr+iks+ibna+ibca+inak+icap+inaca-istim))+Dv;

      //getchar();
//    printf("Vn=%10.9e\n",vt[i][j]);
      }//y-loop end
  }//x-loop end

  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
     v[i][j]=vt[i][j];
   }
  } 

}
//////////////////////////////////
void prtstuff2(double X[Lx][Ly],double tp, int k){
FILE *fout;
int i,j;
char n1[32];
char dtext[11];
char dtext2[11];

strcpy(dtext2,"snaps/");

switch(k){
	  case 1:
                
                strcpy(dtext,"V.dat");
		sprintf(n1,"%s%lf%s",dtext2,tp,dtext);
		break;
         case 2:
                strcpy(dtext,"Cai.dat");
                sprintf(n1,"%s%lf%s",dtext2,tp,dtext);
		break;
         }
 
  fout=fopen(n1,"w");
  for(i=0;i<Lx;i++) {
   for(j=0;j<Ly;j++){
     fprintf(fout,"%10.9e ",X[i][j]);
    }
     fprintf(fout,"\n");
  }
  fclose(fout);
}
//////////////////////////////////
void prtstuff(double X[Lx][Ly],double tp, int k){
 FILE *fout;
 int i,j;
 char datext[20];
 switch(k){
	case 1:
		strcpy(datext,"Volt.dat");
		break;
	case 2:
                strcpy(datext,"Cad.dat");
		break;
	case 3:
		strcpy(datext,"m.dat");
                break;
	case 4:
		strcpy(datext,"h1.dat");
                break;
        case 5:
                strcpy(datext,"h2.dat");
		break;
        case 6:
		strcpy(datext,"dl.dat");
		break;
	case 7:
		strcpy(datext,"fl1.dat");
                break;
	case 8:
		strcpy(datext,"fl2.dat");
                break;
	case 9:
		strcpy(datext,"Nac.dat");
		break;
	case 10:
		strcpy(datext,"Kc.dat");
		break;
	case 11:
		strcpy(datext,"Cac.dat");
	 	break;
	case 12:
		strcpy(datext,"Nai.dat");
		break;
	case 13:
		strcpy(datext,"Ki.dat");
		break;
	case 14:
		strcpy(datext,"Cai.dat");
		break;
	case 15:
		strcpy(datext,"Caup.dat");
		break;
	case 16:
		strcpy(datext,"Carel.dat");
		break;
	case 17:
		strcpy(datext,"x00.dat");
		break;
	case 18:
		strcpy(datext,"x10.dat");
		break;
	case 19:
		strcpy(datext,"x01.dat");
		break;
	case 20:
		strcpy(datext,"x11.dat");
		break;
	case 21:
		strcpy(datext,"F1.dat");
		break;
	case 22:
		strcpy(datext,"F2.dat");
		break;
        case 23:
                strcpy(datext,"Cells.dat");
                break;
	case 24:
		strcpy(datext,"Cas.dat");
                break;
	case 25:
		strcpy(datext,"y00.dat");
		break;
	case 26:
		strcpy(datext,"y10.dat");
		break;
	case 27:
		strcpy(datext,"y01.dat");
		break;
	case 28:
		strcpy(datext,"y11.dat");
                break;
        case 29: 
                strcpy(datext,"Carel2.dat");
                break;


 } 
 
  if(tp==0){fout=fopen(datext,"w");}
  else{fout=fopen(datext,"a");}
 
  if(k!=23){ 
  fprintf(fout,"%10.9e ",tp);
  for(i=0;i<lx;i++) {
   for(j=0;j<ly;j++){
     fprintf(fout,"%10.9e ",X[i][j]);
    }}
  fprintf(fout,"\n");
  }
 
  else{
   for(i=0;i<Lx;i++){for(j=0;j<Ly;j++){if(X[i][j]==1.0){fprintf(fout,"%d %d\n",i,j);}}} 
  }

  fclose(fout);

}
//////////////////////////////////
void inicond2(void){
 int i,j,k;
 int ri,rj;
 double lk;
 double qi,qj,qr,pr,Ro,Rm,R,Rx;
 FILE *fout;

 gsl_rng_env_setup();
 Trnum = gsl_rng_default;
 rnum = gsl_rng_alloc (Trnum); 
 //ri=gsl_rng_uniform_int(rnum,Lx);
 //rj=gsl_rng_uniform_int(rnum,Ly);
 //printf("rx:%d ry=%d\n",ri,rj);
 //getchar();


 //Cell-Type Distribution
 for(i=0;i<Lx;i++){
  for(j=0;j<Ly;j++){
     
    // CT[i][j]=1.0; //1.-NoAlt 2-ALt
      
     qi=(double)((i-110)*(i-110));
     qj=(double)((j-74)*(j-74));
     qr=(double)(Ly); 
     qr*=0.5; 
     qr=qr-2.0;	
     if(((qi+qj)<(qr*qr)) &&(lk<1.0)){CT[i][j]=2.0;}
	//qr=((double)Lx)/2.0;
	//qi=gsl_pow_int((double)i,2);
	//qj=gsl_pow_int((double)qr,2);
        //pr=(qi)/(qi+qj);	
     	//lk=gsl_rng_uniform_pos(rnum);
	//if(lk<=pr){CT[i][j]=2.0;}
	else{CT[i][j]=1.0;} //1.-NoAlt 2-ALt
 }
}
 fout=fopen("Celltypes.dat","w");
 for(i=0;i<Lx;i++){
  for(j=0;j<Ly;j++){
    fprintf(fout,"%d %d %lf\n",i,j,CT[i][j]);
  }
 }

 fclose(fout); 




/* gsl_rng_env_setup();
 Trnum = gsl_rng_default;
 rnum = gsl_rng_alloc (Trnum); 
 for(i=0;i<5;)
 {
  ri=gsl_rng_uniform_int(rnum,Lx);
  rj=gsl_rng_uniform_int(rnum,Ly);
  if(((ri!=0)&&(rj!=0))||((ri!=1)&&(rj!=0))){CT[(int)ri][(int)rj]=1.0; i++;}
 }*/

 //CT[0][0]=1.0;
 //CT[1][0]=2.0;
 /*CT[25][25]=1.0;
 CT[24][25]=1.0;
 CT[26][25]=1.0;
 CT[25][24]=1.0;
 CT[25][26]=1.0;*/ 


/*
  for(i=0;i<Lx;i++){
     for(j=0;j<Ly;j++){

     if(CT[i][j]==2.0){
     v[i][j]=v[1][0];Nac[i][j]=Nac[1][0];Kc[i][j]=Kc[1][0];
     Cac[i][j]=Cac[1][0];Nai[i][j]=Nai[1][0];Ki[i][j]=Ki[1][0];
     Cai[i][j]=Cai[1][0];Cad[i][j]=Cad[1][0];Caup[i][j]=Caup[1][0];
     Carel[i][j]=Carel[1][0];m[i][j]= m[1][0];h1[i][j]=h1[1][0];
     h2[i][j]=h2[1][0];DL[i][j]=DL[1][0];fL1[i][j]=fL1[1][0];
     fL2[i][j]=fL2[1][0];r[i][j]=r[1][0];s[i][j]=s[1][0]; 
     rsus[i][j]=rsus[1][0];ssus[i][j]= ssus[1][0];n[i][j]=n[1][0]; 
     Pa[i][j]=Pa[1][0];OC[i][j]=OC[1][0];OTC[i][j]=OTC[1][0]; 
     OTMgC[i][j]=OTMgC[1][0];OTMgMg[i][j]=OTMgMg[1][0];OCSQ[i][j]=OCSQ[1][0]; 
     Cas[i][j]=Cas[1][0]; Carel2[i][j]=Carel2[1][0];  OCSQ2[i][j]=OCSQ2[1][0];
     X00[i][j]=X00[1][0];X10[i][j]=X10[1][0];X01[i][j]=X01[1][0];X11[i][j]=Y11[1][0];
     Y00[i][j]=Y00[1][0];Y10[i][j]=Y10[1][0];Y01[i][j]=Y01[1][0];Y11[i][j]=Y11[1][0];
     }   
    }	     
  }*/


   for(i=0;i<Lx;i++){
     for(j=0;j<Ly;j++){
//     if(CT[i][j]==1.0){
     v[i][j]=v[0][0]; Nac[i][j]=Nac[0][0]; Kc[i][j]=Kc[0][0];
     Cac[i][j]=Cac[0][0];Nai[i][j]=Nai[0][0];Ki[i][j]=Ki[0][0];
     Cai[i][j]=Cai[0][0];Cad[i][j]=Cad[0][0];Caup[i][j]=Caup[0][0];
     Carel[i][j]=Carel[0][0];m[i][j]= m[0][0];h1[i][j]=h1[0][0];
     h2[i][j]=h2[0][0];DL[i][j]=DL[0][0];fL1[i][j]=fL1[0][0];
     fL2[i][j]=fL2[0][0];r[i][j]=r[0][0];s[i][j]=s[0][0]; 
     rsus[i][j]=rsus[0][0];ssus[i][j]= ssus[0][0];n[i][j]=n[0][0]; 
     Pa[i][j]=Pa[0][0];OC[i][j]=OC[0][0];OTC[i][j]=OTC[0][0]; 
     OTMgC[i][j]=OTMgC[0][0];OTMgMg[i][j]=OTMgMg[0][0];OCSQ[i][j]=OCSQ[0][0]; 
     Cas[i][j]=Cas[0][0]; Carel2[i][j]=Carel2[0][0];  OCSQ2[i][j]=OCSQ2[0][0];
     X00[i][j]=X00[0][0];X10[i][j]=X10[0][0];X01[i][j]=X01[0][0];X11[i][j]=Y11[0][0];
     Y00[i][j]=Y00[0][0];Y10[i][j]=Y10[0][0];Y01[i][j]=Y01[0][0];Y11[i][j]=Y11[0][0];
     //}
   }
 }

}
//////////////////////////////////
void inicond(void){
 int i,j;
 int ri,rj;
 
 //Cell-Type Distribution

 CT[0][0]=1.0;
 CT[1][0]=2.0;

 for(i=0;i<lx;i++){
   for(j=0;j<ly;j++){
     v[i][j]=vo;
     Nac[i][j]=cnaco;
     Kc[i][j]=ckco;
     Cac[i][j]=ccaco;
     Nai[i][j]=cnaio;
     Ki[i][j]=ckio;
     Cai[i][j]=ccaio;
     Cad[i][j]=ccado;
     Caup[i][j]=ccaupo;
     Carel[i][j]=ccarelo;
     m[i][j]= mo;
     h1[i][j]=h1o;
     h2[i][j]=h2o;
     DL[i][j]=dlo;
     fL1[i][j]=fl1o;
     fL2[i][j]=fl2o;
     r[i][j]=ro;
     s[i][j]=so; 
     rsus[i][j]=rsuso; 
     ssus[i][j]= ssuso;
     n[i][j]=no; 
     Pa[i][j]=pao;
     OC[i][j]=OCo; 
     OTC[i][j]=OTCo; 
     OTMgC[i][j]=OTMgCo; 
     OTMgMg[i][j]=OTMgMgo; 
     OCSQ[i][j]=Ocsqo;
     Cas[i][j]=ccaso;
     Carel2[i][j]=ccarelo;
     OCSQ2[i][j]=Ocsqo;
    
    
     if(CT[i][j]==0.0){
     		X00[i][j]=0.0; 
     		X10[i][j]=0.0; 
     		X11[i][j]=0.0; 
     		X01[i][j]=0.0;
                Y00[i][j]=0.0; 
     		Y10[i][j]=0.0; 
     		Y11[i][j]=0.0; 
     		Y01[i][j]=0.0; 
 
     		F1[i][j]=F1o; 
     		F2[i][j]=F2o;
		}
     else{
          	X00[i][j]=x00o; 
     		X10[i][j]=x10o; 
     		X11[i][j]=x11o; 
     		X01[i][j]=x01o; 
        	Y00[i][j]=x00o; 
       		Y10[i][j]=x10o; 
     		Y11[i][j]=x11o; 
     		Y01[i][j]=x01o; 

     		F1[i][j]=0.0; 
     		F2[i][j]=0.0;
	}   
    }
  }
//
}
//////////////////////////////////
int main (void){
 double nm,ncyc,fteff,ftfx,TPR,TPR2,tpr2;
 int qn;
 //printf("HOLA\n");
 //getchar();
 TPR=20000;
 TPR2=6000;
 tpr2=1.0;
 ncyc=floor(tf/Ts); 
 fteff=Ts*ncyc;
 //fteff=10.0;
 //printf("number of cycles:%lf\n",ncyc); 
 //printf("final time:%lf secs,effective final time:%lf secs\n",tf/1000,fteff/1000); 
////////////////////
 volc=0.136*voli;
 vold=0.02*voli;
//NEW SARCOLEMMAL VOLUME
 vols=5*vold;
 volup=voli/12.34; //0.007938;//nL
 volrel2=10.0*volrel; //THEN arel in irel2 is 10*arel
 //volrel2=voli/500.0; //THEN arel in irel2 is 10*arel
//printf("volrel:%lf, volrel2:%lf\n",volrel,volrel2);
//getchar();
////////////////////////// 
 Rtf=(R*T)/F;
///////////////////
//////Phase-I//////
/////////////////// 
 t=0.0;
 inicond();

 ftfx=fteff;
 for(qn=0;qn<2;qn++){
  
  if(qn==1){t=0.0;     //Starts printing after integrating in [0,fteff]
            TPR=0.0;   /*Restarts the time to zero with the closer to 
                        the stationary state conditions*/
            //The new ftfx could be the same as the original if the next line  line is uncommented//                 
            //ftfx=fteff;
            //Or an new interval could be specified, e.g.:
            ftfx=10*Ts; 
           }
 //Time loop start
 while(t<ftfx){
  //Integrating routine
  Eulerupdate(0,qn); 
  t+=dt;

    if((t>TPR)&&(qn==1)){prtstuff(v,t,1);
            prtstuff(Cad,t,2); prtstuff(m,t,3); prtstuff(h1,t,4);
	    prtstuff(h2,t,5); prtstuff(DL,t,6); prtstuff(fL1,t,7);
            prtstuff(fL2,t,8);	prtstuff(Nac,t,9);  prtstuff(Kc,t,10); 
            prtstuff(Cac,t,11); prtstuff(Nai,t,12); prtstuff(Ki,t,13); 
            prtstuff(Cai,t,14); prtstuff(Caup,t,15);  prtstuff(Carel,t,16); 
            prtstuff(X00,t,17); prtstuff(X10,t,18); prtstuff(X01,t,19); 
            prtstuff(X11,t,20); prtstuff(Cas,t,24);
            prtstuff(Y00,t,25); prtstuff(Y10,t,26);
            prtstuff(Y01,t,27); prtstuff(Y11,t,28);
            prtstuff(Carel2,t,29); 
            //printf("%10.9e %10.9e %10.9e\n",t,istim,ina);
            TPR+=1;}
 }
}

/////////////////////
 printf("INTEGRATION-I COMPLETED\n");

/////////////////////
/////Phase-II////////
/////////////////////
 t=0.0;
 inicond2();
 printf("PHASE I completed!\n"); 
 //getchar();
 prtstuff2(v,0.0,1);
 prtstuff2(Cad,0.0,2);	
 printf("Integration in progress !\n"); 
 //getchar(); 
 //Time loop start
 
 while(t<20000){
  //Integrating routine
     Eulerupdate(1,1); 
     t+=dt;
     //printf("Time:%lf\n",t);
     //getchar();

  if(t>TPR2){
            prtstuff2(v,TPR2,1); prtstuff2(Cai,TPR2,2);
	     //printf("Time:%lf\n",t);
             TPR2+=5;
            }
  //if(t>tpr2){printf("%lf\n",tpr2);tpr2+=1.0;}
 }//Time loop end
 
 printf("DONE!\n");
 return 0;
}

