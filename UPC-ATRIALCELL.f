cNOTA: 31-5-13 cmbio en linea 588 en eca ccai por ccas.
c    gcal=0.7*gcal;km2=2*km2; 3-7-13
c Recordar que cuando el tiempo sea largo se tiene que cambiar el maxbeats
c   debo cambiar xcap =xcapt(indexcai))
c  inicializar ccas
c  poner stim   

c Nygren-keener et al. model in a single cell
c steady-state restitution
      program nygren0dssrest
      implicit none
c nx is maximum ring size
      integer nnx
      parameter(nnx=1)
      real*8 v(0:nnx),vt(0:nnx),ccai(0:nnx),ccad(0:nnx),ccaup(0:nnx),
     &     ccarel(0:nnx),xm(0:nnx),xh1(0:nnx),xh2(0:nnx),xdl(0:nnx),
     &     xfl1(0:nnx),xfl2(0:nnx),xr(0:nnx),xs(0:nnx),
     &     rsus(0:nnx),ssus(0:nnx),xn(0:nnx),ccas(0:nnx),
     &     xpa(0:nnx),xf1(0:nnx),xf2(0:nnx),xoc(0:nnx),xotc(0:nnx),
     &     xotmgc(0:nnx),xotmgmg(0:nnx),xocalse(0:nnx),cnai(0:nnx),
     &     cki(0:nnx),cnac(0:nnx),ckc(0:nnx),ccac(0:nnx),
     &     x00(0:nnx),x10(0:nnx),x01(0:nnx),x11(0:nnx),ccarel2(0:nnx),
     &     y00(0:nnx),y10(0:nnx),y01(0:nnx),y11(0:nnx),xocalse2(0:nnx)

c this works for 0d
      real*8 xstimamp,xstimdur
c for xstimdur=3, 2*diastolic thresh=19.0 (for i<5 in 1d, di thresh=9.5)
      parameter(xstimdur=3.0)

      integer count,ierr

c parameters
      real*8 cnab,ckb,ccab,cmgi,ecaapp,xkca,xxr,xxt,xxf,cm,
     &     voli,voldfrac,volrel,volup,tauna,tauk,tauca,taudi,
     &     xinak,xknakk,xknakna,xicap,xkcap,xknaca,gamma,dnaca,
     &     xiup,xkcyca,xksrca,xkxcs,tautr,alpharel,xkreli,xkreld,
     &     rrecov,pna,gcal,gt,gsus,gks,gkr,gk1,gbna,gbca,
     &     phinaenorig,volc,CSQtot,arel,arel2,volrel2,vols
      parameter(cnab=130.0,ckb=5.4,ccab=1.8,cmgi=2.5,ecaapp=60.0,
     &     xkca=0.025,xxr=8314.,xxt=306.15,xxf=96487.,cm=1.0,
     &     voli=0.11768,voldfrac=0.02,volrel=0.000882,volup=0.007938,
     &     tauna=14300.,tauk=10000.,tauca=24700.,taudi=1.,
c taudi=10.,
     &     xinak=1.416506,xknakk=1.,xknakna=11.,xicap=0.08,
     &     xkcap=0.0002,xknaca=0.000749684,gamma=0.45,dnaca=0.0003,
     &     xiup=56.,xkcyca=0.0003,xksrca=0.5,xkxcs=0.4,tautr=10.,
     &     alpharel=4000.,xkreli=0.0003,xkreld=0.003,rrecov=0.815e-3,
     &     pna=3.2e-5,gcal=.7*0.135,gt=0.15,gsus=0.055,
     &     gks=0.02,gkr=0.01,
     &     gk1=0.06,gbna=0.00121198,gbca=0.00157362,phinaenorig=-0.0336,
     &     volc=0.01600448,CSQtot=1.5,arel=0.3,arel2=0.3,
     &     volrel2=10.*volrel,vols=5.*voldfrac*voli)
c Parameters for new compartments and diffussion between them
      real*8 tausi,ld1,ld2,ld3,ld4,ld5,ld6,km2,km6,k2,k6
      parameter(tausi=0.1,ld1=1.0,ld2=0.5,ld3=1.3,ld4=70.0,ld5=1.0/100,
     &          ld6=0.,km2=2.*0.06,km6=1./200,k2=45.,k6=0.3)
c values needed for currents
      real*8 phinaen,rtof,fort,vold,xknakna15,
     &     otauna,otauk,otauca,ovolcf,o2volupf,o2volrelf,o2volsf,
     &     o2voldf,ovolif,otautr,otaudi,xkxcs2srca,xkxcssrca,
     &     dtocm,fca,cnai3,cnac3,cnai15,docdt,dotcdt,dotmgcdt,
     &     dotmgmgdt,docalsedt,ract,rinact,temp,dccareldt,
     &     dccaupdt,dccarel2dt,dodt,docalse2dt,o2volrelf2
c needed for computing current integrals
c      real*8 qna,qbna,qt,qsus,qk1,qks,qkr,qcal,qbca,qnakna,
c     &     qnakk,qcap,qnacana,qnacaca,qnaen,qstim
c currents
      real*8 xna,xcal,xt,xsus,xks,xkr,xk1,xbna,xbca,xnak,xcap,
     &     xnaca,xstim,xup,xtr,xrel,xidi,xisi,dx00,dx10,dx01,
     &     dx11,dy00,dy10,dy01,dy11,xrel2,xtr2
c numerical values
      real*8 dt,endtime,ena,ek,eca

      integer minsteps,imax,idtmin,idtmed,idtmax,istimdur,
     &     iS1delay,istep,BEAT,S1step,nx,nend

c tables
      real*8 dvtable,dcaitable,dvmektable,dvmenatable,dckctable,
     &     dccadtable,odvtable,odcaitable,odvmektable,odvmenatable,
     &     odckctable,odccadtable,vv,ccaii,vmek,ckcc,ccadd,vmena,
     &     dccastable,odccastable,ccass
      integer indexv,indexcai,indexvmek,indexvmena,indexckc,indexccad
      integer  indexccas
      real*8 vlo,vhi
      integer nvt
      parameter(vlo=-90.,vhi=60.,nvt=5000)
      real*8 xmbar(0:nvt),otaum(0:nvt),xhbar(0:nvt),otauh1(0:nvt),
     &       otauh2(0:nvt),dlbar(0:nvt),otaudl(0:nvt),flbar(0:nvt),
     &       otaufl1(0:nvt),otaufl2(0:nvt),rbar(0:nvt),otaur(0:nvt),
     &       sbar(0:nvt),otaus(0:nvt),rsusbar(0:nvt),otaursus(0:nvt),
     &       ssusbar(0:nvt),otaussus(0:nvt),xnbar(0:nvt),otaun(0:nvt),
     &       xpabar(0:nvt),otaupa(0:nvt),xpi(0:nvt),
     &       xnacoeff(0:nvt),xcalcoeff(0:nvt),
     &     xnaca1(0:nvt),xnaca2(0:nvt)
      real*8 cailo,caihi
      integer ncai
      parameter(cailo=.00001,caihi=.01001,ncai=10000)
c      real*8 ecat(0:ncai),xcapt(0:ncai),ccaipow(0:ncai)
      real*8 ecat(0:ncai),ccaipow(0:ncai)
      real*8 vmeklo,vmekhi
      integer nvmekt
      parameter(vmeklo=0.0,vmekhi=150.0,nvmekt=5000)
      real*8 xk1v(0:nvmekt)
      real*8 vmenalo,vmenahi
      integer nvmenat
      parameter(vmenalo=-158.0,vmenahi=-8.0,nvmenat=5000)
      real*8 xna1(0:nvmenat)
      real*8 ckclo,ckchi
      integer nckct
      parameter(ckclo=4.0,ckchi=14.0,nckct=10000)
      real*8 ckcpow(0:nckct)
      real*8 ccadlo,ccadhi
      real*8 ccaslo,ccashi
      integer nccadt
      integer nccast
      parameter(ccadlo=0.0,ccadhi=0.1,nccadt=1000)
      parameter(ccaslo=0.0,ccashi=0.1,nccast=1000)
      real*8 ccadpow(0:nccadt)
      real*8 xcapt(0:nccast)
      real*8 dvdtmax,timefin

c threshold used to measure apd from di
      real*8 vr,vru
c variables for using APD percentage instead of fixed voltage
      real*8 vmin1,vmin2,vmin3,vmin4,vmax1,vmax2,vmax3,vmax4,
     &     vthresh1,vthresh2,vthresh3,vthresh4,vpercent,
     &     dvdt1,dvdt2,dvdt3,dvdt4,dvdtold1,dvdtold2,dvdtold3,dvdtold4,
     &     dvdtolder1,dvdtolder2,dvdtolder3,dvdtolder4,vmaxend
      parameter(vpercent=0.9)
      integer ncl
c      parameter(vr=-69.0,vru=-69.0)
      character*40 voltfile
      character*40 outfile
      parameter(ncl=73)
c      parameter(ncl=100)
      real*8 cls(ncl)
      real*8 dx,diff,ddt_o_dx2,xlap
      integer icl1,icl,icycle,ncycles,i,
     &     intcl
c use this to start from the beginning
c      parameter(icl1=1)
c use this to start from a previously output file
      parameter(icl1=0)
      real*8 maxv,time,cyclelength,apd,di,cl,prevapd,prevdi,prevcl,
     &     cvfront,cvback
      integer istart,ntime,idig1,idig2,idig3,idig4,istp,idecr
      integer nups1,nups2,nups3,nups4,ndowns1,ndowns2,ndowns3,ndowns4,
     &     i1,maxbeats
      integer nvups1,nvups2,nvups3,nvups4,
     &     nvdowns1,nvdowns2,nvdowns3,nvdowns4
      parameter(maxbeats=10000)
      real*8 ups1(maxbeats),ups2(maxbeats),ups3(maxbeats),
     &     ups4(maxbeats),downs1(maxbeats),downs2(maxbeats),
     &     downs3(maxbeats),downs4(maxbeats)
      real*8 vups1(maxbeats),vups2(maxbeats),vups3(maxbeats),
     &     vups4(maxbeats),vdowns1(maxbeats),vdowns2(maxbeats),
     &     vdowns3(maxbeats),vdowns4(maxbeats)

c here specify how to update intracellular and cleft ionic 
c concentrations
c intracellular include cnai and cki
c extracellular include cnac, ckc, and ccac
c 1 = allow to vary freely as originally specified
c 2 = hold fixed as constant (cnai(t+dt)=cnai(t))
      integer ifixed,icourt,iaf
      parameter(ifixed=0)
c use Courtemanche/AF parameters? 0=no,1=yes
      parameter(icourt=0,iaf=0)
c before each new CL, re-initialize all variables to starting values?
      integer iinit
c      parameter(iinit=1)
c   si no queremos reinicializar
      parameter(iinit=0)
c-------------
      dvdtmax=0.0
c-------------

      phinaen=0.0
c these values were taken as 2*diastolic thresh in single cell, 3ms duration
c vr is APD90 measured after 30s pacing at 1000 ms CL
      if(iaf.eq.1.and.icourt.eq.0.and.ifixed.eq.0) then
c AF
         xstimamp=20.0
         vr=-67.03
c         vr=-63.44
      elseif(iaf.eq.1.and.icourt.eq.0.and.ifixed.eq.1) then
c AF-fixed
         xstimamp=19.8
         vr=-66.59
c         vr=-62.99
      elseif(iaf.eq.0.and.icourt.eq.1.and.ifixed.eq.0) then
c Court
         xstimamp=19.8
c         vr=-60.16
         vr=-63.27
      elseif(iaf.eq.0.and.icourt.eq.1.and.ifixed.eq.1) then
c Court-fixed
         xstimamp=19.6
         vr=-60.82
      elseif(ifixed.eq.0) then
c reg
         xstimamp=19.0
c         vr=-62.47
c PROPER ONE
         vr=-62.64
      else
c fixed
         xstimamp=19.6
c         vr=-62.42
c         vr=-66.49
          vr=-63.5
      endif
      vru=vr

      cls(1)=1000.
      cls(2)=950.
      cls(3)=900.
      cls(4)=850.
      cls(5)=800.
      cls(6)=750.
      cls(7)=700.
      cls(8)=650.
      cls(9)=600.
      cls(10)=575.
      cls(11)=550.
      cls(12)=525.
      cls(13)=500.
      cls(14)=480.
      cls(15)=460.
      cls(16)=440.
      cls(17)=420.
      cls(18)=400.
      cls(19)=390.
      cls(20)=380.
      cls(21)=370.
      cls(22)=360.
      cls(23)=350.
      cls(24)=340.
      cls(25)=330.
      cls(26)=320.
      cls(27)=310.
      cls(28)=300.
      cls(29)=290.
      cls(30)=280.
      cls(31)=270.
      cls(32)=260.
      cls(33)=250.
      cls(34)=240.
      cls(35)=230.
      cls(36)=220.
      cls(37)=210.
      cls(38)=200.
      cls(39)=190.
      cls(40)=180.
      cls(41)=170.
      cls(42)=160.
      cls(43)=155.
      cls(44)=150.
      cls(45)=145.
      cls(46)=140.
      cls(47)=135.
      cls(48)=130.
      cls(49)=125.
      cls(50)=120.
      cls(51)=115.
      cls(52)=110.
      cls(53)=105.
      cls(54)=100.
      cls(55)=95.
      cls(56)=90.
      cls(57)=85.
      cls(58)=80.
      cls(59)=75.
      cls(60)=70.
      cls(61)=65.
      cls(62)=60.
      cls(63)=55.
      cls(64)=50.
      cls(65)=45.
      cls(66)=40.
      cls(67)=35.
      cls(68)=30.
      cls(69)=25.
      cls(70)=20.
      cls(71)=15.
      cls(72)=10.
      cls(73)=5.

c      cls(1)=1000.
c      do i=2,ncl
c         cls(i)=cls(i-1)-10.
c      enddo

      open(51,file='apdevolution.dat',form='formatted',status='unknown')
      open(31,file='res0dss.dat',form='formatted',status='unknown')
      open(41,file='res0dssv.dat',form='formatted',status='unknown')

c numerical parameters
      nx=1
      dt=.01
      dx=.025
      diff=0.001
      ddt_o_dx2=dt*diff/(dx*dx)
      
      dvtable=(vhi-vlo)/nvt
      dcaitable=(caihi-cailo)/ncai
      dvmektable=(vmekhi-vmeklo)/nvmekt
      dvmenatable=(vmenahi-vmenalo)/nvmenat
      dckctable=(ckchi-ckclo)/nckct
      dccadtable=(ccadhi-ccadlo)/nccadt
      dccastable=(ccashi-ccaslo)/nccast

      odvtable=1.0/dvtable
      odcaitable=1.0/dcaitable
      odvmektable=1.0/dvmektable
      odvmenatable=1.0/dvmenatable
      odckctable=1.0/dckctable
      odccadtable=1.0/dccadtable
      odccastable=1.0/dccastable

c      endtime=356.0
c      endtime=800000.
      endtime=200000.
      istimdur=nint(xstimdur/dt)
      icl=1
      cyclelength=cls(icl)
      icycle=nint(cyclelength/dt)
c do whole cycles, don't stop in the middle of one
      ncycles=nint(endtime/cyclelength)
      if(ncycles.lt.endtime/cyclelength) then
         ncycles=ncycles+1
      endif
      nend=icycle*ncycles
c      nend=endtime/dt
      write(6,*) 'write n cycles of length x ', ncycles,cyclelength
      write(6,*) 'write nend', nend

c      dvt=(vhi-vlo)/nvt
c      dcait=(caihi-cailo)/ncait
c      write(6,*) 'dvt, dcait = ', dvt, dcait

c initial conditions
      if(icl1.eq.1) then
         do i=0,nx
            v(i)=-74.2525
            vt(i)=v(i)
            cnac(i)=130.0
            ckc(i)=5.4
            ccac(i)=1.8
            cnai(i)=8.55474
            cki(i)=129.435
            ccai(i)=6.72905e-5
            ccas(i)=6.72905e-5
            ccad(i)=7.24947e-5
            ccaup(i)=0.664564
            ccarel(i)=0.646458
            xm(i)=3.20173e-3
            xh1(i)=0.881421
            xh2(i)=0.874233
            xdl(i)=1.30047e-5
            xfl1(i)=0.998639
            xfl2(i)=0.998604
            xr(i)=1.06783e-3
            xs(i)=0.948975
            rsus(i)=1.59491e-4
            ssus(i)=0.991169
            xn(i)=4.83573e-3
            xpa(i)=5.16114e-5
            xf1(i)=0.428396
            xf2(i)=0.0028005
            xoc(i)=0.0274971
            xotc(i)=0.0132801
            xotmgc(i)=0.196085
            xotmgmg(i)=0.709417
            xocalse(i)=0.43686
	    x00(i)=0.998
	    x10(i)=0.
	    x11(i)=0.
	    x01(i)=0.02
	    y00(i)=0.998
	    y10(i)=0.
	    y11(i)=0.
	    y01(i)=0.02	    
         enddo
      else
         open(7,file='check.dat',status='unknown',form='unformatted')
         read(7) v,vt,cnac,ckc,ccac,cnai,cki,ccai,ccad,ccaup,ccarel,
     &     xm,xh1,xh2,xdl,xfl1,xfl2,xr,xs,rsus,ssus,xn,xpa,
     &     xf1,xf2,xoc,xotc,xotmgc,xotmgmg,xocalse,ccas,x00,
     &     x10,x01,x11,ccarel2,y00,y10,y01,y11,xocalse2
      endif
 

c useful values
      rtof=xxr*xxt/xxf  
      fort=1.0/rtof
      vold=voldfrac*voli
      xknakna15=xknakna**1.5
      otauna=1.0/tauna
      otauk=1.0/tauk
      otauca=1.0/tauca
      ovolcf=1.0/(volc*xxf)
      o2volupf=1.0/(2.0*volup*xxf)
      o2volrelf=1.0/(2.0*volrel*xxf)
      o2volrelf2=1.0/(2.0*volrel2*xxf)
      o2voldf=1.0/(2.0*vold*xxf)
      o2volsf=1.0/(2.0*vols*xxf)
      ovolif=1.0/(voli*xxf)
      otautr=1.0/tautr
      otaudi=1.0/taudi
      xkxcs2srca=9.6e-5         ! xkxcs*xkxcs*xkcyca /xksrca
      xkxcssrca=0.00024         ! xkxcs*xkcyca/xksrca

      dtocm=dt/cm
      istimdur=nint(xstimdur/dt)

c make tables
c nvt had better = (vhi-vlo)/dvtable
      do i=0,nvt
        vv=vlo+i*dvtable
        xmbar(i)=1.0/(1.0+exp((vv+27.12)/(-8.21)))
        otaum(i)=1.0/(0.042*exp(-((vv+25.57)/28.8)**2)+0.024)
        xhbar(i)=1.0/(1.0+exp((vv+63.6)/5.3))
        otauh1(i)=1.0/(30.0/(1.0+exp((vv+35.1)/3.2))+0.3)
        otauh2(i)=1.0/(120.0/(1.0+exp((vv+35.1)/3.2))+3.0)
        dlbar(i)=1.0/(1.0+exp((vv+9.0)/(-5.8)))
        otaudl(i)=1.0/(2.7*exp(-((vv+35.0)/30.0)**2)+2.0)
        flbar(i)=1.0/(1.0+exp((vv+27.4)/7.1))
        otaufl1(i)=1.0/(161.0*exp(-((vv+40.0)/14.4)**2)+10.0)
        otaufl2(i)=1.0/(1332.3*exp(-((vv+40.0)/14.2)**2)+62.6)
        rbar(i)=1.0/(1.0+exp((vv-1.0)/(-11.0)))
        otaur(i)=1.0/(3.5*exp(-(vv/30.0)**2)+1.5)
        sbar(i)=1.0/(1.0+exp((vv+40.5)/11.5)) 
        otaus(i)=1.0/(481.2*exp(-((vv+52.45)/14.97)**2)+14.14)
        rsusbar(i)=1.0/(1.0+exp((vv+4.3)/(-8.0)))
        otaursus(i)=1.0/(9.0/(1.0+exp((vv+5.0)/12.0))+0.5)
        ssusbar(i)=0.4/(1.0+exp((vv+20.0)/10.0))+0.6
        otaussus(i)=1.0/(47.0/(1.0+exp((vv+60.0)/10.0))+300.0)
        xnbar(i)=1.0/(1.0+exp((vv-19.9)/(-12.7)))
        otaun(i)=1.0/(700.0+400.0*exp(-((vv-20.0)/20.0)**2))
        xpabar(i)=1.0/(1.0+exp((vv+15.0)/(-6.0)))
        otaupa(i)=1.0/(31.18+217.18*exp(-((vv+20.1376)/22.1996)**2))
        xpi(i)=1.0/(1.0+exp((vv+55.0)/24.0))
        if(vv.ne.0.0) then
          xnacoeff(i)=pna*vv*xxf/rtof/
     &           (exp(vv/rtof)-1.0)
        else
          xnacoeff(i)=pna*vv*xxf/rtof/
     &           (exp(vv+0.001/rtof)-1.0)
        endif
        xcalcoeff(i)=gcal*(vv-ecaapp)
        xnaca1(i)=xknaca*exp(gamma*vv/rtof)
        xnaca2(i)=xknaca*exp((gamma-1.0)*vv/rtof)
      enddo

      do i=0,ncai
	ccaii=cailo+i*dcaitable
c         xcapt(i)=xicap*(ccaii/(ccaii+xkcap))
	ccaipow(i)=(ccaii/(ccaii+xkreli))**4
      enddo

      do i=0,nccast
	 ccass=ccaslo+i*dccastable
         xcapt(i)=xicap*(ccass/(ccass+xkcap))
      enddo

      do i=0,nvmekt
         vmek=vmeklo+i*dvmektable
         xk1v(i)=gk1*vmek/(1.0+exp((1.5*(vmek+3.6))/rtof))
      enddo

      do i=0,nckct
         ckcc=ckclo+i*dckctable
         ckcpow(i)=ckcc**(0.4457)
c         write(19,*) i,ckcc,ckcpow(i)
      enddo

      do i=0,nccadt
         ccadd=ccadlo+i*dccadtable
         ccadpow(i)=(ccadd/(ccadd+xkreld))**4
      enddo

      do i=0,nvmenat
         vmena=vmenalo+i*dvmenatable
         if(abs(vmena).gt.1e-3) then
            xna1(i)=exp(vmena/rtof)-1.0
         else
            xna1(i)=xna1(i-1)
         endif
      enddo

c recording sites
      i1=1

      dvdt1=0.0
      dvdtold1=0.0
      dvdtolder1=0.0
      vmin1=v(i1)
      vmax1=v(i1)

      voltfile='volt.dat.xxxx'
      outfile='output.dat.xxxx'

c      write(6,*) 'starting time loop'
C-_-_-_-_-_-_-_-_ TIME LOOP -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

      do icl=1,ncl
       dvdtmax=0.0
c      do icl=36,ncl
         intcl=nint(cls(icl))
         write(6,*) intcl
         idig1=intcl/1000
         istp=mod(intcl,1000)
         idig2=istp/100
         istp=mod(istp,100)
         idig3=istp/10
         idig4=mod(istp,10)
         voltfile(10:10)=char(ichar('0')+idig1)
         voltfile(11:11)=char(ichar('0')+idig2)
         voltfile(12:12)=char(ichar('0')+idig3)
         voltfile(13:13)=char(ichar('0')+idig4)
         outfile(12:12)=char(ichar('0')+idig1)
         outfile(13:13)=char(ichar('0')+idig2)
         outfile(14:14)=char(ichar('0')+idig3)
         outfile(15:15)=char(ichar('0')+idig4)
         write(6,*) voltfile
         open(80,file=voltfile,status='unknown')

         cyclelength=cls(icl)
         icycle=nint(cyclelength/dt)
         ncycles=nint(endtime/cyclelength)
         if(ncycles.lt.endtime/cyclelength) then
            ncycles=ncycles+1
         endif
         nend=icycle*ncycles
         write(6,*) 'write n cycles of length x ',ncycles,cyclelength
c            timefin=(endtime-2.*cyclelength) 
c             print*,'timefin en if',timefin

         nups1=0
         ndowns1=0
         nups2=0
         ndowns2=0
         nups3=0
         ndowns3=0
         nups4=0
         ndowns4=0
         nvups1=0
         nvdowns1=0
         nvups2=0
         nvdowns2=0
         nvups3=0
         nvdowns3=0
         nvups4=0
         nvdowns4=0

c re-initialize to initial values if desired
         if(iinit.eq.1) then
            do i=1,nx
               v(i)=-74.2525
               vt(i)=v(i)
               cnac(i)=130.0
               ckc(i)=5.4
               ccac(i)=1.8
               cnai(i)=8.55474
               cki(i)=129.435
               ccai(i)=6.72905e-5
               ccas(i)=6.72905e-5
               ccad(i)=7.24947e-5
               ccaup(i)=0.664564
               ccarel(i)=0.646458
               xm(i)=3.20173e-3
               xh1(i)=0.881421
               xh2(i)=0.874233
               xdl(i)=1.30047e-5
               xfl1(i)=0.998639
               xfl2(i)=0.998604
               xr(i)=1.06783e-3
               xs(i)=0.948975
               rsus(i)=1.59491e-4
               ssus(i)=0.991169
               xn(i)=4.83573e-3
               xpa(i)=5.16114e-5
               xf1(i)=0.428396
               xf2(i)=0.0028005
               xoc(i)=0.0274971
               xotc(i)=0.0132801
               xotmgc(i)=0.196085
               xotmgmg(i)=0.709417
               xocalse(i)=0.43686
	       x00(i)=0.998
	       x10(i)=0.
	       x11(i)=0.
	       x01(i)=0.02
	       y00(i)=0.998
	       y10(i)=0.
	       y11(i)=0.
	       y01(i)=0.02
            enddo
         endif

c      do ntime=1,nend

         ntime=0
         time=0.0
         vmaxend=0.0

c 707     if(ntime.ge.nend.and.vmaxend.ge.-40.0) then
 707     continue

         ntime=ntime+1
         time=ntime*dt

         if(mod(ntime,1000000).eq.1) write(6,*) time
c         if(time.ge.45100) then
c            write(6,*) time,ntime,nvdowns1,nvups1
c         endif
c         ierr=4510001

c         if(ntime.eq.10) goto 800

         do i=1,nx

c reversal potentials
            ena=rtof*log(cnac(i)/cnai(i))
            ek=rtof*log(ckc(i)/cki(i))
c            eca=0.5*rtof*log(ccac(i)/ccai(i))
!!!!!!!Cuidado
c            if(ntime.gt.112723800) print*,ccac(i),ccas(i),ccad(i)
            eca=0.5*rtof*log(ccac(i)/ccas(i))
c            if(ntime.gt.112723800) print*,'Eca=',eca


c stimulus
            if(mod(ntime,icycle).le.istimdur.and.i.lt.5) then
c            if(icl1.eq.1.and.ntime.le.istimdur.and.i.lt.5) then
             if(icl.lt.28) then
               xstim=xstimamp
             else
               xstim=xstimamp*1.2
             endif
            else
               xstim=0.0
            endif

c indices into tables
            indexv=nint((v(i)-vlo)*odvtable)
            indexcai=nint((ccai(i)-cailo)*odcaitable)
            indexvmek=nint((v(i)-ek-vmeklo)*odvmektable)
            indexvmena=nint((v(i)-ena-vmenalo)*odvmenatable)
            indexckc=nint((ckc(i)-ckclo)*odckctable)
            indexccad=nint((ccad(i)-ccadlo)*odccadtable)
            indexccas=nint((ccas(i)-ccaslo)*odccastable)
            if (indexv.gt.nvt.or.indexv.lt.0) then
               write(6,*), 'v outside range',
     &           ntime,ntime*dt,v(i)
               goto 800
            endif
            if (indexcai.gt.ncai)write(6,*),'cai outside range',
     &           ntime,ccai(i),ccac(i),log(ccac(i)/ccas(i))
c     &           ntime,ntime*dt,ccai(i)
            if (indexvmek.gt.nvmekt)write(6,*) 'v-ek outside range',
     &           ntime,ntime*dt,v(i)-ek
            if (indexvmena.gt.nvmenat)write(6,*) 'v-ena outside range',
     &           ntime,ntime*dt,v(i)-ena
c            if (indexckc.gt.nckct) write(6,*) 'ckc outside range',
c     &           ntime,ntime*dt,ckc(i)
            if (indexccad.gt.nccadt) write(6,*) 'ccad outside range',
     &           ntime,ntime*dt,ccad(i)
            if (indexccas.gt.nccast) write(6,*) 'ccas outside range',
     &           ntime,ntime*dt,ccas(i)

c sodium current
            xm(i)=xm(i)+dt*(xmbar(indexv)-xm(i))*otaum(indexv)
            xh1(i)=xh1(i)+dt*(xhbar(indexv)-xh1(i))*otauh1(indexv)
            xh2(i)=xh2(i)+dt*(xhbar(indexv)-xh2(i))*otauh2(indexv)

            xna=xnacoeff(indexv)*xm(i)*xm(i)*xm(i)*cnac(i)*
     &           (0.9*xh1(i)+0.1*xh2(i))*xna1(indexvmena)

c L-type calcium current
            fca=ccad(i)/(ccad(i)+xkca)

            xdl(i)=xdl(i)+dt*(dlbar(indexv)-xdl(i))*otaudl(indexv)
            xfl1(i)=xfl1(i)+dt*(flbar(indexv)-xfl1(i))*otaufl1(indexv)
            xfl2(i)=xfl2(i)+dt*(flbar(indexv)-xfl2(i))*otaufl2(indexv)

            xcal=xcalcoeff(indexv)*xdl(i)
     &           *(fca*xfl1(i)+(1.0-fca)*xfl2(i))

c transient and sustained outward k+ currents
            xr(i)=xr(i)+dt*(rbar(indexv)-xr(i))*otaur(indexv)
            xs(i)=xs(i)+dt*(sbar(indexv)-xs(i))*otaus(indexv)

            rsus(i)=rsus(i)
     &           +dt*(rsusbar(indexv)-rsus(i))*otaursus(indexv)
            ssus(i)=ssus(i)
     &           +dt*(ssusbar(indexv)-ssus(i))*otaussus(indexv)

            xt=gt*xr(i)*xs(i)*(v(i)-ek)
            xsus=gsus*rsus(i)*ssus(i)*(v(i)-ek)
          
c delayed rectifier k+ currents
            xn(i)=xn(i)+dt*(xnbar(indexv)-xn(i))*otaun(indexv)
            xpa(i)=xpa(i)+dt*(xpabar(indexv)-xpa(i))*otaupa(indexv)

            xks=gks*xn(i)*(v(i)-ek)
            xkr=gkr*xpa(i)*xpi(indexv)*(v(i)-ek)

c inward rectifier k+ current
            if(ckc(i).lt.0) write(6,*) 'xk1 will crash, ckc negative'
            if(indexckc.ge.0.and.indexckc.le.nckct) then
               xk1=xk1v(indexvmek)*ckcpow(indexckc)
            else
               xk1=xk1v(indexvmek)*ckc(i)**(0.4457)
            endif

c if using Courtemanche parameters, resize currents
            if(icourt.eq.1) then
               xcal=xcal*1.33
               xt=xt*2.0
               xsus=xsus*0.4
               xkr=xkr*3.0
               xks=xks*3.0
            endif

c if using AF parameters, resize currents
            if(iaf.eq.1) then
               xcal=xcal*0.3
               xt=xt*0.5
               xsus=xsus*0.5
            endif

c background inward currents
            xbna=gbna*(v(i)-ena)
            xbca=gbca*(v(i)-eca)

c pump and exchanger currents
            cnai3=cnai(i)*cnai(i)*cnai(i)
            cnac3=cnac(i)*cnac(i)*cnac(i)
            cnai15=sqrt(cnai3)
            xnak=xinak*ckc(i)*cnai15*(v(i)+150.0)/
     &           ((ckc(i)+xknakk)*(cnai15+xknakna15)*(v(i)+200.0))
c          xnak=xinak*(ckc(i)/(ckc(i)+xknakk))
c     &              *(cnai15/(cnai15+xknakna15))
c     &              *(v(i)+150.0)/(v(i)+200.0)

c---- pongo la expresion para probar. Poner despues en tablas
c            xcap=xicap*(ccas(i)/(ccas(i)+xkcap))
            xcap=xcapt(indexccas)

            xnaca=(cnai3*ccac(i)*xnaca1(indexv)
     &           -cnac3*ccas(i)*xnaca2(indexv))/
     &           (1.0+dnaca*(cnac3*ccas(i)+cnai3*ccac(i)))

c intracellular ca2+ buffering

            docdt=200.0*ccai(i)*(1.0-xoc(i))-0.476*xoc(i)
            dotcdt=78.4*ccai(i)*(1.0-xotc(i))-0.392*xotc(i)
            dotmgcdt=200.0*ccai(i)*(1.0-xotmgc(i)-xotmgmg(i))
     &           -0.0066*xotmgc(i)
            dotmgmgdt=2.0*cmgi*(1.0-xotmgc(i)-xotmgmg(i))
     &           -0.666*xotmgmg(i)
            
            xoc(i)=xoc(i)+dt*docdt
            xotc(i)=xotc(i)+dt*dotcdt
            xotmgc(i)=xotmgc(i)+dt*dotmgcdt
            xotmgmg(i)=xotmgmg(i)+dt*dotmgmgdt

c cleft space ion concentrations
            if(ifixed.eq.0) then
               cnac(i)=cnac(i)+dt*((cnab-cnac(i))*otauna
c     &              +(xna+xbna+3.0*xnak+3.0*xnaca+phinaen)
     &              +(xna+xbna+3.0*xnak+3.0*xnaca+phinaen-xstim)
     &              *(ovolcf))
               ckc(i)=ckc(i)+dt*((ckb-ckc(i))*otauk
     &              +(xt+xsus+xk1+xks+xkr-2.0*xnak)*(ovolcf))
               ccac(i)=ccac(i)+dt*((ccab-ccac(i))*otauca
     &              +(xcal+xbca+xcap-2.0*xnaca)*0.5*(ovolcf))
            endif

c ca2+ handling by sarcoplasmic reticulum
            docalsedt=0.48*ccarel(i)*(1.0-xocalse(i))-0.4*xocalse(i)
            docalse2dt=0.48*ccarel2(i)*(1.0-xocalse2(i))-0.4*xocalse2(i)
            xocalse(i)=xocalse(i)+dt*docalsedt
            xocalse2(i)=xocalse2(i)+dt*docalse2dt

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ract=0.2038*(ccaipow(indexcai)+ccadpow(indexccad))
c           rinact=0.03396+0.3396*(ccaipow(indexcai))
            ract=k2*(ld4*ld4)*ccad(i)*ccad(i)
            rinact=k6*ld4*ccad(i)
c            xf1(i)=xf1(i)+dt*(rrecov*(1.0-xf1(i)-xf2(i))-ract*xf1(i))
c            xf2(i)=xf2(i)+dt*(ract*xf1(i)-rinact*xf2(i))

        dx00=km2*X10(i)-ract*X00(i)+km6*X01(i)-rinact*X00(i)
        dx10=ract*X00(i)-km2*X10(i)+km6*X11(i)-rinact*X10(i)
        dx11=rinact*X10(i)-km6*X11(i) + ract*X01(i)-km2*X11(i)
c No deberiamos hacer la resta        
        dx01=km2*X11(i)-ract*X01(i)+rinact*X00(i)-km6*X01(i)
        X00(i)=X00(i)+dt*(dx00)
        X10(i)=X10(i)+dt*(dx10)
        X11(i)=X11(i)+dt*(dx11)
        X01(i)=X01(i)+dt*(dx01)
c        X01(i)=1.-(X00(i)+X10(i)+X11(i))
            ract=k2*(ld4*ld4)*ccai(i)*ccai(i)
            rinact=k6*ld4*ccai(i)
        dy00=km2*Y10(i)-ract*Y00(i)+km6*Y01(i)-rinact*Y00(i)
        dy10=ract*Y00(i)-km2*Y10(i)+km6*Y11(i)-rinact*Y10(i)
        dy11=rinact*Y10(i)-km6*Y11(i) + ract*Y01(i)-km2*Y11(i)
        dy01=km2*Y11(i)-ract*Y01(i)+rinact*Y00(i)-km6*Y01(i)
        Y00(i)=Y00(i)+dt*(dy00)
        Y10(i)=Y10(i)+dt*(dy10)
        Y11(i)=Y11(i)+dt*(dy11)
       Y01(i)=Y01(i)+dt*(dy01)
c        Y01(i)=1.-(Y00(i)+Y10(i)+Y11(i))

            xup=xiup*(ccai(i)-xkxcs2srca*ccaup(i))/
     &           (ccai(i)+xkcyca+xkxcssrca*(ccaup(i)+xksrca))
            xtr=(ccaup(i)-ccarel(i))*2.0*xxf*volrel*otautr
            xtr2=(ccaup(i)-ccarel2(i))*2.0*xxf*volrel2*otautr
c            temp=xf2(i)/(xf2(i)+0.25)
c            xrel=alpharel*(temp*temp)*(ccarel(i)-ccai(i))
          xrel=(ld3*arel/o2volrelf)*X10(i)*(ccarel(i)-ccad(i))
          xrel2=(ld3*arel2/o2volrelf2)*Y10(i)*(ccarel2(i)-ccai(i))
c          if(ntime.gt.112723800) print*,'xrel=',xrel,'xrel2=',xrel2


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            docalsedt=0.48*ccarel(i)*(1.0-xocalse(i))-0.4*xocalse(i)
            docalse2dt=0.48*ccarel2(i)*(1.0-xocalse2(i))-0.4*xocalse2(i)
c            dccareldt=(xtr-xrel)*o2volrelf-31.0*docalsedt
c            dccareldt2=(xtr2-xrel2)*o2volrelf2-31.0*docalsedt2
            dccareldt=(xtr-xrel)*o2volrelf-ld1*docalsedt
            dccarel2dt=(xtr2-xrel2)*o2volrelf2-ld1*docalse2dt
            dccaupdt=(xup-xtr-xtr2)*o2volupf

            ccarel(i)=ccarel(i)+dt*dccareldt
            ccarel2(i)=ccarel2(i)+dt*dccarel2dt
            if (ccarel(i).lt.0.0) write(6,*) 'CALCIUM-REL NEGATIVE ',
     &           ccarel(i),i,ntime,ntime*dt
            ccaup(i)=ccaup(i)+dt*dccaupdt
            if (ccaup(i).lt.0.0) write(6,*) 'CALCIUM-UP NEGATIVE ',
     &           ccaup(i),i,ntime,ntime*dt

c diadic and subsarcolema space components
c            xidi=(ccad(i)-ccai(i))*2.0*xxf*vold*otaudi
            xidi=(ccad(i)-ccas(i))*2.0*xxf*vold*otaudi/ld5
            xisi=(ccas(i)-ccai(i))*2.0*xxf*vols/tausi
            ccad(i)=ccad(i)-dt*(xcal+xidi-xrel)*o2voldf
            ccas(i)=ccas(i)-dt*(xcap+xbca-xidi-2.*xnaca+xisi)*o2volsf
            if (ccad(i).lt.0.0) write(6,*) 'CALCIUM-D NEGATIVE ',
     &           ccad(i),i,ntime,ntime*dt
            if (ccas(i).lt.0.0) write(6,*) 'CALCIUM-S NEGATIVE ',
     &           ccas(i),i,ntime,ntime*dt,xcap,xbca,-xidi,-xnaca,xisi

c may need to update the dotcdt, dotmgcdt, docdt with new O vals
c here now recomputing them
            docdt=200.0*ccai(i)*(1.0-xoc(i))-0.476*xoc(i)
            dotcdt=78.4*ccai(i)*(1.0-xotc(i))-0.392*xotc(i)
            dotmgcdt=200.0*ccai(i)*(1.0-xotmgc(i)-xotmgmg(i))
     &           -0.0066*xotmgc(i)

c trying this here for now
            dodt=0.08*dotcdt+0.16*dotmgcdt+0.045*docdt

c intracellular ion concentrations
c            ccai(i)=ccai(i)-dt*(-xidi+xbca+xcap-2.0*xnaca+xup-xrel)
c     &           *0.5*ovolif-dt*dodt
            ccai(i)=ccai(i)-dt*(-xisi+xup-xrel2)*0.5*ovolif
     &           -dt*((1.0*ld2*voli)/(voli+vols))*dodt
            if(ifixed.eq.0) then
               cnai(i)=cnai(i)
c     &              -dt*(xna+xbna+3.0*xnak+3.0*xnaca+phinaen)
     &              -dt*(xna+xbna+3.0*xnak+3.0*xnaca+phinaen-xstim)
     &              *ovolif
               cki(i)=cki(i)-dt*(xt+xsus+xk1+xks+xkr-2.0*xnak)*ovolif
            endif

c compute laplacian
cc          if(icl1.eq.1.and.ntime.lt.nint(200./dt)) then
c             if(i.eq.0) then
c                xlap=ddt_o_dx2*(-2.*v(0)+2.*v(1))
c             elseif(i.eq.nx) then
c                xlap=ddt_o_dx2*(-2.*v(nx)+2.*v(nx-1))
c             else
c                xlap=ddt_o_dx2*(-2.*v(i)+v(i+1)+v(i-1))
c             endif
cc          else
cc             if(i.eq.0) then
cc                xlap=ddt_o_dx2*(-2.*v(0)+v(1)+v(nx))
cc             elseif(i.eq.nx) then
cc                xlap=ddt_o_dx2*(-2.*v(nx)+v(nx-1)+v(0))
cc             else
cc                xlap=ddt_o_dx2*(-2.*v(i)+v(i+1)+v(i-1))
cc             endif
cc          endif
            xlap=0.0
            vt(i)=v(i) -dtocm*(xna+xcal+xt+xsus+xkr+xks+xk1
     &           +xbna+xbca+xnak+xcap+xnaca-xstim)

c for variable method, compute dvdt to find vmax and dvdt max
             if(i.eq.i1) then
                dvdtolder1=dvdtold1
                dvdtold1=dvdt1
c                dvdt1=(vt(i1)-v(i1))
                dvdt1=(vt(i1)-v(i1))/dt
             endif
       if(dvdt1.gt.dvdtmax) dvdtmax=dvdt1

c record vmin as needed
             if(i.eq.i1.and.dvdt1.ge.0.0.and.dvdtold1.lt.0.0.
     &       and.v(i1).lt.-20.0) then
                vmin1=v(i1)
c                write(6,*) 'vmin1 = ',vmin1,ntime*dt
c                write(6,*) '   ',vt(i1),v(i1),
c     &               dvdt1,dvdtold1
             endif

c record vmax and compute vthresh
             if(i.eq.i1.and.dvdt1.le.0.and.dvdtold1.gt.0.0
     &            .and.abs(dvdt1).gt.1e-2
     &            .and.abs(dvdtold1).gt.1e-2) then
                vmax1=v(i1)
                vthresh1=vmin1+(vmax1-vmin1)*(1.0-vpercent)
c                if(icl.eq.1) then
c                   write(6,*) '   ',ntime*dt,vthresh1,dvdt1,dvdtold1
c                endif
c                write(6,*) 'vmax1 = ',vmax1,ntime*dt
c                write(6,*) '   ', vt(i1),v(i1),
c     &               dvdt1,dvdtold1
c                write(6,*) 'vthresh1 = ',vthresh1,vmin1,vmax1,vpercent
             endif

c check for variable APD threshold crossings
             if(i.eq.i1.and.dvdtold1.gt.dvdtolder1.and.
     &            dvdt1.le.dvdtold1.and.dabs(dvdt1).gt.1.0.and.
     &            dabs(dvdtolder1).gt.1.0.and.nvups1.eq.nvdowns1) then
                nvups1=nvups1+1
                vups1(nvups1)=(ntime-1)*dt
c                write(6,*) 'site 1 upward crossing ',ntime*dt,vt(i1),
c     &               v(i1)
c                write(6,*) '   ',dvdt1,dvdtold1,dvdtolder1
             elseif(i.eq.i1.and.
     &               v(i1).gt.vthresh1.and.vt(i1).le.vthresh1.and.
     &               nvups1.eq.nvdowns1+1) then
                nvdowns1=nvdowns1+1
                vdowns1(nvdowns1)=ntime*dt
     &               +dt*(vthresh1-v(i1))/(vt(i1)-v(i1))
c                write(6,*) 'site 1 downward crossing ',ntime*dt,
c     &               v(i1),vt(i1),vthresh1,nvdowns1
             endif

c          if(mod(ntime,40).eq.1.and.i.eq.1) then
c          if(i.eq.1) then
c             write(6,*) ntime*dt,v(i),vt(i),xstim
c             write(6,*) '   ',xna,xk1,xto,xkur
c             write(6,*) '   ',xkr,xks,xcal,xpca
c             write(6,*) '   ',xnak,xnaca,xbna,xbca
c             write(6,*) '   ',xrel,xtr,xup,xupleak
c             write(6,*) '   ',xu(i),xv(i),xw(i)
c             write(6,*) '     ',xinfu,exptauu,fn,ifn1
c             write(6,*) '        ',vrel,cm,xrel
c             write(6,*) '   ',ccai(i),ccaup(i),ccarel(i)
c          endif
 
c         if(mod(ntime,40).eq.1.and.i.eq.40) then
c            write(81,*) ntime*dt,xina,xik1,xito,xikp
c            write(82,*) ntime*dt,xinab,xiks,xica,xinaca
c            write(83,*) ntime*dt,xipca,xicab,xicak,xinak
c            write(84,*) ntime*dt,xikr,ccai(i),ccasr(i),xlap
cc            write(6,*) vt(0),vt(40),vt(80)
c         endif

        enddo

c check for APD threshold crossings
        if(v(i1).lt.vru.and.vt(i1).ge.vr.and.
     &       nups1.eq.ndowns1) then
           nups1=nups1+1
           ups1(nups1)=ntime*dt+dt*(vr-v(i1))/(vt(i1)-v(i1))
        elseif(v(i1).gt.vr.and.vt(i1).le.vr.and.
     &          nups1.eq.ndowns1+1) then
           ndowns1=ndowns1+1
           downs1(ndowns1)=ntime*dt+dt*(vr-v(i1))/(vt(i1)-v(i1))
	apd=downs1(ndowns1)-ups1(ndowns1)
	write(51,*) ntime*dt,apd
        endif

        do i=1,nx
           v(i)=vt(i)
        enddo

c measure maximum voltage at end of ring where cut would take place        
c        if (icl.le.ncl) then
c           vmaxend=-100.0
c           istart=nx-ndecrease(idecr+1)
c           do i=istart,nx
c              if(vmaxend.lt.v(i)) vmaxend=v(i)
c           enddo
cc           if(mod(ntime,100).eq.1) write(6,*) time,istart,nx,vmaxend
c        endif

        if(mod(ntime,50).eq.1) then
c     Creo que ahora funcionara, he multiplicado ntime*dt
           if(ntime*dt.gt.(endtime-5.*cyclelength)) then
c           if(float(ntime).gt.798000000) then
c           if(float(ntime).gt.250000000) then
           write(80,*) ntime*dt,v(1),ccai(1)
c           write(88,*) ntime*dt,dvdt1
c        if(mod(ntime,50).eq.1) then
c           write(80,*) ntime*dt,v(0),v(nx)
c           write(80,*) ntime*dt,v(1)
           write(95,*) ntime*dt,ccad(1),ccas(1),ccai(1)
c            write(95,*) ntime*dt,ccac(1),ccas(1),ccai(1)
        endif
        endif

c        if (ntime.le.nend.or.vmaxend.ge.-40.0) goto 707
        if(ntime.le.nend) then
           goto 707
        endif

c        endif

c      enddo
c end of time loop for given cycle length

c write out restart information (also useful for S1-S2 protocol)
      open(8,file=outfile,form='unformatted',status='unknown')
      write(8) v,vt,cnac,ckc,ccac,cnai,cki,ccai,ccad,ccaup,ccarel,
     &     xm,xh1,xh2,xdl,xfl1,xfl2,xr,xs,rsus,ssus,xn,xpa,
     &     xf1,xf2,xoc,xotc,xotmgc,xotmgmg,xocalse,ccas,x00,
     &     x10,x01,x11,ccarel2,y00,y10,y01,y11,xocalse2
      close(8)

c find APD's, DI's, and CV's - voltage threshold
c site 1 (i1)
      apd=downs1(ndowns1)-ups1(ndowns1)
      prevapd=downs1(ndowns1-1)-ups1(ndowns1-1)
      di=ups1(ndowns1)-downs1(ndowns1-1)
      prevdi=ups1(ndowns1-1)-downs1(ndowns1-2)
      cl=ups1(ndowns1)-ups1(ndowns1-1)
      prevcl=ups1(ndowns1-1)-ups1(ndowns1-2)
      write(31,*) di,apd,cl
      write(31,*) prevdi,prevapd,prevcl
c      write(31,*) 'dvdtmax = ',dvdtmax
      write(6,*) 'di,apd,cl = ',di,apd,cl
      write(6,*) 'di,apd,cl = ',prevdi,prevapd,prevcl
      write(6,*) 'dvdtmax = ',dvdtmax

c find APD's, DI's, and CV's - variable threshold (true APD90 or whatever)
c site 1 (i1)

      apd=vdowns1(nvdowns1)-vups1(nvdowns1)
      prevapd=vdowns1(nvdowns1-1)-vups1(nvdowns1-1)
      di=vups1(nvdowns1)-vdowns1(nvdowns1-1)
      prevdi=vups1(nvdowns1-1)-vdowns1(nvdowns1-2)
      cl=vups1(nvdowns1)-vups1(nvdowns1-1)
      prevcl=vups1(nvdowns1-1)-vups1(nvdowns1-2)
c      apd=vdowns1(nvups1)-vups1(nvups1)
c      prevapd=vdowns1(nvups1-1)-vups1(nvups1-1)
c      di=vups1(nvups1)-vdowns1(nvups1-1)
c      prevdi=vups1(nvups1-1)-vdowns1(nvups1-2)
c      cl=vups1(nvups1)-vups1(nvups1-1)
c      prevcl=vups1(nvups1-1)-vups1(nvups1-2)
      write(41,*) di,apd,cl,vthresh1
      write(41,*) prevdi,prevapd,prevcl,vthresh1
c      write(41,*) 'dvdtmax = ',dvdtmax
c      write(81,*) cl,vthresh1,vups1(nvdowns1),vdowns1(nvdowns1)
      write(6,*) 'at site 1: ',di,apd,cl,vthresh1

c      maxv=-100.0
c      icl1=0
c      do i=1,nx
c         if(v(i).gt.maxv) maxv=v(i)
c      enddo
c      if(maxv.lt.-40.0) then
c         write(6,*) 'maximum voltage = ', maxv
c         write(6,*) 'ending program'
c         goto 800
c      endif

c      icl=icl+1

c   Cuidado!!! Mejor comentado, podria salir del bucle
c      if(cl.gt.cyclelength*1.5) goto 800

      if(icl.lt.ncl) write(6,*) 'new cl = ',cls(icl+1)

      enddo

 800  continue

      close(31)
c      close(32)
c      close(33)
c      close(34)
c      close(35)
c      close(36)

      end











