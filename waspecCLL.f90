IMPLICIT DOUBLE PRECISION (a-h,o-z)

! Supplementary material of Moretti, Le Losq, Neuville, Geochimica Cosmochimica 2014
! Original code written by Roberto Moretti in F77
! Code converted by Charles Le Losq in F90
! Date: 2014-02-12  Time: 11:12:39

REAL :: k
PARAMETER (ndi=600)
! Declaration of variable in double precision
INTEGER, PARAMETER :: DP = KIND(0.0D0)
REAL(KIND=DP) :: x(49), delta(ndi), tcent(ndi), o(ndi), o1(ndi), o3(ndi), o4(ndi)
real(KIND=DP) :: zcat(13), xcat(ndi,13), aossi(ndi), ersta(ndi), z(13), y(13)
real(KIND=DP) :: tkelv(ndi), pbar(ndi), rappox(ndi), totcat(ndi), totani(ndi)
real(KIND=DP) :: xti(ndi), xmn(ndi), xfe(ndi), xmg(ndi), xca(ndi), xna(ndi)
real(KIND=DP) :: xk(ndi), xsi(ndi), xal(ndi), poss(ndi,13), xossi(ndi), redox(ndi)
real(KIND=DP) ::  acidic(ndi), w(ndi), ah(ndi), somm(ndi), som(ndi), aossiz(ndi)
real(KIND=DP) :: xh(ndi), xwd(ndi), xoh(ndi), xcat0(ndi), xohz(ndi), xw(ndi)
real(KIND=DP) :: ERR(ndi), dv(ndi), delv(ndi), wmol(ndi), watmol(ndi), dx(ndi)
real(KIND=DP) :: xleft(ndi), partw(ndi), wwmol(ndi), wattot(ndi), freeox(ndi)
real(KIND=DP) :: prot(ndi), ohtrue(ndi), pipp(ndi), cost22(ndi), cost20(ndi)
real(KIND=DP) :: xmol(ndi,13), cri(ndi)


OPEN (UNIT=8,FILE='INPUT.txt',STATUS='old')
OPEN (UNIT=19,FILE='OUTPUT.txt',STATUS='UNKNOWN')


! attenzione basicity moderating parameters per il calcolo della
DATA x/2.913128,-1723.12,0.856683,-2033.56,0.,0.,0.,0.,0.,  &
    0.645,-3.06217395577529,872.441431890575,0.,1.4,0.78,0.,0.,  &
    -132.,0.,2.,0.,4.545,2.92,1.37,1.35,2.5,0.,0.5,-1.,-1.5,  &
    0.,1.41600470139180,-3568.68070901824,0.,0.,-2.80718661275843,  &
    0.,4.68206041948470,0.,-1.3,0.,0.,0.,-2.14703618087306,  &
    -2532.05418107778,0.,0.,0.,0./
DATA z/2.09,1.54,2.5,1.67,2.09,1.82,1.354,1.40,1.,1.28,.87, .71,2.5/
DATA y/60.085,79.899,141.94,101.96,159.69,151.99,71.846,  &
    70.937,56.079,40.311,61.979,94.203,18.015/
DATA zcat/1.,1.,2.,2.,2.,2.,1.,1.,1.,1.,2.,2.,2./

! Chemical composition of melt in oxyde wt%
!  (1)  (2)  (3)  (4)   (5)   (6)  (7) (8) (9) (10) (11)(12)(13)
! SiO2 TiO2 P2O5 Al2O3 Fe2O3 Cr2O3 FeO MnO CaO MgO Na2O K2O H2O

READ (8,*) ncomp,kympa
READ (8,*) ! this is to skype the header line with compositions
WRITE (*,*) ncomp
DO  i=1,ncomp
  xossi(i)=-12
  READ (8,*) tcent(i),pbar(i),poss(i,1),poss(i,2),  &
      poss(i,4),poss(i,5),poss(i,6),poss(i,7),poss(i,8),poss(i,10),  &
      poss(i,9),poss(i,11),poss(i,12),poss(i,3), poss(i,13),wwmol(i)
  tkelv(i)=tcent(i)+273.15
  IF (poss(i,5) <= 0.) poss(i,5)=0.00001
  IF (poss(i,7) <= 0.) poss(i,7)=0.00001
  zum=0.
  DO  jk=1,12
    zum=zum+poss(i,jk)
  END DO
  zym=0.
  DO  kj=1,12
    poss(i,kj)=poss(i,kj)*(100.-poss(i,13))/zum
    zym=zym+poss(i,kj)
  END DO
  zym=zym+poss(i,13)
  IF (DABS(zym-100.) >= 0.00001) PAUSE ! consistency check
  
! calcola le proporzioni molari a 100% catione e ossido
  d=0
  dd=0
  DO  j=1,13
    xcat(i,j)=poss(i,j)*zcat(j)/y(j)
    d=d+xcat(i,j)
    dd=dd+poss(i,j)/y(j)
    xmol(i,j)=poss(i,j)/y(j)
  END DO
  som(i)=d
  xcat0(i)=xcat(i,13)
!  considero tutto H2O come OH (o H!)
  dd=dd-poss(i,13)/y(13)
  xw(i)=2.*poss(i,13)/y(13)
  dd=dd+xw(i)
  somm(i)=dd
  DO  j=1,13
    xcat(i,j)=xcat(i,j)/d
    xmol(i,j)=xmol(i,j)/dd
  END DO
  
END DO

DO  i=1,ncomp
  
  k=0
  e=0
  u=0
  v=0
!  Inizia il ricalcolo per H+ e OH-
!  OH-  anione libero di H2O sulla matrice anionica
!  H+ entra nel computo della depolimerizzazione e non ci interessa
!  per niente se poi si rimette con OH- dando H2O molecolare o con O-
!  dando l'OH di Fraser che chiude le terminazioni polimeriche.
!  Di seguito approssimiamo i volumi di H+ e O2- dai raggi ionici...
  
  frit=0.
  nit=50
  kflag=0.
  
  DO  jk=1,nit
    
    xcat(i,1)=poss(i,1)*zcat(1)/y(1)
    xcat(i,3)=poss(i,3)*zcat(3)/y(3)
    xcat(i,4)=poss(i,4)*zcat(4)/y(4)
    xcat(i,5)=poss(i,5)*zcat(5)/y(5)
    xcat(i,6)=poss(i,6)*zcat(6)/y(6)
    xcat(i,7)=poss(i,7)*zcat(7)/y(7)
    xcat(i,2)=poss(i,2)*zcat(2)/y(2)
    xcat(i,13)=poss(i,13)*zcat(13)/y(13)
    xcat(i,8)=poss(i,8)*zcat(8)/y(8)
    xcat(i,9)=poss(i,9)*zcat(9)/y(9)
    xcat(i,10)=poss(i,10)*zcat(10)/y(10)
    xcat(i,11)=poss(i,11)*zcat(11)/y(11)
    xcat(i,12)=poss(i,12)*zcat(12)/y(12)
    xcat(i,1)=xcat(i,1)/som(i)
    xcat(i,3)=xcat(i,3)/som(i)
    xcat(i,4)=xcat(i,4)/som(i)
    xcat(i,5)=xcat(i,5)/som(i)
    xcat(i,6)=xcat(i,6)/som(i)
    xcat(i,2)=xcat(i,2)/som(i)
    xcat(i,7)=xcat(i,7)/som(i)
    xcat(i,8)=xcat(i,8)/som(i)
    xcat(i,9)=xcat(i,9)/som(i)
    xcat(i,10)=xcat(i,10)/som(i)
    xcat(i,11)=xcat(i,11)/som(i)
    xcat(i,12)=xcat(i,12)/som(i)
    xcat(i,13)=xcat(i,13)/som(i)
    xsi(i)=xcat(i,1)*4.
    xal(i)=xcat(i,4)*3.
    xsi(i)=xsi(i)/(xsi(i)+xal(i))
    xal(i)=xal(i)/(xsi(i)+xal(i))
    
    xti(i)=xcat(i,2)*4.
    xna(i)=xcat(i,11)
    xk(i)=xcat(i,12)
    xca(i)=xcat(i,9)*2.
    xmg(i)=xcat(i,10)*2.
    xmn(i)=xcat(i,8)*2.
    xfe(i)=xcat(i,7)*2.
    totale=xti(i)+xna(i)+xk(i)+xca(i)+xmg(i)+xmn(i)+xfe(i)
    xti(i)=xti(i)/totale
    xna(i)=xna(i)/totale
    xk(i)=xk(i)/totale
    xca(i)=xca(i)/totale
    xmg(i)=xmg(i)/totale
    xmn(i)=xmn(i)/totale
    xfe(i)=xfe(i)/totale
    xfe(i)=xfe(i)+xmn(i)
    
    aossiz(i)=aossi(i)
    IF (jk == 1) THEN
      aossiz(i)=1.
      xoh(i)=xcat0(i)
    END IF
    
    xh(i)=xcat0(i)-xoh(i)
    xohz(i)=xoh(i)
    xcat(i,13)=xh(i)/som(i)
    xwd(i)=xoh(i)/som(i)
    
! bilancio cariche creazione complessi MAl4+, ecc.
    w(i)=3.0*xcat(i,6)+2.0*(xcat(i,7)+xcat(i,8)+xcat(i,9)+xcat(i,10))  &
        +xcat(i,11)+xcat(i,12)+xcat(i,13)+xcat(i,3)+4.0*xcat(i,2)
    
    IF (xcat(i,4) > w(i)) GO TO 540
    xfor=xcat(i,1)+xcat(i,4)/2.0
    w(i)=w(i)-xcat(i,4)
    xallu=xcat(i,4)/2.0
    xcat(i,4)=0.
    GO TO 590
    540  u=xcat(i,4)-w(i)
    xfor=xcat(i,1)+w(i)/2.0
    xcat(i,4)=u
    xallu=xcat(i,4)/2.0
    w(i)=0.
!      write (*,*) 'ci casco per il composto ',i
    GO TO 680
    590  IF (xcat(i,5) > w(i)) GO TO 640
    xfor=xfor+xcat(i,5)/2.0
    w(i)=w(i)-xcat(i,5)
    xfe2o3=xcat(i,5)/2.0
    xcat(i,5)=0.
    GO TO 680
    640  v=xcat(i,5)-w(i)
    xfor=xfor+w(i)/2.0
    xcat(i,5)=v
    xfe2o3=w(i)/2.0
    w(i)=0.
!     write (*,*) 'ci casco per il composto ',i
    680  xfor=xfor+xcat(i,3)/2.0
    
!     Compute optic basicity of network former (z prefissati)
    basfor=xcat(i,1)*z(1)/xfor+xcat(i,3)*z(3)/xfor+xallu*z(4)/xfor+  &
        xfe2o3*z(5)/xfor
    
! Compute proportion of network modifiers and the
! associated basicity
    xmod=0.0
    DO  ij=2,13
      xmod=xmod+xcat(i,ij)/zcat(ij)
    END DO
    
    basmod=0.0
    DO  ij=2,12
      xcat(i,ij)=xcat(i,ij)/(zcat(ij)*xmod)
      
!   basicit ottica prefissata
      
      basmod=basmod+xcat(i,ij)*z(ij)
    END DO
    xcat(i,13)=xcat(i,13)/(zcat(13)*xmod)
    basmod=basmod+xcat(i,13)*x(26)
    k=DEXP(((basmod-basfor)/0.2145)-1.1445+0.*(pbar(i)-1)/tkelv(i))
    a=1.0-4.0*k
    xcat(i,1)=xfor/(xfor+xmod)
    acidic(i)=xcat(i,1)
    e=1.0-acidic(i)
    totcat(i)=e
    
    atoop=-a
    btoop=2.0+2.0*acidic(i)
    ctoop=8.0*acidic(i)*(acidic(i)-1.0)
    o(i)=(-btoop+DSQRT(btoop**2.-4.0*atoop*ctoop))/(2.0*atoop)
    
! calcola moli di O2-, Oø e polianioni
! deriva la frazione di O2- sulla matrice anionica
! O1  il numero di O=, O3  il numero di Oø, O(i)  il numero di O-
    
    o1(i)=1.0-acidic(i)-o(i)/2.0
    o3(i)=(4.0*acidic(i)-o(i))/2.0
    o4(i)=o(i)/(o(i)+o3(i)+acidic(i))
    s=DEXP(-1.7165*DLOG(o4(i))+2.8776)
    s1=acidic(i)/s
    o2=o1(i)/(o1(i)+s1+xwd(i))
    totani(i)=o1(i)+s1+xwd(i)
    aossi(i)=o2
    ah(i)=totcat(i)/totani(i)
    
!  I volumi sotto sono in joule/bar ...fattore 0.1...
!  Raggi ionici di Shannon
    voh=(4./3.)*3.14159*((x(14)+x(21)+x(7)/1.d6* (tkelv(i)-298.15))**3)
    voh=voh*0.6022045
    
    vo2=(4./3.)*3.14159*((x(14)+x(7)/1.d6*(tkelv(i)-298.15))**3)
    vo2=vo2*0.6022045
    
    vh=(4./3.)*3.14159*((x(21)+x(16)/1.d6*(tkelv(i)-298.15))**3)
    vh=vh*0.6022045
    
    IF (vh < 0) vh=0.
    
! Calcolo delv assumendo espansione termica = 0.
    delv(i)=vh+vo2-voh
    delv2=delv(i)*0.1/(8.3147*2.303)
    delv2=delv2*(pbar(i)-1.)/tkelv(i)
    
    partna=10.**((1.-xmol(i,13))*  &
        (x(36)+x(12)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
    partk=10.**((1.-xmol(i,13))*  &
        (x(11)+x(12)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
    partti=10.**((1.-xmol(i,13))*  &
        (x(38)+x(39)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
    partca=10.**((1.-xmol(i,13))*  &
        (x(40)+x(41)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
    partmg=10.**((1.-xmol(i,13))*  &
        (x(40)+x(41)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
    partfe=10.**((1.-xmol(i,13))*  &
        (x(44)+x(45)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/  &
        8.3147/tkelv(i)))
!      partsi=10.**((1.-xmol(i,13))*
!     &(x(46)+x(47)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/
!     &8.3147/tkelv(i)))
!      partal=10.**((1.-xmol(i,13))*
!     &(x(48)+x(49)/tkelv(i)-delv2+x(28)*x(27)*pbar(i-1.)/2.303/
!     &8.3147/tkelv(i)))
    
    partw(i)=(partna**(xna(i)))*(partk**(xk(i)))*(partti**(xti(i)))*  &
        (partca**(xca(i)))*(partmg**(xmg(i)))*(partfe**(xfe(i)))
!     &*(partsi**(xsi(i)))*(partal**(xal(i)))
!     &*partsi
    
    ratius=(aossi(i)/ah(i))/partw(i)
! ratius=(aossi(i)/ah(i)*totcat(i))/partw(i)
!     xoh(i)=xcat0(i)*ratius/(ratius+1.)
    xoh(i)=(xcat0(i)/som(i))*ratius/(1.+ratius) ! modified x scalare su 1 mole
    xoh(i)=xoh(i)*som(i)  !inserted x scalare sul tutto mole
    xh(i)=xcat0(i)-xoh(i)
    xoh(i)=DSQRT(xohz(i)*xoh(i))
    
    frit=DABS(aossi(i)-aossiz(i))
    IF (frit < 0.00001.AND.jk > 1) THEN
      kflag = 3
    END IF
!      write (*,*) 'KFLAG ',kflag
    IF (kflag == 3) EXIT
    
    IF (jk == nit) THEN
!       write (*,*) ' STRONZO...SU', i
      PAUSE
      EXIT
      
    END IF
  END DO
  792   CONTINUE
  
! speculazione su H+ e O-
  voh=(4./3.)*3.14159*((x(14)+x(21)+x(7)/1.d6* (tkelv(i)-298.15))**3)
  voh=voh*0.6022045
!  sarebbe x(7) il coeff lineare di espansione
  vo=(4./3.)*3.14159*((x(14)+x(7)/1.d6*(tkelv(i)-298.15))**3)
  vo=vo*0.6022045
!  sarebbe x(16) il coeff lineare di espansione
  vh=(4./3.)*3.14159*((x(21)+x(16)/1.d6*(tkelv(i)-298.15))**3)
  vh=vh*0.6022045
  IF (vh < 0) vh=0.
  dv(i)=-vh-vo+voh
  dv1=dv(i)*0.1/(8.3147*2.303)
  dv1=dv1*(pbar(i)-1)/tkelv(i)
  costx=10.**(x(18)+x(19)/tkelv(i)-dv1)
  aoin=o(i)
  ahin=xh(i)/som(i)
  aspec=costx
  bspec=-(costx*ahin+costx*aoin+totcat(i))
  cspec=costx*ahin*aoin
  root1=(-bspec-DSQRT(bspec**2.-4.*aspec*cspec))/(2.*aspec)
!      root2=(-bspec+dsqrt(bspec**2.-4.*aspec*cspec))/(2.*aspec)
! ATTENZIONE EFFETTUO UN CAMBIO IMPORTANTE:
  o(i)=o(i)-root1
  
! calcola moli di O2-, Oø e polianioni
! deriva la frazione di O2- sulla matrice anionica
! O1  il numero di O=, O3  il numero di Oø, O(i)  il numero di O-
  
  o1(i)=1.0-acidic(i)-o(i)/2.0
  o3(i)=(4.0*acidic(i)-o(i))/2.0
  o4(i)=o(i)/(o(i)+o3(i)+acidic(i)+root1)
!      root(i)=root1
  s=DEXP(-1.7165*DLOG(o4(i))+2.8776)
!      si(i)=s
  s1=acidic(i)/s
  o2=o1(i)/(o1(i)+s1+xwd(i))
  totani(i)=o1(i)+s1+xwd(i)
  aossi(i)=o2
  
  790   cost54=-2.8792+6364.8/tkelv(i)
  cost57=x(1)+x(2)/tkelv(i)
  cost58=x(3)+x(4)/tkelv(i)
  cost59=1.1529-1622.4/tkelv(i)
  
!   VFeO2 usando i CR (assumo struttura vincolata e non free anions)
!   0.63  HS CR di Shannon per Fe3+ (IR sarebbe 0.49) in coord IV
!   0.55  LS IR di Shannon per Fe3+ (HS sarebbe 0.645) in coord VI
!   0.61  LS IR di Shannon per Fe2+ (HS sarebbe 0.780) in coord VI
!   1.24  CR di Shannon per O2- (IR sarebbe 1.38) in coord IV
!   1.26  CR di Shannon per O2- (IR  1.40...) in coord VI
  vfeo2m=4./3.*3.14159*((2*x(14)+0.49+2*x(7)/1.d6)**3)
  vo2m=4./3.*3.14159*((x(14)+x(7)/1.d6*(tkelv(i)-298.15))**3)
  vfe3p=4./3.*3.14159*((x(10)+x(8)/1.d6*(tkelv(i)-298.15))**3)
  vfe2p=4./3.*3.14159*((x(15)+x(8)/1.d6*(tkelv(i)-298.15))**3)
  conv=0.6022045
  vfeo2m=vfeo2m*conv
  vo2m=vo2m*conv
  vfeo2m=vfeo2m-x(20)*vo2m
  vfe3p=vfe3p*conv
  vfe2p=vfe2p*conv
  
! espansione da trovare @ 298K
  vfeo15=21.065+x(22)*0.001*(298-1673)
  vfeo=13.65+x(23)*0.001*(298-1673)
! volumi da trovare @ 298K
  dv57=vfeo2m-vfeo15-0.5*vo2m
  dv58=vfe3p+1.5*vo2m-vfeo15
  dv59=vfe2p+vo2m-vfeo
  
  cost57=cost57-(dv57*0.1)*(pbar(i)-1.)/(tkelv(i)*8.3147*2.303)
  cost57=cost57
  cost58=cost58-(dv58*0.1)*(pbar(i)-1.)/(tkelv(i)*8.3147*2.303)
  cost58=cost58
  cost59=cost59-(dv59*0.1)*(pbar(i)-1.)/(tkelv(i)*8.3147*2.303)
  cost59=cost59
! Therefore:
  cost57=10.**cost57
  cost58=10.**cost58
  cost59=10.**cost59
  
! dai dati di Lange (1994)
  v07=13.65
  v05=42.13
  eoxm7=2.92
  eoxm5=9.09
  coxm7=-0.45
  coxm5=-2.53
  
  v0fe2=v07+eoxm7*0.001*(tkelv(i)-1673)
  v0fe3=v05+eoxm5*0.001*(tkelv(i)-1673)
  
!  Passaggio ai joule/bar
  v0fe2=0.1*v0fe2
  v0fe3=0.1*v0fe3
  
  cfe2=coxm7*0.0001
  cfe3=coxm5*0.0001
  
  cfe2=cfe2*0.1
  cfe3=cfe3*0.1
  
  dvs=(0.5*v0fe3-v0fe2)*(pbar(i)-1)+  &
      (0.5*cfe3-cfe2)*((pbar(i)**2)/2-pbar(i)+0.5)
  
  
  dvs=dvs/(8.3147*tkelv(i)*2.303)
!     dvs=0
!      cost54=cost54+x(27)*0.*(pbar(i)-1.)
  cost54=10.**(cost54-dvs)
  
  xoss =10.0**xossi(i)
  den1=cost54*xoss**0.25
  den2=cost57*aossi(i)**2*totani(i)+cost58*totcat(i)
  ratio=aossi(i)**(0.5)*(cost59)*totcat(i)/(den1*den2)
  
  
  cost20(i)=DEXP(x(32)+x(33)/tkelv(i))*  &
      DEXP(-x(34)*(pbar(i)-1.)/(8.41726*tkelv(i)))*DEXP(x(35))
  cost22(i)=1./partw(i)
  freeox(i)=aossi(i)*totani(i)
  
  wmol(i)=poss(i,13)*0.001
  
  kk=5000
  IF (kympa == 0)  kk=1
  324     DO  jh=1,kk
    
    IF (kympa == 0)  wmol(i)=wwmol(i)/18.015
    IF (jh > 1) THEN
      wmol(i)=(DABS(redox(i)*rappox(i)**11.))**(1./12.)
!        wmol(i)=(dabs(redox(i)*rappox(i)**14.))**(1./15.)
      
    END IF
    rappox(i)=wmol(i)
    wmol(i)=(xcat0(i)/2.-wmol(i))*2. ! moli di OH IR-like
    wattot(i)=xcat0(i)*0.5
    watmol(i)=wattot(i)-0.5*(wmol(i))
    aden=-totcat(i)*wattot(i)+cost22(i)*wattot(i)*freeox(i)  &
        +0.5*totcat(i)*wmol(i)+0.5*cost22(i)*freeox(i)*wmol(i)
    bden=cost22(i)*freeox(i)+totcat(i)
    ohtrue(i)=aden/bden
    IF (ohtrue(i) < 0.) ohtrue(i)=0.00001*wmol(i)
    
    prot(i)=wmol(i)-ohtrue(i)
    
    xleft(i)=(wattot(i)-0.5*wmol(i))/(wattot(i)+0.5*wmol(i))
    dx(i)=(wattot(i)+0.5*wmol(i))/(wattot(i)-0.5*wmol(i))
    
!      redox(i)=1000.*ohtrue(i)*prot(i)/((100.-poss(i,13))*cost20(i))
!      redox(i)=ohtrue(i)*prot(i)/(100.*cost20(i))
    redox(i)=ohtrue(i)*prot(i)/((somm(i)-0.5*xw(i))*cost20(i))
    
    pipp(i)=2.*redox(i)+ohtrue(i)+prot(i)
    ERR(i)=xcat0(i)-pipp(i)
!   e poi va in iterazione...
    
    crit=DABS(redox(i)-rappox(i))/DABS(redox(i))
    cri(jh)=crit
    
    IF (jh > 1.AND.cri(jh) > cri(jh-1)) THEN
      WRITE (*,*) 'ITERATION NOT ATTAINED FOR THIS COMP !!!...'
      PAUSE
      EXIT
    END IF
    IF (crit < 0.00001) THEN
!         write (*,*) 'CHECK POINT!',kk,crit,i
      EXIT
    END IF
  END DO
  8193  CONTINUE
END DO

30      erstan=0.0
xmis=0.
WRITE (19,*) 'Number of composition is ',ncomp
WRITE (19,*) 'KYMPA option is ',kympa
WRITE (19,*)
WRITE (19,*)'H2Ototwt%,nH2Otot,nH2Omolwt%EXP,nH2OmolEXP,1/T,nfreeOH,&
    nfreeH,nOH-IR,nH2OmolEXP,nH2Omolcalc,nHtot,nHtot_after,DISPERSION,&
    XfreeOHcalc,XfreeHcalc,XH2Omolcalc,YfreeOHcalc,YfreeHcalc,YH2Omolcalc'


DO  i=1,ncomp
  xmis=xmis+redox(i)
  delta(i)=rappox(i)-redox(i)
  ersta(i)=delta(i)**2
  ermean=ermean+ABS(delta(i))
  erstan=erstan+ersta(i)
  WRITE(19,451)poss(i,13),wattot(i),xcat0(i),watmol(i)*18.015,watmol  &
      (i),1./tkelv(i),ohtrue(i),prot(i),wmol(i),rappox(i),redox(i),  &
      pipp(i),ERR(i), ohtrue(i)/((somm(i)-0.5*xw(i))),  &
      prot(i)/((somm(i)-0.5*xw(i))), redox(i)/((somm(i)-0.5*xw(i))),  &
      ohtrue(i)/(ohtrue(i)+prot(i)+redox(i)),  &
      prot(i)/(ohtrue(i)+prot(i)+redox(i)), redox(i)/(ohtrue(i)+prot(i)+redox(i))
  
  451   FORMAT (8(f8.5,x),f9.6,1X,f13.8,1X,14(f12.5,x),1X,f9.3,f11.4,x  &
      ,f8.5,1X,f12.5,1X,f12.5,1X,f12.5,1X,f12.5,1X,f12.5,1X,f12.5,1X,  &
      f12.5,1X,f12.5,1X,f12.5)
  
END DO

ermean=ermean/ncomp
erstan=(erstan/ncomp)**(0.5)
xmis=xmis/ncomp
WRITE(19,*)'Number of compositions ',ncomp
WRITE(19,*)'Mean error  = ',ermean
WRITE(19,*)'Std error ',erstan

STOP
END
