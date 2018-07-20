!module esclec is to save and load parameter files and save and output tooth morphology files
!module coreop2d is the model itself. It includes the following subroutines
!      subroutine ciinicial          specifies the default initial conditions
!      subroutine dime               allocates the matrices for the initial conditions and sets them     
!      subroutine redime             not in use
!      subroutine posar              specifies the position of cells in the initial conditions
!      subroutine calculmarges       calculates the shape of each epithelial cell (this is the position of its margins)
!      subroutine reaccio_difusio    calculates diffusion of all the molecules between cells
!      subroutine diferenciacio      updates cells differentiation values  
!      subroutine empu               calculates pushing between cells resulting from epithelial growth and border growth
!      subroutine stelate            calculates pushing between cells resulting from buoyancy
!      subroutine pushing            calculates repulsion between neighboring cells
!      subroutine pushingnovei       checks if non-neighbors cells get too close and applies repulsion between them (it never happens for seals)
!      subroutine biaixbl            applies BMP4 concentration in the buccal and lingual borders of the tooth
!      subroutine promig             calculates the nucleus traction by the cell borders
!      subroutine actualitza         updates cell positions
!      subroutine afegircel          calculates where cell divisions occur and adds new cells accordingly
!      subroutine perextrems         identifies which of the cells added in afegircel are in the border of the tooth
!      subroutine iteracio           determines the other by which the subroutines are called (invariable)

!***************************************************************************
!***************  MODUL ****************************************************
!***************************************************************************
! gfortran -w -fexceptions -fno-underscoring -fbounds-check source_original.f90 -o TM.e

module coreop2d  

implicit none
public :: iteracio,ciinicial,calculmarges,reaccio_difusio,dime

!coreop2d
real*8, public, allocatable  :: malla(:,:)  					! x,y and z positions of each cell
real*8, public, allocatable  :: marge(:,:,:) 
integer, public, allocatable :: vei(:,:)    					!indices of each cells neighbors
integer, public, allocatable :: knots(:)
integer, public, allocatable :: nveins(:)
real*8, public, allocatable  :: q2d(:,:)    
real*8, public, allocatable  :: q3d(:,:,:)  					! concentrations of the different molecules
real*8, public,allocatable   :: difq3d(:),difq2d(:)
real*8, public, allocatable  :: hmalla(:,:),hvmalla(:,:)
real*8, public, allocatable  :: px(:),py(:),pz(:)

integer, public :: ncels   														!number of cells during the simulations
integer, public :: ncals   
integer, public :: ncz     														!thickness of mesenchyme (number of mesenchyme cell below a specific epithelial cell)
integer, parameter, public :: nvmax=30   
integer, public :: radi
integer, public, parameter :: ng=5,ngg=4						  !number of genes (ng -> gen number + 1)F
integer, public :: temps,npas
real*8, public,  parameter :: la=1.      							!distance between cells in the original conditions

!different model parameters (some not in use)
real*8, public :: ud,us  
real*8, public :: tacre 
real*8, public :: tahor
real*8, public :: acac
real*8, public :: acec
real*8, public :: acaca
real*8, public :: ihac
real*8, public :: ih
real*8, public :: elas
real*8, public :: tadi    
real*8, public :: crema   
real*8, public, parameter :: dmax=2.   
real*8, public :: bip,bia,bil,bib   
real*8, public :: ampl 
real*8, public :: mu 
real*8, public :: tazmax 
real*8, public :: radibi 
real*8, public :: radibii 
real*8, public :: fac
real*8, public :: condme 
real*8, public :: tadif 
!semiconstants
real*8, public :: umelas
integer, public :: maxcels

!for implementation
real*8, parameter ::  delta=0.005D1 , vmin=0.015D1 
integer,public :: nca,icentre,centre,ncils,focus
real*8,public:: x,y,xx,yy
real*8, public :: csu,ssu,csd,ssd,cst,sst,csq,ssq,csc,ssc,css,sss
integer, public :: nnous
integer, public :: nmaa,nmap
integer, public,allocatable :: mmaa(:),mmap(:)

!the typical
integer, public :: i,j,k,ii,jj,kk,iii,jjj,kkk,iiii,jjjj,kkkk,iiiii,jjjjj,kkkkk
real*8, public :: a,b,c,d,e,f,g,h,aa,bb,cc,dd,ee,ff,gg,hh,aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,panic

!for visualization
integer, public :: vlinies,vrender,vmarges,vvec,vvecx,vveck,vex,vn
integer, public :: nc
integer, public :: pin,pina 
integer, public :: submenuid
integer, public :: nivell,kko 

!constants
real*8, public, parameter ::  pii = 31.41592653589793D-1

CONTAINS

subroutine ciinicial
  real*8 ua,ub,uc,ux,uy,uz

  !core values
  ncz=4
  temps=0
  npas=1
  !values for visualization
  vlinies=1
  vrender=1
  vvec=0
  vn=0
  vveck=0
  vvecx=0
  vmarges=0
  vex=0
  pin=1
  pina=1
  nivell=1
  panic=0
  allocate(difq3d(ng))
  allocate(difq2d(ngg))
end subroutine ciinicial


subroutine dime
  integer, allocatable :: cv(:,:)
  real*8 , allocatable :: cmalla(:,:)
  integer iit

  umelas=1-elas

  j=0
  do i=1,radi   ; j=j+i ; end do ; ncals=6*j+1 !calculate the number of ncals
  j=0
  do i=1,radi-1 ; j=j+i ; end do ; ncels=6*j+1 !calculate the number of ncels


  a=pii*0.2D1/0.36D3		!a degree value multiplied with this a gives the equivalent rad value (for the cos and sin values)
  !these values added to a center coordinate, will result in the 6 x and y components of the neigbours (in a hexagon shape)
  csu=dsin(0*a)  ; ssu=dcos(0*a)
  csd=dsin(60*a)  ; ssd=dcos(60*a)
  cst=dsin(120*a) ; sst=dcos(120*a)
  csq=dsin(180*a) ; ssq=dcos(180*a)
  csc=dsin(240*a) ; ssc=dcos(240*a)
  css=dsin(300*a) ; sss=dcos(300*a)


  !allocations
  allocate(cv(ncals,nvmax))
  allocate(cmalla(ncals,3))
  allocate(malla(ncals,3))
  allocate(vei(ncals,nvmax))
  allocate(hmalla(ncals,3))
  allocate(hvmalla(ncals,3))
  allocate(marge(ncals,nvmax,8))
  allocate(knots(ncals))
  allocate(nveins(ncals))
  allocate(q2d(ncals,ngg))
  allocate(q3d(ncals,ncz,ng))
  allocate(mmap(radi))
  allocate(mmaa(radi))

  !matrices for visualisation
  allocate(px(ncals)) ; allocate(py(ncals)) ; allocate(pz(ncals))

  ampl=radi*0.75

  !values that are zero
  vei=0. ; nveins=0. ; malla=0. ; q2d=0. ; q3d=0. ; knots=0 ; hmalla=0. ; hvmalla=0.
  
  !initial values
  malla(1,1)=0. ; malla(1,2)=0. ; malla(1,3)=1.
  nca=1
  nveins=6

al: do icentre=1,ncels
    x=malla(icentre,1) ; y=malla(icentre,2)

    !calculate for each center cell (icentre) the 6 neighbouring coordinates. 
    !Posar then checks if the coordinate is already a cell and if it is already a neighbour. If not, it makes a new cell that is a neighbour of the current cell (icentre)
    xx=x+csu*la ; yy=y+ssu*la ; j=1 ; jj=4 ; call posar
    xx=x+csd*la ; yy=y+ssd*la ; j=2 ; jj=5 ; call posar
    xx=x+cst*la ; yy=y+sst*la ; j=3 ; jj=6 ; call posar
    xx=x+csq*la ; yy=y+ssq*la ; j=4 ; jj=1 ; call posar
    xx=x+csc*la ; yy=y+ssc*la ; j=5 ; jj=2 ; call posar
    xx=x+css*la ; yy=y+sss*la ; j=6 ; jj=3 ; call posar

  do i=2,ncels
    do j=1,nvmax
      if (vei(i,j)>ncels) then 	!if a cell has a neighbour, that is not within the ncels
        vei(i,j)=ncals       		!then replace this neigbhour with the highest cell -> its a marked neighbour
      end if 
    end do
  end do

  do k=1,3
    do i=2,ncels
      do j=1,nvmax-1
        if (vei(i,j)==ncals.and.vei(i,j+1)==ncals) then 	!if a cell has two adjacent marked neighbours as defined above
          do jj=j,nvmax-1
            vei(i,jj)=vei(i,jj+1)			      							!then replace this and all the coming neighbours with the next neighbour --> reduce these two marked neighbours to one
          end do
        end if 
      end do
    end do
  end do

  do i=2,ncels
    k=0
    do j=1,nvmax		!search in all neighbours for marked ones
      if (vei(i,j)==ncals.and.k==0) then ; k=1 ; cycle ; end if				!if its the first marked neighbour of the current cell, its still okay
      if (vei(i,j)==ncals.and.k==1) then ; vei(i,j)=0 ; exit ; end if !if its the second one, delete it and go to the next cell
    end do
  end do

  !trim all the values in malla
  malla=dnint(malla*1D14)*1D-14

  !if then a value is negativ, make it zero
  do i=1,ncels
    do j=1,3
      if (abs(malla(i,j))<1D-14) malla(i,j)=0. 
    end do
  end do

!calcul de distancia original entre nodes

!inversio de forma que els primers son als marges

  cv=vei
  cmalla=malla

  do i=ncels,1,-1
    vei(i,:)=cv(ncels-i+1,:)				!This transformes vei to its inverse (-> the order of cells is inverted)
    malla(i,:)=cmalla(ncels-i+1,:)	!This tranformes malla to its inverse
  end do

  cv=vei
  !Invert the values of the matrix (19->1, 18->2, 17->3,...,1->19) 0 and 37 stay -> do this with cv (not with vei)
  do i=ncels,1,-1
    ii=ncels-i+1
    do jj=1,ncels
      do jjj=1,nvmax
        if (cv(jj,jjj)==i)  vei(jj,jjj)=ii
      end do
    end do
  end do

  call calculmarges 

	nveins=3
	nveins(1)=6
	marge(:,:,4:5)=la

	centre=ncels
	ncils=(radi-1)*6+1
	focus=ncils/2-1
	if (radi==2) focus=3
  mmaa=0 ; mmap=0

  do i=1,radi ; mmap(i)=i ; end do
  ii=0

  do i=ncils/2+1,ncils/2+radi ; ii=ii+1 ; mmaa(ii)=i ; end do

  nmaa=radi ; nmap=radi

call calculmarges		!braucht es das??

q3d=0

end subroutine dime

subroutine posar
al: do i=1,nca 
      if (i==icentre) cycle al
      if (dnint(1000000*malla(i,1))==dnint(1000000*xx).and.dnint(1000000*malla(i,2))==dnint(1000000*yy)) then ; !check if the current coordinates are already present in malla
        do ii=1,nvmax 																																													!check if this point is already a neihbour of the current cell
           if (vei(icentre,ii)==i) then 																																				!if this point is already a neighbour of the current cell, return
               return
           end if   
        end do
      vei(icentre,j)=i ; vei(i,jj)=icentre ; nveins(i)=nveins(i)+1 ; nveins(icentre)=nveins(icentre)+1 ; return !if not, declare it as neighbour and return 
														!(it doesn't have to be put into malla, because the point already exists)
      end if
    end do al
    nveins(icentre)=nveins(icentre)+1 										!the current cell will get a new neighbour
    nca=nca+1																							!generate the ID of the new cell		
    vei(icentre,j)=nca																		!the current cell has the new cell as neighbour
    vei(nca,jj)=icentre																		!the new cell will have the current cell as neighbour
    malla(nca,1)=xx ; malla(nca,2)=yy ; malla(nca,3)=1. 	!set the coordinates of the new cell
    nveins(nca)=nveins(nca)+1															!the new cell gets a neighbour (the current cell)
end subroutine posar
 
subroutine calculmarges
  real*8 cont
  integer kl

  marge(:,:,1:3)=0.

  do i=1,ncels	
    aa=0. ; bb=0. ; cc=0. ; kl=0
    do j=1,nvmax
      if (vei(i,j)/=0.) then																												!whenever they is a neighbour vei(i,j)
        a=0. ; b=0. ; c=0. ; cont=0
        iii=i 
        a=malla(i,1) ; b=malla(i,2) ; c=malla(i,3) ; cont=1													!coordinates of i = (a,b,c)
        ii=vei(i,j) 																																!remember the ID of this neighbour
        if (ii>ncels) then																													!if this neighbour is not within ncels
          do jj=j-1,1,-1
            if (vei(i,jj)/=0) then 																									!look if other neighbours of i are present... 
              if (vei(i,jj)<ncels+1) then ; 																				!...that are within ncels..
                ii=vei(i,jj) ; a=a+malla(ii,1) ; b=b+malla(ii,2) ; c=c+malla(ii,3) 	!if yes, then add the coordinates of this neighbour to i
                cont=cont+1 ; goto 77 ; 
              else ; goto 77 ; 																											!if not (j is the only neighbour within ncels of i), then goto 77
              end if
            end if 
          end do 
          do jj=nvmax,j+1,-1 																												!look also at the other part of the neighbours of i, then the same as above
            if (vei(i,jj)/=0) then  
              if (vei(i,jj)<ncels+1) then ; ii=vei(i,jj) ; a=a+malla(ii,1) ; b=b+malla(ii,2) ; c=c+malla(ii,3) 
              cont=cont+1 ; goto 77 ; 
              else ; goto 77 ; 
              end if
            end if 
          end do 
          goto 77 
        end if
66      if (ii==i) goto 77
        kl=kl+1
        if (kl>100) then
          do jj=j-1,1,-1
            if (vei(i,jj)/=0) then  
              if (vei(i,jj)<ncels+1) then ; ii=vei(i,jj) ; a=a+malla(ii,1) ; b=b+malla(ii,2) ; c=c+malla(ii,3) 
                cont=cont+1 ; goto 77 ; else ; goto 77 ; end if
              end if 
           end do 
           do jj=nvmax,j+1,-1
             if (vei(i,jj)/=0) then  
               if (vei(i,jj)<ncels+1) then ; ii=vei(i,jj) ; a=a+malla(ii,1) ; b=b+malla(ii,2) ; c=c+malla(ii,3) 
               cont=cont+1 ; goto 77 ; else ; goto 77 ; end if
             end if 
           end do 
           goto 77 
        end if
        a=a+malla(ii,1) ; b=b+malla(ii,2) ; c=c+malla(ii,3) ; cont=cont+1
        do jj=1,nvmax ; if (vei(ii,jj)==iii) then ; jjj=jj  ; exit ; end if ; end do  !search the correspondent neighbour entry for this
                     								      																						!relationship in the neighbour
        do jj=jjj+1,nvmax  !comencem la gira
          if (vei(ii,jj)/=0) then 																										!for all other neighbours of the neighbour
            if (vei(ii,jj)>ncels) goto 77 
            iii=ii ; ii=vei(iii,jj) 																									!the new ii is now the neighbour of the neighbour
            goto 66 
          end if
        end do
        do jj=1,jjj-1  !comencem la gira
          if (vei(ii,jj)/=0) then ; if (vei(ii,jj)>ncels) goto 77 ; iii=ii ; ii=vei(iii,jj) ; goto 66 ; end if
        end do
      end if
77    marge(i,j,1)=a/cont ; marge(i,j,2)=b/cont ; marge(i,j,3)=c/cont 
    end do
  end do
end subroutine calculmarges

subroutine reaccio_difusio
  real*8 pes(ncals,nvmax)         !contact area between i and vei(i,j)
  real*8 areap(ncals,nvmax)
  real*8 suma,areasota
  real*8 hq3d(ncals,ncz,ng)
  real*8 hq2d(ncals,ngg)
  real*8 ux,uy,uz,dx,dy,dz,ua,ub,uc
  integer primer 

  hq3d=0.
  hq2d=0.

  do i=1,ncels
    pes(i,:)=0. ; areap(i,:)=0.
ui: do j=1,nvmax 																														!for each neighbour of i
      if (vei(i,j)/=0.) then 
        ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)												!Coordinates of i
        do jj=j+1,nvmax																											!Go to the/a next neighbour of i -> jj
          if (vei(i,jj)/=0.) then	
            pes(i,j)=sqrt((marge(i,j,1)-marge(i,jj,1))**2+(marge(i,j,2)-marge(i,jj,2))**2+(marge(i,j,3)-marge(i,jj,3))**2)
            ux=marge(i,j,1)-ua  ; uy=marge(i,j,2)-ub  ;  uz=marge(i,j,3)-uc !Vector u: i->mij
            dx=marge(i,jj,1)-ua ; dy=marge(i,jj,2)-ub ; dz=marge(i,jj,3)-uc	!Vector d: i->mijj
						!Betrag des Kreuzproduktes von Vektor u und d -> Fläche zwischen u und d
            areap(i,j)=0.05D1*sqrt((uy*dz-uz*dy)**2+(uz*dx-ux*dz)**2+(ux*dy-uy*dx)**2)	!Wieso *0.5?? 
            cycle ui 																												!As soon as an area is calculated, the next neighbour of i is taken into account
          end if
        end do

				!if j is the only neighbour of i, take the first neighbour as second reference
        pes(i,j)=sqrt((marge(i,j,1)-marge(i,1,1))**2+(marge(i,j,2)-marge(i,1,2))**2+(marge(i,j,3)-marge(i,1,3))**2) 
        ux=marge(i,j,1)-ua ; uy=marge(i,j,2)-ub ; uz=marge(i,j,3)-uc 
        dx=marge(i,1,1)-ua ; dy=marge(i,1,2)-ub ; dz=marge(i,1,3)-uc 
        areap(i,j)=0.05D1*sqrt((uy*dz-uz*dy)**2+(uz*dx-ux*dz)**2+(ux*dy-uy*dx)**2)
      end if
    end do ui

    areasota=sum(areap(i,:))
    suma=sum(pes(i,:))+2*areasota !suma: Mantel + 2 * Grundfläche = Oberfläche
		areasota=areasota/suma ; pes(i,:)=pes(i,:)/suma !Anteil der Oberfläche (0-1) von Grundfläche (areasota) und Mantel (pes)

		!Mesenchyme-Mesenchyme Diffusion
    do k=1,4 !for each gene
      do kk=2,ncz-1 !within the mesenchyme
        hq3d(i,kk,k)=hq3d(i,kk,k)+areasota*(q3d(i,kk-1,k)-q3d(i,kk,k)) 			 !within a cell (in z-direction)
        hq3d(i,kk,k)=hq3d(i,kk,k)+areasota*(q3d(i,kk+1,k)-q3d(i,kk,k))
        do j=1,nvmax !for each neighbour
          if (vei(i,j)/=0) then 
            ii=vei(i,j)
            if (ii==ncals) then !if the neighbour is a border cell 
              hq3d(i,kk,k)=hq3d(i,kk,k)+pes(i,j)*(-q3d(i,kk,k)*0.044D1)      !sink     
            else !if not
              hq3d(i,kk,k)=hq3d(i,kk,k)+pes(i,j)*(q3d(ii,kk,k)-q3d(i,kk,k))  !normal diffusion between neighbours in horizontal direction
            end if
          end if
        end do
      end do

			!Diffusion in the last mesenchyme cell (in vertical direction)
      hq3d(i,ncz,k)=areasota*(-q3d(i,ncz,k)*0.044D1) 												!sink in the last mesenchyme cell
      hq3d(i,ncz,k)= hq3d(i,ncz,k)+areasota*(q3d(i,ncz-1,k)-q3d(i,ncz,k)) 	!and diffusion into the upper mesenchyme cell
      do j=1,nvmax
        if (vei(i,j)/=0) then 
          ii=vei(i,j)
          if (ii==ncals) then !and in vertical direction
            hq3d(i,ncz,k)=hq3d(i,ncz,k)+pes(i,j)*(-q3d(i,ncz,k)*0.044D1)   	!sink if the neighbour is a border cell  
          else
            hq3d(i,ncz,k)=hq3d(i,ncz,k)+pes(i,j)*(q3d(ii,ncz,k)-q3d(i,ncz,k)) 
          end if
        end if
      end do
    end do

    pes(i,:)=pes(i,:)*suma ; areasota=areasota*suma 
		
		!Epithelial Diffusion
		suma=suma-areasota !The upper area is not part of the diffusion area
		pes(i,:)=pes(i,:)/suma; areasota=areasota/suma
    do k=1,4 !ng 
      hq3d(i,1,k)=areasota*(q3d(i,2,k)-q3d(i,1,k)) !Mesenchyme->Epithel Diffusion (vertical)
      do j=1,nvmax
        if (vei(i,j)/=0) then 
          ii=vei(i,j)
          if (ii==ncals) then
            hq3d(i,1,k)=hq3d(i,1,k)+pes(i,j)*(-q3d(i,1,k)*0.044D1)   !horizontal Diffusion   
          else
            hq3d(i,1,k)=hq3d(i,1,k)+pes(i,j)*(q3d(ii,1,k)-q3d(i,1,k)) 
          end if
        end if
      end do
    end do
  end do

	!Update the q3d Matrix. To reduce the effect: * delta and dependent on diffusion-coefficient difq3d()
  do k=1,4 ! ng 
    q3d(:,:,k)=q3d(:,:,k)+delta*difq3d(k)*hq3d(:,:,k)
  end do

  !REACCIO -> happens only in epithelial cells
  hq3d=0.
  do i=1,ncels
    if (q3d(i,1,1)>1) then 			!if the activator concentration in the epithelial cell is high enough..
      if (i>=ncils) knots(i)=1 	!..and if it is in the centre (within ncils), then it (becomes) a knot cell
    end if
		
		!Reaction and Degradation of Gen1 (Activator)
    a=acac*q3d(i,1,1)-q3d(i,1,4)
    if (a<0) a=0.
    hq3d(i,1,1)=a/(1+ihac*q3d(i,1,2))-mu*q3d(i,1,1)

    if (q2d(i,1)>us) then														!if the differentiation state is higher than the first threshold
      hq3d(i,1,2)=q3d(i,1,1)*q2d(i,1)-mu*q3d(i,1,2)	!then inhibitor is produced (proportionally to diff-state of activator)
    else
      if (knots(i)==1) then													!if the diff-state is not high enough, but the cell is a EK-cell
        hq3d(i,1,2)=q3d(i,1,1)-mu*q3d(i,1,2)				!then inhibitor is produced proportionally to activator concentration
      end if
    end if

    if (q2d(i,1)>ud) then									!if the diff-state  is higher then the second threshold
      a=ih*q2d(i,1)-mu*q3d(i,1,3)					!then Sec1 is produced, proportionally to diff-state and parameter ih
      if(a<0.) a=0.							
      hq3d(i,1,3)=a
    else
      if (knots(i)>ud) then								!if the diff-state is not high enough but the cell is already a EK-cell
        a=ih-mu*q3d(i,1,3)								!then Sec1 is produced, proportionally to parameter ih
        if(a<0.) a=0.
        hq3d(i,1,3)=a
      end if
    end if

    a=acec*q3d(i,1,1)-mu*q3d(i,1,4)-difq3d(ng)*q3d(i,1,3)	
		!Sec2 concentration is proportional to how strong Act acts minus those who are already Sec1
    if(a<0.) a=0.
    hq3d(i,1,4)=a
  end do

	!if Act or Inh concentration is too high, there is an error message
	if (maxval(abs(hq3d(:,1,1:2)))>1D100) then ; 
		print *,"PANIC OVERFLOW" ; 
		panic=1 ; return ; 
	end if

	!Update the q3d Matrix with new values
  do i=1,4 					
    q3d(:,1,i)=q3d(:,1,i)+delta*hq3d(:,1,i)
  end do

  where(q3d<0.) q3d=0.  !remove negativ values

end subroutine reaccio_difusio

subroutine diferenciacio
  do i=1,ncels 
    q2d(i,1)=q2d(i,1)+tadif*(q3d(i,1,3))		!increase the diff-state of each cell by Parameter*Sec1-Concentration
    if (q2d(i,1)>1.) q2d(i,1)=1.						!the diff-state can be maximally 1
  end do
end subroutine diferenciacio

!Epithelial Proliferation
!calculates pushing between cells resulting from epithelial growth and border growth
subroutine empu
  real*8 ux,uy,uz,uux,uuy,uuz,ua,ub,uc,uuuz,uaa,ubb,uuux,uuuy,duux,duuy
  real*8 persu(nvmax,3)
  integer :: indexx(nvmax)

  hvmalla=0.
  hmalla=0

  do i=ncils,ncels					!for all cells in the very center 
    if (knots(i)==1) cycle	!and only for non-EK cells
    ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
    persu=0.
    aa=0 ; bb=0 ; cc=0
    do j=1,nvmax
      k=vei(i,j) 
      if (k==0.or.k>ncels) cycle				!only if the neighbour is not a border cell
      b=uc-malla(k,3) 									!b = difference in z-direction between i and k
      if (b<-1D-4) then									!if the neighbour is a certain amount higher than i
        uux=ua-malla(k,1)      ; uuy=ub-malla(k,2)      ; uuz=uc-malla(k,3) 
        d=sqrt(uux**2+uuy**2+uuz**2) 		!Distance between the cell centres
        d=1/d
        aa=aa-uux*d ; bb=bb-uuy*d ; cc=cc-uuz*d !how strong the neighbours (in total) deviate from i -> inhomogeneity
      end if
    end do
    d=sqrt(aa**2+bb**2+cc**2)						!d = measure for total deviation of all neighbours of i
    if (d>0) then
      d=tacre/d
      a=1-q2d(i,1) ; if (a<0) a=0.
      d=d*a
			!d is higher the more inhomogen the distribution of cells, the less they are differentiated and the higher the higher tacre is
      hmalla(i,1)=aa*d !-> cells drift away from each other
      hmalla(i,2)=bb*d
      hmalla(i,3)=cc*d
    end if
  end do

  do i=1,ncils-1																	!within ncels, but not the centre
    aa=0. ; bb=0. ; a=-0.3 ; b=0. ; c=0. 
    ua=malla(i,1) ; ub=malla(i,2)
    do j=1,nvmax
      k=vei(i,j)
      if (k<1.or.k>ncels) cycle										!z-difference between i and k does not matter as above
      if (k>ncils-1) then 												!when the neighbour is a centre cell
        uux=ua-malla(k,1)   ; uuy=ub-malla(k,2) 	!z-difference between i and k does not matter as above
        d=sqrt(uux**2+uuy**2)
        if (d>0) then
          c=acos(uux/d)														!angle between x-axis and d
          if (uuy<0) c=2*pii-c 										!we want the left (in negative x-direction) angle
        end if
      else 																				!when the neighbour cell is within ncels, but not in the centre (ncils) (same as i)
        uux=ua-malla(k,1)      ; uuy=ub-malla(k,2)
        d=sqrt(uux**2+uuy**2)
        if (d>0) then
        if (a==-0.3) then													!if a is not already calculated
          a=acos(uux/d)														!angle between x-axis and d
          if (uuy<0) a=2*pii-a										!we want the left angle
          if (d>0) then
            dd=1/d
            uuux=-uuy*dd           ; uuuy=uux*dd	!uuux = -uuy/d		uuuy = -uux/d
            uaa=acos(uuux)												!=acos(-uuy/d)
            if (uuuy<0) uaa=2*pii-uaa
          end if
        else																			!if a is already calculated
          b=acos(uux/d)														!make the same as with a but with b
          if (uuy<0) b=2*pii-b
          if (d>0) then
            dd=1/d
            duux=-uuy*dd           ; duuy=uux*dd
            ubb=acos(duux)
            if (duuy<0) ubb=2*pii-ubb
          end if
        end if  
        end if
      end if
    end do

      if (a<b) then ; d=a ; a=b ; b=d ; end if
      if (c<a.and.c>b) then 
        if (uaa<a.and.uaa>b) then ; uuux=-uuux ; uuuy=-uuuy ; end if ! es a la banda de dins i aleshores l'hem d'invertir
        if (ubb<a.and.ubb>b) then ; duux=-duux ; duuy=-duuy ; end if
      else
        if (uaa>a.or.uaa<b) then ; uuux=-uuux ; uuuy=-uuuy ; end if ! es a la banda de dins i aleshores l'hem d'invertir
        if (ubb>a.or.ubb<b) then ; duux=-duux ; duuy=-duuy ; end if 
      end if    
      aa=-uuux-duux ; bb=-uuuy-duuy  

      !ara mirem que sigui cap a fora de la dent a lo cutre
      a=ua+aa ; b=ub+bb
      c=ua-aa ; d=ub-bb
      dd=sqrt(a**2+b**2)
      ddd=sqrt(c**2+d**2)
      if (ddd>dd) then ; aa=-aa ; bb=-bb ; end if

    ! ames tenim la traccio cap abaix deguda a l'adhesio al mesenkima    
		!We get the result due to adehesion on the mesenchyme
    d=sqrt(aa**2+bb**2)
    if (d>0) then         
      d=(d+tahor*q3d(i,1,3))/d
      aa=aa*d  
      bb=bb*d
    end if
    cc=tazmax
    d=sqrt(aa**2+bb**2+cc**2)
    if (d>0) then
      d=tacre/d
      a=1-q2d(i,1) ; if (a<0) a=0.
      d=d*a                           !aixo sembla estar repe i malament
      hmalla(i,1)=aa*d								!This overwrites the hmalle from above??!
      hmalla(i,2)=bb*d
      hmalla(i,3)=cc*d
    end if
  end do

  hvmalla=hmalla

end subroutine empu

!calculates pushing between cells resulting from buoyancy
subroutine stelate
  real*8 ax,ay

  do i=1,ncels													!for all cells
    ax=hmalla(i,1) ; ay=hmalla(i,2) 
    d=sqrt(ax**2+ay**2)									!d= distance to origin in 2D
    if (d/=0) then
      c=hmalla(i,3)
      if (d>0) then											!is this necessary??
        a=sqrt(ax**2+ay**2+c**2) 				!a=distance to origin in 3D
				a=-c/a	
        ax=ax*a ; ay=ay*a
        dd=sqrt(ax**2+ay**2+d**2) ; dd=difq2d(2)*q3d(i,1,3)/dd	!difq2d(d)=buoyancy
        if (dd>0) then
        a=1-q2d(i,1) ; if (a<0) a=0.
        ax=ax*dd*a ; ay=ay*dd*a ; d=d*dd*a
        hmalla(i,1)=hmalla(i,1)-ax ; hmalla(i,2)=hmalla(i,2)-ay ;hmalla(i,3)=hmalla(i,3)-d     
        end if
      end if
    end if
  end do

end subroutine stelate

subroutine pushing
  real*8 ux,uy,uz,ua,ub,uc,hx,hy,hz,dd,d,uux,uuy,sux,suy,dr,rd
  real*8 persu(nvmax,3)
  integer :: indexx(nvmax)

  do i=1,ncels
    ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
    persu=0.
    do j=1,nvmax
      k=vei(i,j)
      if (k>0.and.k<ncels+1) then	!if the neighbour is also in ncels
        ux=malla(k,1)-ua ; uy=malla(k,2)-ub ; uz=malla(k,3)-uc 
        if (abs(ux)<1D-15) ux=0.	!if the neighbour is very close -> maybe rounding error -> give them the same value
        if (abs(uy)<1D-15) uy=0.
        if (abs(uz)<1D-15) uz=0.
        dr=sqrt(ux**2+uy**2+uz**2) 	!distance between i and j
        rd=marge(i,j,5)							!distance between i and j in x/y plane ("theoretical distance")
        if (dr<1D-8) dr=0.
        if (rd<1D-8) rd=0.
        if (knots(i)==1.and.knots(k)==1) then !if both cells are EK-cells
          d=dr-rd 									!deviation from theoretical distance
          dr=d/dr 
          persu(j,1)=ux*dr ; persu(j,2)=uy*dr ; persu(j,3)=uz*dr
        else
          if (dr<rd) then
            d=dr-rd 
            dr=d/dr 
            persu(j,1)=ux*dr ; persu(j,2)=uy*dr ; persu(j,3)=uz*dr
          else
            if (i>ncils-1) then			!if i is within ncils but i and/or k are not EK-cells, then take the default parameter 
              persu(j,1)=ux*crema ; persu(j,2)=uy*crema ; persu(j,3)=uz*crema 
            end if
          end if
        end if
      end if
    end do

    c=elas 
    if (c>1) c=1
    a=0. ; do j=1,nvmax ; a=a+persu(j,1) ; end do ;
		hmalla(i,1)=hmalla(i,1)+a*c
    a=0. ; do j=1,nvmax ; a=a+persu(j,2) ; end do ;
		hmalla(i,2)=hmalla(i,2)+a*c
    a=0. ; do j=1,nvmax ; a=a+persu(j,3) ; end do ;
		hmalla(i,3)=hmalla(i,3)+a*c

  end do

end subroutine pushing

subroutine pushingnovei
  real*8 ux,uy,uz,ua,ub,uc,hx,hy,hz,dd,d,uux,uuy
  real*8, allocatable :: persu(:,:),cpersu(:,:)  
  integer, allocatable :: indexa(:)
  integer :: conta,espai,espaia

  espai=20

  allocate(persu(espai,3))
  allocate(indexa(espai))

  do i=1,ncels																							!for all cells
    ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
    persu=0. ; conta=0
gg: do ii=1,ncels
      if (ii==i) cycle 
      do j=1,nvmax ; if (vei(i,j)==ii) cycle gg ; end do    !if i and ii are neighbours or the same, do nothing     
      ux=malla(ii,1)-ua 
      if (ux>0.14D1) cycle																	!if i and ii are too far away from each other, do nothing
      uy=malla(ii,2)-ub 
      if (uy>0.14D1) cycle
      uz=malla(ii,3)-uc
      if (uz>0.14D1) cycle
      if (abs(ux)<1D-15) ux=0.															!if they are very close, treat them as if they had the same position (-> rounding?)
      if (abs(uy)<1D-15) uy=0.
      if (abs(uz)<1D-15) uz=0.
      d=sqrt(ux**2+uy**2+uz**2)															!d = distance between i and ii
      if (d<0.14D1) then ! ARBITRARY?? RZ										!if they are close enough
        conta=conta+1																				!conta counts how many cells (ii) are close enough to i									
        if (conta>espai) then 
          espaia=espai
          espai=espai+20
          allocate(cpersu(espai,3)) ; cpersu=0.							!make a new matrix that can hold all values 
          cpersu(1:espaia,:)=persu
          deallocate(persu) ; deallocate(indexa)
          allocate(persu(espai,3)) ; allocate(indexa(espai))
          persu=cpersu
          deallocate(cpersu)
        end if 
        dd=1/(d+10D-1)**8 ; d=dd/d ; d=aint(d*1D8)*1D-8			!the smaller d, the even more bigger the force -> **8
        persu(conta,1)=-ux*d ; persu(conta,2)=-uy*d  ; persu(conta,3)=-uz*d !persu: measure of "compression"
      end if
    end do gg

    !versio rapida sense ordenar (possible biaixos per floats)
    c=elas 				!Repulsion between tissues
    if (c>1) c=1	!c is not in use afterwards...
    a=0. ; do j=1,espai ; a=a+persu(j,1) ; end do ;
		hmalla(i,1)=hmalla(i,1)+a*elas
    a=0. ; do j=1,espai ; a=a+persu(j,2) ; end do ;
		hmalla(i,2)=hmalla(i,2)+a*elas
    a=0. ; do j=1,espai ; a=a+persu(j,3) ; end do ;
		hmalla(i,3)=hmalla(i,3)+a*elas
  end do

end subroutine pushingnovei

subroutine biaixbl ! RZ: THIS SEEMS A LITTLE CRAPPY: DISCRETE?!
  do i=1,ncils-1	!for all cells in the very centre (-> primary EK?)
    if (malla(i,2)<-tadi) then !if they are too far away from the long-axis, their Act-Concentration is set to bil/bib
      q3d(i,1,1)=bil
    else
      if (malla(i,2)>tadi) then
        q3d(i,1,1)=bib
      end if
    end if
  end do
end subroutine

!calculates the nucleus traction by the cell borders
subroutine promig  
  real*8 n
  real*8,allocatable :: pmalla(:,:),testdumb(:)
  integer ncalsi

	ncalsi=ncals
	if(allocated(pmalla)) deallocate(pmalla)
	allocate(pmalla(ncals,3)) ! RZ
	pmalla=0d0 ! RZ
  pmalla=malla

  do i=ncils,ncels							!for all centre cells (ncils)
    if (q2d(i,1)==1) cycle			!only if they are not already fully differentiated
    a=0. ; b=0. ; c=0. ; n=0
    do j=1,nvmax
      k=vei(i,j)
      if (k/=0.and.k<ncels+1) then
        a=a+malla(k,1) ; b=b+malla(k,2) ; c=c+malla(k,3)	!sum up the coordinates of all neighbours that are within ncels
        n=n+1
      end if
    end do
    n=1/n
    a=a*n ; b=b*n ; c=c*n 																!and divide it by the number of counted cells -> average position of the neighbours
    a=a-malla(i,1)																				!Difference to the average of neighbours
    b=b-malla(i,2)
    c=c-malla(i,3)
    pmalla(i,1)=malla(i,1)+delta*radibi*a									!radibi: parameter of nuclear traction
    pmalla(i,2)=malla(i,2)+delta*radibi*b
    if (knots(i)==0) then  ! WHY ONLY IN THE EXTREME CASE ? RZ
      a=1-q2d(i,1) ; if (a<0) a=0.												!only if the cell is not a EK cell, the z-position is affected by nuclear traction
      pmalla(i,3)=malla(i,3)+delta*radibi*c*a
    end if
  end do

  do i=1,ncils-1																					!for all cells at the border of ncels (same as with centre cells above)
    if (q2d(i,1)==1) cycle															
    a=0. ; b=0. ; c=0. ; n=0
    do j=1,nvmax
      k=vei(i,j)
      if (k>0.and.k<ncils.and.k<ncels+1) then
        a=a+malla(k,1) ; b=b+malla(k,2) ; c=c+malla(k,3)
        n=n+1
      end if
    end do
    n=1/n
    a=a*n ; b=b*n ; c=c*n 
    a=a-malla(i,1)
    b=b-malla(i,2)
    c=c-malla(i,3)
    pmalla(i,1)=malla(i,1)+delta*radibi*a
    pmalla(i,2)=malla(i,2)+delta*radibi*b
    if (knots(i)==0) then
      a=1-q2d(i,1) ; if (a<0) a=0.
      pmalla(i,3)=malla(i,3)+delta*radibi*c*a
    end if
  end do
  malla=pmalla
end subroutine promig

!updates cell-positions
subroutine actualitza

	!determinem els extrems
  do i=1,ncils-1																							!for all border cells of ncels
    if (abs(malla(i,2))<radibii) then													!radibii = 0.8 -> AP-bias applies only for cells near the x-axis
      if (malla(i,1)>0) then ; hmalla(i,1)=hmalla(i,1)*bia ; 
        hmalla(i,3)=hmalla(i,3)*fac ; 
      end if
      if (malla(i,1)<0) then ; hmalla(i,1)=hmalla(i,1)*bip ; 
        hmalla(i,3)=hmalla(i,3)*fac ; 
      end if
    end if    
  end do

  do i=1,ncels
    if (hmalla(i,3)<0) hmalla(i,3)=0. !there cannot be any force in negativ z-direction due to the pressure of the stelate
  end do

  do i=1,ncels
    if (knots(i)==1) hmalla(i,3)=0.		!there no force at all in z-direction, if the cell is an EK-cell
  end do

  do i=1,ncels
      malla(i,:)=malla(i,:)+delta*hmalla(i,:)  !apply the force vector on the positions (malla)
  end do

end subroutine actualitza

subroutine afegircel

integer primer,segon 

real*8,  allocatable :: cmalla(:,:)  ! les posicions dels nodes x,y,z
integer, allocatable :: cvei(:,:),ccvei(:,:)
integer, allocatable :: cnveins(:)
integer, allocatable :: cknots(:)
real*8,  allocatable :: cq2d(:,:),cmarge(:,:,:),ccmarge(:,:,:)    ! quantitats que son 2d
real*8,  allocatable :: cq3d(:,:,:)    ! quantitats que son 3d act,inh,fgf,ect,p
integer,  allocatable :: scvei(:,:)
integer,  allocatable :: scnveins(:)
real*8 ,  allocatable :: scmalla(:,:)
real*8 ,  allocatable :: scq3d(:,:,:)
real*8 ,  allocatable :: scq2d(:,:)
real*8 ,  allocatable :: scmarge(:,:,:)
integer ,  allocatable :: scknots(:)

integer             up(nvmax),do(nvmax)
real*8 pup(nvmax),pdo(nvmax)
real*8 ua,ub,uc,ux,uy,uz,dx,dy
integer             ord(nvmax)

integer nousnodes(ncels*nvmax,2),externsa(ncels*nvmax)	!nousnodes(ID of new cell, ID of the two mother cells)
integer pillats(nvmax),cpillats(nvmax)
integer nncels,ancels
integer cj,ini,fi,sjj,ji,ij,ijj,jji

	!"nous" = "new"
  nnous=0
  nousnodes=0
  externsa=0

	!first, identify and name new nodes and reset the arrays malla and vei
  do i=1,ncels
    primer=0. ; kkk=0 ; ji=0
    do j=1,nvmax
      if (vei(i,j)>ncels) then ; ji=1 ; exit ; end if				!count (-> ji) how many cells have a border cell (>ncels) as neighbour
    end do
    do j=1,nvmax
      k=vei(i,j)
      ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
			!if k is a neighbour that we do not have looked at it already and that is also within ncels
      if (k/=0.and.k>i.and.k<=ncels) then										
        ux=malla(k,1) ; uy=malla(k,2) ; uz=malla(k,3)
        ux=ux-ua ; uy=uy-ub ; uz=uz-uc
        a=sqrt(ux**2+uy**2+uz**2)														!a = distance between i and k
        a=dnint(a*1D9)*1D-9																	!round a
        if (a>dmax) then  																	!if the distance is big enough (>2) then we add a new cell(=node)
          nnous=nnous+1 ; nousnodes(nnous,1)=i ; nousnodes(nnous,2)=k 
          if (i<ncils.and.k<ncils) then ; externsa(nnous)=1 ; end if	!externsa(x)=1 if the new cell Nr.x is at the border of ncels (not ncils)
        end if
      end if    
    end do 
  end do

  if (nnous>0) then
    do i=1,ncals
      do j=1,nvmax
        if (vei(i,j)==0) then
          do jj=j,nvmax-1
            vei(i,jj)=vei(i,jj+1)		!delete all entries of vei that are 0
          end do
        end if
      end do
    end do

    nncels=ncals+nnous

		!allocate new matrices that are big enough to hold also the new cells
    allocate(cmalla(nncels,3))   ; allocate(cvei(nncels,nvmax))
    allocate(cnveins(nncels))    ; allocate(cq2d(nncels,ngg))   ; allocate(cq3d(nncels,ncz,ng))
    allocate(ccvei(nnous,nvmax)) ; allocate(cknots(nncels))     ; allocate(cmarge(nncels,nvmax,8))

    cmalla=0. ; cvei=0. ; cnveins=0. ; cq2d=0. ; cq3d=0. ; cknots=0 ; cmarge=0.

		!shift all entries in the matrices after the ncels entries by the amount of new cells backwards
    do i=nncels,ncels+nnous+1,-1
      cmalla(i,:)=malla(i-nnous,:)     ; cvei(i,:)=vei(i-nnous,:)
      cnveins(i)=nveins(i-nnous)       ; cq2d(i,:)=q2d(i-nnous,:)
      cq3d(i,:,:)=q3d(i-nnous,:,:)     ; cknots(i)=knots(i-nnous)
      cmarge(i,:,4:8)=marge(i-nnous,:,4:8)
    end do

		!the first entries up to ncels are kept
    cmalla(1:ncels,:)=malla(1:ncels,:)     ; cvei(1:ncels,:)=vei(1:ncels,:) ; cnveins(1:ncels)=nveins(1:ncels)       
    cq2d(1:ncels,:)=q2d(1:ncels,:)
    cq3d(1:ncels,:,:)=q3d(1:ncels,:,:)     ; cknots(1:ncels)=knots(1:ncels) ; cmarge(1:ncels,:,:)=marge(1:ncels,:,:)   

    !update the "marker" of the border cells from the old ncals value to the new one
		do i=1,ncels+nnous ; do j=1,nvmax ; if (cvei(i,j)>ncels) cvei(i,j)=nncels ; end do ; end do

		!set all neighbours of the border cells (beyond ncels) to 0
    cvei(ncels+nnous+1:nncels,:)=0

    do i=1,nnous																									!for all new cells:
      ii=nousnodes(i,1) ; kk=nousnodes(i,2) ; jj=ncels+i					!ii & kk: IDs of mother cells, jj: new ID of the new cell
      cvei(jj,:)=0 ; cvei(jj,1)=ii ; cvei(jj,2)=kk								!Set the mother cells as the first to neighbours of the new cell

      a=cmalla(ii,1)+cmalla(kk,1) ; b=cmalla(ii,2)+cmalla(kk,2)		!a & b: sum of x and y positions of the mother cells
      d=sqrt(a**2+b**2)																						
      a=a/d ; b=b/d
      d=d/0.20D1
      d=dnint(d*1D10)*1D-10
      a=d*a ; b=d*b																								!a and b are divided in half
      cmalla(jj,1)=a ; cmalla(jj,2)=b															!the position of the new cell is in the middle of its mother cells
      cmalla(jj,3)=(malla(ii,3)+malla(kk,3))*0.05D1

			!The concentrations in the new cell are the mean of the mother cells
      cq3d(jj,:,:)=(cq3d(ii,:,:)+cq3d(kk,:,:))*0.05D1 ; cq2d(jj,:)=(cq2d(ii,:)+cq2d(kk,:))*0.05D1

			!Mothercell1 is not a neighbour anymore of mothercell2 -> replace it with the new cell
      do j=1,nvmax ; if (cvei(ii,j)==kk) then ; cvei(ii,j)=jj ; exit ; end if ; end do
      do j=1,nvmax ; if (cvei(kk,j)==ii) then ; cvei(kk,j)=jj ; exit ; end if ; end do
    
      pillats=0
      
			!look which mothercell, that has no border cells at j+1, has to follow
      do j=1,nvmax ; if (vei(ii,j)==kk) then ; jjj=j ; exit ; end if ; end do
			!jjj: mothercell2 is the "jjj"th neighbour of mothercell1 (before the placement of the new cell)
      kkk=0

      do jjjj=jjj+1,nvmax
        if (vei(ii,jjjj)>0) then						!if the next neighbour of mothercell1 
          if (vei(ii,jjjj)<ncels+1) then 		!... is whithin ncels
            ini=ii ; fi=kk ; kkk=1 ; exit	
          else 															!... is not within ncels
            ini=kk ; fi=ii ; kkk=1 ; exit 
          end if
        end if       
      end do

      if (kkk==0) then											!if mothercell1 has no "higher" neighbour, look for lower ones
        do jjjj=1,jjj-1
          if (vei(ii,jjjj)>0) then
            if (vei(ii,jjjj)<ncels+1) then 
              ini=ii ; fi=kk ; exit
            else 
              ini=kk ; fi=ii ; exit 
            end if       
          end if
        end do
      end if

      iii=ini
      cj=1 ; pillats(cj)=iii

      do j=1,nvmax
        if (cvei(iii,j)==jj) then ; jjj=j ; exit ; end if		!the new cell is the jjjth neighbour of a mothercell
      end do

			!Let's reschedule the line next to it
      kkk=0
      do j=jjj+1,nvmax 
				jji=cvei(iii,j) 			!look at the neighbour after the new cell that is within ncels
				if (jji/=0.and.jji<ncels+nnous+1) then 
					iiii=jji 						!iiii: ID of this neighbour
      		kkk=1 ; exit  
				end if 
			end do

      if (kkk==0) then !if we did not find such a "higher" neighbour, look also at the lower ones
        do j=1,jjj-1 
					jji=cvei(iii,j) 
					if (jji/=0.and.jji<ncels+nnous+1) then 
						iiii=jji 
        		kkk=1 ; exit  
					end if 
				end do
      end if

      cj=cj+1 

			if (cj>nvmax) then ; panic=1 ; return ; end if 

			pillats(cj)=iiii

      do j=1,nvmax 
				if (cvei(iiii,j)==iii) then 
					jjjj=j ; exit 							!iiii (neighbour within ncels of newcell) is the jjjjth neighbour of a mothercell
				end if 
			end do

88    iii=iiii ; jjj=jjjj ; kkk=0

      do j=jjj+1,nvmax 
				jji=cvei(iii,j) 
				if (jji/=0.and.jji<ncels+nnous+1) then 	
					iiii=jji 										!iiii is the ID of the neighbour within ncels of a mothercell
      		kkk=1 ; exit  
				end if 
			end do

      if (kkk==0) then 								!if we did not find such a "higher" neighbour, look also at the lower ones
        do j=1,jjj-1 
					jji=cvei(iii,j) 
					if (jji/=0.and.jji<ncels+nnous+1) then 
						iiii=jji 
        		kkk=1 ; exit  
					end if 
				end do
      end if

      cj=cj+1 ; if (cj>nvmax) then ; panic=1 ; return ; end if 

			pillats(cj)=iiii

      do j=1,nvmax 
				if (cvei(iiii,j)==iii) then 
					jjjj=j ; exit 									!the mothercell is the jjjjth neighbour of another neighbour of the mothercell
				end if 
			end do

      if (iiii==fi) then !equinox					!if the neighbour of the mothercell is the other mothercell
        kkk=0
        do kkkk=jjjj+1,nvmax
          if (cvei(iiii,kkkk)/=0.and.kkk==1) then  
            if (cvei(iiii,kkkk)>ncels+nnous) then
              iiii=ini ; kkk=2 ; cj=cj+1 ; exit
            else              
              sjj=kkkk ; kkk=2 ; exit 
            end if
          end if
          if (cvei(iiii,kkkk)/=0.and.kkk==0) then; kkk=1 ; sjj=kkkk ;  end if
        end do

        if (kkk<2) then
          do kkkk=1,jjjj-1
            if (cvei(iiii,kkkk)/=0.and.kkk==1) then
              if (cvei(iiii,kkkk)/=0.and.kkk==1) then  
                if (cvei(iiii,kkkk)>ncels+nnous) then
                  iiii=ini ; cj=cj+1 ; exit
                else              
                  sjj=kkkk ; exit 
                end if
              end if
            end if
            if (cvei(iiii,kkkk)/=0.and.kkk==0) then; kkk=1 ; sjj=kkkk ; end if
          end do
        end if
        jjjj=sjj-1
      end if

      !we have made the whole turn
      if (iiii==ini) then
        cpillats=0
        if (cj>nvmax) then ; panic=1 ; return ; end if ;
        pillats(cj)=0 ; cj=cj-1
        do jjj=1,cj
          cpillats(cj-jjj+1)=pillats(jjj) 
        end do
        pillats=cpillats
				
				!now we look which cells can really exist
        jjj=0
        if (cj>nvmax) then ; panic=1 ; return ; end if ;
        do kkk=1,cj
          kkkk=pillats(kkk)
          if (kkkk>ncels.and.kkkk<=ncels+nnous) jjj=jjj+1
        end do
        if (jjj==0) then  !if we dont have new nodes at the side -> crossing is not possible
          ccvei(i,:)=pillats  
          do j=1,nvmax
            k=ccvei(i,j)
            if (k/=0) then
              ii=nousnodes(i,1) ; iiii=nousnodes(i,2)
              if (k/=ii.and.k/=iiii.and.k<ncels+nnous+1) then ! es un als dels que em de posar conexio
                kkkk=0
uu:             do kk=1,nvmax
uuu:              do kkk=1,nvmax
                    if (cvei(k,kk)/=0.and.cvei(k,kk)==pillats(kkk).and.kkkk==1) then ; ji=kk ; exit uu ; end if
                    if (cvei(k,kk)/=0.and.cvei(k,kk)==pillats(kkk).and.kkkk==0) then ; kkkk=1 ; ij=kk ; exit uuu ; end if
                  end do uuu
                end do uu
                !connect between ij and ji
                if (ji-ij==1) then
                  do kk=nvmax,ji+1,-1 ; cvei(k,kk)=cvei(k,kk-1) ; cmarge(k,kk,4:8)=cmarge(k,kk-1,4:8) ; end do 
                  cvei(k,ji)=jj
                else
                  cvei(k,ji+1)=jj ; 
                end if
              end if
            end if
          end do     
        else                       
          ccvei(i,:)=0
          ccvei(i,1)=ini
          kkkk=1
          if (cj>nvmax) then ; panic=1 ; return ; end if ;
rtt:      do kkk=1,cj
            jjjj=pillats(kkk)
            if (jjjj==fi) then 
              kkkk=kkkk+1 ; ccvei(i,kkkk)=fi 
            else
              if (jjjj>ncels) then
                kkkk=kkkk+1 ; ccvei(i,kkkk)=jjjj
              end if
            end if
          end do rtt
goto 899
          do j=1,nvmax
            k=ccvei(i,j)
            if (k/=0) then
              ii=nousnodes(i,1) ; iiii=nousnodes(i,2)
              if (k/=ii.and.k/=iiii.and.k<ncels+1) then ! es un als dels que em de posar conexio
                kkkk=0
uuuu:           do kk=1,nvmax
uuuuu:            do kkk=1,nvmax
                    if (cvei(k,kk)/=0.and.cvei(k,kk)==pillats(kkk).and.kkkk==1) then ; ji=kk ; exit uuuu ; end if
                    if (cvei(k,kk)/=0.and.cvei(k,kk)==pillats(kkk).and.kkkk==0) then ; kkkk=1 ; ij=kk ; exit uuuuu ; end if
                  end do uuuuu
                end do uuuu
                !connect between ij and ji
                if (ji-ij==1) then
                  do kk=nvmax,ji+1,-1 ; cvei(k,kk)=cvei(k,kk-1) ; cmarge(k,kk,4:8)=cmarge(k,kk-1,4:8) ; end do
                  cvei(k,ji)=jj
                else
                  cvei(k,ji+1)=jj
                end if
              end if
            end if
          end do
899 continue     
        end if

				!add connections to external nodes
        ii=nousnodes(i,1) ; kk=nousnodes(i,2)
        kkk=0 ; jjj=0
        do j=1,nvmax
          if (cvei(ii,j)>ncels+nnous) then ; kkk=1 ; exit ; end if
        end do
        do j=1,nvmax
          if (cvei(kk,j)>ncels+nnous) then ; kkk=kkk+1 ; exit ; end if
        end do
        if (kkk==2) then
          do j=1,nvmax
            if (ccvei(i,j)==ii) then ; ij=j ; exit ; end if
          end do
          do j=1,nvmax
            if (ccvei(i,j)==kk) then ; jjj=j ; exit ; end if
          end do 
          if (ij>jjj) then 
            ji=ij ; ij=jjj
          else
            ji=jjj
          end if
          !connect between ij and ji
          if (ji-ij==1) then
            do kk=nvmax,ji+1,-1 ; ccvei(i,kk)=ccvei(i,kk-1) ; cmarge(i,kk,4:8)=cmarge(i,kk-1,4:8) ; end do 
            ccvei(i,ji)=nncels;
          else
            ccvei(i,ji+1)=nncels 
          end if
        end if 
        cycle 
      end if
      goto 88
    end do	!End of the big loop

    !Update the matrices
    cvei(ncels+1:ncels+nnous,:)=ccvei(1:nnous,:)

    deallocate(malla)
    allocate(malla(nncels,3))
    malla=cmalla

    !calculate the distances between new cells -> marge(i, j, 4:8)
    do i=ncels+1,ncels+nnous
      ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
      do j=1,nvmax
        ii=cvei(i,j)
        if (ii>0.and.ii<ncels+nnous+1) then
          ux=malla(ii,1) ; uy=malla(ii,2) ; uz=malla(ii,3) 
          ux=ux-ua ; uy=uy-ub ; uz=uz-uc
          if (abs(ux)<10D-14) ux=0.
          if (abs(uy)<10D-14) uy=0.
          if (abs(uz)<10D-14) uz=0.
          d=sqrt(ux**2+uy**2)
          cmarge(i,j,5)=d
          d=sqrt(ux**2+uy**2+uz**2)
          cmarge(i,j,4)=d
          cmarge(i,j,6)=ux ;  cmarge(i,j,7)=uy ;  cmarge(i,j,8)=uz 
        end if
      end do
    end do

    !calculate the distances between new and old cells -> marge(i, j, 4:8)
    do i=1,ncels
      ua=malla(i,1) ; ub=malla(i,2) ; uc=malla(i,3)
      do j=1,nvmax
        ii=cvei(i,j)
        if (ii>ncels.and.ii<ncels+nnous+1) then
          ux=malla(ii,1) ; uy=malla(ii,2) ; uz=malla(ii,3) 
          ux=ux-ua ; uy=uy-ub ; uz=uz-uc
          if (abs(ux)<10D-14) ux=0.
          if (abs(uy)<10D-14) uy=0.
          if (abs(uz)<10D-14) uz=0.
          d=sqrt(ux**2+uy**2)
          cmarge(i,j,5)=d
          d=sqrt(ux**2+uy**2+uz**2)
          cmarge(i,j,4)=d
          cmarge(i,j,6)=ux ;  cmarge(i,j,7)=uy ;  cmarge(i,j,8)=uz 
        end if
      end do  
    end do  

    ncals=nncels
    ancels=ncels
    ncels=ncels+nnous

		!Reallocate all the matrices (use the copies made)
    deallocate(vei)    ; deallocate(hmalla)  ; deallocate(hvmalla)
    deallocate(marge) ;  deallocate(nveins) ; deallocate(q2d)     ; deallocate(q3d)
    deallocate(px)    ;  deallocate(py)     ; deallocate(pz)      ; deallocate(knots)

    allocate(vei(nncels,nvmax)) ; allocate(hmalla(nncels,3))  
    allocate(hvmalla(nncels,3))
    allocate(marge(nncels,nvmax,8)) ;  allocate(nveins(nncels))    ; allocate(q2d(nncels,ngg))   
    allocate(q3d(nncels,ncz,ng))
    allocate(px(nncels)) ; allocate(py(nncels)) ; allocate(pz(nncels)) ; allocate (knots(nncels))

    vei=cvei ; nveins=cnveins ; q2d=cq2d ; q3d=cq3d ; knots=cknots ; marge=cmarge ; hmalla=0.
    deallocate(cmalla)  ; deallocate(cvei)    
    deallocate(cnveins) ; deallocate(cq2d)  ; deallocate(cq3d)
    deallocate(ccvei)   ; deallocate(cknots); deallocate(cmarge) 

    !delete the zeros
    do i=1,ncals
      do j=1,nvmax
        if (vei(i,j)==0) then
          do jj=j,nvmax-1
            vei(i,jj)=vei(i,jj+1)
          end do
        end if
      end do
    end do
	
		!update nveins
    do i=1,ncals
      ii=0
      do j=1,nvmax
        if (vei(i,j)>0) ii=ii+1
      end do
      nveins(i)=ii
    end do
  end if   !end of big if loop (if nnous>0)

  do iii=1,nnous
    if (externsa(iii)==1) then  !tenim una nova cel externa
      ii=ancels+iii
      if (ncils==centre) centre=ii
      allocate(scvei(ncels,nvmax))
      allocate(scnveins(ncels))
      allocate(scmalla(ncels,3))
      allocate(scq3d(ncels,ncz,ng))
      allocate(scq2d(ncels,ngg))
      allocate(scmarge(ncels,nvmax,8))
      allocate(scknots(ncels))
      scmalla=malla
      scvei=vei
      scq2d=q2d
      scq3d=q3d
      scnveins=nveins
      scmarge=marge
      scknots=knots
      vei(ii,:)=scvei(ncils,:)
      vei(ncils,:)=scvei(ii,:)
      nveins(ii)=scnveins(ncils)
      nveins(ncils)=scnveins(ii)
      malla(ii,:)=scmalla(ncils,:)
      malla(ncils,:)=scmalla(ii,:)
      q2d(ii,:)=scq2d(ncils,:)
      q2d(ncils,:)=scq2d(ii,:)
      q3d(ii,:,:)=scq3d(ncils,:,:)
      q3d(ncils,:,:)=scq3d(ii,:,:)
      marge(ii,:,:)=scmarge(ncils,:,:)
      marge(ncils,:,:)=scmarge(ii,:,:)
      knots(ii)=scknots(ncils)
      knots(ncils)=scknots(ii)
      scvei=vei

      do i=1,ncels
        do j=1,nvmax
          if (scvei(i,j)==ii) vei(i,j)=ncils
        end do
      end do
      do i=1,ncels
        do j=1,nvmax
         if (scvei(i,j)==ncils) vei(i,j)=ii
        end do
      end do
      deallocate(scvei)
      deallocate(scnveins)
      deallocate(scmalla)
      deallocate(scq3d)
      deallocate(scq2d)
      deallocate(scmarge)
      deallocate(scknots)
      ncils=ncils+1
    end if
  end do

  if (nnous>0) then
    call perextrems
  end if

end subroutine afegircel

! identifies which of the cells added in afegircel are in the border of the tooth
subroutine perextrems
 integer,allocatable :: mau(:)
 integer nousn(ncils)
 integer nunous,nmaaa

     nousn=0
     nunous=0
er:  do i=1,ncils-1
       if (i==3.or.i==6) cycle 
       kk=0
       do ii=1,nmaa
         iii=mmaa(ii)
         if (iii==i) cycle er
       end do
err:   do j=1,nvmax

         k=vei(i,j)
         if (k<ncils) then
           do ii=1,nmaa
             iii=mmaa(ii)
             if (k==iii) then
               if (kk==1) then
                 nunous=nunous+1
                 nousn(nunous)=i    
                 cycle er
               else
                 kk=1
                 cycle err
               end if
             end if
           end do
         end if
       end do err
     end do er

     if (nunous>0) then
       nmaaa=nmaa
       allocate(mau(nmaa))
       mau=mmaa
       deallocate(mmaa)
       nmaa=nmaa+nunous
       allocate(mmaa(nmaa))
       mmaa(1:nmaa-nunous)=mau
       do i=1,nunous
         mmaa(nmaaa+i)=nousn(i)
       end do 
       deallocate(mau)
     end if

     nousn=0
     nunous=0
era:  do i=1,ncils-1
       if (i==3.or.i==6) cycle 
       kk=0
       do ii=1,nmap
         iii=mmap(ii)
         if (iii==i) cycle era
       end do
erra:   do j=1,nvmax
         k=vei(i,j)
         if (k<ncils) then
           do ii=1,nmap
             iii=mmap(ii)
             if (k==iii) then
               if (kk==1) then
                 nunous=nunous+1
                 nousn(nunous)=i    
                 cycle era
               else
                 kk=1
                 cycle erra
               end if
             end if
           end do
         end if
       end do erra
     end do era

     if (nunous>0) then
       nmaaa=nmap
       allocate(mau(nmap))
       mau=mmap
       deallocate(mmap)
       nmap=nmap+nunous
       allocate(mmap(nmap))
       mmap(1:nmap-nunous)=mau
       do i=1,nunous
         mmap(nmaaa+i)=nousn(i)
       end do 
       deallocate(mau)
     end if
end subroutine

subroutine iteracio
    panic=0
    hmalla=0.

    call reaccio_difusio
    if (panic==1) return
    call biaixbl
    call diferenciacio
    call empu
    call stelate
    call pushingnovei
    call pushing
    call biaixbld
    call promig
    call actualitza
    call afegircel
    call calculmarges

    temps=temps+1
   
end subroutine iteracio

end module coreop2d

!***************************************************************************
!***************  MOQUL ***************************************************
!***************************************************************************
module esclec

	use coreop2d
	public:: guardaforma,guardaveins,guardapara,llegirforma,llegirveins,llegirpara,llegir
	character*50, public :: fifr,fivr,fipr,fifw,fivw,fipw 
	integer, public :: nom,map,fora,pass,passs,maptotal,is,maxll
	integer, parameter :: mamax=5000
	real*8,  public :: mallap(1000,3,mamax)  !atencio si ncels>1000 el sistema peta al llegir
	real*8,  public :: parap(30,mamax)
	integer, public :: knotsp(1000,mamax)
	integer, public, allocatable :: veip(:,:,:)
	real*8, public,allocatable :: ma(:)
	character*30, public :: cac
	real*8, public :: vamax,vamin

	contains


subroutine guardaveinsoff(cvei)
  integer cvei(ncals,nvmax),ccvei(nvmax)
  real*8 c(4),mic(4)
    !em de passar de veinatge a faces

  integer ja(ncels*6,4)

  integer face(ncels*20,5)
  integer nfa(ncels*20)
  integer nfaces
  integer pasos(10)
  integer npasos,bi,nop
  real*8 mamx

  integer nfacre

  allocate(ma(ncels))
    call mat

  nfaces=0
  nfacre=0

    do i=1,ncels
ale: do j=1,nvmax
       bi=0
       ii=vei(i,j) ; if (ii==0.or.ii>ncels) cycle
ele:   do k=1,nvmax
         iii=vei(ii,k) ; if (iii==0.or.iii>ncels.or.iii==i) cycle
         do kk=1,nvmax
           iiii=vei(iii,kk) ; if (iiii==0.or.iiii>ncels) cycle
           if (iiii==i) then !triangle trobat
             nfaces=nfaces+1
7891         bi=bi+1
             nop=iii     
             if (bi==1) cycle ele
             cycle ale
           end if
         end do
       end do ele
      
       do k=1,nvmax
         iii=vei(ii,k) ; if (iii==0.or.iii>ncels.or.iii==i.or.iii==nop) cycle
         if (bi==0) cycle
         ! a per els quadrats
         do kk=1,nvmax
           iiii=vei(iii,kk) ; 
           if (iiii==0.or.iiii>ncels.or.iiii==ii.or.iiii==nop) cycle
           do kkk=1,nvmax
             jj=vei(iiii,kkk)
             if (jj==i) then !triangle trobat
               nfaces=nfaces+1
               !write (2,*) nfaces,i,ii,iii,iiii
               cycle ale
             end if
           end do
         end do
       end do

    end do ale
  end do
  
  write (2,*) "COFF"  
  write (2,*) ncels,nfaces,ncels
  write (2,*) " "
  do i=1,ncels ; call get_rainbow(ma(i),vamin,vamax,c) ; write (2,*) malla(i,:),c ; end do

nfaces=0
ja=0
nfacre=0

  write (2,*) " "    
    do i=1,ncels
aale: do j=1,nvmax
       bi=0
       ii=vei(i,j) ; if (ii==0.or.ii>ncels) cycle
aele:   do k=1,nvmax
         iii=vei(ii,k) ; if (iii==0.or.iii>ncels.or.iii==i) cycle
         do kk=1,nvmax
           iiii=vei(iii,kk) ; if (iiii==0.or.iiii>ncels) cycle
           if (iiii==i) then !triangle trobat
               nfaces=nfaces+1
               call get_rainbow(ma(i),vamin,vamax,c)
               mic=c ; mamx=ma(i)
               call get_rainbow(ma(ii),vamin,vamax,c)
               mic=mic+c
               !if (ma(ii)>mamx) then ; mic=c ; mamx=ma(ii) ; end if
               call get_rainbow(ma(iii),vamin,vamax,c)
               mic=mic+c
               !if (ma(iii)>mamx) then ; mic=c; end if 
               mic=mic/3.
             write (2,67) 3,i-1,ii-1,iii-1 !,mic
67 format (4I4,4F10.6)
789          bi=bi+1
             nop=iii          
             if (bi==1) cycle aele
             cycle aale
           end if
         end do
       end do aele
      
       do k=1,nvmax
         iii=vei(ii,k) ; if (iii==0.or.iii>ncels.or.iii==i.or.iii==nop) cycle
         if (bi==0) cycle
         ! a per els quadrats
         do kk=1,nvmax
           iiii=vei(iii,kk) ; 
           if (iiii==0.or.iiii>ncels.or.iiii==ii.or.iiii==nop) cycle
           do kkk=1,nvmax
             jj=vei(iiii,kkk)
             if (jj==i) then !triangle trobat
               nfaces=nfaces+1
               call get_rainbow(ma(i),vamin,vamax,c)
               mic=c ; mamx=ma(i)
               call get_rainbow(ma(ii),vamin,vamax,c)
               mic=mic+c 
               !if (ma(ii)>mamx) then ; mic=c ; mamx=ma(ii) ; end if 
               call get_rainbow(ma(iii),vamin,vamax,c)
               mic=mic+c ; 
               !if (ma(iii)>mamx) then ; mic=c ; mamx=ma(iii) ; end if 
               call get_rainbow(ma(iiii),vamin,vamax,c)
               mic=mic+c ; 
               !if (ma(iiii)>mamx) then ; mic=c ; mamx=ma(iiii) ; end if 
               mic=c/4.
               write (2,68) 4,i-1,ii-1,iii-1,iiii-1 !,mic
68 format (5I4,4F10.6)
               cycle aale
             end if
           end do
         end do
       end do

    end do aale
  end do
deallocate(ma)
end subroutine guardaveinsoff

subroutine mat
  ma=0
      do i=1,ncels
        if (knots(i)==1) then ; ma(i)=1.0
        else
          if (q2d(i,1)>us) ma(i)=0.1
          if (q2d(i,1)>ud) ma(i)=1.0
        end if
      end do
  vamax=maxval(ma)
  vamin=minval(ma)
end subroutine

subroutine get_rainbow(val,minval,maxval,c)
real*8, intent(in) :: val,maxval,minval
real*8, intent(out) :: c(4)

real*8 :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5
endif

if (f < .07) then
   c(1) = 0.6
   c(2) = 0.6
   c(3) = 0.6
   c(4) = 0.8
elseif (f < .2) then
   c(1) = 1.0
   c(2) = f
   c(3) = 0.0
   c(4) = 0.5
elseif (f < 1.0) then
   c(1) = 1.0
   c(2) = f*3
   c(3) = 0.0
   c(4) = 1.0
else
   c(1) = 1.0
   c(2) = 1.0
   c(3) = 0.0
   c(4) = 1.0
endif

end subroutine get_rainbow

subroutine llegirparatxt
  character*20 cf

  do i=3,29
    read (2,*,END=666,ERR=777)  a ; parap(i,map)=a ; print *,i,a 
  end do
	parap(30,map)=0.8

5 format(5F13.6)
  Return
889 print *,"fi";return
999 print *,"fi";return
888 print *,"error";return

777 print *,"error de lectura para" ; fora=1 ; close(2) ; return
666 print *,"fi de fitxer para"     ; print *,parap(1:5,map); fora=1 ; close(2) ; return
end subroutine llegirparatxt

subroutine posarparap(imap)

integer imap

  temps=parap(1,imap) ; ncels=parap(2,imap)
  tacre=parap(3,imap) ; tahor=parap(4,imap) ; elas=parap(5,imap) ; tadi=parap(6,imap)  ; crema=parap(7,imap)
  acac=parap(8,imap)  ; ihac=parap(9,imap) ; acaca=parap(10,imap); ih=parap(11,imap)   ; acec=parap(12,imap)
  do j=1,ng  ; difq3d(j)=parap(12+j,imap)    ; end do
  do j=1,ngg ; difq2d(j)=parap(12+ng+j,imap) ; end do
  us=parap(17,imap) ; ud=parap(18,imap);
  bip=parap(13+ng+ngg,imap) ; bia=parap(14+ng+ngg,imap) ; bib=parap(15+ng+ngg,imap) ; bil=parap(16+ng+ngg,imap)
  radi=parap(17+ng+ngg,imap); mu=parap(18+ng+ngg,imap)  ; tazmax=parap(19+ng+ngg,imap)
  radibi=parap(20+ng+ngg,imap) ; tadif=parap(14+ng+1,imap) 
  fac=parap(15+ng+1,imap) ; radibii=parap(21+ng+ngg,imap) 

end subroutine posarparap

subroutine llegirinicial

	open(2, file=cac, status='old', iostat=i)
	map=1
	call llegirparatxt
	call posarparap(map)
	close(2)

end subroutine llegirinicial

subroutine gnuoutputxyz
	integer row

	open(17, file = "GnuOutput.txt")

	do row = 1, ncels
  	write(17, *) malla(row,:)
	end do
	
	close(17)

end subroutine gnuoutputxyz

subroutine gnuoutputnet
	integer first
	integer second

	open(17, file = "GnuOutput.txt")

	do first = 1, ncels
  	do second = 1, nvmax
    	if (vei(first, second) > 0) then;
    	  write(17, *) malla(first,:)
    	  write(17, *) malla(vei(first, second),:)
    	  write(17, *)
    	end if
  	end do
	end do

	close(17)

end subroutine gnuoutputnet

!produces a screenshot of the tooth and makes output to plot with mesh in gnuplot
subroutine rolandGnuOutput(inputfile, progress)

	character*30 :: inputfile
	integer :: progress

  nfioff=""
  nfioff=trim(inputfile)//".out"//trim(progress)

	open(22,file=nfioff,iostat=i)

	call guardaveinsoff(vei)

  close(22)

  call execute_command_line("./MMFG "//trim(nfioff)//" & ",wait=.true.)

end subroutine rolandGnuOutput

end module esclec

!***************************************************************************
!***************  PROGRAMA           ****************************************
!***************************************************************************

program tresdac

	use coreop2d
	use esclec
	implicit none

	integer :: iteedone,itee,iteestart
	character*10 :: nfi,iterall
	character*60 nfioff!14
	character*4 iteestartc4
	character*5 iteestartc5
	character*6 iteestartc


	! Commandline inputs
	call getarg(1,cac)				!Inputfile
	call getarg(2,iterall)		!number of iterations

	if (cac.eq. "") then; 
  	print *,"you need to indicate an input file after the name of the command" ;        
  	print *,"in the form muscmd.e parameterfile.txt"
  	goto 666 ;
	end if

	!Program
	call ciinicial
	call llegirinicial
	call dime

	temps=0
	pass=0
	maptotal=0

10 continue

	if(len(trim(iterall))>0)then
  	read(iterall,*) itee
	else
		print *,"How many iterations to run; enter -1 to stop the program"
  	read(*,*) itee
  	if (itee==-1) goto 666
  	iterall=""
	endif

	do ite=1, itee
		call iteracio

		!every 1000 iterations, make an outputfile
		if (mod(1000, ite)==0) then
			call rolandGnuOutput(cac, ite)
		end if
	end do

	call gnuoutputnet

	666 print *,"out"
end program tresdac

!***************************************************************************
!***************  FI  PROGRAMA      ****************************************
!***************************************************************************


