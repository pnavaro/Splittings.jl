module splinepp_class
  use used_precision
  use geometry_module
  use clock
  implicit none
  private
  public :: new, interpole
  type, public :: splinepp
     type (geometry) :: geom
     real(wp) :: a1x, a2x, a3x, a4x ! coef de la matrice 2x2 per.
     real(wp) :: a1y, a2y, a3y, a4y ! coef de la matrice 2x2 nat.
     real(wp), dimension(:), pointer :: axd, ayd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: axod, ayod ! termes sous-diag. de ax
     real(wp), dimension(:,:), pointer :: axm1gamma, aym1gamma 
     real(wp), dimension(:,:), pointer :: coef, bcoef ! coefficients des splines
  end type splinepp
  interface new
     module procedure new_splinepp
  end interface
  interface interpole
     module procedure interpole_splinepp,interpole_splineppdep
  end interface
contains
  subroutine new_splinepp(this,geom,iflag)
    type(splinepp), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag
    ! local variables
    integer :: err ! error flag
    integer :: n1,n2 ! dimensions
    integer :: n1p1, n1p2, n1p3, n2p1, n2p2, n2p3
    real(wp) :: aa1x, aa1y, aa2x, aa2y
    integer i,j, info
    ! initialisation des variables locales
    iflag = 0
    n1 = geom%n1
    n2 = geom%n2
    n1p1=n1+1
    n1p2=n1+2
    n1p3=n1+3
    n2p1=n2+1
    n2p2=n2+2
    n2p3=n2+3

    ! initialisation de geom
    this%geom = geom

    ! memory allocation
    allocate(this%axm1gamma(n1p1,2), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%aym1gamma(n2p1,2), stat=err)
    if (err.ne.0) then
       iflag = 15
       return
    end if
    allocate(this%axd(n1p1), stat=err)
    if (err.ne.0) then
       iflag = 20
       return
    end if
    allocate(this%axod(n1), stat=err)
    if (err.ne.0) then
       iflag = 30
       return
    end if
    allocate(this%ayd(n2p1), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%ayod(n2), stat=err)
    if (err.ne.0) then
       iflag = 45
       return
    end if
    allocate(this%coef(n1p3,n2p3), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%bcoef(n1p3,n2), stat=err)
    if (err.ne.0) then
       iflag = 60
       return
    end if

    ! factorize matrices  Ax and Ay
    this%axd = 4_wp
    this%axod = 1_wp
#ifdef _CRAY
    call spttrf(n1p1,this%axd,this%axod,err)
#else
    call dpttrf(n1p1,this%axd,this%axod,err)
#endif
    if (err.ne.0) then
       iflag = 70
       return
    end if
    this%ayd = 4_wp
    this%ayod = 1_wp
#ifdef _CRAY
    call spttrf(n2p1,this%ayd,this%ayod,err)
#else
    call dpttrf(n2p1,this%ayd,this%ayod,err)
#endif
    if (err.ne.0) then
       iflag = 80
       return
    end if
    ! compute Ax-1.gamma
    this%axm1gamma = 0_wp
    this%axm1gamma(1,2) = 1_wp
    this%axm1gamma(n1p1,1) = 1_wp
#ifdef _CRAY
    call spttrs(n1p1,2,this%axd,this%axod,this%axm1gamma,n1p1,err)
#else
    call dpttrs(n1p1,2,this%axd,this%axod,this%axm1gamma,n1p1,err)
#endif
    if (err.ne.0) then
       iflag = 90
       return
    end if
    ! compute Ay-1.gamma
    this%aym1gamma = 0_wp
    this%aym1gamma(1,2) = 1_wp
    this%aym1gamma(n2p1,1) = 1_wp
#ifdef _CRAY
    call spttrs(n2p1,2,this%ayd,this%ayod,this%aym1gamma,n2p1,err)
#else
    call dpttrs(n2p1,2,this%ayd,this%ayod,this%aym1gamma,n2p1,err)
#endif
    if (err.ne.0) then
       iflag = 100
       return
    end if

    aa1x=3_wp/geom%delta1
    aa1y=3_wp/geom%delta2
    aa2x=6_wp/(geom%delta1*geom%delta1)
    aa2y=6_wp/(geom%delta2*geom%delta2)
    ! assemblage de la matrice 2x2 pour la spline dans la direction Ox
       this%a1x = -aa1x*(1. + this%axm1gamma(2,1)+this%axm1gamma(n1,1))
       this%a2x = -aa1x*(1. + this%axm1gamma(2,2)+this%axm1gamma(n1,2))
       this%a3x = aa2x*(-1. + 2*this%axm1gamma(1,1) - this%axm1gamma(2,1) &
            + this%axm1gamma(n1,1) - 2*this%axm1gamma(n1p1,1))
       this%a4x = aa2x*(1. + 2*this%axm1gamma(1,2) - this%axm1gamma(2,2) &
            + this%axm1gamma(n1,2) - 2*this%axm1gamma(n1p1,2))
    ! assemblage de la matrice 2x2 pour spline naturels (direction Oy)
       this%a1y = -aa1y*(1. + this%aym1gamma(2,1)+this%aym1gamma(n2,1))
       this%a2y = -aa1y*(1. + this%aym1gamma(2,2)+this%aym1gamma(n2,2))
       this%a3y = aa2y*(-1. + 2*this%aym1gamma(1,1) - this%aym1gamma(2,1) &
            + this%aym1gamma(n2,1) - 2*this%aym1gamma(n2p1,1))
       this%a4y = aa2y*(1. + 2*this%aym1gamma(1,2) - this%aym1gamma(2,2) &
            + this%aym1gamma(n2,2) - 2*this%aym1gamma(n2p1,2))
 end subroutine new_splinepp

  subroutine interpole_splinepp(this,fin,fout,x,y) 
    type(splinepp), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr

    call per_x(this,fin,ierr)

    if (ierr.ne.0) then
       iflag = 10
       return
    end if

    call per_y(this,ierr)

    if (ierr.ne.0) then
       iflag = 20
       return
    end if

    call evaltab(this,x,y,fout)     

  end subroutine interpole_splinepp

  subroutine interpole_splineppdep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinepp), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    logical  :: aff
    integer ierr, l_a, l_b
    real(wp) vtime(1:4)

    if (aff) then 
       call clck_temps(l_a)
    end if

    call per_x(this,f,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if
    if (aff) then 
       call clck_temps(l_b)
       call clck_diff(l_a,l_b,vtime(1))
       call clck_temps(l_a)
    end if

    call per_y(this,ierr)
    if (ierr.ne.0) then
       iflag = 20
       return
    end if
    if (aff) then 
       call clck_temps(l_b)
       call clck_diff(l_a,l_b,vtime(2))
       call clck_temps(l_a)
    end if

    call evaldep(this,depx,depy,f)     

    if (aff) then 
       call clck_temps(l_b)
       call clck_diff(l_a,l_b,vtime(3))
       write(*,'(A,3(1x,3E14.5))') "splinepp ",vtime(1:3)
    end if

  end subroutine interpole_splineppdep

  !
  ! calcul des splines periodiques
  !
  subroutine per_x(this,gtau,iflag)
    type(splinepp), intent(inout) :: this     ! objet de type splinepp
    real(wp), dimension(:,:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer n1, n2, n1p2, n1p3  !  n1+2, n1+3
    integer :: err ! error flag

    real(wp) :: axm1f(this%geom%n1+1,this%geom%n2)
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    n1=this%geom%n1
    n2=this%geom%n2
    n1p2=n1+2
    n1p3=n1+3
    det=this%a1x*this%a4x - this%a2x*this%a3x

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    do j=1,n2
       do i=1,n1
          axm1f(i,j) = 6*gtau(i,j)
       end do
       axm1f(n1+1,j) = 6*gtau(1,j)  ! calcul par periodicite
    end do

#ifdef _CRAY
    call spttrs(n1+1,n2,this%axd,this%axod,axm1f,n1+1,err)
#else 
    call dpttrs(n1+1,n2,this%axd,this%axod,axm1f,n1+1,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    do  j=1,n2
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  - (3.0/this%geom%delta1)*(axm1f(2,j) +axm1f(this%geom%n1,j))
       gamma2 =  (6.0/(this%geom%delta1)**2)*(2*axm1f(1,j) - axm1f(2,j) &
            + axm1f(this%geom%n1,j) - 2*axm1f(this%geom%n1+1,j))

       this%bcoef(n1p3,j)= (gamma1*this%a4x - gamma2*this%a2x)/det
       this%bcoef(1,j)= (gamma2*this%a1x - gamma1*this%a3x)/det
       do  i=2,n1p2
          this%bcoef(i,j)= axm1f(i-1,j) &
               - this%axm1gamma(i-1,1)*this%bcoef(n1p3,j) &
               - this%axm1gamma(i-1,2)*this%bcoef(1,j)
       end do
    end do
  end subroutine per_x
  !
  subroutine per_y(this,iflag)
    type(splinepp), intent(inout) :: this     ! objet de type spline
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer n2, n1p3, n2p3  !  n1+3, n2+3
    integer :: err ! error flag

    real(wp) :: aym1f(this%geom%n2+1,this%geom%n1+3)
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    n2 = this%geom%n2
    n1p3=this%geom%n1+3
    n2p3=this%geom%n2+3
    det=this%a1y*this%a4y - this%a2y*this%a3y

    ! calcul de coef par resolution de n1p2 systemes lineaires.

    ! Calcul de Ay^-1 f
    ! assemblage du membre de droite pour le calcul de Ay^-1 f
    do i=1,n1p3
       do j=1,n2
          aym1f(j,i) = 6.*this%bcoef(i,j)
       end do
       aym1f(n2+1,i) = 6.*this%bcoef(i,1)
    end do
#ifdef _CRAY
    call spttrs(n2+1,n1p3,this%ayd,this%ayod,aym1f,n2+1,err)
#else
    call dpttrs(n2+1,n1p3,this%ayd,this%ayod,aym1f,n2+1,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    do i=1,n1p3
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  - (3.0/this%geom%delta2)*(aym1f(2,i) +aym1f(this%geom%n2,i))
       gamma2 =  (6.0/(this%geom%delta2)**2)*(2*aym1f(1,i) - aym1f(2,i) &
            + aym1f(this%geom%n2,i) - 2*aym1f(this%geom%n2+1,i))
       ! resolution du syteme lineaire 2x2
       this%coef(i,n2p3) = (gamma1*this%a4y - gamma2*this%a2y)/det
       this%coef(i,1) = (gamma2*this%a1y - gamma1*this%a3y)/det
       do  j=2,this%geom%n2 + 2
          this%coef(i,j)= aym1f(j-1,i)                   &
               - this%aym1gamma(j-1,1)*this%coef(i,n2p3) &
               - this%aym1gamma(j-1,2)*this%coef(i,1)
       end do
    end do
  end subroutine per_y

  subroutine evaltab(this,xd,yd,fout)
    ! ------------------------------------------------------
    ! Evalue la spline en tous les points (xd(i,j), yd(i,j))
    !-------------------------------------------------------
    type(splinepp) :: this
    ! coordonnees du point ou les valeurs sont calculees
    real(wp), dimension(:,:) :: xd, yd 
    ! fout(i,j)contient la valeur de la spline au point (xd(i,j), yd(i,j))
    real(wp), dimension(:,:) :: fout

    ! variables locales
    real(wp) :: sval   ! valeur de la fonction au point (xd,yd)
    real(wp) :: bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3,bvaly4
    real(wp) :: a1,delta1x,delta1xx,delta1xx6,delta2y,delta2yy,delta2yy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4, idelta1, idelta2, lx, ly
    integer i1,j1,i,j
    !
    delta1x=this%geom%delta1*this%geom%delta1
    delta1xx=delta1x*this%geom%delta1
    delta1xx6=1./(6.*delta1xx)
    !
    !
    delta2y=this%geom%delta2*this%geom%delta2
    delta2yy=delta2y*this%geom%delta2
    delta2yy6=1./(6.*delta2yy)

    idelta1 = 1/this%geom%delta1
    idelta2 = 1/this%geom%delta2

    lx = (this%geom%n1-1)*this%geom%delta1
    ly = (this%geom%n2-1)*this%geom%delta2
    !print*,'l ',lx,ly

    ! evaluation des splines pour le calcul des n1-1 lignes 
    ! et n2-1 colonnes de fout
    do j=1,this%geom%n2
       do i=1,this%geom%n1
!!$          do while (xd(i,j)<this%geom%x0)
!!$             !print*,i,j,xd(i,j),this%geom%x0
!!$             xd(i,j)=xd(i,j)+lx
!!$          end do
!!$          do while (yd(i,j)<this%geom%y0)
!!$             !print*,i,j,yd(i,j),this%geom%y0
!!$             yd(i,j)=yd(i,j)+ly
!!$          end do
!!$          do while (xd(i,j)>this%geom%x0+lx)
!!$             !print*,i,j,xd(i,j),this%geom%x0
!!$             xd(i,j)=xd(i,j)-lx
!!$          end do
!!$          do while (yd(i,j)>this%geom%y0+ly)
!!$             !print*,i,j,yd(i,j),this%geom%y0
!!$             yd(i,j)=yd(i,j)-ly
!!$          end do

          i1=(xd(i,j)-this%geom%x0)*idelta1
          j1=(yd(i,j)-this%geom%y0)*idelta2

          if ((i1<0).or.(i1>this%geom%n1-1)) print*,'i1',i1
          if ((j1<0).or.(j1>this%geom%n2-1)) print*,'j1',j1
          
          !i1=mod(i1+this%geom%n1,this%geom%n1)
          !j1=mod(j1+this%geom%n2,this%geom%n2)
          !print*,'rrr ',i,j,xd(i,j),yd(i,j),i1,j1,this%geom%n1,this%geom%n2

          xdp1=this%geom%xgrid(i1+2)-xd(i,j)
          bvalx1=xdp1*xdp1*xdp1
          bvalx2=delta1xx+3.*delta1x*xdp1+3.*this%geom%delta1* &
               &xdp1*xdp1-3.*xdp1*xdp1*xdp1
          xd1=xd(i,j)-this%geom%xgrid(i1+1)
          bvalx3=delta1xx+3.*delta1x*xd1+3.*this%geom%delta1* &
               &xd1*xd1-3.*xd1*xd1*xd1
          bvalx4=xd1*xd1*xd1
          ydp1=this%geom%ygrid(j1+2)-yd(i,j)
          bvaly1=ydp1*ydp1*ydp1
          bvaly2=delta2yy+3.*ydp1*(delta2y+ydp1*(this%geom%delta2-ydp1))          
          yd1=yd(i,j)-this%geom%ygrid(j1+1)
          bvaly3=delta2yy+3.*yd1*(delta2y+yd1*(this%geom%delta2-yd1))
          bvaly4=yd1*yd1*yd1

          sval=0.
          sval1=this%coef(i1+1,j1+1)*bvaly1
          sval1=sval1+this%coef(i1+1,j1+2)*bvaly2
          sval1=sval1+this%coef(i1+1,j1+3)*bvaly3
          sval1=sval1+this%coef(i1+1,j1+4)*bvaly4        
          sval= sval+sval1*bvalx1
          sval2=this%coef(i1+2,j1+1)*bvaly1
          sval2=sval2+this%coef(i1+2,j1+2)*bvaly2
          sval2=sval2+this%coef(i1+2,j1+3)*bvaly3
          sval2=sval2+this%coef(i1+2,j1+4)*bvaly4
          sval=sval+sval2*bvalx2
          sval3=this%coef(i1+3,j1+1)*bvaly1
          sval3=sval3+this%coef(i1+3,j1+2)*bvaly2
          sval3=sval3+this%coef(i1+3,j1+3)*bvaly3
          sval3=sval3+this%coef(i1+3,j1+4)*bvaly4
          sval=sval+sval3*bvalx3
          sval4=this%coef(i1+4,j1+1)*bvaly1
          sval4=sval4+this%coef(i1+4,j1+2)*bvaly2
          sval4=sval4+this%coef(i1+4,j1+3)*bvaly3
          sval4=sval4+this%coef(i1+4,j1+4)*bvaly4
          sval=sval+sval4*bvalx4
         
          fout(i,j) = delta1xx6*delta2yy6*sval
       end do
    end do
 
  end subroutine evaltab

  subroutine evaldep(this,alphax,alphay,fout)
    type(splinepp) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphax,alphay 
    ! fout(i,j) contient en sortie la valeur de la spline au 
    ! point (xi-alphax,yj-alphay), les (xi,yj) etant les points du maillage
    real(wp), dimension(:,:), intent(out) :: fout

    ! variables locales
    real(wp) :: sval   ! valeur de la fonction au point d'evaluation
    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,delta1x,delta1xx,delta1xx6,delta2y,delta2yy,delta2yy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
    integer :: intaxsdelta1, intaysdelta2 
    integer i1,j1,i,j

    ! debut du code
    delta1x=this%geom%delta1*this%geom%delta1
    delta1xx=delta1x*this%geom%delta1
    delta1xx6=1./(6.*delta1xx)

    delta2y=this%geom%delta2*this%geom%delta2
    delta2yy=delta2y*this%geom%delta2
    delta2yy6=1./(6.*delta2yy)

    if (alphax.gt.0) then
       intaxsdelta1=int(-alphax/this%geom%delta1+epsilon)-1  
! intaxsdelta1=int(-alphax/this%geom%delta1)-1 
!       intaxsdelta1 = -1
    else
       intaxsdelta1=int(-alphax/this%geom%delta1)
    end if
    
    if (alphay.gt.0) then
       intaysdelta2=int(-alphay/this%geom%delta2+epsilon)-1
!intaysdelta2=int(-alphax/this%geom%delta1)-1    
!intaysdelta2=-1
    else
       intaysdelta2=int(-alphay/this%geom%delta2)
    end if
    xd1=-alphax-intaxsdelta1*this%geom%delta1
    xdp1=this%geom%delta1-xd1
    yd1=-alphay-intaysdelta2*this%geom%delta2
    ydp1=this%geom%delta2-yd1

    bvalx1=xdp1*xdp1*xdp1
    bvalx2=delta1xx+3.*delta1x*xdp1+3.*this%geom%delta1* &
         &xdp1*xdp1-3.*xdp1*xdp1*xdp1
    bvalx3=delta1xx+3.*delta1x*xd1+3.*this%geom%delta1* &
         &xd1*xd1-3.*xd1*xd1*xd1
    bvalx4=xd1*xd1*xd1
    bvaly1=ydp1*ydp1*ydp1
    bvaly2=delta2yy+3.*ydp1*(delta2y+ydp1*(this%geom%delta2-ydp1))          
    bvaly3=delta2yy+3.*yd1*(delta2y+yd1*(this%geom%delta2-yd1))
    bvaly4=yd1*yd1*yd1

    do j=1,this%geom%n2
       j1=mod(this%geom%n2+j-1+intaysdelta2,this%geom%n2)

       do i=1,this%geom%n1
          i1=mod(this%geom%n1+i-1+intaxsdelta1,this%geom%n1)

          fout(i,j) = delta1xx6*delta2yy6* ( &
               bvalx1*( this%coef(i1+1,j1+1)*bvaly1 &
               +this%coef(i1+1,j1+2)*bvaly2 &
               +this%coef(i1+1,j1+3)*bvaly3 &
               +this%coef(i1+1,j1+4)*bvaly4) &
               + bvalx2* (this%coef(i1+2,j1+1)*bvaly1 &
               +this%coef(i1+2,j1+2)*bvaly2 &
               +this%coef(i1+2,j1+3)*bvaly3 &
               +this%coef(i1+2,j1+4)*bvaly4) &
               + bvalx3* (this%coef(i1+3,j1+1)*bvaly1 &
               +this%coef(i1+3,j1+2)*bvaly2 &
               +this%coef(i1+3,j1+3)*bvaly3 &
               +this%coef(i1+3,j1+4)*bvaly4) &
               + bvalx4* (this%coef(i1+4,j1+1)*bvaly1 &
               +this%coef(i1+4,j1+2)*bvaly2 &
               +this%coef(i1+4,j1+3)*bvaly3 &
               +this%coef(i1+4,j1+4)*bvaly4))
       end do
    end do

  end subroutine evaldep
end module splinepp_class
