module splinenn_class
  use used_precision
  use geometry_module
  use clock
  implicit none
  private
  public :: new, interpole
  type, public :: splinenn
     type (geometry) :: geom
     real(wp) :: a1x, a2x, a3x, a4x ! coef de la matrice 2x2 per.
     real(wp) :: a1y, a2y, a3y, a4y ! coef de la matrice 2x2 nat.
     real(wp), dimension(:), pointer :: axd, ayd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: axod, ayod ! termes sous-diag. de ax
     real(wp), dimension(:), pointer :: aym1gamma1, aym1gamma2
     real(wp), dimension(:), pointer :: axm1gamma1, axm1gamma2
     real(wp), dimension(:,:), pointer :: coef, bcoef ! coefficients des splines
  end type splinenn
  interface new
     module procedure new_splinenn
  end interface
  interface interpole
     module procedure interpole_splinenn, interpole_splinenndep
  end interface
contains
  subroutine new_splinenn(this,geom,iflag)
    type(splinenn), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag
    ! local variables
    integer :: err ! error flag
    integer :: n1,n2 ! dimensions
    integer :: n1p1, n1p2, n2p1, n2p2
    real(wp) :: aa1x, aa1y, aa2x, aa2y
    integer i,j, info
    ! initialisation des variables locales
    iflag = 0
    n1 = geom%n1
    n2 = geom%n2
    n1p1=n1+1
    n1p2=n1+2
    n2p1=n2+1
    n2p2=n2+2

    ! initialisation de geom
    this%geom = geom

    ! memory allocation
    allocate(this%axm1gamma1(n1), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%axm1gamma2(n1), stat=err)
    if (err.ne.0) then
       iflag = 11
       return
    end if
    allocate(this%aym1gamma1(n2), stat=err)
    if (err.ne.0) then
       iflag = 15
       return
    end if
    allocate(this%aym1gamma2(n2), stat=err)
    if (err.ne.0) then
       iflag = 16
       return
    end if
    allocate(this%axd(n1), stat=err)
    if (err.ne.0) then
       iflag = 20
       return
    end if
    allocate(this%axod(n1-1), stat=err)
    if (err.ne.0) then
       iflag = 30
       return
    end if
    allocate(this%ayd(n2), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%ayod(n2-1), stat=err)
    if (err.ne.0) then
       iflag = 45
       return
    end if
    allocate(this%coef(n1p2,n2p2), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%bcoef(n2p2,n1p2), stat=err)
    if (err.ne.0) then
       iflag = 60
       return
    end if

    ! factorize matrices  Ax and Ay
    this%axd = 4_wp
    this%axod = 1_wp
#ifdef _CRAY
    call spttrf(n1,this%axd,this%axod,err)
#else
    call dpttrf(n1,this%axd,this%axod,err)
#endif
    if (err.ne.0) then
       iflag = 70
       return
    end if
    this%ayd = 4_wp
    this%ayod = 1_wp
#ifdef _CRAY
    call spttrf(n2,this%ayd,this%ayod,err)
#else
    call dpttrf(n2,this%ayd,this%ayod,err)
#endif
    if (err.ne.0) then
       iflag = 80
       return
    end if
    ! compute Ax-1.gamma
    this%axm1gamma1 = 0_wp
    this%axm1gamma2 = 0_wp
    this%axm1gamma2(1) = 1_wp
    this%axm1gamma1(n1) = 1_wp
#ifdef _CRAY
    call spttrs(n1,1,this%axd,this%axod,this%axm1gamma1,n1,err)
    call spttrs(n1,1,this%axd,this%axod,this%axm1gamma2,n1,err)
#else
    call dpttrs(n1,1,this%axd,this%axod,this%axm1gamma1,n1,err)
    call dpttrs(n1,1,this%axd,this%axod,this%axm1gamma2,n1,err)
#endif

    if (err.ne.0) then
       iflag = 90
       return
    end if
    ! compute Ay-1.gamma
    this%aym1gamma1 = 0_wp
    this%aym1gamma2 = 0_wp
    this%aym1gamma2(1) = 1_wp
    this%aym1gamma1(n2) = 1_wp
#ifdef _CRAY
    call spttrs(n2,1,this%ayd,this%ayod,this%aym1gamma1,n2,err)
    call spttrs(n2,1,this%ayd,this%ayod,this%aym1gamma2,n2,err)
#else
    call dpttrs(n2,1,this%ayd,this%ayod,this%aym1gamma1,n2,err)
    call dpttrs(n2,1,this%ayd,this%ayod,this%aym1gamma2,n2,err)
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
    this%a1x = aa2x*(1. - this%axm1gamma1(n1-1) + 2*this%axm1gamma1(n1))
    this%a2x = aa2x*(-this%axm1gamma2(n1-1) + 2*this%axm1gamma2(n1))
    this%a3x = aa2x*(2*this%axm1gamma1(1)-this%axm1gamma1(2))
    this%a4x = aa2x*(1. + 2*this%axm1gamma2(1) - this%axm1gamma2(2))
    ! assemblage de la matrice 2x2 pour spline naturels (direction Oy)
    this%a1y = aa2y*(1. - this%aym1gamma1(n2-1) + 2*this%aym1gamma1(n2))
    this%a2y = aa2y*(-this%aym1gamma2(n2-1) + 2*this%aym1gamma2(n2))
    this%a3y = aa2y*(2*this%aym1gamma1(1)-this%aym1gamma1(2))
    this%a4y = aa2y*(1. + 2*this%aym1gamma2(1) - this%aym1gamma2(2))
  end subroutine new_splinenn

  subroutine interpole_splinenn(this,fin,fout,x,y) 
    type(splinenn), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    integer :: iflag ! error flag
    ! variables locales
    integer i,j, ierr

    ! initialisation des variables locales
    iflag = 0


    call nat_x(this,fin,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if

    call nat_y(this,ierr)
    if (ierr.ne.0) then
       iflag = 20
       return
    end if

    call evaltab(this,x,y,fout)     

    return

  end subroutine interpole_splinenn

  subroutine interpole_splinenndep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinenn), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    logical :: aff
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr, l_a, l_b
    real(wp) durat, vtime(1:4)

    if (aff) then 
       call clck_temps(l_a)
    end if
    call nat_x(this,f,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if
    if (aff) then 
       call clck_temps(l_b)
       call clck_diff(l_a,l_b,vtime(1))
       call clck_temps(l_a)
    end if

    call nat_y(this,ierr)
    
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
       write(*,'(A,3(1x,3E14.5))') "splinenn ",vtime(1:3)
    end if

  end subroutine interpole_splinenndep

  !
  ! calcul des "natural splines"
  ! 
  subroutine nat_x(this,gtau,iflag)
    type(splinenn), intent(inout) :: this     ! objet de type spline
    real(wp), dimension(:,:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer n1, n2, n1p1, n1p2  !  n1+1, n1+2
    integer :: err ! error flag

    real(wp) :: axm1f(this%geom%n1,this%geom%n2)
    real(wp) :: det, gamma1, gamma2, coef1, coefnp2

    ! initialisations
    iflag =0
    n1=this%geom%n1
    n2=this%geom%n2
    n1p2=n1+2
    n1p1=n1+1
    det=this%a1x*this%a4x - this%a2x*this%a3x

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    do j=1,n2
       do i=1,n1
          axm1f(i,j) = 6*gtau(i,j)
       end do
    end do
#ifdef _CRAY
    call spttrs(n1,n2,this%axd,this%axod,axm1f,n1,err)
#else
    call dpttrs(n1,n2,this%axd,this%axod,axm1f,n1,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    !print*,'axmf1',axm1f
    do  j=1,n2
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  (6.0/(this%geom%delta1)**2)*(-axm1f(this%geom%n1-1,j) &
            + 2*axm1f(this%geom%n1,j))
       gamma2 = (6.0/(this%geom%delta1)**2)*(2*axm1f(1,j) - axm1f(2,j))

       coefnp2 = (gamma1*this%a4x - gamma2*this%a2x)/det
       coef1 = (gamma2*this%a1x - gamma1*this%a3x)/det
       this%bcoef(j,n1p2)=coefnp2
       this%bcoef(j,1)=coef1

       do  i=2,n1p1
          this%bcoef(j,i)= axm1f(i-1,j) &
               - this%axm1gamma1(i-1)*coefnp2 &
               - this%axm1gamma2(i-1)*coef1
       end do
    end do
  end subroutine nat_x

  subroutine nat_y(this,iflag)
    type(splinenn), intent(inout) :: this     ! objet de type spline
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j          ! indices de boucle
    integer :: n2, n1p2,n2p2 ! n1+2, n2+2
    integer :: err       ! indicateur d erreur

    real(wp) :: aym1f(this%geom%n2,this%geom%n1+2)
    real(wp) :: det, gamma1, gamma2, coef1, coefnp2

    ! initialisations
    iflag =0
    n2 = this%geom%n2
    n1p2=this%geom%n1+2
    n2p2=this%geom%n2+2
    det=this%a1y*this%a4y - this%a2y*this%a3y

    ! calcul de coef par resolution de n1p2 systemes lineaires.

    ! Calcul de Ay^-1 f
    ! assemblage du membre de droite pour le calcul de Ay^-1 f
    do  i=1,n1p2
       do j=1,this%geom%n2
          aym1f(j,i) = 6.*this%bcoef(j,i)
       end do
    end do
#ifdef _CRAY
    call spttrs(n2,n1p2,this%ayd,this%ayod,aym1f,n2,err)
#else
    call dpttrs(n2,n1p2,this%ayd,this%ayod,aym1f,n2,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    ! resolution du syteme lineaire 2x2
    do i=1,n1p2
       gamma1 =  (6.0/(this%geom%delta2)**2)*(-aym1f(this%geom%n2-1,i) &
            + 2*aym1f(this%geom%n2,i))
       gamma2 = (6.0/(this%geom%delta2)**2)*(2*aym1f(1,i) - aym1f(2,i))

       coefnp2 = (gamma1*this%a4y - gamma2*this%a2y)/det
       coef1 = (gamma2*this%a1y - gamma1*this%a3y)/det
       this%coef(i,n2p2) = coefnp2
       this%coef(i,1) = coef1
       do  j=2,this%geom%n2 + 1
          this%coef(i,j)= aym1f(j-1,i)                   &
               - this%aym1gamma1(j-1)*coefnp2 &
               - this%aym1gamma2(j-1)*coef1
       end do
    end do
  end subroutine nat_y
  !

  !
  ! tsvaleur
  !
  subroutine evaltab(this,xd,yd,fout)
    type(splinenn) :: this
    real(wp), dimension(:,:) :: xd, yd ! coordonnees du point ou les valeurs sont calculees

    real(wp) :: sval,idelta1,idelta2   ! valeur de la fonction au point (xd,yd)
    real(wp), dimension(:,:) :: fout

    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,a1,delta1x,delta1xx,delta1xx6,delta2y,delta2yy,delta2yy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
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
    do j=2,this%geom%n2-1
       do i=2,this%geom%n1-1

          i1=(xd(i,j)-this%geom%x0)*idelta1
          j1=(yd(i,j)-this%geom%y0)*idelta2

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

    fout(1,:)=0
    fout(this%geom%n1,:)=0 

    fout(:,1)=0
    fout(:,this%geom%n2)=0 
  end subroutine evaltab

  subroutine evaldep(this,alphax,alphay,fout)
    type(splinenn) :: this
    real(wp) :: alphax,alphay ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp) :: sval   ! valeur de la fonction au point d'evaluation
    real(wp), dimension(:,:) :: fout

    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,delta1x,delta1xx,delta1xx6,delta2y,delta2yy,delta2yy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
    integer :: intaxsdelta1, intaysdelta2
    integer i1,j1,i,j,ideb,ifin,jdeb,jfin
    !
    delta1x=this%geom%delta1*this%geom%delta1
    delta1xx=delta1x*this%geom%delta1
    delta1xx6=1./(6.*delta1xx)
    !
    !
    delta2y=this%geom%delta2*this%geom%delta2
    delta2yy=delta2y*this%geom%delta2
    delta2yy6=1./(6.*delta2yy)

!!$    if (abs(alphax).gt.this%geom%delta1) then
!!$       print*, 'deplacement en x trop grand, alphax=',alphax
!!$       print*, 'delta1=',this%geom%delta1
!!$       stop
!!$    end if
!!$    if (abs(alphay).gt.this%geom%delta2) then
!!$       print*,'deplacement en y trop grand, alphay=',alphay
!!$       print*, 'delta2=',this%geom%delta2
!!$       stop
!!$    end if

    if (alphax.gt.0) then
       intaxsdelta1=int(-alphax/this%geom%delta1+epsilon)-1 
!       i1=this%geom%n1-2
!       intaxsdelta1=-1
       ideb = int(alphax/this%geom%delta1)+2
       ifin = this%geom%n1 - 1
    else
!       intaxsdelta1=0
       intaxsdelta1=int(-alphax/this%geom%delta1)
!       i1=this%geom%n1-1
       ideb = 2
       ifin = -int(-alphax/this%geom%delta1)+this%geom%n1-1
    end if

    if (alphay.gt.0) then
!       intaysdelta2=-1
       intaysdelta2=int(-alphay/this%geom%delta2+epsilon)-1
!       j1=-1
       jdeb = int(alphay/this%geom%delta2)+2
       jfin = this%geom%n2 - 1
    else
!       intaysdelta2=0
       intaysdelta2=int(-alphay/this%geom%delta2)
!       j1=0
        jdeb = 2
        jfin = -int(-alphay/this%geom%delta2)+this%geom%n2-1
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
!print*,'eval ',xd1,yd1,xdp1,ydp1,intaxsdelta1,intaysdelta2
!print*,'eval ',this%geom%delta1,this%geom%delta2
!print*,'eval ',bvalx1,bvalx2,bvalx3,bvalx4
!print*,'eval ',bvaly1,bvaly2,bvaly3,bvaly4
!    do j=2,this%geom%n2-1
    do j=jdeb,jfin
!       j1=j1+1
!       i1=i1-(this%geom%n1-2) ! remise de i1 a sa valeur de depart
       j1=j-1+intaysdelta2
!j1=mod(this%geom%n2+j-2+intaysdelta2,this%geom%n2-1)
!       do i=2,this%geom%n1-1
       do i=ideb,ifin
!          i1=i1+1
          i1=i-1+intaxsdelta1
! i1=mod(this%geom%n1+i-2+intaxsdelta1,this%geom%n1-1)

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

    fout(1:ideb-1,:)=0
    fout(ifin+1:this%geom%n1,:)=0

    fout(:,1:jdeb-1)=0
    fout(:,jfin+1:this%geom%n2)=0

  end subroutine evaldep
end module splinenn_class
