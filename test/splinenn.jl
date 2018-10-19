mutable struct SplineNN

    geom       :: Geometry
    ax         :: Vector{Float64} 
    ay         :: Vector{Float64}
    axd        :: Vector{Float64}
    ayd        :: Vector{Float64}
    axod       :: Vector{Float64}
    ayod       :: Vector{Float64}
    aym1gamma1 :: Vector{Float64}
    aym1gamma2 :: Vector{Float64}
    axm1gamma1 :: Vector{Float64}
    axm1gamma2 :: Vector{Float64}
    coef       :: Vector{Float64}
    bcoef      :: Vector{Float64}

    function SplineNN(geom::Geometry)

        n1   = geom.n1
        n2   = geom.n2
        n1p1 = n1+1
        n1p2 = n1+2
        n2p1 = n2+1
        n2p2 = n2+2
    
        axm1gamma1 = zeros(Float64, n1)
        axm1gamma2 = zeros(Float64, n1)
        aym1gamma1 = zeros(Float64, n2)
        axm1gamma2 = zeros(Float64, n2)
        axd        = 4ones(Float64, n1)
        axod       =  ones(Float64, n1-1)
        ayd        = 4ones(Float64, n2)
        ayod       =  ones(Float64, n2-1)

        coef       = zeros(Float64,(n1p2,n2p2))
        bcoef      = zeros(Float64,(n2p2,n1p2))
    
        # factorize matrices  Ax and Ay
        call dpttrf(n1,this.axd,this.axod,err)
        call dpttrf(n2,this.ayd,this.ayod,err)

        axm1gamma2[1]  = 1
        axm1gamma1[n1] = 1

        call dpttrs(n1,1,this.axd,this.axod,this.axm1gamma1,n1,err)
        call dpttrs(n1,1,this.axd,this.axod,this.axm1gamma2,n1,err)
    
        # compute Ay-1.gamma
        aym1gamma2[1]  = 1
        aym1gamma1[n2] = 1

        call dpttrs(n2,1,this.ayd,this.ayod,this.aym1gamma1,n2,err)
        call dpttrs(n2,1,this.ayd,this.ayod,this.aym1gamma2,n2,err)
    
        aa1x = 3 /  geom.delta1
        aa1y = 3 /  geom.delta2
        aa2x = 6 / (geom.delta1*geom.delta1)
        aa2y = 6 / (geom.delta2*geom.delta2)

        # assemblage de la matrice 2x2 pour la spline dans la direction Ox

        a1x = aa2x*(1. - axm1gamma1[n1-1] + 2*axm1gamma1[n1])
        a2x = aa2x*(   - axm1gamma2[n1-1] + 2*axm1gamma2[n1])
        a3x = aa2x*(   2*axm1gamma1[1]    -   axm1gamma1[2])
        a4x = aa2x*(1. + 2*this.axm1gamma2[1] - this.axm1gamma2[2])

        # assemblage de la matrice 2x2 pour spline naturels (direction Oy)

        this.a1y = aa2y*(1. -   this.aym1gamma1(n2-1) + 2*this.aym1gamma1(n2))
        this.a2y = aa2y*(   -   this.aym1gamma2(n2-1) + 2*this.aym1gamma2(n2))
        this.a3y = aa2y*(     2*this.aym1gamma1(1)    -   this.aym1gamma1(2))
        this.a4y = aa2y*(1. + 2*this.aym1gamma2(1)    -   this.aym1gamma2(2))
    
        new( geom, ax , ay, axd , ayd , axod , ayod , 
             aym1gamma1, aym1gamma2, axm1gamma1, axm1gamma2,
             coef, bcoef )

end 

function interpole_splinenn!(this :: SplineNN, 
                             fin  :: Array{Float64,2},
                             fout :: Array{Float64,2},
                             x    :: Float64,
                             y    :: Float64) 

    nat_x!(this,fin,ierr)
    nat_y!(this,ierr)
    evaltab!(this,x,y,fout)     

end

"""
 interpolation par spline periodique dans les deux directions.
 Les points d'interpolation sont definis grace a depx et depy
 qui definissent le deplacement par rapport au maillage.
 - f contient les valeurs de la fonction de distribution
 - depx et depy : deplacements par rapport au maillage
     des points dans les quels on veut evaluer la spline.
"""
function interpole_splinenndep(this :: SplineNN,
                               f    :: Array{Float64,2},
                               depx :: Float64,
                               depy :: Float64,
                               aff) 


    nat_x!(this,f,ierr)
    nat_y!(this,ierr)
    
    evaldep!(this,depx,depy,f)  

end

"""
    calcul des "natural splines"
""" 
function nat_x( this :: SplineNN, 
                gtau :: Array{Float64,2} )

    axm1f = zeros( Float64, (this.geom.n1,this.geom.n2))

    n1   =  this.geom.n1
    n2   =  this.geom.n2
    n1p2 =  n1+2
    n1p1 =  n1+1
    det  =  this.a1x*this.a4x - this.a2x*this.a3x

    for j=1:n2
       for i=1:n1
          axm1f[i,j] = 6*gtau[i,j]
       end
    end

    call dpttrs(n1,n2,this.axd,this.axod,axm1f,n1,err)

    @assert (err == 0) 

    for j=1:n2
       # assemblage du second membre du systeme 2x2 
       gamma1 =  (6.0/(this.geom.delta1)**2)*(-axm1f[this.geom.n1-1,j] &
            + 2*axm1f[this.geom.n1,j])
       gamma2 = (6.0/(this.geom.delta1)**2)*(2*axm1f[1,j] - axm1f[2,j])

       coefnp2 = (gamma1*this.a4x - gamma2*this.a2x)/det
       coef1 = (gamma2*this.a1x - gamma1*this.a3x)/det
       this.bcoef[j,n1p2]=coefnp2
       this.bcoef[j,1]=coef1

       for  i=2,n1p1
          this.bcoef[j,i]= axm1f[i-1,j] &
               - this.axm1gamma1[i-1]*coefnp2 &
               - this.axm1gamma2[i-1]*coef1
       end
    end
end 

function nat_y(this::SplineNN)

    aym1f = zeros(Float64, (this.geom.n2,this.geom.n1+2))

    n2   = this.geom.n2
    n1p2 = this.geom.n1+2
    n2p2 = this.geom.n2+2
    det  = this.a1y*this.a4y - this.a2y*this.a3y

    for  i=1,n1p2
       for j=1,this.geom.n2
          aym1f[j,i] = 6.*this.bcoef[j,i]
       end
    end

    call dpttrs(n2,n1p2,this.ayd,this.ayod,aym1f,n2,err)

    @assert ( err == 0 )

    ! resolution du syteme lineaire 2x2
    for i=1:n1p2
       gamma1 =  (6.0/(this.geom.delta2)**2)*(-aym1f[this.geom.n2-1,i] &
            + 2*aym1f[this.geom.n2,i])
       gamma2 = (6.0/(this.geom.delta2)**2)*(2*aym1f[1,i] - aym1f[2,i])

       coefnp2 = (gamma1*this.a4y - gamma2*this.a2y)/det
       coef1 = (gamma2*this.a1y - gamma1*this.a3y)/det
       this.coef[i,n2p2] = coefnp2
       this.coef[i,1] = coef1
       for  j=2,this.geom.n2 + 1
          this.coef[i,j]= aym1f[j-1,i] 
               - this.aym1gamma1[j-1]*coefnp2 
               - this.aym1gamma2[j-1]*coef1
       end
    end
end 

function evaltab(this :: SplineNN,
                 xd   :: Array{Float64,2},
                 yd   :: Array{Float64,2},
                 fout :: Array{Float64,2})



    delta1x    = this.geom.delta1*this.geom.delta1
    delta1xx   = delta1x*this.geom.delta1
    delta1xx6  = 1./(6.*delta1xx)
    delta2y    = this.geom.delta2*this.geom.delta2
    delta2yy   = delta2y*this.geom.delta2
    delta2yy6  = 1./(6.*delta2yy)
    idelta1    = 1./this.geom.delta1
    idelta2    = 1./this.geom.delta2

    for j=2:this.geom.n2-1
       for i=2:this.geom.n1-1

          i1  = (xd[i,j]-this.geom.x0)*idelta1
          j1  = (yd[i,j]-this.geom.y0)*idelta2

          xdp1   = this.geom.xgrid[i1+2]-xd[i,j]
          bvalx1 = xdp1*xdp1*xdp1
          bvalx2 = delta1xx+3.*delta1x*xdp1
                 + 3.*this.geom.delta1*xdp1*xdp1-3.*xdp1*xdp1*xdp1
          xd1=xd[i,j]-this.geom.xgrid[i1+1]
          bvalx3=delta1xx+3.*delta1x*xd1
                 +3.*this.geom.delta1*xd1*xd1-3.*xd1*xd1*xd1
          bvalx4=xd1*xd1*xd1
          ydp1=this.geom.ygrid[j1+2]-yd[i,j]
          bvaly1=ydp1*ydp1*ydp1
          bvaly2=delta2yy+3.*ydp1*(delta2y+ydp1*(this.geom.delta2-ydp1))          
          yd1=yd[i,j]-this.geom.ygrid[j1+1]
          bvaly3=delta2yy+3.*yd1*(delta2y+yd1*(this.geom.delta2-yd1))
          bvaly4=yd1*yd1*yd1

          sval  = 0.
          sval1 = this.coef[i1+1,j1+1]*bvaly1
          sval1 = sval1+this.coef[i1+1,j1+2]*bvaly2
          sval1 = sval1+this.coef[i1+1,j1+3]*bvaly3
          sval1 = sval1+this.coef[i1+1,j1+4]*bvaly4        

          sval  =  sval+sval1*bvalx1

          sval2 = this.coef[i1+2,j1+1]*bvaly1
          sval2 = sval2+this.coef[i1+2,j1+2]*bvaly2
          sval2 = sval2+this.coef[i1+2,j1+3]*bvaly3
          sval2 = sval2+this.coef[i1+2,j1+4]*bvaly4

          sval  = sval+sval2*bvalx2

          sval3 = this.coef[i1+3,j1+1]*bvaly1
          sval3 = sval3+this.coef[i1+3,j1+2]*bvaly2
          sval3 = sval3+this.coef[i1+3,j1+3]*bvaly3
          sval3 = sval3+this.coef[i1+3,j1+4]*bvaly4

          sval  = sval+sval3*bvalx3

          sval4 = this.coef[i1+4,j1+1]*bvaly1
          sval4 = sval4+this.coef[i1+4,j1+2]*bvaly2
          sval4 = sval4+this.coef[i1+4,j1+3]*bvaly3
          sval4 = sval4+this.coef[i1+4,j1+4]*bvaly4

          sval  = sval+sval4*bvalx4

          fout[i,j] = delta1xx6*delta2yy6*sval
       end
    end

    fout[1,:]=0
    fout[this.geom.n1,:]=0 

    fout[:,1]=0
    fout[:,this.geom.n2]=0 

end

function evaldep(this::SplineNN,alphax,alphay,fout::Array{Float64,2})

    delta1x   = this.geom.delta1*this.geom.delta1
    delta1xx  = delta1x*this.geom.delta1
    delta1xx6 = 1./(6.*delta1xx)
    
    delta2y   = this.geom.delta2*this.geom.delta2
    delta2yy  = delta2y*this.geom.delta2
    delta2yy6 = 1./(6.*delta2yy)

    if (alphax > 0) 
       intaxsdelta1=trunc(-alphax/this.geom.delta1+epsilon)-1 
       ideb = trunc(alphax/this.geom.delta1)+2
       ifin = this.geom.n1 - 1
    else
       intaxsdelta1=trunc(-alphax/this.geom.delta1)
       ideb = 2
       ifin = -trunc(-alphax/this.geom.delta1)+this.geom.n1-1
    end

    if (alphay>0) 
       intaysdelta2=trunc(-alphay/this.geom.delta2+epsilon)-1
       jdeb = trunc(alphay/this.geom.delta2)+2
       jfin = this.geom.n2 - 1
    else
       intaysdelta2=trunc(-alphay/this.geom.delta2)
       jdeb = 2
       jfin = -trunc(-alphay/this.geom.delta2)+this.geom.n2-1
    end

    xd1    = -alphax-intaxsdelta1*this.geom.delta1
    xdp1   =  this.geom.delta1-xd1
    yd1    = -alphay-intaysdelta2*this.geom.delta2
    ydp1   = this.geom.delta2-yd1
    bvalx1 = xdp1*xdp1*xdp1
    bvalx2 = delta1xx+3.*delta1x*xdp1+3.*this.geom.delta1*xdp1*xdp1-3.*xdp1*xdp1*xdp1
    bvalx3 = delta1xx+3.*delta1x*xd1+3.*this.geom.delta1*xd1*xd1-3.*xd1*xd1*xd1
    bvalx4 = xd1*xd1*xd1
    bvaly1 = ydp1*ydp1*ydp1
    bvaly2 = delta2yy+3.*ydp1*(delta2y+ydp1*(this.geom.delta2-ydp1))          
    bvaly3 = delta2yy+3.*yd1*(delta2y+yd1*(this.geom.delta2-yd1))
    bvaly4 = yd1*yd1*yd1

    for j=jdeb:jfin
       j1=j-1+intaysdelta2
       for i=ideb:ifin
          i1=i-1+intaxsdelta1

          fout[i,j] = delta1xx6*delta2yy6* ( 
               bvalx1*( this.coef(i1+1,j1+1)*bvaly1 
               +this.coef(i1+1,j1+2)*bvaly2 
               +this.coef(i1+1,j1+3)*bvaly3 
               +this.coef(i1+1,j1+4)*bvaly4) 
               + bvalx2* (this.coef(i1+2,j1+1)*bvaly1 
               +this.coef(i1+2,j1+2)*bvaly2 
               +this.coef(i1+2,j1+3)*bvaly3 
               +this.coef(i1+2,j1+4)*bvaly4) 
               + bvalx3* (this.coef(i1+3,j1+1)*bvaly1 
               +this.coef(i1+3,j1+2)*bvaly2 
               +this.coef(i1+3,j1+3)*bvaly3 
               +this.coef(i1+3,j1+4)*bvaly4) 
               + bvalx4* (this.coef(i1+4,j1+1)*bvaly1 
               +this.coef(i1+4,j1+2)*bvaly2 
               +this.coef(i1+4,j1+3)*bvaly3 
               +this.coef(i1+4,j1+4)*bvaly4))
       end
    end

    fout[1:ideb-1,:]=0
    fout[ifin+1:this.geom.n1,:]=0

    fout[:,1:jdeb-1]=0
    fout[:,jfin+1:this.geom.n2]=0

end 
