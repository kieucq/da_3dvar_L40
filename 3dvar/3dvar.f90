!
! NOTE:     
!      This program performs a 3D var assimilation processes
!      for the Lorenz 40-var model, using the full cov matrix B
!
! HISTORY:  
!    - Feb,2 2009: Created  
!
! REFERENCE:
!    - Kalnay, 2005: Book
!
! AUTHOR: 
!      Chanh Q. Kieu, Research Associate
!      Dept. of Atmospheric and Oceanic science
!      Univ. of Maryland, College Park, MD
!      email: kieucq@atmos.umd.edu
!
! COPYRIGHT: (C) 2009
!
!=================================================================
!
  program variational
  use common
  use common_mtx
  implicit none
  integer, parameter          :: nx = 40             ! model domain size
  integer, parameter          :: nv = nx               ! vector 1D
  integer, parameter          :: nxpatch = 1         ! number of local patch in x-direction
  integer, parameter          :: no = nv             ! number of observations
  real, dimension(nx)         :: pf                  ! background (forecast) 
  real, dimension(nx)         :: pa                  ! analysis
  real, dimension(no)         :: po                  ! observation data 
  real, dimension(no)         :: olat,olon           ! observation lat and lon
  real, dimension(nv,nv)      :: B                   ! background error cov matrix
  real, dimension(nv,nv)      :: Binv                ! B inverted
  real, dimension(no,nv)      :: H                   ! observation operator
  real, dimension(no,no)      :: R                   ! obs error cov matrix
  real, dimension(no,no)      :: Rinv                ! R inverted
  real, dimension(no)         :: dx                  ! observation increments       
  real                        :: tem1(nx)            ! temporary var
  real                        :: tem2(nx,nx)         ! temp var
  character*100               :: ifile,ofile         ! I/O files
  integer                     :: i,j,k,i1,j1,k1      ! indexing
  integer                     :: debug               ! debuging
  integer                     :: irec                ! output record
  real                        :: bvar                ! background variance
  real                        :: ovar                ! observation variance
  real                        :: rscale              ! scale of background variance matrix
  integer                     :: iter                ! number of iterations after executing dfpmin
  real(kind=4)                :: fret                ! return value for the minimized cost function
  logical                     :: error               ! error return of dfpmin
  external                    :: fcost,dfcost
  debug          = 0 
  irec           = 1
  bvar           = 3e-1   
  ovar           = 1e-1
  rscale         = 3.
  open(92,file='3dvar.dat',access='direct',form='unformatted',recl=nx*4)
!
! reading observation data
!
  open(90,file='obs.dat')
  i              = 1
2 continue
  read(90,*,end=3)olon(i),po(i)
  if (debug.eq.1) then
   if (i.eq.1) print*,'checking the obs'
   print*,i,olon(i),po(i)
  endif
  i              = i + 1
  if (i.le.no) goto 2
  close(90)
  goto 4
3 print*,'There is not enough obs data as needed by no...stop'
  stop
4 continue
  close(90)
!
! quick check of obs data
!
  if (no.eq.nx) then
   write(92,rec=irec)(po(i),i=1,nx)
   irec          = irec + 1
  endif
!
! reading forecast (background) data
!
  ifile          = 'bgd.dat'
  open(91,file=ifile)
  read(91,*)(pf(i),i=1,nx)
  if (debug.eq.1) then
   print*,'Checking the background data'
   print*,(pf(i),i=1,nx)
  endif
  write(92,rec=irec)(pf(i),i=1,nx)
  irec           = irec + 1
  close(91)
!
! prescribe the background error covariance matrix B
!
  call background_err_cov_mtx(B,nv,bvar,rscale)
  if (debug.eq.1) then
   print*,'Checking the background cov mtx'
   do i           = 1,10
    write(0,'(10E12.3)')(B(i,j),j=1,10)
   enddo
  endif
!
! compute the observational operator
!
  call observation_operator(nv,no,olon,H) 
  if (debug.eq.1) then
   print*,'Checking the H mtx'
   do i           = 1,20
    write(0,'(20F6.3)')(H(i,j),j=1,20)
   enddo
  endif
!
! compute the observational error covariance matrix R
!
  call observational_err_cov_mtx(R,no,ovar)
!
! Invert R and stored in Rinv
!
  call compute_Rinv(R,Rinv,no)
!
! Invert B and stored in Binv
!
  call mtx_inv(nv,B,Binv)
  if (debug.eq.1) then
   print*,'Checking the inversed matrix'
   do i           = 1,10
    write(0,'(10E12.3)')(Binv(i,j),j=1,10)
   enddo
   print*,'Checking matrix inversion'
   tem2          = matmul(B,Binv)
   do i           = 1,20
    write(0,'(20F6.3)')(tem2(i,j),j=1,20)
   enddo
  endif
!
! first guess for analysis
!
  pa             = pf
!
! miniminze the cost fucntion bow
!
  call dfpmin(Binv,Rinv,H,pa,pf,po,nv,no,0.001,300,iter,fret,fcost,dfcost,error)
!
! output analysis files
!
  ofile='ana.dat'
  print*,'Open output file is: ',ofile(1:20)
  open(12,file=ofile,status='unknown')
  write(12,*)(pa(i),i=1,nx)
  write(92,rec=irec)(pa(i),i=1,nx)
  irec       = irec + 1
  print*,'3D-Var is finished safely...!'
  end

  subroutine observation_operator(nx,no,olon,H)
  implicit none
  integer nx,no,i,j,m,n,k,io,jo
  real H(no,nx),olon(no)
  H         = 0.
  do i      = 1,no
   j        = int(olon(i)+0.001)  
   H(i,j)   = 1.
  enddo
  return
  end

  subroutine localiztion_operator(nx,ny,no,nv,olon,olat,lopt)
  implicit none
  integer nx,ny,no,nv,i,j,m,n
  real lopt(nv,no),olat(no),olon(no),rscale,radi
  rscale    = 10.
  do i      = 1,nv
   do j     = 1,no
    m       = mod(i,nx)
    n       = i/nx + 1
    if (m.eq.0) then
     m      = nx
     n      = n - 1
    endif
    radi    = sqrt((olon(j)-m)**2. + (olat(j)-n)**2.)
    lopt(i,j) = exp(-(radi/rscale)**2)
   enddo
  enddo
  print*,'Checking Pat matrix'
  do i       = 1,20
   write(*,'(20F6.2)')(lopt(i,j),j=1,20)
  enddo
  return
  end

  subroutine background_err_cov_mtx(B,nx,bvar,rscale)
  implicit none
  integer i,j,m,n,nx
  real B(nx,nx),rscale,radi,bvar
  do i      = 1,nx
   do j     = 1,nx
    radi    = sqrt((j-i)**2. + (j-i)**2.)
    B(i,j)  = bvar*bvar*exp(-(radi/rscale)**2)
   enddo
  enddo
  return
  end

  subroutine convert_array_vector(p,nx,ny,X,nxp,nyp,nxp0,nyp0,nvp,i1,j1,ip,jp)
  implicit none
  integer nx,ny,ne,nvp,nxp,nyp,i1,j1,nxp0,nyp0
  real ip(nxp),jp(nyp)
  real p(nx,ny),X(nvp),tem
  integer i,j,k,m,n
!
! convert from an ensemble of 2D arrays to 1D vectors
!
  do i     = 1,nvp
!
! compute the local indices of grid points in the local patch first
!
   m       = mod(i,nxp)
   n       = i/nxp + 1
   if (m.eq.0) then
    m      = nxp
    n      = n - 1
   endif
   print*,'patch:',i1,j1,ip(m),jp(n)
!
! now compute the absolute indices w.r.t. to the global grid. Note 
! that the last patch is not necessarily the same size as the
! other patch. So, nxp0 must be used
!p
   m       = (i1-1)*nxp0 + m
   n       = (j1-1)*nyp0 + n
   print*,'patch:',i1,j1,m,n
   read*
   X(i)    = p(m,n)
!  X(i)    = p(ip(m),jp(n)) 
  enddo
  return
  end

  subroutine compute_Htilde(H,no,nv,Xf,ne,Ht)
  implicit none
  integer ne,nv,no
  real H(no,nv),Xf(nv,ne),Ht(no,ne)
  integer i,j,k,m,n
  Ht       = matmul(H,Xf)
  return
  end

  subroutine observational_err_cov_mtx(R,no,ovar)
  implicit none
  integer no,i
  real R(no,no),ovar
  R        = 0
  do i     = 1,no
   R(i,i)  = ovar**2
  enddo
  return
  end


  subroutine compute_Rinv(R,Rinv,no)
  implicit none
  integer no,i
  real R(no,no),Rinv(no,no)
  Rinv     = 0.
  do i     = 1,no
   Rinv(i,i) = 1/R(i,i)
  enddo
  return
  end
  
  subroutine obs_increment(H,no,nv,xfm,po,ne,obs_inc)
  implicit none
  integer no,ne,nv
  real H(no,nv),xfm(nv),po(no),obs_inc(nv)
  real tem
  integer i,j
  do i     = 1,no
   tem     = 0.
   do j    = 1,nv
    tem    = tem + H(i,j)*xfm(j)
   enddo
   obs_inc(i) = po(i) - tem
  enddo
  return
  end

  subroutine analysis_mean(K,lopt,nv,no,xfm,obs_inc,xam)
  implicit none
  integer no,nv
  real K(nv,no),xfm(nv),obs_inc(no),xam(nv),lopt(nv,no)
  integer i,j
  do i     = 1,nv
   xam(i)  = xfm(i)
   do j    = 1,no
    xam(i) = xam(i) + lopt(i,j)*K(i,j)*obs_inc(j)
   enddo
  enddo
  return
  end

  subroutine convert_vector_array(Xa,nx,ny,nv,pa)
  implicit none
  integer nv,nx,ny
  real Xa(nv),pa(nx,ny)
  integer i,j,k,m,n
  do i   = 1,nx
   do j  = 1,ny
    m    = (j-1)*nx + i
    pa(i,j)  = Xa(m)
   enddo
  enddo
  return
  end

  subroutine convert_vector_array1(Xa,nv,pa,nx,ny)
  implicit none
  integer nv,nx,ny
  real Xa(nv),pa(nx,ny)
  integer i,j,k,m,n
  do i       = 1,nx
   do j      = 1,ny
    m        = (j-1)*nx + i
    pa(i,j)  = Xa(m)
   enddo
  enddo
  return
  end

  subroutine mapping_local_global(Xa,nxp,nyp,nvp,ip,jp,nx,ny,nv,Xag)
  implicit none
  integer nv,nx,ny,nxp,nyp,nvp
  integer ip(nxp),jp(nyp),i,j,k,m,n
  real Xag(nv),Xa(nvp)
  do k       = 1,nvp
   m         = mod(k,nxp)
   n         = k/nxp + 1
   if (m.eq.0) then
    m        = nxp
    n        = n - 1
   endif
   i         = (jp(n)-1)*nx + ip(m)
   Xag(i)    = Xa(k)
  enddo
  return
  end

  real function fcost(p,pf,po,Binv,Rinv,H,n,no)
  implicit none
  integer n,no,i,j,debug1
  real H(no,n),Binv(n,n),Rinv(no,no)
  real pf(n),p(n),po(no),dpo(no)
  debug1     = 0
  fcost      = 0
!
! sum over background part
!
  do i       = 1,n
   do j      = 1,n
    fcost    = fcost + (p(i)-pf(i))*Binv(i,j)*(p(j)-pf(j))
   enddo
  enddo
!
! sum over obs increment
!
  dpo        = po - matmul(H,p)
  do i       = 1,no
   do j      = 1,no
    fcost    = fcost + dpo(i)*Rinv(i,j)*dpo(j)
   enddo
  enddo
  fcost      = fcost/2.
  if (debug1.eq.1) print*,'cost function returns',fcost
  return
  end

  subroutine dfcost(p,pf,po,g,Binv,Rinv,H,n,no)
  integer n,no,debug2
  real pf(n),p(n),po(no),dpo(no),g(n),dpb(n)
  real H(no,n),Binv(n,n),Rinv(no,no)
  debug2      = 0 
!
! define some departure
!
  dpb         = p - pf
  dpo         = po - matmul(H,pf)
!
! now gradient
!
  g           = matmul(Binv,dpb)                                   &
              + matmul(transpose(H),matmul(Rinv,matmul(H,dpb)))    &
              - matmul(transpose(H),matmul(Rinv,dpo))  
  if (debug2.eq.1) print*,'Grad cost function return',g(1:n)
  return
  end

     SUBROUTINE dfpmin (Binv,Rinv,H,p,pf,po,n,no,gtol,ITMAX,iter,fret,func,dfunc,error)

! ... Given a starting point p that is a vector of dimension n, the
! ... Broyden-Fletcher-Goldfarb-Shanno (BFGS) variant of the
! ... Davidson-Fletcher-Powell (DFP) minimization algorithm is performed
! ... on a function func, using its gradient as calculated by the routine
! ... dfunc. The convergence requirement on zeroing the gradient is input
! ... as gtol. Returned quantities are p (the location of thr minimum),
! ... iter (the number of iterations that were performed), and fret (the
! ... minimum value of the function). The routine lnsrch is called to
! ... perform approximate line minimizations.
! ... Parameters: ITMAX is the maximum allowed number of iterations;
! ... STPMX is the scaled maximum step length allowed in line searches;
! ... EPS should be the machine precision; TOLX is the convergence values
! ... on x values.

     IMPLICIT NONE
     INTEGER, INTENT(in)                       :: n,ITMAX,no
     INTEGER, INTENT(out)                      :: iter
     REAL(kind=4), DIMENSION(n), INTENT(inout) :: p
     REAL(kind=4), INTENT(in)                  :: Binv(n,n),Rinv(no,no),H(no,n)
     REAL(kind=4), INTENT(in)                  :: gtol,po(no),pf(n)
     REAL(kind=4), INTENT(out)                 :: fret
     REAL(kind=4) func
     EXTERNAL dfunc,func
!
! Local variables:
!
     REAL(kind=4), PARAMETER                   :: STPMX=100.D0
     REAL(kind=4), PARAMETER                   :: EPS=3.0D-8
     REAL(kind=4), PARAMETER                   :: TOLX=4.0D0*EPS
     INTEGER i,its,j
     LOGICAL check,error
     REAL(kind=4) den,fac,fad,fae,fp,stpmax,xsum,sumdg,sumxi
     REAL(kind=4), DIMENSION(n)     :: dg,g,hdg,pnew,xi
     REAL(kind=4), DIMENSION(n,n)   :: hessin

     error=.false.
!
!U    USES dfunc,func,lnsrch
!
     fp = func(p,pf,po,Binv,Rinv,H,n,no)
     CALL dfunc(p,pf,po,g,Binv,Rinv,H,n,no)
     hessin(:,:) = 0.0D0
     DO i=1,n
       hessin(i,i) = 1.0D0
     ENDDO
     xi(:) = -g(:)
     xsum = DOT_PRODUCT(p,p)
     stpmax = STPMX*MAX(SQRT(xsum),DBLE(n))

     DO its=1,ITMAX
       iter=its

! ... The new evaluation occurs in lnsrch; save the function value in fp for
! ... the next line search. It is usually safer to ignore the value of check.
! ...
       CALL lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func,Binv,Rinv,H,no,pf,po)
       fp=fret
       xi(:) = pnew(:) - p(:)                   ! Update the line direction
       p(:)  = pnew(:)                          ! Update the current point.
!
! ... Test for convergence on Delta x:
! ...
       IF (MAXVAL(ABS(xi)/MAX(ABS(p),1.0D0)).LT.TOLX) THEN
         WRITE(*,*) 'Too small space increment'
         error=.true.
         RETURN
       ENDIF

       dg(:) = g(:)                             ! Save the old gradient
       CALL dfunc (p,pf,po,g,Binv,Rinv,H,n,no)
       den = MAX(fret,1.0D0)
!
! ... Test for convergence on zero gradient:
! ...
       IF (MAXVAL(ABS(g)*MAX(ABS(p),1.0D0)/den).LT.gtol) THEN
        RETURN
       ENDIF

       dg(:)  = g(:) - dg(:)                     ! Difference of gradients
       hdg(:) = MATMUL(hessin,dg)
       fac    = DOT_PRODUCT(dg,xi)
       fae    = DOT_PRODUCT(dg,hdg)
       sumdg  = DOT_PRODUCT(dg,dg)
       sumxi  = DOT_PRODUCT(xi,xi)

       IF (fac**2.GT.EPS*sumdg*sumxi) THEN      ! Skip update if not positive enough
         fac   = 1.d0/fac
         fad   = 1.d0/fae
         dg(:) = fac*xi(:) - fad*hdg(:)         ! What makes BFGS differ from DFP
         DO i=1,n
         DO j=1,n
           hessin(i,j) = hessin(i,j)       + &  ! The BFGS updating
                         fac*xi(i)*xi(j)   - &  ! formula
                         fad*hdg(i)*hdg(j) + &
                         fae*dg(i)*dg(j)
         ENDDO
         ENDDO
       ENDIF
       xi(:) = -MATMUL(hessin,g)                ! Next direction to go
     ENDDO
     STOP 'too many iterations in dfpmin'
     END SUBROUTINE dfpmin
!
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
!              
     SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func,Binv,Rinv,H,no,pf,po)
     IMPLICIT NONE
     INTEGER, INTENT(in)                        :: n,no
     REAL(kind=4), INTENT(in)                   :: fold,stpmax
     REAL(kind=4), DIMENSION(n), INTENT(in)     :: xold,g,pf(n),po(no)
     REAL(kind=4), DIMENSION(n), INTENT(inout)  :: p
     REAL(kind=4), INTENT(in)                   :: Binv(n,n),Rinv(no,no),H(no,n)
     REAL(kind=4), INTENT(out)                  :: f
     REAL(kind=4), DIMENSION(n), INTENT(out)    :: x
     LOGICAL, INTENT(out)                       :: check
     REAL(kind=4) func
     EXTERNAL func
     REAL(kind=4), PARAMETER                    :: ALF =1.0D-4
     REAL(kind=4), PARAMETER                    :: TOLX=1.0D-7
     REAL(kind=4) a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2
     REAL(kind=4) slope,xsum,tmplam
     check = .false.
     xsum  = SQRT(DOT_PRODUCT(p,p))
     IF (xsum.GT.stpmax) p(:) = p(:)*stpmax/xsum
     slope = DOT_PRODUCT(g,p)
     IF (slope >= 0.0) THEN
!
! STOP 'Roundoff problem in lnsrch'
!
       WRITE(*,*) '!!!!!!!!!!!!!!! ROUNDOFF PROBLEM IN LNSRCH'
      RETURN
     ENDIF
     alamin = TOLX/MAXVAL(ABS(p(:))/MAX(ABS(xold(:)),1.0D0))       ! Lambda_min
     alam   = 1.0D0   ! Always try full Newton step first
     DO
       x(:) = xold(:) + alam*p(:)
       f = func(x,pf,po,Binv,Rinv,H,n,no)
       IF (alam.LT.alamin) THEN                           ! Convergence on Delta x.
         x(:) = xold(:)
         check = .true.
         RETURN
       ELSE IF (f.LE.fold+ALF*alam*slope) THEN            ! Enough function decrease
         RETURN
       ELSE
         IF (alam.EQ.1.0D0) THEN                          ! First time
           tmplam = - 0.5D0*slope/(f-fold-slope)
         ELSE                                             ! Subsequent times
           rhs1 = f  - fold - alam*slope
           rhs2 = f2 - fold - alam2*slope
           a    = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
           b    = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
           IF (a.EQ.0.D0) THEN
             tmplam = -0.5D0*slope/b
           ELSE
             disc   = b*b - 3.D0*a*slope
             IF (disc.LT.0.0D0) THEN
               tmplam = 0.5D0*alam
             ELSE IF (b.LE.0.0D0) THEN
               tmplam = (-b+SQRT(disc))/(3.D0*a)
             ELSE
               tmplam = -slope/(b+SQRT(disc))
             ENDIF
           ENDIF
           IF (tmplam.GT.0.5D0*alam) tmplam = 0.5D0*alam
         ENDIF
       ENDIF
       alam2 = alam
       f2    = f
       alam  = MAX(tmplam,0.1D0*alam)
     END DO
    END SUBROUTINE lnsrch

