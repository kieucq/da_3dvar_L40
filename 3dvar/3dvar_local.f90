!
! NOTE:     
!      This program performs a 3D var assimilation processes
!      for the barotropic model (2D). The 3Dvar process will be 
!      developed on each local patch to ease for the matrix inversion
!      by my own method. This program is not optimized in computation
!      and memory usage much. Need to re-think for real application.
!      It sounds interesting from discussions with Takemasa that this
!      simple approach has not been done before. Let's see how it works 
!
! HISTORY:  
!    - Jan 19, 2009: Created  
!
! REFERENCE:
!    - Kalnay, 2005: Book
!    - Krishnamurti 1998: Numerical technique
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
  program 3dvar
  use common
  use common_mtx
  implicit none
  integer, parameter          :: nx = 51,ny = 21     ! model domain size
  integer, parameter          :: nv=nx*ny            ! vector 1D
  integer, parameter          :: nxpatch = 10        ! number of local patch in x-direction
  integer, parameter          :: nypatch = 5         ! number of local patch in y-direction
  integer, parameter          :: no = nv             ! number of observations
  real, dimension(nx,ny)      :: pf                  ! background (forecast) 
  real, dimension(nx,ny)      :: pa                  ! analysis
  real, dimension(no)         :: po                  ! observation data 
  real, dimension(no)         :: olat,olon           ! observation lat and lon
  real, dimension(nv)         :: Xag                 ! analysis global vector
  real, dimension(nv,nv)      :: B                   ! background error cov matrix
  real, dimension(nv,nv)      :: Binv                ! B inverted
  real, dimension(no,nv)      :: H                   ! observation operator
  real, dimension(no,no)      :: R                   ! obs error cov matrix
  real, dimension(no,no)      :: Rinv                ! R inverted
  real, allocatable           :: Xf(:)               ! matrix of ensemble pertrubation
  real, allocatable           :: Xa(:)               ! matrix of ensemble analysis perturbation   
  real, dimension(no)         :: dx                  ! observation increments       
  integer, allocatable        :: ip(:),jp(:)         ! local indexing 
  real                        :: tem1(nx,ny)         ! temporate var
  real                        :: tem2(no)            ! temporate var
  real                        :: tem3(no)            ! temporate var
  real                        :: tem4(no)            ! temporate var
  character*100               :: ifile,ofile         ! I/O files
  integer                     :: i,j,i1,j1           ! indexing
  integer                     :: debug               ! debuging
  integer                     :: irec                ! output record
  integer                     :: nxp,nyp,nvp,nop     ! local ranges of each patch
  real                        :: bvar                ! background variance
  real                        :: ovar                ! observation variance
  real                        :: rscale              ! scale of background variance matrix
  integer                     :: iter                ! number of iterations after executing dfpmin
  real(kind=8)                :: fret                ! return value for the minimized cost function
  logical                     :: error               ! error return of dfpmin
  debug          = 2 
  irec           = 1
  bvar           = 1e7   
  ovar           = 1e7
  rscale         = 3.
  open(92,file='3dvar.dat',access='direct',form='unformatted',recl=nx*ny*4)
!
! reading observation data
!
  open(90,file='obs.dat')
  i              = 1
2 continue
  read(90,*,end=3)olon(i),olat(i),po(i)
  if (debug.eq.1) print*,i,olon(i),olat(i),po(i)
  i              = i + 1
  if (i.le.no) goto 2
  close(90)
  goto 4
3 print*,'There is not enough obs data as needed by no...stop'
  stop
4 continue
  close(90)
!
! quick check of obs, only for no=nx*ny
!
  if (no.eq.nx*ny) then
   call convert_vector_array1(po,no,tem1,nx,ny)
   write(92,rec=irec)((tem1(i,j),i=1,nx),j=1,ny)
   irec          = irec + 1
  endif
!
! reading forecast (background) data
!
  ifile          = 'fsc.dat'
  open(91,file=ifile)
  read(91,'(6e13.6)')((pf(i,j),i=1,nx),j=1,ny)
  write(92,rec=irec)((pf(i,j),i=1,nx),j=1,ny)
  irec           = irec + 1
  close(91)
!
! start to loop over the local patchs now.
!
  do i1          = 1,nxpatch
   do j1         = 1,nypatch
!
! define the ranges of the local patch first
!
    if (i1.lt.nxpatch) then
     nxp         = int(nx/nxpatch)
     nxp0        = nxp
    else
     nxp         = int(nx/nxpatch) + mod(nx,nxpatch)
    endif
    if (j1.lt.nypatch) then
     nyp         = int(ny/nypatch)
     nyp0        = nyp
    else
     nyp         = int(ny/nypatch) + mod(ny,nypatch)
    endif
    nvp          = nxp*nyp
    if (debug.eq.1) print*,'Ranges of local patch are',i1,j1,nxp,nyp,nvp
!
! assign the indices for the local patch now
!
    do i         = 1,nxp
     ip(i)       = (i1-1)*nxp + 1
     if (debug.eq.1) print*,i1,i,ip(i)
    enddo
    do j         = 1,nyp
     jp(i)       = (j1-1)*nyp + 1
     if (debug.eq.1) print*,j1,j,jp(j)
    enddo
!
! convert from 2d array to a 1d vector 
!
    allocate(Xf(nvp))
    call convert_array_vector(pf,nx,ny,Xf,nxp,nyp,nxp0,nyp0,nvp,i1,j1,ip,jp)
!
! prescribe the background error covariance matrix B
!
    allocate(B(nvp,nvp))
    call background_err_cov_mtx(B,nvp,nxp,nyp,bvar,rscale,ip,jp)
!
! find the number of observation points within the local patch
!
    nop          = 0
    do i         = 1,nxp-1
     do j        = 1,nyp-1
      do k       = 1,no
       if (ip(i).le.olon(k).and.olon(k).lt.ip(i+1).and. &
           jp(j).le.olat(k).and.olat(k).lt.jp(j+1)) then
        nop      = nop + 1
        tem2(nop)= po(k)
        tem3(nop)= olon(k)
        tem4(nop)= olat(k)
        goto 10
       endif    
      enddo
10    continue
     enddo
    enddo
    allocate(xop(nop),olatp(nop),olonp(nop))
    xop(1:nop)    = tem2(1:nop)
    olonp(1:nop)  = tem3(1:nop)
    olatp(1:nop)  = tem4(1:nop)
!
! compute the observational operator
!
    allocate(H(nop,nvp))
    call observation_operator(nxp,nyp,nxp0,nyp0,nop,nvp,olonp,olatp,H) 
!
! compute the observational error covariance matrix R
!
    call observational_err_cov_mtx(R,nop,ovar)
!
! Invert R and stored in Rinv
!
    call compute_Rinv(R,Rinv,nop)
!
! Invert B and stored in Binv
!
    allocate(Binv(nvp,nvp))
    call mtx_inv(nvp,B,Binv)
!
! Output a temp data file to store Binv and Rinv for
! the external functions fcost and dfcost
!
    open(81,file='temp.dat')
    write(81,*)nop
    do i          = 1,nop
     write(81,*)(Rinv(i,j),j=1,nop)
    enddo
    write(81,*)nvp
    do i          = 1,nvp
     write(81,*)(Binv(i,j),j=1,nvp)
    enddo
    close(81)
!
! first guess for analysis
!
    allocate(Xa(nvp))
    Xa             = Xf
!
! miniminze the cost fucntion bow
!
   call dfpmin(Xa,nvp,0.001,300,iter,fret,fcost,dfcost,error)
!
! map the local analysis to the full range 1d vector
!
    call mapping_local_global(Xa,nxp,nyp,nvp,ip,jp,nx,ny,nv,Xag)
   enddo
  enddo
!
! convert back from 1d vector to 2d array
!
  call convert_vector_array(Xag,nx,ny,nv,pa)
!
! output analysis files
!
  ofile='ana.dat'
  print*,'open output file is: ',ofile(1:20)
  open(12,file=ofile,status='unknown')
  write(12,'(6e13.6)')((pa(i,j,nfile),i=1,nx),j=1,ny)
  write(92,rec=irec)((pa(i,j,nfile),i=1,nx),j=1,ny)
  irec       = irec + 1
  print*,'3D-Var is finished safely...!'
  end

  subroutine observation_operator(nxp,nyp,nxp0,nyp0,nop,nvp,olonp,olatp,H)
  implicit none
  integer nxp,nyp,nop,nvp,i,j,m,n,k
  integer i1,j1,nxp0,nyp0
  real H(nop,nvp),olatp(no),olonp(no)
  io        = (i1-1)*nxp0 + 1
  jo        = (j1-1)*nyp0 + 1
  H         = 0.
  do i      = 1,nop
   m        = int(olonp(i)+0.001) - io 
   n        = int(olatp(i)+0.001) - jo 
   j        = (n-1)*nxp + m
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

  subroutine background_err_cov_mtx(B,nvp,nxp,nyp,bvar,rscale,ip,jp)
  implicit none
  integer nvp,i,j,m,n,nxp,nyp
  integer ip(nxp),jp(nyp)
  real B(nvp,nvp),rscale,radi,bvar
  do i      = 1,nvp
   do j     = 1,nvp
    m       = mod(j,nxp)
    n       = j/nxp + 1
    if (m.eq.0) then
     m      = nxp
     n      = n - 1
    endif
    radi    = sqrt((m-i)**2. + (n-i)**2.)
    B(i,j)  = bvar*exp(-(radi/rscale)**2)
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
  real R(no,no),obserr
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

     SUBROUTINE dfpmin (p,n,gtol,ITMAX,iter,fret,func,dfunc,error)

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
     INTEGER, INTENT(in)                       :: n,ITMAX
     INTEGER, INTENT(out)                      :: iter
     REAL(kind=8), DIMENSION(n), INTENT(inout) :: p
     REAL(kind=8), INTENT(in)                  :: gtol
     REAL(kind=8), INTENT(out)                 :: fret
     REAL(kind=8) func
     EXTERNAL dfunc,func
!
! Local variables:
!
     REAL(kind=8), PARAMETER                   :: STPMX=100.D0
     REAL(kind=8), PARAMETER                   :: EPS=3.0D-8
     REAL(kind=8), PARAMETER                   :: TOLX=4.0D0*EPS
     INTEGER i,its,j
     LOGICAL check,error
     REAL(kind=8) den,fac,fad,fae,fp,stpmax,xsum,sumdg,sumxi
     REAL(kind=8), DIMENSION(n)     :: dg,g,hdg,pnew,xi
     REAL(kind=8), DIMENSION(n,n)   :: hessin

     error=.false.
!
!U    USES dfunc,func,lnsrch
!
     fp = func(p)
     CALL dfunc(p,g)
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
       CALL lnsrch (n,p,fp,g,xi,pnew,fret,stpmax,check,func)
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
       CALL dfunc (p,g)
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
     SUBROUTINE lnsrch (n,xold,fold,g,p,x,f,stpmax,check,func)
     IMPLICIT NONE
     INTEGER, INTENT(in)                        :: n
     REAL(kind=8), INTENT(in)                   :: fold,stpmax
     REAL(kind=8), DIMENSION(n), INTENT(in)     :: xold,g
     REAL(kind=8), DIMENSION(n), INTENT(inout)  :: p
     REAL(kind=8), INTENT(out)                  :: f
     REAL(kind=8), DIMENSION(n), INTENT(out)    :: x
     LOGICAL, INTENT(out)                       :: check
     REAL(kind=8) func
     EXTERNAL func
     REAL(kind=8), PARAMETER                    :: ALF =1.0D-4
     REAL(kind=8), PARAMETER                    :: TOLX=1.0D-7
     REAL(kind=8) a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2
     REAL(kind=8) slope,xsum,tmplam
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
       f = func(x)
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

