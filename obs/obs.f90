!
! This program is for creating an ensemble of  
! observation that are perturbed from the truth
!
  PROGRAM obs
  IMPLICIT NONE
  INTEGER,PARAMETER :: n = 40 
  INTEGER,PARAMETER :: ntime = 21  
  REAL,PARAMETER    :: obs_err = 0.05
  REAL              :: x(n),xp(n)
  REAL(8)           :: rnd(n),xr(n)
  CHARACTER*50      :: ofile,ifile
  INTEGER           :: id,i,j,k,itime,irec
  OPEN(91,file='obs.dat',access='direct',form='unformatted',recl=n*4)
  ifile     = './truth0000.dat'
  ofile     = './obs0000.dat'
  itime     = 1
  id        = 0
  irec      = 1
19 continue
!
! reading the truth
!
  IF (id.lt.10) THEN
   WRITE(ifile(11:11),'(1I1)')id
  ELSEIF (id.lt.100) THEN
   WRITE(ifile(10:11),'(1I2)')id
  ELSEIF (id.lt.1000) THEN
   WRITE(ifile(9:11),'(1I3)')id
  ELSE
   WRITE(ifile(8:11),'(1I4)')id
  ENDIF
  PRINT*,'Open file is:  ',ifile(1:30)
  OPEN(71,file=ifile,status='old')
  READ(71,*)(x(i),i=1,n)
  CLOSE(71)
  WRITE(91,rec=irec)(x(i),i=1,n)
  irec     = irec + 1
  CALL com_randn(n,rnd)
  xr       = rnd*obs_err
  xp       = x + xr
  WRITE(91,rec=irec)(xp(i),i=1,n)
  irec     = irec + 1
!
! Output the observation
!
  IF (id.lt.10) THEN
   WRITE(ofile(9:9),'(1I1)')id
  ELSEIF (id.lt.100) THEN
   WRITE(ofile(8:9),'(1I2)')id
  ELSEIF (id.lt.1000) THEN
   WRITE(ofile(7:9),'(1I3)')id
  ELSE
   WRITE(ofile(6:9),'(1I4)')id
  ENDIF
  PRINT*,'Output file is:  ',ofile(1:30)
  OPEN(72,file=ofile,status='unknown')
  DO i     = 1,n
!    IF (mod(i,3).eq.0.and.mod(j,2).eq.0) WRITE(72,'(2I5,E12.4)')i,j,xp(i,j) 
    WRITE(72,*)i,xp(i)
  ENDDO
  CLOSE(72)
  itime    = itime + 1
  id       = (itime-1)*50 
  IF (itime.le.ntime) GOTO 19
  PRINT*,'Program ends perfectly'
  END 


  SUBROUTINE com_randn(ndim,var)
  USE mt19937
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  REAL(8),INTENT(OUT) :: var(1:ndim)
  REAL(8) :: rnd(2),pi
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.
  pi = 4.*atan(1.)

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_genrand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF
  RETURN
END SUBROUTINE com_randn
