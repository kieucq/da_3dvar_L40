!
! Note:
!         This program is for performing an analysis of the output
!         from the ctl, truth, and assimilation runs.
!
! History: Created Feb 6, 2009
!
! Author: Chanh Q. Kieu
!
!===================================================================
  PROGRAM analysis
  IMPLICIT NONE
  INTEGER,PARAMETER :: n = 40 
  INTEGER,PARAMETER :: ntime = 21  
  REAL,PARAMETER    :: obs_err = 0.05
  REAL              :: xc(n,ntime),xt(n,ntime),xb(n,ntime),xo(n,ntime)
  REAL              :: rmsc(ntime),rmsb(ntime),rmso(ntime),tem
  CHARACTER*50      :: afile,cfile,tfile,ofile
  INTEGER           :: id,i,j,k,itime,irec
  OPEN(91,file='ana.dat',access='direct',form='unformatted',recl=n*4)
  OPEN(92,file='ana.txt')
  tfile     = 'truth0000.dat'
  cfile     = './ctl0000.dat'
  afile     = './ana0000.dat'
  ofile     = './obs0000.dat'
  itime     = 1
  id        = 0
  irec      = 1
19 continue
!
! reading the truth
!
  IF (id.lt.10) THEN
   WRITE(tfile(9:9),'(1I1)')id
   WRITE(cfile(9:9),'(1I1)')id
   WRITE(afile(9:9),'(1I1)')id
   WRITE(ofile(9:9),'(1I1)')id
  ELSEIF (id.lt.100) THEN
   WRITE(tfile(8:9),'(1I2)')id
   WRITE(cfile(8:9),'(1I2)')id
   WRITE(afile(8:9),'(1I2)')id
   WRITE(ofile(8:9),'(1I2)')id
  ELSEIF (id.lt.1000) THEN
   WRITE(tfile(7:9),'(1I3)')id
   WRITE(cfile(7:9),'(1I3)')id
   WRITE(afile(7:9),'(1I3)')id
   WRITE(ofile(7:9),'(1I3)')id
  ELSE
   WRITE(tfile(6:9),'(1I4)')id
   WRITE(cfile(6:9),'(1I4)')id
   WRITE(afile(6:9),'(1I4)')id
   WRITE(ofile(6:9),'(1I4)')id
  ENDIF
  PRINT*,'Open truth file is:  ',tfile(1:30)
  PRINT*,'Open control file is:  ',cfile(1:30)
  PRINT*,'Open background file is:  ',afile(1:30)
  PRINT*,'Open observation file is:  ',ofile(1:30)
  OPEN(71,file=tfile,status='old')
  OPEN(72,file=cfile,status='old')
  OPEN(73,file=afile,status='old')
  OPEN(74,file=ofile,status='old')
  READ(71,*)(xt(i,itime),i=1,n)
  READ(72,*)(xc(i,itime),i=1,n)
  READ(73,*)(xb(i,itime),i=1,n)
  DO i     = 1,n
   READ(74,*)tem,xo(i,itime)
  ENDDO 
  CLOSE(71)
  CLOSE(72)
  CLOSE(73)
  CLOSE(74)

  WRITE(91,rec=irec)(xt(i,itime),i=1,n)
  irec     = irec + 1
  WRITE(91,rec=irec)(xc(i,itime),i=1,n)
  irec     = irec + 1
  WRITE(91,rec=irec)(xb(i,itime),i=1,n)
  irec     = irec + 1
  WRITE(91,rec=irec)(xo(i,itime),i=1,n)
  irec     = irec + 1
!
! Compute the stardard error devidation
!
  rmso(itime)  = 0.
  rmsc(itime)  = 0.
  rmsb(itime)  = 0.
  DO i         = 1,n
   rmsc(itime) = rmsc(itime) + (xc(i,itime)-xt(i,itime))**2
   rmsb(itime) = rmsb(itime) + (xb(i,itime)-xt(i,itime))**2
   rmso(itime) = rmso(itime) + (xo(i,itime)-xt(i,itime))**2
  ENDDO
  rmsc(itime)  = sqrt(rmsc(itime)/n)
  rmsb(itime)  = sqrt(rmsb(itime)/n)
  rmso(itime)  = sqrt(rmso(itime)/n)
  PRINT*,'Checking error now',itime,rmsc(itime),rmsb(itime),rmso(itime)
  if (itime.eq.1) WRITE(92,'(A5,10A16)')'time','rmse ctl','rmse ana','rmse obs'
  WRITE(92,'(I5,10F16.4)')itime,rmsc(itime),rmsb(itime),rmso(itime)
!
! advance and loop now
!
  itime    = itime + 1
  id       = (itime-1)*50 
  IF (itime.le.ntime) GOTO 19
  PRINT*,'Program ends perfectly'
  END
