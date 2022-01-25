      SUBROUTINE inverse_matrix_gauss(N,C,B)
      INTEGER N
      REAL C(N,N),B(N,N),t1
      REAL A(N,2*N)
!
! Start algorithm below
!
      DO i       = 1,N
       DO j      = 1,2*N
        IF (j.le.N) THEN
         A(i,j)  = C(i,j)
        ELSE
         A(i,j)   = 0.
        ENDIF
        IF (j.eq.N+i) A(i,j) = 1.
       ENDDO
      ENDDO
      DO i       = 1,N
       t1        = 0.
       PRINT*,'Inverting at row',i
       DO j      = i,N
        IF (t1.lt.abs(A(j,i))) THEN
         t1      = abs(A(j,i))
         L       = j 
        ENDIF
       ENDDO
       IF (t1.eq.0) GOTO 10
       IF (L.ne.i) THEN
        DO j     = i,2*N
         t1      = A(i,j)
         A(i,j)  = A(L,j)
         A(L,j)  = t1   
        ENDDO
       ENDIF
       t1        = A(i,i)
       DO j      = i,2*N
        A(i,j)   = A(i,j)/t1
       ENDDO
       DO j      = 1,N
       IF (j.ne.i) THEN
         DO k    = i+1,2*N
          A(j,k) = A(j,k) - A(j,i)*A(i,k) 
         ENDDO
       ENDIF
      ENDDO
      ENDDO
      GOTO 11
10    PRINT*,'Degenerate matrix, no inverse matrix..stop'
      STOP
11    B(1:N,1:N) = A(1:N,(N+1):2*N)
      RETURN
      END


