module utils


contains

  subroutine c_set(n,a,b)
    integer i,n
    character*(*)a(n),b
    do i=1,n
       a(i)=b
    enddo
    return
  end subroutine c_set

  subroutine r_set(n,a,b)
    integer i,n
    real*8 a(n),b
    do i=1,n
       a(i)=b
    enddo
    return
  end subroutine r_set

  subroutine i_set(n,a,b)
    integer i,n,a(n),b
    do i=1,n
       a(i)=b
    enddo
    return
  end subroutine i_set

  subroutine dcopy(n,a,ia,b,ib)
    integer ia,ib,n,i
    real*8 a(n),b(n)
    do i=1,n
       b(i)=a(i)
    enddo
    return
  end subroutine dcopy

      SUBROUTINE TQLI(D,E,N,NP,Z)
      implicit real*8(a-h,o-z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      DO 11 I=2,N
        E(I-1)=E(I)
11    CONTINUE
      E(N)=0.d0
      DO 15 L=1,N
        ITER=0
1       DO 12 M=L,N-1
          DD=ABS(D(M))+ABS(D(M+1))
          IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12      CONTINUE
        M=N
2       IF(M.NE.L)THEN
          IF(ITER.EQ.30)PAUSE 'too many iterations'
          ITER=ITER+1
          G=(D(L+1)-D(L))/(2.d0*E(L))
          R=SQRT(G**2+1.d0)
          G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
          S=1.d0
          C=1.d0
          P=0.d0
          DO 14 I=M-1,L,-1
            F=S*E(I)
            B=C*E(I)
            IF(ABS(F).GE.ABS(G))THEN
              C=G/F
              R=SQRT(C**2+1.d0)
              E(I+1)=F*R
              S=1.d0/R
              C=C*S
            ELSE
              S=F/G
              R=SQRT(S**2+1.d0)
              E(I+1)=G*R
              C=1.d0/R
              S=S*C
            ENDIF
            G=D(I+1)-P
            R=(D(I)-G)*S+2.d0*C*B
            P=S*R
            D(I+1)=G+P
            G=C*R-B
!     Omit lines from here ...
            DO 13 K=1,N
              F=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*F
              Z(K,I)=C*Z(K,I)-S*F
13          CONTINUE
!     ... to here when finding only eigenvalues.
14        CONTINUE
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.d0
          GO TO 1
        ENDIF
15    CONTINUE
      RETURN
      END

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit real*8(a-h,o-z)
! Indexes an array ARRIN of length N, i.e. outputs the array INDX 
! such that ARRIN(INDX(J)) is in ascending order for J=1,2,...,N. 
! The input quantities N and ARRIN are not changed.
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
         INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
         IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            Q=ARRIN(INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=ARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
20       IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
               IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(Q.LT.ARRIN(INDX(J)))THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
            GO TO 20
         ENDIF
         INDX(I)=INDXT
      GO TO 10
      RETURN
      END

      subroutine setrn(iseed)
        common /rnyucm/ m(4),l(4),nbit,irnyuc
        integer iseed(4)
        ishft12(ii)=ii/4096
        mask12(ii)=mod(ii,4096)
        do 10 i=1,4
10         l(i)=mod(iseed(i),4096)
           l(4)=2*(l(4)/2)+1
           !
           ! shift everything to the left if not 48 bit
           !
           if (nbit.lt.48) then
              do 20 i=1,48-nbit
                 i1=l(1)*2
                 i2=l(2)*2
                 i3=l(3)*2
                 i4=l(4)*2
                 l(4)=mask12(i4)
                 i3=i3+ishft12(i4)
                 l(3)=mask12(i3)
                 i2=i2+ishft12(i3)
                 l(2)=mask12(i2)
                 l(1)=mask12(i1+ishft12(i2))
20               continue
              endif
              return
            end subroutine setrn

      subroutine setrn2(iseed)
        common /rnyucm2/ m(4),l(4),nbit,irnyuc
        integer iseed(4)
        ishft12(ii)=ii/4096
        mask12(ii)=mod(ii,4096)
        do 10 i=1,4
10         l(i)=mod(iseed(i),4096)
           l(4)=2*(l(4)/2)+1
           !
           ! shift everything to the left if not 48 bit
           !
           if (nbit.lt.48) then
              do 20 i=1,48-nbit
                 i1=l(1)*2
                 i2=l(2)*2
                 i3=l(3)*2
                 i4=l(4)*2
                 l(4)=mask12(i4)
                 i3=i3+ishft12(i4)
                 l(3)=mask12(i3)
                 i2=i2+ishft12(i3)
                 l(2)=mask12(i2)
                 l(1)=mask12(i1+ishft12(i2))
20               continue
              endif
              return
            end subroutine setrn2


end module utils
