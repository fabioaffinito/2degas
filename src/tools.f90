module tools

! double determinant
real(8) :: dtmnt

! complex determinant
complex*16 :: zdtmnt

! must check this
integer, parameter :: lwork=64

!$omp threadprivate(dtmnt,zdtmnt)

CONTAINS
!
! gets inverse and determinant using lapack
!
subroutine dget_inverse(matrix,lda,n,dtmnt,info)
implicit none

! arguments
integer, intent(in) :: lda,n
integer, intent(out) :: info
real(8), intent (INOUT) :: matrix(lda,n)
real(8), intent (out) :: dtmnt


! local vars
integer :: ipiv(n)   !  check this
integer :: i

real(8) :: work(lwork)

#ifdef DEBUG
   write(*,*) 'called dget_inverse'
#endif
! First perform LU decomposition with Gaussian elimination

call dgetrf(lda,n,matrix,lda,ipiv,info)
if (info .ne. 0) then
   return
endif

! Get determinant
dtmnt = 1.0d0
do i=1,n
   if (ipiv(i).ne.i) then
      dtmnt = -dtmnt * matrix(i,i)
   else
      dtmnt = dtmnt * matrix(i,i)
   endif
enddo

! Get inverse matrix via lapack
!
call dgetri(n, matrix, lda,ipiv, work, lwork, info)

if (info .ne. 0) then
   return
endif     

end subroutine dget_inverse

!
! complex version of dget_inverse
!
subroutine zget_inverse(matrix,lda,n,zdtmnt,info)
implicit none

! arguments
integer, intent(in) :: lda,n
integer, intent(out) :: info
complex*16, intent (INOUT) :: matrix(lda,*)
complex*16, intent (OUT) :: zdtmnt

complex*16 :: sgn
! local vars
integer :: ipiv(n)   !  check this
integer :: i

complex*16 :: work(lwork)
COMPLEX*16         ONE
PARAMETER   ( ONE = ( 1.0D+0, 0.0D+0 ) )

! First perform LU decomposition with Gaussian elimination

call zgetrf(lda,n,matrix,lda,ipiv,info)
if (info .ne. 0) then
   return
endif

! Get determinant
  zdtmnt = ONE

 do i=1,n

    zdtmnt=zdtmnt*matrix(i,i)

 enddo


! do i=1,n
!    if (ipiv(i).ne.i) then
!       zdtmnt = -zdtmnt * matrix(i,i)
!    else
!       zdtmnt = zdtmnt * matrix(i,i)
!    endif
! enddo

!  This code should be quicker

  sgn=ONE

  do i=1, N
     if(ipiv(i)/=i)then
        sgn=-sgn
     end if
  end do
  zdtmnt=sgn*zdtmnt

#ifdef DEBUG
    write(*,*) 'in zget_inverse: n=',n,' zdtmnt=',zdtmnt
#endif
!

! Get inverse matrix via lapack
!
call zgetri(n, matrix, lda,ipiv, work, lwork, info)

! if info .ne.0 there is a problem
return

end subroutine zget_inverse

complex(kind=8) function get_det(N, mat)

    implicit none

    integer(kind=8), intent(in) :: N 

    complex(kind=8), intent(inout), dimension(:,:) :: mat

    integer(kind=8) :: i, info

    integer, allocatable :: ipiv(:)

    real(kind=8) :: sgn


    allocate(ipiv(N))


    ipiv = 0


    call zgetrf(N, N, mat, N, ipiv, info)

    get_det = 1

    do i = 1, N

        get_det = get_det*mat(i, i)

    end do

    sgn = 1

    do i = 1, N

        if(ipiv(i) /= i) then

            sgn = -sgn

        end if

    end do

    get_det = sgn*get_det   

end function get_det
end module
