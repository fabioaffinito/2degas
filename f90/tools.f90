module tools

! double determinant
real(8) :: dtmnt

! must check this
integer, parameter :: lwork=64

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

end module
