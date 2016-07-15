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


end module utils
