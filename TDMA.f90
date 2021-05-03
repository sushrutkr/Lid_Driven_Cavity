subroutine TDMA(n, b, Coeffd, Coeffl, Coeffu)
    implicit none
    
    integer :: i, n
    real :: Coeffd, Coeffl, Coeffu
    real, dimension(n-2) :: b, d
    real,dimension(n-3) :: l,u
    
    d = Coeffd
    u = Coeffu
    l = Coeffl
    
    b(1) = b(1)/d(1)
    u(1) = u(1)/d(1)
    d(1) = 1
    
    do i=2,n-3
        b(i) = b(i) - b(i-1)*l(i-1)
        d(i) = d(i) - u(i-1)*l(i-1)
        l(i-1) = 0
        
        b(i) = b(i)/d(i)
        u(i) = u(i)/d(i)
        d(i) = 1
    end do 
    
    b(n-2) = b(n-2) - b(n-3)*l(n-3)
    d(n-2) = d(n-2) - u(n-3)*l(n-3)
    l(n-3) = 0
    
    b(n-2) = b(n-2)/d(n-2)
    
    !backsubsitution
    do i=n-3,1,-1
        b(i) = b(i) - b(i+1)*u(i)
    end do     
    
end subroutine
    