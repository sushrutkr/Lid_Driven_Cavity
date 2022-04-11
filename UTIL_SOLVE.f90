subroutine PPESolver()
    use global_variables
    
    
    errorppe = 1
    pk(:,:) = p(:,:)

    ! Calculating face velocities 
    DO j=1,ny-2
        DO i=2,nx-2
            uf(i,j) = (u(i+1,j+1) + u(i,j+1))/2
        END DO
    END DO

    DO j = 2,ny-2
        DO i = 1,nx-2
            vf(i,j) = (v(i+1,j+1) + v(i+1,j))/2
        END DO 
    END DO

    IF (solvetype_ppe .EQ. 2) THEN
        DO WHILE(errorppe>errormax)
            DO j=2,ny-1
                bp = 0
                DO i=1,nx-2
                    i2d = i + 1
                    bp(i) = (1/dt)*((uf(i2d,j-1) - uf(i2d-1,j-1))/dx  + (vf(i2d-1,j) - vf(i2d-1,j-1))/(dy))
                    bp(i) = bp(i) - (p(i2d,j+1) - 2*p(i2d,j) + p(i2d,j-1))/(dy**2)
                END DO
                
                bp(1) = bp(1) - BPPE*p(1,j)
                bp(nx-2) = bp(nx-2) - BPPE*p(nx,j)

                CALL TDMA(nx,bp,APPE,BPPE,BPPE)
                p(2:nx-1,j) = bp(:) 
                p(:,j) = pk(:,j)*(1-w_PPE) + w_PPE*p(:,j)
            END DO

            errorppe = MAXVAL(ABS(p) - ABS(pk))
            pk(:,:) = p(:,:)        
        END DO
    
    ELSE
        !DO iter = 1,PPE_itermax
        DO WHILE(errorppe>errormax)
            DO j=2,ny-1
                DO i=2,nx-1
                    p(i,j) = (1/dt)*(((uf(i,j-1) - uf(i-1,j-1))/dx) + ((vf(i-1,j) - vf(i-1,j-1))/dy))
                    p(i,j) = p(i,j) - ((p(i+1,j) + p(i-1,j))/(dx**2)) - ((p(i,j+1) + p(i,j-1))/(dy**2))
                    p(i,j) = -p(i,j)/((2/(dx**2)) + (2/(dy**2)))
                    p(i,j) = (1-w_PPE)*pk(i,j) + w_PPE*p(i,j)
                END DO
            END DO
            errorppe = MAXVAL(ABS(p) - ABS(pk))
            pk(:,:) = p(:,:)
            !call set_neumann_bc()
        END DO
    END IF
end subroutine

subroutine ADSolver()
    use global_variables
    REAL, DIMENSION(2) :: err
    REAL, DIMENSION(nx,ny) :: sx, sy 

    sx = 0
    sy = 0 

    error = 1
    errorx = 1
    errory = 1
    uk(:,:) = u(:,:)
    vk(:,:) = v(:,:)
    ukp1(:,:) = u(:,:)
    vkp1(:,:) = v(:,:) 

    DO j = 2,ny-1
        DO i = 2,nx-1
            i2d = i
            sx(i,j) = u(i2d,j) - (dt/(2*dx))*(u(i2d+1,j)*uf(i2d+1,j-1) - u(i2d-1,j)*uf(i2d-1,j-1)) &
                                - (dt/(2*dy))*(u(i2d,j+1)*vf(i2d-1,j+1) - u(i2d,j-1)*vf(i2d-1,j-1))

            sy(i,j) = v(i2d,j) - (dt/(2*dx))*(v(i2d+1,j)*uf(i2d+1,j-1) - v(i2d-1,j)*uf(i2d-1,j-1)) &
                                 - (dt/(2*dy))*(v(i2d,j+1)*vf(i2d-1,j+1) - v(i2d,j-1)*vf(i2d-1,j-1))
        END DO
    END DO
    
    
    DO WHILE(error>errormax)
        DO j=2,ny-1
            ! x-Advection Difussion Iterations 
            un = 0
            ukx = 0
            bx = 0 
            bcx = 0

            DO i=1,nx-2
                i2d = i+1
                un(i) = sx(i2d,j)
                ukx(i) = (dt/(Re*(dy**2)))*(ukp1(i2d,j+1) - 2*ukp1(i2d,j) + ukp1(i2d,j-1))
            END DO

            bcx(1) = -BAD*u(1,j) 
            bcx(nx-2) = -BAD*u(nx,j)
            bx = un + ukx + bcx
            
            CALL TDMA(nx, bx, AAD, BAD, BAD)
            
            ukp1(2:nx-1,j) = bx(:)
            ukp1(:,j) = uk(:,j)*(1-w_AD) + w_AD*ukp1(:,j)
            
            ! y-Advection Difussion Iterations
            vn = 0
            vky = 0
            by = 0 
            bcy = 0

            DO i=1,nx-2
                i2d = i+1
                vn(i) = sy(i2d,j)
                vky(i) = (dt/(Re*(dy**2)))*(vkp1(i2d,j+1) - 2*vkp1(i2d,j) + vkp1(i2d,j-1))
            END DO

            bcy(1) = -BAD*v(1,j) 
            bcy(nx-2) = -BAD*v(nx, j)
            by = vn + vky + bcy
            
            CALL TDMA(nx, by, AAD, BAD, BAD)
            
            vkp1(2:nx-1,j) = by(:)
            vkp1(:,j) = vk(:,j)*(1-w_AD) + w_AD*vkp1(:,j) 
        END DO
        errorx = MAXVAL(ABS(ukp1) - ABS(uk))
        uk(:,:) = ukp1(:,:)
        errory = MAXVAL(ABS(vkp1) - ABS(vk))
        vk(:,:) = vkp1(:,:)
        err(1) = errorx
        err(2) = errory
        error = MAXVAL(err) 
        ! call set_ad_bc()
    END DO
    
    u(:,:) = ukp1(:,:)
    v(:,:) = vkp1(:,:)
end subroutine

subroutine vel_correct()
    use global_variables
    !Cell Center Velocity Correction
    DO j=2,ny-1
        DO i=2,nx-1
            u(i,j) = u(i,j) - (dt/(2*dx))*(p(i+1,j) - p(i-1,j))
            v(i,j) = v(i,j) - (dt/(2*dy))*(p(i,j+1) - p(i,j-1))
        END DO 
    END DO

    !Face Center Velocity Correction
    DO j=1,ny-2
        DO i=2,nx-2
            uf(i,j) = uf(i,j) - (dt/dx)*(p(i+1,j+1) - p(i,j+1))
        END DO
    END DO

    DO j = 2,ny-2
        DO i = 1,nx-2
            vf(i,j) = vf(i,j) - (dt/dy)*(p(i+1,j+1) - p(i+1,j))
        END DO 
    END DO
end subroutine
