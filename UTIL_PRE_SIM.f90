MODULE global_variables
    INTEGER:: i,j,i2d,iter,nx,ny,AD_itermax,PPE_itermax,solvetype_ppe,solvetype_AD,iterx,itery
    REAL :: errorx, errory, errorppe, errormax
    REAL :: dx, dy, lx, ly
    REAL :: w_AD, w_PPE, tmax, dt, Re, mu
    REAL :: cflx, cfly, rx, ry, AAD, BAD, APPE, BPPE
    REAL :: x0, y0, r0, r, vt
    REAL :: start, finish
    REAL, ALLOCATABLE, DIMENSION(:) :: x, y, bx, un, ukx, bcx, by, vn, vky, bcy, bp
    REAL, ALLOCATABLE, DIMENSION(:,:) :: u, v, uk, vk, ukp1, vkp1, vor, p, pk, uf, vf,velmag
END MODULE global_variables

MODULE boundary_conditions
    REAL, PARAMETER :: u_bc_w = 0, u_bc_e = 0, u_bc_n = 1, u_bc_s = 0
    REAL, PARAMETER :: v_bc_w = 0, v_bc_e = 0, v_bc_n = 0, v_bc_s = 0 
    REAL, PARAMETER :: p_bc_w = 0, p_bc_e = 0, p_bc_n = 0, p_bc_s = 0
END MODULE boundary_conditions

SUBROUTINE readdata()
    USE global_variables
    ! READing Input Data
    OPEN(2, file='input.dat', status='old')
    READ(2,*) 
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) nx, ny
    READ(2,*)
    READ(2,*)
    READ(2,*) lx, ly
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) w_AD, w_PPE, AD_itermax, PPE_itermax, solvetype
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*)
    READ(2,*) errormax, tmax, dt, Re, mu
    CLOSE(2)

    dx = lx/nx
    dy = lx/ny
END SUBROUTINE

SUBROUTINE domain_init()
    USE global_variables
    REAL :: dxc, dyc
    dxc = lx/(nx-2)
    dyc = ly/(ny-2)

    x(1) = -dxc/2
    y(1) = -dyc/2
    
    DO i = 2,nx
        x(i) = x(i-1) + dxc
    END DO
    
    DO j = 2,ny
        y(j) = y(j-1) + dyc
    END DO


END SUBROUTINE

SUBROUTINE flow_init()

    USE global_variables    

    u = 0
    v = 0
    p = 0

    CALL set_bc()

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
    
END SUBROUTINE

SUBROUTINE set_bc()
    USE global_variables
    USE boundary_conditions

    !for cell center velocities
    u(1,:) = -u(2,:) + u_bc_w*2
    u(nx,:) = -u(nx-1,:) + u_bc_e*2
    u(:,1) = -u(:,2) + u_bc_s*2
    u(:,ny) = -u(:,ny-1) + u_bc_n*2

    v(1,:) = -v(2,:) + v_bc_w*2
    v(nx,:) = -v(nx-1,:) + v_bc_e*2
    v(:,1) = -v(:,2) + v_bc_s*2
    v(:,ny) = -v(:,ny-1) + v_bc_n*2

    p(nx,:) = p(nx-1,:)
    p(1,:) = p(2,:)
    p(:,ny) = p(:,ny-1)
    p(:,1) = p(:,2)

    !for face velocities
    uf(nx-1,:) = u_bc_e
    uf(1,:) = u_bc_w
    vf(:,1) = v_bc_s
    vf(:,ny-1) = v_bc_n

    uf(nx,:) = 2*u(nx,2:ny-1) - uf(nx-1,:)
    vf(:,ny) = 2*v(2:nx-1,ny) - vf(:,ny-1)

END SUBROUTINE
    
    