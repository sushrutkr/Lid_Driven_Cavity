subroutine writepostproc()
    
    use global_variables
    use boundary_conditions

    !modifying domain data
    x(1) = 0 
    x(nx) = lx
    y(1) = 0
    y(ny) = ly

    !modifying the velocity data
    j = 1
    DO i = 1,nx
        u(i,j) = u_bc_s
        v(i,j) = v_bc_s
    END DO

    !top wall data
    j = ny
    DO i=1,nx
        u(i,j) = u_bc_n
        v(i,j) = v_bc_n
    END DO
    
    !left wall data
    i=1
    DO j=2,ny-1
        u(i,j) = u_bc_w
        v(i,j) = v_bc_w
    END DO

    !right wall data
    i=nx
    DO j=2,ny-1
        u(i,j) = u_bc_e
        v(i,j) = v_bc_e
    END DO


    vor = 0
    DO j=2,ny-1
        DO i = 2,nx-1
            vor(i,j) = ((v(i+1,j) - v(i-1,j))/(2*dx)) - ((u(i,j+1) - u(i,j-1))/(2*dy)) 
        END DO
    END DO

    DO j=1,ny
        DO i = 1,nx
            velmag(i,j) = (u(i,j)**2 + v(i,j)**2)**0.5 
        END DO
    END DO
    
    
    !Writing Files For Post Processing
    open(12, file='data.dat', status='unknown')
    WRITE(12,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(12,*) 'VARIABLES = "X", "Y", "u", "v","Velocity Magnitude", "P", "Vorticity"'
    WRITE(12,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    DO j=1,ny
        DO i = 1,nx
            WRITE(12,*) x(i), y(j), U(i,j), V(i,j),velmag(i,j), P(i,j), vor(i,j) 
        END DO
    END DO
    close(12)
end subroutine