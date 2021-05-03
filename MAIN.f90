! My CFD code 
! Compile - gfortran UTIL_PRE_SIM.f90 UTIL_POST_PROC.f90 TDMA.f90 UTIL_SOLVE.f90 MAIN.f90

program CFDCode
    use global_variables
    CALL cpu_time(start)
    CALL readdata()

    open(11, file='log.dat', status='old')
    write(11,*) " CODE BEGUN.........................."
    write(11,*)
    write(11,*) " Simulation Parameters are nx, ny, dx, dy, w, errormax,tmax, dt, Re "
    write(11,*) nx, ny, dx, dy, w_PPE, errormax,tmax, dt, Re
    write(11,*)
    write(11,*) "Initialising the domain .........................."
    write(11,*)  
    allocate(x(nx), y(ny))
    allocate(u(nx,ny), uk(nx,ny),ukp1(nx,ny)) 
    allocate(v(nx,ny), vk(nx,ny), vkp1(nx,ny)) 
    allocate(p(nx,ny), pk(nx,ny), uf(nx,ny-2), vf(nx-2,ny),bp(nx-2))
    allocate(vor(nx,ny), velmag(nx,ny))
    allocate(bx(nx-2), by(nx-2), un(nx-2), vn(nx-2), ukx(nx-2), vky(nx-2), bcx(nx-2), bcy(nx-2))
    
    ! Defining Coefficient for [A]x=b   
    AAD = 1 + ((2*dt)/(Re*(dx**2)))
    BAD = -dt/(Re*(dx**2))

    APPE = -2/(dx**2)
    BPPE = 1/(dx**2)

    CALL domain_init()    
    CALL flow_init()
    write(11,*) "Initiating Calculations .........................."
    write(11,*)
    write(11,*) "time,  Error_x-AD,  Error_y-AD,  Error_PPE"
    
    t = 0    
    DO WHILE (t<tmax)
        CALL set_bc()

        CALL ADSolver()
        
        CALL PPESolver()

        CALL vel_correct() 

        write(11,*) t, errorx, errory, errorppe
        t = t+dt
    END DO

    write(11,*)
    write(11,*) "Calculations Done  .............................."
    write(11,*)
    write(11,*) "Writing file for Post-Processing  ................"
    write(11,*)
      
    CALL writepostproc()
    
    deallocate(x,y)
    deallocate(u,uk,ukp1)
    deallocate(v,vk,vkp1)
    deallocate(vor,velmag)
    deallocate(p, pk, bp, vf, uf)
    deallocate(bx, by, un, vn, ukx, vky, bcx, bcy)
    call cpu_time(finish)

    write(11,*) "Total Simulation Time : ", (finish-start)/60 , "mins"
    close(11)

end program CFDCode
