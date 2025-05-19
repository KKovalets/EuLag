!  EuLag.f90 
!
!  FUNCTIONS:
!  EuLag - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: EuLag
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program EuLag

    implicit none
    
    real(kind=8):: Sp0
    integer:: i, j, k, n, m, nt, ni
    integer:: n_data1, n_data2, count, vartime
    real(kind=8):: k1, k2, k3, k4, m1, m2, m3, m4, x, y, dlim, dmin, dmax, Cm, ratio, dstart, dend
    real (kind=8):: time, dt, null, T1, O1,  T2, O2,  T3, O3, T4, O4, gamma1, gamma2, gamma3, gamma4 
    real (kind=4):: a, b, Tref, Q10, Cr, Cw, K_O2, sigma, eta, epsilon, M0, gamma0,  M0xCm 
    real(kind=8)::  wp2, dwp2,  H, dz, temp, dd0, Z_measur, Sp_measur, Fd_measur
    real (kind=8), allocatable:: O2_conc(:), T(:), d2Tdt(:), d2O2dt(:), d2wpdt(:), d3wpdt(:), d2gammadt(:), d2ddt(:), T_data(:), O2_data(:)
    real (kind=8), allocatable:: gamma_z(:), d(:), Sp(:), Sp_norm(:), Fd_norm(:), wp(:), dwp(:), gamma_interp(:),  wp_interp(:), dwp_interp(:), d_interp(:), Total_Sp(:), Total_Fd(:)
    real (kind = 8), allocatable:: d0(:), z(:),Zgrid(:), Z_data1(:), Z_data2(:)
    common /cm1/a, b,Cr, Tref, Q10, Cw, K_O2, sigma, eta, n_data1, n_data2, vartime
   
    Q10=0. 
    Cw=0.
    Cm=0.
    namelist /Inp1D/ H, n, ni, n_data1, n_data2, vartime, a, b, Cr, Q10, Cw, K_O2, sigma, eta, epsilon, dstart, dend, Tref, Z_measur, Sp_measur, Fd_measur
 
    open(1111,file='Input.nml',form='formatted')
    read(1111,Inp1D)
    
    time = 0.
    
    allocate(Zgrid(n+1), Sp(n+1), gamma_interp(n+1), wp_interp(n+1), dwp_interp(n+1), d_interp(n+1), Total_Sp(n+1), Total_Fd(n+1))
    allocate( O2_conc(n+1), T(n+1), d2Tdt(n_data1), O2_data(n_data2), d2O2dt(n_data2), d2wpdt(2*n), d2ddt(2*n), d3wpdt(2*n), d2gammadt(2*n), T_data(n_data1))   
    allocate(d(2*n), z(2*n), gamma_z(2*n),wp(2*n), dwp(2*n), Sp_norm(n+1), Fd_norm(n+1), Z_data1(n_data1), Z_data2(n_data2))
    dz = (H-Z_measur)/n
    dd0 = (dend-dstart)/ni
    
    
    allocate(d0(ni+1))
    
    d0 = 0.     !particle size spectrum from d_min to d_max
    
    do i = 1, ni+1
        d0(i) = dstart + dd0*(i-1)
    enddo    
    
    dlim = (0.01*dstart)
    null = 0.
    
    Zgrid = 0. !regular Z coordinate, m
    d = 0.     !diameter d_i of some size class i, m
    z = 0.     !irregular z-coordinate z*, m
    wp = 0.    !particle settling velocity, wp, m/s
    Sp = 0.    !Particulate organic matter concentration of some size class i, Sp_i, kg/m^3
    gamma_z = 0. !degradation rate gamma, 1/s
    !Normalized POM concentration Sp and flux Fd (kg/(m^2*s))
    Sp_norm = 0.
    Fd_norm = 0.
    
    O2_conc = 0. !Oxygen concentration, kg/m^3 
    T = 0.       !Temperature, C
    !Measurements data
    O2_data = 0.
    T_data = 0.
    !Derivatives
    d2Tdt = 0. 
    d2O2dt = 0.
    dwp = 0.
    d2wpdt = 0.
    d2ddt = 0. 
    d3wpdt = 0.
    d2gammadt = 0.
    !interpolated over regular Z-grid values
    wp_interp = 0. 
    dwp_interp = 0.
    d_interp = 0.
    gamma_interp = 0.
    !Total POM concentration and flux (calculated using numerical integration)
    Total_Sp = 0. !Total POM concentration Sp, kg/m^3
    Total_Fd = 0. !Total POM flux Fd, kg/(m^2*s)
    open(1, file = "Temperature.dat", status='old', form='formatted', action='read')

    do i = 1, n_data1
        read(1, *) Z_data1(i), T_data(i)
    enddo
    close(1)
    

!    
    open(1, file = "O2 concentration.dat", status='old', form='formatted', action='read')
    do i = 1, n_data2
        read(1, *) Z_data2(i), O2_data(i)
    enddo
    close(1)
   
    call spline(Z_data1, T_data, n_data1, null, null,d2Tdt)
    call spline(Z_data2, O2_data, n_data2, null, null,d2O2dt)
    
    do i = 1, n+1
        Zgrid(i) = Z_measur+(i-1)*dz
        call splint(Z_data1, T_data, d2Tdt, n_data1, Zgrid(i), T(i))
        call splint(Z_data2, O2_data, d2O2dt, n_data2, Zgrid(i), O2_conc(i))
    enddo
     
    Total_Sp = 0.
    Total_Fd = 0.
    
    do j = 1, ni+1, 1
    
    gamma_interp = 0.
        z(1) = Zgrid(1)
        d(1) = d0(j)
        wp(1) = Cw * d0(j)**eta
        Sp(1) = d0(j)**(sigma-epsilon)
        dt = dz/wp(1)
        call gamma_func(time, T(1), O2_conc(1), gamma_z(1))
    i = 0
    do while ((z(i+1).le.H).and.(d(i+1).gt.dlim))
        i = i+1
       
        call splint(Z_data1, T_data, d2Tdt, n_data1,z(i), T1)
        call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i), O1)
        call gamma_func(time, T1, O1, gamma_z(i))
        
        k1 = -gamma_z(i)*d(i)/sigma
        m1 = Cw*d(i)**eta
        
        call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m1*0.5, T2)
        call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m1*0.5, O2)
        call gamma_func(time+dt*0.5, T2, O2, gamma2)
        
        k2 = -gamma2*(d(i)+dt*k1*0.5)/sigma
        m2 = Cw*(d(i)+dt*k1*0.5)**eta
        
        call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m2*0.5, T3)
        call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m2*0.5, O3)
        call gamma_func(time+dt*0.5, T3, O3, gamma3)
        
        k3 = -gamma3*(d(i)+dt*k2*0.5)/sigma
        m3 = Cw*(d(i)+dt*k2*0.5)**eta
        
        
        call splint(Z_data1, T_data, d2Tdt, n_data1, z(i)+dt*m3, T4)
        call splint(Z_data2, O2_data, d2O2dt, n_data2, z(i)+dt*m3, O4)
        call gamma_func(time+dt, T4, O4, gamma4)
        
        k4 = -gamma4*(d(i)+dt*k3)/sigma
        m4 = Cw*(d(i)+dt*k3)**eta
        
        d(i+1) = d(i) + dt*(k1 + 2.*k2 + 2.*k3 + k4)/6.
        z(i+1) = z(i) + dt*(m1 + 2.*m2 + 2.*m3 + m4)/6.
        time = time + dt
        wp(i+1) = Cw * d(i+1)**eta
        x = Cw * d(i+1)**eta
        dwp(i+1) = (wp(i+1)-wp(i))/(z(i+1)-z(i))
        dt = dz/wp(i+1)

    enddo
    m = i
    dwp(1) =(wp(2)-wp(1))/dz
    time = 0 

     if (m.gt.n+1) then
        count = n+1
    else
        count = m
    endif
    
    call spline(z, wp, m, null, null,d2wpdt)
    call spline(z, dwp, m,null, null,d3wpdt)
    call spline(z, gamma_z, m, null, null,d2gammadt)
    call spline(z, d, m, null, null,d2ddt)
    
    do i = 1, count, 1
    
        call splint(z, wp, d2wpdt, m, Zgrid(i), wp_interp(i))
        call splint(z, dwp, d3wpdt, m, Zgrid(i), dwp_interp(i))
        call splint(z, d, d2ddt, m, Zgrid(i), d_interp(i))
        call splint(z, gamma_z, d2gammadt, m, Zgrid(i), gamma_interp(i))

        if (Zgrid(i).gt.z(m)) then
            exit
        endif
    enddo
    if (i.ne.(count+1)) then
        count = i 
    endif
      
    if (j.eq.1.or.j.eq.ni+1) then
        Total_Sp(1) = Total_Sp(1) + Sp(1)
        Total_Fd(1) = Total_Fd(1) + Cw*d0(j)**(eta+sigma-epsilon)
    else 
        if (mod(j,2).eq.0) then
            Total_Sp(1) = Total_Sp(1) + 4*Sp(1)
            Total_Fd(1) = Total_Fd(1) + 4*Cw*d0(j)**(eta+sigma-epsilon)
        else
            Total_Sp(1) = Total_Sp(1) + 2*Sp(1)
            Total_Fd(1) = Total_Fd(1) + 2*Cw*d0(j)**(eta+sigma-epsilon)
        endif
    endif
    
        do i = 1, count-1, 1
        
            k1 = -(dwp_interp(i)+gamma_interp(i))*Sp(i)/wp_interp(i)
            call splint(z, wp, d2wpdt, m, Zgrid(i) + dz/2, wp2)
            call splint(z, dwp, d3wpdt, m, Zgrid(i) + dz/2, dwp2)
            call splint(z, gamma_z, d2gammadt, m, Zgrid(i) + dz/2, gamma1)
           
            
            k2 = -(dwp2+gamma1)*(Sp(i) + dz*k1/2)/wp2
            k3 = -(dwp2+gamma1)*(Sp(i) + dz*k2/2)/wp2
            k4 = -(dwp_interp(i+1)+gamma_interp(i+1))*(Sp(i) + dz*k3)/wp_interp(i+1)
           
            Sp(i+1)= Sp(i)+dz*(k1+2.*k2+2.*k3+k4)/6.
        
            if (j.eq.1.or.j.eq.ni+1) then
                    Total_Sp(i+1) = Total_Sp(i+1) + Sp(i+1)
                    Total_Fd(i+1) = Total_Fd(i+1) + Sp(i+1)* wp(i+1)
            else 
                if (mod(j,2).eq.0) then
                    Total_Sp(i+1) = Total_Sp(i+1) + 4*Sp(i+1)
                    Total_Fd(i+1) = Total_Fd(i+1) + 4*Sp(i+1)* wp(i+1)
                else
                    Total_Sp(i+1) = Total_Sp(i+1) + 2*Sp(i+1)
                    Total_Fd(i+1) = Total_Fd(i+1) + 2*Sp(i+1)* wp(i+1)
                endif 
            endif
        enddo
     enddo
     
    do i = 1, count, 1  
        Total_Sp(i) = Total_Sp(i)*dd0/3
        Total_Fd(i) = Total_Fd(i)*dd0/3
    enddo
    if (Sp_measur.ne.0.) then
        M0xCm = Sp_measur/Total_Sp(1)
    else 
        M0xCm = Fd_measur/Total_Fd(1)
    endif
        
    do i = 1, count, 1  
        Total_Sp(i) = Total_Sp(i)*M0xCm
        Total_Fd(i) = Total_Fd(i)*M0xCm
    enddo
    
    do i = 1, count, 1
        Sp_norm(i) = Total_Sp(i)/Total_Sp(1)
        Fd_norm(i) = Total_Fd(i)/Total_Fd(1)
    enddo   

    OPEN(17,FILE="M0xCm.dat")
        write(17,*) M0xCm
    CLOSE(17)
    
    do i = 1, count, 1
        Zgrid(i) = Zgrid(i)-Z_measur
    enddo   
    if (vartime.eq.1) then
        call write_arr(Zgrid(:),Total_Sp(:),count,'Sp(z)_a=',8,int(a/3600/24))
        call write_arr(Zgrid(:),Total_Fd(:),count,'Fd(z)_a=',8,int(a/3600/24))
        call write_arr(Zgrid(:),gamma_interp(:),count,'gamma_a=',8,int(a/3600/24))
        call write_arr(Zgrid(:),Sp_norm(:),count,'Sp(z)_norm_a=',13,int(a/3600/24))
        call write_arr(Zgrid(:),Fd_norm(:),count,'Fd(z)_norm_a=',13,int(a/3600/24))
    else
        
        call write_arr(Zgrid(:),Total_Sp(:),count,'Sp(z)_',6,0)
        call write_arr(Zgrid(:),Total_Fd(:),count,'Fd(z)_',6,0)
        call write_arr(Zgrid(:),gamma_interp(:),count,'gamma_',6,0)
        call write_arr(Zgrid(:),Sp_norm(:),count,'Sp(z)_norm_',11,0)
        call write_arr(Zgrid(:),Fd_norm(:),count,'Fd(z)_norm_',11,0)
    endif

    deallocate(Zgrid, Sp, gamma_interp, wp_interp, dwp_interp, d_interp, Total_Sp, Total_Fd)
    deallocate( O2_conc)
    deallocate(T)
    deallocate(d2Tdt)
    deallocate(O2_data) 
    deallocate(d2O2dt)
    deallocate(d2wpdt)
    deallocate(d2ddt)
    deallocate(d3wpdt)
    deallocate(d2gammadt)
    deallocate(T_data)   
    deallocate(d)
    deallocate(z, gamma_z,wp, dwp, Sp_norm, Fd_norm, Z_data1, Z_data2)    
    deallocate(d0)
    


    end program EuLag


    subroutine gamma_func(time, T, O2_conc, gamma_z)
        real (kind = 8):: time, gamma_z, T, O2_conc
        real (kind = 4):: a, b, Cr, Tref, Q10, Cw, K_O2, sigma, eta
        integer n_data1,n_data2, vartime
        common /cm1/a, b, Cr, Tref, Q10, Cw, K_O2, sigma, eta, n_data1, n_data2, vartime
        
        if (vartime.eq.0) then
            gamma_z = Cr
        else
            gamma_z = b/(a + time)
        endif
        
        gamma_z = gamma_z * (Q10**((T-Tref)/10.))
        if ((O2_conc.eq.0.).or.(K_O2.eq.0.)) then
            gamma_z = gamma_z
        else
            gamma_z = gamma_z * O2_conc/(K_O2 + O2_conc)
        endif

    end subroutine
    
    subroutine spline(x,y,n,yp1,ypn,y2)
    
        integer n,NMAX
        real(kind = 8)::x(n), yp1,ypn,y(n),y2(n)
        integer i,k
        real(kind = 8):: p,qn,sig,un,u(n)
       
        if (yp1.gt..99e30) then 
            y2(1)=0.
            u(1)=0.
        else 
            y2(1)=-0.5
            u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        
        do i=2,n-1 
            sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
            p=sig*y2(i-1)+2.
            y2(i)=(sig-1.)/p
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo 
    
        
        if (ypn.gt..99e30) then
            qn=0.
            un=0.
        else 
            qn=0.5
            un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
        do k=n-1,1,-1 
            y2(k)=y2(k)*y2(k+1)+u(k)
        enddo 
    end

   SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER n
    REAL(kind = 8):: x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    REAL(kind = 8)::  a,b,h
    klo=1 
    khi=n
    1 if (khi-klo.gt.1)then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
    khi=k
    else
    klo=k
    endif
    goto 1
    endif 
    h=xa(khi)-xa(klo)
    if (h.eq.0.) pause 
    a=(xa(khi)-x)/h 
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

    return
    END
   
       
    subroutine write_arr(arr1,arr2,N,fname1,len,tlev)
    implicit none
    ! writes arr1 and arr2 as columns in file
    integer N,tlev,len
    real(kind = 8) arr1(N)
    real(kind = 8) arr2(N)
    character(len) fname1
    integer len1,curlen,i
    character*1000 fname
    character*10 str1,tstr
    
    str1= ''
    tstr = ''
    fname = ''
    len1=len_trim(fname1)
    fname(1:len1)=fname1(1:len1)
    curlen=len1
    write(str1,'(I10)')tlev
    tstr=ADJUSTL(str1)
    len1=len_trim(tstr)
    fname(curlen+1:curlen+len1)=tstr(1:len1)
    curlen=curlen+len1
    fname(curlen+1:curlen+4)='.dat'

    OPEN(2024,FILE=fname)
    do i=1,N
    write(2024,*)arr1(i),arr2(i)
    end do
    CLOSE(2024)
    end subroutine

