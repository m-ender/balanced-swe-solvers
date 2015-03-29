! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! Solve Riemann problems for the non-dimensionalised x-split 2D
! shallow water equations with source terms:
!   (h)_t + (u h)_x = 0
!   (uh)_t + ( uuh + .5*h^2 )_x = Khv - hB_x
!   (vh)_t + (uvh)_x = -Khu + KhU
! using Roe's approximate Riemann solver with entropy fix for
! transonic rarefractions.

! This is a modified solver, which solves the source terms in
! a balanced manner, following LeVeque (1998) "Balancing Source
! Terms and Flux Gradients in High-Resolution Godunov Methods:
! The Quasi-Steady Wave-Propagation Algorithm".

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x-momentum
!       3 y-momentum

! This function solves the Riemann problem at all interfaces in one call

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routine step1, rp is called with ql = qr = q.


    implicit double precision (a-h,o-z)

    dimension   ql(meqn,           1-mbc:maxmx+mbc)
    dimension   qr(meqn,           1-mbc:maxmx+mbc)
    dimension    s(mwaves,         1-mbc:maxmx+mbc)
    dimension wave(meqn,   mwaves, 1-mbc:maxmx+mbc)
    dimension amdq(meqn,           1-mbc:maxmx+mbc)
    dimension apdq(meqn,           1-mbc:maxmx+mbc)
    dimension auxl(2,              1-mbc:maxmx+mbc)
    dimension auxr(2,              1-mbc:maxmx+mbc)

!     # Local storage
!     ---------------
    dimension delta(3)
    logical :: efix, newton

    data efix /.true./    !# Use entropy fix for transonic rarefactions
    data newton /.false./ !# Use Newton-Raphson to find cubic roots

    common /cparam/ grav, K, U, dx

!   # Arrays for cell halves
    dimension hl(1-mbc:maxmx+mbc), hr(1-mbc:maxmx+mbc)
    dimension hvl(1-mbc:maxmx+mbc), hvr(1-mbc:maxmx+mbc)

    complex*16 :: CC, u1, u2, u3, x1, x2, x3

    do 10 i=1-mbc,mx+mbc
        h  = ql(1,i)
        hu = ql(2,i)
        hv = ql(3,i)
        DB = auxl(1,i)

        if (newton) then
    !       # delta h for the case hu=0:
            delh = 0.5d0*(-DB + K*hv/h*dx)

    !       # Newton iteration to improve delh:
            do 5 iter=1,5
                hp = ql(1,i) + delh
                hm = ql(1,i) - delh
                F = hu*hu*(1.d0/hp - 1.d0/hm) + 2.0d0 * delh * ql(1,i) + (h*DB - K*hv*dx)
                Fprime = -hu*hu*(1.d0/hp**2 + 1.d0/hm**2) + 2.0d0 * ql(1,i)
                dnewton = F/Fprime
                delh = delh - dnewton
                if (dabs(dnewton).lt.1d-6) go to 8
            5 continue

            write(6,*) 'nonconvergence of newton in rp1swt'
            write(6,*) '   dnewton =',dnewton

            8 continue

        else

    !       # delta h for the case hu=0:
            if (dabs(hu) < 1d-10) then
                delh = 0.5d0*(-DB + K*hv/h*dx)
            else
                a = - 2 * h
                b = (- h * DB + K*hv*dx)
                c = 2*(h**3 - hu**2)
                d = -h*h*b

                !write(6,*) '    h            =', h
                !write(6,*) '    (a, b, c, d) =', a, b, c, d

                Del = 18*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27*a*a*d*d
                Del0 = b*b - 3*a*c
                Del1 = 2*b*b*b - 9*a*b*c + 27*a*a*d
                CC = ((Del1 + sqrt(cmplx(-27*a*a*Del,0,8)))/2)**(1.0d0/3.0d0)
                u1 = 1
                u2 = complex(-1.d0/2.d0, sqrt(3.d0)/2.d0)
                u3 = complex(-1.d0/2.d0, -sqrt(3.d0)/2.d0)

                x1 = -1/(3*a) * (b + u1*CC + Del0/(u1*CC))
                x2 = -1/(3*a) * (b + u2*CC + Del0/(u2*CC))
                x3 = -1/(3*a) * (b + u3*CC + Del0/(u3*CC))

                if (Del > 0) then
                    ! Discriminant greater 0, three real roots
                    ax1 = dabs(real(x1))
                    ax2 = dabs(real(x2))
                    ax3 = dabs(real(x3))
                    if (ax1 <= ax2 .AND. ax1 <= ax3) then
                        delh = real(x1)
                    else if (ax2 <= ax1 .AND. ax2 <= ax3) then
                        delh = real(x2)
                    else
                        delh = real(x3)
                    endif
                else
                    ! Only one real root. Need to find the right one.
                    if (dabs(real(x1)/aimag(x1)) > 10) then
                        delh = real(x1)
                    else if (dabs(real(x2)/aimag(x2)) > 10) then
                        delh = real(x2)
                    else
                        delh = real(x3)
                    endif
                endif

                !write(6,*) '    delh           =', delh
            endif
        endif

        hl(i) = ql(1,i) - delh
        hr(i) = ql(1,i) + delh

        if (U == 0.d0) then
            delhv = 0.5d0*K*(-1/h)*dx*(h*h - delh*delh) + hv*delh/h
        else
            delhv = 0.5d0*K*(U/hu - 1/h)*dx*(h*h - delh*delh) + hv*delh/h
        endif

        hvl(i) = ql(3,i) - delhv
        hvr(i) = ql(3,i) + delhv
    10 continue

!     # Main loop of the Riemann solver.
    do 30 i=2-mbc,mx+mbc

    !    write(6,*) '   i =',i
    !    write(6,*) '   hl =',hr(i-1)
    !    write(6,*) '   hr =',hl(i)

    !     # compute  Roe-averaged quantities:
        hsqrtl = dsqrt(hr(i-1))
        hsqrtr = dsqrt(hl(i))
        hsq2 = hsqrtl + hsqrtr
        ubar = (qr(2,i-1)/hsqrtl + ql(2,i)/hsqrtr) / hsq2
        vbar = (hvr(i-1)/hsqrtl + hvl(i)/hsqrtr) / hsq2
        cbar=dsqrt(0.5d0*(hr(i-1) + hl(i)))

    !     # delta(1)=h(i)-h(i-1), delta(2)=hu(i)-hu(i-1), delta(3)=hv(i)-hv(i-1)
        delta(1) = hl(i) - hr(i-1)
        delta(2) = ql(2,i) - qr(2,i-1)
        delta(3) = hvl(i) - hvr(i-1)

    !     # Compute coeffs in the vector expansion of delta(1),delta(2)
        a1 = 0.5d0*(-delta(2) + (ubar + cbar) * delta(1))/cbar
        a2 = -vbar * delta(1) + delta(3)
        a3 = 0.5d0*( delta(2) - (ubar - cbar) * delta(1))/cbar

    !     # Finally, compute the waves.
        wave(1,1,i) = a1
        wave(2,1,i) = a1*(ubar - cbar)
        wave(3,1,i) = a1*vbar
        s(1,i) = ubar - cbar

        wave(1,2,i) = 0.d0
        wave(2,2,i) = 0.d0
        wave(3,2,i) = a2
        s(2,i) = ubar

        wave(1,3,i) = a3
        wave(2,3,i) = a3*(ubar + cbar)
        wave(3,3,i) = a3*vbar
        s(3,i) = ubar + cbar

    !    write(6,*) '   s1 =',s(1,i)
    !    write(6,*) '   s2 =',s(2,i)
    !    write(6,*) '   s3 =',s(3,i)
    30 END DO

!     # Compute Godunov flux f0:
!     --------------------------

    if (efix) go to 110

!     # No entropy fix
!     ----------------------------------------------
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do 100 m=1, 3
        do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            90 END DO
    100 END DO
    go to 900

!    -----------------------------------------------


    110 continue

!     # With entropy fix
!     ------------------

!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.

    do 200 i=2-mbc,mx+mbc

    ! ------------------------------------------------------
    !        # check 1-wave:
    !        ---------------

    !        # u-c in left state (cell i-1)
        s0 = qr(2,i-1)/hr(i-1) - dsqrt(hr(i-1))

    !        # check for fully supersonic case:
        if (s0 > 0.d0 .AND. s(1,i) > 0.d0)  then
        !            # everything is right-going
            do 60 m=1,3
                amdq(m,i) = 0.d0
            60 END DO
            go to 200
        endif

    !        # u-c to right of 1-wave
        hr1  = hr(i-1) + wave(1,1,i)
        uhr1 = qr(2,i-1) + wave(2,1,i)
        s1 =  uhr1/hr1 - dsqrt(hr1)

        if (s0 < 0.d0 .AND. s1 > 0.d0) then
        !            # transonic rarefaction in the 1-wave
            sfract = s0 * (s1-s(1,i)) / (s1-s0)
        else if (s(1,i) < 0.d0) then
        !	     # 1-wave is leftgoing
            sfract = s(1,i)
        else
        !	     # 1-wave is rightgoing
            sfract = 0.d0   !# this shouldn't happen since s0 < 0
        endif

        do 120 m=1,3
            amdq(m,i) = sfract*wave(m,1,i)
        120 END DO

    ! -------------------------------------------------------
    !        # check 2-wave:
    !        ---------------

        if (s(2,i) > 0.0d0) then
    !       #2 and 3 waves are right-going
            go to 200
        endif

        do 140 m=1,3
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
        140 END DO
    ! -------------------------------------------------------
    !        # check 3-wave:
    !        ---------------
    !        # u+c in right state  (cell i)
        s03 = ql(2,i)/hl(i) + dsqrt(hl(i))

    !        # u+c to left of 3-wave
        hl3  = hl(i) - wave(1,3,i)
        uhl3 = ql(2,i) - wave(2,3,i)
        s3 = uhl3/hl3 + dsqrt(hl3)

        if (s3 < 0.d0 .AND. s03 > 0.d0) then
        !            # transonic rarefaction in the 3-wave
            sfract = s3 * (s03-s(3,i)) / (s03-s3)
        else if (s(3,i) < 0.d0) then
        !            # 3-wave is leftgoing
            sfract = s(3,i)
        else
        !            # 3-wave is rightgoing
            go to 200
        endif

        do 160 m=1,3
            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
        160 END DO
    200 END DO


!     # compute the rightgoing flux differences:
!     # df = SUM s*wave   is the total flux difference and apdq = df - amdq

    do 220 m=1,3
        do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
                df = df + s(mw,i)*wave(m,mw,i)
            210 END DO
            apdq(m,i) = df - amdq(m,i)
    220 END DO

    900 continue
    return
    end subroutine rp1
