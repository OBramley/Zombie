!*********************************************************************
!** This MODULE contains the random number generator routines.      **
!** It is based on rangen.f written (26/8/2003) by                  **
!** Richard Chandler richard@stats.ucl.ac.uk                        **
!** Paul Northrop northrop@stats.ox.ac.uk                           **
!** and taken from a modified version (25/4/2008) by                **
!** Stewart Reed S.K.Reed@physics.org                               ** 
!** who added the function to convert day month to a day of the     **
!** year. This version has been made to remove the use of common    **
!** blocks and to make it a module.                                 **
!** Oliver Bramley 12/2023                                          **
!*********************************************************************
!** The orginal version can be found at:                            **
!** https://www.ucl.ac.uk/~ucakarc/work/software/randgen.f          **
!** This is true as of 14/12/2023                                   **
!** For ease of reading old comments have been removed.             **
!*********************************************************************



MODULE randgen 

use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, input_unit, output_unit, error_unit

implicit none
real(real64),private::B,C
real(real64),private,dimension(43)::ZBQLIX
integer,private::CURPOS,ID22,ID43

integer,private::INIT
integer,private::STATUS
real(real64)::SPARE
real(real64),parameter::PI = 4.0D0*DATAN(1.0D0)
real(real64),parameter::RLN2P = 0.5D0*DLOG(2.0D0*PI)
real(real64),private,parameter,dimension(0:6)::zbqllg_c=[1.000000000190015D0,76.18009172947146D0,-86.50532032941677D0,&
24.01409824083091D0,-1.231739572450155D0,0.1208650973866179D-2,-0.5395239384953D-5]

contains

    subroutine ZBQLBD01()

        implicit none 
        
        ZBQLIX=[8.001441D7,5.5321801D8,1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,7.9952507D7,3.81202335D8,3.11575334D8,&
        4.02878631D8,2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,&
        2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,3.20029657D8,&
        3.26369301D8,9.441177D6,3.53244738D8,2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,3.75423178D8,7.0893571D7,&
        4.26059785D8,3.95854390D8,2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,2.63576576D8]

        B = 4.294967291D9
        C = 0.0D0 
        CURPOS = 1
        ID22 = 22
        ID43 = 43
        STATUS = -1
    
    end subroutine ZBQLBD01

    function yearday(day,month,year)
    !*****************************************************************
    !  25/4/2008 SKR
    !  a simple function to convert day, month to the day of the year
    !  i.e. a number from 1 to 365 or 366.
    !*****************************************************************
        implicit none
        integer:: day, month, year
        integer:: yearday

        yearday = day
        if (month .ge. 2) then
            yearday = yearday + 31
        end if
        
        if (month .ge. 3) then
            ! account for leap years
            if(dmod(dble(year),4.0d0) .lt. 1.0d-14) then
            if (dmod(dble(year),100d0).lt. 1.0d-14)  then
                yearday= yearday+28
            else
                yearday = yearday + 29
            end if
            else
            yearday = yearday + 28
            end if
        endif
        if (month .ge. 4) then
            yearday = yearday + 31
        end if
        if (month .ge. 5) then
            yearday = yearday + 30
        end if
        if (month .ge. 6) then
            yearday = yearday + 31 ! May
        end if
        if (month .ge. 7) then
            yearday = yearday + 30 ! June
        end if
        if (month .ge. 8) then
            yearday = yearday + 31 !July
        end if
        if (month .ge. 9) then
            yearday = yearday + 31 !August
        end if
        if (month .ge. 10) then
            yearday = yearday + 30 ! September
        end if
        if (month .ge. 11) then
            yearday = yearday + 31 ! October
        end if
        if (month .ge. 12) then
            yearday = yearday + 30 ! November
        endif   
    end function yearday

    SUBROUTINE ZBQLINI(SEED,FU)
        !*****************************************************************
        !       To initialize the random number generator - either
        !       repeatably or nonrepeatably. Need double precision
        !       variables because integer storage can't handle the
        !       numbers involved
        !*****************************************************************
        !	ARGUMENTS
        !	=========
        !	SEED	(integer, input). User-input number which generates
        !		elements of the array ZBQLIX, which is subsequently used 
        !		in the random number generation algorithm. If SEED=0,
        !		the array is seeded using the system clock if the 
        !		FORTRAN implementation allows it.
        !       FU      optional. File unit to write messages to.
        !               Now compulsary as the portland group compiler does
        !               not like it if it is ommited. However if set to 0
        !               then subroutine will not print anything.
        !*****************************************************************
        !	VARIABLES
        !	=========
        !	SEED	See above
        !	ZBQLIX	Seed array for the random number generator. Defined
        !		in ZBQLBD01
        !	B,C	Used in congruential initialisation of ZBQLIX
        !	SS,MM,HH,DD	System clock secs, mins, hours and days
        !	INIT	Indicates whether generator has already been initialised
        !
        integer(int64),intent(inout) :: SEED
        integer :: fu
        integer :: SS,MM,HH,DD,counter
        real(real64) :: TMPVAR1,DSS,DMM,DHH,DDD
        integer,dimension(8) :: values
       
        !	Ensure we don't call this more than once in a program
        IF (INIT.GE.1) THEN
            IF(INIT.EQ.1) THEN
                WRITE(output_unit,"(a)") "****WARNING**** You have called routine ZBQLINI more than once &
                & ignoring any subsequent calls"
                INIT = 2
            END IF
               RETURN
        ELSE
            call ZBQLBD01()
            INIT = 1
        END IF
        
        IF (SEED.EQ.0) THEN
            call date_and_time(VALUES=values)
    
            SS = values(7)   !DATE_AND_TIME used so ifort will be happier
            MM = values(6)   !May need to be changed for other compilers
            HH = values(5)
        
            if (values(1) < 100) values(1) = 2000+values(1)
        
            DD = yearday(values(3), values(2), values(1))
        
            DSS = DINT((DBLE(SS)/6.0D1) * B)
            DMM = DINT((DBLE(MM)/6.0D1) * B)
            DHH = DINT((DBLE(HH)/2.4D1) * B)
            DDD = DINT((DBLE(DD)/3.65D2) * B)
            TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
     
            SEED = IDINT(DSS+DMM+DHH+DDD)
            if(fu .gt. 0) then
                write(fu,*) "Seed from clock ",SEED
            else 
                write(output_unit,*) "Seed from clock ",SEED
            end if 
        else
            if(fu .gt. 0) then
                write(fu,*) "Seed ", SEED
            else 
                write(output_unit,*) "Seed  ",SEED             
            end if
                TMPVAR1 = DMOD(DBLE(SEED),B)
            end if
            ZBQLIX(1) = TMPVAR1
            DO counter = 2,43
                 TMPVAR1 = ZBQLIX(counter-1)*3.0269D4
                 TMPVAR1 = DMOD(TMPVAR1,B)       
                 ZBQLIX(counter) = TMPVAR1
            ENDDO
        
    END SUBROUTINE ZBQLINI

    !*****************************************************************
    ! Returns a uniform random number between 0 & 1, using
    ! a Marsaglia-Zaman type subtract-with-borrow generator.
    ! Uses double precision, rather than integer, arithmetic 
    ! throughout because MZ's integer constants overflow
    ! 32-bit integer storage (which goes from -2^31 to 2^31).
    ! Ideally, we would explicitly truncate all integer 
    ! quantities at each stage to ensure that the double
    ! precision representations do not accumulate approximation
    ! error; however, on some machines the use of DNINT to
    ! accomplish this is *seriously* slow (run-time increased
    ! by a factor of about 3). This double precision version 
    ! has been tested against an integer implementation that
    ! uses long integers (non-standard and, again, slow) -
    ! the output was identical up to the 16th decimal place
    ! after 10^10 calls, so we're probably OK ...
    !*****************************************************************
    FUNCTION ZBQLU01(dummy)
        implicit none
        real(real64):: ZBQLU01,X,B2,BINV
        integer::dummy

        !$omp critical
        B2 = B
        BINV = 1.0D0/B

        X = ZBQLU01_X_value()
        ZBQLIX(ID43) = X
        call ZBQLU01_pointer_check()
       
        ! The integer arithmetic there can yield X=0, which can cause 
        ! problems in subsequent routines (e.g. ZBQLEXP). The problem
        ! is simply that X is discrete whereas U is supposed to 
        ! be continuous - hence if X is 0, go back and generate another
        ! X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
        do while(X.LT.BINV)
            B2 = B2*B
            X = ZBQLU01_X_value()
            ZBQLIX(ID43) = X
            call ZBQLU01_pointer_check()
        end do 
       
        ZBQLU01 = X/B2
        !$omp end critical
        
    END FUNCTION ZBQLU01

    function ZBQLU01_X_value()
        implicit none
        real(real64):: ZBQLU01_X_value,X

        X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
        IF (X.LT.0.0D0) THEN
            X = X + B
            C = 1.0D0
        ELSE
            C = 0.0D0
        END IF

        ZBQLU01_X_value = X
        return

    end function ZBQLU01_X_value
    
    ! Update array pointers. Do explicit check for bounds of each to
    ! avoid expense of modular arithmetic. If one of them is 0 the others
    ! won't be    
    subroutine ZBQLU01_pointer_check()
        implicit none
       
        CURPOS = CURPOS - 1
        ID22 = ID22 - 1
        ID43 = ID43 - 1
        IF (CURPOS.EQ.0) THEN
            CURPOS=43
        ELSE IF (ID22.EQ.0) THEN
            ID22 = 43
        ELSE IF (ID43.EQ.0) THEN
            ID43 = 43
        END IF

    end subroutine ZBQLU01_pointer_check

    !*****************************************************************
    ! Returns a random number uniformly distributed on (X1,X2)
    FUNCTION ZBQLUAB(X1,X2)
        implicit none
        real(real64)::X1,X2,ZBQLUAB
              
        ! Even if X1 > X2, this will work as X2-X1 will then be -ve
        if (abs(X1-X2).lt.1.0d-15) THEN
            ZBQLUAB = X1 + ( (B-X1)*ZBQLU01(1) )
        ELSE
            ZBQLUAB = X1
            WRITE(output_unit,"(a)") "****WARNING**** (function ZBQLUAB) Upper and lower limits on uniform &
            & distribution are identical"
        ENDIF
    
    END function ZBQLUAB
    !*****************************************************************

    ! Returns a random number exponentially distributed with mean MU
    FUNCTION ZBQLEXP(MU)
        implicit none
        real(real64)::MU,ZBQLEXP
    
        ZBQLEXP = 0.0D0
    
        IF (MU.LT.0.0D0) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLEXP"
            RETURN
        ENDIF
    
        ZBQLEXP = -DLOG(ZBQLU01(1))*MU
    
    END function ZBQLEXP
    !*****************************************************************
    ! Returns a random number Normally distributed with mean
    ! MU and standard deviation |SIGMA|, using the Box-Muller
    ! algorithm
    FUNCTION ZBQLNOR(MU,SIGMA)
        implicit none
        real(real64)::THETA,R,ZBQLNOR,MU,SIGMA
        
        !$omp critical 
        IF (STATUS.LE.0) THEN
            THETA = 2.0D0*PI*ZBQLU01(1)
            R = DSQRT( -2.0D0*DLOG(ZBQLU01(1)) )
            ZBQLNOR = (R*DCOS(THETA))
            SPARE = (R*DSIN(THETA))
            STATUS = 1
        ELSE
            ZBQLNOR = SPARE
            STATUS = 0
        ENDIF
        !$omp end critical
        ZBQLNOR = MU + (SIGMA*ZBQLNOR)
    
    END function ZBQLNOR
    !*****************************************************************
    ! Returns a random number binomially distributed (N,P)
    FUNCTION ZBQLBIN(N,P)
        implicit none
        real(real64)::P
        real(real64)::PP,PPP,G,Y,TINY
        INTEGER N,ZBQLBIN,IZ,NN

        TINY = 1.0D-8
        ZBQLBIN = 0

        IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
            WRITE(error_unit,"(a)")"****ERROR**** Illegal parameter value in ZBQLBIN"
            RETURN
        ELSEIF(N.LE.0) THEN
            WRITE(error_unit,"(a)")"****ERROR**** Illegal parameter value in ZBQLBIN"
            RETURN
        ENDIF
    
        ! First step: if NP > 10, say, things will be expensive, and 
        ! we can get into the right ballpark by guessing a value for
        ! ZBQLBIN (IZ, say), and simulating Y from a Beta distribution 
        ! with parameters IZ and NN-IZ+1 (NN starts off equal to N).
        ! If Y is less than PP (which starts off as P) then the IZth order 
        ! statistic from NN U(0,1) variates is less than P, and we know 
        ! that there are at least IZ successes. In this case we focus on 
        ! the remaining (NN-IZ) order statistics and count how many are
        ! less than PP, which is binomial (NN-IZ,(PP-Y)/(1-Y)). 
        ! Otherwise, if Y is greater than PP there must be less 
        ! than IZ successes, so we can count the number of order statistics
        ! under PP, which is binomial (IZ-1,P/Y). When we've got NN*PP
        ! small enough, we go to the next stage of the algorithm and 
        ! generate the final bits directly.
    
        NN = N
        PP = P

        IZ = INT(DBLE(NN)*PP) + 1
        do while((IZ.GT.10).AND.(IZ.LT.NN-10))
            Y = ZBQLBET1(DBLE(IZ),DBLE(NN-IZ+1))
            IF (Y.LT.PP) THEN
                ZBQLBIN = ZBQLBIN + IZ
                NN = NN - IZ
                PP = (PP-Y) / (1.0D0-Y)
            ELSE
                NN = IZ-1
                PP = PP/Y
            ENDIF
            IZ = INT(DBLE(NN)*PP) + 1
        end do
        
        ! PP is the probability of the binomial we're currently
        ! simulating from. For the final part, we simulate either number 
        ! of failures or number of success, depending which is cheaper.
             
        IF (PP.GT.0.5) THEN
            PPP = 1.0D0-PP
        ELSE
            PPP = PP
        ENDIF
        
        G = 0
        IZ = 0
    
        ! ZBQLGEO falls over for miniscule values of PPP, so ignore these
        ! (tiny probability of any successes in this case, anyway)
    
        IF (PPP.GT.TINY) THEN
            G = G + ZBQLGEO(PPP)
            do while(G.LE.NN)
                IZ = IZ + 1
                G = G + ZBQLGEO(PPP)
            end do 
        ENDIF
        
        IF (PP.GT.0.5) IZ = NN - IZ
        ZBQLBIN = ZBQLBIN + IZ
        
       
    END function ZBQLBIN
    !*****************************************************************
    ! Returns a random number geometrically distributed with parameter P ie. mean 1/P
    FUNCTION ZBQLGEO(P)
        implicit none
        real(real64)::P,U,TINY
        INTEGER ZBQLGEO

        TINY = 1.0D-12
        ZBQLGEO = 0
    
        IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLGEO"
            RETURN
        ENDIF

        IF (P.GT.0.9D0) THEN
            ZBQLGEO = ZBQLGEO + 1 
            U = ZBQLU01(1)
            do while(U.GT.P)
                ZBQLGEO = ZBQLGEO + 1
                U = ZBQLU01(1)
            end do
        ELSE
            U = ZBQLU01(1)
            ! For tiny P, 1-p will be stored inaccurately and log(1-p) may
            ! be zero. In this case approximate log(1-p) by -p
            IF (P.GT.TINY) THEN
                ZBQLGEO = 1 + INT( DLOG(U)/DLOG(1.0D0-P) )
            ELSE
                ZBQLGEO = 1 + INT(-DLOG(U)/P)
            ENDIF
        ENDIF

    END function ZBQLGEO
    !*****************************************************************
    ! Returns a random number Poisson distributed with mean MU
    FUNCTION ZBQLPOI(MU)
        implicit none
        real(real64)::X,Y,MU,MU1,TMP1,TMP2,T
        INTEGER ZBQLPOI,K
    
        ZBQLPOI = 0
    
        IF (MU.LT.0.0D0) THEN
            WRITE(error_unit,"(a)") "****ERROR****Illegal parameter value in ZBQLPOI"
            RETURN
        ENDIF

        ! For small MU, generate exponentials till their sum exceeds 1
        ! (equivalently, uniforms till their product falls below e^-MU)

        IF (MU.LE.1.0D3) THEN
            MU1 = MU

            ! For values of MU less than 1000, use order statistics - the Kth
            ! event in a Poisson process of rate MU has a Gamma distribution
            ! with parameters (MU,K); if it's greater than 1 we know that there 
            ! are less than K events in (0,1) (and the exact number is binomial)
            ! and otherwise the remaining number is another Poisson. Choose K so
            ! that we'll get pretty close to 1 in the first go but are unlikely
            ! to overshoot it.
        
            do while (MU1.GT.1.0D1)        
                K = INT(MU1-DSQRT(MU1))
                Y = ZBQLGAM(DBLE(K),MU1)
                IF (Y.GT.1.0D0) THEN
                    ZBQLPOI = ZBQLPOI + ZBQLBIN(K-1,(1.0D0/Y))
                    RETURN
                ENDIF
                ZBQLPOI = ZBQLPOI + K
                MU1 = MU1  * (1.0D0-Y)
            end do 
            Y = DEXP(-MU1)
            X = 1.0D0
            X = X*ZBQLU01(1)

            do while (X.GT.Y)
                ZBQLPOI = ZBQLPOI + 1
                X = X*ZBQLU01(1)
            end do 
            ! For really huge values of MU, use rejection sampling as in 
            ! Press et al (1992) - large numbers mean some accuracy may be
            ! lost, but it caps the execution time.
        ELSE
            TMP1 = DSQRT(2.0D0*MU)
            TMP2 = ZBQLLG(MU+1.0D0)-(MU*DLOG(MU))
            call ZBQLPOI_y_value(Y,MU,TMP1,ZBQLPOI)
            X = DBLE(ZBQLPOI)
            T = (X*DLOG(MU)-ZBQLLG(X+1.0D0)) + TMP2
            do
                IF (DABS(T).LT.1.0D2) THEN
                    T = 0.9D0*(1.0D0+(Y*Y))*DEXP(T)
                    IF (ZBQLU01(1).GT.T) then 
                        call ZBQLPOI_x_y_T_value(X,Y,MU,TMP1,TMP2,ZBQLPOI)
                    else 
                        exit
                    end if
                else 
                    T = DLOG(0.9D0*(1.0D0+(Y*Y))) + T
                    IF (DLOG(ZBQLU01(1)).GT.T) then 
                        call ZBQLPOI_x_y_T_value(X,Y,MU,TMP1,TMP2,ZBQLPOI)
                    else 
                        exit
                    end if
                END IF
            end do
        END IF 
        
    END function ZBQLPOI

    subroutine ZBQLPOI_y_value(Y,MU,TMP1,ZBQLPOI_val)
        implicit none
        real(real64)::Y,MU,TMP1
        INTEGER ZBQLPOI_val
        Y = DTAN(PI*ZBQLU01(1))
        ZBQLPOI_val = INT(MU + (TMP1*Y))
        do while(ZBQLPOI_val.LT.0)
            Y = DTAN(PI*ZBQLU01(1))
            ZBQLPOI_val = INT(MU + (TMP1*Y))
        end do
    end subroutine ZBQLPOI_y_value

    subroutine ZBQLPOI_x_y_T_value(X,Y,MU,TMP1,TMP2,ZBQLPOI_val)
        implicit none
        real(real64)::X,Y,T,MU,TMP1,TMP2
        INTEGER ZBQLPOI_val

        call ZBQLPOI_y_value(Y,MU,TMP1,ZBQLPOI_val)
        X = DBLE(ZBQLPOI_val)
        T = (X*DLOG(MU)-ZBQLLG(X+1.0D0)) + TMP2

    end subroutine ZBQLPOI_x_y_T_value

    !*****************************************************************
    ! Returns a random number with a gamma distribution with mean
    ! G/H and variance G/(H^2). (ie. shape parameter G & scale
    ! parameter H)
    FUNCTION ZBQLGAM(G,H)
        implicit none
        real(real64)::D,R,ZBQLGAM,G,H,A,z1,z2,B1,B2,M
        real(real64)::U1,U2,U,V,TEST,X
        real(real64)::c1,c2,c3,c4,c5,w
        logical::flag=.true.
        ZBQLGAM = 0.0D0
        
        IF ( (G.LE.0.0D0).OR.(H.LT.0.0D0) ) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLGAM, both parameters must be positive"
            RETURN
        ENDIF
        !$omp critical 
        IF (G.LT.1.0D0) THEN
            do 
                u=ZBQLU01(1)
                v=ZBQLU01(1)
                if (u.gt.exp(1.0d0)/(g+exp(1.0d0)))then 
                    ZBQLGAM=-log((g+exp(1.0d0))*(1.0d0-u)/(g*exp(1.0d0)))
                    if (v.gt.ZBQLGAM**(g-1.0))then 
                        cycle 
                    else 
                        exit
                    end if 
                end if
                ZBQLGAM=((g+exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g)
                if (v.gt.exp(-ZBQLGAM)) then
                    cycle
                else
                    exit
                end if
            end do 
            ZBQLGAM=ZBQLGAM/h
            flag=.false.
        ELSE IF (G.LT.2.0D0) THEN
            M = 0.0D0
        ELSE IF (G.gt.10.0d0) then
            c1=g-1.0d0
            c2=(g-1.0d0/(6.0d0*g))/c1
            c3=2.0d0/c1
            c4=c3+2.0d0
            c5=1.0d0/sqrt(g)
            do 
                u=ZBQLU01(1)
                v=ZBQLU01(1)
                if (g.gt.2.50d0) then
                    u=v+c5*(1.0d0-1.860d0*u)
                end if 
                if (u.le.0.0d0.or.u.ge.1.0d0) cycle
                w=c2*v/u 
                if (c3*u+w+1.0d0/w.le.c4) exit
                if (c3*log(u)-log(w)+w.ge.1.0d0) then
                    cycle
                else 
                    exit
                end if
            end do 
            ZBQLGAM=c1*w/h 
            flag=.false.
        ELSE
            M = -(G-2.0D0) 
        ENDIF
        if(flag.eqv..true.)then
            R = 0.50D0
            a = ((g-1.0d0)/exp(1.0d0))**((g-1.0d0)/(r+1.0d0))
            C = (R*(M+G)+1.0D0)/(2.0D0*R)
            D = M*(R+1.0D0)/R
            z1 = C-DSQRT(C*C-D)
    
            ! On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
            ! is not exactly zero - this needs trapping if negative.
    
            IF ((Z1-M.LT.0.0D0).AND.(Z1-M.GT.-1.0D-12)) Z1 = M
            z2 = C+DSQRT(C*C-D)
            B1=(z1*(z1-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z1-M)/(R+1.0D0))
            B2=(z2*(z2-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z2-M)/(R+1.0D0))
            do 
                U1=ZBQLU01(1)
                U2=ZBQLU01(1)
                U=A*U1
                V=B1+(B2-B1)*U2
                X=V/(U**R)
                IF (X.LE.M) cycle
                TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0D0))
                IF (U.LE.TEST) THEN
                    ZBQLGAM = (X-M)/H
                    exit
                ENDIF
            end do 
        end if
        !$omp end critical    
        
    END function ZBQLGAM
    !**************************************************************
    ! Returns a random number, beta distributed with degrees of freedom NU1 and NU2.
    FUNCTION ZBQLBET1(NU1,NU2)
        implicit none
        real(real64)::NU1,NU2,ZBQLBET1,X1,X2
        
        ZBQLBET1 = 0.0D0
        
        IF ( (NU1.LE.0.0).OR.(NU2.LE.0.0) ) THEN
            WRITE(error_unit,"(a)") "****ERROR**** llegal parameter value in ZBQLBET1, both degrees of freedom must be positive"
            RETURN
        ENDIF
            
        ! If parameters are too small, gamma subroutine tends to return zero
        ! as all the probability goes to the origin and we get rounding
        ! errors, even with double precision. In this case, we use Johnk's
        ! method, suitably scaled to avoid rounding errors as much as possible.
    
              
        IF ( (NU1.LT.0.9D0).AND.(NU2.LT.0.9D0) ) THEN
            do
                X1 = ZBQLU01(1)
                X2 = ZBQLU01(1)
                IF ( (X1**(1.0D0/NU1))+(X2**(1.0D0/NU2)).GT.1.0D0) cycle    
                X1 = (DLOG(X2)/NU2) - (DLOG(X1)/NU1)
                ZBQLBET1 = (1.0D0 + DEXP(X1))**(-1)
                IF (ZBQLBET1.GT.1.0D0) then 
                        cycle 
                else 
                        exit
                end if
            end do
        ELSE
            X1 = ZBQLGAM(NU1,1.0D0)
            X2 = ZBQLGAM(NU2,1.0D0)
            ZBQLBET1 = X1/(X1+X2)
        ENDIF
               
        RETURN
        
    
        
    END function ZBQLBET1
    !**************************************************************
    ! Returns a random number, Weibull distributed with shape parameter
    ! A and location parameter B, i.e. density is
    ! f(x) = ( H/(G**H) ) * x**(H-1) * EXP( -(x/G)**H )
    FUNCTION ZBQLWEI(H,G)
        implicit none
        real(real64)::H,G,ZBQLWEI,U
        
        ZBQLWEI = 0.0D0
        
        IF ( (H.LE.0.0).OR.(G.LE.0.0) ) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLWEI, both parameters must be positive"
            RETURN
        ENDIF
    
        U = ZBQLU01(1)
        ZBQLWEI = G * ( (-DLOG(U))**(1.0D0/H) )

    END function ZBQLWEI

    !**************************************************************
    ! Returns a pseudo-random number according to a Negative
    ! Binomial distribution with parameters (R,P). NB these are
    ! both DOUBLE - it copes with non-integer R as well. The
    ! form of the distribution is *not* the no. of trials to 
    ! the Rth success - see documentation for full spec.
    FUNCTION ZBQLNB(R,P)
        implicit none
        real(real64)::R,P,Y
        INTEGER ZBQLNB

        ZBQLNB = 0

        IF ( (R.LE.0.0D0).OR.(P.LE.0.0D0).OR.(P.GE.1.0D0) ) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLNB"
            RETURN
        ENDIF

        Y = ZBQLGAM(R,1.0D0)
        Y = Y*P/(1.0D0-P)
        ZBQLNB = ZBQLPOI(Y)

 
    END function ZBQLNB
    !**************************************************************
    ! Returns a random number, Pareto distributed with parameters
    ! H and G. The density is H*(G**H) / (G+X)**(H+1) for X > 0.
    ! (this is slightly nonstandard - see documentation in 
    ! randgen.txt). The algorithm is straightforward - it uses the
    ! inverse CDF method.
    FUNCTION ZBQLPAR(H,G)
        implicit none
        real(real64)::H,G,ZBQLPAR,U
        
        ZBQLPAR = 0.0D0
        
        IF ( (H.LE.0.0D0).OR.(G.LE.0.0D0) ) THEN
            WRITE(error_unit,"(a)") "****ERROR**** Illegal parameter value in ZBQLPAR, both parameters must be positive"
            RETURN
        ENDIF
         
        U = ZBQLU01(1)
        ZBQLPAR = G * (U**(-1.0D0/H)-1.0D0)
        
       
    END function ZBQLPAR
    !**************************************************************
    ! Returns log(G(X)) where G is the Gamma function. The algorithm is
    ! that given in Press et al (1992), Section 6.1, although this
    ! version also allows for arguments less than 1.
    FUNCTION ZBQLLG(X)
        implicit none
        real(real64)::X,Z,Z2,ZBQLLG,TMP,TOT
        integer::j
    
        ! Compute for x > 1, then use transformation if necessary. Z is
        ! our working argument.
        
        IF (X.GE.1.0D0) THEN
               Z = X
        ELSE 
            Z = 2.0D0-X
            Z2 = 1.0D0-X
        ENDIF
        
        IF (DABS(Z-1.0D0).LT.1.0D-12) THEN
            ZBQLLG = 0.0D0
            RETURN
        ENDIF
        
        TMP = Z + 4.5D0
        TMP = ( (Z-0.5D0)*DLOG(TMP) ) - TMP + RLN2P
        TOT = zbqllg_c(0)
        DO j=1,6
            TOT = TOT + (zbqllg_c(j)/(Z+DBLE(j-1)))
        end do 

        ZBQLLG = TMP + DLOG(TOT)
        
        !Transformation required if X<1
        
        IF (X.LT.1.0D0) THEN
            TMP = PI*Z2
            ZBQLLG = DLOG(TMP/DSIN(TMP)) - ZBQLLG
        ENDIF
        
    END function ZBQLLG
        


END MODULE randgen
