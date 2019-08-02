PROGRAM plainTrajectory

    USE mtmod
    
    IMPLICIT NONE
    
    INTEGER	::			I, J, &
                    Openstatus1, Openstatus2 !, Openstatus3
    
    INTEGER, PARAMETER ::		NumMT = 10, NumCal = 2000, &
                    DP = SELECTED_REAL_KIND(14), &
                    Range15 = SELECTED_INT_KIND(15), &
                    seed = 7358 ! 4397, 242, 6990, 52804, 2072544
    
    INTEGER(KIND = Range15) ::	R_Start, R, R_Next
    
    REAL(KIND = DP) ::		UR, UR1, UR2, NR1, NR2, X, Y, Dis, Angle, dA
    
    REAL(KIND = DP), PARAMETER :: 	pi = 3.14159265358979_DP, &
            dt = 0.1_DP, &
            V_avg = 7.12534_DP, D_v_flu = 0.002_DP, Lp = 15.0_DP
    
    CALL sgrnd(seed)
    
    !OPEN OUTPUT FILE-----------------------------------------------------
    
    OPEN (UNIT = 13, FILE = "mtpaths.txt", STATUS = "NEW", &
            ACTION = "WRITE", POSITION = "REWIND", IOSTAT = Openstatus2)
    IF (Openstatus2 > 0) STOP "*** Can't open file ***"
    
    !INITIAL CONDITION----------------------------------------------------
    
    R = 777137
    
    loopFilament: DO I=1, NumMT
    
    R_Start = R
    X = 0.0_DP
    Y = 0.0_DP
    Angle = 0.0_DP
    
    WRITE(13,'(2I8,3F20.8)') I, 1, X, Y, Angle
    
    !MT TRAJECTORY--------------------------------------------------------
    
    loopTime: DO J=2, NumCal
    
        UR1 = grnd()
        UR2 = grnd()
    
    
        NR1 = DSQRT(-2.0_DP*DLOG(UR1)) * DCOS(2.0_DP*pi*(UR2))
        NR2 = DSQRT(-2.0_DP*DLOG(UR1)) * DSIN(2.0_DP*pi*(UR2) )
    
    
        Dis = V_avg * dt + NR1 * DSQRT(2.0_DP*D_v_flu*dt)
        dA = NR2 * DSQRT(V_avg * dt / Lp)
    
    
        X = X + Dis * DCOS(Angle + dA)
        Y = Y + Dis * DSIN(Angle + dA)
        Angle = Angle + dA
    
    
        WRITE(13,'(2I8,3F20.8)') I, J, X, Y, Angle
    
    END DO loopTime
    
    END DO loopFilament
    
    END PROGRAM plainTrajectory