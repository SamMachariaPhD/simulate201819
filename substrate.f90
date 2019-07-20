module substrateContraction
    implicit none
    contains
    subroutine substrateDeformation(Ts, xP0Ip, yP0Ip, xPIp, yPIp)
        use parameters, only : range15, DP, timeStepEquil, maxStrain, poissonRatio, strainRate, dt
        integer(kind=range15), intent(in) :: Ts
        real(kind=DP), intent(in) :: xP0Ip, yP0Ip
        real(kind=DP), intent(out) :: xPIp, yPIp
        if(Ts<timeStepEquil) then
            xPIp = xP0Ip + maxStrain*xP0Ip
            yPIp = yP0Ip - poissonRatio*maxStrain*yP0Ip
        else if(maxStrain - strainRate*(Ts - timeStepEquil)*dt>=0.0_DP) then
            xPIp = xP0Ip + (maxStrain - strainRate*(Ts - timeStepEquil)*dt)*xP0Ip
            yPIp = yP0Ip - poissonRatio*(maxStrain - strainRate*(Ts - timeStepEquil)*dt)*yP0Ip
        else
            xPIp = xP0Ip
            yPIp = yP0Ip
        end if
    end subroutine substrateDeformation
    subroutine mapping(Ts, xi, yi, xx, yy)
        use parameters, only : range15, DP, timeStepEquil, maxStrain, poissonRatio, strainRate, dt
        integer(kind=range15), intent(in) :: Ts
        real(kind=DP), intent(in) :: xi, yi
        real(kind=DP), intent(out) :: xx, yy
        if(Ts<timeStepEquil) then
            xx = xi/(1.0_DP + maxStrain)
            yy = yi/(1.0_DP - poissonRatio * maxStrain)
        else if(maxStrain - strainRate*(Ts - timeStepEquil)*dt>=0.0_DP) then
            xx = xi/(1.0_DP + (maxStrain - strainRate*(Ts - timeStepEquil)*dt))
            yy = yi/(1.0_DP - poissonRatio * (maxStrain - strainRate*(Ts - timeStepEquil)*dt))
        else
            xx = xi
            yy = yi
        end if
    end subroutine mapping
    end module substrateContraction
    
    module substrateStretch
    implicit none
    contains
    subroutine substrateDeformation(Ts, xP0Ip, yP0Ip, xPIp, yPIp)
        use parameters, only : range15, DP, timeStepEquil, maxStrain, poissonRatio, strainRate, dt
        integer(kind=range15), intent(in) :: Ts
        real(kind=DP), intent(in) :: xP0Ip, yP0Ip
        real(kind=DP), intent(out) :: xPIp, yPIp
        if(Ts<timeStepEquil) then
            xPIp = xP0Ip
            yPIp = yP0Ip
        else if(maxStrain - strainRate*(Ts - timeStepEquil)*dt >= 0.0_DP) then
            xPIp = xP0Ip + strainRate*(Ts - timeStepEquil)*dt*xP0Ip
            yPIp = yP0Ip - poissonRatio*strainRate*(Ts - timeStepEquil)*dt*yP0Ip
        else
            xPIp = xP0Ip + maxStrain*xP0Ip
            yPIp = yP0Ip - poissonRatio*maxStrain*yP0Ip
        end if
    end subroutine substrateDeformation
    
    subroutine mapping(Ts,xi,yi,xx,yy)
        use parameters, only : range15, DP, timeStepEquil, maxStrain, poissonRatio, strainRate, dt
        integer(kind=range15), intent(in) :: Ts
        real(kind=DP), intent(in) :: xi, yi
        real(kind=DP), intent(out) :: xx, yy
        if(Ts<timeStepEquil) then
            xx = xi
            yy = yi
        else if(maxStrain - strainRate*(Ts - timeStepEquil)*dt >= 0.0_DP) then
            xx = xi/(1.0_DP + strainRate*(Ts - timeStepEquil)*dt)
            yy = yi/(1.0_DP - poissonRatio * strainRate * (Ts - timeStepEquil)*dt)
        else
            xx = xi/(1.0_DP + maxStrain)
            yy = yi/(1.0_DP - poissonRatio * maxStrain)
        end if
    end subroutine mapping
    end module substrateStretch
    
    module substratePeriodicStretch
    implicit none
    contains
    subroutine substrateDeformation(Ts, xP0Ip, yP0Ip, xPIp, yPIp)
    use parameters, only : range15, DP, pi, maxStrain, freq, poissonRatio, dt
    integer(kind=range15), intent(in) :: Ts
    real(kind=DP), intent(in) :: xP0Ip, yP0Ip
    real(kind=DP), intent(out) :: xPIp, yPIp
    xPIp = xP0Ip + maxStrain*xP0Ip*(0.5_DP - 0.5_DP * dcos(2.0_DP*pi*freq*Ts*dt))
    yPIp = yP0Ip - poissonRatio*maxStrain*yP0Ip*(0.5_DP-0.5_DP*dcos(2.0_DP*pi*freq*Ts*dt))
    end subroutine substrateDeformation
    
    subroutine mapping(Ts, xi, yi, xx, yy)
    use parameters, only : range15, DP, pi, maxStrain, freq, poissonRatio, dt
    integer(kind=range15), intent(in) :: Ts
    real(kind=DP), intent(in) :: xi, yi
    real(kind=DP), intent(out) :: xx, yy
    xx = xi/(1.0_DP + maxStrain*(0.5_DP - 0.5_DP*dcos(2.0_DP*pi*freq*Ts*dt)))
    yy = yi/(1.0_DP - poissonRatio*maxStrain*(0.5_DP - 0.5_DP*dcos(2.0_DP*pi*freq*Ts*dt)))
    end subroutine mapping
    end module substratePeriodicStretch