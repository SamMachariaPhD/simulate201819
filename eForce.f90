module uniformForceX
    use parameters
implicit none
contains
function extFx(Ts,x,y,z)
    use parameters, only : DP, range15, numMol, numBeads, bondLength, extForceDensity0, timeForceOn
    real(kind=DP), dimension(numBeads, numMol) :: extFx
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: x, y, z
    integer(kind=range15), intent(in) :: Ts
    integer :: i,j
    extFx = 0.0_DP
    if(Ts>=timeForceOn) then
        do i=1, numMol
            extFx(1,i) = 0.5_DP*extForceDensity0*bondLength
            do j=2, numBeads-1
                extFx(j,i) = extForceDensity0*bondLength
            end do
            extFx(numBeads,i) = 0.5_DP*extForceDensity0*bondLength
        end do
    end if
end function extFx
function extFy(Ts,x,y,z)
    use parameters, only : DP, range15, numMol, numBeads, bondLength, extForceDensity0, timeForceOn
    real(kind=DP), dimension(numBeads, numMol) :: extFy
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: x, y, z
    integer(kind = range15), intent(in) :: Ts
    integer :: i,j
    extFy = 0.0_DP
end function extFy
function extFz(Ts,x,y,z)
    use parameters, only : DP, range15, numMol, numBeads, bondLength, extForceDensity0, timeForceOn
    real(kind=DP), dimension(numBeads, numMol) :: extFz
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: x,y,z
    integer(kind=range15), intent(in) :: Ts
    integer :: i,j
    extFz = 0.0_DP
end function extFz
end module uniformForceX