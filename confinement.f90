module planarTrackConfinement
    use parameters
implicit none
contains
subroutine outputBoundary
    use parameters, only : DP, horizontalLength, verticalLength
    integer :: openStatus333
    open (unit=333, file="chamberBoundary.vtk", status="new", action="write", position="rewind", iostat=openStatus333)
    if(openStatus333>0) stop "*** can't open file ***"
    write(333,'(A)') "# vtk dataFile version 2.0"
    write(333,'(A)') "chamber boundary"
    write(333,'(A)') "ASCII"
    write(333,'(A)') "dataset unstructured grid"
    write(333,'(A, I10, A)') "points", 4, " float"
    write(333,'(3F10.5)') (0.5_DP)*horizontalLength, (-0.5_DP)*verticalLength, 0.0_DP
    write(333,'(3F10.5)') (0.5_DP)*horizontalLength, (0.5_DP)*verticalLength, 0.0_DP
    write(333,'(3F10.5)') (-0.5_DP)*horizontalLength, (0.5_DP)*verticalLength, 0.0_DP
    write(333,'(3F10.5)') (-0.5_DP)*horizontalLength, (-0.5_DP)*verticalLength, 0.0_DP
    write(333,'(A,I10,I10)') "cells", 1, 5
    write(333,'(4I5)') 4, 0, 1, 2, 3
    write(333,'(A, I10)') "cell types ", 1
    write(333,'(I5)') 7
    close(333)
end subroutine outputBoundary
subroutine confinement(tempXCoordinate, tempYCoordinate, tempZcoordinate, xCoordinate, yCoordinate, zCoordinate, confinementStatus)
    use parameters, only : DP, horizontalLength, verticalLength
    real(kind=DP), intent(inout) :: tempXCoordinate, tempYCoordinate, tempZCoordinate
    real(kind=DP), intent(in) :: xCoordinate, yCoordinate, zCoordinate
    logical, intent(inout) :: confinementStatus
    if(tempZCoordinate<0.0_DP) then
        confinementStatus = .false.
        tempZcoordinate = 0.0_DP
    end if
end subroutine confinement
end module planarTrackConfinement
