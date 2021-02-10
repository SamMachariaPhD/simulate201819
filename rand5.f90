program rnd

use mtmod

implicit none

real :: Rx,Ry,randX,randY, probX, probY,generate
real :: density,axisX, axisY, counting
density = 10
axisX = 30
axisY = 30
counting = 0.0

open(1, file="rnd51.txt")
open(2, file="rnd52.txt")

do generate = 1,density*axisX*axisY

randX = grnd()
randY = grnd()
Rx = randX*axisX
Ry = randY*axisY

probX = 1.0-(1.0/axisX)*Rx
!print *,probX
!probX = 1-(1/density*axisX*axisY)*generate
!probX = generate/(density*axisX*axisY)
probY = 1.0!-(1/axisX)*Ry
!probY = 1-generate/(density*axisX*axisY)

if (randX < probX) then
write(1,*) Rx
write(2,*) Ry
end if
counting = counting + 1.0
end do

close(1)

end program rnd
