program rnd

use mtmod

implicit none

real :: Rx,Ry,randX,randY, probX, probY,generate
real :: density,axisX, axisY,a,b,d,part1X,part2X,Fc, counting
density = 10.0
axisX = 30.0
axisY = 30.0
counting = 0.0

open(1, file="rnd61.txt")
open(2, file="rnd62.txt")

do generate = 1,density*axisX*axisY

randX = grnd()
randY = grnd()
Rx = randX*axisX
Ry = randY*axisY
a = 5.0
b = axisX
d = (1/2)*axisX
Fc = (d-a)/(b-a)
part1X = a+sqrt(Rx*(b-a)*(d-a)) !Triangular distribution
part2X = b-sqrt((axisX-Rx)*(b-a)*(b-d))

if (Rx < Fc .and. Rx > 0) then
write(1,*) part1X

elseif (Rx < axisX .and. Rx >= Fc) then
write(1,*) part2X
end if
counting = counting + 1.0
end do

do generate = 1,counting
randY = grnd()
Ry = randY*axisY

if (Ry < axisY .and. Ry > 0) then
write(2,*) Ry
end if
end do 

close(1)

end program rnd
