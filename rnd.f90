program rnd

use random

use mtmod

implicit none

real :: m,nx
!real :: random_exponential

open(1, file="rnd.txt")

!write(1,*)"Here are some data:"

do m = 1,10000
!call random_number(n)
!nx = random_exponential()
!nx = random_chisq(3,.true.)
!nx = random_Weibull(3.0)
!nx = grnd()
nx = random_beta(1.0,2.0,.true.)
write(1,*) nx
end do

close(1)

end program rnd
