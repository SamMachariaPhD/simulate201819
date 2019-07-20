program beadRodPolymer
    use parameters
    use mtmod
    use planarTrackConfinement
    use uniformForceX
implicit none

!==========================================INITIALIZATION======================================!

integer :: seedAssay, iAssay, i, j, ii, jj, jjj, jContact, areaCounter, iArea, eraseCounter, iErase, &
numIteration, outputFileCounter, iInt, counterBuffer, iSpecies
integer(kind=range15) :: Ts, Ip, numParticles
integer(kind=range15), dimension(maxAreaNum) :: addedParticleNum
integer, dimension(maxNumParticles) :: motorSpecies
integer, dimension(numMol*numBeads, maxNumParticles) :: candidateList
integer, dimension(maxNumParticles) :: releaseADP
logical :: stateConstraint, stateConfinement, stateConfinementParticle
real(kind=DP) :: initialAngle, ur1, ur2, ur3, normRandParticles, interceptContact, xContact, yContact, zContact, &
disHeadTail, projHeadTail, gammaBead, dBead, gammaParticle, vMotor, intercept, xcm, ycm
real(kind=DP), dimension(maxAreaNum) :: areaOriginX, areaOriginY, areaOriginUx, areaOriginUy
real(kind=DP), dimension(numBeads) :: xi, yi, zi, xiTemp, yiTemp, ziTemp, fiX, fiY, fiZ, &
normRandVectorBeads, xiNormRandVectorBeads, yiNormRandVectorBeads, ziNormRandVectorBeads
real(kind=DP), dimension(maxNumParticles) :: xp, yp, zp, normRandVectorParticles, xpBuffer, &
ypBuffer, zpBuffer, elongation
real(kind=DP), dimension(numMol, maxNumParticles) :: contactState, tempContact
real(kind=DP), dimension(maxNumParticles) :: fMotorX, fMotorY, fMotorZ
real(kind=DP), dimension(numBeads, numMol) :: x, y, z, fX, fY, fZ, forceBeadX, forceBeadY, forceBeadZ, &
xNormRandVectorBeads, yNormRandVectorBeads, zNormRandVectorBeads
character(len=30) :: outFileName, outFileNameF, outFileNameP
seedAssay = seed

!===============================INITIAL CALCULATIONS AND OUTPUT================================!

gammaBead = 3.0_DP*pi*0.001_DP*bondLength/DLOG(bondLength/0.006_DP)
dBead = kBT/gammaBead
call outputBoundary
open(10,file='initialCondition.txt')
loopAssay: do iAssay=1, numAssay
call sgrnd(seedAssay)
Ts = 0
numParticles = 0
areaCounter = 1
eraseCounter = 0
outputFileCounter = 0
xp = -100.0_DP
yp = 5.0_DP
initialAngle = 0.0_DP*pi
write(10,'(I3,I7,F15.10)') iAssay, seedAssay, initialAngle
seedAssay = seedAssay+1
write(outFileName,'(A,I3.3,A)') 'conformationA', iAssay, '.txt'
open(11,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'plusEndPositionsA', iAssay, '.txt'
open(15,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'tipXYA', iAssay, '.txt'
open(16,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'trajectoriesA', iAssay, '.txt'
open(17,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'headTailA', iAssay, '.txt'
open(20,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'allA', iAssay, '.txt'
open(21,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'forceA', iAssay, '.txt'
open(22,file=outFileName)
write(outFileName,'(A,I3.3,A)') 'areaEraseCounterA', iAssay, '.txt'
open(23,file=outFileName)

!=====================================INITIAL ORIENTATION===================================!

do i=1, numMol
    x(numBeads,i)=dcos(initialAngle)
    y(numBeads,i)=dsin(initialAngle)
    z=0.0125_DP
    do j=numBeads-1, 1, -1
        x(j,i)=x(j+1,i)+bondLength*dcos(initialAngle)
        y(j,i)=y(j+1,i)+bondLength*dsin(initialAngle)
    end do
end do

xcm=sum(x(:,1))/dble(numBeads)
ycm=sum(y(:,1))/dble(numBeads)
areaOriginUx(areaCounter)=dcos(initialAngle)
areaOriginUy(areaCounter)=dsin(initialAngle)
areaOriginX(areaCounter)=xcm
areaOriginY(areaCounter)=ycm

do Ip=1, int(motorDensity*horizontalLength*verticalLength, range15)
    ur1=grnd()
    ur2=grnd()
    xp(Ip)=xcm+horizontalLength*(ur1-0.5_DP)*areaOriginUx(areaCounter)-verticalLength*(ur2-0.5_DP)*areaOriginUy(areaCounter)
    yp(Ip)=ycm+horizontalLength*(ur1-0.5_DP)*areaOriginUy(areaCounter)+verticalLength*(ur2-0.5_DP)*areaOriginUx(areaCounter)
    write(21,'(3F15.10)') xp(Ip), yp(Ip)
    ur1=grnd()
    if(ur1<=speciesRatio)then
        motorSpecies(Ip)=1
        ur2=grnd()
        if(ur2<=0.091)then
            contactState(:,Ip)=-1.0_DP
        else
            contactState(:,Ip)=0.0_DP
        end if
    else
        motorSpecies(Ip)=2
        contactState(:,Ip)=0.0_DP
    end if
end do

numParticles = numParticles+int(motorDensity*horizontalLength*verticalLength,range15)
addedParticleNum(areaCounter) = int(motorDensity*horizontalLength*verticalLength,range15)
releaseADP=0
candidateList(:,1:numParticles)=1
candidateList(:,(numParticles+1):maxNumParticles)=0
tempContact=0.0_DP
call checkMotorFilamentContact(x,y,z,xp,yp,zp,numParticles,candidateList,contactState,tempContact)
do Ip=1, numParticles
    do ii=1, numMol
        if((tempContact(ii,Ip)>=1.0_DP) .and. (tempContact(ii,Ip)<=dble(numBeads))) then
            ur1=grnd()
            if(ur1<=0.1_DP)then
                contactState(ii,Ip)=tempContact(ii,Ip)
                tempContact(ii,Ip)=0.0_DP
            end if
        end if
    end do
end do
call calculateForceMotor(x,y,z,xp,yp,zp,contactState,fMotorX,fMotorY,fMotorZ,numParticles,elongation)
do Ip=1,numParticles
    call motorForcedDetachment(fMotorX(Ip),fMotorY(Ip),fMotorZ(Ip),contactState(:,Ip),elongation(Ip))
end do
call calculateForceBead(contactState,forceBeadX,forceBeadY,forceBeadZ,fMotorX,fMotorY,fMotorZ,numParticles)
fX=forceBending(x)+forceBeadX
fY=forceBending(y)+forceBeadY
fZ=forceBending(z)+forceBeadZ
fX=forceBending(x)+forceBeadX+extFx(Ts,x,y,z)
fY=forceBending(y)+forceBeadY+extFy(Ts,x,y,z)
fZ=forceBending(z)+forceBeadZ+extFz(Ts,x,y,z)
do j=1, numBeads
    write(11,'(I10)',advance="no")Ts
    do i=1, numMol-1
        write(11,'(3F25.15)',advance="no") x(j,i),y(j,i),z(j,i)
    end do
    write(11,'(3F25.15)') x(j,numMol), y(j,numMol), z(j,numMol)
end do
do i=1, numMol-1
    write(15,'(3F14.6)',advance="no") x(numBeads,i), y(numBeads,i), z(numBeads,i)
end do
write(15,'(3F14.6)') x(numBeads,numMol), y(numBeads,numMol), z(numBeads,numMol)
do i=1, numMol-1
    write(16,'(I10,2F14.6)',advance="no") Ts, x(1,i), y(1,i)
end do
write(16,'(I10,2F14.6)') Ts, x(1,numMol), y(1,numMol)
do Ip=1, numParticles-1
    write(17,'(F14.6)',advance="no") contactState(1,Ip)*bondLength
end do
write(17,'(F14.6)') contactState(1,numParticles)*bondLength
do Ip=1, numParticles
    if(contactState(1,Ip)>=1.0_DP) then
        jContact=int(contactState(1,Ip))
        interceptContact=contactState(1,Ip)-dble(jContact)
        if(jContact>=numBeads) then
            xContact=x(jContact,1)+interceptContact*(x(jContact,1)-x(jContact-1,1))
            yContact=y(jContact,1)+interceptContact*(y(jContact,1)-y(jContact-1,1))
            zContact=z(jContact,1)+interceptContact*(z(jContact,1)-z(jContact-1,1))
        else
            xContact=x(jContact,1)+interceptContact*(x(jContact+1,1)-x(jContact,1))
            yContact=y(jContact,1)+interceptContact*(y(jContact+1,1)-y(jContact,1))
            zContact=z(jContact,1)+interceptContact*(z(jContact+1,1)-z(jContact,1))
        end if
        write(20,'(2I10,7F14.6,I2)') Ts, Ip, contactState(1,Ip), xContact, yContact, zContact, xp(Ip), yp(Ip), zp(Ip), motorSpecies(Ip)
    end if
end do
do Ip=1, numParticles
    if(contactState(1,Ip)>=1.0_DP)then
        write(22,'(2I10,7F14.6)') Ts, Ip, contactState(1,Ip), fMotorX(Ip), fMotorY(Ip), fMotorZ(Ip)
    end if
end do
write(23,*) Ts, areaCounter, eraseCounter
write(outFileNameF,'(A,I3.3,A,I7.7,A)') 'filamentA', iAssay, 'T', outputFileCounter, '.vtk'
open(33,file=outFileNameF)
write(33,'(A)') "# vtk dataFile version 2.0"
write(33,'(A, I7)') "filament: outputFileCounter = ", outputFileCounter
write(33,'(A)') "ASCII"
write(33,'(A)') "dataset unstructured grid"
write(33,'(A, I10, A)') "points", numMol*numBeads, "float"
do i=1, numMol
    do j=1, numBeads
        write(33,'(3F10.5)') x(j,i), y(j,i), z(j,i)
    end do
end do
write(33,'(A, I10, I10)') "cells", numMol, (1+numBeads)*numMol
do i=1, numMol
    write(33,*) numBeads, (j+(i-1)*numBeads-1, j=1, numBeads)
end do
write(33,'(A, I10)') "cellTypes", numMol
do i=1, numMol
    write(33,'(I3)') 4
end do
close(33)
write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'particleA', iAssay, 'T', outputFileCounter, '.vtk'
open(34,file=outFileNameP)
write(34,'(A)') "# vtk dataFile version 2.0"
write(34,'(A, I7)') "particle: outputFileCounter = ", outputFileCounter
write(34,'(A)') "ASCII"
write(34,'(A)') "dataset unstructured grid"
write(34,'(A, I10, A)') "points", numParticles, "float"
do Ip=1, numParticles
    write(34,'(3F10.5)') xp(Ip), yp(Ip), zp(Ip)
end do
write(34,'(A, I10, I10)') "cells", numParticles, 2*numParticles
do Ip=1, numParticles
    write(34,'(I3, I10)') 1, Ip-1
end do
write(34,'(A, I10)') "cell types", numParticles
do Ip=1, numParticles
    write(34,'(I3)') 1
end do
close(34)
write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'fPlusEndA', iAssay, 'T', outputFileCounter, '.vtk'
open(35, file=outFileNameP)
write(35,'(A)') "# vtk dataFile version 2.0"
write(35,'(A, I7)') "filament plus ends: outputFileCounter = ", outputFileCounter
write(35,'(A)') "ASCII"
write(35,'(A)') "dataset unstructured grid"
write(35,'(A, I10, A)') "points", numMol, "float"
do i=1, numMol
    write(35,'(3F10.5)') x(numBeads,i), y(numBeads,i), z(numBeads,i)
end do
write(35,'(A, I10, I10)') "cells", numMol, 2*numMol
do i=1, numMol
    write(35,'(I3, I10)') 1, i-1
end do
write(35,'(A, I10)') "cell types ", numMol
do i=1, numMol
    write(35,'(I3)') 1
end do
close(35)
write(outFileName,'(A,I3.3,A,I7.7,A)') 'contactStatesA', iAssay, 'T', outputFileCounter, '.txt'
open(36,file=outFileName)
do Ip=1, numParticles
    do ii=1,numMol-1
        write(36,'(F10.5)',advance="no") contactState(ii,Ip)
    end do
    write(36,'(F10.5)') contactState(numMol,Ip)
end do
close(36)
write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'intParticleA', iAssay, 'T', outputFileCounter, '.vtk'
open(37,file=outFileNameP)
counterBuffer=0
do Ip=1, numParticles
    if(sum(contactState(:,Ip))>1) then
        counterBuffer=counterBuffer+1
        xpBuffer(counterBuffer)=xp(Ip)
        ypBuffer(counterBuffer)=yp(Ip)
        zpBuffer(counterBuffer)=zp(Ip)
    end if
end do
write(37,'(A)') "# vtk dataFile version 2.0"
write(37,'(A, I7)') "particle: outputFileCounter = ", outputFileCounter
write(37,'(A)') "ASCII"
write(37,'(A)') "dataset unstructured grid"
write(37,'(A, I10, A)') "points", counterBuffer, "float"
do Ip=1,counterBuffer
    write(37,'(3F10.5)') xpBuffer(Ip), ypBuffer(Ip), zpBuffer(Ip)
end do
write(37,'(A, I10, I10)') "cells", counterBuffer, 2*counterBuffer
do Ip=1, counterBuffer
    write(37, '(I3,I10)') 1, Ip-1
end do
write(37,'(A,I10)') "cell types", counterBuffer
do Ip=1, counterBuffer
    write(37,'(I3)') 1
end do
close(37)
do iSpecies=1, numSpecies
    write(outFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'motorSpecie', iSpecies, 'A', iAssay, 'T', outputFileCounter, '.vtk'
    open(38,file=outFileName)
    counterBuffer=0
    do Ip=1, numParticles
        if(motorSpecies(Ip)==iSpecies) then
            counterBuffer=counterBuffer+1
            xpBuffer(counterBuffer)=xp(Ip)
            ypBuffer(counterBuffer)=yp(Ip)
            zpBuffer(counterBuffer)=zp(Ip)
        end if
    end do
    write(38,'(A)') "# vtk dataFile version 2.0"
    write(38,'(A, I7)') "particle: outputFileCounter = ", outputFileCounter
    write(38,'(A)') "ASCII"
    write(38,'(A)') "dataset unstructured grid"
    write(38,'(A,I10,A)') "points", counterBuffer, "float"
    do Ip=1, counterBuffer
        write(38,'(3F10.5)') xpBuffer(Ip), ypBuffer(Ip), zpBuffer(Ip)
    end do
    write(38,'(A,I10,I10)') "cells", counterBuffer, 2*counterBuffer
    do Ip=1, counterBuffer
        write(38,'(I3, I10)') 1, Ip-1
    end do
    write(38,'(A, I10)') "cell types", counterBuffer
    do Ip=1, counterBuffer
        write(38,'(I3)') 1
    end do
    close(38)
end do
loopTime: do Ts=1, numTimeStep
do i=1, numMol
    call normRNDVector(numBeads, normRandVectorBeads)
    do j=1,numBeads
        xNormRandVectorBeads(j,i)=normRandVectorBeads(j)
    end do
end do
do i=1, numMol
    call normRNDVector(numBeads,normRandVectorBeads)
    do j=1,numBeads
        yNormRandVectorBeads(j,i)=normRandVectorBeads(j)
    end do
end do
do i=1, numMol
    call normRNDVector(numBeads,normRandVectorBeads)
    do j=1,numBeads
        zNormRandVectorBeads(j,i)=normRandVectorBeads(j)
    end do
end do
loopMol: do i=1, numMol
do j=1, numBeads
    xi(j)=x(j,i)
end do
do j=1, numBeads
    yi(j)=y(j,i)
end do
do j=1, numBeads
    zi(j)=z(j,i)
end do
do j=1, numBeads
    fiX(j)=fX(j,i)
end do
do j=1, numBeads
    fiY(j)=fY(j,i)
end do
do j=1, numBeads
    fiZ(j)=fZ(j,i)
end do
do j=1, numBeads
    xiNormRandVectorBeads(j)=xNormRandVectorBeads(j,i)
end do
do j=1, numBeads
    yiNormRandVectorBeads(j)=yNormRandVectorBeads(j,i)
end do
do j=1, numBeads
    ziNormRandVectorBeads(j)=zNormRandVectorBeads(j,i)
end do
call normRNDVector(numBeads,normRandVectorBeads)
xiTemp=xi+fiX/gammaBead*dt+xiNormRandVectorBeads*dsqrt(2.0_DP*dBead*dt)
call normRNDVector(numBeads,normRandVectorBeads)
yiTemp=yi+fiY/gammaBead*dt+yiNormRandVectorBeads*dsqrt(2.0_DP*dBead*dt)
call normRNDVector(numBeads,normRandVectorBeads)
ziTemp=zi+fiZ/gammaBead*dt+ziNormRandVectorBeads*dsqrt(2.0_DP*dBead*dt)
numIteration=0
stateConstraint=.FALSE.
stateConfinement=.FALSE.
iterationConstraintConfinement: do
if(stateConstraint .and. stateConfinement) exit iterationConstraintConfinement
if(numIteration > maxNumIteration) then
    write(*,*) "too many iterations!"
    exit iterationConstraintConfinement
end if
numIteration=numIteration+1
stateConfinement=.TRUE.
loopBeadConfinement: do j=1, numBeads
call confinement(xiTemp(j), yiTemp(j), ziTemp(j), xi(j), yi(j), zi(j), stateConfinement)
end do loopBeadConfinement
stateConstraint=.TRUE.
loopBondConstraint: do j=1, numBeads-1
call constraint(xiTemp(j), yiTemp(j), ziTemp(j), xiTemp(j+1), yiTemp(j+1), ziTemp(j+1), xi(j), yi(j), zi(j), xi(j+1), yi(j+1), zi(j+1), stateConstraint)
end do loopBondConstraint
end do iterationConstraintConfinement
do j=1, numBeads
    x(j,i)=xiTemp(j)
end do
do j=1, numBeads
    y(j,i)=yiTemp(j)
end do
do j=1, numBeads
    z(j,i)=ziTemp(j)
end do
end do loopMol
do Ip=1, numParticles
    if(motorSpecies(Ip) == 1) then
        call motorStateConv(x,y,z,fMotorX(Ip),fMotorY(Ip),fMotorZ(Ip),contactState(:,Ip),releaseADP(Ip))
    end if
end do
if(mod(Ts,areaRenewDiv) == 0) then
    call renewMotorPopulation(x,y,xp,yp,zp,contactState,candidateList,addedParticleNum,areaCounter,eraseCounter,areaOriginX,areaOriginY,areaOriginUx,areaOriginUy,numParticles,motorSpecies)
end if
tempContact=0.0_DP
call checkMotorFilamentContact(x,y,z,xp,yp,zp,numParticles,candidateList,contactState,tempContact)
call motorBinding(numParticles,contactState,tempContact)
do Ip=1, numParticles
    if(motorSpecies(Ip)==1) then
        call motorMove(contactState(:,Ip),tempContact(:,Ip))
    end if
end do
call calculateForceMotor(x,y,z,xp,yp,zp,contactState,fMotorX,fMotorY,fMotorZ,numParticles,elongation)
do Ip=1, numParticles
    call motorForcedDetachment(fMotorX(Ip),fMotorY(Ip), fMotorZ(Ip), contactState(:,Ip), elongation(Ip))
end do
call calculateForceBead(contactState, forceBeadX, forceBeadY, forceBeadZ, fMotorX, fMotorY, fMotorZ, numParticles)
fX=forceBending(x)+forceBeadX
fY=forceBending(y)+forceBeadY
fZ=forceBending(z)+forceBeadZ
fX=forceBending(x)+forceBeadX+extFx(Ts,x,y,z)
fY=forceBending(y)+forceBeadY+extFy(Ts,x,y,z)
fZ=forceBending(z)+forceBeadZ+extFz(Ts,x,y,z)
if(mod(Ts,outputDiv)==0) then
    outputFileCounter=outputFileCounter+1
    do j=1, numBeads
        write(11, '(I10)', advance="no") Ts
        do i=1, numMol-1
            write(11,'(3F25.15)', advance="no") x(j,i), y(j,i), z(j,i)
        end do
        write(11,'(3F25.15)') x(j,numMol), y(j,numMol), z(j,numMol)
    end do
    do i=1, numMol-1
        write(15,'(3F14.6)', advance="no") x(numBeads,i), y(numBeads,i), z(numBeads,i)
    end do
    write(15,'(3F14.6)') x(numBeads,numMol), y(numBeads,numMol), z(numBeads,numMol)
    do i=1, numMol-1
        write(16,'(I10, 2F14.6)', advance="no") Ts, x(1,i), y(1,i)
    end do
    write(16,'(I10, 2F14.6)') Ts, x(1,numMol), y(1,numMol)
    do Ip=1, numParticles-1
        write(17,'(F14.6)',advance="no") contactState(1,Ip)*bondLength
    end do
    write(17,'(F14.6)') contactState(1,numParticles)*bondLength
    do Ip=1, numParticles
        if(contactState(1,Ip)>=1.0_DP) then
            jContact=int(contactState(1,Ip))
            interceptContact=contactState(1,Ip)-dble(jContact)
            if(jContact>=numBeads) then
                xContact=x(jContact,1)+interceptContact*(x(jContact,1)-x(jContact-1,1))
                yContact=y(jContact,1)+interceptContact*(y(jContact,1)-y(jContact-1,1))
                zContact=z(jContact,1)+interceptContact*(z(jContact,1)-z(jContact-1,1))
            else
                xContact=x(jContact,1)+interceptContact*(x(jContact+1,1)-x(jContact,1))
                yContact=y(jContact,1)+interceptContact*(y(jContact+1,1)-y(jContact,1))
                zContact=z(jContact,1)+interceptContact*(z(jContact+1,1)-z(jContact,1))
            end if
            write(20,'(2I10, 7F14.6, I2)') Ts, Ip, contactState(1,Ip), xContact, yContact, zContact, xp(Ip), yp(Ip), zp(Ip), motorSpecies(Ip)
        end if
    end do
    do Ip=1, numParticles
        if(contactState(1,Ip)>=1.0_DP) then
            write(22,'(2I10, 7F14.6)') Ts, Ip, contactState(1,Ip), fMotorX(Ip), fMotorY(Ip), fMotorZ(Ip)
        end if
    end do
    write(23,*) Ts, areaCounter, eraseCounter
    write(outFileNameF,'(A,I3.3,A,I7.7,A)') 'filamentA', iAssay, 'T', outputFileCounter, '.vtk'
    open(33,file=outFileNameF)
    write(33,'(A)') "# vtk dataFile version 2.0"
    write(33,'(A, I7)') "Filament: outputFileCounter = ", outputFileCounter
    write(33,'(A)') "ASCII"
    write(33,'(A)') "dataset unstructured grid"
    write(33,'(A, I10, A)') "points", numMol*numBeads, "float"
    do I=1, numMol
        do j=1, numBeads
            write(33,'(3F10.5)') x(j,i), y(j,i), z(j,i)
        end do
    end do
    write(33,'(A, I10, I10)') "cells", numMol, (1+numBeads)*numMol
    do I=1, numMol
        write(33,*) numBeads, (j+(i-1)*numBeads-1,j=1,numBeads)
    end do
    write(33,'(A, I10)') "cell types", numMol
    do i=1, numMol
        write(33,'(I3)') 4
    end do
    close(33)
    write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'particleA', iAssay, 'T', outputFileCounter, '.vtk'
    open(34,file=outFileNameP)
    write(34,'(A)') "# vtk DataFile version 2.0"
    write(34,'(A,I7)') "particle: outputFileCounter = ", outputFileCounter
    write(34,'(A)') "ASCII"
    write(34,'(A)') "dataset unstructured grid"
    write(34,'(A, I10, A)') "points", numParticles, "float"
    do Ip=1, numParticles
        write(34,'(3F10.5)') xp(Ip), yp(Ip), zp(Ip)
    end do
    write(34,'(A, I10, I10)') "cells", numParticles, 2*numParticles
    do Ip=1, numParticles
        write(34,'(I3, I10)') 1, Ip-1
    end do
    write(34,'(A, I10)') "cell types", numParticles
    do Ip=1, numParticles
        write(34,'(I3)') 1
    end do
    close(34)
    write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'mtPlusEndA', iAssay, 'T', outputFileCounter, '.vtk'
    open(35,file=outFileNameP)
    write(35,'(A)') "# vtk dataFile version 2.0"
    write(35,'(A, I7)') "MT plus ends: outputFileCounter = ", outputFileCounter
    write(35,'(A)') "ASCII"
    write(35,'(A)') "dataset unstructured grid"
    write(35,'(A, I10, A)') "points", numMol, " float"
    do i=1, numMol
        write(35,'(3F10.5)') x(numBeads,i), y(numBeads,i), z(numBeads,i)
    end do
    write(35,'(A,I10,I10)') "cells", numMol, 2*numMol
    do i=1, numMol
        write(35,'(I3, I10)') 1, i-1
    end do
    write(35,'(A,I10)') "cell types", numMol
    do i=1, numMol
        write(35,'(I3)') 1
    end do
    close(35)
    write(outFileName,'(A,I3.3,A,I7.7,A)') 'contactStatesA', iAssay, 'T', outputFileCounter, '.txt'
    open(36,file=outFileName)
    do Ip=1, numParticles
        do ii=1, numMol-1
            write(36,'(F10.5)',advance="no") contactState(ii,Ip)
        end do
        write(36,'(F10.5)') contactState(numMol,Ip)
    end do
    close(36)
    write(outFileNameP,'(A,I3.3,A,I7.7,A)') 'intParticleA', iAssay, 'T', outputFileCounter, '.vtk'
    open(37, file=outFileNameP)
    counterBuffer=0
    do Ip=1, numParticles
        if (sum(contactState(:,Ip))>1) then
            counterBuffer=counterBuffer+1
            xpBuffer(counterBuffer)=xp(Ip)
            ypBuffer(counterBuffer)=yp(Ip)
            zpBuffer(counterBuffer)=zp(Ip)
        end if
    end do
    write(37,'(A)') "# vtk dataFile version 2.0"
    write(37,'(A,I7)') "particle: outputFileCounter = ", outputFileCounter
    write(37,'(A)') "ASCII"
    write(37,'(A)') "dataSet unstructured grid"
    write(37,'(A, I10, A)') "points", counterBuffer, " float"
    do Ip=1, counterBuffer
        write(37,'(3F10.5)') xpBuffer(Ip), ypBuffer(Ip), zpBuffer(Ip)
    end do
    write(37,'(A, I10, I10)') "cells", counterBuffer, 2*counterBuffer
    do Ip=1, counterBuffer
        write(37,'(I3, I10)') 1, Ip-1
    end do
    write(37,'(A, I10)') "cell types", counterBuffer
    do Ip=1, counterBuffer
        write(37,'(I3)') 1
    end do
    close(37)
    do iSpecies=1, numSpecies
        write(outFileName,'(A,I1.1,A,I3.3,A,I7.7,A)') 'motorSpecie', iSpecies, 'A', iAssay, 'T', outputFileCounter, '.vtk'
        open(38,file=outFileName)
        counterBuffer=0
        do Ip=1, numParticles
            if(motorSpecies(Ip)==iSpecies) then
                counterBuffer=counterBuffer+1
                xpBuffer(counterBuffer)=xp(Ip)
                ypBuffer(counterBuffer)=yp(Ip)
                zpBuffer(counterBuffer)=zp(Ip)
            end if
        end do
        write(38,'(A)') "# vtk dataFile version 2.0"
        write(38,'(A, I7)') "particle: outputFileCounter = ", outputFileCounter
        write(38,'(A)') "ASCII"
        write(38,'(A)') "dataSet unstructured grid"
        write(38,'(A, I10, A)') "points", counterBuffer, " float"
        do Ip=1, counterBuffer
            write(38,'(3F10.5)') xpBuffer(Ip), ypBuffer(Ip), zpBuffer(Ip)
        end do
        write(38,'(A, I10, I10)') "cells", counterBuffer, 2*counterBuffer
        do Ip=1, counterBuffer
            write(38,'(I3, I10)') 1, Ip-1
        end do
        write(38,'(A,I10)') "cell types", counterBuffer
        do Ip=1, counterBuffer
            write(38,'(I3)') 1
        end do
        close(38)
    end do
end if
if(x(1,1)>xLimit) exit loopTime
end do loopTime
close(11)
close(15)
close(16)
close(17)
close(20)
end do loopAssay
close(10)
contains
subroutine normRNDNum(NRN)
    use parameters, only : DP
    use mtmod, only : grnd
    real(kind=DP), intent(out) :: NRN
    real(kind=DP) :: ur1, ur2
    ur1=grnd()
    if(ur1==0.0_DP) ur1=grnd()
    ur2=grnd()
    NRN=dsqrt(-2.0_DP*dlog(ur1)) * dcos(2.0_DP*pi*(ur2))
end subroutine normRNDNum
subroutine normRNDVector(nDim,NRV)
    use parameters, only : DP
    use mtmod, only : grnd
    integer, intent(in) :: nDim
    real(kind=DP), dimension(nDim), intent(out) :: NRV
    integer :: i
    real(kind=DP) :: ur1, ur2, NR
    do i=1, nDim
        ur1=grnd()
        if(ur1==0.0_DP) ur1=grnd()
        ur2=grnd()
        NR=dsqrt(-2.0_DP*dlog(ur1)) * dcos(2.0_DP*pi*(ur2))
        NRV(i)=NR
    end do
end subroutine normRNDVector
subroutine makeList(x,y,z,xp,yp,zp,contactState,candidateList)
    use parameters, only : DP, maxNumParticles, numMol, numBeads, bondLength, merginList
    real(kind=DP), dimension(numBeads,numMol),intent(in) :: x,y,z
    real(kind=DP), dimension(maxNumParticles),intent(in) :: xp,yp,zp
    real(kind=DP), dimension(numMol,maxNumParticles), intent(in) :: contactState
    integer, dimension(numMol*numBeads, maxNumParticles), intent(out) :: candidateList
    real(kind=DP) :: intercept, sqDistance
    candidateList=0
    do Ip=1, numParticles
        loopFilamentMakelist: do ii=1, numMol
        if(contactState(ii,Ip)>=1.0_DP) cycle loopFilamentMakelist
        do jj=1, numBeads-1
            intercept=(xp(Ip)-x(jj,ii))*(x(jj+1,ii)-x(jj,ii))+(yp(Ip)-y(jj,ii))*(y(jj+1,ii)-y(jj,ii))+(zp(Ip)-z(jj,ii))*(z(jj+1,ii)-z(jj,ii))
            intercept=intercept/bondLength
            if((intercept>bondLength+merginList) .or. (intercept<0.0_DP - merginList)) cycle
            sqDistance=(xp(Ip)-x(jj,ii))*(xp(Ip)-x(jj,ii))+(yp(Ip)-y(jj,ii))*(yp(Ip)-y(jj,ii))+(zp(Ip)-z(jj,ii))*(zp(Ip)-z(jj,ii))-intercept**2
            if(sqDistance<=merginList**2) then
                candidateList(numBeads*(ii-1)+jj,Ip)=1
            end if
        end do
    end do loopFilamentMakelist
end do
end subroutine makeList
subroutine checkMotorFilamentContact(x,y,z,xp,yp,zp,numParticles,candidateList,contactState,tempContact)
    use parameters, only : DP, maxNumParticles, numMol, numBeads, bondLength, merginList, stepSize
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: x,y,z
    real(kind=DP), dimension(maxNumParticles), intent(in) :: xp,yp,zp
    real(kind=DP), dimension(numMol, maxNumParticles), intent(inout) :: contactState, tempContact
    integer(kind=range15), intent(in) :: numParticles
    integer, dimension(numMol*numBeads, maxNumParticles), intent(in) :: candidateList
    real(kind=DP) :: intercept, sqDistance
    integer(kind=range15) :: Ip
    integer :: ii,jj
    loopParticleContactCheck: do Ip=1, numParticles
    if(sum(candidateList(:,Ip))==0) cycle loopParticleContactCheck
    loopFilamentContactCheck: do ii=1, numMol
    if((candidateList(ii,Ip)>=1.0_DP) .OR. (contactState(ii,Ip) == -1.0_DP)) cycle loopFilamentContactCheck
    if(sum(candidateList((numBeads*(ii-1)+1) : (numBeads*ii),Ip))==0) cycle loopFilamentContactCheck
    do jj=1, numBeads-1
        if(candidateList(numBeads*(ii-1)+jj, Ip)==1) then
            intercept = (xp(Ip)-x(jj,ii))*(x(jj+1,ii) - x(jj,ii)) + (yp(Ip)-y(jj,ii))*(y(jj+1,ii)-y(jj,ii))+(zp(Ip)-z(jj,ii))*(z(jj+1,ii)-z(jj,ii))
            intercept = intercept/bondLength
            if((intercept>bondLength) .OR. (intercept<0.0_DP)) cycle
            sqDistance = (xp(Ip) - x(jj,ii))*(xp(Ip)-x(jj,ii))+(yp(Ip)-y(jj,ii))*(yp(Ip)-y(jj,ii))+(zp(Ip)-z(jj,ii))*(zp(Ip)-z(jj,ii))-intercept**2
            if(sqDistance<=captureRadius**2) then
                tempContact(ii,Ip)=dble(jj)+intercept/bondLength
            end if
        end if
    end do
end do loopFilamentContactCheck
end do loopParticleContactCheck
end subroutine checkMotorFilamentContact
subroutine motorBinding(numParticles, contactState, tempContact)
    use parameters, only : DP, numMol, dt, ka
    use mtmod, only : grnd
    real(kind=DP), dimension(numMol, maxNumParticles), intent(inout) :: contactState, tempContact
    integer(kind=range15), intent(in) :: numParticles
    real(kind=DP), dimension(maxNumParticles) :: urArray
    integer(kind=range15) :: Ip
    integer :: ii
    do Ip=1, numParticles
        urArray(Ip) = grnd()
    end do
    loopMotorContactCheck: do Ip=1, numParticles
    loopFilamentContactCheck: do ii=1, numMol
    if((tempContact(ii,Ip)>=1.0_DP) .and. (int(contactState(ii,Ip))==0)) then
        if(urArray(Ip)>ka*dt) then
            tempContact(ii,Ip)=0.0_DP
        end if
    end if
end do loopFilamentContactCheck
end do loopMotorContactCheck
end subroutine motorBinding
subroutine motorMove(contactStateIp, tempContactIp)
    use parameters, only : DP, numMol, dt, ka
    use mtmod, only : grnd
    real(kind=DP), dimension(numMol), intent(inout) :: contactStateIp, tempContactIp
    integer :: ii
    loopFilamentContactCheck: do ii=1, numMol
    if(tempContactIp(ii)>=1.0_DP) then
        contactStateIp(ii)=tempContactIp(ii)+stepSize/bondLength
    end if
end do loopFilamentContactCheck
end subroutine motorMove
subroutine calculateForceMotor(x,y,z,xp,yp,zp,contactState,fMotorX,fMotorY,fMotorZ,numParticles,elongation)
    use parameters, only : DP, maxNumParticles, numMol, captureRadius, k, fMotorDetach
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: x,y,z
    real(kind=DP), dimension(maxNumParticles), intent(in) :: xp,yp,zp
    real(kind=DP), dimension(numMol, maxNumParticles), intent(in) :: contactState
    real(kind=DP), dimension(maxNumParticles), intent(out) :: fMotorX, fMotorY, fMotorZ
    integer(kind=range15), intent(in) :: numParticles
    real(kind=DP) :: interceptContact, xContact, yContact, zContact
    real(kind=DP), dimension(maxNumParticles), intent(out) :: elongation
    integer(kind=range15) :: Ip
    integer :: ii, jContact
    fMotorX=0.0_DP
    fMotorY=0.0_DP
    fMotorZ=0.0_DP
    do Ip=1, numParticles
        scanContactFilament: do ii=1, numMol
        if(contactState(ii,Ip)>=1.0_DP) then
            jContact=int(contactState(ii,Ip))
            interceptContact=contactState(ii,Ip) - dble(jContact)
            if(jContact>=numBeads) then
                xContact=x(numBeads,ii)+interceptContact*(x(numBeads,ii)-x(numBeads-1,ii))
                yContact=y(numBeads,ii)+interceptContact*(y(numBeads,ii)-y(numBeads-1,ii))
                zContact=z(numBeads,ii)+interceptContact*(z(numBeads,ii)-z(numBeads-1,ii))
            else
                xContact=x(jContact,ii)+interceptContact*(x(jContact+1,ii)-x(jContact,ii))
                yContact=y(jContact,ii)+interceptContact*(y(jContact+1,ii)-y(jContact,ii))
                zContact=z(jContact,ii)+interceptContact*(z(jContact+1,ii)-z(jContact,ii))
            end if
            elongation(Ip)=(xp(Ip)-xContact)**2 + (yp(Ip) - yContact)**2 + (zp(Ip) - zContact)**2
            elongation(Ip) = dsqrt(elongation(Ip))
            fMotorX(Ip) = k*elongation(Ip)*(xp(Ip)-xContact)/elongation(Ip)
            fMotorY(Ip) = k*elongation(Ip)*(yp(Ip)-yContact)/elongation(Ip)
            fMotorZ(Ip) = k*elongation(Ip)*(zp(Ip)-zContact)/elongation(Ip)
        end if
    end do scanContactFilament
end do
end subroutine calculateForceMotor
subroutine motorForcedDetachment(fMotorXiP, fMotorYiP, fMotorZiP, contactStateIp, elongationIp)
    use parameters, only : DP, maxNumParticles, numMol, k, fMotorDetach
    real(kind=DP), intent(inout) :: fMotorXiP, fMotorYiP, fMotorZiP
    real(kind=DP), dimension(numMol), intent(inout) :: contactStateIp
    real(kind=DP), intent(in) :: elongationIp
    integer :: ii
    if(k*elongationIp>=fMotorDetach) then
        fMotorXiP=0.0_DP
        fMotorYiP=0.0_DP
        fMotorZiP=0.0_DP
        scanContactFilament: do ii=1, numMol
        if(contactStateIp(ii)>=1.0_DP) then
            contactStateIp(ii) = -1.0_DP
        end if
    end do scanContactFilament
end if
end subroutine motorForcedDetachment
subroutine calculateForceBead(contactState, forceBeadX, forceBeadY, forceBeadZ, fMotorX, fMotorY, fMotorZ, numParticles)
    use parameters, only : DP, maxNumParticles, numMol
    real(kind=DP), dimension(numMol, maxNumParticles), intent(in) :: contactState
    real(kind=DP), dimension(numBeads, numMol), intent(out) :: forceBeadX, forceBeadY, forceBeadZ
    real(kind=DP), dimension(maxNumParticles), intent(in) :: fMotorX, fMotorY, fMotorZ
    integer(kind=range15), intent(in) :: numParticles
    real(kind=DP) :: interceptContact
    integer(kind=range15) :: Ip
    integer :: ii, jContact
    forceBeadX = 0.0_DP
    forceBeadY = 0.0_DP
    forceBeadZ = 0.0_DP
    do Ip=1, numParticles
        scanContactFilament: do ii=1, numMol
        if(contactState(ii,Ip)>=1.0_DP) then
            jContact=int(contactState(ii,Ip))
            interceptContact=contactState(ii,Ip) - dble(jContact)
            if(jContact>=numBeads) then
                forceBeadX(numBeads,ii)=forceBeadX(numBeads,ii)+fMotorX(Ip)
                forceBeadY(numBeads,ii)=forceBeadY(numBeads,ii)+fMotorY(Ip)
                forceBeadZ(numBeads,ii)=forceBeadZ(numBeads,ii)+fMotorZ(Ip)
            else
                forceBeadX(jContact,ii)=forceBeadX(jContact,ii)+(1.0_DP - interceptContact)*fMotorX(Ip)
                forceBeadY(jContact,ii)=forceBeadY(jContact,ii)+(1.0_DP - interceptContact)*fMotorY(Ip)
                forceBeadZ(jContact,ii)=forceBeadZ(jContact,ii)+(1.0_DP - interceptContact)*fMotorZ(Ip)
                forceBeadX(jContact+1,ii)=forceBeadX(jContact+1,ii)+interceptContact*fMotorX(Ip)
                forceBeadY(jContact+1,ii)=forceBeadY(jContact+1,ii)+interceptContact*fMotorY(Ip)
                forceBeadZ(jContact+1,ii)=forceBeadZ(jContact+1,ii)+interceptContact*fMotorZ(Ip)
            end if
        end if
    end do scanContactFilament
end do
end subroutine calculateForceBead
subroutine constraint(xiTempJ0, yiTempJ0, ziTempJ0, xiTempJ1, yiTempJ1, ziTempJ1, xiJ0, yiJ0, ziJ0, xiJ1, yiJ1, ziJ1, stateConstraint)
    use parameters, only : DP, tol, bondLength
    real(kind=DP), intent(inout) :: xiTempJ0, xiTempJ1, yiTempJ0, yiTempJ1, ziTempJ0, ziTempJ1
    real(kind=DP), intent(in) :: xiJ0, xiJ1, yiJ0, yiJ1, ziJ0, ziJ1
    logical, intent(inout) :: stateConstraint
    real(kind=DP) :: xiTempAB, yiTempAB, ziTempAB, xiAB, yiAB, ziAB, beadsDisSq, diffSq, gAB, dX, dY, dZ
    xiTempAB=xiTempJ1-xiTempJ0
    yiTempAB=yiTempJ1-yiTempJ0
    ziTempAB=ziTempJ1-ziTempJ0
    xiAB=xiJ1-xiJ0
    yiAB=yiJ1-yiJ0
    ziAB=ziJ1-ziJ0
    beadsDisSq=xiTempAB**2 + yiTempAB**2 + ziTempAB**2
    diffSq=bondLength**2 - beadsDisSq
    if(dabs(diffSq)>2.0_DP*tol*bondLength**2) then
        stateConstraint= .false.
        gAB=diffSq/4.0_DP/(xiTempAB*xiAB+yiTempAB*yiAB+ziTempAB*ziAB)
        dX=gAB*xiAB
        dY=gAB*yiAB
        dZ=gAB*ziAB
        xiTempJ0=xiTempJ0-dX
        yiTempJ0=yiTempJ0-dY
        ziTempJ0=ziTempJ0-dZ
        xiTempJ1=xiTempJ1+dX
        yiTempJ1=yiTempJ1+dY
        ziTempJ1=ziTempJ1+dZ
    end if
end subroutine constraint
subroutine renewMotorPopulation(x,y,xp,yp,zp,contactState,candidateList,addedParticleNum,areaCounter,eraseCounter,areaOriginX,areaOriginY,areaOriginUx,areaOriginUy,numParticles,motorSpecies)
    use parameters, only : range15, DP, maxNumParticles, numMol, numBeads, bondLength, horizontalLength, verticalLength, motorDensity
    integer, dimension(maxNumParticles), intent(inout) :: motorSpecies
    real(kind=DP), dimension(numBeads,numMol), intent(in) :: x,y
    real(kind=DP), dimension(maxNumParticles), intent(inout) :: xp, yp, zp
    real(kind=DP), dimension(numMol, maxNumParticles), intent(inout) :: contactState
    integer, dimension(numMol*numBeads,maxNumParticles), intent(inout) :: candidateList
    integer(kind=range15), dimension(maxAreaNum), intent(inout) :: addedParticleNum
    integer, intent(inout) :: areaCounter, eraseCounter
    real(kind=DP), dimension(maxAreaNum), intent(inout) :: areaOriginX, areaOriginY, areaOriginUx, areaOriginUy
    integer(kind=range15), intent(inout) :: numParticles
    integer(kind=range15) :: Ip, addedParticleCounter
    integer :: jj,jjj
    real(kind=DP) :: xcm,ycm,xpNew,ypNew,ur1,ur2,ur3
    logical :: nearBoundary, inNewArea, outOldArea
    outOldArea=.true.
    do jjj=1, numBeads
        if((dabs((x(jjj,1)-areaOriginX(eraseCounter+1))*areaOriginUx(eraseCounter+1)+(y(jjj,1)-areaOriginY(eraseCounter+1))*areaOriginUy(eraseCounter+1)) <= 0.5_DP*horizontalLength + 0.5_DP) .and. &
        (dabs(-(x(jjj,1)-areaOriginX(eraseCounter+1))*areaOriginUy(eraseCounter+1)+(y(jjj,1)-areaOriginY(eraseCounter+1))*areaOriginUx(eraseCounter+1)) <= 0.5_DP*verticalLength+0.5_DP)) then 
            outOldArea = .false.
        end if
    end do
    if(outOldArea) then
        eraseCounter=eraseCounter+1
        do Ip=1, numParticles - addedParticleNum(eraseCounter)
            candidateList(:,Ip) = candidateList(:,Ip+addedParticleNum(eraseCounter))
        end do
        do Ip=1, numParticles-addedParticleNum(eraseCounter)
            contactState(:,Ip)=contactState(:,Ip+addedParticleNum(eraseCounter))
        end do
        do Ip=1, numParticles - addedParticleNum(eraseCounter)
            xp(Ip)=xp(Ip+addedParticleNum(eraseCounter))
        end do
        do Ip=1, numParticles - addedParticleNum(eraseCounter)
            yp(Ip)=yp(Ip+addedParticleNum(eraseCounter))
        end do
        do Ip=1, numParticles - addedParticleNum(eraseCounter)
            zp(Ip)=zp(Ip+addedParticleNum(eraseCounter))
        end do
        do Ip=1, numParticles-addedParticleNum(eraseCounter)
            motorSpecies(Ip)=motorSpecies(Ip+addedParticleNum(eraseCounter))
        end do
        numParticles=numParticles-addedParticleNum(eraseCounter)
    end if
    xcm=sum(x(:,1))/dble(numBeads)
    ycm=sum(y(:,1))/dble(numBeads)
    nearBoundary= .false.
    do jj=1, numBeads
        if(dabs((x(jj,1)-areaOriginX(areaCounter))*areaOriginUx(areaCounter)+(y(jj,1)-areaOriginY(areaCounter))*areaOriginUy(areaCounter))>=0.5_DP*horizontalLength-0.5_DP*dble(numBeads-1)*bondLength-1.0_DP) nearBoundary=.true.
        if(dabs(-(x(jj,1)-areaOriginX(areaCounter))*areaOriginUy(areaCounter)+(y(jj,1)-areaOriginY(areaCounter))*areaOriginUx(areaCounter))>=0.5_DP*verticalLength-0.5_DP*dble(numBeads-1)*bondLength-1.0_DP) nearBoundary=.true.
    end do
    if(nearBoundary == .true.) then
        areaCounter=areaCounter+1
        areaOriginX(areaCounter)=xcm
        areaOriginY(areaCounter)=ycm
        areaOriginUx(areaCounter)=(x(1,1)-x(2,1))/dsqrt((x(1,1)-x(2,1))*(x(1,1)-x(2,1))+(y(1,1)-y(2,1))*(y(1,1)-y(2,1)))
        areaOriginUy(areaCounter)=(y(1,1)-y(2,1))/dsqrt((x(1,1)-x(2,1))*(x(1,1)-x(2,1))+(y(1,1)-y(2,1))*(y(1,1)-y(2,1)))
        addedParticleCounter=0
        do Ip=1, int(motorDensity*horizontalLength*verticalLength,range15)
            ur1=grnd()
            ur2=grnd()
            xpNew=xcm+horizontalLength*(ur1-0.5_DP)*areaOriginUx(areaCounter)-verticalLength*(ur2-0.5_DP)*areaOriginUy(areaCounter)
            ypNew=ycm+horizontalLength*(ur1-0.5_DP)*areaOriginUy(areaCounter)+verticalLength*(ur2-0.5_DP)*areaOriginUx(areaCounter)
            inNewArea=.true.
            if(areaCounter==1) then
                if((dabs((xpNew-areaOriginX(areaCounter))*areaOriginUx(areaCounter)+(ypNew-areaOriginY(areaCounter))*areaOriginUy(areaCounter))<=0.5_DP*horizontalLength) .and. &
                (dabs(-(xpNew-areaOriginX(areaCounter))*areaOriginUy(areaCounter)+(ypNew-areaOriginY(areaCounter))*areaOriginUx(areaCounter))<=0.5_DP*verticalLength)) then
                    inNewArea = .false.
                end if
            else
                do iArea = eraseCounter+1, areaCounter-1
                    if((dabs((xpNew-areaOriginX(iArea))*areaOriginUx(iArea)+(ypNew-areaOriginY(iArea))*areaOriginUy(iArea))<=0.5_DP*horizontalLength) .and. &
                    (dabs(-(xpNew-areaOriginX(iArea))*areaOriginUy(iArea)+(ypNew-areaOriginY(iArea))*areaOriginUx(iArea))<=0.5_DP*verticalLength)) then
                        inNewArea = .false.
                    end if
                end do
            end if
            if(inNewArea == .true.) then
                 addedParticleCounter = addedParticleCounter+1
                xp(numParticles+addedParticleCounter)=xpNew
                yp(numParticles+addedParticleCounter)=ypNew
                zp(numParticles+addedParticleCounter)=0.0_DP
                candidateList(:,numParticles+addedParticleCounter)=1
                ur1 = grnd()
                if(ur1<=speciesRatio) then
                    motorSpecies(numParticles+addedParticleCounter)=1
                    ur2 = grnd()
                    if(ur2 <= 0.091) then
                        contactState(:,numParticles+addedParticleCounter)=-1.0_DP
                    else
                        contactState(:,numParticles+addedParticleCounter)=0.0_DP
                    end if
                else
                    motorSpecies(numParticles+addedParticleCounter)=2
                    contactState(:,numParticles+addedParticleCounter)=0.0_DP
                end if
            end if
        end do
        addedParticleNum(areaCounter)=addedParticleCounter
        numParticles = numParticles+addedParticleCounter
    end if
end subroutine renewMotorPopulation
subroutine motorStateConv(x,y,z,fMotorXiP,fMotorYiP,fMotorZiP,contactStateIp,releaseADPiP)
    use parameters
    use mtmod, only : grnd
    real(kind=DP), dimension(numBeads,numMol), intent(in) :: x,y,z
    integer, intent(inout) :: releaseADPiP
    real(kind=DP), dimension(numMol), intent(inout) :: contactStateIp
    real(kind=DP), dimension(numMol), intent(in) :: fMotorXiP, fMotorYiP, fMotorZiP
    real(kind=DP) :: ur, fMotorTan, kd
    integer(kind=range15) :: Ip
    integer :: ii, jContact
    do ii=1, numMol
        jContact=int(contactStateIp(ii))
        ur=grnd()
    select case(jContact)
    case(-1)
        if(ur<=khp*dt) then
            contactStateIp(ii)=0.0_DP
        end if
    case(0)
        if(ur<=khm*dt) then
            contactStateIp(ii)=-0.0_DP
        end if
    case(1:)
        if(releaseADPiP == 0) then
            if(jContact>=numBeads) then
                fMotorTan=(fMotorXiP(ii)*(x(jContact,ii)-x(jContact-1,ii))+fMotorYiP(ii)*(y(jContact,ii)-y(jContact-1,ii))+fMotorZiP(ii)*(z(jContact,ii)-z(jContact-1,ii)))/bondLength
            else
                fMotorTan=(fMotorXiP(ii)*(x(jContact+1,ii) - x(jContact,ii))+fMotorYiP(ii)*(y(jContact+1,ii) - y(jContact,ii))+fMotorZiP(ii)*(z(jContact+1,ii)-z(jContact,ii)))/bondLength
            end if
            kd=kd0*dexp(-fMotorTan*deltaX/kBT)
            if(ur<=kd*dt) then
                releaseADPiP=1
            end if
        else
            if(ur<=kt*ATP*dt) then
                releaseADPiP=0
                contactStateIp(ii)=-1.0_DP
            end if
        end if
    end select
end do
end subroutine motorStateConv
function forceBending(coordinate)
    use parameters, only : DP, numMol, numBeads, bondLength, EI
    real(kind=DP), dimension(numBeads,numMol) :: forceBending
    real(kind=DP), dimension(numBeads,numMol), intent(in) :: coordinate
    integer :: i,j
    real(kind=DP) :: F
    forceBending=0.0_DP
    do i=1, numMol
        do j=2, numBeads-1
            f=coordinate(j+1,i) - 2.0_DP*coordinate(j,i)+coordinate(j-1,i)
            forceBending(j-1,i)=forceBending(j-1,i)+(-1.0_DP)*f
            forceBending(j,i)=forceBending(j,i)+2.0_DP*f
            forceBending(j+1,i)=forceBending(j+1,i)+(-1.0_DP)*f
        end do
    end do
    forceBending = forceBending*EI/(bondLength**3)
end function forceBending
function boundParticlePosition(Ip,coordinate,contactState)
    use parameters,only : DP, numMol, numBeads, maxNumParticles
    real(kind=DP) :: boundParticlePosition, interceptContact
    integer :: ii, bondCounter, jContact
    integer, intent(in) :: Ip
    real(kind=DP), dimension(numBeads, numMol), intent(in) :: coordinate
    real(kind=DP), dimension(numMol, maxNumParticles) :: contactState
    boundParticlePosition=0.0_DP
    bondCounter=0
    do ii=1, numMol
        if(contactState(ii,Ip)>1.0_DP) then
            bondCounter=bondCounter+1
            jContact=int(contactState(ii,Ip))
            interceptContact=contactState(ii,Ip) - dble(jContact)
            boundParticlePosition = boundParticlePosition + coordinate(jContact,ii)+interceptContact*(coordinate(jContact+1,ii)-coordinate(jContact,ii))
        end if
    end do
    boundParticlePosition = boundParticlePosition/bondCounter
end function boundParticlePosition
end program beadRodPolymer
