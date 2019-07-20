module parameters
    implicit none
    integer, parameter :: range15 = selected_Int_Kind(15), DP = selected_Real_Kind(14), & 
    numAssay = 1, numMol = 1, numBeads = 13, numSpecies = 2, maxNumIteration = 1E6, & 
    areaRenewDiv = 2E4, outputDiv = 2E4, maxInteractinNum = 1, seed = 733
    integer(kind=range15), parameter :: numTimeStep = 2.0E6, timeStepEquil = 1E9, maxNumParticles = 1E6, timeForceOn = 1E9
    integer, parameter :: maxAreaNum = numTimeStep/areaRenewDiv + 1
    real(kind=DP), parameter :: pi = 3.14159265358979_DP, kBT = 0.0041_DP, dt = 5.0E-7_DP, tol = 1.0E-6_DP, &
    bondLength = 0.25_DP, EI = 0.073_DP, motorDensity = 3000.0_DP, k = 300.0_DP, stepSize = 1.0E-2_DP, &
    ka = 40.0_DP, kd0 = 350.0_DP, kt = 2.0_DP, ATP = 2000.0_DP, khp = 100.0_DP, khm = 10.0_DP, deltaX = -1.86E-3_DP, &
    fMotorDetach = 9.2_DP, captureRadius = 0.020_DP, merginList = bondLength*1.5_DP, trackWidth = 4.0_DP, & 
    speciesRatio = 0.60_DP, extForceDensity0 = 0.0_DP, xLimit = 500.0_DP, freq = 1.0_DP, maxStrain = 0.0_DP, & 
    strainRate = 0.05_DP, poissonRatio = 0.5_DP, horizontalLength = dble(numBeads-1)*bondLength + 1.0_DP, &
    verticalLength = dble(numBeads-1)*bondLength + 1.0_DP
end module parameters
