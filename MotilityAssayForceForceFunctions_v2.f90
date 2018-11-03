

MODULE UNIFORM_FORCE_X

USE PARAMETERS

IMPLICIT NONE

CONTAINS


!---External Force Calculation (X-coordinate)-------------------------

FUNCTION ExtF_X(TS, X, Y, Z)


	USE PARAMETERS, ONLY : DP, Range15, NumMol, NumBeads, BondLength, ExtForceDensity0, TimeForceON


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol) :: ExtF_X

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	INTEGER(KIND = Range15), INTENT(IN) :: TS

	INTEGER :: I, J

!	REAL(KIND = DP) :: F


	ExtF_X = 0.0_DP

	IF (TS >= TimeForceON) THEN
		DO I=1, NumMol
			ExtF_X(1,I) = 0.5_DP * ExtForceDensity0*BondLength
			DO J=2, NumBeads-1
				ExtF_X(J,I) = ExtForceDensity0 * BondLength
			END DO
			ExtF_X(NumBeads,I) = 0.5_DP * ExtForceDensity0 * BondLength
		END DO
	END IF

END FUNCTION ExtF_X


!---External Force Calculation (Y-coordinate)-------------------------

FUNCTION ExtF_Y(TS, X, Y, Z)


	USE PARAMETERS, ONLY : DP, Range15, NumMol, NumBeads, BondLength, ExtForceDensity0, TimeForceON


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol) :: ExtF_Y

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	INTEGER(KIND = Range15), INTENT(IN) :: TS

	INTEGER :: I, J

!	REAL(KIND = DP) :: F


	ExtF_Y = 0.0_DP


END FUNCTION ExtF_Y


!---External Force Calculation (Z-coordinate)-------------------------

FUNCTION ExtF_Z(TS, X, Y, Z)


	USE PARAMETERS, ONLY : DP, Range15, NumMol, NumBeads, BondLength, ExtForceDensity0, TimeForceON


	REAL(KIND = DP), DIMENSION(NumBeads, NumMol) :: ExtF_Z

	REAL(KIND = DP), DIMENSION(NumBeads, NumMol), INTENT(IN) :: X, Y, Z

	INTEGER(KIND = Range15), INTENT(IN) :: TS

	INTEGER :: I, J

!	REAL(KIND = DP) :: F


	ExtF_Z = 0.0_DP


END FUNCTION ExtF_Z



END MODULE UNIFORM_FORCE_X


!---------------------------------------------------------------------------------------------------


