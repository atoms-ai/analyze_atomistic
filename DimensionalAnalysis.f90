!     Calculation of particle width and height while impacting a flat rigid substrate
!     sumit.suresh@uconn.edu

PROGRAM PD

        IMPLICIT REAL*8(A-H,O-Z)
        INTEGER, PARAMETER:: KREAL=SELECTED_REAL_KIND(14,99) !14 sig figs, 10+-99
        INTEGER, PARAMETER:: LPMX = 50000000        !Maximum number of particles


        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: XX,YY,ZZ,CSP        !Coordinates

        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE  !Type & History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDGG,CNA

        INTEGER:: NAN, timenew                               !Number of particles
        REAL(KREAL):: TIME                          !Starting time
        REAL(KREAL):: XL,YL,ZL       !Size of the computational celL
        REAL(KREAL):: XCENTR,YCENTR,ZCENTR        !CENTRE of the computational cel

        CHARACTER *15 :: file_name

!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!	Use this to run
! ifort -shared-intel -mcmodel=large  Analysis2DBIN.f90


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEVEL OF COARSENING !!!!!!!!!!!!!!!!!!!!!!!!!!

	LOC = VAR_ACG
	NCG = LOC**3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Filename
	OPEN (UNIT = 14,FILE='/gpfs/scratchfs1/sumit/Coldspray/Rigid_impact/VAR_PSIZE/Impact/VAR_LOC/dump.imp1000mps.timestep')


        REWIND 14
        READ(14,*)
	      READ(14,*) TIMEST
	      READ(14,*)
	      READ(14,*) NAN
	      READ(14,*)

ALLOCATE (XX(NAN),YY(NAN),ZZ(NAN))
ALLOCATE (CSP(NAN))
ALLOCATE (KTYPE(NAN))
ALLOCATE (CNA(NAN))

        READ(14,*) XLO,XHI
        READ(14,*) YLO,YHI
        READ(14,*) ZLO,ZHI
        READ(14,*)
        DO J=1,NAN
         READ(14,*) IDdummy,KTYPE(J),XX(J),YY(J),ZZ(J),Q1X,Q1Y,Q1Z,CNA(J),CSP(J),StrX,StrY,StrZ     !!!!! LOOP 1
        end do

!       ##############################	CHANGE ts to t (ps)  ###########################

	dt = VAR_TS	!	delta t of integration
	timed = timest*dt
	timenew=INT(timed)

	! Write the integer into a string:

  WRITE(file_name, '(I3)')  timenew

	! Open the file with this name
	open(unit = 70, file = '../Dimensions/'//trim(adjustl(file_name))//'_Particledimensions.dat')

!       ##############################	SIMULATION CELL SIZE AND CENTERS  ###########################

	XL=XHI-XLO
	YL=YHI-YLO
	ZL=ZHI-ZLO

	XCENTR=(XHI+XLO)/2.0d0
	YCENTR=(YHI+YLO)/2.0d0
	ZCENTR=(ZHI+ZLO)/2.0d0

!       ##############################	PARTICLE DIMENSIONS, WIDTHS AND HEIGHTS  ###########################
!	2 lattice units square cross section (10A) and ~1 lattice unit (5A) in height for MD, scale as needed

	wcut=10.0*VAR_ACG
	hcut=5.0*VAR_ACG

	xlim1=xcentr-(wcut/2.0)
	xlim2=xcentr+(wcut/2.0)

	ylim1=ycentr-(wcut/2.0)
	ylim2=ycentr+(wcut/2.0)

	ParticleHmax = ZCENTR
	ParticleHmin = ZCENTR

	ParticleWmax = YCENTR
	ParticleWmin = YCENTR


  DO I=1,NAN

    IF ((CNA(I).EQ.5).AND.(CSP(I).EQ.0.0)) THEN
      CYCLE
    END IF

    IF((KTYPE(I).EQ.1).AND.(XX(I).LT.xlim2).AND.(XX(I).GT.xlim1).AND.(YY(I).LT.ylim2).AND.(YY(I).GT.ylim1)) THEN      !!!!! LOOP 2
		    IF(ParticleHmax.LE.ZZ(I)) ParticleHmax = ZZ(I)
		    IF(ParticleHmin.GE.ZZ(I)) ParticleHmin = ZZ(I)
    END IF

  END DO

  WRITE(70,*) "Particle Height in Angstroms ="
  WRITE(70,*) (ParticleHmax-ParticleHmin)

  zlim = ParticleHmin + hcut

  DO I=1,NAN                                                                                                                      !!!!! LOOP 3

    IF ((CNA(I).EQ.5).AND.(CSP(I).EQ.0)) THEN
      CYCLE
    END IF

	   IF((KTYPE(I).EQ.1).AND.(XX(I).LT.xlim2).AND.(XX(I).GT.xlim1).AND.(ZZ(I).LT.zlim)) THEN
		   IF(ParticleWmax.LE.YY(I)) ParticleWmax = YY(I)
		   IF(ParticleWmin.GE.YY(I)) ParticleWmin = YY(I)
     END IF

	END DO

  WRITE(70,*) "Particle Width in Angstroms ="
  WRITE(70,*) (ParticleWmax-ParticleWmin)


222 FORMAT(1x,I4,1x,I4,1x,F12.5,1x,F12.5,1x,I7,1x,F12.5,1x,F12.5,x,F12.5,1x,F12.5,1x,F12.5)
666	FORMAT(1x,F12.5,1x,F12.5,1x,F12.5,1x,I4,1x,I7,1x,F12.5,1x,F12.5,1x,F12.5,1x,F12.5,1x,F12.5)

798 FORMAT(1x,F11.4,1x,F11.4,1x,F11.4,1x,I7)

END PROGRAM PD
