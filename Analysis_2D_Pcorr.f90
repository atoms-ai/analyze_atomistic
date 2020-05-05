!     Calculation of Pressure/Temperature Distribution in two-dimensional space
!     sumit.suresh@uconn.edu

PROGRAM PT

        IMPLICIT REAL*8(A-H,O-Z)
        INTEGER, PARAMETER:: KREAL=SELECTED_REAL_KIND(14,99) !14 sig figs, 10+-99
        INTEGER, PARAMETER:: LPMX = 50000000        !Maximum number of particles
        INTEGER, PARAMETER:: NVD = 1000

        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: XX,YY,ZZ,TKE          !Coordinates
        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: StrX,StrY,StrZ,StrXY,StrYZ,StrZX    !STRAIN
        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: Q1X,Q1Y,Q1Z,VX,VY,VZ,CSP  !Velocities
        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: TPRESS,TEMP  !Velocities

        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE, KHIST !Type & History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDGG,KNCGG,CNA

        REAL(KREAL), DIMENSION(NVD,NVD):: ZPRESS,VMS,STRESS_X,STRESS_Y,STRESS_Z,STRESS_XY,STRESS_YZ,STRESS_ZX,ZVEL,YVEL,XVEL,VZAVE,VXAVE,VYAVE,VZSEC,VXSEC,VYSEC,TTKE,ZTEMP,CYMIN,CYMAX,YLEN,AvY
        REAL(KREAL), DIMENSION(NVD):: CZMIN,CZMAX,DYL
        REAL(KREAL), DIMENSION(NVD)::XStr,YStr,ZLength,YZMAX,YZMIN,AvZ,YWID

        INTEGER, DIMENSION(NVD,NVD):: KAT

        INTEGER:: NAN,KZ,NAS,KY,NY,NZ, timenew                               !Number of particles
        REAL(KREAL):: TIME,XKE,YKE,ZKE,mass,KB                          !Starting time
        REAL(KREAL):: XL,YL,ZL,DZL        !Size of the computational celL
        REAL(KREAL):: XCENTR,YXCENTR,ZXCENTR        !CENTRE of the computational cel
        REAL(KREAL):: AvPress, AvEN,AvMises,rvel,rpress

	CHARACTER *15 :: file_name

        INTEGER NBRMAX,inx,iny,inz,ncx,ncy,ncz
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!	Use this to run
! ifort -shared-intel -mcmodel=large  Analysis2DBIN.f90


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEVEL OF COARSENING !!!!!!!!!!!!!!!!!!!!!!!!!!

	LOC = VAR_ACG
	NCG = LOC**3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINE CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	AMUTOKG = 1.6605402D-27
	EVTOJOU = 1.60219D-19
	XJOUTOEV = 1.d0/EVTOJOU

	XMASS=4.480D-26 !in kg of Al - scale this with level of coarsening to correct for pressures
	XMASS=XMASS*NCG

  mass=4.480d0 ! For calculating temperatures from KE
  KB=1.380d0  ! Boltzmann constant

	ENUNIT = AMUTOKG*1.d4*XJOUTOEV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Filename
	OPEN (UNIT = 14,FILE='/gpfs/scratchfs1/sumit/Coldspray/Rigid_impact/VAR_PSIZE/Impact/VAR_LOC/dump.imp1000mps.timestep')


        REWIND 14
        READ(14,*)
	      READ(14,*) TIMEST
	      READ(14,*)
	      READ(14,*) NAN
	      READ(14,*)

ALLOCATE (XX(NAN),YY(NAN),ZZ(NAN),TKE(NAN))
ALLOCATE (StrX(NAN),StrY(NAN),StrZ(NAN),StrXY(NAN),StrYZ(NAN),StrZX(NAN))
ALLOCATE (Q1X(NAN),Q1Y(NAN),Q1Z(NAN),VX(NAN),VY(NAN),VZ(NAN),CSP(NAN))
ALLOCATE (TPRESS(NAN),TEMP(NAN))
ALLOCATE (KTYPE(NAN), KHIST(NAN))
ALLOCATE (IDGG(NAN),KNCGG(NAN),CNA(NAN))

        READ(14,*) XLO,XHI
        READ(14,*) YLO,YHI
        READ(14,*) ZLO,ZHI
        READ(14,*)
        DO J=1,NAN
         READ(14,*) IDdummy,KTYPE(J),XX(J),YY(J),ZZ(J),Q1X(J),Q1Y(J),Q1Z(J),CNA(J),CSP(J),StrX(J),StrY(J),StrZ(J)     !!!!! LOOP 1
        end do

!       ##############################	CHANGE ts to t (ps)  ###########################

	dt = VAR_TS	!	delta t of integration
	timed = timest*dt
	timenew=INT(timed)

	! Write the integer into a string:

  WRITE(file_name, '(I3)')  timenew

	! Open the file with this name
	open(unit = 50, file = '../Contour/'//trim(adjustl(file_name))//'_Contour.dat')
	open(unit = 60, file = '../Scatter/'//trim(adjustl(file_name))//'_Scatter.dat')
	open(unit = 70, file = '../Dimensions/'//trim(adjustl(file_name))//'_Particledimensions.dat')

!       ##############################	SIMULATION CELL SIZE AND CENTERS  ###########################

	XL=XHI-XLO
	YL=YHI-YLO
	ZL=ZHI-ZLO

	XCENTR=(XHI+XLO)/2.0d0
	YCENTR=(YHI+YLO)/2.0d0
	ZCENTR=(ZHI+ZLO)/2.0d0

!       ##############################	PARTICLE DIMENSIONS, WIDTHS AND HEIGHTS  ###########################
! 2 lattice units - for L8 that is ~60A square in xy plane through the centre to determine max and min z coordinate and thus height
! 2 lattice units - for L8 that is ~60A slice in x through the centre to determine max and min z coordinate and thus height

	ParticleHmax = ZCENTR
	ParticleHmin = ZCENTR

	ParticleWmax = YCENTR
	ParticleWmin = YCENTR


  DO I=1,NAN

    IF((KTYPE(I).EQ.1).AND.((CNA(I).NE.0).OR.(CSP(I).NE.0.0)).AND.(XX(I).LT.(XCENTR+30.0)).AND.(XX(I).GT.(XCENTR-30.0)).AND.(YY(I).LT.(YCENTR+30.0)).AND.(YY(I).LT.(YCENTR-30.0))) THEN      !!!!! LOOP 2
		    IF(ParticleHmax.LE.ZZ(I)) ParticleHmax = ZZ(I)
		    IF(ParticleHmin.GE.ZZ(I)) ParticleHmin = ZZ(I)
    END IF

	END DO

  WRITE(70,*) "Particle Height in Angstroms ="
  WRITE(70,*) (ParticleHmax-ParticleHmin)

	ZSMIN = ParticleHmin
	ZSMAX = ParticleHmax
  ZSL=ZSMAX-ZSMIN

  DO I=1,NAN                                                                                                                      !!!!! LOOP 3

	   IF((KTYPE(I).EQ.1).AND.((CNA(I).NE.0).OR.(CSP(I).NE.0.0)).AND.(XX(I).LT.(XCENTR+30.0)).AND.(XX(I).GT.(XCENTR-30.0)).AND.(ZZ(I).LT.(ParticleHmin+30.0))) THEN
		   IF(ParticleWmax.LE.YY(I)) ParticleWmax = YY(I)
		   IF(ParticleWmin.GE.YY(I)) ParticleWmin = YY(I)
     END IF

	END DO

  WRITE(70,*) "Particle Width in Angstroms ="
  WRITE(70,*) (ParticleWmax-ParticleWmin)


  XSMIN = XCENTR
  XSMAX = XSMIN
  YSMIN = YCENTR
  YSMAX = YSMIN
  ZSMIN = ZCENTR
  ZSMAX = ZSMIN

  DO I=1,NAN                                                      !!!!! LOOP 4
    IF((KTYPE(I).EQ.1).AND.((CNA(I).NE.0).OR.(CSP(I).NE.0.0))) THEN
      IF(XSMIN.GE.XX(I)) XSMIN = XX(I)
      IF(XSMAX.LE.XX(I)) XSMAX = XX(I)
      IF(YSMIN.GE.YY(I)) YSMIN = YY(I)
      IF(YSMAX.LE.YY(I)) YSMAX = YY(I)
      IF(ZSMIN.GE.ZZ(I)) ZSMIN = ZZ(I)
		  IF(ZSMAX.LE.ZZ(I)) ZSMAX = ZZ(I)
    END IF

 ENDDO

 print*,"XSMAX=", XSMAX
 print*,"XSMIN=", XSMIN
 print*,"YSMAX=", YSMAX
 print*,"YSMIN=", YSMIN
 print*,"ZSMAX=", ZSMAX
 print*,"ZSMIN=", ZSMIN


XSL=XSMAX-XSMIN
YSL=YSMAX-YSMIN
ZSL=ZSMAX-ZSMIN


!Changing to account for jetted atoms
YMID=YCENTR
XMID=XCENTR

!	Taking a 15 nm thick section in X direction,

XWID = VAR_WID

XLIM1=XMID-(XWID/2.0)
XLIM2=XMID+(XWID/2.0)


!	Number of bins in Y and Z
NZ = VAR_BINZ
NY = VAR_BINY

DZL = ZSL/NZ   ! Divide by the same number of cells as used earlier


XKE=0
YKE=0
ZKE=0
VZAVE(1:NVD,1:NVD)=0
VYAVE(1:NVD,1:NVD)=0
VXAVE(1:NVD,1:NVD)=0
VZSEC(1:NVD,1:NVD)=0
VYSEC(1:NVD,1:NVD)=0
VXSEC(1:NVD,1:NVD)=0
ZTEMP(1:NVD,1:NVD)=0
VMS(1:NVD,1:NVD)=0
ZVEL(1:NVD,1:NVD)=0
YVEL(1:NVD,1:NVD)=0
XVEL(1:NVD,1:NVD)=0
TKE(1:NAN)=0
TTKE(1:NVD,1:NVD)=0
YZMAX(1:NVD)=0
YZMIN(1:NVD)=YMID
YWID(1:NVD)=0

KAT(1:NVD,1:NVD)=0
ZPRESS(1:NVD,1:NVD)=0

STRESS_X(1:NVD,1:NVD)=0
STRESS_Y(1:NVD,1:NVD)=0
STRESS_Z(1:NVD,1:NVD)=0
STRESS_XY(1:NVD,1:NVD)=0
STRESS_YZ(1:NVD,1:NVD)=0
STRESS_ZX(1:NVD,1:NVD)=0


!!!!!!!!!!!!!!!!!! BINNING FOR THE SYSTEM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	CALCULATING YLENGTH AND YDIVISION OF EACH BIN FOR EACH Z (YLENGTH IS CALCULATED BY THE MAXIMUM Y LENGTH)
SECTION:Do KNZ=1,NZ
   CZMIN(KNZ) = ZSMIN+(KNZ-1)*DZL
   CZMAX(KNZ) = ZSMIN+KNZ*DZL

   YDIM:DO I=1,NAN                        !!!!! LOOP 5

      IF((ZZ(I).GE.(CZMIN(KNZ))).AND.(ZZ(I).LT.(CZMAX(KNZ))).AND.(XX(I).GE.XLIM1).AND.(XX(I).LT.XLIM2).AND.((CNA(I).NE.0).OR.(CSP(I).NE.0.0))) THEN

       	IF(YZMIN(KNZ).GE.YY(I)) YZMIN(KNZ) = YY(I)
       	IF(YZMAX(KNZ).LE.YY(I)) YZMAX(KNZ) = YY(I)

        GO TO 33
      END IF

33    CONTINUE
   END DO YDIM

!		Changing Ybin width as the the Ylength changes (i.e. moving from substrate region to particle region)
		YWID(KNZ) = YZMAX(KNZ)-YZMIN(KNZ)
		DYL(KNZ) = YWID(KNZ)/NY

END Do SECTION

    Print*,"  Binning over"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!! Calculate center of mass velocities in each bin for correction !!!!!!!!!!!!!!!!!!!

    SECTIONZ:Do KNZ=1,NZ
       CZMIN(KNZ) = ZSMIN+(KNZ-1)*DZL
       CZMAX(KNZ) = ZSMIN+KNZ*DZL

		   SECTIONY:Do KNY=1,NY
					CYMIN(KNZ,KNY) = YZMIN(KNZ)+(KNY-1)*DYL(KNZ)
					CYMAX(KNZ,KNY) = YZMIN(KNZ)+KNY*DYL(KNZ)

					NATOMS:DO I=1,NAN                !!!!! LOOP 6
						IF((ZZ(I).GE.(CZMIN(KNZ))).AND.(ZZ(I).LT.(CZMAX(KNZ))).AND.(YY(I).GE.(CYMIN(KNZ,KNY))).AND.(YY(I).LT.(CYMAX(KNZ,KNY))).AND.(XX(I).GE.XLIM1).AND.(XX(I).LT.XLIM2)) THEN
								KAT(KNZ,KNY)=KAT(KNZ,KNY)+1

								VZSEC(KNZ,KNY)=VZSEC(KNZ,KNY) + Q1Z(I)
								VYSEC(KNZ,KNY)=VYSEC(KNZ,KNY) + Q1Y(I)
								VXSEC(KNZ,KNY)=VXSEC(KNZ,KNY) + Q1X(I)

								GO TO 89
						END IF
89					CONTINUE
					END DO NATOMS

					VZAVE(KNZ,KNY)=VZSEC(KNZ,KNY)/KAT(KNZ,KNY)
					VYAVE(KNZ,KNY)=VYSEC(KNZ,KNY)/KAT(KNZ,KNY)
					VXAVE(KNZ,KNY)=VXSEC(KNZ,KNY)/KAT(KNZ,KNY)

		END Do SECTIONY
	END Do SECTIONZ

  Print*,"Average centre of mass Velocities calculated"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ATOM:DO I=1,NAN                       !!!!! LOOP 7

        Do KNZ=1,NZ
          CZMIN(KNZ) = ZSMIN+(KNZ-1)*DZL
          CZMAX(KNZ) = ZSMIN+KNZ*DZL

			    ZLength(KNZ)=CZMAX(KNZ)-CZMIN(KNZ)

				      Do KNY=1,NY
		              CYMIN(KNZ,KNY) = YZMIN(KNZ)+(KNY-1)*DYL(KNZ)
	                CYMAX(KNZ,KNY) = YZMIN(KNZ)+KNY*DYL(KNZ)

					        YLEN(KNZ,KNY) = CYMAX(KNZ,KNY)-CYMIN(KNZ,KNY)

					        !Calculating volume of each bin, the width in x dimension is 800nm
					        RVOL=XWID*YLEN(KNZ,KNY)*ZLength(KNZ)
						      IF((ZZ(I).GE.(CZMIN(KNZ))).AND.(ZZ(I).LT.(CZMAX(KNZ))).AND.(YY(I).GE.(CYMIN(KNZ,KNY))).AND.(YY(I).LT.(CYMAX(KNZ,KNY))).AND.(XX(I).GE.XLIM1).AND.(XX(I).LT.XLIM2)) THEN

                  !Removing bias from velocities to calculate temperature componenet of kinetic energies
							           VZ(I)=Q1Z(I)-VZAVE(KNZ,KNY)
							           VY(I)=Q1Y(I)-VYAVE(KNZ,KNY)
							           VX(I)=Q1X(I)-VXAVE(KNZ,KNY)

							           !MULTIPLYING BY 10^4 TO ACCOUNT FOR A/PS TO M/S
							           XKE = 0.5d0*mass*1E+04*VX(I)*VX(I)
							           YKE = 0.5d0*mass*1E+04*VY(I)*VY(I)
							           ZKE = 0.5d0*mass*1E+04*VZ(I)*VZ(I)

							           TKE(I) = (XKE+YKE+ZKE)
							           TTKE(KNZ,KNY)=TTKE(KNZ,KNY)+TKE(I)

                         !!!!!!!!! SWRITE LOOP THAT CALCULATES PRESSURE BEGINS

                         !	VIRIAL COMPONENT CORRECTION (REFER to https://en.wikipedia.org/wiki/Virial_stress)
                         !	Correcting for the term Tau_ij in terms of  m(vi - vi_avg)(vj-vj_avg)

                         ! Convert bar.A3 to Joules , multiply 10^-25


						            StrX(I) = StrX(I)*1D-25
          	            StrY(I) = StrY(I)*1D-25
          	            StrZ(I) = StrZ(I)*1D-25

                        !Velocities are multiplied by 100 to convert from A/ps to m/s
						            StrX(I) = StrX(I)+2*XMASS*Q1X(I)*VXAVE(KNZ,KNY)*1.0D04-XMASS*VXAVE(KNZ,KNY)*VXAVE(KNZ,KNY)*1.0D04
          	            StrY(I) = StrY(I)+2*XMASS*Q1Y(I)*VYAVE(KNZ,KNY)*1.0D04-XMASS*VYAVE(KNZ,KNY)*VYAVE(KNZ,KNY)*1.0D04
          	            StrZ(I) = StrZ(I)+2*XMASS*Q1Z(I)*VZAVE(KNZ,KNY)*1.0D04-XMASS*VZAVE(KNZ,KNY)*VZAVE(KNZ,KNY)*1.0D04

							          XP = StrX(I)*(KAT(KNZ,KNY)/RVOL)/(1.0d+04)
							          YP = StrY(I)*(KAT(KNZ,KNY)/RVOL)/(1.0d+04)
							          ZP = StrZ(I)*(KAT(KNZ,KNY)/RVOL)/(1.0d+04)

							          !CONVERT BACK TO REGULAR UNITS
							          XP = XP*1D+25
							          YP = YP*1D+25
							          ZP = ZP*1D+25

							          TPRESS(I) = -(XP+YP+ZP)/3

							          ZPRESS(KNZ,KNY)=ZPRESS(KNZ,KNY)+ TPRESS(I)

							          Stress_X(KNZ,KNY)=Stress_X(KNZ,KNY) + XP
							          Stress_Y(KNZ,KNY)=Stress_Y(KNZ,KNY) + YP
							          Stress_Z(KNZ,KNY)=Stress_Z(KNZ,KNY) + ZP

							          !CONVERTING THE ANGSTROMS/PS TO M/S
							          ZVEL(KNZ,KNY)=VZAVE(KNZ,KNY)*100.0d0
							          YVEL(KNZ,KNY)=VYAVE(KNZ,KNY)*100.0d0
							          XVEL(KNZ,KNY)=VXAVE(KNZ,KNY)*100.0d0

							          GO TO 34
						      END IF
				      END Do
        END Do

34      CONTINUE

  END DO ATOM


!	Writing output - in fort.50
  WRITE(50,*) "#KZ KY AvgZ AVGY Natoms TEMP Press SX SY SZ"

  Do KNZ = 1, NZ

		AvZ(KNZ) = (CZMAX(KNZ)+CZMIN(KNZ))/2.0d0

			Do KNY = 1, NY

				AvY(KNZ,KNY) = (CYMAX(KNZ,KNY)+CYMIN(KNZ,KNY))/2.0d0

				!CONVERT BACK TO GPa

				ZPRESS(KNZ,KNY)=ZPRESS(KNZ,KNY)/KAT(KNZ,KNY)

				Stress_X(KNZ,KNY)=Stress_X(KNZ,KNY)/KAT(KNZ,KNY)
				Stress_Y(KNZ,KNY)=Stress_Y(KNZ,KNY)/KAT(KNZ,KNY)
				Stress_Z(KNZ,KNY)=Stress_Z(KNZ,KNY)/KAT(KNZ,KNY)

				ZTEMP(KNZ,KNY)=(2.0d0*TTKE(KNZ,KNY))/(3000.0d0*KB*KAT(KNZ,KNY))

        !!!! Write only if the bins have atoms
				IF(KAT(KNZ,KNY).GT.1) THEN

			      WRITE(50,222) KNZ ,KNY, AvZ(KNZ), AvY(KNZ,KNY), KAT(KNZ,KNY),ZTEMP(KNZ,KNY),ZPRESS(KNZ,KNY),Stress_X(KNZ,KNY),Stress_Y(KNZ,KNY),Stress_Z(KNZ,KNY)

				END IF




			END Do
  END Do

  DO KNZ=1,NZ
	   DO KNY=1,NY
		     DO I=1,NAN                     !!!!! LOOP 8
			        IF((KAT(KNZ,KNY).GT.1).AND.(ZZ(I).GE.(CZMIN(KNZ))).AND.(ZZ(I).LT.(CZMAX(KNZ))).AND.(YY(I).GE.(CYMIN(KNZ,KNY))).AND.(YY(I).LT.(CYMAX(KNZ,KNY))).AND.(XX(I).GE.XLIM1).AND.(XX(I).LT.XLIM2)) THEN

				            WRITE(60,666) XX(I),YY(I),ZZ(I),KTYPE(I),KAT(KNZ,KNY),ZTEMP(KNZ,KNY),ZPRESS(KNZ,KNY),Stress_X(KNZ,KNY),Stress_Y(KNZ,KNY),Stress_Z(KNZ,KNY)

			        END IF
		     END DO
	   END DO
  END DO

	DO I=1,NAN                           !!!!! LOOP 9
		IF((KTYPE(I).EQ.2).AND.(XX(I).GE.XLIM1).AND.(XX(I).LT.XLIM2)) THEN
				WRITE(60,666) XX(I),YY(I),ZZ(I),KTYPE(I),0,0.0,0.0,0.0,0.0,0.0
		END IF
	END DO



222 FORMAT(1x,I4,1x,I4,1x,F12.5,1x,F12.5,1x,I7,1x,F12.5,1x,F12.5,x,F12.5,1x,F12.5,1x,F12.5)
666	FORMAT(1x,F12.5,1x,F12.5,1x,F12.5,1x,I4,1x,I7,1x,F12.5,1x,F12.5,1x,F12.5,1x,F12.5,1x,F12.5)

798 FORMAT(1x,F11.4,1x,F11.4,1x,F11.4,1x,I7)

END PROGRAM PT
