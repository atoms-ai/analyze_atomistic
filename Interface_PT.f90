!     Calculation of Pressure/Temperature Distribution at the interface
!     sumit.suresh@uconn.edu

PROGRAM PT

        IMPLICIT REAL*8(A-H,O-Z)
        INTEGER, PARAMETER:: KREAL=SELECTED_REAL_KIND(14,99) !14 sig figs, 10+-99


        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: XX,YY,ZZ,YC          !Coordinates
        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: StrX,StrY,StrZ,SX,SY,SZ,Vol,PRESS   !Stress
        REAL(KREAL), DIMENSION(:), ALLOCATABLE :: Q1X,Q1Y,Q1Z,VX,VY,VZ,VX_avg,VY_avg,VZ_avg,V_avg,CSP,KE,TEMPB  !Velocities


        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE !Type & History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: CNA, BINY, BIN_COUNT


        INTEGER:: NAN, count, NBIN, LOC, NCG                               !Number of particles
        REAL(KREAL):: TIME,mass,KB, LATT, CSP_MAX                          !Starting time
        REAL(KREAL):: XL,YL,ZL,DZL, height, thicc, bwidth, bin_width, intf_width        !Size of the computational celL
        REAL(KREAL):: XCENTR,YCENTR,ZCENTR        !CENTRE of the computational cel

	      CHARACTER *15 :: file_name

!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!	Use this to run
! ifort -shared-intel -mcmodel=large  Analysis2DBIN.f90


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEVEL OF COARSENING !!!!!!!!!!!!!!!!!!!!!!!!!!

!	LOC = VAR_ACG
  LOC = VAR_ACG
	NCG = LOC**3

  LATT = 4.05d0*LOC
  CSP_MAX = 3.0*(LATT**2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINE CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
AMUTOKG = 1.6605402D-27
bar2gpa = 1.0D-4
EVTOJOU = 1.60219D-19
XJOUTOEV = 1.d0/EVTOJOU

mass_amu = 26.98        ! Species = aluminum
mass_amu = mass_amu * NCG
r = 1.43          ! particle radius in angstroms
at_vol = (4/3)*3.14156*(r)**3

mass = mass_amu*AMUTOKG

amu_by_kb = 1.2027D-4       ! amu to kg divided by boltzmann constant in SI units: for temperature calculation
aps_2_mps = 100.0           ! Angstroms/ps to m/s

XMASS=4.480D-26 !in kg of Al - scale this with level of coarsening to correct for pressures
XMASS=XMASS*NCG

mass=4.480d0 ! For calculating temperatures from KE
KB=1.380d0  ! Boltzmann constant

ENUNIT = AMUTOKG*1.d4*XJOUTOEV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Filename
	!OPEN (UNIT = 14,FILE='/gpfs/scratchfs1/sumit/Coldspray/Rigid_impact/VAR_PSIZE/Impact/VAR_LOC/dump.imp1000mps.timestep')
OPEN (UNIT = 14,FILE='/gpfs/scratchfs1/sumit/Coldspray/Rigid_impact/VAR_PSIZE/Impact/VAR_LOC/Jetting_fine/dump.imp1000mps.timestep')

        REWIND 14
        READ(14,*)
	      READ(14,*) TIMEST
	      READ(14,*)
	      READ(14,*) NAN
	      READ(14,*)

ALLOCATE (XX(NAN),YY(NAN),ZZ(NAN))
ALLOCATE (StrX(NAN),StrY(NAN),StrZ(NAN))
ALLOCATE (Q1X(NAN),Q1Y(NAN),Q1Z(NAN),CSP(NAN))
ALLOCATE (KTYPE(NAN))
ALLOCATE (CNA(NAN))

        READ(14,*) XLO,XHI
        READ(14,*) YLO,YHI
        READ(14,*) ZLO,ZHI
        READ(14,*)


        XL=XHI-XLO
      	YL=YHI-YLO
      	ZL=ZHI-ZLO

      	XCENTR=(XHI+XLO)/2.0d0
      	YCENTR=(YHI+YLO)/2.0d0
      	ZCENTR=(ZHI+ZLO)/2.0d0

        SubstrateHmax=0.0

        ParticleWmax = 0.0
        Ymax = 0.0
        ParticleWmin = YCENTR

        !       ##############################	CARVE OUT THE INTERFACE  ###########################
        !	(20A cube) and ~1 lattice unit (5A) in height for MD, scale as needed

          bwidth=20.0
          height=20.0
          h=5.0*LOC       ! height for accurate prediction of Y max (keep it small ~ 5 angstrom)
          thicc=20.0

          print*, "x dimension thickness in Angstroms = ", thicc
          print*, "z dimension height in Angstroms = ", height

        DO J=1,NAN

         READ(14,*) IDdummy,KTYPE(J),XX(J),YY(J),ZZ(J),Q1X(J),Q1Y(J),Q1Z(J),CNA(J),CSP(J),StrX(J),StrY(J),StrZ(J)     !!!!! LOOP to assign jetted atoms and find max sub height

          IF ((CNA(J).EQ.5).AND.((CSP(J).EQ.0.0).OR.(CSP(J).GT.CSP_MAX))) THEN     !Assign different CNA values to jetted atoms
           CNA(J)=10
          END IF

          IF((KTYPE(J).EQ.2).AND.(CNA(J).NE.10)) THEN      !!!!! Substrate maximum height
              IF(SubstrateHmax.LE.ZZ(J)) SubstrateHmax = ZZ(J)
          END IF

        END DO

        pxlim1=xcentr-(thicc/2.0)
        pxlim2=xcentr+(thicc/2.0)

        pzlim1=SubstrateHmax
        pzlim2=SubstrateHmax+height

        pzlim_int=SubstrateHmax + h

        count=0   !counter for number of atoms in the interface for binning

        INTERF:DO J=1,NAN                                                                                                                      !!!!! LOOP to find width of interface and atom counting

           IF((KTYPE(J).EQ.1).AND.(XX(J).LT.pxlim2).AND.(XX(J).GT.pxlim1).AND.(ZZ(J).LT.pzlim2).AND.(CNA(J).NE.10)) THEN
             IF(ParticleWmax.LE.YY(J)) ParticleWmax = YY(J)
             IF(ParticleWmin.GE.YY(J)) ParticleWmin = YY(J)

             count=count+1

             XX(count)=XX(J)
             YY(count)=YY(J)
             ZZ(count)=ZZ(J)

             Q1X(count)=Q1X(J)
             Q1Y(count)=Q1Y(J)
             Q1Z(count)=Q1Z(J)

             StrX(count)=StrX(J)
             StrY(count)=StrY(J)
             StrZ(count)=StrZ(J)

           END IF

           IF((KTYPE(J).EQ.1).AND.(XX(J).LT.pxlim2).AND.(XX(J).GT.pxlim1).AND.(ZZ(J).LT.pzlim_int).AND.(CNA(J).NE.10)) THEN
             IF(Ymax.LE.YY(J)) Ymax = YY(J)
           END IF

        END DO INTERF

        intf_width = (ParticleWmax-ParticleWmin)
        print *, "Y max in Angstroms =", Ymax
        print *, "Interface Width in Angstroms =", intf_width
        print *, "Number of atoms in the binning interface =", count

        ALLOCATE (BINY(count))


!       ##############################	Binning  ###########################

        NBIN = NINT(intf_width/bwidth)     !Rounding to the nearest whole number
        bin_width=(intf_width/NBIN)

        print *, "Bin Width in Y dimension in Angstroms =", bin_width
        print *, "No of bins =", NBIN

        !Minimum bin value should be 1

        IF (NBIN.LT.1) THEN
          NBIN=1
        END IF

        ALLOCATE (bin_count(NBIN),YC(NBIN))
        ALLOCATE (SX(NBIN),SY(NBIN),SZ(NBIN),PRESS(NBIN), VOL(NBIN))
        ALLOCATE (VX(NBIN),VY(NBIN),VZ(NBIN),VX_avg(NBIN),VY_avg(NBIN),VZ_avg(NBIN),V_avg(NBIN),KE(NBIN),TEMPB(NBIN))


        bin_count = 0
        YC = 0.0
        SX = 0.0
        SY = 0.0
        SZ = 0.0
        VOL = 0.0
        PRESS = 0.0
        VX = 0.0
        VX_avg = 0.0
        VY = 0.0
        VY_avg = 0.0
        VZ = 0.0
        VZ_avg = 0.0
        KE = 0.0
        TEMPB = 0.0

        !	Assign atoms to their specific bin

        COUNTING:DO J=1,count                        !!!!! Assign bin atoms and averages (J for atom count, K for bin count, I for binatom count)

          DO K=1,NBIN

            YMIN = ParticleWmin
            YMAX = ParticleWmax

            IF((YY(J).GE.(YMIN+(K-1)*bin_width)).AND.(YY(J).LT.(YMAX-(NBIN-K)*bin_width))) THEN

              BINY(J) = K                 !Bin assigned
              bin_count(K)=bin_count(K)+1 !Bin counter updated
              YC(K)= YMIN + (2*K-1)*(bin_width/2.0)                      !Bin y center coordinate

              VZ(K)=VZ(K) + Q1Z(J)
              VY(K)=VY(K) + Q1Y(J)
              VX(K)=VX(K) + Q1X(J)

              SX(K)=SX(K) + StrX(J)
              SY(K)=SY(K) + StrY(J)
              SZ(K)=SZ(K) + StrZ(J)

            END IF

          END DO

        END DO COUNTING


        pressure:DO K=1,NBIN

            VZ_avg(K)=VZ(K)/bin_count(K)
            VY_avg(K)=VY(K)/bin_count(K)
            VX_avg(K)=VX(K)/bin_count(K)

            VOL(K) = height*thicc*bin_width     !Initial volume

            !Volume of side bins
            IF (NBIN.EQ.1) then
              VOL(K) = bin_count(K)*at_vol
            ELSE
              VOL(1) = bin_count(1)*at_vol
              VOL(NBIN) = bin_count(NBIN)*at_vol
            END IF

            !Calculate bin stresses and convert from bars to gpa
            SX(K)=SX(K)/VOL(K) * bar2gpa
            SY(K)=SY(K)/VOL(K) * bar2gpa
            SZ(K)=SZ(K)/VOL(K) * bar2gpa

            !Calculate bin pressure
            PRESS(K)=(-0.3333333)*(SX(K)+SY(K)+SZ(K))

        END DO pressure

        print*, " Pressure calculation over"


        kin_en:DO J=1,count

            DO K=1,NBIN

              IF(BINY(J).EQ.K) then         ! If atom belongs to a specific bin

                Q1Z(J) = Q1Z(J) - VZ_avg(K)
                Q1Y(J) = Q1Y(J) - VY_avg(K)
                Q1X(J) = Q1X(J) - VX_avg(K)

                !KE in amu units
                XKE = 0.5d0*mass_amu*Q1X(J)*Q1X(J)*(aps_2_mps)**2
                YKE = 0.5d0*mass_amu*Q1Y(J)*Q1Y(J)*(aps_2_mps)**2
                ZKE = 0.5d0*mass_amu*Q1Z(J)*Q1Z(J)*(aps_2_mps)**2

                KE(K) = KE(K) + (XKE+YKE+ZKE)

                !TEMPB(K) = TEMPB(K) + 2.0d0*KE(K)*amu_by_kb/(3.0d0*bin_count(K))

              END IF

            END DO

        END DO kin_en

        temperature:DO K=1,NBIN

          TEMPB(K) = 2.0d0*KE(K)*amu_by_kb/(3.0d0*bin_count(K))

        END DO temperature

        Print*, " Temperature calculation over"



!       ##############################	CHANGE ts to t (ps)  ###########################

!	dt = VAR_TS	!	delta t of integration
  dt = VAR_TS
	timed = timest*dt
	timenew=INT(timed)

	! Write the integer into a string:

  WRITE(file_name, '(I0)')  timenew

	! Open the file with this name
	!open(unit = 50, file = trim(adjustl(file_name))//'fs_PT_v_Y.dat')
  open(unit = 50, file = '../../IntP/PT_v_Y_realtime_fs.dat')

  DO K=1, NBIN

    WRITE(50,798) K, YC(K), bin_count(K), PRESS(K), TEMPB(K), VOL(K), Vy_avg(K)*aps_2_mps

  END DO

798 FORMAT(I0,1x,F0.4,1x,I0,1x,F0.4,1x,F0.4,1x,F0.4,1x,F0.4)

END PROGRAM PT
