#!/bin/bash
#SBATCH -p debug
#SBATCH -n 1

declare -a StringArray=("MD" "L2" "L4" "L6" "L8" )
declare -a pd=("100nm" )

##rm -rf MD L2 L4 L6 L8

for i in ${StringArray[@]}; do

	mkdir $i
	mkdir $i/Contour
	mkdir $i/Scatter
	mkdir $i/Dimensions

	cp Analysis_2D_PCorr.f90 job_2DA.sh $i/

	cd $i/

	if [ "$i" = "MD" ]; then

		for j in `seq 0 1000 50000`; do
		mkdir $j
		cp Analysis_2D_PCorr.f90 job_2DA.sh $j/

		cd $j/

		sed -i "s/VAR_ACG/1/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_PSIZE/${pd}/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_LOC/$i/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_TS/0.002/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_WID/150.0/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINZ/40/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINY/40/g" Analysis_2D_PCorr.f90

		sed -i "s/timestep/$j/g" Analysis_2D_PCorr.f90

		sbatch job_2DA.sh

		cd ../

		done


	elif [ "$i" = "L2" ]; then

		for j in `seq 0 500 25000`; do
		mkdir $j
		cp Analysis_2D_PCorr.f90 job_2DA.sh $j/
		cd $j/

		sed -i "s/VAR_ACG/2/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_PSIZE/${pd}/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_LOC/$i/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_TS/0.004/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_WID/150.0/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINZ/40/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINY/40/g" Analysis_2D_PCorr.f90

		sed -i "s/timestep/$j/g" Analysis_2D_PCorr.f90

		sbatch job_2DA.sh

		cd ../

		done

	elif [ "$i" = "L4" ]; then

		for j in `seq 0 250 12500`; do
		mkdir $j
		cp Analysis_2D_PCorr.f90 job_2DA.sh $j/
		cd $j/

		sed -i "s/VAR_ACG/4/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_PSIZE/${pd}/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_LOC/$i/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_TS/0.008/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_WID/150.0/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINZ/40/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINY/40/g" Analysis_2D_PCorr.f90

				sed -i "s/timestep/$j/g" Analysis_2D_PCorr.f90

		sbatch job_2DA.sh

		cd ../

		done

	elif [ "$i" = "L6" ]; then

		for j in `seq 0 167 8183`; do
		mkdir $j
		cp Analysis_2D_PCorr.f90 job_2DA.sh $j/
		cd $j/

		sed -i "s/VAR_ACG/6/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_PSIZE/${pd}/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_LOC/$i/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_TS/0.012/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_WID/150.0/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINZ/40/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINY/40/g" Analysis_2D_PCorr.f90

				sed -i "s/timestep/$j/g" Analysis_2D_PCorr.f90

		sbatch job_2DA.sh

		cd ../

		done

	elif [ "$i" = "L8" ]; then

		for j in `seq 0 125 6250`; do
		mkdir $j
		cp Analysis_2D_PCorr.f90 job_2DA.sh $j/
		cd $j/

		sed -i "s/VAR_ACG/8/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_PSIZE/${pd}/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_LOC/$i/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_TS/0.016/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_WID/150.0/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINZ/40/g" Analysis_2D_PCorr.f90
		sed -i "s/VAR_BINY/40/g" Analysis_2D_PCorr.f90

				sed -i "s/timestep/$j/g" Analysis_2D_PCorr.f90

		sbatch job_2DA.sh

		cd ../

		done

	fi

	cd ../

done
