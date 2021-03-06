rm -r output
mkdir output
cd output
mkdir datafile distribution
cd ..
rm nohup.out
rm -r ifortfile
mkdir ifortfile
ifort -O3 -ipo -parallel -traceback -CB -gen_interfaces -warn all -o HZT.exe \
program/parameter.f90 \
program/global.f90 \
program/mersenne.f90 \
program/initialcondition.f90 \
program/lhs.f90 \
program/electron_fluid.f90 \
program/xecollision.f90 \
program/particle.f90 \
program/main.f90
mv *.mod ifortfile/
mv *__genmod.f90 ifortfile/
./HZT.exe

