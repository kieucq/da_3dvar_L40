#!/bin/sh
echo "enter the time cycle to plot"
read timing
ln -sf ../fsc/fsc01_$timing.dat ./baro01.dat
ln -sf ../fsc/fsc02_$timing.dat ./baro02.dat
ln -sf ../fsc/fsc03_$timing.dat ./baro03.dat
ln -sf ../fsc/fsc04_$timing.dat ./baro04.dat
ln -sf ../fsc/fsc05_$timing.dat ./baro05.dat
ln -sf ../fsc/fsc06_$timing.dat ./baro06.dat
ln -sf ../fsc/fsc07_$timing.dat ./baro07.dat
ln -sf ../fsc/fsc08_$timing.dat ./baro08.dat
ln -sf ../fsc/fsc09_$timing.dat ./baro09.dat
ln -sf ../fsc/fsc10_$timing.dat ./baro10.dat
ln -sf ../fsc/fsc11_$timing.dat ./baro11.dat
