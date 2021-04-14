#!/bin/bash
for dsp in *dsp
do
    dsp2=${dsp}2
    ps=$dsp.ps
    scale=`gmt gmtinfo -I1/0.1 $dsp`	
    gmt psxy $dsp $scale -Ba4f2:"period(s)":/a1f0.5:"v(km/s)"::."$dsp":WSne -JX4i/6i -W1p,red -Y3i -K -P >$ps
    gmt psxy $dsp2 $scale -Ba4f2:"period(s)":/a1f0.5:"v(km/s)":WSne -JX -O -W1p,black -P >>$ps
done
