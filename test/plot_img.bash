#!/bin/bash
par=loc_par
out=`awk 'NR==4{print $1}' $par`
range=`minmax -I0.1 $out`
echo $range
ps=img.ps
#pscoast $range -JM5i -Ba5f1/a1f1WSNE -K -P -W1p>img.ps
evlo=`awk 'NR==3{print $1}' $par`
evla=`awk 'NR==3{print $2}' $par`
#pscoast $range -JE$evlo/$evla/1/7i -Df -Ba1f1/a1f1 -W0.1p -Ggray -K -P -Y3i -X4i >$ps
range1=-R-54/-50/71/72.3
pscoast $range1 -Jl$evlo/$evla/70/73/5.0i -Df -Ba1f1/a1f1 -W0.1p -Ggray -K -P -Y3i -X3i >$ps
psclip -R  -J -O -K>>$ps<<eof
-53 71
-53 72
-51 72
-51 71
-53 71
eof
tscal=`minmax -T1/2 $out `
makecpt $tscal -Chaxby -I>a.cpt
surface $range -I0.001 $out -Ga.grd
grdimage a.grd $range1 -J -O -K -Ca.cpt >>$ps
psclip -C -O -K>>$ps
pscoast $range1 -J -B -Df -O -K -P -W1p>>img.ps
echo $evlo $evla | psxy -R -J -O -K -Sa0.2i -Gblack >>$ps
echo -53.1996  71.5384 | psxy -R -J -O -K -St0.2i -Gblack >>$ps
echo -53.1996 71.5384
min=`minmax -Ip -C $out | awk '{print $5}'`
echo $min
evlo=`grep $min $out | awk '{print $1}'`
evla=`grep $min $out | awk '{print $2}'`
t000=`grep $min $out | awk '{print $4}'`
echo $evlo $evla $t000

echo $evlo  $evla | psxy -R -J -O -Sa0.2i -Gred >>$ps
