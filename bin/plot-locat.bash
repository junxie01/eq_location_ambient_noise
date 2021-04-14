#!/bin/bash
if [ $# -ne 1 ];then
	echo Usage: plot-imag.bash locat.par
	exit
fi
if [ $(uname) != 'Linux' ];then
	out=`awk '{print $2}' $1`
        min=`gmt gmtinfo -Ip -C $out | awk '{print $5}'`
        scale=`gmt gmtinfo -I0.5/0.5 $out`
        echo hello1
        gmt xyz2grd $scale $out -Gou.grd -I0.025
        tscale=`gmt gmtinfo -T1/2 $out`
        gmt makecpt -Cseis $tscale > ou.cpt
echo hello2
        gmt grdimage -Cou.cpt $scale ou.grd -JM5i -Ba0.5f0.25g0.5 -K -P -Y3i >img.ps
        gmt grdcontour $scale ou.grd -JM -Ba0.5f0.25g0.5 -W0.5 -A5 -O -K -P >>img.ps
        stlo=`grep $min $out | awk -v m=$min '$3==m{print $1}'`
        stla=`grep $min $out | awk -v m=$min '$3==m{print $2}'`
echo $stla $stlo
        gmt psxy -R -JM -B -O -K -Ss0.2i -Gblue >>img.ps<<eof
$stlo $stla
eof
        echo `distance -19.79 134.38 $stla $stlo`
        gmt psxy -R -JM -B -O -K -Sa0.2i -Gblack >>img.ps<<eof
134.38  -19.79  
eof
        gmt pstext -R -JM -B -O -K -W1p>>img.ps <<eof
134.38  -19.79 10 0 4 BL WB9 
eof
        gmt psxy -R -JM -B -O -K -Ss0.2i -Gblack >>img.ps<<eof
134.35  -19.96  
eof
        gmt pstext -R -JM -B -O -K -W1p>>img.ps <<eof
134.35  -19.96 10 0 4 BL WB1 
eof



else
        sta=`awk 'NR==1{print $1}' $1`
        ref=`awk 'NR==2{print $1}' $1`
        tgt=`awk 'NR==3{print $1}' $1`
        out=`awk 'NR==4{print $1}' $1`
        ps=$out.ps
        stlo_ref=`awk -v a=$ref '{if($1==a)print $2}' $sta`
        stla_ref=`awk -v a=$ref '{if($1==a)print $3}' $sta`
        stlo_tgt=`awk -v a=$tgt '{if($1==a)print $2}' $sta`
        stla_tgt=`awk -v a=$tgt '{if($1==a)print $3}' $sta`
        min=`minmax -Ip -C $out | awk '{print $5}'`
        scale=`minmax -I0.1/0.1 $out`
        xyz2grd $scale $out -Gou.grd -I0.01
        tscale=`minmax -T5000/2 $out`
        makecpt -Cseis $tscale > ou.cpt
        grdimage -Cou.cpt $scale ou.grd -JM4i -Ba0.5f0.25g0.5 -K -P -Y3i >$ps
        grdcontour $scale ou.grd -JM4i -A50000 -Ba0.5f0.25g0.5 -W0.5 -O -K -P >>$ps
        stlo=`grep $min $out | awk -v m=$min '$3==m{print $1}'`
        stla=`grep $min $out | awk -v m=$min '$3==m{print $2}'`
        echo $stlo $stla | psxy -R -JM -B -O -K -Sa0.2i -Gblue >>$ps
        dist1=`distance $stla $stlo $stla_tgt $stlo_tgt`
        dist0=`distance $stla_ref $stlo_ref $stla_tgt $stlo_tgt`
        echo $dist0 $dist1
        dist0=`printf "%5.2f\n" $dist0`
        dist1=`printf "%5.2f\n" $dist1`
        echo $stlo_ref $stla_ref | psxy -R -JM -B -O -K -St0.2i -Gblack >>$ps
        echo $stlo_tgt $stla_tgt | psxy -R -JM -B -O -K -Sa0.2i -Gblack >>$ps
        echo $stlo_ref $stla_ref 15 0 4 BL $ref | pstext -R -JM -B -O -K >>$ps 
        echo $stlo_tgt $stla_tgt 15 0 4 BL $tgt | pstext -R -JM -B -O -K >>$ps 
        echo 134.7 -19.42 10 0 4 BL REF: ${dist0}km | pstext -R -J -B -O -K  >>$ps
        echo 134.7 -19.45 10 0 4 BL ERR: ${dist1}km | pstext -R -J -B -O >>$ps
fi
