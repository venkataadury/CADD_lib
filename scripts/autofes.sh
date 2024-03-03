#!/bin/bash
SHIFT_BLOCKS=3
file=$1
if test -n "$2"; then
	SHIFT_BLOCKS=$2
fi

if test -n "$3"
then
	shift
	shift
fi

creations=""
for i in `seq 1 $SHIFT_BLOCKS`
do
	if test -n "$3"
	then
		python3 do_fes_mod.py $file 1 $i $SHIFT_BLOCKS $@ "fes_2cv_s"$i".dat"
	else
		python3 do_fes_mod.py $file 1 $i $SHIFT_BLOCKS 0.6 5.0 128 2.49 fes_2cv_s"$i".dat
	fi
	creations=$creations" fes_2cv_s"$i".dat "
done
str="awk '{print \$1\" \""
p=`expr $SHIFT_BLOCKS "*" 2 "+" 1`
echo $p
for i in `seq 2 2 $p`
do
	str=$str"\$"$i"\" \""
done
str=$str"}'"

extr="\$"2
end=`expr $SHIFT_BLOCKS "+" 1`
for i in `seq 3 $end`
do
	extr=$extr"\" \"\$"$i
done

echo $str
echo $extr

paste -d" " $creations |eval $str > fes_2cv_sall.dat
comm=`echo awk -e \''{printf($1" ");system("stats.mean "'$extr'"|tr \"\n\"  \" \" && stats.stddev "'$extr');}'\' fes_2cv_sall.dat`
eval $comm > fes_final_auto.dat
