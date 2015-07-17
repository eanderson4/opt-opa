#!/bin/bash
NAME=nocon
for mod in jcc oj
do
    echo "N c r rl0 rl T LS SD SE min max" > S$mod.out
    for i in {100..2900..100}
    do
	echo "${mod}S$NAME$i.out"
	sed '1d' data/${mod}S$NAME$i.out >> S$mod.out
    done

    rm -f D$mod.out
    touch D$mod.out
    for i in {100..2900..100}
    do
	echo "${mod}S$NAME$i.out"
	cat data/${mod}D$NAME$i.out >> D$mod.out
    done
done

#    echo "pow case/2383wp.db 1.03 1.5 1 .015 .95 .01 .5 .9 .005 .25 10 .9 .5 10 50 nocon$i $old $i"
#    ( nohup pow case/2383wp.db 1.03 1.5 1 .015 .95 .01 .5 .9 .005 .25 150 .9 .5 10 50 nocon$i $old $i > nocon$i.log & )
#    old=$i
