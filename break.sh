#!/bin/bash
old=0
for i in {150..3000..150}
do
    echo "$i"
#    echo "pow case/2383wp.db 1.03 1.5 1 .015 .95 .01 .5 .9 .005 .25 10 .9 .5 10 50 nocon$i $old $i"
#    ( nohup pow case/2383wp.db 1.03 1.5 1 .015 .95 .01 .5 .9 .005 .25 150 .9 .5 10 50 newdes$i $old $i > newdes$i.log & )

    echo "pow case/2383wp.db 1.03 1.5 1 .0013 .5 .2 .5 .99 .005 .5 10 .9 .5 10 50 olddes$i $old $i"
    ( nohup pow case/2383wp.db 1.03 1.5 1 .0013 .5 .2 .5 .99 .005 .5 150 .9 .5 .5 18.5 old18$i old $old $i > old18$i.log & )

    old=$i
done
