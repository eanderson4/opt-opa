#!/bin/bash

rm -f rts96.db

#bus.txt
sqlite3 rts96.db "create table bus (num INTEGER, name VARCHAR(20), type INTEGER, pd DOUBLE, qd DOUBLE, gs DOUBLE, bs DOUBLE, area INTEGER, basekv DOUBLE, zone INTEGER);"
cat bus.txt | head -n -9 | sed -e '1,6d' | sed 's/ \+/,/g' > tmpbus
sqlite3 -separator ',' rts96.db ".import tmpbus bus"

sqlite3 rts96.db "select * from bus"

#branch.txt
sqlite3 rts96.db "create table branch (name VARCHAR(20), fbus INTEGER, tbus INTEGER, l DOUBLE, permlam DOUBLE, duration DOUBLE, translam DOUBLE, br_r DOUBLE, br_x DOUBLE, br_b DOUBLE, rate_a DOUBLE, rate_b DOUBLE, rate_C DOUBLE, tap DOUBLE);"
cat branch.txt | head -n -2 | sed -e '1,17d' | sed 's/ \+/,/g' > tmpbranch
sqlite3 -separator ',' rts96.db ".import tmpbranch branch"

sqlite3 rts96.db "select * from branch"

#busload.txt
sqlite3 rts96.db "create table busload (bus1 INTEGER, bus2 INTEGER, bus3 INTEGER, percentload DOUBLE, pd DOUBLE, qd DOUBLE, peakpd DOUBLE, peakqd DOUBLE);"
cat busload.txt | head -n -2 | sed -e '1,5d' | sed 's/ \+/,/g' > tmpload
sqlite3 -separator ',' rts96.db ".import tmpload busload"

sqlite3 rts96.db "select * from busload"


#daily.txt
sqlite3 rts96.db "create table daily (daynum INTEGER, day VARCHAR(20), load DOUBLE)"
cat daily.txt | head -n -0 | sed -e '1,4d' | nl -n ln -s ' '  | sed 's/ \+/,/g' >tmpdaily
sqlite3 -separator ',' rts96.db ".import tmpdaily daily"

sqlite3 rts96.db "select * from daily"

#gencost
sqlite3 rts96.db "create table gencost (gid VARCHAR(20), size DOUBLE, type VARCHAR(20), ramp DOUBLE, c2 DOUBLE, c1 DOUBLE);"
cat gencost.txt | head -n -0 | sed -e '1,3d' | sed 's/ \+/,/g' | sed 's/\t//g' >tmpcost
sqlite3 -separator ',' rts96.db ".import tmpcost gencost"

sqlite3 rts96.db "select * from gencost"


#gendynamic
sqlite3 rts96.db "create table gendynamic (gid VARCHAR(20), size DOUBLE, type VARCHAR(20), mva DOUBLE, react1 DOUBLE, react2 DOUBLE, inertia DOUBLE, damp DOUBLE);"
cat gendynamic.txt | head -n -6 | sed -e '1,6d' | sed 's/ \+/,/g' > tmpdyn
sqlite3 -separator ',' rts96.db ".import tmpdyn gendynamic"

sqlite3 rts96.db "select * from gendynamic"

#genparm.txt
sqlite3 rts96.db "create table genparm (bus INTEGER, gid VARCHAR(20), id INTEGER, pg DOUBLE, qg DOUBLE, qmax DOUBLE, qmin DOUBLE, vs DOUBLE);"
cat genparm.txt | head -n -4 | sed -e '1,5d' | sed 's/ \+/,/g' >tmpparm
sqlite3 -separator ',' rts96.db ".import tmpparm genparm"

sqlite3 rts96.db "select * from genparm"

#genramp.txt
sqlite3 rts96.db "create table genramp (gid VARCHAR(20), size DOUBLE, name VARCHAR(20), mindown DOUBLE, minup DOUBLE, starthot DOUBLE, startcold DOUBLE, warmstart DOUBLE, ramp DOUBLE);"
cat genramp.txt | head -n -0 | sed -e '1,5d' | sed 's/ \+/,/g' >tmpramp
sqlite3 -separator ',' rts96.db ".import tmpramp genramp"

sqlite3 rts96.db "select * from genramp"

#genrel.txt
sqlite3 rts96.db "create table genrel (gid VARCHAR(20), size DOUBLE, name VARCHAR(20), outage DOUBLE, mttf DOUBLE, mttr DOUBLE, maint DOUBLE);"
cat genrel.txt | head -n -0 | sed -e '1,6d' | sed 's/ \+/,/g' >tmprel
sqlite3 -separator ',' rts96.db ".import tmprel genrel"

sqlite3 rts96.db "select * from genrel"


#genstart.txt
sqlite3 rts96.db "create table genstart (gid VARCHAR(20), size DOUBLE, name VARCHAR(20), hot DOUBLE, cold DOUBLE);"
cat genstart.txt | head -n -0 | sed -e '1,6d' | sed 's/ \+/,/g' >tmpstart
sqlite3 -separator ',' rts96.db ".import tmpstart genstart"

sqlite3 rts96.db "select * from genstart"

#hourly.txt
sqlite3 rts96.db "create table hourly (time INTEGER, tname VARCHAR(20), winterwkdy DOUBLE, winterwknd DOUBLE, summerwkdy DOUBLE, summerwknd DOUBLE, springwkdy DOUBLE, springwknd DOUBLE);"
cat hourly.txt | head -n -1 | sed -e '1,6d' | nl -n ln -s ' ' | sed 's/ \+/,/g' >tmphourly
sqlite3 -separator ',' rts96.db ".import tmphourly hourly"

sqlite3 rts96.db "select * from hourly"

#weekly.txt
sqlite3 rts96.db "create table weekly (week INTEGER, peak DOUBLE, season INTEGER);"
cat weekly.txt | head -n -0 | sed -e '1,4d' | sed 's/ \+/,/g' >tmpweekly
sqlite3 -separator ',' rts96.db ".import tmpweekly weekly"

sqlite3 rts96.db "select * from weekly"

rm tmp*
