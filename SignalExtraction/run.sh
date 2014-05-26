#!/bin/bash

#root -l -b -q newBoxCuts.C+\(\"merged.conf\",0\) > all.out
root -l -b -q newBoxCuts.C+\(\"new.conf\",1\) > dijet45.out # dijet
root -l -b -q newBoxCuts.C+\(\"new.conf\",2\) > jetmu.out # jetmu
root -l -b -q newBoxCuts.C+\(\"new.conf\",3\) > jetele.out # jetele
root -l -b -q newBoxCuts.C+\(\"new.conf\",4\) > muele.out # muele
root -l -b -q newBoxCuts.C+\(\"new.conf\",5\) > dijet30.out # dijet