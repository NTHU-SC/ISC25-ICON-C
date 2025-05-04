#!/bin/bash

function check(){
    cdo -diffv $1 $2
    testok=$(cdo -diffv $1 $2 | sed -n '/records differ$/p' | tr -s " " | cut -f2 -d " ")
   
    if [[ $testok -eq 0 ]] ; then
        echo "[OK] Bit-identical results"
    else
        echo "[FAILED] Not bit-identical results. See diffs below:"
        cdo infon -sub $1 $2
    fi 
}
