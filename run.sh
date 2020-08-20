#!/bin/bash

timer=$SECONDS

nextflow run main.nf -resume

timer=$(($SECONDS-timer))
printf "Time used: %02d:%02d:%02d\n" "$((timer/3600))" "$((timer/60%60))" "$((timer%60))"

