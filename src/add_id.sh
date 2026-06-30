#!/bin/bash

./plink/plink2 --pfile "$1" --set-all-var-ids @:# --recode vcf --out "$2"
