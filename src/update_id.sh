#!/bin/bash

./plink/plink2 --vcf "$1" --update-name "$2" --recode vcf --out "$3"