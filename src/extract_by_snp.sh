#!/bin/bash

./plink/plink2 --vcf "$1" --snps "$2" --recode vcf --out "$3"
