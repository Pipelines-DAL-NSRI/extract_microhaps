#!/bin/bash

./plink/plink2 --vcf "$1" --chr "$2" --from-bp "$3" --to-bp "$4" --allow-extra-chr --make-pgen --out "$5"