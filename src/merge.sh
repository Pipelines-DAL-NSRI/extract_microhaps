#!/bin/bash

bcftools concat -f "$1" -Oz > "$2"