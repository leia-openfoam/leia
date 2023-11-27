#! /usr/bin/env bash

mkdir -p plots
cp $(realpath $1)/study_*.json plots/
study_cases.py $1 | xargs -I{} bash -c 'case={} && cp -v $case/*.png plots/'
