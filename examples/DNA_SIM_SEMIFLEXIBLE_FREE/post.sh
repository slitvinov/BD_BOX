#! /bin/bash

for dcd in *.dcd; do
     echo $dcd
done | parallel octave -qf ~/google-svn/octave/matdcd-scripts/post-new.m
