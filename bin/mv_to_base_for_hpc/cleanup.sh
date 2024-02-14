#!/bin/bash

rm -rf pgt_tmp
rm -rf pst_tmp

find data/* -type f -name wait* 2>/dev/null | while read file
do
rm file
done

pgt_primers='data/primers/pgt_primers.csv'
pst_primers='data/primers/pst_primers.csv'

head -n 1 -q data/primers/pgt/* | uniq > $pgt_primers
tail -n 1 -q data/primers/pgt/* >> $pgt_primers
sort -k 2 -n -o $pgt_primers $pgt_primers

head -n 1 -q data/primers/pst/* | uniq > $pst_primers
tail -n 1 -q data/primers/pst/* >> $pst_primers
sort -k 2 -n -o $pst_primers $pst_primers

rm -rf data/primers/pgt data/primers/pst

