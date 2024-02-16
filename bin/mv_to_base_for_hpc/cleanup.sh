#!/bin/bash

rm -rf pgt-tmp
rm -rf pst-tmp

find data/* -type f -name wait* 2>/dev/null | while read file
do
rm file
done

pgt_primers='data/primers/pgt_primers.csv'
pst_primers='data/primers/pst_primers.csv'

(head -n 1 -q data/primers/pgt/* | uniq && tail -n 1 -q data/primers/pgt/* | sort -k 2 -rn) > $pgt_primers
(head -n 1 -q data/primers/pst/* | uniq && tail -n 1 -q data/primers/pst/* | sort -k 2 -rn) > $pst_primers

rm -rf data/primers/pgt data/primers/pst
