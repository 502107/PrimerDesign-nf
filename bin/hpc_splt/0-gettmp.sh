#!/bin/bash
rust=$1
conserved_sites=$2
basedir=$3

mkdir -p $basedir/$rust-tmp
sed 's/>//g' $conserved_sites > $basedir/$rust-tmp/conserved.fna
tmpc=${basedir}/$rust-tmp/conserved.fna
# Reading two lines at a time
exec 3< $tmpc
while read -u 3 ID && read -u 3 SEQ
do
    echo ">$ID" > $basedir/$rust-tmp/$ID.fna
    echo "$SEQ" >> $basedir/$rust-tmp/$ID.fna
done
exec 3<&-
rm $tmpc
