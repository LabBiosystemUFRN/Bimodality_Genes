#!/bin/bash

#this script join two parts of the dataset and restore it to the original format

cd $1"/data/parts"

cat KIRC_FPKM_TP.tsv.tar.bz.part.* > KIRC_FPKM_TP.tsv.tar.bz

tar -xvf KIRC_FPKM_TP.tsv.tar.bz -C ../


