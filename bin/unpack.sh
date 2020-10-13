#!/bin/bash

#this script join two parts of the dataset and restore it to the original format

cd $1"/data/"

cat KIRC_FPKM_TP.tsv.tar.bz.part.aa KIRC_FPKM_TP.tsv.tar.bz.part.ab > KIRC_FPKM_TP.tsv.tar.bz

tar -xjvf KIRC_FPKM_TP.tsv.tar.bz -C ./


