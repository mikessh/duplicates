#!/bin/bash

ME=`dirname $0`
#SAME_BARCODE_SPACE="1"
PROJECT_NAME="project"

export JAVA_OPTS="$JAVA_OPTS -Xmx64G"

groovy $ME/groovy/Duplicates.groovy `ls *R1*fastq*` $PROJECT_NAME
Rscript $ME/R/plot_duplicates.R $PROJECT_NAME