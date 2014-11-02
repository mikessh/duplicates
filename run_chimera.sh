ME=`dirname $0`
PROJECT_NAME="project"

export JAVA_OPTS="$JAVA_OPTS -Xmx64G"

groovy $ME/groovy/DuplicatesChimera.groovy `ls *R2*.fastq*` $PROJECT_NAME
Rscript $ME/R/plot_duplicates_chimera.R $PROJECT_NAME