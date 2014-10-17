duplicates
==========

Analysis of duplicate reads in Illumina sequencing using unique molecular identifiers (UMIs)

![alt text](https://github.com/mikessh/duplicates/blob/master/example/ex.png "Example")

How to use

* Install Java/Groovy and R, in R call ```install.packages("ggplot2")``` to resolve a dependency, Groovy dependency (Gpars) will be automatically resolved via Grapes

* Clone the repository with ```git clone <repo_url>```

* Make sure your FASTQ files were prepared using [MiGEC/Checkout](https://github.com/mikessh/migec#1-checkout) utility with ```-u``` option that stores UMI tags in read headers

* From the folder containing *.fastq[.gz]* files run ```$bash path/to/duplicates/run.sh```

* Enjoy fancy pdf output