duplicates
==========

Analysis of duplicate reads in Illumina sequencing using unique molecular identifiers (UMIs)

![alt text](https://github.com/mikessh/duplicates/blob/master/example/ex.png "Example")

### How to use

* Install [Java/Groovy](http://groovy.codehaus.org/) and R, call ```install.packages("ggplot2")``` in R to resolve dependency. Groovy dependency (Gpars) will be automatically resolved via [Grape](http://groovy.codehaus.org/Grape), additional configuration is only required if youre behind firewall/using proxy

* Clone the repository with ```git clone <repo_url>```

* Make sure your **fastq** files were prepared using [MiGEC/Checkout](https://github.com/mikessh/migec#1-checkout) utility with ```-u``` option that stores UMI tags in read headers

* Change dir to the folder containing **.fastq[.gz]** files you wish to analyze and run ```$bash path/to/duplicates/run.sh``` (only R1 will be considered)

* Enjoy fancy pdf output