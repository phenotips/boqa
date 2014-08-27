#!/usr/bin/env bash

if [[ $1 == "-h" ]]; then
	echo "usage: $0 hpo_dir out_dir"
	exit
fi

hpo=$1
out=$2

cd /filer/tools/boqa/boqa-test/boqa-dist/boqa

ant clean

ant build

java -Xmx16G -cp bin:jars/commons-cli-1.2.jar sonumina.boqa.BOQABenchmark -o data/hp.obo.gz -a  /filer/tools/boqa/boqa-test/boqa-dist/boqa/data/new_phenotype.gz -p $hpo -d $out

