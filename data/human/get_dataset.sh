#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

if [ ! -f "chm13v2.0.fa.gz" ]
then
	wget 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
	gunzip -c chm13v2.0.fa.gz | \
		awk '{if (substr($1, 0, 1) != ">") {print}}' | \
		tr "[:lower:]" "[:upper:]" | \
		cat <(echo ">chm13v2.0_concat") - | \
		gzip -c \
		> chm13v2.0_concat.fa.gz
fi

if [ ! -f "sample_100k.fa.gz" ]
then
	wget 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64004_210224_230828.hifi_reads.fastq.gz'
	seqtk sample -s 1 m64004_210224_230828.hifi_reads.fastq.gz 100000 | seqtk seq -A - | gzip -c > sample_100k.fa.gz
fi
