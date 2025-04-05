#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

if (( $# != 2 ))
then
	>&2 echo "usage: run_chainx_test.sh {MEM|MUM} length"
	exit 1
fi


chainx="../../chainX"
ref="chm13v2.0_concat.fa.gz"
query="sample_100k.fasta"
mode="-m sg"
usrbintime="/usr/bin/time -f"
options="%e total time elapsed (s)\n%M maxresident k"
anchors="-a $1 -l $2"

# run original ChainX
echo -n > original
$usrbintime "$options" $chainx $anchors $mode -t $ref -q $query \
	--originalmagicnumbers >> original 2>> original

# run optimal ChainX with non-constant initial guess
echo -n > optimal_nonconstant
$usrbintime "$options" $chainx $anchors $mode -t $ref -q $query \
	--optimal >> optimal_nonconstant 2>> optimal_nonconstant

echo -n "Checking (sanity) that the optimal chaining cost improves or stays the same..."
check=$(paste \
	<(grep "distance =" original            | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "distance =" optimal_nonconstant | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") | \
	awk '{if ($1 < $2) {print}}')
if [ "$check" != "" ]
then
	echo " it increases instead!"
	exit 1
fi
echo " done (it improves or stays the same)."

echo -n $(paste \
	<(grep "distance =" original            | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "distance =" optimal_nonconstant | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") | \
	awk '{if ($1 > $2) {improved += 1}} END {print improved}') "chaining costs improve"
echo " with average improvement ratio of" $(paste \
	<(grep "distance =" original            | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "distance =" optimal_nonconstant | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") | \
	awk '{if ($1 > $2) {ratiosum += $1 / $2 ; n += 1}} END {print ratiosum / n}')

# time, memory, iterations
for t in original optimal_nonconstant
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	echo "avg iters" > stats_iters_$t
done

for t in original optimal_nonconstant
do
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" $t | cut -d' ' -f1 >> stats_space_$t
	grep "iterations" $t | cut -d'=' -f2 | cut -d'(' -f2 | cut -d' ' -f1 | awk '{n += 1; sum += $1} END {print sum / n}' >> stats_iters_$t
done

paste -d'$' \
	stats_time_original             stats_space_original            stats_iters_original \
	stats_time_optimal_nonconstant  stats_space_optimal_nonconstant stats_iters_optimal_nonconstant \
	| cat <(echo -e "ChainX\$\$\$ChainX-opt*") - | column -t -s'$'

rm original optimal_nonconstant stats_time_original stats_space_original stats_iters_original stats_time_optimal_nonconstant stats_space_optimal_nonconstant stats_iters_optimal_nonconstant
