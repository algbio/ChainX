#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

chainx="../../chainX"
ref="Chromosome_2890043_3890042_0.fasta"
mode="-m g"
usrbintime="/usr/bin/time -f"
options="%M maxresident k"

rm -f original optimal optimal_nonconstant

# run original ChainX
for m in 70 80 90 94 97 99
do
	$usrbintime "$options" $chainx $mode -t $ref -q "mutated_${m}_perc.fasta" \
		--originalmagicnumbers >> original 2>> original
done

# run optimal ChainX
for m in 70 80 90 94 97 99
do
	$usrbintime "$options" $chainx $mode -t $ref -q "mutated_${m}_perc.fasta" \
		--originalmagicnumbers --optimal >> optimal 2>> optimal
done

# run optimal ChainX with non-constant initial guess
for m in 70 80 90 94 97 99
do
	$usrbintime "$options" $chainx $mode -t $ref -q "mutated_${m}_perc.fasta" \
		--optimal >> optimal_nonconstant 2>> optimal_nonconstant
done

echo -n "Checking if the optimal chaining cost improves..."
for t in optimal optimal_nonconstant
do
	check=$(diff \
		<(grep "distance =" original | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
		<(grep "distance =" $t| cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
	if [ "$check" != "" ]
	then
		echo " it improves!"
		exit 1
	fi
done
echo " done (no improvement)."

# number of MUMs, should be the same as ChainX paper
echo "MUMs" > stats_mums
grep "count of anchors" optimal | cut -d' ' -f9 | tr -d "," | tac >> stats_mums

# time, memory, iterations
for t in original optimal optimal_nonconstant
do
	echo "time (s)" > stats_time_$t
	grep "distance computation finished" $t | cut -d'(' -f2 | cut -d' ' -f1 | tac >> stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "maxresident" $t | cut -d' ' -f1 | tac >> stats_space_$t
	echo "iters" > stats_iters_$t
	grep "iterations" $t | cut -d'=' -f2 | cut -d'(' -f2 | cut -d' ' -f1 | tac >> stats_iters_$t
done

paste -d'$' stats_mums \
	stats_time_original             stats_space_original            stats_iters_original \
	stats_time_optimal              stats_space_optimal             stats_iters_optimal \
	stats_time_optimal_nonconstant  stats_space_optimal_nonconstant stats_iters_optimal_nonconstant \
	| cat <(echo -e "\$ChainX\$\$\$ChainX-opt\$\$\$ChainX-opt*") - | column -t -s'$'

rm -f original optimal optimal_nonconstant stats_*
