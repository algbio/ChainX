#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

chainx="../../chainX"
ref="Chromosome_2890043_3890042_0.fasta"
mode="-m g"
usrbintime="/usr/bin/time -f"
options="%e total time elapsed (s)\n%M maxresident k"

# run original ChainX
for m in "75_80" "80_90" "90_100" 
do
	echo -n > original_$m
	$usrbintime "$options" $chainx $mode -t $ref -q "count100_mutated_$m.fa.gz" \
		--originalmagicnumbers >> original_$m 2>> original_$m
done

# run optimal ChainX
for m in "75_80" "80_90" "90_100" 
do
	echo -n > optimal_$m
	$usrbintime "$options" $chainx $mode -t $ref -q "count100_mutated_$m.fa.gz" \
		--originalmagicnumbers --optimal >> optimal_$m 2>> optimal_$m
done

# run optimal ChainX with non-constant initial guess
for m in "75_80" "80_90" "90_100" 
do
	echo -n > optimal_nonconstant_$m
	$usrbintime "$options" $chainx $mode -t $ref -q "count100_mutated_$m.fa.gz" \
		--optimal >> optimal_nonconstant_$m 2>> optimal_nonconstant_$m
done

echo -n "Checking if the optimal chaining cost improves..."
for t in optimal optimal_nonconstant
do
	for m in "75_80" "80_90" "90_100" 
	do
		check=$(diff \
			<(grep "distance =" original_$m | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
			<(grep "distance =" ${t}_${m} | cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
		if [ "$check" != "" ]
		then
			echo " it improves!"
			exit 1
		fi
	done
done
echo " done (no improvement)."

# time, memory, iterations
echo "similarity" > stats_headers

for t in original optimal optimal_nonconstant
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	echo "avg iters" > stats_iters_$t
done

for m in "75_80" "80_90" "90_100" 
do
	echo "$m" >> stats_headers
	for t in original optimal optimal_nonconstant
	do
		grep "total time elapsed" ${t}_${m} | cut -d' ' -f1 >> stats_time_$t
		grep "maxresident" ${t}_${m} | cut -d' ' -f1 >> stats_space_$t
		grep "iterations" ${t}_${m} | cut -d'=' -f2 | cut -d'(' -f2 | cut -d' ' -f1 | awk '{n += 1; sum += $1} END {print sum / n}' >> stats_iters_$t
	done
done

paste -d'$' stats_headers \
	stats_time_original             stats_space_original            stats_iters_original \
	stats_time_optimal              stats_space_optimal             stats_iters_optimal \
	stats_time_optimal_nonconstant  stats_space_optimal_nonconstant stats_iters_optimal_nonconstant \
	| cat <(echo -e "\$ChainX\$\$\$ChainX-opt\$\$\$ChainX-opt*") - | column -t -s'$'

rm optimal_* original_* stats_*
