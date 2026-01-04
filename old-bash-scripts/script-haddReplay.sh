
nFiles=$(ls -l replay/partfiles/replay.$1.* | wc -l) 


numMB=$((nFiles * 1100))
numMB=$(($numMB + 256))

netInputBytes=$(ls -l replay/partfiles/replay.$1.* | awk '{sum += $5} END{ print sum}')

sizeMB=1000000

netInputMB=$(( $netInputBytes / $sizeMB   +  512 ))

echo "size of files to hadd = $netInputMB MB"

sbatch --job-name="hadd_run_$1" --mem-per-cpu=$netInputMB script-haddReplay $1

echo "number of flies = $nFiles  memory requested=$netInputMB MB"
