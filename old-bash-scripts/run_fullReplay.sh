ls  /cache/mss/halla/apex/raw/apex_$1.dat.* > listRaw.$1.lis

$num_part=0

while IFS= read -r file; do
    
    let num_part++
    
    echo "processing $file ($num_part)..." 
    
    sbatch --job-name="replay_$1_$(( $num_part -1 ))" script-replay-partFile $1 $(( $num_part-1 )) 
    
done < "listRaw.$1.lis"

rm listRaw.$1.lis

