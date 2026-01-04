rm replay/test/test.$1.root

ls  /cache/mss/halla/apex/raw/apex_$1.dat.* > listRaw.$1.lis

$num_part=0

while IFS= read -r file; do
    let num_part++
    echo "processing $file ($num_part)..." 
    analyzer -l -b -q 'setup_forTest.C('$1', "'$file'", "replay/test/test_'$num_part'.'$1'.root", '$2')' 
done < "listRaw.$1.lis"

rm listRaw.$1.lis

hadd replay/test/test.$1.root replay/test/test_*.$1.root
rm replay/test/test_*.$1.root
