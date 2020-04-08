#!/bin/bash -l


outdir="results/macs2/reproducible_peaks"
BEDs="results/macs2"

final_outfile=$outdir"/"all.reproduciblepeaks.merged.bed

echo "the outdir is $outdir"
echo "the final file is $final_outfile"

mkdir -p $outdir


declare -a INTersect_sample


##### ctrl


sample="ctrl"

pairs=("P11457_101_peaks.bed,P11457_102_peaks.bed" "P11457_102_peaks.bed,P11457_103_peaks.bed" "P11457_101_peaks.bed,P11457_103_peaks.bed")


declare -a INTersect_pairs

smpl_outfile=$outdir"/"$sample".reproducible_peaks.bed"
INTersect_sample=("${INTersect_sample[@]}" "$smpl_outfile")

for pair in "${pairs[@]}"; do
    #Create an array in itself, of the two values.
    tmpArr=(`echo $pair | tr ',' ' '`)

    #Then either use Array Indexing on the array ie. 
    #Or place the offsets in their own variables. ie Run=${tmpArr[0] ...;
    echo "intersecting peaks in files"
    echo "$BEDs/${tmpArr[0]} $BEDs/${tmpArr[1]}"

    f1=${tmpArr[0]}
    f2=${tmpArr[1]}

	smplID1=${f1#P11457_}
    smplID1=${smplID1%_peaks.bed}

	smplID2=${f2#P11457_}
    smplID2=${smplID2%_peaks.bed}

    outfile=$outdir"/"$smplID1.$smplID2.$sample.intersect.bed
    echo "$outfile"

    INTersect_pairs=("${INTersect_pairs[@]}" "$outfile")

    intersectBed -a $BEDs"/"${tmpArr[0]} -b $BEDs"/"${tmpArr[1]} -f 0.5 -r > $outfile

done

pairwise_intersect=`echo ${INTersect_pairs[@]}`

bedops --merge $pairwise_intersect > $smpl_outfile


##### LPS 6h


sample="LPS6h"

pairs=("P11457_104_peaks.bed,P11457_105_peaks.bed" "P11457_105_peaks.bed,P11457_106_peaks.bed" "P11457_104_peaks.bed,P11457_106_peaks.bed")


declare -a INTersect_pairs

smpl_outfile=$outdir"/"$sample".reproducible_peaks.bed"
INTersect_sample=("${INTersect_sample[@]}" "$smpl_outfile")

for pair in "${pairs[@]}"; do
    #Create an array in itself, of the two values.
    tmpArr=(`echo $pair | tr ',' ' '`)

    #Then either use Array Indexing on the array ie. 
    #Or place the offsets in their own variables. ie Run=${tmpArr[0] ...;
    echo "$BEDs/${tmpArr[0]} $BEDs/${tmpArr[1]}"

    f1=${tmpArr[0]}
    f2=${tmpArr[1]}

	smplID1=${f1#P11457_}
    smplID1=${smplID1%_peaks.bed}

	smplID2=${f2#P11457_}
    smplID2=${smplID2%_peaks.bed}

    outfile=$outdir"/"$smplID1.$smplID2.$sample.intersect.bed

    INTersect_pairs=("${INTersect_pairs[@]}" "$outfile")

    intersectBed -a $BEDs"/"${tmpArr[0]} -b $BEDs"/"${tmpArr[1]} -f 0.5 -r > $outfile

done

pairwise_intersect=`echo ${INTersect_pairs[@]}`

bedops --merge $pairwise_intersect > $smpl_outfile



##### LPS 24h


sample="LPS24h"

pairs=("P11457_107_peaks.bed,P11457_108_peaks.bed" "P11457_108_peaks.bed,P11457_109_peaks.bed" "P11457_107_peaks.bed,P11457_109_peaks.bed")


declare -a INTersect_pairs

smpl_outfile=$outdir"/"$sample".reproducible_peaks.bed"
INTersect_sample=("${INTersect_sample[@]}" "$smpl_outfile")

for pair in "${pairs[@]}"; do
    #Create an array in itself, of the two values.
    tmpArr=(`echo $pair | tr ',' ' '`)

    #Then either use Array Indexing on the array ie. 
    #Or place the offsets in their own variables. ie Run=${tmpArr[0] ...;
    echo "$BEDs/${tmpArr[0]} $BEDs/${tmpArr[1]}"

    f1=${tmpArr[0]}
    f2=${tmpArr[1]}

	smplID1=${f1#P11457_}
    smplID1=${smplID1%_peaks.bed}

	smplID2=${f2#P11457_}
    smplID2=${smplID2%_peaks.bed}

    outfile=$outdir"/"$smplID1.$smplID2.$sample.intersect.bed

    INTersect_pairs=("${INTersect_pairs[@]}" "$outfile")

    intersectBed -a $BEDs"/"${tmpArr[0]} -b $BEDs"/"${tmpArr[1]} -f 0.5 -r > $outfile

done

pairwise_intersect=`echo ${INTersect_pairs[@]}`

bedops --merge $pairwise_intersect > $smpl_outfile

### merged reproducible

smpl_reprodicible=`echo ${INTersect_sample[@]}`

bedops --merge $smpl_reprodicible > $final_outfile
