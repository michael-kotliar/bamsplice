#!/bin/bash
echo "Searching for all sam files"
ls | grep .*.sam | while read -r SAM ; do
    EXT="${SAM##*.}"
    NAME="${SAM%.*}"
    BAM=${NAME}.bam
    echo "Convert ${SAM} to ${BAM}"
    samtools view -b -S -o ${BAM} ${SAM}
    echo "Sort ${BAM}"
    samtools sort ${BAM} temp_sorted
    rm ${BAM}
    mv temp_sorted.bam ${BAM}
done
