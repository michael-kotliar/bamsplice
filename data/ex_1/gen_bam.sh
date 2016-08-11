#!/bin/bash
echo "Searching for all sam files"
ls | grep .*.sam | while read -r SAM ; do
    EXT="${SAM##*.}"
    NAME="${SAM%.*}"
    BAM=${NAME}.bam
    echo "Convert ${SAM} to ${BAM}"
    samtools view -b -S -o ${BAM} ${SAM}
    echo "Sort ${BAM}"
    samtools sort ${BAM} -o temp_sorted
    echo "Remove temporal files"
    rm ${BAM}
    mv temp_sorted ${BAM}
    echo "Finish"
done
