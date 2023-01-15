#!/bin/bash
for file in ./*.txt; do
    if [[ $file =~ ^.*B[1-9][0-9]*0+\.txt$ ]]; then
        base=${file%.txt}
        zeros=`echo "$base" | grep -o '0' | wc -l`
        sci_num=`echo "1$(printf %0.${zeros}d 0 | tr 0 '*10^-')" | bc -l`
        mv "$file" "${base%[0-9]*}B$sci_num.txt"
    fi
done
