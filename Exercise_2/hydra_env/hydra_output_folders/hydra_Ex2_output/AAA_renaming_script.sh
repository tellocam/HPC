#!/bin/bash
for file in ./*.txt; do
    base=${file%.txt}
    if [[ $base =~ ^.*B[1-9][0-9]*0+$ ]]; then
        zeros=`echo "$base" | grep -o '0' | wc -l`
        sci_num=$(printf "1e%d" $((zeros-1)))
        mv "$file" "${base%B[0-9]*}B$sci_num.txt"
    fi
done
