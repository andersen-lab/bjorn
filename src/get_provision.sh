#!/usr/bin/env bash

cd $1
for filename in *; do
    if [[ $filename > `date -d "$2" '+%Y-%m-%d'` ]]; then
        cat $filename
    fi
done
cd - > /dev/null
