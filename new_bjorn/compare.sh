comm -3 $1 $2 | tee >(grep -v $'^\t' > $3) | grep $'^\t' | cut -f2- > $4
