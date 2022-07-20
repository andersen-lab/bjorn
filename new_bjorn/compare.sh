comm -3 $1 $2 | tee >(grep -v $'^\t' | cut -f1 > $6) | grep $'^\t' | cut -f2- > $4
grep -vfF <(cut -f1 $4) $6 > $3
