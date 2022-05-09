BEGIN{FS="\t"; OFS="\t"}{ $1 = sprintf("%09i", substr($1, 9)); $5 = substr($5, 14)}1
