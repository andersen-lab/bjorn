mkdir -p /dev/shm/parallel && \
parallel -j12 --tmpdir /dev/shm/parallel --linebuffered \
  " xz -dc $1 | sed -n '{}~12p' | buffer -s128K -m128M | \
    parallel --pipe --roundrobin -j5 --linebuffered \
     ' jq --unbuffered -crf ./filterextract_gisaid.jq ' | \
    python ./compresslines.py --level 3 --dict ./dictionary.zstd | \
    awk -f ./trim_bytes.awk " \
  ::: 0 1 2 3 4 5 6 7 8 9 10 11 | \
parsort -T /dev/shm/parallel
