select( (.covv_host|ascii_downcase == "human") and
        ( .sequence|length > 24576 ) and 
        ( (.sequence|split("N")|length) < (.sequence|length / 8) ) ) |
"\(.covv_accession_id)\t\(.covv_collection_date)\t\(.covv_subm_date)\t\(.covv_location)\t\(.sequence|split("\n")|join(""))"
