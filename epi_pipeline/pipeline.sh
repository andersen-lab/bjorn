#!/usr/bin/env bash

rundir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
workdir=$1
mkdir -p "$workdir"

cd "$workdir"
>&2 echo "downloading global cases from JHU"
wget -Nq "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
>&2 echo "downloading global losses from JHU"
wget -Nq "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
>&2 echo "downloading county-level cases from JHU"
wget -Nq "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
>&2 echo "downloading county-level losses from JHU"
wget -Nq "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv"
>&2 echo "downloading loc lookup table from JHU"
wget -Nq "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv"

##>&2 echo "downloading country-level geoJSON"
##wget -Nq "https://raw.githubusercontent.com/outbreak-info/biothings_covid19/master/geo/countries.json"
##jq -r '.features[]|"\(.properties.location_id)\t\(.properties.NAME)\t\(.geometry)"' "countries.json" > "countries_geo.tsv"
##>&2 echo "downloading US-state-level geoJSON"
##wget -Nq "https://raw.githubusercontent.com/outbreak-info/biothings_covid19/master/geo/US_states.json"
##jq -r '.features[]|"\(.properties.location_id)\t\(.properties.state_id)\t\(.geometry)"' "US_states.json" > "states_geo.tsv"
##>&2 echo "downloading US-county-level geoJSON"
##wget -Nq "https://raw.githubusercontent.com/outbreak-info/biothings_covid19/master/geo/US_counties.json"
##jq -r '.features[]|"\(.properties.location_id)\t\(.properties.state_id)\t\(.geometry)"' "US_counties.json" > "counties_geo.tsv"

>&2 echo "starting pipeline"
python3 $rundir/metamatcher.py $workdir $rundir | parallel --pipe --lb -N32 -j10 $rundir/exploder.py > "$workdir/epi_data.jsonl"
>&2 echo "processing complete"

>&2 echo "casting jsonl to csv"
header=' "date","location_id","iso3","admin1","admin2","name","admin_level","confirmed_per_100k","confirmed_rolling","confirmed_rolling_per_100k","confirmed_rolling_14days_ago_diff","confirmed_rolling_14days_ago_diff_per_100k","dead_per_100k","dead_rolling","dead_rolling_per_100k","dead_rolling_14days_ago_diff","dead_rolling_14days_ago_diff_per_100k"'
getter='[.date, .location_id, .iso3, .admin1, .admin2, .name, .admin_level, .confirmed_per_100k, .confirmed_rolling, .confirmed_rolling_per_100k, .confirmed_rolling_14days_ago_diff, .confirmed_rolling_14days_ago_diff_per_100k, .dead_per_100k, .dead_rolling, .dead_rolling_per_100k, .dead_rolling_14days_ago_diff, .dead_rolling_14days_ago_diff_per_100k]'
echo $header > "$workdir/counties_data.csv"
parallel --pipepart -j10 --quote jq -cr "select( .admin_level==2 and .iso3=="'"'"USA"'"'" ) | $getter | @csv" :::: "$workdir/epi_data.jsonl" >> "$workdir/counties_data.csv"
echo $header > "$workdir/metros_data.csv"
parallel --pipepart -j10 --quote jq -cr "select( .admin_level==1.5 and .iso3=="'"'"USA"'"'" ) | $getter | @csv" :::: "$workdir/epi_data.jsonl" >> "$workdir/metros_data.csv"
echo $header > "$workdir/states_data.csv"
parallel --pipepart -j10 --quote jq -cr "select( .admin_level==1 and .iso3=="'"'"USA"'"'" ) | $getter | @csv" :::: "$workdir/epi_data.jsonl" >> "$workdir/states_data.csv"
echo $header > "$workdir/countries_data.csv"
parallel --pipepart -j10 --quote jq -cr "select( .admin_level==0 ) | $getter | @csv" :::: "$workdir/epi_data.jsonl" >> "$workdir/countries_data.csv"
>&2 echo "calculating breaks"
Rscript $rundir/breaks_and_gifs.R $workdir/
>&2 echo "integrating breaks into jsonl"
jq -c --slurpfile breaks <(jq 'map(del(.id, .location))' "$workdir/breaks.json") '. + $breaks[][.admin_level]' "$workdir/epi_data.jsonl" | gzip > $1/out.jsonl.gz
>&2 echo "output data written to $1/out.jsonl.gz"
