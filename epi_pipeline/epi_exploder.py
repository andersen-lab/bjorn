#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import sys

# load from upstream
print("reading into parallel", file=sys.stderr)
all_data = pd.read_json("/dev/stdin", lines=True)

# explode by dates:
print("exploding by dates", file=sys.stderr)
def Exploder():
    data = []
    def f(row):
        nonlocal data
        ndates = min(len(row["cases"]), len(row["losses"]))
        c = row["cases"]
        l = row["losses"]
        d = row["dates_list"]
        y = list(row.drop(labels=["cases", "losses", "dates_list"]))
        data += [(y+[c[n-21:n+1], l[n-21:n+1], d[n], n == ndates-1]) for n in range(21, ndates)]
        return None
    f.data = lambda: data
    return f
exploder = Exploder()
col_header = list(all_data.drop(columns=["cases", "losses", "dates_list"]).columns)
all_data.apply(exploder, axis=1)
all_data = pd.DataFrame(exploder.data(), columns=col_header + ["cases", "losses", "date", "mostRecent"])

# compute basic statistics for all rows
print("calculating stats", file=sys.stderr)
def calc_stats(data):
   a = data.apply(lambda x: x[-1]).to_numpy()
   b = data.apply(lambda x: x[-2]).to_numpy()
   r = data.apply(lambda x: np.sum(x[-7:])).to_numpy()
   q = data.apply(lambda x: np.sum(x[-21:-14])).to_numpy()
   return pd.DataFrame([ a, a-b, (a-b)/(b+1.), r, q, r-q ]).transpose()
all_data[["confirmed", "confirmed_numIncrease", "confirmed_pctIncrease", "confirmed_rolling",
          "confirmed_rolling_14days_ago", "confirmed_rolling_14days_ago_diff"]] = calc_stats(all_data["cases"])
all_data[["dead", "dead_numIncrease", "dead_pctIncrease", "dead_rolling",
          "dead_rolling_14days_ago", "dead_rolling_14days_ago_diff"]] = calc_stats(all_data["losses"])
all_data = all_data.drop(columns=["cases", "losses"])
# compute population-relative statistics for all rows
data_pop = ~pd.isna(all_data["population"]) & all_data["population"] > 0
all_data.loc[data_pop, ["confirmed_per_100k", "confirmed_numIncrease_per_100k",
           "confirmed_rolling_per_100k", "confirmed_rolling_14days_ago_per_100k",
           "confirmed_rolling_14days_ago_diff_per_100k",
           "dead_per_100k", "dead_numIncrease_per_100k", "dead_rolling_per_100k",
           "dead_rolling_14days_ago_per_100k", "dead_rolling_14days_ago_diff_per_100k"]] = all_data.loc[data_pop, [
                "confirmed", "confirmed_numIncrease", "confirmed_rolling", "confirmed_rolling_14days_ago",
                "confirmed_rolling_14days_ago_diff", "dead", "dead_numIncrease", "dead_rolling",
                "dead_rolling_14days_ago", "dead_rolling_14days_ago_diff" ]].to_numpy() / \
                np.tile(all_data.loc[data_pop, ["population"]].to_numpy(), (1, 10))
# doubling rates
np.seterr(divide='ignore')
all_data["confirmed_doublingRate"] = 14 * np.log(2) / np.log(all_data["confirmed_rolling"] / all_data["confirmed_rolling_14days_ago"])
all_data["dead_doublingRate"] = 14 * np.log(2) / np.log(all_data["dead_rolling"] / all_data["dead_rolling_14days_ago"])

# fill blank with None and export
print("exporting", file=sys.stderr)
all_data = all_data.fillna(np.nan).replace([np.nan], [None])
all_data.to_json("/dev/stdout", orient="records", lines=True)