import pandas as pd
d=pd.read_csv("data/attpc_all_clusters.csv")
d["gid"]=d.run*100000+d.event
d=d.rename(columns={"cluster":"label"})
d=d[["run","event","gid","x","y","z","q","label"]]
d.to_parquet("data/attpc_all.parquet")
print(f"attpc_all.parquet: {d.gid.nunique()} events, {len(d)} pts, runs {sorted(d.run.unique())}")
print("tracks/event med:", int(d.groupby('gid').label.nunique().median()))
