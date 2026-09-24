#!/usr/bin/env python3
"""Real-weather hindcast data for one season -> js/hind.bin + js/hind.json
Setup:  py -3.11 -m pip install xarray h5netcdf h5py scipy requests numpy pandas
Run:    py -3.11 make_hindcast.py 2019                       (Jun 1 - Nov 30 of that year)
        py -3.11 make_hindcast.py 2019 --start 2019-08-01 --days 60
Weather: NCEP/NCAR Reanalysis daily means at 850/500/200 hPa (NOAA PSL), regridded to 1 deg. Each year file is ~100-200 MB (cached).
Storms:  IBTrACS best tracks (fixes at 6 h, wind >= 25 kt) for storms active in the window.
bin layout (int16, 0.01 m/s): [day][level 850,500,200][u,v][lat 50->-5 (56)][lon 90->180 (91)]"""
import sys, json, argparse, numpy as np, pandas as pd, xarray as xr, requests
from make_stats import HERE, JS, load_ibtracs, storms, download
ap=argparse.ArgumentParser();ap.add_argument("year",type=int);ap.add_argument("--start");ap.add_argument("--days",type=int,default=183);a=ap.parse_args()
start=pd.Timestamp(a.start or f"{a.year}-06-01");year=start.year
CACHE=HERE/"ncep_cache";CACHE.mkdir(exist_ok=True)
BASE="https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/"
def load(v):
    f=CACHE/f"{v}.{year}.nc"
    if not f.exists() or f.stat().st_size<1e6:
        download(BASE+f.name,f)
    errs=[]
    for eng in ("h5netcdf","netcdf4","scipy"):
        try: return xr.open_dataset(f,engine=eng,decode_times=False)[v]
        except Exception as e: errs.append(f"{eng}: {str(e)[:120]}")
    sys.exit("Could not open "+str(f)+"\n  "+"\n  ".join(errs))
u,v=load("uwnd"),load("vwnd")
dates=pd.to_datetime("1800-01-01")+pd.to_timedelta(u["time"].values,unit="h")
sel=np.where((dates>=start)&(dates<start+pd.Timedelta(days=a.days)))[0]
if len(sel)<10: sys.exit("not enough days in that window for this year file")
days=len(sel);lat=np.arange(50,-6,-1.0);lon=np.arange(90,181,1.0);out=[]
print("reading",days,"days from the file (disk/CPU work, no network - can take a minute) ...",flush=True)
u=u.isel(time=sel).sel(level=[850,500,200],lat=slice(90,-20),lon=slice(60,200)).load();v=v.isel(time=sel).sel(level=[850,500,200],lat=slice(90,-20),lon=slice(60,200)).load()
for n in range(days):
    if n%20==0: print(f"  regridding day {n}/{days}",flush=True)
    for lev in (850,500,200):
        for da in (u,v): out.append(da.isel(time=n).sel(level=lev).interp(lat=lat,lon=lon).values)
arr=np.array(out).reshape(days,3,2,len(lat),len(lon));assert not np.isnan(arr).any(),"NaN after regrid"
JS.mkdir(exist_ok=True);(np.round(arr*100).astype("<i2")).tofile(JS/"hind.bin")
end=start+pd.Timedelta(days=days);obs=[]
for s in storms(load_ibtracs()):
    T=pd.to_datetime(s["T"])
    if T[-1]<start or T[0]>=end: continue
    nm=s["name"] if s["name"] not in ("NOT_NAMED","nan","UNNAMED") else s["sid"][-5:]
    obs.append(dict(name=nm.title(),sid=s["sid"],trk=[[round((t-start)/pd.Timedelta(hours=1),1),round(float(lo),2),round(float(la),2),float(w)] for t,lo,la,w in zip(T,s["lon"],s["lat"],s["w"])]))
(JS/"hind.json").write_text(json.dumps(dict(start=str(start.date()),days=days,storms=obs)))
print(f"wrote hind.bin ({days} days) and hind.json ({len(obs)} storms: "+", ".join(o["name"] for o in obs)+")")