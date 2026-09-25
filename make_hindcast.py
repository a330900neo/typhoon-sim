#!/usr/bin/env python3
"""Real-weather hindcast data for one season -> js/hind.bin + js/hind_sst.bin + js/hind.json
Setup:  py -3.11 -m pip install xarray h5netcdf h5py scipy requests numpy pandas
Run:    py -3.11 make_hindcast.py 2019                       (Jun 1 - Nov 30 of that year)
        py -3.11 make_hindcast.py 2019 --start 2019-08-01 --days 60
Weather: NCEP/NCAR Reanalysis daily means at 850/500/200 hPa (NOAA PSL), regridded to 1 deg. Each year file is ~100-200 MB (cached).
SST:     NOAA ERSSTv5 monthly means (NOAA PSL, single ~small file covering 1854-present) for the exact calendar
         months the window spans - gives the sim the real ENSO/SST state for that year instead of a generic
         climatology. Sub-monthly variability (the storm's own cold wake) is still handled dynamically in sim.js.
Storms:  IBTrACS best tracks (fixes at 6 h, wind >= 25 kt) for storms active in the window.
bin layout (int16, 0.01 m/s): [day][level 850,500,200][u,v][lat 50->-5 (56)][lon 90->180 (91)]
hind_sst.bin layout (int16, 0.01 degC): [month][lat 50->-5 (56)][lon 90->180 (91)]; month0/year0/n in hind.json"""
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
end=start+pd.Timedelta(days=days)

# --- real SST (ERSSTv5 monthly, 2 deg): tiny download, gives the exact ENSO/SST state for this year+season ---
SST_URL="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
def load_sst():
    f=CACHE/"sst.mnmean.nc"
    if not f.exists() or f.stat().st_size<1e6: download(SST_URL,f)
    errs=[]
    for eng in ("h5netcdf","netcdf4","scipy"):
        try: return xr.open_dataset(f,engine=eng)["sst"]
        except Exception as e: errs.append(f"{eng}: {str(e)[:120]}")
    print("Could not open "+str(f)+" for SST ("+"; ".join(errs)+") - continuing without real SST, sim will use analytic climatology")
    return None
print("fetching real SST (ERSSTv5 monthly, small download the first time)...",flush=True)
sda=load_sst()
if sda is not None:
    if "zlev" in sda.dims: sda=sda.isel(zlev=0)
    sda=sda.fillna(27.0)   # land/ice cells -> a plain ocean placeholder (never actually read; sim masks by its own terrain, not this)
    m0=pd.Timestamp(year=start.year,month=start.month,day=1);m1=pd.Timestamp(year=end.year,month=end.month,day=1)
    months=pd.date_range(m0,m1,freq="MS")
    sst_fields=[]
    for ts in months:
        try: sst_fields.append(sda.sel(time=ts,method="nearest").interp(lat=lat,lon=lon).values)
        except Exception: sst_fields.append(np.full((len(lat),len(lon)),27.0))
    sst_arr=np.nan_to_num(np.array(sst_fields),nan=27.0)
    (np.round(sst_arr*100).astype("<i2")).tofile(JS/"hind_sst.bin")
    sst_meta=dict(month0=int(months[0].month-1),year0=int(months[0].year),n=len(months))
    print("wrote",JS/"hind_sst.bin",sst_arr.shape,"months",str(months[0].date()),"..",str(months[-1].date()))
else:
    sst_meta=None

obs=[]
for s in storms(load_ibtracs()):
    T=pd.to_datetime(s["T"])
    if T[-1]<start or T[0]>=end: continue
    nm=s["name"] if s["name"] not in ("NOT_NAMED","nan","UNNAMED") else s["sid"][-5:]
    obs.append(dict(name=nm.title(),sid=s["sid"],trk=[[round((t-start)/pd.Timedelta(hours=1),1),round(float(lo),2),round(float(la),2),float(w)] for t,lo,la,w in zip(T,s["lon"],s["lat"],s["w"])]))
(JS/"hind.json").write_text(json.dumps(dict(start=str(start.date()),days=days,storms=obs,sst=sst_meta)))
print(f"wrote hind.bin ({days} days) and hind.json ({len(obs)} storms: "+", ".join(o["name"] for o in obs)+")")