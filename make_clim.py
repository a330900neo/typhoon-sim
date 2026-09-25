#!/usr/bin/env python3
"""Build js/clim.bin (monthly mean winds at 850/500/200 hPa, 1 deg grid, 90-180E, 50N..-5S).
Source: NCEP/NCAR Reanalysis long-term monthly means (NOAA PSL, no login needed).
Setup:  pip install xarray h5netcdf h5py scipy requests numpy
Run:    python make_clim.py          (works from any folder; finds the js/ folder next to it)
Files are cached in ./ncep_cache so a re-run does not download again.
Layout (int16, 0.01 m/s): [month 12][level 850,500,200][u,v][lat 50->-5 (56)][lon 90->180 (91)]"""
import sys, pathlib, numpy as np, xarray as xr, requests
HERE=pathlib.Path(__file__).resolve().parent
JS=HERE/"js" if (HERE/"js").is_dir() else HERE.parent/"js"
CACHE=HERE/"ncep_cache";CACHE.mkdir(exist_ok=True)
BASE="https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/"
def load(v):
    f=CACHE/(v+".mon.ltm.nc")
    if not f.exists() or f.stat().st_size<1e5:
        print("downloading",v,"...");r=requests.get(BASE+v+".mon.ltm.nc",timeout=600);r.raise_for_status();f.write_bytes(r.content)
    errs=[]
    for eng in ("h5netcdf","netcdf4","scipy"):
        try: return xr.open_dataset(f,engine=eng,decode_times=False)[v].load()   # climatology time axis is year 0001 -> can't be decoded; months are in order
        except Exception as e: errs.append(f"{eng}: {type(e).__name__}: {str(e)[:150]}")
    sys.exit(f"Could not open {f}:\n  "+"\n  ".join(errs)+"\nTry: python -m pip install h5netcdf h5py")
lat=np.arange(50,-6,-1.0);lon=np.arange(90,181,1.0);out=[]
u,v=load("uwnd"),load("vwnd")
for m in range(12):
    for lev in (850,500,200):
        for da in (u,v):
            out.append(da.isel(time=m).sel(level=lev).interp(lat=lat,lon=lon).values)
a=np.array(out).reshape(12,3,2,len(lat),len(lon))
assert not np.isnan(a).any(),"NaN in interpolated field"
JS.mkdir(exist_ok=True);(np.round(a*100).astype("<i2")).tofile(JS/"clim.bin");print("wrote",JS/"clim.bin",a.shape)