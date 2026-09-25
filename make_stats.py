#!/usr/bin/env python3
"""IBTrACS West-Pacific best tracks -> js/stats.json (used by the sim for genesis + validation).
Setup:  py -3.11 -m pip install pandas numpy requests
Run:    py -3.11 make_stats.py
If the download fails, get ibtracs.WP.list.v04r01.csv from
https://www.ncei.noaa.gov/products/international-best-track-archive and save it next to this script as ibtracs.WP.csv
Definitions (the SAME ones validate.js uses on the simulator):
  storm = max 1-min wind >= 34 kt; track = fixes with wind >= 25 kt; genesis = first such fix;
  recurved = (last lon - westmost lon > 6 deg) and (last lat > lat at westmost point + 4 deg); needs >= 36 h."""
import sys, pathlib, json, numpy as np, pandas as pd, requests
HERE=pathlib.Path(__file__).resolve().parent
JS=HERE/"js" if (HERE/"js").is_dir() else HERE.parent/"js"
URL="https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.WP.list.v04r01.csv"
import shutil, subprocess, time, os
def download(url,dest,tries=6):
    """Download with visible progress, resume, and stall detection (uses curl.exe if present - Windows 10/11 have it)."""
    dest=pathlib.Path(dest);part=dest.with_name(dest.name+".part");curl=None if os.environ.get("NO_CURL") else shutil.which("curl")
    for k in range(1,tries+1):
        print(f"[{k}/{tries}] {url}\n   -> {dest.name}",flush=True)
        if curl:
            rc=subprocess.call([curl,"-L","--fail","--connect-timeout","20","--speed-limit","20000","--speed-time","45","-C","-","-A","Mozilla/5.0","-o",str(part),url])
            if rc==0: part.replace(dest);return
            if rc in (33,) and part.exists(): part.unlink()          # server can't resume -> restart
        else:
            try:
                have=part.stat().st_size if part.exists() else 0;h={"User-Agent":"Mozilla/5.0"}
                if have: h["Range"]=f"bytes={have}-"
                with requests.get(url,headers=h,stream=True,timeout=(20,30)) as r:
                    if r.status_code==416: part.replace(dest);return
                    r.raise_for_status();mode="ab" if r.status_code==206 else "wb";t0=time.time();n=0;tot=int(r.headers.get("content-length",0))+(have if r.status_code==206 else 0)
                    with open(part,mode) as o:
                        for c in r.iter_content(1<<18):
                            o.write(c);n+=len(c)
                            if time.time()-t0>2: print(f"\r   {(have+n)/1e6:.1f} / {tot/1e6:.1f} MB  ({n/1e6/(time.time()-t0):.2f} MB/s)  ",end="",flush=True);t0=time.time()-1e-9 if False else t0
                part.replace(dest);print();return
            except Exception as e: print("\n   stalled/failed:",str(e)[:100],"- retrying (resumes)",flush=True)
        time.sleep(3)
    sys.exit(f"\nDownload failed. Download it in your browser instead and save it as:\n  {dest}\nURL: {url}")
def load_ibtracs():
    f=HERE/"ibtracs.WP.csv"
    if not f.exists() or f.stat().st_size<1e5:
        print("downloading IBTrACS WP (about 30-40 MB) ...");download(URL,f)
    want={"SID","SEASON","NAME","ISO_TIME","LAT","LON","WMO_WIND","USA_WIND","DIST2LAND","TRACK_TYPE"}
    d=pd.read_csv(f,skiprows=[1],usecols=lambda c:c in want,low_memory=False)
    for c in ("LAT","LON","USA_WIND","WMO_WIND","DIST2LAND","SEASON"): d[c]=pd.to_numeric(d[c],errors="coerce")
    d["T"]=pd.to_datetime(d["ISO_TIME"],errors="coerce");d=d.dropna(subset=["T","LAT","LON"])
    # Keep MAIN tracks (final, reanalyzed) and PROVISIONAL/US-PROVISIONAL tracks (real-time
    # data for storms not yet reanalyzed - typically anything within ~2 years of occurring,
    # which is exactly what recent hindcast seasons are). Only drop the "spur" variants
    # (secondary/alternate tracks), which IBTrACS says shouldn't be counted as storms.
    if "TRACK_TYPE" in d: d=d[~d.TRACK_TYPE.astype(str).str.contains("spur",case=False,na=False)]
    d=d[d["T"].dt.hour%6==0].copy();d["LON"]=d["LON"]%360;d["W"]=d["USA_WIND"].fillna(d["WMO_WIND"])
    return d.sort_values(["SID","T"])
def storms(d):
    """yield dict per storm: fixes restricted to first..last wind>=25 kt; only storms reaching 34 kt"""
    for sid,g in d.groupby("SID"):
        w=g["W"].values
        if not np.nanmax(np.where(np.isnan(w),0,w))>=34: continue
        ok=np.where(np.nan_to_num(w)>=25)[0]
        if len(ok)<6: continue
        g=g.iloc[ok[0]:ok[-1]+1]
        yield dict(sid=sid,name=str(g["NAME"].iloc[0]),season=int(g["SEASON"].iloc[0]),T=g["T"].values,lat=g["LAT"].values,lon=g["LON"].values,
                   w=np.nan_to_num(g["W"].values),land=(g["DIST2LAND"].values==0))
def recurved(lon,lat):
    i=int(np.argmin(lon));return bool(lon[-1]-lon[i]>6 and lat[-1]>lat[i]+4)
if __name__=="__main__":
    d=load_ibtracs();Y0,Y1=1991,2020;S=[s for s in storms(d) if Y0<=s["season"]<=Y1];ny=Y1-Y0+1
    cnt=np.zeros(12);rc=np.zeros(12);rn=np.zeros(12);gen=[];spd={};land=0;pk=0;life=[]
    for s in S:
        t0=pd.Timestamp(s["T"][0]);m=t0.month-1;cnt[m]+=1;gen.append([int(t0.dayofyear),round(float(s["lon"][0]),1),round(float(s["lat"][0]),1)])
        dur=(s["T"][-1]-s["T"][0])/np.timedelta64(1,"h")
        if dur>=36: rn[m]+=1;rc[m]+=recurved(s["lon"],s["lat"])
        life.append(dur/24);land+=bool(s["land"].any());pk+=s["w"].max()>=64
        for i in range(1,len(s["lat"])):
            la=(s["lat"][i]+s["lat"][i-1])/2;dx=(s["lon"][i]-s["lon"][i-1])*111*np.cos(np.radians(la));dy=(s["lat"][i]-s["lat"][i-1])*111
            b=int(la//5*5)
            if 5<=b<=45: spd.setdefault(b,[]).append(np.hypot(dx,dy)/6)
    out=dict(years=[Y0,Y1],n_storms=len(S),rate=[round(c/ny,2) for c in cnt],gen=gen,
        recurve_by_month=[round(a/b,3) if b>=5 else None for a,b in zip(rc,rn)],n_by_month=[int(x) for x in rn],
        speed_by_lat={str(k):round(float(np.mean(v)),1) for k,v in sorted(spd.items())},land_share=round(land/len(S),3),
        peak64_share=round(pk/len(S),3),mean_life_days=round(float(np.mean(life)),2))
    JS.mkdir(exist_ok=True);(JS/"stats.json").write_text(json.dumps(out));print("wrote",JS/"stats.json");
    print({k:v for k,v in out.items() if k!="gen"})