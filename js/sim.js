// v3: MEAN flow = monthly climatology (js/clim.bin from reanalysis if present, else analytic). ANOMALIES evolve with a
// barotropic vorticity equation on a 1° spherical grid (beta effect, advection by mean+eddy flow, eddy advection of mean
// vorticity, damping, sponge), stirred by random vorticity sources -> waves, troughs, ridge breaks arise dynamically.
// Storms are steered by the evolving total flow (one-way coupling). Older notes:
// Streamfunction version: each level's eddy flow = curl of a streamfunction (chain of subtropical-high cells that
// drift/weaken/break, eastward Rossby wave train, transient troughs, monsoon gyre, TUTT) + zonal-mean flow.
// Seasonal grid typhoon model. Three atmospheric levels are modelled separately:
//  850 hPa (low-level: trades, monsoon trough, ridge gyre), 500 hPa (mid: subtropical ridge
//  steering), 200 hPa (upper: subtropical jet, tropical easterly jet, TUTT). The surface wind is
//  derived from 850 hPa (friction: ~0.78x speed, backed 20° across the isobars, and terrain drag).
// Horse latitudes = subtropical high belt (light winds, high pressure) at the ridge latitude,
// which migrates seasonally (~19°N in winter -> ~30°N in Aug). Storm motion uses a deep-layer mean:
// weak storms follow the low levels, intense storms feel the deep layer (so they recurve into the
// westerlies). Vertical shear = |V200 - V850| drives intensity. Date sets SST, ridge, jet, monsoon
// trough, trough frequency and genesis rate (peak Aug-Oct, near zero Jan-Mar).
(function(){
const C=TS.cfg,GW=C.gw,GH=C.gh,N=GW*GH,D=(C.lon1-C.lon0)/GW,R=Math.PI/180,
 clamp=(x,a,b)=>x<a?a:x>b?b:x, rn=()=>(Math.random()+Math.random()+Math.random()-1.5)*2, sg=x=>1/(1+Math.exp(-x));
const S=TS.sim={storms:[],hist:[],t:0,t0:Date.UTC(2026,7,10),running:false,nextSpawn:0,uid:0,E:null,wx:'synth'};
// F.lv = [sfcU,sfcV,u850,v850,u500,v500,u200,v200]
const F=S.F={lv:Array.from({length:8},()=>new Float32Array(N)),p:new Float32Array(N),q:new Float32Array(N),
 q2:new Float32Array(N),sst0:new Float32Array(N),sa:new Float32Array(N),pr:new Float32Array(N),el:null};
F.u=F.lv[0];F.v=F.lv[1];
const NAMES=['Haishen','Noul','Dolphin','Kujira','Chan-hom','Linfa','Nangka','Saudel','Molave','Goni','Vamco','Krovanh','Dujuan','Surigae','Choi-wan','Koguma','Champi','In-fa','Cempaka','Nepartak','Lupit','Mirinae','Nida','Omais','Conson','Chanthu','Dianmu','Mindulle','Lionrock','Kompasu'];
// Genesis starts as an unnamed invest (JTWC-style 90W..99W numbering); a real name is only assigned once the
// system actually reaches tropical-storm strength (34 kt), same as real operational practice.
let investN=90,nameSeq=0;
function nameStorm(s){s.name=NAMES[nameSeq%NAMES.length];nameSeq++;s.named=true}
const GS=.25/D,RR=Math.round(6*GS),SC=16*GS*GS;   // grid-resolution scale factors (tuned at 0.25 deg)
const TH=[34,48,64,85,105],catOf=k=>{let c=0;while(c<5&&k>=TH[c])c++;return c};
// intensity-change node with +-2 kt hysteresis so a storm hovering on a threshold doesn't spam nodes
function catEv(s){const k=s.v*1.944;s.ev=s.ev||[];if(s.cat==null){s.cat=catOf(k);if(s.cat>=1&&!s.named)nameStorm(s);return}let c=s.cat;
  while(c<5&&k>=TH[c]+2)c++;while(c>0&&k<TH[c-1]-2)c--;
  if(c!==s.cat){if(c>=1&&!s.named)nameStorm(s);s.ev.push({t:'i',lon:s.lon,lat:s.lat,c,up:c>s.cat,h:S.t});s.cat=c}}
const pOf=v=>1010-Math.pow(v*1.944/6.7,1/.644);
const cellIdx=S.cellIdx=(lon,lat)=>{const i=Math.floor((lon-C.lon0)/D),j=Math.floor((C.lat1-lat)/D);return i<0||j<0||i>=GW||j>=GH?-1:j*GW+i};

// ---------------- calendar / climatology ----------------
const doyOf=t=>{const d=new Date(S.t0+t*36e5);return(d-Date.UTC(d.getUTCFullYear(),0,0))/864e5};
const seas=t=>.5+.5*Math.cos(2*Math.PI*(doyOf(t)-205)/365);            // 0 = winter, 1 = late-July peak
const clim=(sn,en=0)=>({rx:150-10*sn+5*en,ry:19+11*sn,ra:1.12-.14*sn-.08*en,jl:30+8*sn,mtx:132,mty:6+7*sn});
const act=()=>{const d=doyOf(S.t),x=Math.min(Math.abs(d-240),365-Math.abs(d-240));return .05+.95*Math.exp(-((x/65)**2))};
function setSST(sn){const useReal=hindOn()&&TS.hindSST,lc=6+6*sn,b=29.1+1.1*sn,k=.0125-.005*sn,en=S.E.enso;
  for(let j=0;j<GH;j++){const lat=C.lat1-(j+.5)*D,dl=lat-lc,base=b-k*dl*dl;
    for(let i=0;i<GW;i++){const lon=C.lon0+(i+.5)*D,idx=j*GW+i;
      if(useReal){const r=realSST(lon,lat);F.sst0[idx]=r!=null?r:base}
      else{ // analytic fallback: ENSO shifts the warm pool - El Nino relaxes/spreads it east with a weaker
            // cold tongue; La Nina compresses it west of ~150E with a cooler equatorial tongue further east
        const east=clamp((lon-150)/25,-1,1),tongue=clamp(1-Math.abs(lat)/4,0,1)*(.6-.35*en),
          scs=.8*clamp((124-lon)/6,0,1);   // South China Sea slightly cooler than open Pacific - smooth ramp, not a hard step, so it never shows as a line
        F.sst0[idx]=base-scs+en*.7*east-tongue}}}}

// ---------------- background = climatological mean + barotropic anomaly dynamics ----------------
const LV=new Float64Array(8),CA=Math.cos(20*R),CS=Math.sin(20*R),AE=6.371e6,OM=7.292e-5,A2H=(AE*R)**2,LF=[.75,1,1.15];let bp=0;
const CW=91,CH=56,NC=CW*CH,CG=Array.from({length:8},()=>new Float32Array(NC)),CP=new Float32Array(NC);
const bil=(a,fx,fy)=>{fx=clamp(fx,0,CW-1.001);fy=clamp(fy,0,CH-1.001);const x=fx|0,y=fy|0,tx=fx-x,ty=fy-y,k=y*CW+x;
  return(a[k]*(1-tx)+a[k+1]*tx)*(1-ty)+(a[k+CW]*(1-tx)+a[k+CW+1]*tx)*ty};
const cs=(m,lon,lat)=>bil(CG[m],lon-C.lon0,C.lat1-lat);
const MU=[0,1,2].map(()=>new Float32Array(NC)),MV=[0,1,2].map(()=>new Float32Array(NC)),
  CRH=[0,1].map(()=>new Float32Array(NC)),   // climatological RH background (850/700 hPa, coarse grid), from js/clim_rh.bin - level 0=850 (low-level), 1=700 (mid-level, Gray's genesis parameter)
  PM=new Float32Array(NC),ZX=new Float64Array(NC),ZY=new Float64Array(NC),
  ZP=new Float64Array(NC),PP=new Float64Array(NC),T1=new Float64Array(NC),Z1=new Float64Array(NC),Z2=new Float64Array(NC),
  rc=new Float64Array(CH),rN=new Float64Array(CH),rS=new Float64Array(CH),den=new Float64Array(CH);
for(let j=0;j<CH;j++){const lat=C.lat1-j;rc[j]=Math.cos(lat*R);rN[j]=Math.cos((lat+.5)*R);rS[j]=Math.cos((lat-.5)*R);den[j]=2/(rc[j]*rc[j])+(rN[j]+rS[j])/rc[j]}
TS.loadClim=async url=>{try{const r=await fetch(url,{cache:'no-store'});if(!r.ok){console.warn('loadClim: HTTP',r.status,url);return false}
  const b=new Int16Array(await r.arrayBuffer());const want=12*3*2*NC;
  if(b.length!==want){console.warn('loadClim: size mismatch, got',b.length,'ints, expected',want,'(check js/clim.bin was rebuilt for this grid)');return false}
  TS.clim=b;return true}catch(e){console.warn('loadClim: error',e);return false}};   // layout: [month12][level 850/500/200][u,v][CH][CW], 0.01 m/s
TS.loadClimRH=async url=>{try{const r=await fetch(url,{cache:'no-store'});if(!r.ok){console.warn('loadClimRH: HTTP',r.status,url);return false}
  const b=new Int16Array(await r.arrayBuffer());const want=12*2*NC;
  if(b.length!==want){console.warn('loadClimRH: size mismatch, got',b.length,'ints, expected',want,'(check js/clim_rh.bin was rebuilt for this grid)');return false}
  TS.climRH=b;return true}catch(e){console.warn('loadClimRH: error',e);return false}};   // layout: [month12][level 850/700][CH][CW], RH % x100 - real background moisture (monsoon trough/ridge pattern), see tools/make_clim.py
TS.tune=TS.tune||{rate:1.4};            // rate: multiplier on observed named-storm genesis rate (sim storms often fail to reach 34 kt)
TS.loadStats=async url=>{try{const r=await fetch(url,{cache:'no-store'});if(!r.ok){console.warn('loadStats: HTTP',r.status,url);return false}TS.stats=await r.json();return true}catch(e){console.warn('loadStats: error',e);return false}};   // js/stats.json from tools/make_stats.py
TS.loadHind=async(bu,ju)=>{try{const[a,b]=await Promise.all([fetch(bu,{cache:'no-store'}),fetch(ju,{cache:'no-store'})]);if(!a.ok||!b.ok){console.warn('loadHind: HTTP',a.status,bu,'/',b.status,ju);return false}
  const j=await b.json(),d=new Int16Array(await a.arrayBuffer()),want=j.days*6*NC;
  if(d.length!==want){console.warn('loadHind: size mismatch, hind.bin has',d.length,'ints, hind.json says days=',j.days,'-> expected',want,'(regenerate both together with tools/make_hindcast.py)');return false}
  TS.hind={...j,t0:Date.parse(j.start+'T00:00:00Z'),data:d};return true}catch(e){console.warn('loadHind: error',e);return false}};   // js/hind.bin+json from tools/make_hindcast.py
TS.loadHindSST=async url=>{try{const r=await fetch(url,{cache:'no-store'});if(!r.ok){console.warn('loadHindSST: HTTP',r.status,url);return false}
  TS.hindSST={data:new Int16Array(await r.arrayBuffer())};return true}catch(e){console.warn('loadHindSST: error',e);return false}};   // js/hind_sst.bin: real monthly-mean SST (ERSSTv5); month0/year0/n metadata lives in hind.json's "sst" field
TS.loadHindRH=async url=>{try{const r=await fetch(url,{cache:'no-store'});if(!r.ok){console.warn('loadHindRH: HTTP',r.status,url);return false}
  TS.hindRH={data:new Int16Array(await r.arrayBuffer())};return true}catch(e){console.warn('loadHindRH: error',e);return false}};   // js/hind_rh.bin: real daily-mean RH (NCEP/NCAR) at 850/700 hPa, same day axis as hind.bin/hind.json
// Real-weather mode is only active while the simulated date lies inside the reanalysis window that was loaded. If the chosen start date
// is outside it (e.g. a future date), or the data runs out mid-run, the sim falls back to climatology + stochastic genesis instead of
// sitting there with genesis disabled.
S.hm=false;S.hOff=false;
const winEnd=()=>TS.hind.t0+(TS.hind.days-1)*864e5,hindOn=()=>S.hm&&!!TS.hind&&S.t0+S.t*36e5<=winEnd();
S.mode=()=>S.wx==='synth'||!TS.hind?'synthetic weather':!S.hm?'no reanalysis for this date: climatology + stochastic genesis':hindOn()?(S.wx==='replay'?'real weather (replay)':'real weather (free-run)'):'reanalysis ended: climatology + stochastic genesis';
const rateNow=()=>{const st=TS.stats;if(!st)return null;const x=doyOf(S.t)/30.4375-.5,m0=((Math.floor(x)%12)+12)%12,f=x-Math.floor(x);return st.rate[m0]*(1-f)+st.rate[(m0+1)%12]*f};
const gap=first=>{const rt=rateNow(),mact=clamp(1+.7*S.E.mjo,.2,2.5);   // MJO-like active/quiet spell modulation of genesis frequency
  if(rt==null)return Math.min(720,(first?Math.random()*60:30+Math.random()*110)/act()/mact);
  return Math.min(1500,-Math.log(1-Math.random())*720/Math.max(.05,rt*TS.tune.rate*mact)*(first?Math.random():1))};
const AU=[0,1,2].map(()=>new Float32Array(NC)),AV=[0,1,2].map(()=>new Float32Array(NC)),ZA=new Float64Array(NC);
function loadAna(){const H=TS.hind,x=clamp((S.t0+S.t*36e5-H.t0)/864e5-.5,0,H.days-1.001),d0=x|0,f=x-d0,D=H.data;   // analysis fields, time-interpolated between daily means
  for(let l=0;l<3;l++)for(let c=0;c<2;c++){const o0=((d0*3+l)*2+c)*NC,o1=o0+6*NC,A=c?AV[l]:AU[l];for(let k=0;k<NC;k++)A[k]=.01*(D[o0+k]*(1-f)+D[o1+k]*f)}
  const ua=AU[1],va=AV[1];ZA.fill(0);                                   // 500 hPa anomaly (analysis - climatology) relative vorticity
  for(let j=1;j<CH-1;j++)for(let i=1;i<CW-1;i++){const k=j*CW+i;
    ZA[k]=(((va[k+1]-MV[1][k+1])-(va[k-1]-MV[1][k-1]))-((ua[k-CW]-MU[1][k-CW])*rc[j-1]-(ua[k+CW]-MU[1][k+CW])*rc[j+1]))/(2*R*AE*rc[j])}}
function hindInit(){const H=TS.hind,off=(S.t0-H.t0)/36e5;for(const o of H.storms){o.done=false;o.T=o.trk.map(p=>p[0]-off)}}
// Real SST for the exact simulated year/month (ERSSTv5 monthly means, time-interpolated), used in place of the
// analytic climatology while hindcast free-run/replay is active. The fast storm-induced cold wake is still
// handled separately and dynamically by F.sa, so a monthly-mean background is all this needs to provide.
function realSST(lon,lat){const M=TS.hind&&TS.hind.sst,D=TS.hindSST&&TS.hindSST.data;if(!M||!D||!M.n)return null;
  const d=new Date(S.t0+S.t*36e5),mIdx=(d.getUTCFullYear()-M.year0)*12+(d.getUTCMonth()-M.month0),
    dim=new Date(Date.UTC(d.getUTCFullYear(),d.getUTCMonth()+1,0)).getUTCDate(),f=(d.getUTCDate()-1+d.getUTCHours()/24)/dim,
    i0=clamp(mIdx,0,M.n-1),i1=clamp(mIdx+1,0,M.n-1),fx=lon-C.lon0,fy=C.lat1-lat,
    a=bil(D.subarray(i0*NC,i0*NC+NC),fx,fy)*.01,b=bil(D.subarray(i1*NC,i1*NC+NC),fx,fy)*.01;
  return a+(b-a)*f}
// Real RH for the exact simulated day (NCEP/NCAR daily means, time-interpolated), used as the background moisture
// field in place of climatology while hindcast free-run/replay is active - the actual monsoon trough position and
// actual dry intrusions for that real day, not just a monthly-mean or synthetic proxy. lvl: 0=850 hPa, 1=700 hPa.
function realRH(lon,lat,lvl){const H=TS.hind,D=TS.hindRH&&TS.hindRH.data;if(!H||!D)return null;
  const x=clamp((S.t0+S.t*36e5-H.t0)/864e5-.5,0,H.days-1.001),d0=x|0,f=x-d0,fx=lon-C.lon0,fy=C.lat1-lat,
    o0=(d0*2+lvl)*NC,o1=(clamp(d0+1,0,H.days-1)*2+lvl)*NC,
    a=bil(D.subarray(o0,o0+NC),fx,fy)*.0001,b=bil(D.subarray(o1,o1+NC),fx,fy)*.0001;
  return a+(b-a)*f}
function obsAt(o,t){const T=o.T,n=T.length;if(t<T[0]||t>T[n-1])return null;let i=1;while(i<n-1&&T[i]<t)i++;
  const f=(t-T[i-1])/(T[i]-T[i-1]||1),p=o.trk[i-1],q=o.trk[i];return[p[1]+(q[1]-p[1])*f,p[2]+(q[2]-p[2])*f,p[3]+(q[3]-p[3])*f]}
function hindSpawn(){for(const o of TS.hind.storms){if(o.done||o.T[0]>S.t||o.T[o.T.length-1]<S.t)continue;o.done=true;
  const p=obsAt(o,S.t),a=obsAt(o,S.t+6),b=a?p:obsAt(o,S.t-6),c=a||p,cl=Math.cos(p[1]*R),ok=b&&c!==b,v=Math.max(10,p[2]/1.944);
  S.storms.push({id:++S.uid,name:o.name,named:true,lon:p[0],lat:p[1],v,pmin:pOf(v),rm:60,mu:ok?(c[0]-b[0])*111*cl/21.6:-3,mv:ok?(c[1]-b[1])*111/21.6:1,nu:0,nv:0,age:0,max:v,land:false,spd:5,sv:[0,0],obs:o,trk:[[p[0],p[1],v]],ev:[]})}}
// analytic climatology (fallback; also gives the schematic mean pressure)
function psiM(lon,lat,o){const E=S.E,cl=Math.cos(lat*R);let a=0,b=0,c=0,g;
  for(const h of[[0,0,1],[-24,-2,.45],[22,1,.5]]){const dx=(lon-E.rx-h[0])*cl/17,dy=(lat-E.ry-h[1])/9;g=110*E.ra*h[2]*Math.exp(-(dx*dx+dy*dy));a+=.9*g;b+=1.25*g;c+=.55*g}
  {const dx=(lon-E.mtx)*cl/8,dy=(lat-E.mty)/6;g=-45*(.15+.85*E.sn)*Math.exp(-(dx*dx+dy*dy));a+=g;b+=.3*g;c-=.5*g}
  {const dx=(lon-145)*cl/8,dy=(lat-18)/6;c-=55*E.sn*Math.exp(-(dx*dx+dy*dy))}
  o[0]=a;o[1]=b;o[2]=c}
const P0=new Float64Array(3),PE=new Float64Array(3),PW=new Float64Array(3),PN=new Float64Array(3),PS=new Float64Array(3);
function evalMean(lon,lat){const h=.25,cl=Math.cos(lat*R),E=S.E,sn=E.sn;
  psiM(lon,lat,P0);psiM(lon+h,lat,PE);psiM(lon-h,lat,PW);psiM(lon,lat+h,PN);psiM(lon,lat-h,PS);
  const tr=sg((E.ry-4-lat)/2.5)*sg((lat-4)/2),wl=sg((lat-E.ry-4)/3),jet=Math.exp(-(((lat-E.jl)/9)**2)),
    Z=[-2.5*tr+5*wl,-tr+11*wl,-9*sn*Math.exp(-(((lat-11)/7)**2))+(62-32*sn)*jet];
  for(let m=0;m<3;m++){LV[2+2*m]=-(PN[m]-PS[m])/(2*h)+Z[m];LV[3+2*m]=(PE[m]-PW[m])/(2*h)/cl}
  bp=1010+5*Math.exp(-(((lat-E.ry)/6)**2))+.06*P0[0]-3*Math.exp(-(((lat-E.mty)/5)**2))-5*sg((lat-45)/5)}
function meanField(){const E=S.E,c=clim(E.sn,E.enso);E.rx=c.rx;E.ry=c.ry;E.ra=c.ra;E.jl=c.jl;E.mtx=c.mtx+E.mo;E.mty=c.mty;
  const x=(doyOf(S.t)-1)/365.25*12-.5,m0=((Math.floor(x)%12)+12)%12,m1=(m0+1)%12,fr=x-Math.floor(x);   // calendar month index, used by both the wind and RH climatologies below
  for(let gy=0;gy<CH;gy++)for(let gx=0;gx<CW;gx++){const k=gy*CW+gx;evalMean(C.lon0+gx,C.lat1-gy);PM[k]=bp;
    for(let l=0;l<3;l++){if(TS.clim){const B=TS.clim,o0=((m0*3+l)*2)*NC+k,o1=((m1*3+l)*2)*NC+k;
        MU[l][k]=.01*(B[o0]*(1-fr)+B[o1]*fr);MV[l][k]=.01*(B[o0+NC]*(1-fr)+B[o1+NC]*fr)}
      else{MU[l][k]=LV[2+2*l];MV[l][k]=LV[3+2*l]}}
    for(let l=0;l<2;l++){if(TS.climRH){const B=TS.climRH,o0=(m0*2+l)*NC+k,o1=(m1*2+l)*NC+k;CRH[l][k]=.0001*(B[o0]*(1-fr)+B[o1]*fr)}
      else CRH[l][k]=clamp(.5+.15*Math.exp(-((( C.lat1-gy-(E.mty+2))/9)**2)),.25,.85)}}   // fallback if clim_rh.bin isn't loaded: a mild moist band near the mean monsoon-trough latitude instead of a flat value
  const Zb=new Float64Array(NC);                                  // mean 500 hPa relative vorticity -> its gradient
  for(let j=1;j<CH-1;j++)for(let i=1;i<CW-1;i++){const k=j*CW+i,c=rc[j];
    Zb[k]=((MV[1][k+1]-MV[1][k-1])-(MU[1][k-CW]*rc[j-1]-MU[1][k+CW]*rc[j+1]))/(2*R*AE*c)}
  for(let j=2;j<CH-2;j++)for(let i=2;i<CW-2;i++){const k=j*CW+i;ZX[k]=(Zb[k+1]-Zb[k-1])/(2*R*AE*rc[j]);ZY[k]=(Zb[k-CW]-Zb[k+CW])/(2*R*AE)}}
function solve(z,ps){for(let it=0;it<15;it++)for(let j=1;j<CH-1;j++){const c=rc[j],cN=rN[j],cS=rS[j],d=den[j];
  for(let i=1;i<CW-1;i++){const k=j*CW+i;ps[k]+=1.8*(((ps[k+1]+ps[k-1])/(c*c)+(cN*ps[k-CW]+cS*ps[k+CW])/c-z[k]*A2H)/d-ps[k])}}}
const NU=2e5,DAMP=1/5.2e5,SPG=1/21600;
function tend(z,o){solve(z,PP);
  for(let j=1;j<CH-1;j++){const c=rc[j],cN=rN[j],cS=rS[j],f=2*OM*c/AE,ax=1/(2*R*AE*c),ay=1/(2*R*AE);
    for(let i=1;i<CW-1;i++){const k=j*CW+i,u=-(PP[k-CW]-PP[k+CW])*ay,v=(PP[k+1]-PP[k-1])*ax,zx=(z[k+1]-z[k-1])*ax,zy=(z[k-CW]-z[k+CW])*ay,
      lap=((z[k+1]+z[k-1]-2*z[k])/(c*c)+(cN*(z[k-CW]-z[k])-cS*(z[k]-z[k+CW]))/c)/A2H,e=Math.min(i,j,CW-1-i,CH-1-j),sp=e<5?SPG*(1-e/5)**2:0;
      o[k]=-((MU[1][k]+u)*zx+(MV[1][k]+v)*zy+u*ZX[k]+v*ZY[k]+f*v)-(DAMP+sp)*z[k]+NU*lap}}}
function dyn(){const dt=1800;                                     // Wicker-Skamarock RK3, 2 steps = 1 h
  for(let n=0;n<2;n++){tend(ZP,T1);for(let k=0;k<NC;k++)Z1[k]=ZP[k]+dt/3*T1[k];tend(Z1,T1);for(let k=0;k<NC;k++)Z2[k]=ZP[k]+dt/2*T1[k];
    tend(Z2,T1);for(let k=0;k<NC;k++)ZP[k]+=dt*T1[k]}}
function blob(lon,lat,r0,amp){for(let j=1;j<CH-1;j++)for(let i=1;i<CW-1;i++){const dx=(C.lon0+i-lon)*rc[j]/r0,dy=(C.lat1-j-lat)/r0,q=dx*dx+dy*dy;if(q<9)ZP[j*CW+i]+=amp*Math.exp(-q)}}
function force(amp){amp=amp==null?1:amp;const E=S.E,sn=E.sn,sgn=()=>Math.random()<.5?-1:1;       // random vorticity stirring: mid-lat waves (more in winter), tropical disturbances (more in summer)
  if(Math.random()<(.5+.5*(1-sn))/30)blob(94+Math.random()*30,E.jl+(Math.random()-.5)*16,8,sgn()*(2+3*Math.random())*1e-5*amp);
  if(Math.random()<(.4+.6*sn)/40)blob(110+Math.random()*60,7+Math.random()*16,5,(Math.random()<.65?1:-1)*(1.2+2*Math.random())*1e-5*amp)}
function composeCG(){solve(ZP,PP);const ay=1/(2*R*AE);
  for(let j=0;j<CH;j++){const ax=ay/rc[j];for(let i=0;i<CW;i++){const k=j*CW+i,ed=i==0||j==0||i==CW-1||j==CH-1,u=ed?0:-(PP[k-CW]-PP[k+CW])*ay,v=ed?0:(PP[k+1]-PP[k-1])*ax;
    for(let l=0;l<3;l++){const rp=S.wx==='replay'&&hindOn();CG[2+2*l][k]=rp?AU[l][k]:MU[l][k]+LF[l]*u;CG[3+2*l][k]=rp?AV[l][k]:MV[l][k]+LF[l]*v}
    const u8=CG[2][k],v8=CG[3][k];CG[0][k]=.78*(u8*CA-v8*CS);CG[1][k]=.78*(u8*CS+v8*CA);CP[k]=PM[k]+4e-7*PP[k]}}}
function env(){const E=S.E;E.sn=seas(S.t);if(S.t%24==0){meanField();E.mjo+=-E.mjo*.05+rn()*.15}
  if(!hindOn())force();else if(S.wx!=='replay')force(.5)}   // free-run: keep new troughs/disturbances forming (at half strength, since we started from a real, already-perturbed state) instead of letting the anomaly relax back to pure climatology
// Deep-layer mean steering on a 3° ring. Weak/sheared storms follow 850/500; intense storms add 200 hPa.
const RING=[[0,0],[3,0],[-3,0],[0,3],[0,-3],[2.2,2.2],[-2.2,2.2],[2.2,-2.2],[-2.2,-2.2]];
function steer(lon,lat,vs){const A=[0,0,0,0,0,0];for(const[dx,dy]of RING)for(let m=0;m<6;m++)A[m]+=cs(m+2,lon+dx,lat+dy)/RING.length;
  const[a,b,c,d,e,f]=A,sh=.8*Math.hypot(e-a,f-b),
    ridge=clamp((lat-(S.E.ry-8))/12,0,1),          // approaching/poleward of the ridge nose -> feel the deep layer even if not intense, so recurvature isn't gated on strength alone
    w=Math.max(clamp((vs-8)/40,0,1),ridge)*(1-clamp(sh/30,0,.6));
  return{u:(1-w)*(.55*a+.45*c)+w*(.25*a+.4*c+.35*e),v:(1-w)*(.55*b+.45*d)+w*(.25*b+.4*d+.35*f),shear:sh}}

const BL=Array.from({length:8},()=>new Float32Array(N)),BP=new Float32Array(N),FR=new Float32Array(N),
  BQ=[0,1].map(()=>new Float32Array(N));   // background RH target on the fine grid (0=850hPa,1=700hPa): real hindcast day > climatology+synoptic-anomaly nudge > flat fallback
let bgOk=false;
const ELS=new Float32Array(N);                                          // 1.75°-box-smoothed elevation for terrain blocking
S.init=function(el){F.el=el;const t=new Float32Array(N),r=Math.round(3*GS),w=2*r+1;
  for(let j=0;j<GH;j++)for(let i=0;i<GW;i++){let a=0;for(let d=-r;d<=r;d++)a+=Math.max(0,el[j*GW+clamp(i+d,0,GW-1)]);t[j*GW+i]=a/w}
  for(let j=0;j<GH;j++)for(let i=0;i<GW;i++){let a=0;for(let d=-r;d<=r;d++)a+=t[clamp(j+d,0,GH-1)*GW+i];ELS[j*GW+i]=a/w}
  S.reset()};
S.reset=function(){S.storms=[];S.hist=[];S.t=0;S.running=false;S.hm=S.wx!=='synth'&&!!TS.hind&&S.t0>=TS.hind.t0&&S.t0<=winEnd();S.hOff=false;bgOk=false;F.sa.fill(0);ZP.fill(0);PP.fill(0);investN=90;nameSeq=0;rPN=0;
  const sn=seas(0),en=clamp(rn()*.6,-1,1);        // en>0 El Niño (ridge retreats east, genesis shifts east), <0 La Niña
  S.E={sn,enso:en,mo:rn()*4,rx:0,ry:0,ra:1,jl:0,mtx:0,mty:0,mjo:0};meanField();
  if(hindOn()){hindInit();loadAna();ZP.set(ZA)}                        // start from the observed 500 hPa anomaly
  else{for(let n=0;n<8;n++)force();for(let n=0;n<4;n++){blob(94+Math.random()*60,12+Math.random()*30,7,(Math.random()<.5?-1:1)*3e-5)}
    for(let h=0;h<96;h++){S.E.sn=sn;force();dyn()}}                     // spin-up so waves/troughs already exist at T+0
  composeCG();setSST(sn);bgField();S.nextSpawn=hindOn()?1e9:S.t+gap(true);if(hindOn())hindSpawn();
  for(let k=0;k<N;k++)F.q[k]=BQ[0][k];fields();precip()};
// Genesis: rejection-sampled from SST, low shear, latitude, monsoon-trough proximity and ENSO-shifted longitude.
// Local relative vorticity of the actual simulated 850 hPa flow (climatological mean + evolving anomaly) at a point.
// Real cyclogenesis happens inside pre-existing tropical waves/monsoon-trough disturbances, not on the open, undisturbed
// ocean - this lets genesis "see" the same troughs/waves the isobars and wind layers are already drawing on screen.
function vort850(lon,lat){const h=1,cl=Math.max(.15,Math.cos(lat*R));
  return((cs(3,lon+h,lat)-cs(3,lon-h,lat))/cl-(cs(2,lon,lat+h)*Math.cos((lat+h)*R)-cs(2,lon,lat-h)*Math.cos((lat-h)*R))/cl)/(2*h*R*AE)}
const MINSEP=850;   // km - independent (non-companion) genesis events are rejected within this range of an existing storm, so systems don't spawn on top of each other
S.spawn=function(near){const E=S.E,st=TS.stats;
  for(let n=0;n<300;n++){let lon,lat;
    if(near){lon=near.lon+(Math.random()<.5?-1:1)*(8+6*Math.random());lat=near.lat+(Math.random()-.5)*6}
    else if(st){const d=doyOf(S.t);let g;for(let q=0;q<80;q++){g=st.gen[Math.random()*st.gen.length|0];const dd=Math.abs(g[0]-d);if(Math.min(dd,365-dd)<=20)break}
      lon=g[1]+rn()*1.5+6*E.enso;lat=g[2]+rn()*1.2}                     // observed genesis points within +-20 days of the calendar date
    else{lon=108+Math.random()*68;lat=5+Math.random()*22}
    const k=cellIdx(lon,lat);if(k<0||F.el[k]>0)continue;
    if(!near){let tooClose=false;for(const o of S.storms){if(Math.hypot((lon-o.lon)*Math.cos(lat*R)*111,(lat-o.lat)*111)<MINSEP){tooClose=true;break}}if(tooClose)continue}
    const sh=.8*Math.hypot(cs(6,lon,lat)-cs(2,lon,lat),cs(7,lon,lat)-cs(3,lon,lat)),
      vf=clamp(1+vort850(lon,lat)*4e4,.2,3),                            // favor spots with a live cyclonic disturbance; disfavor (but don't zero out) ridges
      rh=clamp((F.q[k]-.35)/.3,.12,1),                                  // moist mid-troposphere (Gray's 4th genesis parameter): dry air aloft suppresses convection and genesis
      w=clamp((F.sst0[k]-26.5)/1.5,0,1)*clamp(1-sh/16,0,1)*clamp((lat-5)/5,0,1)*vf*rh*(st&&!near?1:(.25+1.6*Math.exp(-(((lon-E.mtx)/22)**2+((lat-E.mty-2)/6)**2)))*Math.exp(-(((lon-142-14*E.enso)/30)**2)));
    if(Math.random()>=w)continue;const v=12+Math.random()*3,inv=investN;investN=investN>=99?90:investN+1;
    S.storms.push({id:++S.uid,name:'Invest '+inv+'W',named:false,lon,lat,v,pmin:pOf(v),rm:60,mu:-3+rn()*1.5,mv:1+rn(),nu:0,nv:0,age:0,max:v,land:false,spd:5,sv:[0,0],trk:[[lon,lat,v]],ev:[]});return S.storms[S.storms.length-1]}};
// Manual placement (D + click in the UI): drop a fresh TD-strength invest at an exact point, bypassing the
// genesis favorability sampling in S.spawn (the user is choosing the spot on purpose). Still refuses land and
// open water outside the grid, since a cyclone can't spin up there. Returns null (and a reason string) on refusal.
S.spawnAt=function(lon,lat){
  if(!F.el)return{ok:false,why:'terrain not loaded yet'};
  const k=cellIdx(lon,lat);if(k<0)return{ok:false,why:'outside the modeled domain'};
  if(F.el[k]>0)return{ok:false,why:'over land'};
  const v=12+Math.random()*3,inv=investN;investN=investN>=99?90:investN+1;
  const s={id:++S.uid,name:'Invest '+inv+'W',named:false,lon,lat,v,pmin:pOf(v),rm:60,mu:-3+rn()*1.5,mv:1+rn(),nu:0,nv:0,age:0,max:v,land:false,spd:5,sv:[0,0],trk:[[lon,lat,v]],ev:[]};
  S.storms.push(s);return{ok:true,storm:s}};
// Fujiwhara: each storm is advected by the (depth-averaged) flow induced by the others, counter-clockwise about the neighbour (NH), plus a
// weak inward drift at close range. Depth-averaging over the storm's own scale means the steering push is only a small fraction of the
// neighbour's peak tangential wind: ~4 m/s max near 300 km, tapering to roughly 1 m/s by ~800-1000 km and fading out by ~1400 km,
// matching the separation at which real binary-typhoon interaction is typically described as becoming noticeable.
function fuji(s){let U=0,V=0;for(const o of S.storms){if(o===s)continue;
  const dx=(s.lon-o.lon)*Math.cos(s.lat*R)*111,dy=(s.lat-o.lat)*111,r=Math.hypot(dx,dy)+1e-3;
  if(r<1400){const x=Math.pow(o.rm/r,1.4),vt=Math.min(4,.2*o.v*Math.sqrt(x*Math.exp(1-x))*Math.exp(-((r/750)**2))),inw=r<600?.12:0;
    U+=-dy/r*vt-dx/r*vt*inw;V+=dx/r*vt-dy/r*vt*inw}}
  return[U,V]}
S.fuji=fuji;
// ---------------- Forecast: extrapolate the CURRENT flow forward, same physics as upd() ----------------
// A short-range track forecast is "where does this storm go if today's steering pattern verifies" - steering + beta
// drift + storm-storm interaction, run forward WITHOUT the hourly wobble (nu/nv), since that fine-scale turbulence is
// unknowable ahead of time (the same simplification real steady-state/CLIPER-style guidance makes). The flow field
// itself (CG) is held at its current snapshot rather than re-run through dyn() - a short-lead persistence assumption,
// same one real 1-3 day steering nowcasts implicitly make. Cone radius follows the approximate SHAPE of published
// NHC/JTWC mean track-error-by-lead-time curves (present-day skill, km).
const CONE=[[0,0],[24,40],[48,80],[72,125],[96,175],[120,225]];
const coneR=h=>{if(h<=0)return 0;let i=1;while(i<CONE.length-1&&CONE[i][0]<h)i++;
  const a=CONE[i-1],b=CONE[i],t=(h-a[0])/(b[0]-a[0]||1);return a[1]+(b[1]-a[1])*t};
function fujiF(f,set){let U=0,V=0;for(const o of set){if(o===f||o.dead)continue;
  const dx=(f.lon-o.lon)*Math.cos(f.lat*R)*111,dy=(f.lat-o.lat)*111,r=Math.hypot(dx,dy)+1e-3;
  if(r<1400){const x=Math.pow(o.rm/r,1.4),vt=Math.min(4,.2*o.v*Math.sqrt(x*Math.exp(1-x))*Math.exp(-((r/750)**2))),inw=r<600?.12:0;
    U+=-dy/r*vt-dx/r*vt*inw;V+=dx/r*vt-dy/r*vt*inw}}
  return[U,V]}
function nbShearF(f,set){let u8=0,v8=0,u2=0,v2=0;for(const o of set){if(o===f||o.dead)continue;
  const dx=(f.lon-o.lon)*Math.cos(f.lat*R)*111,dy=(f.lat-o.lat)*111,r=Math.hypot(dx,dy)+1e-3;if(r>1400)continue;
  const x=Math.pow(o.rm/r,1.4),vt=o.v*Math.sqrt(x*Math.exp(1-x))*Math.exp(-((r/900)**2)),tu=-vt*dy/r,tv=vt*dx/r;
  u8+=tu;v8+=tv;u2-=tu*.2;v2-=tv*.2}
  return .8*Math.hypot(u2-u8,v2-v8)}
S.forecast=function(hours){hours=hours||120;
  const set=S.storms.map(s=>({id:s.id,name:s.name,lon:s.lon,lat:s.lat,mu:s.mu,mv:s.mv,v:s.v,rm:s.rm,dead:false,trk:[[s.lon,s.lat,s.v,0]]}));
  for(let h=1;h<=hours;h++)for(const f of set){if(f.dead)continue;const k=cellIdx(f.lon,f.lat);if(k<0){f.dead=true;continue}
    const el=F.el[k],land=el>0,sst=F.sst0[k]+F.sa[k],els=ELS[k],gi=k%GW;let tx=0,ty=0;
    if(gi>0&&gi<GW-1&&k>=GW&&k<N-GW){const gx=(ELS[k+1]-ELS[k-1])/2*GS,gy=(ELS[k-GW]-ELS[k+GW])/2*GS,m=Math.hypot(gx,gy);if(m>1){const ff=Math.min(2.5,m/100)/m;tx=-gx*ff;ty=-gy*ff}}
    const slow=1-.3*Math.min(1,els/1000)-(land?.1:0),st=steer(f.lon,f.lat,f.v),
      betaMag=(.6+1.5*clamp(f.rm/90,0,1))*clamp(1.3-f.lat/60,.5,1.3),betaU=-betaMag*.7,betaV=betaMag*.85,
      [fjU,fjV]=fujiF(f,set);
    f.mu+=((st.u+betaU+fjU)*.95*slow+tx-f.mu)/13;f.mv+=((st.v+betaV+fjV)*.95*slow+ty-f.mv)/13;
    f.lon+=f.mu*3.6/(111*Math.cos(f.lat*R));f.lat+=f.mv*3.6/111;
    const k2=Math.max(0,cellIdx(f.lon,f.lat)),shear=st.shear+nbShearF(f,set),rel=clamp((F.q[k2]-.35)/.4,0,1),
      mpi=12+80*clamp((sst-26)/4,0,1),tgt=mpi*clamp(1-shear/32,.12,1)*(.65+.35*rel)*clamp((f.lat-5)/8,0,1),
      rmf=clamp(65/f.rm,.7,1.6),dv=land?-(f.v-8)*(.07+.06*Math.min(1,el/1500)):(tgt-f.v)*(tgt>f.v?.024:.032)*rmf;
    f.v=Math.max(3,f.v+dv);f.rm+=(clamp(65-.55*f.v+(f.lat-12),20,100)-f.rm)*.03;
    if(h%3===0)f.trk.push([f.lon,f.lat,f.v,coneR(h)]);
    if(f.v<10&&(land||f.lat>36))f.dead=true}
  return set.map(f=>({id:f.id,name:f.name,trk:f.trk}))};
// Vertical shear a storm feels from a NEIGHBOUR's own circulation: the exact same Rankine profile fields() uses to
// paint each storm's wind footprint (850 hPa cyclonic core, weak 200 hPa outflow anticyclone), evaluated at this
// storm's centre from every other storm and read back as u200-u850. Real close-packed systems ventilate/shear each
// other through outflow this way - this isn't a separate distance penalty, it's the model's own wind field read
// back at a neighbour's location instead of just this storm's own.
function nbShear(s){let u8=0,v8=0,u2=0,v2=0;for(const o of S.storms){if(o===s)continue;
  const dx=(s.lon-o.lon)*Math.cos(s.lat*R)*111,dy=(s.lat-o.lat)*111,r=Math.hypot(dx,dy)+1e-3;if(r>1400)continue;
  const x=Math.pow(o.rm/r,1.4),vt=o.v*Math.sqrt(x*Math.exp(1-x))*Math.exp(-((r/900)**2)),tu=-vt*dy/r,tv=vt*dx/r;
  u8+=tu;v8+=tv;u2-=tu*.2;v2-=tv*.2}
  return .8*Math.hypot(u2-u8,v2-v8)}
// Merger: centres closer than ~0.9x the sum of the two RMWs (80-170 km, i.e. eyewalls essentially overlapping). Similar intensity -> complete merger (new centre = intensity-weighted, slightly stronger, larger);
// otherwise the stronger storm absorbs the weaker one, gaining a little.
function merge(){for(let i=0;i<S.storms.length;i++)for(let j=i+1;j<S.storms.length;j++){const a=S.storms[i],b=S.storms[j];
  if(a.land||b.land||Math.hypot((a.lon-b.lon)*Math.cos(a.lat*R)*111,(a.lat-b.lat)*111)>Math.min(170,Math.max(80,.9*(a.rm+b.rm))))continue;
  const[big,sm]=a.v>=b.v?[a,b]:[b,a],q=sm.v/big.v;
  if(q>.8){const t=big.v+sm.v;big.lon=(big.lon*big.v+sm.lon*sm.v)/t;big.lat=(big.lat*big.v+sm.lat*sm.v)/t;big.mu=(big.mu*big.v+sm.mu*sm.v)/t;big.mv=(big.mv*big.v+sm.mv*sm.v)/t;
    big.v*=1+.08*q;big.rm=Math.max(big.rm,sm.rm)+10}
  else{big.v*=1+.05*q;big.rm+=4*q}
  big.max=Math.max(big.max,big.v);sm.by=big.name;kill(sm);return true}return false}
function upd(s){const k=cellIdx(s.lon,s.lat);if(k<0)return kill(s);
  const el=F.el[k],land=el>0,sst=F.sst0[k]+F.sa[k];s.land=land;if(land)s.landed=true;
  s.ev=s.ev||[];if(land!==!!s.ls){if(s.cand===land){if(++s.cn>=3){s.ev.push({t:land?'in':'out',lon:s.cpos[0],lat:s.cpos[1],v:s.v,h:S.t});if(land)s.lf=(s.lf||0)+1;s.ls=land;s.cand=null}}else{s.cand=land;s.cn=1;s.cpos=[s.lon,s.lat]}}else s.cand=null;
  const els=ELS[k],gi=k%GW;let tx=0,ty=0;                                // terrain: pushed away from high ground, slower over/near mountains
  if(gi>0&&gi<GW-1&&k>=GW&&k<N-GW){const gx=(ELS[k+1]-ELS[k-1])/2*GS,gy=(ELS[k-GW]-ELS[k+GW])/2*GS,m=Math.hypot(gx,gy);if(m>1){const f=Math.min(2.5,m/100)/m;tx=-gx*f;ty=-gy*f}}
  const slow=1-.3*Math.min(1,els/1000)-(land?.1:0);
  const st=steer(s.lon,s.lat,s.v);s.sv=[st.u,st.v];
  const betaMag=(.6+1.5*clamp(s.rm/90,0,1))*clamp(1.3-s.lat/60,.5,1.3),betaU=-betaMag*.7,betaV=betaMag*.85;
  const[fjU,fjV]=fuji(s);
  const nz=1+clamp((4-Math.hypot(st.u,st.v))/4,0,1);   // persistent (~3 day) wobble, larger when steering is weak
  s.nu=s.nu*.99+rn()*.055*nz;s.nv=s.nv*.99+rn()*.055*nz;   // amplitude cut ~4x (was .22) - real tracks are only faintly trochoidal, not visibly looping
  s.mu+=((st.u+betaU+fjU)*.95*slow+tx+s.nu-s.mu)/13;s.mv+=((st.v+betaV+fjV)*.95*slow+ty+s.nv-s.mv)/13;
  s.lon+=s.mu*3.6/(111*Math.cos(s.lat*R));s.lat+=s.mv*3.6/111;
  const shear=st.shear+nbShear(s),rel=clamp((F.q[k]-.35)/.4,0,1),mpi=12+80*clamp((sst-26)/4,0,1);
  let tgt=mpi*clamp(1-shear/32,.12,1)*(.55+.45*rel)*clamp((s.lat-5)/8,0,1);
  // eyewall replacement cycle: a mature, near-peak storm occasionally spends ~a day reorganising its inner
  // core, capping its target intensity a bit below MPI before it's free to re-strengthen - gives the
  // plateau-then-wobble near peak intensity real best-tracks show, instead of a smooth climb straight to MPI.
  if(!s.erc&&s.v>47&&s.v>s.max*.92&&Math.random()<1/260)s.erc=16+Math.random()*16;
  if(s.erc>0){tgt*=.82;s.erc--}
  const rmf=clamp(65/s.rm,.7,1.6);   // compact (small-RMW) storms intensify/weaken faster than sprawling ones
  const dv=land?-(s.v-8)*(.07+.06*Math.min(1,el/1500)):(tgt-s.v)*(tgt>s.v?.024:.032)*rmf+rn()*.05;
  // Dissipation: once the kill condition is met the storm doesn't vanish, it spins down over ~18h - intensity
  // and drawn opacity both ramp toward zero (s.fade) so a weakening TD/extratropical system visibly winds down,
  // still drifting with the flow, rather than disappearing between two frames.
  if(s.decay){s.v=Math.max(0,s.decayV0*(1-s.decay/18));s.fade=clamp(1-s.decay/18,0,1)}
  else{s.v=Math.max(3,s.v+dv);s.fade=1}
  s.max=Math.max(s.max,s.v);s.pmin=pOf(s.v);catEv(s);
  s.rm+=(clamp(65-.55*s.v+(s.lat-12),20,100)-s.rm)*.03;s.age++;s.spd=Math.hypot(s.mu,s.mv);
  if(s.age%3==0)s.trk.push([s.lon,s.lat,s.v]);
  if(!land&&!s.decay){const ci=Math.floor((s.lon-C.lon0)/D),cj=Math.floor((C.lat1-s.lat)/D),
      mldf=clamp(1.3-clamp((s.lat-10)/20,0,1)*.7,.5,1.3);   // shallower thermocline away from the deep warm pool -> more self-cooling there
    for(let dj=-RR;dj<=RR;dj++)for(let di=-RR;di<=RR;di++){const i=ci+di,j=cj+dj;if(i>=0&&j>=0&&i<GW&&j<GH)
      F.sa[j*GW+i]=Math.max(-4,F.sa[j*GW+i]-.0007*s.v*mldf*Math.exp(-(di*di+dj*dj)/SC))}}
  if(s.obs){const p=obsAt(s.obs,S.t);s.err=p?Math.hypot((s.lon-p[0])*Math.cos(s.lat*R)*111,(s.lat-p[1])*111):null}
  if(!s.decay&&((s.v<10&&(land||s.lat>36||s.age>120))||s.age>650)){s.decay=1;s.decayV0=s.v}
  else if(s.decay){s.decay++;if(s.decay>=18)kill(s)}}
function kill(s){S.storms=S.storms.filter(x=>x!==s);S.hist.push(s);if(S.hist.length>12)S.hist.shift()}
function bgField(){
    const el=F.el,useRealRH=hindOn()&&TS.hindRH;
    for(let j=0;j<GH;j++)for(let i=0;i<GW;i++){const k=j*GW+i,fx=(i+.5)*D,fy=(j+.5)*D,x=fx|0,y=fy|0,tx=fx-x,ty=fy-y,a=y*CW+x,
    w0=(1-tx)*(1-ty),w1=tx*(1-ty),w2=(1-tx)*ty,w3=tx*ty,lon=C.lon0+fx,lat=C.lat1-fy;
    FR[k]=el[k]>0?Math.max(.45,.8-el[k]/12000):1;
    for(let m=0;m<8;m++){const c=CG[m];BL[m][k]=(c[a]*w0+c[a+1]*w1+c[a+CW]*w2+c[a+CW+1]*w3)*(m<2?FR[k]:1)}
    BP[k]=(CP[a]*w0+CP[a+1]*w1+CP[a+CW]*w2+CP[a+CW+1]*w3)-(el[k]>0?2.5+(lon<125&&lat>20?4:0):0);
    // background RH target: real hindcast day when available; otherwise climatology nudged by the evolving 500 hPa
    // vorticity anomaly (a transient trough/ridge moistens/dries on top of the climatological mean, the same
    // anomaly field the winds and isobars already show) - either way this replaces the old SST-only proxy.
    const zv=ZP[a]*w0+ZP[a+1]*w1+ZP[a+CW]*w2+ZP[a+CW+1]*w3;
    for(let l=0;l<2;l++){const rv=useRealRH?realRH(lon,lat,l):null,cl=CRH[l],base=cl[a]*w0+cl[a+1]*w1+cl[a+CW]*w2+cl[a+CW+1]*w3;
      BQ[l][k]=clamp(rv!=null?rv:base+clamp(zv*6e4,-.12,.12),.15,.97)}}
    bgOk=true}
function fields(){const{u,v,p,q,lv}=F;if(!bgOk||S.t%6===0)bgField();for(let m=0;m<8;m++)lv[m].set(BL[m]);p.set(BP);
  for(const s of S.storms){const cl=Math.cos(s.lat*R),dl=12.6/cl,i0=Math.max(0,Math.floor((s.lon-dl-C.lon0)/D)),i1=Math.min(GW-1,Math.ceil((s.lon+dl-C.lon0)/D)),
    j0=Math.max(0,Math.floor((C.lat1-s.lat-12.6)/D)),j1=Math.min(GH-1,Math.ceil((C.lat1-s.lat+12.6)/D));
   for(let j=j0;j<=j1;j++){const lat=C.lat1-(j+.5)*D,c2=Math.cos(lat*R);for(let i=i0;i<=i1;i++){const k=j*GW+i,lon=C.lon0+(i+.5)*D,
    dx=(lon-s.lon)*111*c2,dy=(lat-s.lat)*111,r=Math.max(1,Math.hypot(dx,dy));if(r>1400)continue;
    const x=Math.pow(s.rm/r,1.4),vt=s.v*Math.sqrt(x*Math.exp(1-x))*Math.exp(-((r/900)**2)),tr=.4*Math.exp(-((r/500)**2)),f=FR[k],tu=-vt*dy/r,tv=vt*dx/r;
    u[k]+=(tu-.3*vt*dx/r+tr*s.mu)*f;v[k]+=(tv-.3*vt*dy/r+tr*s.mv)*f;p[k]-=(1010-s.pmin)*(1-Math.exp(-x));
    lv[2][k]+=tu;lv[3][k]+=tv;lv[4][k]+=tu*.45;lv[5][k]+=tv*.45;lv[6][k]-=tu*.2;lv[7][k]-=tv*.2;   // cyclone at 850, weak at 500, outflow anticyclone at 200
    // moisture: convective core moistens (eyewall + spiral rainbands); the compensating subsidence in the storm's
    // OWN secondary circulation dries an annulus around it (real TC structure - the dry moat just outside the
    // eyewall). That annulus reaches out to ~800-1000+ km, so two storms at typical binary/companion separations
    // will each dry out the shared air the other would otherwise draw on - genuine moisture competition arising
    // from each storm's own physics, not a bolt-on distance penalty.
    // Once a storm is organized/intense enough to have earned a real eye, that eye is itself dry and subsiding -
    // real TCs are saturated in the eyewall ring, not moist all the way to the center. eyeStrength ramps up from
    // ~TS strength to a clear eye by the low end of typhoon strength; weak/disorganized systems keep the old
    // filled (moist all the way to the center) profile, since they don't have a real eye to hollow out.
    const fill=Math.exp(-((r/(1.4*s.rm))**2)),eyeStrength=clamp((s.v-23)/27,0,1),eye=eyeStrength*Math.exp(-((r/(.32*s.rm))**2)),
      dry=Math.max(0,Math.exp(-((r/900)**2))-Math.exp(-((r/220)**2)));
    q[k]+=(.95-q[k])*.05*Math.max(0,fill-.9*eye)-q[k]*(.03*dry+.05*eye)}}}}
function humidity(){const{u,v,q,q2,el}=F;
  for(let j=0;j<GH;j++){const lat=C.lat1-(j+.5)*D,ck=D*111,cx=ck*Math.cos(lat*R);
   for(let i=0;i<GW;i++){const k=j*GW+i,x=clamp(i-u[k]*3.6/cx,0,GW-2.001),y=clamp(j+v[k]*3.6/ck,0,GH-2.001),
     x0=x|0,y0=y|0,fx=x-x0,fy=y-y0,a=y0*GW+x0;
    let val=(q[a]*(1-fx)+q[a+1]*fx)*(1-fy)+(q[a+GW]*(1-fx)+q[a+GW+1]*fx)*fy;
    // Background target is now BQ[0] (850 hPa): real hindcast RH when in real-weather mode, else RH climatology
    // nudged by the evolving vorticity anomaly (see bgField()) - carries the real monsoon-trough/ridge and
    // land moisture pattern instead of the old SST-only proxy. Relaxation is slower than before so storm-advected
    // moisture plumes and dry moats persist for hours instead of snapping back to the background immediately.
    const tgt=BQ[0][k];
    if(el[k]>0){val+=(tgt-val)*.03;if(el[k]>300)val-=.02*clamp(el[k]/2000,0,1)*Math.hypot(u[k],v[k])/10}   // still extra-dry on high terrain (thinner boundary layer, rain-shadow)
    else val+=(tgt-val)*.035;
    q2[k]=clamp(val,.05,1)}}
  F.q=q2;F.q2=q}
// ---------------- precipitation (diagnostic, mm/h) ----------------
// Background convective rain from how saturated the (already-advected) moisture field is - this is what paints
// the monsoon trough/ITCZ as a broad rainy band and the subtropical ridge as dry, same q field humidity() just
// updated. Orographic term adds the extra rain wrung out where the low-level flow is forced up a slope (and by
// omission leaves the leeward side a rain shadow) - same terrain data terrain.js/upd() already use.
//
// A storm's OWN rain is no longer a fixed radial mask re-painted from scratch every frame. Instead it's a swarm of
// Lagrangian rain parcels: they are BORN near the eyewall/spiral bands (rejection-sampled with the same shape real
// TC rainfall has - a compact core peaking at the radius of max wind, a dry eye, 2-armed log-spiral bands, a
// downshear-left bias), then ridden by the ACTUAL wind - the average of the surface flow and 850 hPa (the same
// fields the storm's own spiral inflow/tangential wind/outflow are already painted into) - so the rain visibly
// streams, curls and lags the flow instead of being re-stamped in place every step. A parcel fades (dry air, old
// age) or, having spiralled back in on the low-level inflow, is finally swallowed once it reaches the core again -
// born at the storm, carried by its wind, reabsorbed at the storm: a genuine recycling loop, not a static mask.
const RPMAX=20000,RP_R=2,RP_S=4;                        // deposit-kernel radius (cells) & spread pulled well back down - the band's width/curl is now shaped explicitly above, so the old wide kernel + blur was piling on top of that and ballooning everything (core included) far past its intended size
const rLon=new Float32Array(RPMAX),rLat=new Float32Array(RPMAX),rAmt=new Float32Array(RPMAX),rAge=new Float32Array(RPMAX),PRT=new Float32Array(N);
let rPN=0;
const sampG=(a,lon,lat)=>{let fx=(lon-C.lon0)/D-.5,fy=(C.lat1-lat)/D-.5;fx=clamp(fx,0,GW-1.001);fy=clamp(fy,0,GH-1.001);
  const x=fx|0,y=fy|0,tx=fx-x,ty=fy-y,k=y*GW+x;
  return(a[k]*(1-tx)+a[k+1]*tx)*(1-ty)+(a[k+GW]*(1-tx)+a[k+GW+1]*tx)*ty};
function rAdd(lon,lat,amt){const i=rPN<RPMAX?rPN++:(Math.random()*RPMAX)|0;   // pool full -> recycle a random slot rather than growing further
  rLon[i]=lon;rLat[i]=lat;rAmt[i]=amt;rAge[i]=0}
// Birth: the core/eyewall ring is still a light rejection-sampled fill (a handful of random points reads fine
// over a small, roughly-circular area). Feeder bands are NOT rejection-sampled any more - randomly accepting or
// rejecting points scattered over the whole disc left gaps of well over 100 km between hits along any one arm,
// which is exactly what read as "broken/scattered" instead of a band, no matter how the accept odds or the
// after-the-fact blur radius were tuned. Instead, each band is walked DETERMINISTICALLY along its own exact
// log-spiral curve (ang as a function of r - the same relation the old code only used to test closeness): a
// point is dropped every DR km the whole way out, so consecutive seeds are always a short, fixed distance apart
// and the arm is unbroken by construction. Intensity (core/far/downshear-bias) still shapes how much rain each
// step deposits, and a small jitter plus the wind advection in the ride phase below keep the result looking like
// organic weather rather than a drafting-compass curve, but the LINE ITSELF no longer depends on random luck.
function spawnRain(){for(const s of S.storms){if(s.decay)continue;
  const eyeStrength=clamp((s.v-18)/22,0,1),peak=9+100*clamp((s.v-15)/50,0,1),
    // The dynamical RMW (s.rm) barely contracts with intensity (it's shared with steering/wind-field code that
    // needs its own tuning), so anchoring the rain ring directly at s.rm left a huge, essentially rain-free gap
    // between the small fixed-radius centerRain patch and the eyewall for any strong storm. rw contracts the RAIN
    // profile only, in step with eye maturity (eyeStrength, already 0->1 as the storm develops a real eye) - a
    // mature intense storm's visible eyewall pulls in toward the centre the way a real one does; s.rm itself,
    // and everything else that reads it, is untouched.
    rw=s.rm*(1-.55*eyeStrength),
    // Stronger storms get MORE distinct feeder bands (real major typhoons commonly show 2-3 spiral arms wrapped
    // around a well-organized core, not just one), and a thicker/more robust core convective ring (a wider CDO),
    // while weak/moderate storms keep a single thin arm and a tight core - matStr keeps climbing well past where
    // eyeStrength saturates (v=40) so the strongest super typhoons still gain bands other storms don't.
    matStr=clamp((s.v-20)/40,0,1),nBands=1+Math.round(2*matStr),
    wallFrac=.22+.24*eyeStrength,                        // core/eyewall ring thickness grows with intensity: .22*rw (weak) up to .46*rw (mature intense)
    shu=cs(6,s.lon,s.lat)-cs(2,s.lon,s.lat),shv=cs(7,s.lon,s.lat)-cs(3,s.lon,s.lat),
    biasAng=Math.atan2(shv,shu)+Math.PI/2,          // downshear-left quadrant: climatological TC rainfall maximum (NH)
    cl=Math.cos(s.lat*R),
    nCore=Math.max(4,Math.round(peak/3));
  for(let a=0;a<nCore;a++)for(let tries=0;tries<6;tries++){
    const r=Math.max(1,Math.random()*1.15*rw),ang=Math.random()*2*Math.PI,x=Math.pow(rw/r,1.4),core=Math.sqrt(x*Math.exp(1-x)),
      wallR=(r-rw)/(wallFrac*rw),wall=Math.exp(-(wallR*wallR)),
      w=core*.9+eyeStrength*wall*1.3;
    if(Math.random()<w/1.1){rAdd(s.lon+r*Math.cos(ang)/(111*cl),s.lat+r*Math.sin(ang)/111,peak*Math.max(.3,w)*(.7+.5*Math.random()));break}}
  const RMIN=Math.max(8,rw*.85),RMAX=850,DR=20;   // km spacing along each arm's curve - reach pulled in from 1350 to 850 km, a more realistic total rain-shield radius
  for(let bnd=0;bnd<nBands;bnd++){const phase=bnd*(2*Math.PI/nBands);
    for(let r=RMIN;r<=RMAX;r+=DR){
      // Winding rate is no longer one constant for the whole arm. Real tangential wind is close to solid-body
      // rotation inside the RMW and decays outside it, so the ANGULAR speed a parcel of air (and the band riding
      // it) sweeps through is fast near the core and slack far out - a single fixed pitch made every band spin as
      // one rigid rod, same rate at 50 km and 1000 km, which is why it read as an oversized, too-uniform sheet.
      // wind(r) sets the spiral's pitch (tight near the core, loose far out); omega(r) makes that same falloff
      // apply to how fast the pattern rotates in TIME, so the inner curl visibly winds up faster hour to hour
      // than the outer reach of the same arm.
      // Sign matters too: this storm's tangential flow is counterclockwise (spin advances the angle forward in
      // time), and air spiraling INTO a counterclockwise vortex sweeps CLOCKWISE as you trace it outward from the
      // centre - the classic NH pinwheel look. That trailing shape needs a NEGATIVE pitch here; the old constant
      // +1.15 wound the opposite way, which is what made the whole thing look mirrored/off.
      const wind=-(0.5+1.3*rw/(rw+r*.55)),
        omega=-.17*clamp(1.7*rw/(rw+r),.2,1.7),
        angC=phase-omega*S.t+wind*Math.log(r/rw),
        // Width has to be pinned to a roughly constant distance in km, NOT a fixed number of degrees - a fixed
        // degree width times a growing radius is exactly a wedge that gets wider forever, which is what was
        // actually making the whole thing balloon into an oversized pie-slice at long range. Capping the km width
        // (and keeping the cap itself modest) instead keeps the arm reading as a slim ribbon all the way out.
        halfWkm=Math.min(38,12+26*r/500),
        halfW=halfWkm/r;
      for(let sub=0;sub<3;sub++){
        const off=(sub/2-.5)*2*halfW,taper=Math.cos(clamp(off/halfW,-1,1)*Math.PI/2)**2,
          ang=angC+off,
          far=Math.exp(-((r/380)**2))+.3*Math.exp(-((r/850)**2)),   // tighter e-folding (was 650/1350) so intensity actually fades out within the smaller RMAX above, instead of staying strong out to the old edge
          az=.35+1*Math.max(0,Math.cos(ang-biasAng)),
          w=far*az*taper;
        if(w<.05)continue;
        const clon=s.lon+r*Math.cos(ang)/(111*cl),clat=s.lat+r*Math.sin(ang)/111,
          jr=Math.random()*4,ja=Math.random()*2*Math.PI;   // small organic jitter only - width and curl are now shaped explicitly, not left to randomness
        rAdd(clon+jr*Math.cos(ja)/(111*cl),clat+jr*Math.sin(ja)/111,peak*w*(.75+.5*Math.random()))}}}}}
// Ride: each parcel is advected one hour mostly by the surface wind, with a smaller 850 hPa contribution for a
// touch of steering depth (both already carry the storm's spiral inflow/tangential/outflow, so parcels naturally
// curve inward near the core and stream outward along the bands, exactly tracking the evolving flow). It dries
// out faster over low-humidity air, and is removed once it ages out, fades below a visible threshold, drifts off
// the grid, or - closing the cycle - spirals back within a small radius of ANY storm's centre.
// A single one-hour Euler step is far too coarse near the eyewall: at r~rm an intense storm's tangential wind
// (~s.v) gives an orbital period of only 1.5-2 h, so one hourly jump covers upwards of half an orbit and flings
// the parcel outward along the tangent instead of letting it curve - eyewall rain self-ejects into the far field
// after a step or two instead of building up a ring. Sub-stepping (resampling the wind at each intermediate
// position) keeps each leg's turning angle small so tight inner-core orbits stay closed loops; outer parcels in
// weaker, slower-turning flow are unaffected since their per-substep displacement is proportionally tiny anyway.
const ADV_SUB=6,ADV_DT=3.6/ADV_SUB;
function advectRain(){const{u,v,lv,q}=F,u8=lv[2],v8=lv[3];let i=0;
  while(i<rPN){let lon=rLon[i],lat=rLat[i];
    for(let sub=0;sub<ADV_SUB;sub++){const cl=Math.cos(lat*R),
        uu=.9*sampG(u,lon,lat)+.1*sampG(u8,lon,lat),vv=.9*sampG(v,lon,lat)+.1*sampG(v8,lon,lat);
      lon+=uu*ADV_DT/(111*cl);lat+=vv*ADV_DT/111}
    const nl=lon,na=lat;
    rLon[i]=nl;rLat[i]=na;rAge[i]++;
    rAmt[i]*=.94-.06*clamp(1-sampG(q,nl,na),0,1);      // dry air along the way -> the parcel evaporates faster
    let dead=rAge[i]>30||rAmt[i]<.6||nl<C.lon0||nl>C.lon1||na<C.lat0||na>C.lat1;
    // Skip the reabsorption check on a parcel's first tick: without this, a parcel born right at the center
    // (from spawnRain's new centerFloor term) got reabsorbed here before ever reaching the deposit loop below,
    // so the middle never actually got painted no matter how much rain was spawned there.
    if(!dead&&rAge[i]>1)for(const s of S.storms){if(Math.hypot((nl-s.lon)*111*Math.cos(s.lat*R),(na-s.lat)*111)<Math.min(18,.14*s.rm)){dead=true;break}}
    if(dead){rPN--;rLon[i]=rLon[rPN];rLat[i]=rLat[rPN];rAmt[i]=rAmt[rPN];rAge[i]=rAge[rPN]}else i++}}
// Fill-in: a parcel swarm, however dense, still leaves gaps between individual footprints at any single instant.
// A cheap local blur (two box-blur passes, restricted to each storm's own footprint so the open-ocean background
// keeps its texture) closes those gaps into a continuous shield/band look without hiding the underlying motion -
// next frame's parcel positions still drive where the smoothed field ends up.
function smoothPr(pr){for(const s of S.storms){const cl=Math.cos(s.lat*R),RCH=8.5,dl=RCH/cl,   // reach just needs to cover where spawnRain seeds (now out to ~850 km, was 1350); NOT relied on any more to control band width - that's shaped explicitly in spawnRain
    i0=Math.max(1,Math.floor((s.lon-dl-C.lon0)/D)),i1=Math.min(GW-2,Math.ceil((s.lon+dl-C.lon0)/D)),
    j0=Math.max(1,Math.floor((C.lat1-s.lat-RCH)/D)),j1=Math.min(GH-2,Math.ceil((C.lat1-s.lat+RCH)/D));
  for(let pass=0;pass<1;pass++){   // one pass is enough now that the raw deposits are already the right size/shape - a second pass was extra spreading the band didn't need
    for(let j=j0;j<=j1;j++)for(let i=i0;i<=i1;i++){const k=j*GW+i;PRT[k]=(pr[k]*2+pr[k-1]+pr[k+1])*.25}
    for(let j=j0;j<=j1;j++)for(let i=i0;i<=i1;i++){const k=j*GW+i;pr[k]=(PRT[k]*2+PRT[k-GW]+PRT[k+GW])*.25}}}}
// Carve the dry eye back out AFTER the parcel deposits and box-blur. Deposits + smoothPr already leave the eye
// mostly rain-free at birth (spawnRain suppresses it via eyeHole), but for small, intense storms the eye radius
// (~0.3*rm) can be well under one grid cell (0.25 deg, ~27 km) or under the blur radius, so the surrounding
// eyewall's max-intensity ring gets smeared straight across the middle - the center reads as the MOST intense
// spot instead of the driest one. This runs last, after smoothPr, so nothing re-fills what it clears, and the
// radius has a floor so the hole stays visually resolvable even at coarse grid spacing.
function carveEyes(pr){for(const s of S.storms){if(s.decay)continue;
  const eyeStrength=clamp((s.v-18)/22,0,1);if(eyeStrength<=0)continue;
  const cl=Math.cos(s.lat*R),rad=Math.max(D*111*.55,.24*s.rm),   // km; floored to just over half a grid cell
    dlon=(rad/111)/cl,dlat=rad/111,
    i0=Math.max(0,Math.floor((s.lon-dlon-C.lon0)/D)),i1=Math.min(GW-1,Math.ceil((s.lon+dlon-C.lon0)/D)),
    j0=Math.max(0,Math.floor((C.lat1-s.lat-dlat)/D)),j1=Math.min(GH-1,Math.ceil((C.lat1-s.lat+dlat)/D));
  for(let j=j0;j<=j1;j++){const lat=C.lat1-(j+.5)*D;
    for(let i=i0;i<=i1;i++){const lon=C.lon0+(i+.5)*D,km=Math.hypot((lon-s.lon)*111*cl,(lat-s.lat)*111);
      if(km>rad)continue;const t=km/rad,supp=1-eyeStrength*.92*(1-t*t);pr[j*GW+i]*=Math.max(.06,supp)}}}}
// Deterministic center rain: painted directly onto the grid (not through the random parcel swarm, which has no
// reliable way to keep a dense, stable feature confined to a handful of cells - spawn variance and the
// reabsorption/advection cycle either starved the center or let it blow up unpredictably). Iterates the same
// small per-storm bounding box carveEyes uses, so cost stays proportional to storm count, not grid size. Each
// cell's amount is modulated by the ACTUAL local wind vector (already fully updated by fields() for this step) -
// both its magnitude (calmer center -> lighter rain, consistent with the storm's own dynamics rather than an
// arbitrary radial number) and its direction (feeds a low-order angular wobble so neighboring cells vary with
// the real flow instead of every cell in the disc getting an identical, flat value) - so the patch reads as
// textured convection tied to the storm's own wind field, not a uniform painted circle.
function centerRain(pr){const{u,v}=F;
  for(const s of S.storms){if(s.decay)continue;
    const eyeStrength=clamp((s.v-18)/22,0,1),
      rad=16+10*eyeStrength,                    // km - grows 16->26 with intensity (bigger core for stronger storms), still kept well inside rw so it never reaches the eyewall ring
      cl=Math.cos(s.lat*R),dlon=(rad/111)/cl,dlat=rad/111,
      i0=Math.max(0,Math.floor((s.lon-dlon-C.lon0)/D)),i1=Math.min(GW-1,Math.ceil((s.lon+dlon-C.lon0)/D)),
      j0=Math.max(0,Math.floor((C.lat1-s.lat-dlat)/D)),j1=Math.min(GH-1,Math.ceil((C.lat1-s.lat+dlat)/D));
    for(let j=j0;j<=j1;j++){const lat=C.lat1-(j+.5)*D;
      for(let i=i0;i<=i1;i++){const k=j*GW+i,lon=C.lon0+(i+.5)*D,
          dx=(lon-s.lon)*111*cl,dy=(lat-s.lat)*111,km=Math.hypot(dx,dy);
        if(km>rad)continue;
        const t=km/rad,base=Math.exp(-((t*2.1)**2)),               // tight Gaussian, ~zero by the time it reaches rad
          spd=Math.hypot(u[k],v[k]),mag=.5+.5*Math.min(1,spd/25),  // scales with this cell's own wind speed
          ang=Math.atan2(v[k],u[k]),wobble=.75+.25*Math.cos(3*ang+2*Math.atan2(dy,dx)); // cell-to-cell texture from local flow direction
        pr[k]+=42*base*mag*wobble}}}}
function precip(){const{pr}=F;
  // Precipitation is now purely the storm's own convection: no terrain/orographic term at all, even inside a
  // storm's own footprint. Gating it to nearby-storm-only (previous version) still let the ambient wind hitting
  // terrain show up as rain whenever a storm happened to be nearby, which still read as unwanted "environment"
  // rain. Now the ENTIRE precip field comes from the Lagrangian rain swarm (spawnRain/advectRain/smoothPr) and
  // the explicit center patch (centerRain) below - both are storm-only constructs with no terrain input.
  pr.fill(0);
  spawnRain();advectRain();
  for(let i=0;i<rPN;i++){const fx=(rLon[i]-C.lon0)/D-.5,fy=(C.lat1-rLat[i])/D-.5,ci=Math.round(fx),cj=Math.round(fy);
    for(let dj=-RP_R;dj<=RP_R;dj++){const j=cj+dj;if(j<0||j>=GH)continue;
      for(let di=-RP_R;di<=RP_R;di++){const ii=ci+di;if(ii<0||ii>=GW)continue;
        pr[j*GW+ii]+=rAmt[i]*Math.exp(-(di*di+dj*dj)/RP_S)}}}
  smoothPr(pr);centerRain(pr)/*,carveEyes(pr) - disabled for now, punched a hole bigger than the storm's real eye; revisit later*/}
S.step=function(){S.t++;if(S.hm&&!S.hOff&&!hindOn()){S.hOff=true;S.nextSpawn=S.t+gap(true)}env();if(hindOn()){loadAna();if(S.wx==='replay')ZP.set(ZA);else dyn();hindSpawn()}else dyn();composeCG();if(S.t%24==0)setSST(S.E.sn);
  if(S.t>=S.nextSpawn){if(S.storms.length<6){const f=S.spawn();if(f&&Math.random()<.3)S.spawn(f)}   /* ~30% of genesis events get a companion 900-1500 km away */ S.nextSpawn=S.t+gap()}
  for(let k=0;k<N;k++)F.sa[k]*=.998;   // cold-wake recovery, ~21 d e-folding (real re-stratification is roughly 1-3 weeks)
  for(const s of S.storms.slice())upd(s);while(merge());fields();humidity();precip()};
})();