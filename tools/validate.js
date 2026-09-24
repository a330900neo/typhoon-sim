// Headless calibration check: simulator vs IBTrACS statistics (js/stats.json).
// usage: node tools/validate.js [--months 6,7,8,9] [--reps 2] [--days 45] [--stats js/stats.json] [--clim js/clim.bin] [--land tools/land.bin]
const fs=require('fs'),path=require('path'),arg=(k,d)=>{const i=process.argv.indexOf('--'+k);return i>0?process.argv[i+1]:d},root=path.resolve(__dirname,'..');
global.TS={cfg:{lon0:90,lon1:180,lat0:-5,lat1:50,gw:720,gh:440,bw:3600,bh:2200,z:6}};
require(path.join(root,'js/sim.js'));const S=TS.sim,GW=TS.cfg.gw,GH=TS.cfg.gh,N=GW*GH;
const rd=(p,d)=>{const f=path.resolve(root,arg(p,d));return fs.existsSync(f)?fs.readFileSync(f):null};
const sj=rd('stats','js/stats.json'),cb=rd('clim','js/clim.bin'),lb=rd('land','tools/land.bin');
if(sj)TS.stats=JSON.parse(sj);if(cb)TS.clim=new Int16Array(cb.buffer.slice(cb.byteOffset,cb.byteOffset+cb.byteLength));
const el=new Float32Array(N);if(lb){const lw=lb.length===N?GW:360,lh=lb.length===N?GH:220;for(let j=0;j<GH;j++)for(let i=0;i<GW;i++)el[j*GW+i]=lb[((j*lh/GH)|0)*lw+((i*lw/GW)|0)]?300:-100}   // land.bin may be 360x220 (upsampled) or full-reselse console.log('(no land mask: landfall/terrain not exercised)');
const months=arg('months','0,1,2,3,4,5,6,7,8,9,10,11').split(',').map(Number),reps=+arg('reps',2),days=+arg('days',45);
const R=Math.PI/180,rows=[],all=[];
for(const m of months){const c=[];for(let r=0;r<reps;r++){S.t0=Date.UTC(2026,m,1);S.uid=0;S.init(el);for(let i=0;i<days*24;i++)S.step();
  for(const s of S.hist.concat(S.storms))if(s.max*1.944>=34&&s.trk.length>=12){c.push(s);all.push(s)}}
  let rc=0;for(const s of c){const t=s.trk,e=t[t.length-1];let mn=t[0];for(const p of t)if(p[0]<mn[0])mn=p;if(e[0]-mn[0]>6&&e[1]>mn[1]+4)rc++}
  rows.push({m,n:c.length,rate:c.length/(reps*days/30.4375),rec:c.length>=5?rc/c.length:null})}
const O=TS.stats,f=(x,d=1)=>x==null?'  - ':x.toFixed(d).padStart(5);
console.log('\nmonth | obs rate sim rate | obs recurve sim recurve   (rate = named storms/month; sim n in brackets)');
for(const r of rows)console.log(String(r.m+1).padStart(5),'|',f(O&&O.rate[r.m]),'   ',f(r.rate),'  |',f(O&&O.recurve_by_month[r.m]&&O.recurve_by_month[r.m],2),'     ',f(r.rec,2),' (n='+r.n+')');
const sp={};for(const s of all)for(let i=1;i<s.trk.length;i++){const a=s.trk[i-1],b=s.trk[i],la=(a[1]+b[1])/2,d=Math.hypot((b[0]-a[0])*111*Math.cos(la*R),(b[1]-a[1])*111)/3,k=Math.floor(la/5)*5;if(k>=5&&k<=45)(sp[k]=sp[k]||[]).push(d)}
console.log('\nlat band | obs km/h  sim km/h');for(const k of Object.keys(sp).sort((a,b)=>a-b))console.log(String(k).padStart(8),'|',f(O&&O.speed_by_lat[k]),'   ',f(sp[k].reduce((a,b)=>a+b)/sp[k].length));
const n=all.length||1,L=all.reduce((a,s)=>a+s.trk.length*3/24,0)/n;
console.log(`\nlandfall share: obs ${O?O.land_share:'-'}  sim ${(all.filter(s=>s.landed).length/n).toFixed(2)}`);
console.log(`reach 64 kt:    obs ${O?O.peak64_share:'-'}  sim ${(all.filter(s=>s.max*1.944>=64).length/n).toFixed(2)}`);
console.log(`mean life (d):  obs ${O?O.mean_life_days:'-'}  sim ${L.toFixed(2)}`);