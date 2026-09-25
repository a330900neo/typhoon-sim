(function(){
const C=TS.cfg,S=TS.sim,F=S.F,GW=C.gw,GH=C.gh,D=(C.lon1-C.lon0)/GW,KX=.9,$=id=>document.getElementById(id),cv=$('map'),g=cv.getContext('2d');
const clamp=(x,a,b)=>x<a?a:x>b?b:x;let W=0,H=0,dpr=1;const V={lon:135,lat:22,s:10},R=Math.PI/180;let fcData=null,fcT=-999;
const minS=()=>Math.max(W/(90*KX),H/55);
function fixV(){V.s=clamp(V.s,minS(),600);const hw=W/(2*V.s*KX),hh=H/(2*V.s);V.lon=clamp(V.lon,C.lon0+hw,C.lon1-hw);V.lat=clamp(V.lat,C.lat0+hh,C.lat1-hh)}
function resize(){dpr=Math.min(devicePixelRatio||1,2);W=innerWidth;H=innerHeight;cv.width=W*dpr;cv.height=H*dpr;g.setTransform(dpr,0,0,dpr,0,0);fixV()}
const X=lon=>(lon-V.lon)*V.s*KX+W/2,Y=lat=>(V.lat-lat)*V.s+H/2,LON=x=>V.lon+(x-W/2)/(V.s*KX),LAT=y=>V.lat-(y-H/2)/V.s;
function zoomAt(x,y,f){const lo=LON(x),la=LAT(y);V.s=clamp(V.s*f,minS(),600);V.lon=lo-(x-W/2)/(V.s*KX);V.lat=la+(y-H/2)/V.s;fixV()}
const ptr=new Map(),pinch=()=>{if(ptr.size<2)return null;const[a,b]=[...ptr.values()];return{d:Math.hypot(a[0]-b[0],a[1]-b[1]),x:(a[0]+b[0])/2,y:(a[1]+b[1])/2}};let pd=null;
// hold D + click: drop a manual TD invest at the clicked point (see S.spawnAt in sim.js)
let keyD=false,downPos=null;
const setKeyD=v=>{keyD=v;cv.classList.toggle('dspawn',v)};
addEventListener('keydown',e=>{if((e.key==='d'||e.key==='D')&&!e.repeat&&!/input|select|textarea/i.test(e.target.tagName))setKeyD(true)});
addEventListener('keyup',e=>{if(e.key==='d'||e.key==='D')setKeyD(false)});
addEventListener('blur',()=>setKeyD(false));
let toastT=0;function toast(msg,bad){const t=$('toast');t.textContent=msg;t.classList.toggle('bad',!!bad);t.classList.add('show');
  clearTimeout(toastT);toastT=setTimeout(()=>t.classList.remove('show'),1500)}
cv.onpointerdown=e=>{cv.setPointerCapture(e.pointerId);downPos=[e.clientX,e.clientY];
  if(keyD)return;   // D-held click is a spawn gesture only, never a pan/zoom drag
  ptr.set(e.pointerId,[e.clientX,e.clientY]);pd=pinch()};
cv.onpointermove=e=>{hover(e);const o=ptr.get(e.pointerId);if(!o)return;
  if(ptr.size==1){V.lon-=(e.clientX-o[0])/(V.s*KX);V.lat+=(e.clientY-o[1])/V.s;fixV()}
  ptr.set(e.pointerId,[e.clientX,e.clientY]);const p=pinch();if(p&&pd){zoomAt(p.x,p.y,p.d/pd.d);V.lon-=(p.x-pd.x)/(V.s*KX);V.lat+=(p.y-pd.y)/V.s;fixV()}pd=p};
cv.onpointerup=cv.onpointercancel=e=>{
  if(keyD&&downPos&&ptr.size===0&&Math.hypot(e.clientX-downPos[0],e.clientY-downPos[1])<6){
    const r=S.spawnAt(LON(e.clientX),LAT(e.clientY));
    if(r&&r.ok){toast(r.storm.name+' spawned');dirty=true;ui()}else toast("Can't spawn — "+((r&&r.why)||'try open ocean'),true)}
  ptr.delete(e.pointerId);pd=null;downPos=null;
  try{cv.releasePointerCapture(e.pointerId)}catch(_){}};
cv.addEventListener('wheel',e=>{e.preventDefault();zoomAt(e.clientX,e.clientY,Math.exp(-e.deltaY*.0015))},{passive:false});
$('zin').onclick=()=>zoomAt(W/2,H/2,1.5);$('zout').onclick=()=>zoomAt(W/2,H/2,1/1.5);$('zfit').onclick=()=>{V.s=minS();V.lon=135;V.lat=22;fixV()};
// colormaps
const build=st=>{const a=new Uint8Array(768);for(let k=0;k<256;k++){const f=k/255;let n=1;while(n<st.length-1&&st[n][0]<f)n++;
  const A=st[n-1],B=st[n],t=clamp((f-A[0])/(B[0]-A[0]),0,1);for(let c=0;c<3;c++)a[k*3+c]=A[1][c]+(B[1][c]-A[1][c])*t}return a};
const L={wind:{lo:0,hi:90,u:'m/s',st:[[0,[8,16,60]],[.19,[0,170,225]],[.27,[35,205,95]],[.37,[255,235,0]],[.49,[255,120,0]],[.6,[225,15,15]],[.8,[190,0,120]],[1,[255,255,255]]]},
 pressure:{lo:900,hi:1025,u:'hPa',st:[[0,[200,60,200]],[.2,[220,40,60]],[.5,[240,200,60]],[.8,[80,200,120]],[1,[30,60,140]]]},
 humidity:{lo:0,hi:1,u:'rel. humidity',st:[[0,[190,150,90]],[.5,[110,190,190]],[.8,[40,120,200]],[1,[10,40,150]]]},
 sst:{lo:15,hi:31,u:'°C',st:[[0,[40,60,160]],[.5,[60,180,200]],[.7,[240,220,80]],[.87,[240,130,40]],[1,[200,30,40]]]},
};
for(const k in L)L[k].lut=build(L[k].st);
const CAT=[[105,'Violent Typhoon','#ff2bd6'],[85,'Very Strong Typhoon','#ff3b3b'],[64,'Typhoon','#ff8a2b'],[48,'Severe Tropical Storm','#ffd23b'],[34,'Tropical Storm','#7fe06a'],[0,'Tropical Depression','#5ec8ff']];
const cat=v=>CAT.find(c=>v*1.944>=c[0]);
let lvU=F.u,lvV=F.v,base,layer='wind',speed=6,acc=0,last=performance.now(),fr=0,dirty=true,isoT=0,iso=[],SM;
const announced=new WeakSet();
function announceEvents(){for(const s of S.storms.concat(S.hist))for(const e of s.ev||[]){if(announced.has(e))continue;announced.add(e);
  const name=s.name||'Storm';
  if(e.t==='i')toast(`${name} ${e.up?'intensifying':'weakening'}: ${CAB[e.c]}`);
  else if(e.t==='in')toast(`${name} entered land`,true);
  else if(e.t==='out')toast(`${name} returned to sea`)} }
const ov=document.createElement('canvas');ov.width=GW;ov.height=GH;const og=ov.getContext('2d'),oi=og.createImageData(GW,GH);
function legend(){const l=L[layer],c=$('legend').getContext('2d');c.clearRect(0,0,220,10);if(!l){$('legtxt').textContent='';return}
  for(let x=0;x<220;x++){const k=((x/219*255)|0)*3;c.fillStyle=`rgb(${l.lut[k]},${l.lut[k+1]},${l.lut[k+2]})`;c.fillRect(x,0,1,10)}$('legtxt').textContent=`${l.lo} → ${l.hi} ${l.u}`}
function updOverlay(){const l=L[layer];if(!l)return;const d=oi.data;
  for(let k=0;k<GW*GH;k++){const val=layer==='wind'?Math.hypot(lvU[k],lvV[k]):layer==='sst'?F.sst0[k]+F.sa[k]:layer==='pressure'?F.p[k]:F.q[k],
    i=clamp(((val-l.lo)/(l.hi-l.lo)*255)|0,0,255)*3;d[k*4]=l.lut[i];d[k*4+1]=l.lut[i+1];d[k*4+2]=l.lut[i+2];
    d[k*4+3]=layer==='sst'&&F.el[k]>0?0:150}
  og.putImageData(oi,0,0)}
// isobars via marching squares (4 hPa), segments in grid coords
function calcIso(){const p=F.p;
  // Light 5-point smoothing purely for contouring: marching squares point-samples the field, so a single noisy
  // grid cell (leftover computational noise from the barotropic solver) can otherwise draw its own tiny, unlabeled
  // closed loop that has nothing to do with the actual synoptic pattern.
  if(!SM)SM=new Float32Array(GW*GH);SM.set(p);
  for(let j=1;j<GH-1;j++)for(let i=1;i<GW-1;i++){const k=j*GW+i;SM[k]=(p[k]*4+p[k-1]+p[k+1]+p[k-GW]+p[k+GW])*.125}
  const s=[],st=Math.max(1,Math.round(.25/D)),G=GW*st;   // isobars on a 0.25 deg subsample (smooth field, keeps it cheap at 0.125 deg)
  for(let j=0;j<GH-st;j+=st)for(let i=0;i<GW-st;i+=st){const k=j*GW+i,a=SM[k],b=SM[k+st],c=SM[k+G+st],d=SM[k+G],mn=Math.min(a,b,c,d),mx=Math.max(a,b,c,d);
    for(let lv=Math.ceil(mn/4)*4;lv<=mx;lv+=4){const q=[],e=(v1,v2,x1,y1,x2,y2)=>{if((v1<lv)!==(v2<lv)){const t=(lv-v1)/(v2-v1);q.push(x1+(x2-x1)*t,y1+(y2-y1)*t)}};
      e(a,b,i,j,i+st,j);e(b,c,i+st,j,i+st,j+st);e(d,c,i,j+st,i+st,j+st);e(a,d,i,j,i,j+st);
      if(q.length>=4)s.push(lv,q[0],q[1],q[2],q[3]);if(q.length==8)s.push(lv,q[4],q[5],q[6],q[7])}}
  iso=s}
const gx=x=>X(C.lon0+(x+.5)*D),gy=y=>Y(C.lat1-(y+.5)*D);
function drawIso(){for(let pass=0;pass<2;pass++){g.beginPath();for(let n=0;n<iso.length;n+=5){if((iso[n]%20===0)!==(pass===1))continue;
    g.moveTo(gx(iso[n+1]),gy(iso[n+2]));g.lineTo(gx(iso[n+3]),gy(iso[n+4]))}
    g.strokeStyle=pass?'rgba(255,255,255,.9)':'rgba(255,255,255,.4)';g.lineWidth=pass?1.6:.8;g.stroke()}
  g.font='bold 11px sans-serif';g.strokeStyle='#000';g.fillStyle='#fff';g.lineWidth=3;
  for(let n=0;n<iso.length;n+=5)if(iso[n]%20===0&&(n/5)%70===0){const x=gx(iso[n+1]),y=gy(iso[n+2]);if(x>0&&y>0&&x<W&&y<H){g.strokeText(iso[n],x,y);g.fillText(iso[n],x,y)}}}
// wind sampling (bilinear)
let su=0,sv=0;function samp(lon,lat){const x=(lon-C.lon0)/D-.5,y=(C.lat1-lat)/D-.5;if(x<0||y<0||x>=GW-1||y>=GH-1){su=sv=0;return false}
  const x0=x|0,y0=y|0,fx=x-x0,fy=y-y0,a=y0*GW+x0,U=lvU,Vv=lvV;
  su=(U[a]*(1-fx)+U[a+1]*fx)*(1-fy)+(U[a+GW]*(1-fx)+U[a+GW+1]*fx)*fy;sv=(Vv[a]*(1-fx)+Vv[a+1]*fx)*(1-fy)+(Vv[a+GW]*(1-fx)+Vv[a+GW+1]*fx)*fy;return true}
const wcol=w=>`hsl(${Math.max(0,200-w*4)},90%,75%)`;
const NP=3500,TL=7,pl=new Float32Array(NP),pa=new Float32Array(NP),pg=new Uint16Array(NP),pt=new Float32Array(NP*TL*2);
function seed(i){let t=0;do{const x=Math.random()*W,y=Math.random()*H;pl[i]=LON(x);pa[i]=LAT(y)}while(!samp(pl[i],pa[i])&&t++<5);
  pg[i]=(Math.random()*-70)|0;for(let m=0;m<TL;m++){pt[(i*TL+m)*2]=pl[i];pt[(i*TL+m)*2+1]=pa[i]}}
for(let i=0;i<NP;i++){seed(i);pg[i]=Math.random()*90|0}
function particles(){const bins=[[],[],[]],slot=fr%TL;
  for(let i=0;i<NP;i++){if(!samp(pl[i],pa[i])||++pg[i]>100||pl[i]<LON(0)||pl[i]>LON(W)||pa[i]>LAT(0)||pa[i]<LAT(H)){seed(i);continue}
    const w=Math.hypot(su,sv);pl[i]+=su*.16/(V.s*KX);pa[i]+=sv*.16/V.s;pt[(i*TL+slot)*2]=pl[i];pt[(i*TL+slot)*2+1]=pa[i];bins[w<10?0:w<25?1:2].push(i)}
  const cols=['rgba(255,255,255,.45)','rgba(255,240,150,.7)','rgba(255,150,90,.9)'];g.lineWidth=1.1;
  bins.forEach((b,bi)=>{g.beginPath();for(const i of b){for(let m=1;m<=TL;m++){const s=(slot+m)%TL,x=X(pt[(i*TL+s)*2]),y=Y(pt[(i*TL+s)*2+1]);m==1?g.moveTo(x,y):g.lineTo(x,y)}}g.strokeStyle=cols[bi];g.stroke()})}
function arrows(){const sp=36;g.lineWidth=1.4;g.lineCap='round';
  for(let y=sp/2;y<H;y+=sp)for(let x=sp/2;x<W;x+=sp){if(!samp(LON(x),LAT(y)))continue;const w=Math.hypot(su,sv);if(w<.3)continue;
    const len=Math.min(32,6+w*1.1),a=Math.atan2(-sv,su),c=Math.cos(a),s=Math.sin(a),x1=x-c*len/2,y1=y-s*len/2,x2=x+c*len/2,y2=y+s*len/2,h=Math.min(8,len*.4);
    g.strokeStyle=wcol(w);g.beginPath();g.moveTo(x1,y1);g.lineTo(x2,y2);g.lineTo(x2-h*Math.cos(a-.5),y2-h*Math.sin(a-.5));g.moveTo(x2,y2);g.lineTo(x2-h*Math.cos(a+.5),y2-h*Math.sin(a+.5));g.stroke()}}
// track nodes: circle = intensity category change (filled = strengthened, hollow = weakened), triangle = entered terrain, diamond = exited to sea
const catOf=k=>{let c=0;for(const t of[34,48,64,85,105])if(k>=t)c++;return c};
const CCOL=['#5ec8ff','#7fe06a','#ffd23b','#ff8a2b','#ff3b3b','#ff2bd6'],CAB=['TD','TS','STS','TY','VSTY','VTY'];
function nodes(){const lab=V.s>=14;g.font='bold 10px sans-serif';
  for(const s of S.hist.concat(S.storms)){if(!s.ev)continue;g.globalAlpha=S.hist.includes(s)?.55:(s.fade==null?1:s.fade);
    for(const e of s.ev){const x=X(e.lon),y=Y(e.lat);if(x<-20||y<-20||x>W+20||y>H+20)continue;g.beginPath();g.lineWidth=2;
      if(e.t==='i'){g.arc(x,y,5,0,7);g.fillStyle=e.up?CCOL[e.c]:'#101828';g.strokeStyle=e.up?'#fff':CCOL[e.c]}
      else if(e.t==='in'){g.moveTo(x,y-7.5);g.lineTo(x+7,y+5);g.lineTo(x-7,y+5);g.closePath();g.fillStyle='#ffb347';g.strokeStyle='#fff'}
      else{g.moveTo(x,y-7.5);g.lineTo(x+6.5,y);g.lineTo(x,y+7.5);g.lineTo(x-6.5,y);g.closePath();g.fillStyle='#56e0ff';g.strokeStyle='#fff'}
      g.fill();g.stroke();
      if(lab){const t=e.t==='i'?(e.up?'↑ ':'↓ ')+CAB[e.c]:e.t==='in'?`LANDFALL ${(e.v*1.944)|0} kt`:'exit to sea';g.lineWidth=3;g.strokeStyle='#000';g.fillStyle='#fff';g.strokeText(t,x+8,y-6);g.fillText(t,x+8,y-6)}}}
  g.globalAlpha=1}
function tracks(){if(S.hm&&TS.hind)for(const o of TS.hind.storms){if(!o.done)continue;g.save();g.setLineDash([5,4]);g.strokeStyle='rgba(255,255,255,.75)';g.lineWidth=1.6;g.beginPath();o.trk.forEach((p,i)=>i?g.lineTo(X(p[1]),Y(p[2])):g.moveTo(X(p[1]),Y(p[2])));g.stroke();g.restore()}
  for(const s of S.hist.concat(S.storms)){const isHist=S.hist.includes(s);g.globalAlpha=isHist?.35:(s.fade==null?1:s.fade);g.lineWidth=2.2;
  for(let i=1;i<s.trk.length;i++){g.strokeStyle=cat(s.trk[i][2])[2];g.beginPath();g.moveTo(X(s.trk[i-1][0]),Y(s.trk[i-1][1]));g.lineTo(X(s.trk[i][0]),Y(s.trk[i][1]));g.stroke()}}
  nodes();g.globalAlpha=1;g.font='bold 12px sans-serif';
  for(const s of S.storms){g.globalAlpha=(s.named?1:.55)*(s.fade==null?1:s.fade);const x=X(s.lon),y=Y(s.lat),c=cat(s.v)[2],r=8+s.v/8,a=fr*clamp(.02+s.v*.0016,.02,.13);g.fillStyle=c+'44';g.strokeStyle=c;g.lineWidth=2.5;g.beginPath();g.arc(x,y,r,0,7);g.fill();
    for(let k=0;k<2;k++){g.beginPath();g.arc(x,y,r*.65,a+k*Math.PI,a+k*Math.PI+2.2);g.stroke()}
    if(s.sv){g.strokeStyle='#fff';g.lineWidth=2;g.beginPath();g.moveTo(x,y);g.lineTo(x+s.sv[0]*7,y-s.sv[1]*7);g.stroke()}
    const kt=s.v*1.944,sp=s.spd*1.944,hd=(Math.atan2(s.mu,s.mv)*180/Math.PI+360)%360,dir=['N','NE','E','SE','S','SW','W','NW'][Math.round(hd/45)%8],
      L1=`${s.name} · ${CAB[catOf(kt)]}`,L2=`${kt|0} kt · ${s.pmin|0} hPa`,L3=`→${dir} ${sp.toFixed(0)} kt`,L4=s.land?'⚠ LANDFALL / over land':null,tx=x+r+4;let ty=y-r;
    g.fillStyle='#fff';g.strokeStyle='#000';g.lineWidth=3;
    for(const[t,col]of[[L1,'#fff'],[L2,c],[L3,'#cfe8ff'],[L4,'#ffb347']]){if(!t)continue;g.fillStyle=col;g.strokeText(t,tx,ty);g.fillText(t,tx,ty);ty+=13}}g.globalAlpha=1}
function drawForecast(){for(const f of fcData){const trk=f.trk;if(trk.length<2)continue;
    const L=[],Rt=[];
    for(let i=0;i<trk.length;i++){const[lo,la,,km]=trk[i],p=trk[Math.max(0,i-1)],n=trk[Math.min(trk.length-1,i+1)],
      dx=(n[0]-p[0])*Math.cos(la*R)*111,dy=(n[1]-p[1])*111,ang=Math.atan2(dy,dx)+Math.PI/2,
      dLon=km*Math.cos(ang)/(111*Math.cos(la*R)),dLat=km*Math.sin(ang)/111;
      L.push([X(lo+dLon),Y(la+dLat)]);Rt.push([X(lo-dLon),Y(la-dLat)])}
    g.beginPath();L.forEach((p,i)=>i?g.lineTo(p[0],p[1]):g.moveTo(p[0],p[1]));for(let i=Rt.length-1;i>=0;i--)g.lineTo(Rt[i][0],Rt[i][1]);g.closePath();
    g.fillStyle='rgba(255,255,255,.10)';g.strokeStyle='rgba(255,255,255,.3)';g.lineWidth=1;g.fill();g.stroke();
    g.save();g.setLineDash([6,4]);g.strokeStyle='rgba(255,255,255,.85)';g.lineWidth=2;g.beginPath();
    trk.forEach((p,i)=>i?g.lineTo(X(p[0]),Y(p[1])):g.moveTo(X(p[0]),Y(p[1])));g.stroke();g.restore();
    g.font='bold 10px sans-serif';
    for(let i=0;i<trk.length;i++){const[lo,la,v]=trk[i],x=X(lo),y=Y(la),c=cat(v)[2],h=i*3;
      g.beginPath();g.arc(x,y,4,0,7);g.fillStyle=c;g.strokeStyle='#fff';g.lineWidth=1.5;g.fill();g.stroke();
      if(h%24===0&&h>0){g.lineWidth=3;g.strokeStyle='#000';g.fillStyle='#fff';g.strokeText('+'+h+'h',x+6,y-6);g.fillText('+'+h+'h',x+6,y-6)}}}}
function ui(){const d=new Date(S.t0+S.t*36e5);$('clock').textContent=`${d.toISOString().slice(0,16).replace('T',' ')}Z  (T+${S.t}h) · ${S.E.enso>.3?'El Niño':S.E.enso<-.3?'La Niña':'ENSO-neutral'} · ${S.E.mjo>.25?'active spell':S.E.mjo<-.25?'quiet spell':'near-avg activity'} · ${S.mode()}`;
  let h=S.storms.map(s=>{const c=cat(s.v);return `<div class="st" style="border-color:${c[2]}"><b>${s.name}</b> — ${c[1]}<br>${(s.v*1.944)|0} kt (${s.v.toFixed(0)} m/s) · ${s.pmin.toFixed(0)} hPa · RMW ${s.rm.toFixed(0)} km<br>${s.lat.toFixed(1)}°N ${s.lon.toFixed(1)}°E · ${(s.spd*1.944).toFixed(0)} kt${s.land?' · <b style="color:#ffb347">LANDFALL</b>':''}${s.lf?' · '+s.lf+' landfall'+(s.lf>1?'s':''):''}${s.err!=null?' · obs err '+(s.err|0)+' km':''}</div>`}).join('');
  h+=S.hist.slice(-5).reverse().map(s=>`<div class="st" style="border-color:#556;opacity:.7">${s.name} ${s.by?'absorbed by '+s.by:'dissipated'} · peak ${(s.max*1.944)|0} kt</div>`).join('');
  $('storms').innerHTML=h||'<em>None yet — press Start Sim</em>'}
function hover(e){if(!TS.elev0)return;const lon=LON(e.clientX),lat=LAT(e.clientY),k=S.cellIdx(lon,lat);if(k<0)return;
  const ix=Math.min(C.bw-1,Math.max(0,((lon-C.lon0)*C.bw/90)|0)),iy=Math.min(C.bh-1,Math.max(0,((C.lat1-lat)*C.bh/55)|0)),el=TS.elev0[iy*C.bw+ix],w=Math.hypot(F.u[k],F.v[k]);
  $('hover').innerHTML=`${lat.toFixed(2)}°N ${lon.toFixed(2)}°E<br>${el>0?'Elev':'Depth'}: ${Math.abs(el)|0} m<br>Surface wind: ${w.toFixed(1)} m/s (${(w*1.944)|0} kt)<br>${[['850',2],['500',4],['200',6]].map(([n,m])=>{const u=F.lv[m][k],v=F.lv[m+1][k];return `${n} hPa: ${Math.hypot(u,v).toFixed(0)} m/s from ${((270-Math.atan2(v,u)*180/Math.PI)%360+360)%360|0}°`}).join('<br>')}<br>Pressure: ${F.p[k].toFixed(1)} hPa<br>Humidity: ${(F.q[k]*100)|0}%<br>SST: ${(F.sst0[k]+F.sa[k]).toFixed(1)} °C`}
function loop(now){const dt=Math.min(.1,(now-last)/1000);last=now;
  if(S.running){acc+=dt*speed;const t0=performance.now();let n=0;
    while(acc>=1&&n<48&&performance.now()-t0<14){
      S.step();
      acc--;n++;dirty=true}
    if(acc>3)acc=3}
  announceEvents();if(dirty){updOverlay();if(now-isoT>250){calcIso();isoT=now;dirty=false}}
  g.drawImage(base,0,0,C.bw,C.bh,X(C.lon0),Y(C.lat1),90*V.s*KX,55*V.s);
  if(L[layer]){g.imageSmoothingEnabled=true;g.drawImage(ov,0,0,GW,GH,X(C.lon0),Y(C.lat1),90*V.s*KX,55*V.s)}
  if($('iso').checked)drawIso();const wm=$('wmode').value;if(wm==='particles')particles();else if(wm==='arrows')arrows();
  tracks();if($('fcst').checked){if(!fcData||S.t-fcT>=3){fcData=S.forecast(120);fcT=S.t}drawForecast()}
  if(++fr%8==0)ui();requestAnimationFrame(loop)}
$('start').onclick=()=>{S.running=!S.running;$('start').textContent=S.running?'⏸ Pause':'▶ Start Sim'};
$('reset').onclick=()=>{S.reset();$('start').textContent='▶ Start Sim';dirty=true;ui();fcData=null;fcT=-999};
$('speed').onchange=e=>speed=+e.target.value;$('layer').onchange=e=>{layer=e.target.value;legend();dirty=true};
$('lvl').onchange=e=>{const m={sfc:0,850:2,500:4,200:6}[e.target.value];lvU=F.lv[m];lvV=F.lv[m+1];dirty=true};
$('date').onchange=e=>{let t=Date.parse(e.target.value+'T00:00:00Z');if(isNaN(t))return;S.t0=t;$('reset').onclick()};
$('wx').onchange=e=>{S.wx=e.target.value;$('date').onchange({target:$('date')})};
$('ptoggle').onclick=()=>$('panel').classList.toggle('hide');
$('ttoggle').onclick=()=>{document.body.classList.add('top-hidden');$('top').classList.add('hide')};
$('tshow').onclick=()=>{document.body.classList.remove('top-hidden');$('top').classList.remove('hide')};addEventListener('resize',resize);
S.t0=Date.parse($('date').value+'T00:00:00Z');resize();V.s=minS();fixV();
TS.loadTerrain(p=>$('loading').textContent=`Loading terrain… ${(p*100)|0}%`).then(async el=>{await TS.loadClim('js/clim.bin');await TS.loadStats('js/stats.json');await TS.loadHind('js/hind.bin','js/hind.json');if(TS.hind)await TS.loadHindSST('js/hind_sst.bin');$('wx').disabled=!TS.hind;document.querySelector('.note').insertAdjacentHTML('beforeend','<br>Mean flow: '+(TS.clim?'reanalysis climatology':'analytic (run tools/make_clim.py)')+'; anomalies: barotropic model.'+(TS.stats?' Genesis from IBTrACS.':'')+(TS.hind?' Hindcast window: '+TS.hind.start+' + '+TS.hind.days+' d (dates outside it fall back to climatology + stochastic genesis).':' No hindcast loaded — Weather stays on Synthetic (see console for why; regenerate with tools/make_hindcast.py).')+' SST: '+(TS.hind&&TS.hind.sst&&TS.hindSST?'observed monthly means (ERSSTv5) while in real weather mode':'analytic climatology')+'.');
  const RX=C.bw/GW,ge=new Float32Array(GW*GH);for(let j=0;j<GH;j++)for(let i=0;i<GW;i++)ge[j*GW+i]=el[((j*RX+RX/2)|0)*C.bw+((i*RX+RX/2)|0)];
  base=TS.buildBase(el);S.init(ge);legend();ui();
  if(TS.tilesOk)$('loading').style.display='none';else $('loading').textContent='Could not load elevation tiles (offline?). Reload to retry.';
  $('start').disabled=$('reset').disabled=false;requestAnimationFrame(loop)});
})();