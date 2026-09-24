// Loads real elevation/bathymetry (Terrarium PNG tiles) and builds the base map.
window.TS={cfg:{lon0:90,lon1:180,lat0:-5,lat1:50,gw:720,gh:440,bw:3600,bh:2200,z:6}};
TS.loadTerrain=async function(onProg){
  // Tiles are decoded one at a time straight into the elevation array (no giant mosaic canvas), so 3600x2200 stays light on memory.
  const c=TS.cfg,n=2**c.z,mx=l=>(l+180)/360*n,
   my=a=>{const r=a*Math.PI/180;return(1-Math.log(Math.tan(r)+1/Math.cos(r))/Math.PI)/2*n};
  const x0=Math.floor(mx(c.lon0)),x1=Math.floor(mx(c.lon1-1e-6)),y0=Math.floor(my(c.lat1)),y1=Math.floor(my(c.lat0));
  const el=new Float32Array(c.bw*c.bh).fill(-2000),cx=new Float64Array(c.bw),ry=new Float64Array(c.bh),cr={},rr={};
  for(let i=0;i<c.bw;i++){cx[i]=mx(c.lon0+(i+.5)/c.bw*(c.lon1-c.lon0));const t=Math.floor(cx[i]);cr[t]?cr[t][1]=i:cr[t]=[i,i]}
  for(let j=0;j<c.bh;j++){ry[j]=my(c.lat1-(j+.5)/c.bh*(c.lat1-c.lat0));const t=Math.floor(ry[j]);rr[t]?rr[t][1]=j:rr[t]=[j,j]}
  const tc=document.createElement('canvas');tc.width=tc.height=256;const g=tc.getContext('2d',{willReadFrequently:true});
  const jobs=[];for(let x=x0;x<=x1;x++)for(let y=y0;y<=y1;y++)jobs.push([x,y]);
  const total=jobs.length;let done=0,ok=0;
  const one=([x,y])=>new Promise(res=>{const im=new Image();im.crossOrigin='anonymous';
    im.onload=()=>{g.clearRect(0,0,256,256);g.drawImage(im,0,0);const d=g.getImageData(0,0,256,256).data,ci=cr[x],rj=rr[y];
      if(ci&&rj)for(let j=rj[0];j<=rj[1];j++){const py=Math.min(255,((ry[j]-y)*256)|0);
        for(let i=ci[0];i<=ci[1];i++){const px=Math.min(255,((cx[i]-x)*256)|0),k=(py*256+px)*4;
          if(d[k+3])el[j*c.bw+i]=d[k]*256+d[k+1]+d[k+2]/256-32768}}
      ok++;onProg(++done/total);res()};
    im.onerror=()=>{onProg(++done/total);res()};
    im.src=`https://s3.amazonaws.com/elevation-tiles-prod/terrarium/${c.z}/${x}/${y}.png`});
  const worker=async()=>{while(jobs.length)await one(jobs.shift())};
  await Promise.all(Array.from({length:12},worker));
  TS.elev0=el;TS.tilesOk=ok;TS.tilesTotal=total;return el;
};
TS.buildBase=function(el){
  const c=TS.cfg,W=c.bw,H=c.bh,cv=document.createElement('canvas');cv.width=W;cv.height=H;
  const g=cv.getContext('2d'),im=g.createImageData(W,H),o=im.data,mix=(a,b,t)=>a+(b-a)*t;
  for(let j=0;j<H;j++)for(let i=0;i<W;i++){const k=j*W+i,e=el[k];let r,gg,b;
    if(e<=0){const t=Math.min(1,-e/6000);r=mix(40,8,t);gg=mix(95,22,t);b=mix(150,55,t)}
    else{const t=Math.min(1,e/4500),a=t<.4?t/.4:(t-.4)/.6;
      [r,gg,b]=t<.4?[mix(75,200,a),mix(135,180,a),mix(75,120,a)]:[mix(200,245,a),mix(180,245,a),mix(120,245,a)];
      const dz=(el[k-(i>0)]-el[k+(i<W-1)]+el[k+(j<H-1?W:0)]-el[k-(j>0?W:0)])*.5,s=Math.max(.65,Math.min(1.35,1+dz*.0009*(1800/W)));r*=s;gg*=s;b*=s}
    o[k*4]=r;o[k*4+1]=gg;o[k*4+2]=b;o[k*4+3]=255}
  g.putImageData(im,0,0);const q=W/1800;g.strokeStyle='rgba(255,255,255,.18)';g.fillStyle='rgba(255,255,255,.6)';g.font=(11*q)+'px sans-serif';g.lineWidth=q;
  for(let lo=90;lo<=180;lo+=10){const x=(lo-90)*W/90;g.beginPath();g.moveTo(x,0);g.lineTo(x,H);g.stroke();g.fillText(lo+'°E',x+3*q,H-4*q)}
  for(let la=0;la<=50;la+=10){const y=(50-la)*H/55;g.beginPath();g.moveTo(0,y);g.lineTo(W,y);g.stroke();g.fillText(la+'°N',3*q,y-3*q)}
  return cv;
};