/* External spur HPSTC and finite surface patch. mm, N, N*mm internally. */
(function(root){
  'use strict';
  const G=typeof module!=='undefined'&&module.exports?require('./gear'):root.GearGeom;
  const F=typeof module!=='undefined'&&module.exports?require('./fea'):root.FEASolver;
  function hpstc(p,{mateTeeth=40,mateShift=0,centerDistance=null}={}){
    const g=G.basic(p),mate=G.basic({...p,z:mateTeeth,x:mateShift}),aw=centerDistance===null?g.r+mate.r:centerDistance;
    if(!(aw>g.rb+mate.rb))throw Error('Invalid operating centre distance');
    const alphaW=Math.acos((g.rb+mate.rb)/aw),rw=g.rb/Math.cos(alphaW),rw2=mate.rb/Math.cos(alphaW),pb=Math.PI*g.m*Math.cos(g.alphaRad);
    const sA=rw2*Math.sin(alphaW)-Math.sqrt(mate.ra**2-mate.rb**2),sE=Math.sqrt(g.ra**2-g.rb**2)-rw*Math.sin(alphaW),ratio=(sE-sA)/pb;
    if(!(ratio>1&&ratio<2))throw Error('HPSTC requires 1 < contact ratio < 2');
    const s=sA+pb,R=Math.hypot(g.rb,rw*Math.sin(alphaW)+s),theta=G.involuteHalfAngle(p,R),alphaR=Math.acos(g.rb/R);
    const point=[R*Math.cos(theta),R*Math.sin(theta),0],direction=[-Math.sin(alphaR-theta),-Math.cos(alphaR-theta),0];
    return {R,theta,point,direction,s,ratio,alphaW,rw,rb:g.rb};
  }
  function resultant(mesh,force){const total=[0,0,0],moment=[0,0,0];mesh.nodes.forEach(([x,y,z],i)=>{const [fx,fy,fz]=force.slice(3*i,3*i+3);total[0]+=fx;total[1]+=fy;total[2]+=fz;moment[0]+=y*fz-z*fy;moment[1]+=z*fx-x*fz;moment[2]+=x*fy-y*fx;});return {force:total,moment};}
  function patch(mesh,p,{torqueNm=40,patchWidth=.1*G.basic(p).m,KHbeta=1,...pair}={}){
    if(!(patchWidth>0)||!(KHbeta>=1&&KHbeta<=2)||!(torqueNm>0))throw Error('Invalid patch load');
    const g=G.basic(p),contact=hpstc(p,pair),faces=mesh.faces.filter(f=>f.group==='tooth:0:plus:flank');
    if(!faces.length)throw Error('Loaded flank missing');
    const force=F.traction(mesh,faces,xyz=>{
      const R=Math.hypot(xyz[0],xyz[1]),distance=(R*R-contact.R*contact.R)/(2*g.rb);
      const profile=Math.abs(distance)<patchWidth/2?(1+Math.cos(2*Math.PI*distance/patchWidth))/patchWidth:0;
      const width=1+(KHbeta-1)*2*xyz[2]/g.b;
      return contact.direction.map(x=>x*profile*width/g.b);
    });
    const raw=resultant(mesh,force),target=-1000*torqueNm;
    if(!(raw.moment[2]<0))throw Error('Patch unresolved: refine loaded flank');
    const scale=target/raw.moment[2];for(let i=0;i<force.length;i++)force[i]*=scale;
    const balance=resultant(mesh,force),expectedFn=1000*torqueNm/g.rb,forceError=Math.abs(Math.hypot(...balance.force)/expectedFn-1);
    if(forceError>.005)throw Error('Patch force/moment mismatch: refine or reduce patch width');
    return {force,contact,resultant:balance,expectedFn,forceError,patchWidth,KHbeta};
  }
  const api={hpstc,resultant,patch};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.FEALoading=api;
})(typeof globalThis!=='undefined'?globalThis:this);
