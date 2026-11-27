/* Conforming straight-sided TET10 extrusion. mm; edge order 01,12,20,03,13,23. */
(function (root) {
  'use strict';
  const edges = [[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]];
  function determinant(a,b,c,d) {
    const u=b.map((x,i)=>x-a[i]), v=c.map((x,i)=>x-a[i]), w=d.map((x,i)=>x-a[i]);
    return u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0]);
  }
  function elevate(vertices, tets) {
    const nodes=vertices.map(p=>p.slice()), shared=new Map();
    const elements=tets.map(input=>{
      const t=input.slice();
      let det=determinant(...t.map(i=>nodes[i]));
      if (!Number.isFinite(det)||Math.abs(det)<1e-14) throw Error('Degenerate tetrahedron');
      if(det<0) [t[1],t[2]]=[t[2],t[1]];
      return t.concat(edges.map(([a,b])=>{
        const i=t[a],j=t[b],key=Math.min(i,j)+':'+Math.max(i,j);
        if(!shared.has(key)) { shared.set(key,nodes.length); nodes.push(nodes[i].map((x,k)=>(x+nodes[j][k])/2)); }
        return shared.get(key);
      }));
    });
    return {nodes,elements,type:'TET10',edgeOrder:edges};
  }
  function extrude(points, triangles, width, layers=6, mirrorAxial=false) {
    if(!(width>0)||!Number.isInteger(layers)||layers<1) throw Error('Invalid extrusion');
    if(points.some(p=>p.length!==2||p.some(x=>!Number.isFinite(x)))) throw Error('Invalid 2D coordinates');
    const vertices=[],tets=[],n=points.length;
    if(mirrorAxial&&layers%2!==0)throw Error('Mirrored extrusion requires an even layer count');
    for(let k=0;k<=layers;k++) for(const p of points) vertices.push([p[0],p[1],width*(k/layers-.5)]);
    for(let k=0;k<layers;k++) for(const tri of triangles) {
      if(tri.length!==3||tri.some(i=>!Number.isInteger(i)||i<0||i>=n)||new Set(tri).size!==3) throw Error('Invalid triangle');
      // Global vertex ordering gives the same diagonal on every shared prism face.
      const flip=mirrorAxial&&k<layers/2,first=flip?k+1:k,step=flip?-n:n;
      const [a,b,c]=tri.slice().sort((x,y)=>x-y).map(i=>i+first*n), A=a+step,B=b+step,C=c+step;
      tets.push([a,b,c,C],[a,b,B,C],[a,A,B,C]);
    }
    return elevate(vertices,tets);
  }
  function box(length,height,width,nx=10,ny=2,nz=2) {
    if(![length,height,width].every(x=>Number.isFinite(x)&&x>0)||![nx,ny,nz].every(x=>Number.isInteger(x)&&x>0)) throw Error('Invalid box');
    const p=[],tri=[],id=(i,j)=>i*(ny+1)+j;
    for(let i=0;i<=nx;i++) for(let j=0;j<=ny;j++) p.push([length*i/nx,height*j/ny]);
    for(let i=0;i<nx;i++) for(let j=0;j<ny;j++) {const a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);tri.push([a,b,c],[a,c,d]);}
    return extrude(p,tri,width,nz);
  }
  function boundary(mesh, classify=()=> 'surface') {
    const faces=new Map(), local=[[0,2,1],[0,1,3],[1,2,3],[2,0,3]];
    mesh.elements.forEach((t,element)=>local.forEach(c=>{
      const v=c.map(i=>t[i]),key=v.slice().sort((a,b)=>a-b).join(':');
      if(faces.has(key)) { faces.get(key).count++; return; }
      const mids=c.map((a,i)=>{const b=c[(i+1)%3];return t[4+edges.findIndex(e=>e.includes(a)&&e.includes(b))];});
      faces.set(key,{nodes:v.concat(mids),element,count:1});
    }));
    const out=[];
    for(const f of faces.values()) {
      if(f.count>2) throw Error('Non-manifold tetrahedron face');
      if(f.count===1) { f.group=classify(f.nodes.map(i=>mesh.nodes[i]),f);out.push(f); }
    }
    return out;
  }
  function gear(section,{radialLayers=4,widthLayers=6,sectorTeeth=null,boundaryRadius=null,radialGrade=1,flankKnots=[],coarsenBore=false,mirrorAxial=false,resolveCurvature=false}={}){
    if(!section.analyzable)throw Error('Geometry not analyzable');
    if(section.holes.length!==1)throw Error('Only single circular bore supported by radial mesher');
    if(!Number.isInteger(radialLayers)||radialLayers<1)throw Error('Invalid radial layers');
    if(!Number.isFinite(radialGrade)||radialGrade<1)throw Error('Invalid radial grade');
    const hole=section.holes[0],radius=Math.hypot(hole[0].x,hole[0].y),tol=section.g.m*1e-8;
    if(!(radius>0))throw Error('Invalid bore radius');
    const keyed=hole.some(p=>Math.abs(Math.hypot(p.x,p.y)-radius)>tol);
    if(keyed&&(sectorTeeth!==null||coarsenBore))throw Error('Keyway topology requires a full gear and an uncoarsened bore');
    const angle=-(Math.PI/2+Math.PI/section.g.z),c=Math.cos(angle),s=Math.sin(angle);
    let innerAt=()=>radius;
    if(keyed){
      // GearGeom.boreHole appends the two roof corners after the round arc.
      const left=hole.at(-2),right=hole.at(-1),half=right.x,top=right.y;
      if(!(half>0&&half<radius&&top>Math.sqrt(radius*radius-half*half))||
        Math.abs(left.x+half)>tol||Math.abs(left.y-top)>tol||
        hole.slice(0,-2).some(p=>Math.abs(Math.hypot(p.x,p.y)-radius)>tol)){
        throw Error('Unrecognized keyway outline');
      }
      innerAt=theta=>{
        const x=Math.cos(theta-angle),y=Math.sin(theta-angle);
        return y>0?Math.max(radius,Math.min(top/y,half/Math.abs(x))):radius;
      };
    }
    let outer=section.outer.map(p=>({x:c*p.x-s*p.y,y:s*p.x+c*p.y,part:p.part}));
    if(outer.length>1&&Math.hypot(outer[0].x-outer[outer.length-1].x,outer[0].y-outer[outer.length-1].y)<tol)outer.pop();
    // Monotone polar outline is required: no undercut or radial reentry.
    let sweep=0;
    for(let i=0;i<outer.length;i++){const p=outer[i],q=outer[(i+1)%outer.length];let d=Math.atan2(q.y,q.x)-Math.atan2(p.y,p.x);while(d<0)d+=2*Math.PI;if(!(d>1e-12&&d<Math.PI))throw Error('Non-star-shaped or duplicate outline');sweep+=d;}
    if(Math.abs(sweep-2*Math.PI)>1e-6)throw Error('Invalid polar outline winding');
    if(sectorTeeth!==null){
      if(!Number.isInteger(sectorTeeth)||sectorTeeth<3||sectorTeeth%2!==1||sectorTeeth>=section.g.z)throw Error('Sector must have odd tooth count below gear count');
      const limit=sectorTeeth*Math.PI/section.g.z;
      const sorted=outer.map(p=>({...p,theta:Math.atan2(p.y,p.x)})).sort((a,b)=>a.theta-b.theta);
      const cut=theta=>({x:section.g.rf*Math.cos(theta),y:section.g.rf*Math.sin(theta),part:'root',theta});
      outer=[cut(-limit),...sorted.filter(p=>p.theta> -limit+1e-10&&p.theta<limit-1e-10),cut(limit)];
    }
    if(flankKnots.length){
      if(!boundaryRadius)throw Error('Flank knots require exact boundary radius');
      const min=flankKnots[0],max=flankKnots.at(-1);
      outer=outer.map(p=>({...p,theta:Math.atan2(p.y,p.x)})).filter(p=>p.theta<min-1e-10||p.theta>max+1e-10);
      for(const theta of flankKnots){const r=boundaryRadius(theta);outer.push({x:r*Math.cos(theta),y:r*Math.sin(theta),part:'involute',theta});}
      outer.sort((a,b)=>a.theta-b.theta);
    }
    if(keyed){
      if(!boundaryRadius)throw Error('Keyway topology requires exact tooth boundary radius');
      const corners=[hole[0],...hole.slice(-3)];
      outer=outer.map(p=>({...p,theta:Math.atan2(p.y,p.x)}));
      for(const p of corners){
        const theta=Math.atan2(s*p.x+c*p.y,c*p.x-s*p.y);
        if(outer.some(q=>Math.abs(Math.atan2(Math.sin(q.theta-theta),Math.cos(q.theta-theta)))<1e-10))continue;
        const r=boundaryRadius(theta);
        const before=outer.reduce((best,q)=>{
          const delta=(theta-q.theta+2*Math.PI)%(2*Math.PI);
          return delta<best.delta?{q,delta}:best;
        },{delta:Infinity}).q;
        outer.push({x:r*Math.cos(theta),y:r*Math.sin(theta),part:before.part,theta});
      }
      outer.sort((a,b)=>a.theta-b.theta);
    }
    if(resolveCurvature&&boundaryRadius){
      // Limit midside movement relative to the wall thickness of the nearest
      // radial layer, especially around the steep fillet-to-flank transition.
      const step=radialLayers**(-radialGrade);
      function subdivide(a,b,depth){
        const theta=Math.atan2(a.y+b.y,a.x+b.x);
        const chord=Math.hypot(a.x+b.x,a.y+b.y)/2,r=boundaryRadius(theta);
        const thickness=r-innerAt(theta);
        if(!(thickness>0))throw Error('Bore intersects gear outline');
        if(Math.abs(r-chord)<.15*thickness*step)return [a];
        if(depth===10)throw Error('Surface crosses an unresolved profile transition');
        const parts=[a.part,b.part];
        const part=parts.includes('tip')?'tip':parts.includes('involute')?'involute':parts.includes('fillet')?'fillet':'root';
        const mid={x:r*Math.cos(theta),y:r*Math.sin(theta),part};
        return [...subdivide(a,mid,depth+1),...subdivide(mid,b,depth+1)];
      }
      const limit=sectorTeeth?outer.length-1:outer.length,refined=[];
      for(let i=0;i<limit;i++)refined.push(...subdivide(outer[i],outer[(i+1)%outer.length],0));
      if(sectorTeeth)refined.push(outer.at(-1));
      outer=refined;
    }
    const points=[],tri=[],n=outer.length,rings=[],pointTags=[];
    for(let k=0;k<=radialLayers;k++){
      const ring=new Map();
      const stride=coarsenBore?2**Math.min(2,Math.max(0,radialLayers-k-1)):1;
      for(let j=0;j<n;j++){
        if(j%stride!==0&&!(sectorTeeth&&j===n-1))continue;
        const p=outer[j],R=Math.hypot(p.x,p.y);
        const inner=innerAt(Math.atan2(p.y,p.x));
        if(!(R>inner))throw Error('Bore intersects gear outline');
        const r=R-(R-inner)*(1-k/radialLayers)**radialGrade;
        ring.set(j,points.length);points.push([p.x*r/R,p.y*r/R]);pointTags.push({ring:k,column:j});
      }
      rings.push(ring);
    }
    const area=(i,j,k)=>{
      const [a,b,c]=[points[i],points[j],points[k]];
      return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
    };
    for(let k=0;k<radialLayers;k++){
      const inner=rings[k],outerRing=rings[k+1],keys=[...inner.keys()];
      if(!sectorTeeth)keys.push(n);
      for(let i=0;i<keys.length-1;i++){
        const first=keys[i],last=keys[i+1],a=inner.get(first),d=inner.get(last%n);
        const columns=[];
        for(let j=first;j<=last;j++)if(outerRing.has(j%n))columns.push(outerRing.get(j%n));
        let best=null;
        // Coarse-to-fine angular transitions can cut across a steep flank.
        // Try the admissible bridges in the local strip and retain positive area.
        for(let bridge=columns.length-1;bridge>=0;bridge--){
          const candidate=[];
          for(let j=0;j<bridge;j++)candidate.push([a,columns[j],columns[j+1]]);
          candidate.push([a,columns[bridge],d]);
          for(let j=bridge;j<columns.length-1;j++)candidate.push([d,columns[j],columns[j+1]]);
          const min=Math.min(...candidate.map(t=>area(...t)));
          if(min>0&&(!best||min>best.min))best={candidate,min};
          if(best&&columns.length===2)break;
        }
        if(!best)throw Error('Radial triangulation overlaps or degenerates');
        tri.push(...best.candidate);
      }
    }
    const mesh=extrude(points,tri,section.g.b,widthLayers,mirrorAxial),count=points.length;
    mesh.faces=boundary(mesh,(coords,f)=>{
      const tags=f.nodes.slice(0,3).map(i=>pointTags[i%count]);
      if(tags.every(p=>p.ring===0)){
        const theta=Math.atan2(coords[0][1]+coords[1][1]+coords[2][1],coords[0][0]+coords[1][0]+coords[2][0]);
        return innerAt(theta)>radius+tol?'keyway':'bore';
      }
      if(sectorTeeth&&tags.every(p=>p.column===0))return 'cut:minus';
      if(sectorTeeth&&tags.every(p=>p.column===n-1))return 'cut:plus';
      if(tags.every(p=>p.ring===radialLayers)){
        const mid=coords.slice(0,3).reduce((p,q)=>[p[0]+q[0]/3,p[1]+q[1]/3],[0,0]),theta=Math.atan2(mid[1],mid[0]),pitch=2*Math.PI/section.g.z,tooth=Math.round(theta/pitch),local=theta-tooth*pitch;
        const parts=tags.map(p=>outer[p.column].part);
        const part=parts.includes('tip')?'tip':parts.includes('involute')?'flank':parts.includes('fillet')?'fillet':'root';
        return 'tooth:'+((tooth+section.g.z)%section.g.z)+':'+(local<0?'minus':'plus')+':'+part;
      }
      return 'end';
    });
    if(boundaryRadius){
      const updated=new Set();
      for(const f of mesh.faces)if(f.group==='bore'||f.group==='keyway'||f.group.startsWith('tooth:'))for(const i of f.nodes.slice(3)){
        if(updated.has(i))continue;updated.add(i);
        const p=mesh.nodes[i],theta=Math.atan2(p[1],p[0]),r=f.group.startsWith('tooth:')?boundaryRadius(theta):innerAt(theta);
        if(!(r>0)||!Number.isFinite(r))throw Error('Invalid boundary mapping');p[0]=r*Math.cos(theta);p[1]=r*Math.sin(theta);
      }
    }
    mesh.metadata={rotation:angle,geometry:boundaryRadius?'curved TET10':'radial straight TET10',boreRadius:radius,radialLayers,radialGrade,widthLayers,sectorTeeth};return mesh;
  }
  const api={edges,determinant,elevate,extrude,box,boundary,gear};
  if(typeof module!=='undefined'&&module.exports) module.exports=api; else root.FEAMesh=api;
})(typeof globalThis!=='undefined'?globalThis:this);
