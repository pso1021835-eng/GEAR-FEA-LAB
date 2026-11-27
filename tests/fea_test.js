'use strict';
const assert=require('assert'),M=require('../src/mesh'),F=require('../src/fea');
// B: exact uniaxial stress with Poisson contraction, no transverse end clamp.
const mesh=M.box(40,10,10,4,1,1),faces=M.boundary(mesh),end=faces.filter(f=>f.nodes.every(i=>Math.abs(mesh.nodes[i][0]-40)<1e-9));
const force=F.traction(mesh,end,[100,0,0]),fixed=[];
mesh.nodes.forEach(([x,y,z],i)=>{if(x===0)fixed.push(3*i);if(x===0&&y===0&&z===-5)fixed.push(3*i+1,3*i+2);if(x===0&&y===10&&z===-5)fixed.push(3*i+2);});
const result=F.solve(F.assemble(mesh,206000,.3),force,fixed);
const stress=F.stress(mesh,result.displacements,206000,.3);
for(const s of stress){assert(Math.abs(s.stress[0]-100)<.1);assert(s.stress.slice(1).every(x=>Math.abs(x)<.1));assert(Math.abs(s.vonMises-100)<.1);}
mesh.nodes.forEach(([x,y,z],i)=>{assert(Math.abs(result.displacements[3*i]-100*x/206000)<1e-8);assert(Math.abs(result.displacements[3*i+1]+.3*100*y/206000)<1e-8);assert(Math.abs(result.displacements[3*i+2]+.3*100*(z+5)/206000)<1e-8);});
assert(Math.abs(fixed.filter(i=>i%3===0).reduce((s,i)=>s+result.reactions[i],0)+10000)<.001);
assert(faces.length===36);
console.log('FEA B patch passed',result.iterations,result.relativeResidual);
// A: end shear, clamp; displacement plus force and moment balance.
function beam(nx,ny){
  const m=M.box(100,10,10,nx,ny,ny),f=F.traction(m,M.boundary(m).filter(s=>s.nodes.every(i=>m.nodes[i][0]===100)),[0,1,0]),lock=[];
  m.nodes.forEach((p,i)=>{if(p[0]===0)lock.push(3*i,3*i+1,3*i+2);});
  const r=F.solve(F.assemble(m,206000,.3),f,lock),end=m.nodes.map((p,i)=>p[0]===100?r.displacements[3*i+1]:null).filter(x=>x!==null);
  const mean=end.reduce((a,b)=>a+b)/end.length;
  let fy=0,mz=0;for(let i=0;i<m.nodes.length;i++)if(m.nodes[i][0]===0){fy+=r.reactions[3*i+1];mz+=m.nodes[i][0]*r.reactions[3*i+1]-m.nodes[i][1]*r.reactions[3*i];}
  assert(Math.abs(fy+100)<.005);assert(Math.abs(mz+10000)<.5);
  const s=F.stress(m,r.displacements,206000,.3);
  assert(s.every(v=>Number.isFinite(v.vonMises)));
  if(nx===10){const probes=F.probe(m,r.displacements,206000,.3,[20,0,0]);const meanSxx=probes.reduce((sum,s)=>sum+s.stress[0],0)/probes.length;assert(Math.abs(meanSxx/48-1)<.05);}
  return mean;
}
const coarse=beam(5,1),fine=beam(10,2);
assert(Math.abs(fine/.1941747572815534-1)<.03);
assert(Math.abs(fine/.1941747572815534-1)<Math.abs(coarse/.1941747572815534-1));
console.log('FEA A displacement refinement passed',coarse,fine);
assert.deepStrictEqual(F.invariants([10,10,10,0,0,0]).principal,[10,10,10]);
assert(Math.abs(F.invariants([0,0,0,100,0,0]).maxPrincipal-100)<1e-10);
const G=require('../src/gear');
const gear=M.gear(G.section({z:20},{bore:20},{counts:{root:1,fillet:3,involute:4},nTip:2}),{radialLayers:2,widthLayers:2});
assert(gear.faces.some(f=>f.group==='tooth:0:minus:flank'));
assert(gear.faces.some(f=>f.group==='tooth:0:plus:fillet'));
assert(gear.faces.some(f=>f.group==='bore'));
assert(gear.elements.every(t=>M.determinant(...t.slice(0,4).map(i=>gear.nodes[i]))>0));
assert.throws(()=>M.gear({analyzable:false}),/not analyzable/);
assert.throws(()=>M.gear(G.section({z:20},{bore:20,key:{w:4,t:2}})),/Keyways/);
console.log('Gear radial topology and face ID checks passed');
for(const sectorTeeth of [5,7]){
  const sector=M.gear(G.section({z:20},{bore:20},{counts:{root:1,fillet:3,involute:4},nTip:2}),{sectorTeeth,radialLayers:2,widthLayers:2});
  for(const sign of ['minus','plus'])assert(sector.faces.some(f=>f.group==='cut:'+sign));
  assert(sector.elements.length<gear.elements.length);
}
// Curved element must reproduce an affine physical displacement exactly.
const curved=M.elevate([[0,0,0],[2,0,0],[0,2,0],[0,0,2]],[[0,1,2,3]]);
curved.nodes[4][1]=-.1;
const affine=Float64Array.from(curved.nodes.flatMap(([x,y,z])=>[.001*x,-.0003*y,-.0003*z]));
assert(F.isCurved(curved,curved.elements[0]));
for(const L of [[.25,.25,.25,.25],[.7,.1,.1,.1]]){
  const s=F.elementStress(curved,curved.elements[0],affine,206000,.3,L);
  assert(Math.abs(s.stress[0]-206)<1e-9);assert(s.stress.slice(1).every(x=>Math.abs(x)<1e-9));
}
const K=F.assemble(curved,206000,.3),rigid=new Float64Array(30);for(let i=0;i<10;i++)rigid[3*i]=1;
assert(Math.max(...F.multiply(K,rigid).map(Math.abs))<1e-8);
console.log('Sector boundaries and curved affine strain passed');
const Load=require('../src/loading'),contact=Load.hpstc({z:20});
assert(Math.abs(contact.R-20.304798860015772)<1e-10);
assert(Math.abs(contact.s-.8456863940860675)<1e-10);
assert(Math.abs(contact.ratio-1.6351859635714583)<1e-10);
const curvedGear=M.gear(G.section({z:20},{bore:20},{counts:{root:1,fillet:4,involute:40},nTip:2}),{sectorTeeth:5,radialLayers:1,widthLayers:2,boundaryRadius:t=>G.boundaryRadius({z:20},t)});
const load=Load.patch(curvedGear,{z:20});assert(Math.abs(load.resultant.moment[2]+40000)<1e-8);assert(load.forceError<.005);
for(const p of G.outline({z:20}).points){const theta=Math.atan2(p.y,p.x)-Math.PI/2-Math.PI/20;assert(Math.abs(Math.hypot(p.x,p.y)-G.boundaryRadius({z:20},theta))<1e-8);}
console.log('Exact rack boundary and C HPSTC finite patch balance passed');
// C: physical equilibrium and finite root response. The hand comparison is diagnostic.
const A=require('../src/analysis');
const c=A.run({z:20},{bore:20},{widthLayers:2,counts:{root:1,fillet:4,involute:8}});
assert(c.residual<1e-8&&c.forceError<.005&&c.momentError<.005);
assert(c.root.central.maxPrincipal>0&&c.root.full.maxPrincipal>=c.root.central.maxPrincipal);
assert(c.loadDisplacement>0&&Number.isFinite(c.loadDisplacement));
const comparison=c.root.central.maxPrincipal/144.67427765268215-1;
console.log('C equilibrium passed, root MPa and hand-reference deviation:',c.root.central.maxPrincipal,comparison);
