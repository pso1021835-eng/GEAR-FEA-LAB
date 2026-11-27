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
  return mean;
}
const coarse=beam(5,1),fine=beam(10,2);
assert(Math.abs(fine/.1941747572815534-1)<.03);
assert(Math.abs(fine/.1941747572815534-1)<Math.abs(coarse/.1941747572815534-1));
console.log('FEA A displacement refinement passed',coarse,fine);
