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
  function extrude(points, triangles, width, layers=6) {
    if(!(width>0)||!Number.isInteger(layers)||layers<1) throw Error('Invalid extrusion');
    if(points.some(p=>p.length!==2||p.some(x=>!Number.isFinite(x)))) throw Error('Invalid 2D coordinates');
    const vertices=[],tets=[],n=points.length;
    for(let k=0;k<=layers;k++) for(const p of points) vertices.push([p[0],p[1],width*(k/layers-.5)]);
    for(let k=0;k<layers;k++) for(const tri of triangles) {
      if(tri.length!==3||tri.some(i=>!Number.isInteger(i)||i<0||i>=n)||new Set(tri).size!==3) throw Error('Invalid triangle');
      // Global vertex ordering gives the same diagonal on every shared prism face.
      const [a,b,c]=tri.slice().sort((x,y)=>x-y).map(i=>i+k*n), A=a+n,B=b+n,C=c+n;
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
  const api={edges,determinant,elevate,extrude,box,boundary};
  if(typeof module!=='undefined'&&module.exports) module.exports=api; else root.FEAMesh=api;
})(typeof globalThis!=='undefined'?globalThis:this);
