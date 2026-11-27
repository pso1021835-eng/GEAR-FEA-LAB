/* Straight-sided quadratic tetrahedral elasticity. Stress: xx yy zz xy yz zx. */
(function(root){
  'use strict';
  const edge=[[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]];
  const a=.5854101966249685,b=.1381966011250105;
  const quadrature=[[a,b,b,b],[b,a,b,b],[b,b,a,b],[b,b,b,a]];
  function material(E,nu){
    if(!(E>0)||!Number.isFinite(E)||!(nu> -1&&nu<.5))throw Error('Invalid elastic material');
    const mu=E/(2*(1+nu)),lambda=E*nu/((1+nu)*(1-2*nu));
    return Array.from({length:6},(_,i)=>Float64Array.from({length:6},(_,j)=>i<3&&j<3?lambda+(i===j?2*mu:0):(i===j?mu:0)));
  }
  function geometry(mesh,t){
    const p=t.slice(0,4).map(i=>mesh.nodes[i]),u=p[1].map((x,i)=>x-p[0][i]),v=p[2].map((x,i)=>x-p[0][i]),w=p[3].map((x,i)=>x-p[0][i]);
    const cross=(x,y)=>[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]],dot=(x,y)=>x.reduce((s,z,i)=>s+z*y[i],0);
    const det=dot(u,cross(v,w));if(!(det>0))throw Error('Nonpositive Jacobian');
    const g=[cross(v,w),cross(w,u),cross(u,v)].map(x=>x.map(y=>y/det));
    g.unshift([0,1,2].map(i=>-g[0][i]-g[1][i]-g[2][i]));
    // This implementation explicitly rejects curved midside nodes.
    edge.forEach(([i,j],k)=>{if(mesh.nodes[t[k+4]].some((x,d)=>Math.abs(x-(p[i][d]+p[j][d])/2)>1e-9))throw Error('Curved TET10 mapping not yet supported');});
    return {g,det};
  }
  function strainMatrix(g,L){
    const grad=L.map((x,i)=>g[i].map(y=>(4*x-1)*y));
    edge.forEach(([i,j])=>grad.push(g[i].map((x,d)=>4*(L[j]*x+L[i]*g[j][d]))));
    const B=Array.from({length:6},()=>new Float64Array(30));
    grad.forEach(([x,y,z],i)=>{const k=3*i;B[0][k]=x;B[1][k+1]=y;B[2][k+2]=z;B[3][k]=y;B[3][k+1]=x;B[4][k+1]=z;B[4][k+2]=y;B[5][k]=z;B[5][k+2]=x;});return B;
  }
  function assemble(mesh,E,nu){
    const D=material(E,nu),rows=Array.from({length:mesh.nodes.length*3},()=>new Map());
    for(const t of mesh.elements){const {g,det}=geometry(mesh,t),K=new Float64Array(900);
      for(const L of quadrature){const B=strainMatrix(g,L),DB=D.map(row=>Float64Array.from({length:30},(_,j)=>row.reduce((s,x,i)=>s+x*B[i][j],0)));
        for(let i=0;i<30;i++)for(let j=0;j<30;j++){let s=0;for(let k=0;k<6;k++)s+=B[k][i]*DB[k][j];K[i*30+j]+=s*det/24;}}
      const ids=t.flatMap(i=>[3*i,3*i+1,3*i+2]);
      for(let i=0;i<30;i++)for(let j=0;j<30;j++){const row=rows[ids[i]],id=ids[j];row.set(id,(row.get(id)||0)+K[i*30+j]);}
    }
    const ptr=new Uint32Array(rows.length+1),cols=[],values=[];
    rows.forEach((r,i)=>{for(const [j,v]of r){cols.push(j);values.push(v);}ptr[i+1]=cols.length;});
    return {ptr,cols:Uint32Array.from(cols),values:Float64Array.from(values),size:rows.length};
  }
  function multiply(A,x){const y=new Float64Array(A.size);for(let i=0;i<A.size;i++)for(let k=A.ptr[i];k<A.ptr[i+1];k++)y[i]+=A.values[k]*x[A.cols[k]];return y;}
  const dot=(x,y)=>x.reduce((s,v,i)=>s+v*y[i],0);
  function solve(A,force,fixed,{tolerance=1e-8,maxIterations=10000}={}){
    const n=A.size,locked=new Set(fixed),x=new Float64Array(n),r=Float64Array.from(force),diag=new Float64Array(n);
    if(force.length!==n||force.some(v=>!Number.isFinite(v)))throw Error('Invalid load vector');
    for(let i=0;i<n;i++){if(locked.has(i))r[i]=0;for(let k=A.ptr[i];k<A.ptr[i+1];k++)if(A.cols[k]===i)diag[i]=A.values[k];if(!(diag[i]>0))throw Error('Missing stiffness');}
    let z=Float64Array.from(r,(v,i)=>v/diag[i]),p=z.slice(),rz=dot(r,z),norm=Math.sqrt(dot(r,r)),relativeResidual=norm?1:0,iterations=0;
    while(relativeResidual>tolerance&&iterations<maxIterations){const q=multiply(A,p);for(const i of locked)q[i]=0;const pq=dot(p,q);if(!(pq>0))throw Error('Nonpositive stiffness / rigid mode');const alpha=rz/pq;
      for(let i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*q[i];}iterations++;relativeResidual=Math.sqrt(dot(r,r))/norm;if(relativeResidual<=tolerance)break;
      z=Float64Array.from(r,(v,i)=>v/diag[i]);const next=dot(r,z),beta=next/rz;for(let i=0;i<n;i++)p[i]=z[i]+beta*p[i];rz=next;
    }
    if(relativeResidual>tolerance)throw Error('PCG did not converge: '+relativeResidual);
    const reactions=multiply(A,x).map((v,i)=>v-force[i]);return {displacements:x,reactions,iterations,relativeResidual};
  }
  function stress(mesh,u,E,nu){const D=material(E,nu);return mesh.elements.map(t=>{const {g}=geometry(mesh,t),B=strainMatrix(g,[.25,.25,.25,.25]),v=t.flatMap(i=>[u[3*i],u[3*i+1],u[3*i+2]]),eps=B.map(row=>dot(row,v)),s=D.map(row=>dot(row,eps));const [x,y,z,xy,yz,zx]=s;return {stress:s,vonMises:Math.sqrt(((x-y)**2+(y-z)**2+(z-x)**2)/2+3*(xy*xy+yz*yz+zx*zx))};});}
  function traction(mesh,faces,vector){
    const f=new Float64Array(3*mesh.nodes.length);
    for(const face of faces){const p=face.nodes.slice(0,3).map(i=>mesh.nodes[i]),a=p[1].map((x,i)=>x-p[0][i]),b=p[2].map((x,i)=>x-p[0][i]);
      const area=Math.hypot(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])/2;
      // Constant traction on TRI6: vertex integrals vanish, each edge node gets A/3.
      for(const i of face.nodes.slice(3))for(let d=0;d<3;d++)f[3*i+d]+=area*vector[d]/3;
    }return f;
  }
  const api={material,geometry,strainMatrix,assemble,multiply,solve,stress,traction};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.FEASolver=api;
})(typeof globalThis!=='undefined'?globalThis:this);
