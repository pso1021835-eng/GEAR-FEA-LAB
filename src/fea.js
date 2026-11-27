/* Straight-sided quadratic tetrahedral elasticity. Stress: xx yy zz xy yz zx. */
(function(root){
  'use strict';
  const edge=[[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]];
  const a=.5854101966249685,b=.1381966011250105;
  const quadrature=[[a,b,b,b],[b,a,b,b],[b,b,a,b],[b,b,b,a]];
  const curvedQuadrature=[];
  const gp=[(1-Math.sqrt(3/5))/2,.5,(1+Math.sqrt(3/5))/2],gw=[5/18,4/9,5/18];
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++){
    const r=gp[i],s=gp[j],t=gp[k],L=[(1-r)*(1-s)*(1-t),r,(1-r)*s,(1-r)*(1-s)*t];
    curvedQuadrature.push({L,weight:gw[i]*gw[j]*gw[k]*(1-r)**2*(1-s)});
  }
  function mapping(mesh,t,L){
    const natural=[[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1]],grad=L.map((x,i)=>natural[i].map(y=>(4*x-1)*y));
    edge.forEach(([i,j])=>grad.push(natural[i].map((x,d)=>4*(L[j]*x+L[i]*natural[j][d]))));
    const J=Array.from({length:3},(_,d)=>[0,1,2].map(k=>grad.reduce((s,g,i)=>s+mesh.nodes[t[i]][d]*g[k],0)));
    const u=J.map(r=>r[0]),v=J.map(r=>r[1]),w=J.map(r=>r[2]),cross=(x,y)=>[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]],det=dot(u,cross(v,w));
    if(!(det>0)||!Number.isFinite(det))throw Error('Nonpositive curved Jacobian');
    const g=[cross(v,w),cross(w,u),cross(u,v)].map(x=>x.map(y=>y/det));g.unshift([0,1,2].map(i=>-g[0][i]-g[1][i]-g[2][i]));return {g,det};
  }
  function isCurved(mesh,t){return edge.some(([i,j],k)=>mesh.nodes[t[k+4]].some((x,d)=>Math.abs(x-(mesh.nodes[t[i]][d]+mesh.nodes[t[j]][d])/2)>1e-9));}
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
  function shapeGradients(g,L){
    const grad=L.map((x,i)=>g[i].map(y=>(4*x-1)*y));
    edge.forEach(([i,j])=>grad.push(g[i].map((x,d)=>4*(L[j]*x+L[i]*g[j][d]))));
    return grad;
  }
  function strainMatrix(g,L){
    const grad=shapeGradients(g,L);
    const B=Array.from({length:6},()=>new Float64Array(30));
    grad.forEach(([x,y,z],i)=>{const k=3*i;B[0][k]=x;B[1][k+1]=y;B[2][k+2]=z;B[3][k]=y;B[3][k+1]=x;B[4][k+1]=z;B[4][k+2]=y;B[5][k]=z;B[5][k+2]=x;});return B;
  }
  function assemble(mesh,E,nu){
    const D=material(E,nu),mu=D[3][3],lambda=D[0][1],rows=Array.from({length:mesh.nodes.length*3},()=>new Map());
    for(const t of mesh.elements){const curved=isCurved(mesh,t),base=curved?null:geometry(mesh,t),K=new Float64Array(900);
      const rule=curved?curvedQuadrature:quadrature.map(L=>({L,weight:1/24}));
      if(curved)for(const L of [[.25,.25,.25,.25],[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])mapping(mesh,t,L);
      for(const {L,weight} of rule){
        const {g,det}=base||mapping(mesh,t,L),grad=shapeGradients(g,L),volume=det*weight;
        for(let i=0;i<10;i++)for(let j=i;j<10;j++){
          const a=grad[i],b=grad[j],diagonal=mu*dot(a,b);
          for(let x=0;x<3;x++)for(let y=0;y<3;y++){
            const row=3*i+x,col=3*j+y;
            if(col<row)continue;
            const value=volume*(lambda*a[x]*b[y]+mu*a[y]*b[x]+(x===y?diagonal:0));
            K[row*30+col]+=value;if(row!==col)K[col*30+row]+=value;
          }
        }
      }
      const ids=t.flatMap(i=>[3*i,3*i+1,3*i+2]);
      for(let i=0;i<30;i++)for(let j=0;j<30;j++){const row=rows[ids[i]],id=ids[j];row.set(id,(row.get(id)||0)+K[i*30+j]);}
    }
    const ptr=new Uint32Array(rows.length+1),cols=[],values=[];
    rows.forEach((r,i)=>{for(const [j,v]of r){cols.push(j);values.push(v);}ptr[i+1]=cols.length;});
    // All three components of a node were inserted together for whenever that
    // node appeared in an element. Their CSR rows share the column ordering.
    return {ptr,cols:Uint32Array.from(cols),values:Float64Array.from(values),size:rows.length,triplets:true};
  }
  function multiply(A,x,y=new Float64Array(A.size)){
    const {ptr,cols,values}=A;
    if(A.triplets){
      for(let i=0;i<A.size;i+=3){
        const start=ptr[i],end=ptr[i+1],length=end-start;
        let a=0,b=0,c=0;
        for(let k=start;k<end;k+=3){
          const index=cols[k],u=x[index],v=x[index+1],w=x[index+2];
          const p=k+length,q=p+length;
          a+=values[k]*u+values[k+1]*v+values[k+2]*w;
          b+=values[p]*u+values[p+1]*v+values[p+2]*w;
          c+=values[q]*u+values[q+1]*v+values[q+2]*w;
        }
        y[i]=a;y[i+1]=b;y[i+2]=c;
      }
      return y;
    }
    for(let i=0;i<A.size;i++){
      let sum=0;for(let k=ptr[i],end=ptr[i+1];k<end;k++)sum+=values[k]*x[cols[k]];
      y[i]=sum;
    }
    return y;
  }
  function dot(x,y){let sum=0;for(let i=0;i<x.length;i++)sum+=x[i]*y[i];return sum;}
  function checkConstraints(mesh,fixed){
    const count=mesh.nodes.length,parent=Int32Array.from({length:count},(_,i)=>i);
    function root(i){while(parent[i]!==i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
    for(const t of mesh.elements)for(const i of t)parent[root(i)]=root(t[0]);
    const components=new Map();
    for(let i=0;i<count;i++){const j=root(i);if(!components.has(j))components.set(j,[]);}
    for(const dof of fixed){
      if(!Number.isInteger(dof)||dof<0||dof>=3*count)throw Error('Invalid restraint DOF');
      const i=Math.floor(dof/3),j=root(i),axis=dof%3;
      const p=mesh.nodes[i].map((x,d)=>x-mesh.nodes[j][d]);
      components.get(j).push({axis,p});
    }
    for(const rows of components.values()){
      const scale=rows.reduce((max,{p})=>Math.max(max,Math.hypot(...p)),0)||1;
      const basis=[];
      for(const {axis,p} of rows){
        const [x,y,z]=p.map(v=>v/scale);
        const v=axis===0?[1,0,0,0,z,-y]:axis===1?[0,1,0,-z,0,x]:[0,0,1,y,-x,0];
        for(const q of basis){
          const projection=dot(v,q);
          for(let k=0;k<6;k++)v[k]-=projection*q[k];
        }
        const magnitude=Math.hypot(...v);
        if(magnitude>1e-9)basis.push(v.map(x=>x/magnitude));
        if(basis.length===6)break;
      }
      if(basis.length<6)throw Error('A mesh component has an unrestrained rigid-body mode');
    }
  }
  // Each elastic node has a 3x3 diagonal block. Constrained rows are replaced
  // by identity in this preconditioner, so partial constraints remain valid.
  function blockJacobi(A,locked,diag){
    const blocks=new Float64Array(6*(A.size/3));
    for(let i=0,k=0;i<A.size;i+=3,k+=6){
      const a=locked.has(i)?1:diag[i],d=locked.has(i+1)?1:diag[i+1],f=locked.has(i+2)?1:diag[i+2];
      const entry=(row,col)=>{
        if(locked.has(row)||locked.has(col))return 0;
        for(let j=A.ptr[row];j<A.ptr[row+1];j++)if(A.cols[j]===col)return A.values[j];
        return 0;
      };
      const b=entry(i,i+1),c=entry(i,i+2),e=entry(i+1,i+2);
      const cof=[d*f-e*e,c*e-b*f,b*e-c*d,a*f-c*c,b*c-a*e,a*d-b*b];
      const det=a*cof[0]+b*cof[1]+c*cof[2];
      if(!(det>0)||!Number.isFinite(det))throw Error('Nonpositive Jacobi block');
      for(let j=0;j<6;j++)blocks[k+j]=cof[j]/det;
    }
    return (r,z)=>{
      for(let i=0,k=0;i<A.size;i+=3,k+=6){
        const x=r[i],y=r[i+1],w=r[i+2];
        z[i]=blocks[k]*x+blocks[k+1]*y+blocks[k+2]*w;
        z[i+1]=blocks[k+1]*x+blocks[k+3]*y+blocks[k+4]*w;
        z[i+2]=blocks[k+2]*x+blocks[k+4]*y+blocks[k+5]*w;
      }
    };
  }
  function solve(A,force,fixed,{tolerance=1e-8,maxIterations=10000,preconditioner='jacobi'}={}){
    const n=A.size,locked=new Set(fixed),x=new Float64Array(n),r=Float64Array.from(force),diag=new Float64Array(n);
    if(force.length!==n||force.some(v=>!Number.isFinite(v)))throw Error('Invalid load vector');
    for(let i=0;i<n;i++){if(locked.has(i))r[i]=0;for(let k=A.ptr[i];k<A.ptr[i+1];k++)if(A.cols[k]===i)diag[i]=A.values[k];if(!(diag[i]>0))throw Error('Missing stiffness');}
    let precondition;
    switch(preconditioner){
      case 'jacobi':precondition=(r,z)=>{for(let i=0;i<n;i++)z[i]=r[i]/diag[i];};break;
      case 'block-jacobi':precondition=blockJacobi(A,locked,diag);break;
      default:throw Error('Unknown preconditioner');
    }
    const z=new Float64Array(n),q=new Float64Array(n);
    precondition(r,z);
    const p=z.slice();
    let rz=dot(r,z),norm=Math.sqrt(dot(r,r)),relativeResidual=norm?1:0,iterations=0;
    while(relativeResidual>tolerance&&iterations<maxIterations){multiply(A,p,q);for(const i of locked)q[i]=0;const pq=dot(p,q);if(!(pq>0))throw Error('Nonpositive stiffness / rigid mode');const alpha=rz/pq;
      for(let i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*q[i];}iterations++;relativeResidual=Math.sqrt(dot(r,r))/norm;if(relativeResidual<=tolerance)break;
      precondition(r,z);const next=dot(r,z),beta=next/rz;for(let i=0;i<n;i++)p[i]=z[i]+beta*p[i];rz=next;
    }
    if(relativeResidual>tolerance)throw Error('PCG did not converge: '+relativeResidual);
    const reactions=multiply(A,x).map((v,i)=>v-force[i]);
    let freeResidual=0;
    for(let i=0;i<n;i++)if(!locked.has(i))freeResidual+=reactions[i]**2;
    relativeResidual=norm?Math.sqrt(freeResidual)/norm:0;
    if(relativeResidual>tolerance)throw Error('Recomputed PCG residual exceeds target: '+relativeResidual);
    return {displacements:x,reactions,iterations,relativeResidual};
  }
  function invariants(s){
    const [x,y,z,xy,yz,zx]=s,mean=(x+y+z)/3,p=Math.sqrt(((x-mean)**2+(y-mean)**2+(z-mean)**2+2*(xy*xy+yz*yz+zx*zx))/6);
    let principal;
    if(p<1e-14*Math.max(1,...s.map(Math.abs)))principal=[mean,mean,mean];
    else {const a=(x-mean)/p,b=(y-mean)/p,c=(z-mean)/p,d=xy/p,e=yz/p,f=zx/p,r=Math.max(-1,Math.min(1,(a*b*c+2*d*e*f-a*e*e-b*f*f-c*d*d)/2)),theta=Math.acos(r)/3;
      const high=mean+2*p*Math.cos(theta),low=mean+2*p*Math.cos(theta+2*Math.PI/3);principal=[high,3*mean-high-low,low];}
    return {stress:s,principal,maxPrincipal:principal[0],vonMises:Math.sqrt(((x-y)**2+(y-z)**2+(z-x)**2)/2+3*(xy*xy+yz*yz+zx*zx))};
  }
  function elementStress(mesh,t,u,E,nu,L=[.25,.25,.25,.25]){
    const {g}=isCurved(mesh,t)?mapping(mesh,t,L):geometry(mesh,t),B=strainMatrix(g,L),v=t.flatMap(i=>[u[3*i],u[3*i+1],u[3*i+2]]),eps=B.map(row=>dot(row,v));return invariants(material(E,nu).map(row=>dot(row,eps)));
  }
  function stress(mesh,u,E,nu){return mesh.elements.map(t=>elementStress(mesh,t,u,E,nu));}
  function probe(mesh,u,E,nu,point){
    const candidates=[];
    mesh.elements.forEach((t,element)=>{const {g}=geometry(mesh,t),origin=mesh.nodes[t[0]],delta=point.map((x,i)=>x-origin[i]),L=[0,...g.slice(1).map(row=>dot(row,delta))];L[0]=1-L[1]-L[2]-L[3];
      if(L.every(x=>x>=-1e-8&&x<=1+1e-8))candidates.push({element,...elementStress(mesh,t,u,E,nu,L)});
    });
    if(!candidates.length)throw Error('Stress probe outside mesh');
    // Preserve element-side values. Caller chooses a consistent recovery convention.
    return candidates;
  }
  function traction(mesh,faces,vector){
    const f=new Float64Array(3*mesh.nodes.length);
    for(const face of faces)for(let i=0;i<3;i++)for(let j=0;j<3;j++){
      const r=gp[i],s=gp[j],L=[(1-r)*(1-s),r,(1-r)*s],N=L.map(x=>x*(2*x-1)),ng=[[-1,-1],[1,0],[0,1]],dN=L.map((x,k)=>ng[k].map(y=>(4*x-1)*y));
      for(const [a,b]of [[0,1],[1,2],[2,0]]){N.push(4*L[a]*L[b]);dN.push(ng[a].map((x,k)=>4*(L[b]*x+L[a]*ng[b][k])));}
      const xyz=[0,1,2].map(d=>N.reduce((v,n,k)=>v+n*mesh.nodes[face.nodes[k]][d],0));
      const tangent=[0,1].map(k=>[0,1,2].map(d=>dN.reduce((v,n,a)=>v+n[k]*mesh.nodes[face.nodes[a]][d],0))),[a,b]=tangent;
      const cross=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],jac=Math.hypot(...cross);
      if(!(jac>0))throw Error('Degenerate load face');
      const q=typeof vector==='function'?vector(xyz,face,cross.map(x=>x/jac)):vector;
      if(q.length!==3||q.some(x=>!Number.isFinite(x)))throw Error('Invalid traction');
      const weight=gw[i]*gw[j]*(1-r)*jac;
      for(let k=0;k<6;k++)for(let d=0;d<3;d++)f[3*face.nodes[k]+d]+=N[k]*q[d]*weight;
    }return f;
  }
  const api={material,geometry,mapping,isCurved,strainMatrix,assemble,multiply,solve,stress,traction,invariants,elementStress,probe,checkConstraints};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.FEASolver=api;
})(typeof globalThis!=='undefined'?globalThis:this);
