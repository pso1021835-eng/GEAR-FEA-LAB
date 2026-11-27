/* Browser and Node entry point for a bore-clamped external spur sector. */
(function(root){
  'use strict';
  const node=typeof module!=='undefined'&&module.exports;
  const G=node?require('./gear'):root.GearGeom;
  const M=node?require('./mesh'):root.FEAMesh;
  const F=node?require('./fea'):root.FEASolver;
  const L=node?require('./loading'):root.FEALoading;

  // Find the mesh surface point in the HPSTC radial plane at mid-width.
  function contactDisplacement(mesh,displacements,contact){
    const normal=[-Math.sin(contact.theta),Math.cos(contact.theta)];
    for(const face of mesh.faces){
      if(face.group!=='tooth:0:plus:flank')continue;
      const samples=face.nodes.map(i=>{
        const [x,y,z]=mesh.nodes[i];
        return [normal[0]*x+normal[1]*y,z];
      });
      const [p0,p1,p2]=samples;
      const a=p1.map((x,i)=>x-p0[i]),b=p2.map((x,i)=>x-p0[i]);
      const det=a[0]*b[1]-a[1]*b[0];
      if(Math.abs(det)<1e-14)continue;
      let u=(-p0[0]*b[1]+p0[1]*b[0])/det,v=(-a[0]*p0[1]+a[1]*p0[0])/det;
      if(Math.min(u,v,1-u-v)<-.02)continue;
      let shape;
      for(let iter=0;iter<10;iter++){
        const l=1-u-v;
        shape=[l*(2*l-1),u*(2*u-1),v*(2*v-1),4*l*u,4*u*v,4*v*l];
        const du=[1-4*l,4*u-1,0,4*(l-u),4*v,-4*v];
        const dv=[1-4*l,0,4*v-1,-4*u,4*u,4*(l-v)];
        const project=(w,i)=>w.reduce((sum,x,k)=>sum+x*samples[k][i],0);
        const x=project(shape,0),z=project(shape,1);
        if(Math.hypot(x,z)<1e-11)break;
        const xu=project(du,0),xv=project(dv,0),zu=project(du,1),zv=project(dv,1);
        const j=xu*zv-xv*zu;
        if(Math.abs(j)<1e-14){shape=null;break;}
        u-=(x*zv-z*xv)/j;v-=(z*xu-x*zu)/j;
        if(iter===9)shape=null;
      }
      if(!shape||Math.min(u,v,1-u-v)<-1e-8)continue;
      const interpolate=get=>[0,1,2].map(d=>shape.reduce((sum,n,k)=>sum+n*get(face.nodes[k],d),0));
      const point=interpolate((i,d)=>mesh.nodes[i][d]),displacement=interpolate((i,d)=>displacements[3*i+d]);
      return {point,directional:displacement.reduce((sum,x,d)=>sum+x*contact.direction[d],0)};
    }
    throw Error('HPSTC displacement probe missed the loaded flank');
  }

  function run(p,body,options={}){
    const {E=206000,nu=.3,sectorTeeth=5,radialLayers=2,widthLayers=6,radialGrade=2,coarsenBore=false,
      counts={root:1,fillet:4,involute:6},backgroundCounts=null,cutFixed=false,preconditioner='jacobi',axialSymmetry=false,mirrorAxial=false,...loading}=options;
    const g=G.basic(p);
    if(!Number.isInteger(g.z)||g.z<3||g.z>300||![g.m,g.b,g.ha,g.hf].every(v=>Number.isFinite(v)&&v>0)
      ||!(g.alpha>0&&g.alpha<90)||!Number.isFinite(g.x)||!(g.rho>0)||!(body.bore>0))
      throw Error('잘못된 기어 또는 보어 치수입니다.');
    const section=G.section(p,body,{counts,nTip:3});
    if(!section.analyzable)throw Error('Geometry not analyzable');
    if(backgroundCounts){
      const coarse=G.outline(p,{counts:backgroundCounts,nTip:3}).points;
      const start=Math.PI/2,end=start+2*Math.PI/g.z;
      const angle=q=>{
        const theta=Math.atan2(q.y,q.x);
        return theta<start-1e-10?theta+2*Math.PI:theta;
      };
      // Both tooth profiles share their exact root-circle endpoints.
      section.outer=[
        ...section.outer.filter(q=>angle(q)>start+1e-10&&angle(q)<end-1e-10),
        ...coarse.filter(q=>angle(q)<=start+1e-10||angle(q)>=end-1e-10)
      ].sort((a,b)=>angle(a)-angle(b));
      section.outer=section.outer.filter((q,i,all)=>i===0||Math.hypot(q.x-all[i-1].x,q.y-all[i-1].y)>1e-9*g.m);
    }
    if(axialSymmetry&&((loading.KHbeta??1)!==1||widthLayers%2!==0)){
      throw Error('Axial symmetry requires uniform face load and an even full-width layer count');
    }
    const meshSection=axialSymmetry?{...section,g:{...g,b:g.b/2}}:section;
    const mesh=M.gear(meshSection,{sectorTeeth,radialLayers,widthLayers:axialSymmetry?widthLayers/2:widthLayers,radialGrade,coarsenBore,mirrorAxial:!axialSymmetry&&mirrorAxial,flankKnots:L.patchKnots(p,loading),boundaryRadius:theta=>G.boundaryRadius(p,theta)});
    if(axialSymmetry)for(const point of mesh.nodes)point[2]+=g.b/4;
    if(axialSymmetry)loading.torqueNm=(loading.torqueNm??40)/2;
    const load=L.patch(mesh,p,loading),fixedNodes=new Set(),boreNodes=new Set();
    for(const face of mesh.faces){
      if(face.group==='bore')for(const i of face.nodes)boreNodes.add(i);
      if(face.group==='bore'||(cutFixed&&face.group.startsWith('cut:')))for(const i of face.nodes)fixedNodes.add(i);
    }
    if(!boreNodes.size)throw Error('No bore restraint');
    const fixed=[...fixedNodes].flatMap(i=>[3*i,3*i+1,3*i+2]);
    if(axialSymmetry)mesh.nodes.forEach((point,i)=>{
      if(Math.abs(point[2])<1e-10&&!fixedNodes.has(i))fixed.push(3*i+2);
    });
    F.checkConstraints(mesh,fixed);
    const solution=F.solve(F.assemble(mesh,E,nu),load.force,fixed,{preconditioner});
    const reactionVector=new Float64Array(load.force.length);
    for(const i of fixed)reactionVector[i]=solution.reactions[i];
    const reaction=L.resultant(mesh,reactionVector);
    let cutTorque=0;
    for(const i of fixedNodes)if(!boreNodes.has(i)){
      const [x,y]=mesh.nodes[i];cutTorque+=x*reactionVector[3*i+1]-y*reactionVector[3*i];
    }
    const error=(a,b,scale)=>Math.hypot(...a.map((x,i)=>x+b[i]))/scale;
    const forceError=error(reaction.force,load.resultant.force,load.expectedFn);
    const momentError=error(reaction.moment,load.resultant.moment,Math.abs(load.resultant.moment[2]));
    if(forceError>.005||momentError>.005)throw Error('Reaction balance exceeded 0.5%');
    const root={central:null,full:null};
    let filletSurfaceInterval=0,filletElementSpan=0;
    for(const face of mesh.faces){
      const match=/^tooth:0:(minus|plus):fillet$/.exec(face.group);
      if(!match)continue;
      const tetra=mesh.elements[face.element],bary=[0,0,0,0];
      const inPlaneSpan=(ids)=>Math.max(...ids.flatMap((a,i)=>ids.slice(i+1).map(b=>{
        const [x,y]=mesh.nodes[a],[u,v]=mesh.nodes[b];
        return Math.hypot(x-u,y-v);
      })));
      filletSurfaceInterval=Math.max(filletSurfaceInterval,inPlaneSpan(face.nodes.slice(0,3)));
      filletElementSpan=Math.max(filletElementSpan,inPlaneSpan(tetra.slice(0,4)));
      for(const id of face.nodes.slice(0,3))bary[tetra.indexOf(id)]=1/3;
      const point=[0,1,2].map(d=>(-face.nodes.slice(0,3).reduce((s,i)=>s+mesh.nodes[i][d],0)+4*face.nodes.slice(3).reduce((s,i)=>s+mesh.nodes[i][d],0))/9);
      const stress=F.elementStress(mesh,tetra,solution.displacements,E,nu,bary);
      const sample={...stress,point,side:match[1]};
      if(!root.full||sample.maxPrincipal>root.full.maxPrincipal)root.full=sample;
      if(Math.abs(point[2])<=section.g.b/4&&(!root.central||sample.maxPrincipal>root.central.maxPrincipal))root.central=sample;
    }
    if(!root.central)throw Error('No central fillet samples');
    let work=0;for(let i=0;i<load.force.length;i++)work+=load.force[i]*solution.displacements[i];
    const contactProbe=contactDisplacement(mesh,solution.displacements,load.contact);
    const maximumVonMises=F.stress(mesh,solution.displacements,E,nu).reduce((max,s)=>Math.max(max,s.vonMises),0);
    const filletMinRadius=G.filletMinRadius(p);
    return {root,maximumVonMises,load,reaction,forceError,momentError,cutTorque,
      meshMetrics:{filletMinRadius,filletSurfaceInterval,filletElementSpan,target:Math.min(filletMinRadius/4,g.m/8)},
      loadDisplacement:contactProbe.directional,contactProbe,axialSymmetry,
      patchMeanDisplacement:work/Math.hypot(...load.resultant.force),dofs:load.force.length,
      iterations:solution.iterations,residual:solution.relativeResidual};
  }
  function compare(coarse,fine,{stressTolerance=.03,displacementTolerance=.02}={}){
    const relative=(a,b)=>b===0?Infinity:Math.abs(a-b)/Math.abs(b);
    const stressChange=relative(coarse.root.central.maxPrincipal,fine.root.central.maxPrincipal);
    const displacementChange=relative(coarse.loadDisplacement,fine.loadDisplacement);
    const sameFillet=coarse.root.central.side===fine.root.central.side;
    return {stressChange,displacementChange,sameFillet,
      passes:stressChange<stressTolerance&&displacementChange<displacementTolerance&&sameFillet};
  }
  const api={run,compare};if(node)module.exports=api;else root.FEAAnalysis=api;
})(typeof globalThis!=='undefined'?globalThis:this);
