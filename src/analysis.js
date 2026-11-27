/* Browser and Node entry point for a bore-clamped external spur sector. */
(function(root){
  'use strict';
  const node=typeof module!=='undefined'&&module.exports;
  const G=node?require('./gear'):root.GearGeom;
  const M=node?require('./mesh'):root.FEAMesh;
  const F=node?require('./fea'):root.FEASolver;
  const L=node?require('./loading'):root.FEALoading;

  function run(p,body,options={}){
    const {E=206000,nu=.3,sectorTeeth=5,radialLayers=2,widthLayers=6,radialGrade=2,coarsenBore=false,
      counts={root:1,fillet:4,involute:6},cutFixed=false,...loading}=options;
    const g=G.basic(p);
    if(!Number.isInteger(g.z)||g.z<3||g.z>300||![g.m,g.b,g.ha,g.hf].every(v=>Number.isFinite(v)&&v>0)
      ||!(g.alpha>0&&g.alpha<90)||!Number.isFinite(g.x)||!(g.rho>0)||!(body.bore>0))
      throw Error('잘못된 기어 또는 보어 치수입니다.');
    const section=G.section(p,body,{counts,nTip:3});
    if(!section.analyzable)throw Error('Geometry not analyzable');
    const mesh=M.gear(section,{sectorTeeth,radialLayers,widthLayers,radialGrade,coarsenBore,flankKnots:L.patchKnots(p,loading),boundaryRadius:theta=>G.boundaryRadius(p,theta)});
    const load=L.patch(mesh,p,loading),fixedNodes=new Set(),boreNodes=new Set();
    for(const face of mesh.faces){
      if(face.group==='bore')for(const i of face.nodes)boreNodes.add(i);
      if(face.group==='bore'||(cutFixed&&face.group.startsWith('cut:')))for(const i of face.nodes)fixedNodes.add(i);
    }
    if(!boreNodes.size)throw Error('No bore restraint');
    const fixed=[...fixedNodes].flatMap(i=>[3*i,3*i+1,3*i+2]);
    const solution=F.solve(F.assemble(mesh,E,nu),load.force,fixed);
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
    for(const face of mesh.faces){
      const match=/^tooth:0:(minus|plus):fillet$/.exec(face.group);
      if(!match)continue;
      const tetra=mesh.elements[face.element],bary=[0,0,0,0];
      for(const id of face.nodes.slice(0,3))bary[tetra.indexOf(id)]=1/3;
      const point=[0,1,2].map(d=>(-face.nodes.slice(0,3).reduce((s,i)=>s+mesh.nodes[i][d],0)+4*face.nodes.slice(3).reduce((s,i)=>s+mesh.nodes[i][d],0))/9);
      const stress=F.elementStress(mesh,tetra,solution.displacements,E,nu,bary);
      const sample={...stress,point,side:match[1]};
      if(!root.full||sample.maxPrincipal>root.full.maxPrincipal)root.full=sample;
      if(Math.abs(point[2])<=section.g.b/4&&(!root.central||sample.maxPrincipal>root.central.maxPrincipal))root.central=sample;
    }
    if(!root.central)throw Error('No central fillet samples');
    let work=0;for(let i=0;i<load.force.length;i++)work+=load.force[i]*solution.displacements[i];
    const maximumVonMises=F.stress(mesh,solution.displacements,E,nu).reduce((max,s)=>Math.max(max,s.vonMises),0);
    return {root,maximumVonMises,load,reaction,forceError,momentError,cutTorque,
      loadDisplacement:work/Math.hypot(...load.resultant.force),dofs:load.force.length,
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
