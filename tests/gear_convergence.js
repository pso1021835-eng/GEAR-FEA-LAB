'use strict';
// Deliberately separate from fast CI. Run from any working directory with Node.
const {run,compare}=require('../src/analysis');
const refined=process.argv.includes('--refined');
const spatial=process.argv.includes('--3d'),meshOnly=process.argv.includes('--mesh-only');
const meshOptions=spatial?{
  axialSymmetry:true,coarsenBore:true,resolveCurvature:true,radialGrade:2,
  backgroundCounts:{root:1,fillet:8,involute:12},preconditioner:'block-jacobi'
}:refined?{
  axialSymmetry:true,coarsenBore:true,radialGrade:1.8,widthLayers:8,
  backgroundCounts:{root:1,fillet:7,involute:10},preconditioner:'block-jacobi'
}:{};
function sample([radialLayers,fillet,involute,widthLayers=6],sectorTeeth=spatial?7:5,cutFixed=false){
  const started=Date.now();
  const r=run({z:20},{bore:20},{...meshOptions,widthLayers,radialLayers,counts:{root:1,fillet,involute},sectorTeeth,cutFixed});
  console.log(JSON.stringify({radialLayers,fillet,involute,widthLayers,sectorTeeth,cutFixed,dofs:r.dofs,
    durationMs:Date.now()-started,iterations:r.iterations,residual:r.residual,
    rootMPa:r.root.central.maxPrincipal,filletSide:r.root.central.side,
    displacementMm:r.loadDisplacement,patchMeanDisplacementMm:r.patchMeanDisplacement,
    forceError:r.forceError,momentError:r.momentError,
    cutTorqueNmm:r.cutTorque,meshMetrics:r.meshMetrics}));
  return r;
}
const levels=spatial?[[9,20,12,8],[11,24,14,10],[13,30,16,12]]:
  refined?[[6,18,12,8],[7,18,12,8],[8,18,12,8]]:[[2,3,4],[2,4,6],[3,6,9]];
const coarse=sample(levels[0]);
const middle=sample(levels[1]),fine=sample(levels[2]);
const boundary={stressTolerance:.05,displacementTolerance:.05};
const {filletElementSpan,target}=fine.meshMetrics;
const result={mesh:compare(middle,fine),filletSize:{filletElementSpan,target,passes:filletElementSpan<=target}};
if(!meshOnly){
  // The spatial profile already meets the fillet size target at its first
  // level. Keep the sector comparisons at that common resolution.
  const comparable=spatial?coarse:refined?fine:middle;
  const sectorLevel=spatial?levels[0]:refined?levels[2]:levels[1];
  if(spatial){
    result.sector=compare(sample(sectorLevel,5),comparable,boundary);
    result.extendedSector=compare(comparable,sample(sectorLevel,9),boundary);
  }else result.sector=compare(comparable,sample(sectorLevel,7),boundary);
  result.cut=compare(comparable,sample(sectorLevel,spatial?7:5,true),boundary);
}
console.log(JSON.stringify(result,null,2));
if(!Object.values(result).every(check=>check.passes))process.exitCode=1;
