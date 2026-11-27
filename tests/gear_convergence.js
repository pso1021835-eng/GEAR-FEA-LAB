'use strict';
// Deliberately separate from fast CI. Run from any working directory with Node.
const {run,compare}=require('../src/analysis');
const refined=process.argv.includes('--refined');
const meshOptions=refined?{
  axialSymmetry:true,coarsenBore:true,radialGrade:1.8,widthLayers:8,
  backgroundCounts:{root:1,fillet:7,involute:10},preconditioner:'block-jacobi'
}:{};
function sample(radialLayers,fillet,involute,sectorTeeth=5,cutFixed=false){
  const started=Date.now();
  const r=run({z:20},{bore:20},{...meshOptions,radialLayers,counts:{root:1,fillet,involute},sectorTeeth,cutFixed});
  console.log(JSON.stringify({radialLayers,fillet,involute,sectorTeeth,cutFixed,dofs:r.dofs,
    durationMs:Date.now()-started,iterations:r.iterations,residual:r.residual,
    rootMPa:r.root.central.maxPrincipal,filletSide:r.root.central.side,
    displacementMm:r.loadDisplacement,patchMeanDisplacementMm:r.patchMeanDisplacement,
    forceError:r.forceError,momentError:r.momentError,
    cutTorqueNmm:r.cutTorque,meshMetrics:r.meshMetrics}));
  return r;
}
const levels=refined?[[6,18,12],[7,18,12],[8,18,12]]:[[2,3,4],[2,4,6],[3,6,9]];
sample(...levels[0]);
const middle=sample(...levels[1]),fine=sample(...levels[2]);
const comparable=refined?fine:middle,sectorLevel=refined?levels[2]:levels[1];
const seven=sample(...sectorLevel,7),fixed=sample(...sectorLevel,5,true);
const boundary={stressTolerance:.05,displacementTolerance:.05};
const {filletElementSpan,target}=fine.meshMetrics;
const result={mesh:compare(middle,fine),sector:compare(comparable,seven,boundary),
  cut:compare(comparable,fixed,boundary),filletSize:{filletElementSpan,target,passes:filletElementSpan<=target}};
console.log(JSON.stringify(result,null,2));
if(!Object.values(result).every(check=>check.passes))process.exitCode=1;
