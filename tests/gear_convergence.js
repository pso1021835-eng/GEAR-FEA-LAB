'use strict';
// Deliberately separate from fast CI. Run from any working directory with Node.
const {run,compare}=require('../src/analysis');
function sample(radialLayers,fillet,involute,sectorTeeth=5,cutFixed=false){
  const started=Date.now();
  const r=run({z:20},{bore:20},{radialLayers,counts:{root:1,fillet,involute},sectorTeeth,cutFixed});
  console.log(JSON.stringify({radialLayers,fillet,involute,sectorTeeth,cutFixed,dofs:r.dofs,
    durationMs:Date.now()-started,iterations:r.iterations,residual:r.residual,
    rootMPa:r.root.central.maxPrincipal,filletSide:r.root.central.side,
    displacementMm:r.loadDisplacement,patchMeanDisplacementMm:r.patchMeanDisplacement,
    forceError:r.forceError,momentError:r.momentError,
    cutTorqueNmm:r.cutTorque}));
  return r;
}
sample(2,3,4);
const middle=sample(2,4,6),fine=sample(3,6,9);
const seven=sample(2,4,6,7),fixed=sample(2,4,6,5,true);
const boundary={stressTolerance:.05,displacementTolerance:.05};
const result={mesh:compare(middle,fine),sector:compare(middle,seven,boundary),
  cut:compare(middle,fixed,boundary)};
console.log(JSON.stringify(result,null,2));
if(!Object.values(result).every(check=>check.passes))process.exitCode=1;
