/* global importScripts, FEAAnalysis */
importScripts('gear.js','mesh.js','fea.js','loading.js','analysis.js');
self.onmessage=({data})=>{
  try{
    const result=FEAAnalysis.run(data.gear,data.body,data.options);
    self.postMessage({rootMPa:result.root.central.maxPrincipal,vmMPa:result.maximumVonMises});
  }catch(error){self.postMessage({error:error.message});}
};
