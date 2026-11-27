/* app.js — 입력 패널 · 3D 기어 표시 (FEA-F1). 해석(FEA)은 F2 이후 단계에서 추가 */
(function () {
  'use strict';
  const G = window.GearGeom;
  const $ = id => document.getElementById(id);
  const FIELDS = [
    // [id, 라벨, 기본값, 단위, step]
    ['m', '모듈 m', 2, 'mm', 0.25], ['z', '잇수 z', 30, '', 1], ['alpha', '압력각 α', 20, '°', 0.5],
    ['x', '전위계수 x', 0, '', 0.05], ['rho', '공구 이끝 둥글기 ρf*', 0.38, '×m', 0.01], ['hf', '이뿌리 높이 hf*', 1.25, '×m', 0.05],
    ['b', '치폭 b', 20, 'mm', 1],
    ['bore', '보어 지름', 20, 'mm', 1], ['keyW', '키홈 폭', 6, 'mm', 0.5], ['keyT', '키홈 깊이(보어 기준)', 2.8, 'mm', 0.1],
    ['webN', '웨브 구멍 수', 0, '개', 1], ['webPcd', '웨브 구멍 PCD', 36, 'mm', 1], ['webD', '웨브 구멍 지름', 8, 'mm', 0.5]
  ];
  const GROUPS = [['치형', ['m', 'z', 'alpha', 'x', 'rho', 'hf', 'b']], ['몸체', ['bore', 'keyW', 'keyT', 'webN', 'webPcd', 'webD']]];

  function buildForm() {
    const host = $('inputs'); host.innerHTML = '';
    for (const [title, ids] of GROUPS) {
      const fs = document.createElement('fieldset'); fs.innerHTML = `<legend>${title}</legend>`;
      for (const id of ids) {
        const [, label, val, unit, step] = FIELDS.find(f => f[0] === id);
        const row = document.createElement('label'); row.className = 'row';
        row.innerHTML = `<span>${label}</span><input id="in-${id}" type="number" value="${val}" step="${step}"><em>${unit}</em>`;
        fs.appendChild(row);
      }
      host.appendChild(fs);
    }
    host.addEventListener('input', () => { clearTimeout(buildForm.t); buildForm.t = setTimeout(update, 120); });
  }
  const val = id => parseFloat($('in-' + id).value);

  function params() {
    return {
      gear: { m: val('m'), z: Math.round(val('z')), alpha: val('alpha'), x: val('x'), rho: val('rho'), hf: val('hf'), b: val('b') },
      body: { bore: val('bore'), key: { w: val('keyW'), t: val('keyT') }, web: { n: Math.round(val('webN')), pcd: val('webPcd'), d: val('webD') } }
    };
  }

  // ---------------------------------------------------------------- 3D
  let renderer, scene, camera, mesh, edges, cam = { th: 0.6, ph: 0.9, R: 200 };
  function init3D() {
    const el = $('view');
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    el.appendChild(renderer.domElement);
    scene = new THREE.Scene(); scene.background = new THREE.Color(0x0d1118);
    camera = new THREE.PerspectiveCamera(35, 1, 0.1, 5000);
    scene.add(new THREE.HemisphereLight(0xdfe8f5, 0x202530, 0.8));
    const d1 = new THREE.DirectionalLight(0xffffff, 0.8); d1.position.set(1, 2, 3); scene.add(d1);
    const d2 = new THREE.DirectionalLight(0xffffff, 0.3); d2.position.set(-2, -1, -1); scene.add(d2);
    // 궤도 조작: 드래그 회전, 휠 확대
    let drag = null;
    el.addEventListener('pointerdown', e => { drag = { x: e.clientX, y: e.clientY }; el.setPointerCapture(e.pointerId); });
    el.addEventListener('pointerup', () => drag = null);
    el.addEventListener('pointermove', e => {
      if (!drag) return;
      cam.th -= (e.clientX - drag.x) * 0.008; cam.ph = Math.min(3.1, Math.max(0.05, cam.ph - (e.clientY - drag.y) * 0.008));
      drag = { x: e.clientX, y: e.clientY }; render();
    });
    el.addEventListener('wheel', e => { e.preventDefault(); cam.R *= Math.exp(e.deltaY * 0.001); render(); }, { passive: false });
    window.addEventListener('resize', resize); resize();
  }
  function resize() {
    const el = $('view'), w = el.clientWidth, h = el.clientHeight;
    renderer.setSize(w, h); camera.aspect = w / Math.max(h, 1); camera.updateProjectionMatrix(); render();
  }
  function render() {
    if (!renderer) return;
    camera.position.set(cam.R * Math.sin(cam.ph) * Math.cos(cam.th), cam.R * Math.cos(cam.ph), cam.R * Math.sin(cam.ph) * Math.sin(cam.th));
    camera.lookAt(0, 0, 0); renderer.render(scene, camera);
  }

  function toShape(sec) {
    const v = p => new THREE.Vector2(p.x, p.y);
    const shape = new THREE.Shape(sec.outer.map(v));
    for (const h of sec.holes) shape.holes.push(new THREE.Path(h.map(v)));
    return shape;
  }

  function showGear(sec, b) {
    if (mesh) {   // 이전 형상의 GPU 자원(기하·재질) 해제
      scene.remove(mesh); mesh.geometry.dispose(); mesh.material.dispose();
      scene.remove(edges); edges.geometry.dispose(); edges.material.dispose();
    }
    const geo = new THREE.ExtrudeGeometry(toShape(sec), { depth: b, bevelEnabled: false, curveSegments: 1 });
    geo.translate(0, 0, -b / 2); geo.rotateX(-Math.PI / 2);          // 기어 축 = 화면 위(+Y)
    geo.computeVertexNormals();
    mesh = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: 0x9aa7b8, metalness: 0.55, roughness: 0.38 }));
    edges = new THREE.LineSegments(new THREE.EdgesGeometry(geo, 25), new THREE.LineBasicMaterial({ color: 0x2b3442 }));
    scene.add(mesh); scene.add(edges);
    cam.R = fit(sec.g.ra * 2, b);
    render();
  }
  let lastSize = null;
  function fit(da, b) { const s = Math.max(da, b); if (lastSize === s) return cam.R; lastSize = s; return s * 1.9; }

  // ---------------------------------------------------------------- 결과 패널
  const f = (v, d = 3) => Number.isFinite(v) ? v.toFixed(d) : '—';
  function update() {
    const { gear, body } = params();
    const out = $('result');
    try {
      const sec = G.section(gear, body), inf = G.info(gear);
      showGear(sec, gear.b);
      const rows = [
        ['피치원 지름 d', f(inf.d), 'mm'], ['이끝원 지름 da', f(inf.da), 'mm'], ['이뿌리원 지름 df', f(inf.df), 'mm'],
        ['기초원 지름 db', f(inf.db), 'mm'], ['피치원 원호 이두께', f(inf.sPitch), 'mm'], ['기초 피치 pb', f(inf.pb), 'mm'],
        ['림 두께 (이뿌리원~구멍)', sec.rim == null ? '—' : f(sec.rim, 2), 'mm'], ['언더컷 한계 잇수(근사)', f(inf.zUndercutApprox, 1), '']
      ];
      out.innerHTML = `<table>${rows.map(r => `<tr><th>${r[0]}</th><td>${r[1]}</td><td>${r[2]}</td></tr>`).join('')}</table>` +
        `<p class="${sec.analyzable ? 'ok' : 'err'}">FEA 해석 대상: ${sec.analyzable ? '가능' : '불가 (아래 형상 문제 해결 필요)'}</p>` +
        (sec.warns.length ? `<ul class="warn">${sec.warns.map(w => `<li>${w}</li>`).join('')}</ul>` : '') +
        '<p class="note">랙 공구로 창성한 외접 스퍼 기어입니다. 해석하려면 키홈 폭과 웨브 구멍 수를 0으로 설정하세요.</p>';
    } catch (e) {
      out.innerHTML = `<p class="err">형상 오류: ${e.message}</p>`;
    }
  }

  function bindFEA(){
    const button=$('run-fea'),output=$('fea-result');
    let worker;
    function clear(){if(worker){worker.terminate();worker=null;}button.disabled=false;output.textContent='';}
    for(const id of ['inputs','mate-teeth','torque'])$(id).addEventListener('input',clear);
    button.addEventListener('click',()=>{
      clear();
      const {gear,body}=params(),mateTeeth=Number($('mate-teeth').value),torqueNm=Number($('torque').value);
      if(!Number.isInteger(mateTeeth)||mateTeeth<3||!(torqueNm>0)){output.textContent='맞물림 잇수와 양의 토크를 입력하세요.';return;}
      button.disabled=true;output.textContent='해석 중…';
      try{worker=new Worker('src/analysis_worker.js');}catch(error){clear();output.textContent='해석 실행에는 HTTP로 연 페이지가 필요합니다.';return;}
      worker.onmessage=({data})=>{
        clear();
        if(data.error){output.textContent=data.error;return;}
        output.textContent='치근 중앙폭 최대주응력 '+f(data.rootMPa,2)+' MPa · 요소 중심 최대 von Mises '+f(data.vmMPa,2)+' MPa. 5치 섹터, 보어 고정, 자유 절단면, E=206 GPa, ν=0.30. 메시 수렴 미확인. 공학적 추정. 시험을 통한 확인 필요.';
      };
      worker.onerror=()=>{clear();output.textContent='해석 실행 환경을 불러오지 못했습니다.';};
      worker.postMessage({gear,body,options:{mateTeeth,torqueNm}});
    });
  }
  window.addEventListener('DOMContentLoaded', () => { buildForm(); bindFEA(); init3D(); update(); });
})();
