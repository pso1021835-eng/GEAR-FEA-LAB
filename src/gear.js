/* gear.js — 스퍼 기어 2D 단면 형상 (랙 공구 창성: 인벌류트 + 트로코이드 필렛)
 * 단위: mm, rad. 브라우저(window.GearGeom)·Node(module.exports) 공용.
 *
 * 창성 모델: 기어 중심 = 원점, 기어 피치원(r)이 랙의 구름선 y = r 위를 미끄럼 없이 구른다.
 * 랙 기준선 y = r + x·m (양의 전위 = 공구가 바깥으로), 공구 이끝 y = rf = r + x·m − hf*·m,
 * 공구 이끝 둥글기 ρ = ρf*·m. 랙 위의 점 P(법선 n)가 기어를 깎는 순간은 법선이 피치점 I=(0,r)을 지날 때:
 *   (P + s·ex − I) × n = 0  →  s = (Py − r)·nx/ny − Px,   기어 좌표 = Rot(+s/r)·(Px + s, Py)
 */
(function (root) {
  'use strict';
  const deg = Math.PI / 180;
  const inv = a => Math.tan(a) - a;
  const rot = (p, a) => [p[0] * Math.cos(a) - p[1] * Math.sin(a), p[0] * Math.sin(a) + p[1] * Math.cos(a)];

  function defaults(p) {
    return Object.assign({ m: 2, z: 30, alpha: 20, x: 0, ha: 1, hf: 1.25, rho: 0.38, b: 20 }, p);
  }

  function basic(p) {
    const g = defaults(p), a = g.alpha * deg, r = g.m * g.z / 2;
    return {
      ...g, alphaRad: a, r, rb: r * Math.cos(a),
      ra: r + g.m * (g.ha + g.x), rf: r + g.m * (g.x - g.hf),
      psi: (Math.PI / 2 + 2 * g.x * Math.tan(a)) / g.z   // 피치원 반 이두께 각
    };
  }

  // 랙 공구 오른쪽 옆면의 프로파일(이끝선 → 이끝 둥글기 → 직선 옆면), 각 점과 바깥 법선
  function rackRight(g, nArc, nFlank) {
    const m = g.m, a = g.alphaRad, yRef = g.r + g.x * m, yTip = yRef - g.hf * m, rho = g.rho * m;
    const xFl = y => Math.PI * m / 4 - (yRef - y) * Math.tan(a);
    const yc = yTip + rho, cx = xFl(yc) - rho / Math.cos(a);
    if (cx < 0) {
      const rhoMax = (Math.PI / 4 - g.hf * Math.tan(a)) / (1 / Math.cos(a) - Math.tan(a));
      throw new Error(`공구 이끝 둥글기 ρf* ${g.rho}가 너무 큽니다 (이 압력각·이뿌리 높이에서 최대 ${rhoMax.toFixed(3)})`);
    }
    const pts = [];
    // 이끝선: x = 0 … cx  (법선 0,−1)
    const nRoot = 12;
    for (let i = 0; i <= nRoot; i++) pts.push({ p: [cx * i / nRoot, yTip], n: [0, -1], part: 'root' });
    // 이끝 둥글기: θ = −90° … −α
    for (let i = 1; i <= nArc; i++) {
      const t = -Math.PI / 2 + (Math.PI / 2 - a) * i / nArc, c = Math.cos(t), s = Math.sin(t);
      pts.push({ p: [cx + rho * c, yc + rho * s], n: [c, s], part: 'fillet' });
    }
    // 직선 옆면: 접점 → 기어 이끝원을 넘을 때까지
    const t0 = pts[pts.length - 1].p, n = [Math.cos(a), -Math.sin(a)], dir = [Math.sin(a), Math.cos(a)];
    const len = (g.ha + g.x + g.hf + 1) * m / Math.cos(a);
    for (let i = 1; i <= nFlank; i++) {
      const l = len * i / nFlank;
      pts.push({ p: [t0[0] + dir[0] * l, t0[1] + dir[1] * l], n, part: 'involute' });
    }
    return pts;
  }

  function envelope(g, q) {
    const [px, py] = q.p, [nx, ny] = q.n;
    const s = (py - g.r) * nx / ny - px;
    return rot([px + s, py], s / g.r);
  }

  /** 반 이홈(이홈 중심 = +y축) 오른쪽 경계: 이뿌리원 → 필렛 → 인벌류트 → 이끝원 */
  function halfSpace(p, opt = {}) {
    const g = basic(p);
    const raw = rackRight(g, opt.nArc || 200, opt.nFlank || 400).map(q => ({ ...q, g: envelope(g, q) }));
    const out = [];
    let undercut = false, lastR = -Infinity;
    for (const q of raw) {
      const R = Math.hypot(q.g[0], q.g[1]);
      if (R < lastR - 1e-9) undercut = true;
      if (R > g.ra) {   // 이끝원에서 정확히 자르기
        const prev = out[out.length - 1], Rp = Math.hypot(prev.g[0], prev.g[1]), t = (g.ra - Rp) / (R - Rp);
        const x = prev.g[0] + t * (q.g[0] - prev.g[0]), y = prev.g[1] + t * (q.g[1] - prev.g[1]), k = g.ra / Math.hypot(x, y);
        out.push({ g: [x * k, y * k], part: 'involute' });
        break;
      }
      lastR = Math.max(lastR, R);
      out.push(q);
    }
    // 연속 중복점 제거 (공구 이끝 둥글기 한계에서 이끝선 길이 0 → 같은 점 반복)
    const pts = [];
    for (const q of out) { const l = pts[pts.length - 1]; if (!l || Math.hypot(q.g[0] - l.x, q.g[1] - l.y) > 1e-9 * g.m) pts.push({ x: q.g[0], y: q.g[1], part: q.part }); }
    return { g, pts, undercut };
  }

  /** 반 이홈 점열을 부위별 호 길이 등간격으로 줄임 (표시·메시용). counts: {root, fillet, involute} */
  function resampleHalf(pts, counts) {
    const out = [];
    for (const part of ['root', 'fillet', 'involute']) {
      const own = pts.filter(q => q.part === part);
      if (!own.length) continue;
      const seg = out.length ? [out[out.length - 1], ...own] : own;   // 앞 부위 끝점(=전이 접점)에서 시작
      if (seg.length < 2) { out.push(...own); continue; }
      const L = [0]; for (let i = 1; i < seg.length; i++) L.push(L[i - 1] + Math.hypot(seg[i].x - seg[i - 1].x, seg[i].y - seg[i - 1].y));
      const n = Math.max(1, counts[part]), tot = L[L.length - 1];
      let j = 0;
      for (let k = 0; k <= n; k++) {
        if (out.length && k === 0) continue;   // 부위 경계점 중복 방지
        const s = tot * k / n;
        while (j < seg.length - 2 && L[j + 1] < s) j++;
        const t = (s - L[j]) / Math.max(L[j + 1] - L[j], 1e-15);
        out.push({ x: seg[j].x + t * (seg[j + 1].x - seg[j].x), y: seg[j].y + t * (seg[j + 1].y - seg[j].y), part });
      }
    }
    return out;
  }

  /** 전체 기어 외곽(반시계, 닫힌 다각형) — 이 하나당: 이홈 오른쪽 반 → 이끝 호 → 다음 이홈 왼쪽 반
   *  opt.counts 지정 시 부위별 점 수로 재샘플(예: {root:3, fillet:12, involute:20}) */
  function outline(p, opt = {}) {
    const hs = halfSpace(p, opt), g = hs.g, undercut = hs.undercut;
    const pts = opt.counts ? resampleHalf(hs.pts, opt.counts) : hs.pts;
    const z = g.z, pitch = 2 * Math.PI / z;
    // 이홈 중심을 각 0 기준으로: 반 이홈 점들의 각 = atan2 − 90°
    const polar = pts.map(q => ({ R: Math.hypot(q.x, q.y), t: Math.atan2(q.y, q.x) - Math.PI / 2, part: q.part }));
    // 오른쪽 반 이홈은 t < 0 (시계 방향 쪽). 반시계로 돌기 위해: 왼쪽 반(거울: −t)은 뿌리→끝, 이끝 호, 다음 이홈 오른쪽 반은 끝→뿌리
    const right = polar;                                   // t ≤ 0, 뿌리(0) → 끝(−ta)
    const left = polar.map(q => ({ ...q, t: -q.t }));      // t ≥ 0
    const tipAngle = -right[right.length - 1].t;           // 이끝 점의 이홈 중심 기준 각
    const nTip = opt.nTip || 8, poly = [];
    for (let k = 0; k < z; k++) {
      const c = k * pitch;
      for (const q of left) poly.push({ ...q, t: c + q.t });                         // 이홈 왼쪽: 뿌리 → 끝
      const t1 = c + tipAngle, t2 = c + pitch - tipAngle;                             // 이끝 호
      for (let i = 1; i < nTip; i++) poly.push({ R: g.ra, t: t1 + (t2 - t1) * i / nTip, part: 'tip' });
      for (let i = right.length - 1; i >= 0; i--) poly.push({ ...right[i], t: c + pitch + right[i].t });  // 다음 이홈 오른쪽: 끝 → 뿌리
    }
    // 이홈 중심의 중복점 제거
    const clean = poly.filter((q, i) => { const n = poly[(i + 1) % poly.length]; return Math.abs(q.R - n.R) + Math.abs(q.t - n.t) > 1e-9; });
    return { g, undercut, points: clean.map(q => ({ x: q.R * Math.cos(q.t + Math.PI / 2), y: q.R * Math.sin(q.t + Math.PI / 2), part: q.part })) };
  }

  /** 해석 기준값 (검증·화면 표시용) */
  function info(p) {
    const g = basic(p);
    const zmin = 2 * (g.hf - g.rho * (1 - Math.sin(g.alphaRad)) - g.x) / Math.sin(g.alphaRad) ** 2; // 이뿌리 둥글기 고려 언더컷 한계(근사)
    return {
      m: g.m, z: g.z, alpha: g.alpha, x: g.x, d: 2 * g.r, db: 2 * g.rb, da: 2 * g.ra, df: 2 * g.rf,
      sPitch: 2 * g.r * g.psi,                        // 피치원 원호 이두께 [mm]
      pb: Math.PI * g.m * Math.cos(g.alphaRad),       // 기초 피치
      zUndercutApprox: zmin
    };
  }

  /** 반경 R에서 인벌류트 이론 반 이두께 각 (검증용) */
  function involuteHalfAngle(p, R) {
    const g = basic(p);
    return g.psi + inv(g.alphaRad) - inv(Math.acos(g.rb / R));
  }

  /** 원 다각형 (cw=true면 시계 — 구멍용) */
  function circle(cx, cy, R, n, cw) {
    const pts = [];
    for (let i = 0; i < n; i++) { const t = (cw ? -1 : 1) * 2 * Math.PI * i / n; pts.push({ x: cx + R * Math.cos(t), y: cy + R * Math.sin(t) }); }
    return pts;
  }

  /** 보어 + 키홈(보어 위쪽 +y 방향) 구멍 경계, 시계 방향
   *  bore: 보어 지름, key: { w: 키홈 폭, t: 보어 표면에서 키홈 깊이 } (KS B 1311 등 규격값은 사용자가 입력) */
  function boreHole(bore, key, n = 96) {
    const R = bore / 2;
    if (!key || !(key.w > 0) || !(key.t > 0)) return circle(0, 0, R, n, true);
    if (key.w >= bore) throw new Error(`키홈 폭 ${key.w} mm가 보어 지름 ${bore} mm 이상입니다`);
    const hw = key.w / 2, a0 = Math.asin(Math.min(1, hw / R)), top = Math.sqrt(R * R - hw * hw) + key.t;
    // 시계 방향: 키홈 오른쪽 아래 모서리(보어와 만나는 점)에서 출발해 보어 호를 시계로 돌아 왼쪽 모서리 → 키홈 위로
    const pts = [], t1 = Math.PI / 2 - a0, t2 = -(3 * Math.PI / 2 - a0);   // 오른쪽 모서리 각 → 왼쪽 모서리 각(시계로 감소)
    const m = Math.max(16, Math.round(n * (t1 - t2) / (2 * Math.PI)));
    for (let i = 0; i <= m; i++) { const t = t1 + (t2 - t1) * i / m; pts.push({ x: R * Math.cos(t), y: R * Math.sin(t) }); }
    pts.push({ x: -hw, y: top }, { x: hw, y: top });
    return pts;
  }

  /** 기어 단면 전체: 외곽 + 구멍(보어·키홈·웨브 경량 구멍) + 형상 경고
   *  body: { bore, key:{w,t}, web:{n, pcd, d} } [mm] */
  function section(p, body = {}, opt = {}) {
    const o = outline(p, { counts: opt.counts || { root: 3, fillet: 12, involute: 20 }, nTip: opt.nTip || 6 });
    const g = o.g, holes = [], warns = [];
    if (o.undercut) warns.push(`언더컷 발생 (z ${g.z} < 한계 약 ${info(p).zUndercutApprox.toFixed(1)}) — 이뿌리가 깎여 강도가 낮아집니다. 양의 전위를 검토하세요`);
    const h = (g.ha + g.hf) * g.m;              // 전체 이 높이
    let inner = 0;                               // 구멍이 차지하는 최대 반경
    if (body.bore > 0) {
      holes.push(boreHole(body.bore, body.key));
      inner = body.bore / 2 + (body.key && body.key.t > 0 ? body.key.t : 0);
    }
    if (body.web && body.web.n > 0 && body.web.d > 0) {
      const { n, pcd, d } = body.web;
      for (let i = 0; i < n; i++) { const t = Math.PI / n + 2 * Math.PI * i / n; holes.push(circle(pcd / 2 * Math.cos(t), pcd / 2 * Math.sin(t), d / 2, 48, true)); }
      if (pcd / 2 - d / 2 <= (body.bore || 0) / 2 + (body.key && body.key.t || 0)) warns.push('웨브 구멍이 보어·키홈과 겹칩니다');
      if (n > 1 && 2 * (pcd / 2) * Math.sin(Math.PI / n) <= d) warns.push('웨브 구멍끼리 겹칩니다 (형상 불가)');
      inner = Math.max(inner, pcd / 2 + d / 2);
    }
    const rim = g.rf - inner;                    // 림 두께 (이뿌리원 ~ 가장 바깥 구멍)
    if (inner > 0 && rim <= 0) warns.push('구멍이 이뿌리원을 넘습니다 — 형상 불가');
    else if (inner > 0 && rim < 1.2 * h) warns.push(`림 두께 ${rim.toFixed(2)} mm < 1.2 × 이 높이(${(1.2 * h).toFixed(2)} mm) — 얇은 림은 치근 응력을 키웁니다 (ISO 6336-3 YB 영역, 기준값은 원문 확인 필요)`);
    const blocking = o.undercut || warns.some(w => /형상 불가|겹칩니다/.test(w));
    return { g, outer: o.points, holes, warns, rim: inner > 0 ? rim : null, undercut: o.undercut,
      analyzable: !blocking };   // false면 F2 FEA에 넘기지 않음 (언더컷 경계는 포락선 교점 재구성 전까지 해석 대상 아님)
  }

  // Exact rack-envelope boundary radius at an angle measured from a tooth centre.
  function boundaryRadius(p,theta){
    const g=basic(p),pitch=2*Math.PI/g.z;
    let q=theta+pitch/2;q=((q+pitch/2)%pitch+pitch)%pitch-pitch/2;
    const delta=Math.abs(q),rack=rackRight(g,1,1),corner=rack.find(v=>v.part==='fillet');
    const rho=g.rho*g.m,yc=g.rf+rho,cx=corner.p[0]-rho*Math.cos(-g.alphaRad);
    const arc=t=>envelope(g,{p:[cx+rho*Math.cos(t),yc+rho*Math.sin(t)],n:[Math.cos(t),Math.sin(t)]});
    const angle=v=>Math.abs(Math.atan2(v[0],v[1])),start=arc(-Math.PI/2),end=arc(-g.alphaRad);
    if(delta<=angle(start))return g.rf;
    if(delta<angle(end)){
      let lo=-Math.PI/2,hi=-g.alphaRad;
      for(let i=0;i<50;i++){const mid=(lo+hi)/2;if(angle(arc(mid))<delta)lo=mid;else hi=mid;}
      return Math.hypot(...arc((lo+hi)/2));
    }
    const tip=pitch/2-involuteHalfAngle(g,g.ra);if(delta>=tip)return g.ra;
    let lo=Math.max(g.rb,Math.hypot(...end)),hi=g.ra;
    for(let i=0;i<50;i++){const mid=(lo+hi)/2;if(pitch/2-involuteHalfAngle(g,mid)<delta)lo=mid;else hi=mid;}return (lo+hi)/2;
  }
  const api = { basic, halfSpace, resampleHalf, outline, info, involuteHalfAngle, boreHole, section, boundaryRadius };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.GearGeom = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
