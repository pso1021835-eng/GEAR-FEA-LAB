// gear_test.js — 창성 치형이 이론 인벌류트·치수와 맞는지
const G = require('../src/gear.js');
let fail = 0;
const ok = (c, msg) => { if (!c) { fail++; console.log('FAIL', msg); } };

const cases = [
  { m: 2, z: 30, alpha: 20, x: 0 },
  { m: 1.5, z: 18, alpha: 20, x: 0.3 },
  { m: 3, z: 45, alpha: 25, x: -0.2, rho: 0.25 },
  { m: 1, z: 60, alpha: 20, x: 0.5, rho: 0.25 }
];
for (const p of cases) {
  const tag = JSON.stringify(p);
  const { g, pts, undercut } = G.halfSpace(p);
  ok(!undercut, tag + ' undercut 오판');
  const R = pts.map(q => Math.hypot(q.x, q.y));
  ok(Math.abs(R[0] - g.rf) < 1e-9, tag + ` rf ${R[0]} vs ${g.rf}`);
  ok(Math.abs(R[R.length - 1] - g.ra) < 1e-9, tag + ' ra');
  // 인벌류트 구간: 이론 반 이두께 각과 비교 (이 중심 = 90° − π/z)
  let maxErr = 0, n = 0;
  for (const q of pts) {
    if (q.part !== 'involute') continue;
    const Rq = Math.hypot(q.x, q.y); if (Rq < g.rb * 1.001) continue;
    const half = Math.atan2(q.y, q.x) - (Math.PI / 2 - Math.PI / g.z);
    maxErr = Math.max(maxErr, Math.abs(half - G.involuteHalfAngle(p, Rq)) * Rq * 1000); n++;   // µm 원호 오차
  }
  ok(n > 50 && maxErr < 0.01, tag + ` 인벌류트 오차 ${maxErr.toFixed(4)} µm (n=${n})`);
  // 필렛은 이뿌리원과 인벌류트 사이에서 연속(점 간 간격이 m의 5% 이하)
  let gap = 0; for (let i = 1; i < pts.length; i++) gap = Math.max(gap, Math.hypot(pts[i].x - pts[i - 1].x, pts[i].y - pts[i - 1].y));
  ok(gap < 0.05 * p.m, tag + ` 경계 불연속 ${gap}`);
  // 전체 외곽: 점 수, 닫힘, 면적 > 0(반시계)
  const o = G.outline(p).points;
  let A = 0; for (let i = 0; i < o.length; i++) { const a = o[i], b = o[(i + 1) % o.length]; A += a.x * b.y - b.x * a.y; }
  ok(A > 0, tag + ' 외곽 방향');
  const Aring = Math.PI * (g.ra ** 2 - g.rf ** 2), Aroot = Math.PI * g.rf ** 2;
  ok(A / 2 > Aroot && A / 2 < Aroot + Aring, tag + ' 외곽 면적 범위');
}
// ρf* 한계 초과는 명확한 오류
try { G.halfSpace({ m: 3, z: 45, alpha: 25 }); ok(false, 'ρ 한계 미검출'); } catch (e) { ok(/최대/.test(e.message), 'ρ 한계 메시지'); }
// 언더컷: z = 8, x = 0 은 언더컷이 나야 함
ok(G.halfSpace({ m: 2, z: 8, x: 0 }).undercut, 'z=8 언더컷 미검출');
ok(G.info({ m: 2, z: 8 }).zUndercutApprox > 8, '언더컷 한계 잇수 근사');


// 단면: 보어·키홈·웨브 구멍, 경고
{
  const area = P => { let A = 0; for (let i = 0; i < P.length; i++) { const a = P[i], b = P[(i + 1) % P.length]; A += a.x * b.y - b.x * a.y; } return A / 2; };
  const s1 = G.section({ m: 2, z: 40 }, { bore: 20, key: { w: 6, t: 2.8 }, web: { n: 4, pcd: 46, d: 10 } });
  ok(s1.holes.length === 5, '구멍 수');
  ok(s1.holes.every(h => area(h) < 0), '구멍은 시계 방향');
  const kb = Math.PI * 100 + 6 * 2.8;   // 보어 면적 + 키홈 근사
  ok(Math.abs(-area(s1.holes[0]) - kb) / kb < 0.02, '보어+키홈 면적 ' + (-area(s1.holes[0])).toFixed(2));
  ok(s1.outer.length === 40 * (2 * (3 + 12 + 20) + 5) || s1.outer.length > 40 * 60, '재샘플 점 수 ' + s1.outer.length);
  ok(s1.warns.length === 0, '정상 형상 경고 없음: ' + s1.warns.join('/'));
  const s2 = G.section({ m: 2, z: 40 }, { bore: 20, web: { n: 4, pcd: 60, d: 14 } });
  ok(s2.warns.some(w => /림|이뿌리원/.test(w)), '얇은 림 경고');
  const s3 = G.section({ m: 2, z: 8 }, {});
  ok(s3.warns.some(w => /언더컷/.test(w)), '언더컷 경고');
}

// F1-REV 회귀: 구멍 현 간격, 키홈 폭 오류, 공구 둥글기 한계 중복점, 재샘플 전이점 보존, 해석 가능 플래그
{
  const w = G.section({ m: 2, z: 30 }, { web: { n: 4, pcd: 20, d: 15 } });
  ok(w.warns.some(x => /웨브 구멍끼리/.test(x)) && !w.analyzable, '웨브 구멍 현 간격 겹침 검출');
  try { G.boreHole(10, { w: 11, t: 2 }); ok(false, '키홈 폭 오류 미검출'); } catch (e) { ok(/키홈 폭/.test(e.message), '키홈 폭 메시지'); }
  const a = 20 * Math.PI / 180, rhoMax = (Math.PI / 4 - 1.25 * Math.tan(a)) / (1 / Math.cos(a) - Math.tan(a));
  const hs = G.halfSpace({ m: 2, z: 30, rho: rhoMax * (1 - 1e-12) }).pts;
  let dup = 0; for (let i = 1; i < hs.length; i++) if (Math.hypot(hs[i].x - hs[i - 1].x, hs[i].y - hs[i - 1].y) < 1e-9) dup++;
  ok(dup === 0, '둥글기 한계 중복점 ' + dup);
  const raw = G.halfSpace({ m: 2, z: 30 }).pts, rs = G.resampleHalf(raw, { root: 3, fillet: 12, involute: 20 });
  const lastOf = part => { const q = raw.filter(p => p.part === part); return q[q.length - 1]; };
  for (const part of ['root', 'fillet']) { const t = lastOf(part); ok(rs.some(p => Math.hypot(p.x - t.x, p.y - t.y) < 1e-9), part + ' 전이점 보존'); }
  ok(rs.length === 1 + 3 + 12 + 20, '재샘플 점 수 ' + rs.length);
  ok(!G.section({ m: 2, z: 8 }).analyzable, '언더컷은 해석 불가');
  ok(G.section({ m: 2, z: 30 }, { bore: 20, key: { w: 6, t: 2.8 } }).analyzable, '정상 형상 해석 가능');
}

console.log(fail ? `gear_test: ${fail} 실패` : 'gear_test: 모두 통과');
process.exit(fail ? 1 : 0);
