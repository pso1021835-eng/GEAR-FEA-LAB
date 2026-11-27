// 전체 테스트 실행 (CI)
const { execFileSync } = require('child_process');
const path = require('path');
const tests = ['gear_test.js'];
let failed = 0;
for (const t of tests) {
  try { process.stdout.write(execFileSync(process.execPath, [path.join(__dirname, t)], { encoding: 'utf8' })); }
  catch (e) { failed++; process.stdout.write(e.stdout || ''); console.error(`${t} 실패`); }
}
if (failed) { console.error(`${failed}개 테스트 실패`); process.exit(1); }
console.log('모든 테스트 통과');
