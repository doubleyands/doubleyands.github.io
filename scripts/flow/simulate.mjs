// Flow past a cylinder in a channel: the DFG 2D-2 benchmark setup
// (Schäfer & Turek, 1996), Re = 100.
//
// Method: D2Q9 lattice Boltzmann, TRT collision, halfway bounce-back on the
// channel walls and the cylinder, non-equilibrium extrapolation at the inlet
// (velocity) and the outlet (pressure).
//
// One thing differs from the benchmark: the inflow is ramped up smoothly from
// rest over RAMP seconds, so the start-up can be shown. After the ramp the
// setup is the benchmark's.
//
// Usage: node scripts/flow/simulate.mjs [--t=10] [--snap=0.02] [--out=scripts/flow/out]
// Writes velocity snapshots (Float32 u then v, row-major) and forces.csv.

import fs from 'node:fs';
import path from 'node:path';

const args = Object.fromEntries(
	process.argv.slice(2).map((arg) => arg.replace(/^--/, '').split('=')),
);
const T_END = Number(args.t ?? 10);
const SNAP_DT = Number(args.snap ?? 0.02);
const OUT = args.out ?? 'scripts/flow/out';

// Physical setup (SI units, as in the benchmark)
const H = 0.41;
const L = 2.2;
const D = 0.1;
const CX = 0.2;
const CY = 0.2;
const U_MAX = 1.5;
const U_MEAN = (2 / 3) * U_MAX;
const RE = 100;
const RAMP = 1.0;

// Lattice
const DX = 0.0025;
const NX = Math.round(L / DX);
const NY = Math.round(H / DX);
const N = NX * NY;
const D_LAT = D / DX;
const U_MEAN_LAT = 0.05;
const U_MAX_LAT = 1.5 * U_MEAN_LAT;
const NU_LAT = (U_MEAN_LAT * D_LAT) / RE;
const TAU = 3 * NU_LAT + 0.5;
const DT = (DX * U_MEAN_LAT) / U_MEAN;

// TRT relaxation rates, magic parameter 3/16
const W_PLUS = 1 / TAU;
const W_MINUS = 1 / (0.5 + 3 / 16 / (TAU - 0.5));

const EX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const EY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const WT = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// Node (x, y) sits at ((x + 0.5) DX, (y + 0.5) DX), so the walls lie half a cell outside
const solid = new Uint8Array(N);
for (let y = 0; y < NY; y++) {
	for (let x = 0; x < NX; x++) {
		const dx = (x + 0.5) * DX - CX;
		const dy = (y + 0.5) * DX - CY;
		if (dx * dx + dy * dy < (D / 2) ** 2) solid[x + y * NX] = 1;
	}
}

// Parabolic inflow profile in lattice units
const inflow = new Float64Array(NY);
for (let y = 0; y < NY; y++) {
	const yy = (y + 0.5) * DX;
	inflow[y] = (4 * U_MAX_LAT * yy * (H - yy)) / (H * H);
}

let f = new Float64Array(9 * N);
let g = new Float64Array(9 * N);
for (let i = 0; i < 9; i++) f.fill(WT[i], i * N, (i + 1) * N);

const feq = new Float64Array(9);
function equilibrium(rho, ux, uy) {
	const usq = 1.5 * (ux * ux + uy * uy);
	for (let i = 0; i < 9; i++) {
		const cu = 3 * (EX[i] * ux + EY[i] * uy);
		feq[i] = WT[i] * rho * (1 + cu + 0.5 * cu * cu - usq);
	}
}

const fs_ = new Float64Array(9); // post-collision populations of one node
let forceX = 0;
let forceY = 0;

function step(time) {
	forceX = 0;
	forceY = 0;

	// Collide, then push to the neighbours
	for (let y = 0; y < NY; y++) {
		for (let x = 0; x < NX; x++) {
			const idx = x + y * NX;
			if (solid[idx]) continue;

			let rho = 0;
			let ux = 0;
			let uy = 0;
			for (let i = 0; i < 9; i++) {
				const fi = f[i * N + idx];
				rho += fi;
				ux += fi * EX[i];
				uy += fi * EY[i];
			}
			ux /= rho;
			uy /= rho;
			equilibrium(rho, ux, uy);

			fs_[0] = f[idx] - W_PLUS * (f[idx] - feq[0]);
			for (let a = 1; a < 9; a++) {
				const b = OPP[a];
				if (b < a) continue;
				const fa = f[a * N + idx];
				const fb = f[b * N + idx];
				const sym = W_PLUS * (0.5 * (fa + fb) - 0.5 * (feq[a] + feq[b]));
				const asym = W_MINUS * (0.5 * (fa - fb) - 0.5 * (feq[a] - feq[b]));
				fs_[a] = fa - sym - asym;
				fs_[b] = fb - sym + asym;
			}

			g[idx] = fs_[0];
			for (let i = 1; i < 9; i++) {
				const xn = x + EX[i];
				const yn = y + EY[i];
				if (xn < 0 || xn >= NX) continue; // inlet and outlet columns are set below
				if (yn < 0 || yn >= NY) {
					g[OPP[i] * N + idx] = fs_[i]; // channel wall
					continue;
				}
				const nidx = xn + yn * NX;
				if (solid[nidx]) {
					g[OPP[i] * N + idx] = fs_[i]; // cylinder
					forceX += 2 * fs_[i] * EX[i];
					forceY += 2 * fs_[i] * EY[i];
				} else {
					g[i * N + nidx] = fs_[i];
				}
			}
		}
	}

	// Inlet: prescribed velocity, density from the neighbour column
	const ramp = time >= RAMP ? 1 : 0.5 * (1 - Math.cos((Math.PI * time) / RAMP));
	for (let y = 0; y < NY; y++) {
		setBoundaryNode(0 + y * NX, 1 + y * NX, null, inflow[y] * ramp);
	}
	// Outlet: prescribed density, velocity from the neighbour column
	for (let y = 0; y < NY; y++) {
		setBoundaryNode(NX - 1 + y * NX, NX - 2 + y * NX, 1, null);
	}

	[f, g] = [g, f];
}

// Non-equilibrium extrapolation (Guo et al., 2002)
const neq = new Float64Array(9);
function setBoundaryNode(idx, nidx, rhoFixed, uxFixed) {
	let rho = 0;
	let ux = 0;
	let uy = 0;
	for (let i = 0; i < 9; i++) {
		const gi = g[i * N + nidx];
		rho += gi;
		ux += gi * EX[i];
		uy += gi * EY[i];
	}
	ux /= rho;
	uy /= rho;
	equilibrium(rho, ux, uy);
	for (let i = 0; i < 9; i++) neq[i] = g[i * N + nidx] - feq[i];

	if (uxFixed === null) equilibrium(rhoFixed, ux, uy);
	else equilibrium(rho, uxFixed, 0);
	for (let i = 0; i < 9; i++) g[i * N + idx] = feq[i] + neq[i];
}

function velocity() {
	const out = new Float32Array(2 * N);
	const scale = U_MEAN / U_MEAN_LAT; // lattice velocity -> m/s
	for (let idx = 0; idx < N; idx++) {
		if (solid[idx]) continue;
		let rho = 0;
		let ux = 0;
		let uy = 0;
		for (let i = 0; i < 9; i++) {
			const fi = f[i * N + idx];
			rho += fi;
			ux += fi * EX[i];
			uy += fi * EY[i];
		}
		out[idx] = (ux / rho) * scale;
		out[N + idx] = (uy / rho) * scale;
	}
	return out;
}

fs.mkdirSync(OUT, { recursive: true });
const steps = Math.round(T_END / DT);
const snapEvery = Math.round(SNAP_DT / DT);
const forceEvery = 20;
const forces = ['t,cd,cl'];
const coeff = 2 / (U_MEAN_LAT * U_MEAN_LAT * D_LAT);
let snapshots = 0;

fs.writeFileSync(
	path.join(OUT, 'meta.json'),
	JSON.stringify({ nx: NX, ny: NY, dx: DX, dt: DT, snapDt: snapEvery * DT, tau: TAU, re: RE, ramp: RAMP }, null, 2),
);
console.log(`grid ${NX}x${NY}, tau ${TAU.toFixed(3)}, dt ${DT} s, ${steps} steps`);

const started = Date.now();
for (let n = 0; n <= steps; n++) {
	const time = n * DT;
	if (n % snapEvery === 0) {
		const name = `uv_${String(snapshots).padStart(4, '0')}.f32`;
		fs.writeFileSync(path.join(OUT, name), Buffer.from(velocity().buffer));
		snapshots++;
	}
	step(time);
	if (n % forceEvery === 0) {
		forces.push(`${time.toFixed(5)},${(forceX * coeff).toFixed(5)},${(forceY * coeff).toFixed(5)}`);
	}
	if (n % 4000 === 0) {
		const elapsed = (Date.now() - started) / 1000;
		if (!Number.isFinite(forceX)) throw new Error(`diverged at step ${n}`);
		console.log(`t=${time.toFixed(2)} s  step ${n}/${steps}  ${elapsed.toFixed(0)} s elapsed`);
		fs.writeFileSync(path.join(OUT, 'forces.csv'), forces.join('\n'));
	}
}
fs.writeFileSync(path.join(OUT, 'forces.csv'), forces.join('\n'));
console.log(`done: ${snapshots} snapshots`);
