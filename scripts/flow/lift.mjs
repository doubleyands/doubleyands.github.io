// Samples the lift coefficient from simulate.mjs at the times of the home page
// frames, and finds the two moments the page marks on its lift chart:
// when vortex shedding starts, and when enough steady cycles have passed to
// read off the Strouhal number.
//
// Usage: node scripts/flow/lift.mjs [--from=0] [--to=4.5] [--frames=271]
//          [--in=scripts/flow/out] [--out=src/data/flow-lift.json]
//
// Use the same --from, --to and --frames as render.mjs.

import fs from 'node:fs';
import path from 'node:path';

const args = Object.fromEntries(
	process.argv.slice(2).map((arg) => arg.replace(/^--/, '').split('=')),
);
const IN = args.in ?? 'scripts/flow/out';
const OUT = args.out ?? 'src/data/flow-lift.json';
const T_FROM = Number(args.from ?? 0);
const T_TO = Number(args.to ?? 4.5);
const FRAMES = Number(args.frames ?? 271);

const D = 0.1; // cylinder diameter
const U = 1.0; // mean inflow speed, as in the DFG definition St = f D / U
const SHEDDING = 0.1; // |C_L| above this counts as the start of shedding
const STEADY = 0.3; // a cycle whose |C_L| peak is above this counts as steady
const CYCLES = 2; // steady cycles to wait for before showing St

const rows = fs
	.readFileSync(path.join(IN, 'forces.csv'), 'utf8')
	.trim()
	.split('\n')
	.slice(1)
	.map((line) => line.split(',').map(Number));
const ts = rows.map((row) => row[0]);
const cls = rows.map((row) => row[2]);

// Linear interpolation of C_L at time t
const liftAt = (t) => {
	let hi = ts.findIndex((time) => time >= t);
	if (hi <= 0) return cls[Math.max(hi, 0)];
	const lo = hi - 1;
	const w = (t - ts[lo]) / (ts[hi] - ts[lo]);
	return (1 - w) * cls[lo] + w * cls[hi];
};

const times = Array.from({ length: FRAMES }, (_, i) => T_FROM + ((T_TO - T_FROM) * i) / (FRAMES - 1));
const cl = times.map((t) => Number(liftAt(t).toFixed(3)));

// Shedding starts where the lift first leaves the band around zero
const start = ts.findIndex((t, i) => t >= T_FROM && Math.abs(cls[i]) > SHEDDING);
const sheddingStart = ts[start];

// Upward zero crossings, interpolated between samples
const ups = [];
for (let i = 1; i < ts.length && ts[i] <= T_TO; i++) {
	if (cls[i - 1] < 0 && cls[i] >= 0) {
		ups.push(ts[i - 1] + ((ts[i] - ts[i - 1]) * -cls[i - 1]) / (cls[i] - cls[i - 1]));
	}
}

// The first crossing that opens a steady cycle, then CYCLES cycles after it.
// St uses only those cycles, so the number shown never relies on later data.
const peak = (a, b) => Math.max(...cls.filter((_, i) => ts[i] >= a && ts[i] <= b).map(Math.abs));
const first = ups.findIndex((t, k) => k + 1 < ups.length && peak(t, ups[k + 1]) > STEADY);
if (first < 0 || first + CYCLES >= ups.length) {
	throw new Error('Not enough steady cycles inside the window to measure St');
}
const period = (ups[first + CYCLES] - ups[first]) / CYCLES;
const strouhal = {
	value: Number(((D / period) / U).toFixed(3)),
	from: Number(ups[first].toFixed(3)),
	shownAt: Number(ups[first + CYCLES].toFixed(3)),
	cycles: CYCLES,
};

const result = {
	from: T_FROM,
	to: T_TO,
	frames: FRAMES,
	sheddingStart: Number(sheddingStart.toFixed(3)),
	strouhal,
	peak: Number(peak(T_FROM, T_TO).toFixed(3)),
	cl,
};
fs.writeFileSync(OUT, JSON.stringify(result) + '\n');
console.log(
	`shedding starts at t = ${result.sheddingStart} s; St = ${strouhal.value} ` +
		`from ${CYCLES} cycles, ${strouhal.from}-${strouhal.shownAt} s; peak |C_L| = ${result.peak}`,
);
