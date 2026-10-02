// Compares the simulated drag, lift and shedding frequency with the reference
// ranges of the DFG 2D-2 benchmark (Schäfer & Turek, 1996).
//
// Usage: node scripts/flow/validate.mjs [--in=scripts/flow/out] [--from=8]

import fs from 'node:fs';
import path from 'node:path';

const args = Object.fromEntries(
	process.argv.slice(2).map((arg) => arg.replace(/^--/, '').split('=')),
);
const IN = args.in ?? 'scripts/flow/out';
const T_FROM = Number(args.from ?? 8);

const rows = fs
	.readFileSync(path.join(IN, 'forces.csv'), 'utf8')
	.trim()
	.split('\n')
	.slice(1)
	.map((line) => line.split(',').map(Number))
	.filter(([t]) => t >= T_FROM);

const cd = rows.map((row) => row[1]);
const cl = rows.map((row) => row[2]);
const mean = (values) => values.reduce((sum, value) => sum + value, 0) / values.length;

// Shedding period from upward zero crossings of the lift about its mean
const clMean = mean(cl);
const crossings = [];
for (let i = 1; i < rows.length; i++) {
	const a = cl[i - 1] - clMean;
	const b = cl[i] - clMean;
	if (a < 0 && b >= 0) crossings.push(rows[i - 1][0] + ((rows[i][0] - rows[i - 1][0]) * -a) / (b - a));
}
const period = (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1);
const strouhal = 0.1 / period / 1.0; // St = f D / U_mean

const report = [
	['max drag coefficient', Math.max(...cd), '3.22 - 3.24'],
	['max lift coefficient', Math.max(...cl), '0.99 - 1.01'],
	['min lift coefficient', Math.min(...cl), ''],
	['Strouhal number', strouhal, '0.295 - 0.305'],
];
console.log(`t >= ${T_FROM} s, ${crossings.length - 1} shedding periods, period ${period.toFixed(4)} s`);
for (const [name, value, reference] of report) {
	console.log(`${name.padEnd(22)} ${value.toFixed(4).padStart(8)}   reference ${reference}`);
}
