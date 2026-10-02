// Turns the velocity snapshots from simulate.mjs into the image sequence the
// home page plays on scroll.
//
// Shown quantity: vorticity. Counter-clockwise rotation is blue, clockwise is
// white, and fluid that does not rotate is black.
//
// Usage: node scripts/flow/render.mjs [--from=0] [--to=6] [--frames=120]
//          [--in=scripts/flow/out] [--out=public/flow] [--scale=40]

import fs from 'node:fs';
import path from 'node:path';
import sharp from 'sharp';

const args = Object.fromEntries(
	process.argv.slice(2).map((arg) => arg.replace(/^--/, '').split('=')),
);
const IN = args.in ?? 'scripts/flow/out';
const OUT = args.out ?? 'public/flow';
const T_FROM = Number(args.from ?? 0);
const T_TO = Number(args.to ?? 6);
const FRAMES = Number(args.frames ?? 120);
const SCALE = Number(args.scale ?? 40); // vorticity (1/s) at which the color is ~76% saturated
const GAMMA = Number(args.gamma ?? 1.5); // above 1 dims the weak shear along the walls

const meta = JSON.parse(fs.readFileSync(path.join(IN, 'meta.json'), 'utf8'));
const { nx: NX, ny: NY, dx: DX, snapDt } = meta;
const N = NX * NY;

const H = 0.41;
const L = 2.2;
const D = 0.1;
const CX = 0.2;
const CY = 0.2;

const SIZES = [
	{ name: 'lg', width: NX * 2, height: NY * 2, quality: 70 },
	{ name: 'sm', width: NX, height: NY, quality: 68 },
];
const BLUE = [41, 151, 255]; // counter-clockwise
const WHITE = [245, 245, 247]; // clockwise

function vorticity(uv) {
	const u = uv.subarray(0, N);
	const v = uv.subarray(N);
	const w = new Float32Array(N);
	for (let y = 0; y < NY; y++) {
		for (let x = 0; x < NX; x++) {
			const idx = x + y * NX;
			// a mirrored ghost value puts u = 0 on the wall, half a cell away
			const uUp = y + 1 < NY ? u[idx + NX] : -u[idx];
			const uDown = y > 0 ? u[idx - NX] : -u[idx];
			const xr = Math.min(x + 1, NX - 1);
			const xl = Math.max(x - 1, 0);
			const dvdx = (v[xr + y * NX] - v[xl + y * NX]) / ((xr - xl) * DX);
			const dudy = (uUp - uDown) / (2 * DX);
			w[idx] = dvdx - dudy;
		}
	}
	return w;
}

function toPixels(w) {
	const rgb = Buffer.alloc(N * 3);
	for (let y = 0; y < NY; y++) {
		for (let x = 0; x < NX; x++) {
			const value = w[x + y * NX];
			const strength = Math.tanh((Math.abs(value) / SCALE) ** GAMMA);
			const color = value > 0 ? BLUE : WHITE;
			const p = (x + (NY - 1 - y) * NX) * 3; // image rows run top to bottom
			rgb[p] = color[0] * strength;
			rgb[p + 1] = color[1] * strength;
			rgb[p + 2] = color[2] * strength;
		}
	}
	return rgb;
}

function cylinder(width, height) {
	const cx = (CX / L) * width;
	const cy = ((H - CY) / H) * height;
	const r = (D / 2 / L) * width + width / NX; // one cell wider, to cover the staircase edge
	return Buffer.from(
		`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}">` +
			`<circle cx="${cx}" cy="${cy}" r="${r}" fill="#2a2a2c"/></svg>`,
	);
}

for (const size of SIZES) fs.mkdirSync(path.join(OUT, size.name), { recursive: true });

const times = [];
let bytes = { lg: 0, sm: 0 };
for (let k = 0; k < FRAMES; k++) {
	const time = T_FROM + ((T_TO - T_FROM) * k) / (FRAMES - 1);
	const snapshot = Math.round(time / snapDt);
	const file = path.join(IN, `uv_${String(snapshot).padStart(4, '0')}.f32`);
	const raw = fs.readFileSync(file);
	const uv = new Float32Array(raw.buffer, raw.byteOffset, 2 * N);
	const pixels = toPixels(vorticity(uv));
	times.push(Number((snapshot * snapDt).toFixed(3)));

	for (const size of SIZES) {
		const target = path.join(OUT, size.name, `${String(k).padStart(3, '0')}.webp`);
		const info = await sharp(pixels, { raw: { width: NX, height: NY, channels: 3 } })
			.resize(size.width, size.height, { kernel: 'cubic' })
			.composite([{ input: cylinder(size.width, size.height) }])
			.webp({ quality: size.quality, effort: 5 })
			.toFile(target);
		bytes[size.name] += info.size;
	}
}

// The page reads the frame count and the time of each frame from here
fs.writeFileSync(
	'src/data/flow.json',
	JSON.stringify({ frames: FRAMES, width: NX, height: NY, times }, null, '\t') + '\n',
);
console.log(
	`${FRAMES} frames, t = ${times[0]} .. ${times[times.length - 1]} s; ` +
		`lg ${(bytes.lg / 1e6).toFixed(2)} MB, sm ${(bytes.sm / 1e6).toFixed(2)} MB`,
);
