// Prints the CV page (dist/cv-print/) to dist/cv.pdf with headless Chrome,
// then removes the page from dist so only the PDF is published.
// Run it after `astro build`: `npm run cv:pdf`. The deploy workflow runs it.
// Chrome is found in the usual places. Set CHROME_PATH to use another one.

import { spawn } from 'node:child_process';
import { existsSync } from 'node:fs';
import { mkdtemp, readFile, rm, stat } from 'node:fs/promises';
import { createServer } from 'node:http';
import { tmpdir } from 'node:os';
import { extname, join, normalize, resolve, sep } from 'node:path';

const DIST = resolve('dist');
const PAGE_DIR = join(DIST, 'cv-print');
const OUT = join(DIST, 'cv.pdf');

const TYPES = {
	'.html': 'text/html; charset=utf-8',
	'.css': 'text/css',
	'.js': 'text/javascript',
	'.svg': 'image/svg+xml',
	'.png': 'image/png',
	'.woff2': 'font/woff2',
};

function findChrome() {
	if (process.env.CHROME_PATH) return process.env.CHROME_PATH;
	const env = process.env;
	const candidates =
		{
			win32: [
				`${env.PROGRAMFILES}\\Google\\Chrome\\Application\\chrome.exe`,
				`${env['PROGRAMFILES(X86)']}\\Google\\Chrome\\Application\\chrome.exe`,
				`${env.LOCALAPPDATA}\\Google\\Chrome\\Application\\chrome.exe`,
				`${env['PROGRAMFILES(X86)']}\\Microsoft\\Edge\\Application\\msedge.exe`,
			],
			darwin: ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'],
			linux: ['/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/usr/bin/chromium-browser'],
		}[process.platform] ?? [];
	return candidates.find((path) => existsSync(path));
}

function fail(message) {
	console.error(`cv:pdf: ${message}`);
	process.exit(1);
}

if (!existsSync(join(PAGE_DIR, 'index.html'))) fail('dist/cv-print/index.html is missing. Run `npm run build` first.');
const chrome = findChrome();
if (!chrome) fail('Chrome was not found. Set CHROME_PATH.');

// The page uses absolute paths, so serve dist over HTTP instead of opening the file.
const server = createServer(async (req, res) => {
	let path = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
	if (path.endsWith('/')) path += 'index.html';
	const file = normalize(join(DIST, path));
	if (!file.startsWith(DIST + sep)) return res.writeHead(403).end();
	try {
		const body = await readFile(file);
		res.writeHead(200, { 'content-type': TYPES[extname(file)] ?? 'application/octet-stream' }).end(body);
	} catch {
		res.writeHead(404).end();
	}
});
await new Promise((done) => server.listen(0, '127.0.0.1', done));
const url = `http://127.0.0.1:${server.address().port}/cv-print/`;

const profileDir = await mkdtemp(join(tmpdir(), 'cv-chrome-'));
const args = [
	'--headless',
	'--disable-gpu',
	'--no-first-run',
	'--no-default-browser-check',
	`--user-data-dir=${profileDir}`,
	'--no-pdf-header-footer',
	// Wait for the web font before printing
	'--virtual-time-budget=15000',
	`--print-to-pdf=${OUT}`,
	url,
];
// The CI runner (Ubuntu) blocks Chrome's sandbox. The page is our own, so turn it off there.
if (process.platform === 'linux') args.unshift('--no-sandbox');

await rm(OUT, { force: true });
const code = await new Promise((done) => {
	const child = spawn(chrome, args, { stdio: ['ignore', 'ignore', 'pipe'] });
	let stderr = '';
	child.stderr.on('data', (chunk) => (stderr += chunk));
	const timer = setTimeout(() => child.kill(), 90_000);
	child.on('close', (exitCode) => {
		clearTimeout(timer);
		if (exitCode !== 0) console.error(stderr);
		done(exitCode);
	});
	child.on('error', (error) => {
		clearTimeout(timer);
		console.error(error);
		done(1);
	});
});
server.close();
await rm(profileDir, { recursive: true, force: true }).catch(() => {});

if (code !== 0 || !existsSync(OUT)) fail(`Chrome did not write the PDF (exit code ${code}).`);
const { size } = await stat(OUT);
if (size < 10_000) fail(`the PDF is too small (${size} bytes).`);
const pages = ((await readFile(OUT, 'latin1')).match(/\/Type\s*\/Page(?!s)/g) ?? []).length;

await rm(PAGE_DIR, { recursive: true, force: true });
console.log(`cv:pdf: wrote dist/cv.pdf (${pages || '?'} pages, ${Math.round(size / 1024)} KB)`);
