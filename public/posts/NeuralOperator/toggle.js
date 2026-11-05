// On/off widget for the FNO post: pick a layout, click conductors to switch them on (1 V) or off (0 V).
// Every pattern was computed in advance (toggle.py in the project): one finite-difference solve and one model output
// per pattern. toggle.bin holds uint16 planes z = 0, reference then prediction, value = -0.25 + 1.5 q / 65535.
(() => {
	const fig = document.getElementById('onoff');
	if (!fig) return;
	const BASE = '/posts/NeuralOperator/';
	const SEQ = ['#ffffff', '#9ec5f4', '#2a78d6', '#0d366b'];
	const DIV = ['#1c5cab', '#86b6ef', '#f0efec', '#f0a3a2', '#b42b2a'];
	const lut = (stops) => {
		const rgb = stops.map((h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)));
		const out = new Uint8ClampedArray(256 * 3);
		for (let k = 0; k < 256; k++) {
			const t = (k / 255) * (rgb.length - 1), i = Math.min(Math.floor(t), rgb.length - 2), f = t - i;
			for (let c = 0; c < 3; c++) out[3 * k + c] = rgb[i][c] + f * (rgb[i + 1][c] - rgb[i][c]);
		}
		return out;
	};
	const seq = lut(SEQ), div = lut(DIV);
	const el = (tag, cls, text) => {
		const e = document.createElement(tag);
		if (cls) e.className = cls;
		if (text != null) e.textContent = text;
		return e;
	};
	const svgEl = (tag, attrs) => {
		const e = document.createElementNS('http://www.w3.org/2000/svg', tag);
		for (const k in attrs) e.setAttribute(k, attrs[k]);
		return e;
	};
	const pct = (v) => (100 * v).toFixed(2) + '%';
	const val = (q) => -0.25 + (1.5 * q) / 65535;
	const names = (ks) => (ks.length ? ks.join(', ') : 'none');

	let meta, planes, n, offsets, L = 0, mask = 0;
	const alternating = (K) => [...Array(K).keys()].reduce((m, k) => (k % 2 ? m : m | (1 << k)), 0);

	// Layout
	const layoutRow = el('div', 'onoff-row');
	layoutRow.setAttribute('role', 'group');
	layoutRow.setAttribute('aria-label', 'Layout');
	const panels = el('div', 'onoff-panels');
	const frames = ['Reference (FDM)', 'FNO prediction', 'Prediction − reference'].map((title) => {
		const p = el('div');
		const f = el('div', 'onoff-frame');
		const c = el('canvas');
		c.width = c.height = 65;
		c.setAttribute('role', 'img');
		const s = svgEl('svg', { preserveAspectRatio: 'none' });
		f.append(c, s);
		p.append(el('div', 'onoff-title', title), f);
		panels.append(p);
		return { c, s, f };
	});
	const bars = el('div', 'onoff-bars');
	const bar = (stops, ticks, wide) => {
		const w = el('div', wide ? 'onoff-wide' : '');
		const b = el('div', 'onoff-bar');
		b.style.background = `linear-gradient(to right, ${stops.join(', ')})`;
		const t = el('div', 'onoff-ticks');
		ticks.forEach((x) => t.append(el('span', '', x)));
		w.append(b, t);
		bars.append(w);
		return t;
	};
	bar(SEQ, ['0 V', '0.5 V', '1 V'], true);
	const diffTicks = bar(DIV, ['', '0', ''], false);
	const hover = el('div', 'onoff-hover');
	const switchRow = el('div', 'onoff-row');
	switchRow.setAttribute('role', 'group');
	switchRow.setAttribute('aria-label', 'Conductors');
	const readout = el('div', 'onoff-readout', 'Loading the precomputed fields (1.6 MB)…');
	readout.setAttribute('aria-live', 'polite');
	const details = el('details', 'viz-table');
	details.innerHTML = '<summary>Show the charges</summary><table><thead><tr><th>conductor</th><th>state</th>' +
		'<th>charge, solver</th><th>charge, model</th></tr></thead><tbody></tbody></table>';
	const cap = fig.querySelector('figcaption');
	[layoutRow, panels, bars, hover, switchRow, readout, details].forEach((e) => fig.insertBefore(e, cap));
	const HINT = 'Hover over a panel to read the values. Click a conductor, or its button, to switch it.';
	hover.textContent = HINT;

	function toggle(k) {
		if (mask === 1 << k) {
			readout.textContent = 'At least one conductor has to be on (with all off, the potential is 0 everywhere).';
			return;
		}
		mask ^= 1 << k;
		draw();
	}

	function pills(row, labels, pressed, onClick) {
		row.replaceChildren(...labels.map((label, i) => {
			const b = el('button', 'pill pill--sm' + (pressed(i) ? '' : ' pill--ghost'), label);
			b.type = 'button';
			b.setAttribute('aria-pressed', pressed(i));
			b.addEventListener('click', () => onClick(i));
			return b;
		}));
	}

	function paint(canvas, lut, index) {
		const ctx = canvas.getContext('2d'), img = ctx.createImageData(n, n);
		for (let ix = 0; ix < n; ix++) {
			for (let iy = 0; iy < n; iy++) {
				const k = 3 * index(ix * n + iy), p = 4 * ((n - 1 - iy) * n + ix);
				img.data[p] = lut[k];
				img.data[p + 1] = lut[k + 1];
				img.data[p + 2] = lut[k + 2];
				img.data[p + 3] = 255;
			}
		}
		ctx.putImageData(img, 0, 0);
	}

	function draw() {
		const lay = meta.layouts[L], K = lay.K, N2 = n * n, pat = lay.patterns[mask - 1];
		const idx = offsets[L] + 2 * (mask - 1);
		const ref = planes.subarray(idx * N2, (idx + 1) * N2), pred = planes.subarray((idx + 1) * N2, (idx + 2) * N2);
		let lim = 0;
		for (let k = 0; k < N2; k++) lim = Math.max(lim, Math.abs(val(pred[k]) - val(ref[k])));
		lim = lim || 1e-9;
		const clamp = (v, a, b) => Math.min(Math.max(v, a), b);
		const field = (a) => (k) => Math.round(clamp(val(a[k]), 0, 1) * 255);
		paint(frames[0].c, seq, field(ref));
		paint(frames[1].c, seq, field(pred));
		paint(frames[2].c, div, (k) => Math.round(((clamp((val(pred[k]) - val(ref[k])) / lim, -1, 1) + 1) / 2) * 255));
		diffTicks.children[0].textContent = '−' + lim.toFixed(3) + ' V';
		diffTicks.children[2].textContent = '+' + lim.toFixed(3) + ' V';
		const on = [], off = [];
		for (let k = 0; k < K; k++) (mask >> k & 1 ? on : off).push(k + 1);
		// Conductor outlines (clickable) and numbers on all three panels
		const h = 1 / (n - 1), lo = -0.5 - h / 2, size = 1 + h;
		frames.forEach(({ c, s }, panel) => {
			s.setAttribute('viewBox', `${lo} ${lo} ${size} ${size}`);
			s.replaceChildren();
			lay.boxes.forEach(([x0, x1, y0, y1], k) => {
				const isOn = mask >> k & 1;
				const r = svgEl('rect', { x: x0 * h - 0.5, y: 0.5 - y1 * h, width: (x1 - x0) * h, height: (y1 - y0) * h, class: 'onoff-box' });
				r.addEventListener('click', () => toggle(k));
				const t = svgEl('text', { x: ((x0 + x1) / 2) * h - 0.5, y: 0.5 - ((y0 + y1) / 2) * h, class: 'onoff-num' });
				t.setAttribute('fill', isOn && panel < 2 ? '#ffffff' : '#1d1d1f');
				t.textContent = k + 1;
				s.append(r, t);
			});
			c.setAttribute('aria-label', `${['Reference', 'Prediction', 'Difference'][panel]} on the plane z = 0, conductors ${names(on)} on`);
		});
		pills(layoutRow, meta.layouts.map((l) => l.name), (i) => i === L, (i) => {
			L = i;
			mask = alternating(meta.layouts[i].K);
			draw();
		});
		pills(switchRow, lay.boxes.map((_, k) => `${k + 1} ${mask >> k & 1 ? 'on' : 'off'}`), (k) => mask >> k & 1, toggle);
		readout.textContent = `Conductors ${names(on)} on (1 V), ${names(off)} off (0 V). ` +
			`Field error ${pct(pat.err)} over the whole box, charge error ${pct(pat.qerr)}.`;
		details.querySelector('tbody').replaceChildren(...pat.Q.map((q, k) => {
			const tr = el('tr');
			[k + 1, mask >> k & 1 ? 'on' : 'off', q.toFixed(4), pat.Qp[k].toFixed(4)].forEach((v) => tr.append(el('td', '', v)));
			return tr;
		}));
		fig._read = (ix, iy) => [val(ref[ix * n + iy]), val(pred[ix * n + iy])];
	}

	frames.forEach(({ f }) => {
		f.addEventListener('mousemove', (e) => {
			if (!fig._read) return;
			const b = f.getBoundingClientRect();
			const ix = Math.min(n - 1, Math.max(0, Math.floor(((e.clientX - b.left) / b.width) * n)));
			const iy = n - 1 - Math.min(n - 1, Math.max(0, Math.floor(((e.clientY - b.top) / b.height) * n)));
			const [r, p] = fig._read(ix, iy), d = p - r;
			hover.textContent = `x = ${(ix / (n - 1) - 0.5).toFixed(3)}, y = ${(iy / (n - 1) - 0.5).toFixed(3)}: ` +
				`solver ${r.toFixed(4)} V, model ${p.toFixed(4)} V, difference ${d >= 0 ? '+' : '−'}${Math.abs(d).toFixed(4)} V`;
		});
		f.addEventListener('mouseleave', () => (hover.textContent = HINT));
	});

	let started = false;
	const load = () => {
		if (started) return;
		started = true;
		Promise.all([
			fetch(BASE + 'toggle.json').then((r) => r.json()),
			fetch(BASE + 'toggle.bin').then((r) => r.arrayBuffer()),
		]).then(([m, buf]) => {
			meta = m;
			n = m.n;
			planes = new Uint16Array(buf);
			offsets = [];
			let o = 0;
			for (const l of m.layouts) {
				offsets.push(o);
				o += 2 * l.patterns.length;
			}
			mask = alternating(m.layouts[0].K);
			draw();
		}).catch(() => (readout.textContent = 'Could not load the data for this figure.'));
	};
	if ('IntersectionObserver' in window) {
		new IntersectionObserver((es, ob) => {
			if (es.some((e) => e.isIntersecting)) {
				ob.disconnect();
				load();
			}
		}, { rootMargin: '600px' }).observe(fig);
	} else load();
})();
