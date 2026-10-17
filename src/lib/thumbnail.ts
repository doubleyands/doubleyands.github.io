import type { CollectionEntry } from 'astro:content';
import sharp from 'sharp';

type Post = CollectionEntry<'blog'>;

// 링크 미리보기 그림의 표준 크기 (1.91 : 1)
export const THUMB_WIDTH = 1200;
export const THUMB_HEIGHT = 630;

const BACKGROUND = '#f5f5f7';
const GRID = '#dcdce1';
const PAD_X = 96;
const PAD_Y = 72;

// 글의 첫 .viz 그래프에서 선과 격자만 남긴다. 그래프가 없으면 undefined.
// 눈금과 라벨 글자는 썸네일 크기에서 읽히지 않으므로 뺀다.
function plotShapes(body: string) {
	const figure = body.match(/<figure class="viz">[\s\S]*?<\/figure>/)?.[0];
	const svg = figure?.match(/<svg[\s\S]*?<\/svg>/)?.[0];
	if (!svg) return undefined;
	return svg
		.replace(/^<svg[^>]*>/, '')
		.replace(/<\/svg>$/, '')
		.replace(/<text[\s\S]*?<\/text>/g, '')
		.replace(/<[a-z]+\b[^>]*(?:visibility="hidden"|class="viz-(?:crosshair|leader)")[^>]*\/>/g, '')
		.replace(/class="viz-grid"/g, `stroke="${GRID}" stroke-width="1"`)
		// 카드 크기로 줄여도 선이 보이도록 데이터 선만 굵게 한다
		.replace(/(<path\b[^>]*?)stroke-width="2"/g, '$1stroke-width="3"');
}

// 남은 도형의 좌표로 그림이 차지하는 범위를 구한다.
function bounds(shapes: string) {
	const xs: number[] = [];
	const ys: number[] = [];
	for (const [, x, y] of shapes.matchAll(/(-?[\d.]+),(-?[\d.]+)/g)) {
		xs.push(Number(x));
		ys.push(Number(y));
	}
	for (const [, axis, value] of shapes.matchAll(/\b(?:c?([xy])[12]?)="(-?[\d.]+)"/g)) {
		(axis === 'x' ? xs : ys).push(Number(value));
	}
	if (xs.length === 0 || ys.length === 0) return undefined;
	const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
	const [y0, y1] = [Math.min(...ys), Math.max(...ys)];
	return { x: x0, y: y0, width: Math.max(x1 - x0, 1), height: Math.max(y1 - y0, 1) };
}

function thumbnailSvg(body: string) {
	const shapes = plotShapes(body);
	if (!shapes) return undefined;
	const box = bounds(shapes);
	if (!box) return undefined;
	const scale = Math.min(
		(THUMB_WIDTH - 2 * PAD_X) / box.width,
		(THUMB_HEIGHT - 2 * PAD_Y) / box.height,
	);
	const dx = (THUMB_WIDTH - box.width * scale) / 2 - box.x * scale;
	const dy = (THUMB_HEIGHT - box.height * scale) / 2 - box.y * scale;
	return `<svg xmlns="http://www.w3.org/2000/svg" width="${THUMB_WIDTH}" height="${THUMB_HEIGHT}" viewBox="0 0 ${THUMB_WIDTH} ${THUMB_HEIGHT}"><rect width="100%" height="100%" fill="${BACKGROUND}"/><g transform="translate(${dx.toFixed(1)} ${dy.toFixed(1)}) scale(${scale.toFixed(4)})">${shapes}</g></svg>`;
}

// 썸네일이 있는 글이면 그 주소, 없으면 undefined
export function thumbnailPath(post: Post) {
	return post.body && plotShapes(post.body) ? `/thumbnails/${post.id}.png` : undefined;
}

export async function renderThumbnail(post: Post) {
	const svg = post.body && thumbnailSvg(post.body);
	if (!svg) return undefined;
	return sharp(Buffer.from(svg)).png({ compressionLevel: 9, palette: true }).toBuffer();
}
