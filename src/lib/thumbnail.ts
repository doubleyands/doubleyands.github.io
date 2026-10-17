import type { CollectionEntry } from 'astro:content';
import sharp from 'sharp';

type Post = CollectionEntry<'blog'>;

// 링크 미리보기 그림의 표준 크기 (1.91 : 1)
export const THUMB_WIDTH = 1200;
export const THUMB_HEIGHT = 630;

// 글마다 손으로 정한 주제 그림. 파일 이름이 글 id와 같다.
// 그림 원본은 scripts/thumbnails/make.py 가 만든다.
const sources = import.meta.glob<string>('../thumbnails/*.svg', {
	query: '?raw',
	import: 'default',
	eager: true,
});

function source(post: Post) {
	return sources[`../thumbnails/${post.id}.svg`];
}

// 썸네일이 있는 글이면 그 주소, 없으면 undefined
export function thumbnailPath(post: Post) {
	return source(post) ? `/thumbnails/${post.id}.png` : undefined;
}

export async function renderThumbnail(post: Post) {
	const svg = source(post);
	if (!svg) return undefined;
	return sharp(Buffer.from(svg)).png({ compressionLevel: 9, palette: true }).toBuffer();
}
