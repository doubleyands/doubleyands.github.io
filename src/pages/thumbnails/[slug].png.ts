import type { APIRoute, GetStaticPaths } from 'astro';
import type { CollectionEntry } from 'astro:content';
import { getPosts } from '../../lib/posts';
import { renderThumbnail, thumbnailPath } from '../../lib/thumbnail';

export const getStaticPaths = (async () => {
	const posts = await getPosts();
	return posts
		.filter((post) => thumbnailPath(post))
		.map((post) => ({ params: { slug: post.id }, props: { post } }));
}) satisfies GetStaticPaths;

export const GET: APIRoute = async ({ props }) => {
	const png = await renderThumbnail((props as { post: CollectionEntry<'blog'> }).post);
	if (!png) return new Response(null, { status: 404 });
	return new Response(new Uint8Array(png), { headers: { 'Content-Type': 'image/png' } });
};
