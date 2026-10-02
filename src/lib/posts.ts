import { type CollectionEntry, getCollection } from 'astro:content';

type Post = CollectionEntry<'blog'>;

export async function getPosts() {
	const posts = await getCollection('blog');
	return posts.sort((a, b) => b.data.pubDate.valueOf() - a.data.pubDate.valueOf());
}

// posts 는 최신순이어야 한다. 처음 만나는 글이 그 프로젝트의 최근 글이다.
export function getProjects(posts: Post[]) {
	const projects = new Map<string, { name: string; count: number; latest: Date }>();
	for (const post of posts) {
		const entry = projects.get(post.data.project);
		if (entry) entry.count += 1;
		else
			projects.set(post.data.project, {
				name: post.data.project,
				count: 1,
				latest: post.data.pubDate,
			});
	}
	return [...projects.values()];
}

const seoulDate = new Intl.DateTimeFormat('sv-SE', { timeZone: 'Asia/Seoul' });

// 2026-10-03
export function isoDate(date: Date) {
	return seoulDate.format(date);
}
