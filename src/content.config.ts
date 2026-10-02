import { defineCollection } from 'astro:content';
import { glob } from 'astro/loaders';
import { z } from 'astro/zod';

const blog = defineCollection({
	loader: glob({ base: './src/content/blog', pattern: '**/*.md' }),
	schema: z.object({
		title: z.string(),
		description: z.string(),
		pubDate: z.coerce.date(),
		updatedDate: z.coerce.date().optional(),
		// 이 글이 속한 프로젝트. 커밋할 때 공개 허락 여부를 이 값으로 확인한다.
		project: z.string().regex(/^[A-Za-z0-9._-]+$/),
		tags: z.array(z.string()).default([]),
	}),
});

export const collections = { blog };
