// @ts-check

import { unified } from '@astrojs/markdown-remark';
import sitemap from '@astrojs/sitemap';
import { defineConfig } from 'astro/config';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';

// https://astro.build/config
export default defineConfig({
	site: 'https://doubleyands.github.io',
	integrations: [sitemap()],
	markdown: {
		// 수식을 빌드할 때 KaTeX로 그린다
		processor: unified({
			remarkPlugins: [remarkMath],
			rehypePlugins: [rehypeKatex],
		}),
		// 코드 블록은 어두운 타일 위에 놓인다
		shikiConfig: {
			theme: 'github-dark',
		},
	},
});
