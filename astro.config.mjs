// @ts-check

import { unified } from '@astrojs/markdown-remark';
import sitemap from '@astrojs/sitemap';
import { defineConfig } from 'astro/config';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';

// https://astro.build/config
export default defineConfig({
	site: 'https://doubleyands.github.io',
	// /cv-print/는 cv.pdf를 만드는 데만 쓰고 배포 전에 지운다 (scripts/cv/pdf.mjs)
	integrations: [sitemap({ filter: (page) => !page.includes('/cv-print/') })],
	// CV 페이지는 About으로 합쳤다 (2026-11-06). 잠깐 있던 /cv/ 링크는 About으로 보낸다.
	redirects: {
		'/cv': '/about/',
	},
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
