export const profile = {
	name: 'Youngsoo Yang',
	role: 'M.S. student in Computational Science and Engineering at Yonsei University.',
	affiliation: 'Department of Computational Science and Engineering, Yonsei University',
	intro: [
		'I work on neural-network methods for solving partial differential equations.',
		'This site is where I keep notes on my experiments.',
	],
	interests: [
		{
			name: 'Physics-Informed Neural Networks (PINNs)',
			description: 'Neural networks trained to satisfy the governing equations of a physical system.',
		},
		{
			name: 'Neural Operators',
			description: 'Models that learn maps between function spaces, such as from the inputs of a PDE to its solution.',
		},
		{
			name: 'Numerical Analysis',
			description: 'The design and analysis of algorithms that solve mathematical problems numerically.',
		},
	],
	education: [
		{
			degree: 'M.S. student in Computational Science and Engineering',
			school: 'Yonsei University',
			period: '2026 – present',
			note: '',
		},
		{
			degree: 'B.Eng. in Data Science and B.S. in Mathematics',
			school: 'The Catholic University of Korea',
			period: '2022 – 2026',
			note: 'Primary major in Data Science, double major in Mathematics',
		},
	],
	publications: [
		{
			title:
				'Physics-Informed Perceiver IO: scalable mesh-free solving of partial differential equations without supervision',
			authors: 'Y. Yang, E. Lee',
			venue: 'Journal of Computational Physics, 2026',
			url: 'https://doi.org/10.1016/j.jcp.2026.115455',
		},
	],
};
