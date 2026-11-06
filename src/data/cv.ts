// Content of the About page (/about/), in the order the page shows it.
// Write short, simple sentences. Unpublished work stays out until it is public.
// The name and the short intro come from profile.ts.

export interface Entry {
	period: string;
	title: string;
	place?: string;
	url?: string;
	lines?: string[];
}

export interface Publication {
	authors: string[];
	title: string;
	venue: string;
	url?: string;
}

export const ME = 'Y. Yang';

export const cv = {
	// The profile card on the left
	position: ['M.S. student', 'Computational Science and Engineering', 'Yonsei University'],
	email: 'youngsooyang@yonsei.ac.kr',
	links: [
		{ label: 'GitHub', icon: 'github', url: 'https://github.com/doubleyands' },
		{ label: 'Google Scholar', icon: 'scholar', url: 'https://scholar.google.com/citations?user=PYyqO6AAAAAJ' },
		{ label: 'ORCID', icon: 'orcid', url: 'https://orcid.org/0009-0007-4676-120X' },
	],

	// Each card opens the posts on that topic. The pictures are post thumbnails.
	interests: [
		{
			name: 'Physics-informed neural networks',
			description: 'Neural networks trained to satisfy the governing equations of a physical system.',
			image: '/thumbnails/2026-11-06-pinn-01-decay.png',
			url: '/projects/PINN/',
		},
		{
			name: 'Neural operators',
			description: 'Models that learn maps between function spaces, such as from the inputs of a PDE to its solution.',
			image: '/thumbnails/2026-11-05-neural-operator-fno.png',
			url: '/projects/NeuralOperator/',
		},
		{
			name: 'Numerical analysis',
			description: 'The design and analysis of algorithms that solve mathematical problems numerically.',
			image: '/thumbnails/2026-10-03-numerical-analysis-week-1.png',
			url: '/projects/numerical-analysis/',
		},
	],
	moreInterests: 'Also: physics foundation models and finite element methods.',

	// Newest first. Only things that are already public.
	news: [
		{
			date: 'Nov 2026',
			text: 'Started a series on this site: one PINN setup on 30 problems.',
			url: '/projects/PINN/',
			label: 'Posts',
		},
		{
			date: 'Nov 2026',
			text: 'New post: a Fourier neural operator for 3D capacitance.',
			url: '/blog/2026-11-05-neural-operator-fno/',
			label: 'Read',
		},
		{
			date: 'Oct 2026',
			text: 'Physics-Informed Perceiver IO was published in the Journal of Computational Physics.',
			url: 'https://doi.org/10.1016/j.jcp.2026.115455',
			label: 'Paper',
		},
		{ date: 'Sep 2026', text: 'Started my M.S. at Yonsei University.', url: '', label: '' },
		{
			date: 'Jul 2026',
			text: 'Preprint on arXiv: hierarchical collocation point refinement for PINNs.',
			url: 'https://arxiv.org/abs/2607.14665',
			label: 'arXiv',
		},
	],

	publications: [
		{
			heading: 'Journal articles',
			items: [
				{
					authors: ['Y. Yang', 'E. Lee*'],
					title:
						'Physics-Informed Perceiver IO: scalable mesh-free solving of partial differential equations without supervision',
					venue: 'Journal of Computational Physics 569 (2027) 115455',
					url: 'https://doi.org/10.1016/j.jcp.2026.115455',
				},
			],
		},
		{
			heading: 'Under revision',
			items: [
				{
					authors: ['D. Shin', 'M. Choi', 'Y. Yang', 'E. Lee*'],
					title:
						'Physics-Informed Transformer Networks with Slope Limiter Embedding for PDEs with Sharp Transitions',
					venue: 'Engineering with Computers, under revision (2026)',
				},
			],
		},
		{
			heading: 'Preprints',
			items: [
				{
					authors: ['M. Choi', 'D. Shin', 'Y. Yang', 'E. Lee'],
					title:
						'A new strategy for physics-informed neural networks based on hierarchical collocation point refinement',
					venue: 'arXiv:2607.14665 [math.NA] (2026)',
					url: 'https://arxiv.org/abs/2607.14665',
				},
			],
		},
	] as { heading: string; items: Publication[] }[],

	research: [
		{
			period: 'Sep 2026 – present',
			title: 'Graduate Researcher',
			place: 'Applied Mathematics and Computing (AMC) Lab, Yonsei University',
			lines: ['Neural network methods for PDEs: PINNs and neural operators.'],
		},
		{
			period: 'Jan 2025 – Aug 2026',
			title: 'Undergraduate Research Intern',
			place: 'Applied Mathematics and Computing (AMC) Lab, Yonsei University',
			lines: [
				'Proposed Physics-Informed Perceiver IO (PIP-IO). It solves PDEs without labeled data and without a mesh. Published as first author in the Journal of Computational Physics.',
				'Cut the attention cost of Transformer-based PINNs from O(N²) to O(mN) with latent cross-attention.',
				'Compared PIP-IO with baselines in a fair way: the same number of parameters and the same input embedding. Reported the error of the final model, not the best error during training.',
				'Co-developed a physics-informed Transformer with slope limiter embedding for PDEs with sharp transitions (under revision).',
			],
		},
	] as Entry[],

	projects: [
		{
			period: 'Nov 2026 – present',
			title: 'PINN benchmark on 30 problems',
			url: '/projects/PINN/',
			lines: [
				'One PINN setup on 30 problems, from simple ODEs to 2D PDEs.',
				'Three ways to weight the loss (NTK, RBA, BRDR), 5 seeds each.',
				'One blog post for each problem, with the code and the results.',
			],
		},
		{
			period: 'Nov 2026',
			title: 'Neural operator for 3D capacitance',
			url: '/projects/NeuralOperator/',
			lines: [
				'Generated 3D Laplace data with my own finite-difference solver.',
				'Trained a Fourier neural operator that maps a conductor layout to the electric potential, and computed the capacitance matrix from its output.',
				'Mean field error 0.3%, mean error of the total capacitance 0.3%.',
			],
		},
		{
			period: 'Oct 2026 – present',
			title: 'This website',
			url: '/',
			lines: [
				'Built with Astro. The home page plays a flow simulation (lattice Boltzmann) as you scroll.',
				'Weekly notes on numerical analysis and numerical PDEs.',
			],
		},
		{
			period: 'Mar 2025 – Dec 2025',
			title: 'PINN guided by FEM for electric field analysis',
			place: 'Capstone Design, Department of Data Science, The Catholic University of Korea',
			lines: [
				'Used coarse-grid solutions from my CUDA C++ FEM solver to guide PINN training.',
				'Got stable solutions where plain PINNs did not converge. Checked the results with ANSYS.',
				'Lower relative L2 error than plain PINNs, SA-PINNs, and NTK-weighted PINNs.',
			],
		},
	] as Entry[],

	education: [
		{
			period: 'Sep 2026 – present',
			title: 'M.S. in Computational Science and Engineering',
			place: 'Yonsei University',
			lines: ['Applied Mathematics and Computing (AMC) Lab. Advisor: Prof. Eunjung Lee.'],
		},
		{
			period: 'Mar 2022 – Aug 2026',
			title: 'B.Eng. in Data Science and B.S. in Mathematics',
			place: 'The Catholic University of Korea',
			lines: ['Primary major in Data Science, double major in Mathematics.'],
		},
	] as Entry[],

	awards: [
		{
			period: '2025',
			title: 'Excellence Award, AI Mathematics Competition',
			place: 'Department of Mathematics, The Catholic University of Korea',
			lines: [
				'Non-convex optimization of the Rastrigin function: Monte Carlo and differential evolution on the GPU in CUDA C++.',
			],
		},
		{
			period: '2025',
			title: 'Excellence Award, Mathematics Academic Festival',
			place: 'Department of Mathematics, The Catholic University of Korea',
			lines: ['A PINN for electric potential, with FEM solutions as a training guide.'],
		},
		{
			period: '2024',
			title: 'First Prize, Mathematics Academic Festival',
			place: 'Department of Mathematics, The Catholic University of Korea',
			lines: ['Adaptive sampling for PINNs, based on the idea of adaptive mesh refinement.'],
		},
		{
			period: '2023',
			title: 'First Prize, Mathematics Academic Festival',
			place: 'Department of Mathematics, The Catholic University of Korea',
			lines: ['A Python tool that solves university math problems.'],
		},
	] as Entry[],

	skills: [
		{ label: 'Programming', value: 'Python, CUDA C++' },
		{ label: 'Libraries and tools', value: 'PyTorch, NumPy, SciPy, Matplotlib, ANSYS, Git, LaTeX' },
		{
			label: 'Methods',
			value: 'PINNs, neural operators, Transformers and attention, FEM, finite differences, GPU parallel computing',
		},
		{ label: 'Languages', value: 'Korean (native), English (OPIc IH, 2026)' },
	],

	military: [
		{
			period: 'Oct 2019 – May 2021',
			title: 'Republic of Korea Army',
			lines: ['Sergeant. Completed mandatory military service.'],
		},
	] as Entry[],
};
