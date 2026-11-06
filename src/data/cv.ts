// Content of the CV page (/cv/), in the order the page shows it.
// Write short, simple sentences. Unpublished work stays out until it is public.
// The About page has its own short list in profile.ts: update both when a paper changes.

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
	updated: 'November 2026',
	headline: 'M.S. student in Computational Science and Engineering, Yonsei University',
	// Shown only when it is filled in
	email: 'youngsooyang@yonsei.ac.kr',
	links: [
		{ label: 'GitHub', url: 'https://github.com/doubleyands' },
		{ label: 'Google Scholar', url: 'https://scholar.google.com/citations?user=PYyqO6AAAAAJ' },
		{ label: 'ORCID', url: 'https://orcid.org/0009-0007-4676-120X' },
	],
	interests: [
		'Physics-informed neural networks (PINNs)',
		'Neural operators',
		'Physics foundation models',
		'Finite element methods',
		'Numerical analysis',
	],

	education: [
		{
			period: 'Sep 2026 – present',
			title: 'M.S. in Computational Science and Engineering',
			place: 'Yonsei University',
			lines: ['Numerical Analysis Lab. Advisor: Prof. Eunjung Lee.'],
		},
		{
			period: 'Mar 2022 – Aug 2026',
			title: 'B.Eng. in Data Science and B.S. in Mathematics',
			place: 'The Catholic University of Korea',
			lines: ['Primary major in Data Science, double major in Mathematics.'],
		},
	] as Entry[],

	research: [
		{
			period: 'Sep 2026 – present',
			title: 'Graduate Researcher',
			place: 'Numerical Analysis Lab, Yonsei University',
			lines: ['Neural network methods for PDEs: PINNs and neural operators.'],
		},
		{
			period: 'Jan 2025 – Aug 2026',
			title: 'Undergraduate Research Intern',
			place: 'Numerical Analysis Lab, Yonsei University',
			lines: [
				'Proposed Physics-Informed Perceiver IO (PIP-IO). It solves PDEs without labeled data and without a mesh. Published as first author in the Journal of Computational Physics.',
				'Cut the attention cost of Transformer-based PINNs from O(N²) to O(mN) with latent cross-attention.',
				'Compared PIP-IO with baselines in a fair way: the same number of parameters and the same input embedding. Reported the error of the final model, not the best error during training.',
				'Co-developed a physics-informed Transformer with slope limiter embedding for PDEs with sharp transitions (under revision).',
			],
		},
	] as Entry[],

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
				'Built with Astro. The home page plays a flow simulation (lattice Boltzmann) as you scroll. The About page shows a PINN while it learns.',
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
