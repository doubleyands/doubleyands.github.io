// A view-transition-name shared by a card and the page it opens.
// It must be a valid CSS identifier and unique on each page.
export const transitionName = (kind: 'post' | 'project', id: string) =>
	`view-transition-name: ${kind}-${id.replace(/[^A-Za-z0-9-]/g, '-')}; view-transition-class: title`;
