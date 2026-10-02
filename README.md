# doubleyands.github.io

실험 기록 블로그. <https://doubleyands.github.io>

[Astro](https://astro.build)로 만들고, `main`에 푸시하면 GitHub Actions가 빌드해서 GitHub Pages에 올린다.

## 처음 받은 뒤

```
npm install
```

설치할 때 커밋 검사(`.githooks/pre-commit`)가 함께 켜진다.

## 글 쓰기

`src/content/blog/YYYY-MM-DD-제목.md` 파일을 만든다.

```markdown
---
title: "글 제목"
description: "목록과 글 머리에 보이는 한 줄 설명"
pubDate: 2026-10-03
project: 프로젝트이름
tags: [태그1, 태그2]
---

본문
```

- 글과 화면의 글귀는 영어로 쓴다. 짧고 쉬운 문장으로.
- `title`, `description`, `pubDate`, `project`는 필수다.
- 그림은 `src/assets/posts/<프로젝트이름>/`에 두고 본문에서 상대 경로로 넣는다: `![설명](../../assets/posts/프로젝트이름/그림.png)`
- 수식은 `$ ... $`(문장 안)과 `$$ ... $$`(별도 줄)로 쓴다.

## 소개 고치기

이름, 소속, 연구 분야, 학력, 논문은 `src/data/profile.ts`에 있다. 첫 화면과 `/about/`이 이 파일을 읽는다. 논문을 추가하려면 `publications`에 항목을 하나 더 넣는다.

## 공개 허락

커밋할 때 각 글의 `project`가 공개를 허락받은 프로젝트인지 확인하고, 아니면 커밋을 막는다. 허락 여부는 이 저장소 밖(`~/Desktop/ME/projects/<프로젝트이름>/README.md` 머리말의 `blog: true`)에 있다. 확인할 수 없을 때도 막는다.

## 로컬에서 보기

```
npm run dev
```

## 디자인

`src/styles/global.css` 맨 위의 변수(색, 반경, 간격)가 기준이다.

- 강조색은 파란색 하나만 쓴다. 누를 수 있는 것은 모두 이 색이다.
- 구역은 흰색, 옅은 회색, 어두운 타일을 번갈아 놓아 나눈다. 테두리나 그림자로 나누지 않는다.
- 주요 버튼은 알약 모양, 카드는 18px 반경에 얇은 테두리만 둔다.
- 글자 굵기는 300, 400, 600, 700만 쓴다.
