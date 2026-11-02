# Style rules for chapters of the ATTPCROOTv2 LaTeX manual

- File: `docs/manual/chapters/<name>.tex`, starting with `\chapter{...}` and `\label{ch:<name>}`. No preamble, no
  `\begin{document}`. The shared preamble (`docs/manual/preamble.tex`) provides: amsmath, graphicx, booktabs, array,
  url, hyperref, and the macros `\code{}` (inline code/class/file names; escape `_` as `\_`), `\brho`, `\dedx`,
  `\sqdedx`. Use nothing else (no colours, no tcolorbox, no listings, no siunitx, no custom fonts).
- Plain default LaTeX look: `\section`, `\subsection`, `\paragraph`; `itemize`/`enumerate`; tables with booktabs
  (`\toprule`, `\midrule`, `\bottomrule`) and `p{}` columns wide enough that nothing overflows the margin; code examples
  in the `verbatim` environment (keep lines under ~95 characters); file paths with `\url{...}` when long.
- Audience: physicists using and extending the framework. Explain WHAT each component does, HOW to configure it
  (tasks, setters with defaults, branch names in/out), and WHY a choice matters; give a short real example taken from
  a macro in the repository (cite its path); list known pitfalls.
- Every statement must be backed by the code or by existing docs in the repository. Read the headers and sources;
  quote defaults from the code. Do not invent classes, setters, defaults or behaviour. If something is unclear, say so
  in the text ("not documented in the code") rather than guessing.
- Write in prose, plain technical English, active voice, no marketing language, no "simply"/"just".
- Before finishing, compile the chapter alone with `docs/manual/test_chapter.sh <name>` and fix every error and every
  `Overfull \hbox` larger than 5pt.
