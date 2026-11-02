#!/bin/bash
# test_chapter.sh <name>: compile docs/manual/chapters/<name>.tex alone with the shared preamble; prints errors and overfull boxes
set -u
D=$(cd "$(dirname "$0")" && pwd); N=$1; T=$(mktemp -d)
cat > $T/t.tex <<EOT
\documentclass[11pt,a4paper]{report}
\input{$D/preamble}
\begin{document}
\input{$D/chapters/$N}
\end{document}
EOT
# run from docs/manual so figures/... paths resolve
cd $D && pdflatex -interaction=nonstopmode -output-directory=$T $T/t.tex >/dev/null 2>&1; pdflatex -interaction=nonstopmode -output-directory=$T $T/t.tex >/dev/null 2>&1
grep -E "^!|Overfull|Undefined control" -A2 $T/t.log | head -40
echo "pages: $(grep -o 'Output written.*(\([0-9]*\) page' $T/t.log | grep -o '[0-9]* page')"
rm -rf $T
