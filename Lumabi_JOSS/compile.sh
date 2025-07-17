pandoc paper.md \
  --from markdown \
  --to pdf \
  --output paper.pdf \
  --citeproc \
  --pdf-engine=xelatex
