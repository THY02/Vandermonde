# 论文与源码

- `articlee_final.pdf` 与 `article.tex`：完整英文稿。
- `chinese/article_final.pdf` 与 `chinese/article.tex`：完整中文稿。
- `figures/` 和 `chinese/figures/`：正文使用的三张图，采用用户2026-10-28完整运行的输出。

英文在本目录执行 `pdflatex -interaction=nonstopmode -halt-on-error -jobname=articlee_final article.tex`；中文在 `chinese/` 中执行 `xelatex -interaction=nonstopmode -halt-on-error -jobname=article_final article.tex`。重复至交叉引用稳定，通常三遍。中文需要 Windows 中文字体；其他平台可将 fontset=windows 改为 fontset=fandol。

作者信息按用户要求留待后续填写。代码及复现命令见 `../code/README.md`。

当前主记录为 `20261028T112716905285Z`：Table 2、Figure 3、正文计时及对应原始记录已同步采用此次运行。早期计时CSV、环境和图3归档在 `../reproducibility/previous_records/`，不再用于当前正文。
