A study of the codec2 700C Rate K resampler

# Basic setup:

$ sudo apt install texmaker texlive-bibtex-extra
$ texmaker ratek_resampler.tex &

Options -> Config Texmaker -> Quick Build -> select 2nd option: "PDFLatex+BibLatex+PDFLatex2+ViewPDF"

# EPS Latex figures:

Watch out for any control characters like & or _ in figure text, this will cause the build step to choke.

Comment out any font size defaults in `.octaverc`:
```
#set(0, "defaulttextfontsize", 24)  % title
#set(0, "defaultaxesfontsize", 24)  % axes labels
#set(0, "defaultlinelinewidth", 2)
```
Create your figure then:

octave:1> print("testepslatex","-depslatex","-S300,300");

In `/path/to/ratek_resampler` Optionally clean out any older .eps files of same name (e.g. previously rendered with encapsulated Postscript).

$ cp testepslatex.tex testepslatex-inc.eps /path/to/ratek_resampler

Latex code to insert figure:

\begin{figure}[h]
\caption{Test Octave -depslatex}
\begin{center}
\input testepslatex.tex
\end{center}
\end{figure}

In texmaker, build PDF from ratek_resampler.tex, not testepslatex.tex, you might end up in the latter if an error was encountered.
