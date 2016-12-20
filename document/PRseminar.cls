%========================== PRseminar.cls ============================
%
% for PR seminar in Komurasaki Lab.
% Version 1.0
% Jun. 11 2013 ver 1.0
% Rei Kawashima
% Modified jjap3.cls by The Japan Society of Applied Physics
%
%=================================================================

\typeout{***** Document Class =<<<jjap3.cls>>>= ver. 2.0b ******}

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{jjap3}[2011/05/11 ver 2.0b JJAP document class]

\RequirePackage{amsmath}[1997/03/20]
\RequirePackage{amssymb}
\RequirePackage[dvipdfmx]{graphicx}
\RequirePackage[dvipdfmx]{color}
\RequirePackage{url}
%\RequirePackage{cite}
\RequirePackage{overcite}
\renewcommand\citeright{)}
\def\@citess#1{\textsuperscript{#1\citeright}}

\def\jjap3@warning{\ClassWarningNoLine{jjap3}}
\def\jjapjversion{2.0b}

\newif\if@ap
\DeclareOption{apex}{\@aptrue}
\newif\if@ip
\DeclareOption{ip}{\@iptrue}
\newif\if@hr
\DeclareOption{hr}{\@hrtrue}
\newif\if@st
\DeclareOption{stap}{\@sttrue}
\newif\if@rc
\DeclareOption{rc}{\@rctrue}
\newif\if@bn
\DeclareOption{bn}{\@bntrue}
\newif\if@cr
\DeclareOption{cr}{\@crtrue}
\newif\if@er
\DeclareOption{er}{\@ertrue}

\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}

\newif\if@lp
\DeclareOption{letterpaper}
   {\@lptrue%
    \setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}

\newif\if@seceq
\DeclareOption{seceq}{%
\@seceqtrue
}

\DeclareOption{twocolumn}{\@twocolumntrue
% To balance the final page,
% following code is borrowed from nicetwocolumn style
\newdimen\dimen@@

\def\balancecolumns{\setbox0\vbox{\unvbox\@outputbox} \dimen@@=\ht0
   \advance\dimen@@ by\topskip \advance\dimen@@ by-\baselineskip
   \divide\dimen@@ by2 \splittopskip=\topskip
   {\vbadness=10000 \loop \global\setbox3=\copy0
     \global\setbox1=\vsplit3 to\dimen@@
     \ifdim\ht3>\dimen@@ \global\advance\dimen@@ by1pt\repeat}
   \ifdim\dimen@@>\textheight \dimen@@=\textheight\fi
   \ifdim\dimen@@<\textheight \@lastpagebalancingfalse
   \setbox0=\vbox to \dimen@@{\unvbox1} \setbox2=\vbox to \dimen@@{\unvbox3}
   \if\ht0>\ht2 \setbox2=\vbox to\ht0{\unvbox2}\dimen@@=\ht0\else
   \if\ht2>\ht0 \setbox0=\vbox to\ht2{\unvbox0}\dimen@@=\ht2\fi\fi
  \wd0=\hsize \wd2=\hsize
%\message{in balancecolumns dimen@@=\the\dimen@@, dp0=\the\dp0, dp2=\the\dp2}
\ifdim\dp0>\dp2\advance\dimen@@ by\dp0\else\advance\dimen@@ by\dp2\fi
\setbox\@outputbox\vbox to\dimen@@{\hbox to\textwidth{\box0\hfil\box2}\vfill}}

% you can refrain from balancing by falsifying following switch!
\newif\if@lastpagebalancing \@lastpagebalancingfalse % default is to balance

% incorporate balanced output to \enddocument
\def\enddocument{\if@lastpagebalancing \global\let\@outputdblcol=\balancing@outputdblcol\fi%
\@checkend{document}\clearpage
\begingroup  
\if@filesw \immediate\closeout\@mainaux 
\def\global\@namedef##1##2{}\def\newlabel{\@testdef r}%
\def\bibcite{\@testdef b}\@tempswafalse \makeatletter\input \jobname.aux
\if@tempswa \@warning{Label(s) may have changed.  Rerun to get
cross-references right}\fi\fi\endgroup
\deadcycles\z@\@@end}

\newlength{\@interlinegap}

% to change \interlinegap appropriately at appropriate time
\def\@outputdblcol{\if@firstcolumn \global\@firstcolumnfalse
    \global\setbox\@leftcolumn\box\@outputbox
% Following interlinegap value may be wrong!!
% It is set to gap+depth for normalsize, and depth(2.031999Q) for footnotesize.
% This setting is inconsistent, but seems to work properly now.
%\ifx \@currsize\normalsize \global\@interlinegap 8.54815pt\else% setup
% reset to 10Q after final page footnote hack!
\ifx \@currsize\normalsize \global\@interlinegap \z@\else %12\p@% setup
 \ifx \@currsize\footnotesize \global\@interlinegap \z@\else% appropriate
  \ifx \@currsize\small \global\@interlinegap \z@\fi\fi\fi% gap here!
  \else \global\@firstcolumntrue
\if@leftfootnoteexist\setbox\@leftcolumn\vbox to\@colht{\boxmaxdepth\maxdepth
\@texttop\dimen128=\dp\@leftcolumn
\unvbox\@leftcolumn\vskip\skip\footins\footnoterule%
%\message{in outputdblcol leftcolumnfootnote=ht:\the\ht\@leftcolumnfootnote, dp:\the\dp\@leftcolumnfootnote, wd:\the\wd\@leftcolumnfootnote}
\unvbox\@leftcolumnfootnote
\vskip-\dimen128\@textbottom}
\global\maxdepth\@maxdepth\else
\setbox\@leftcolumn\vbox to\@colht{\boxmaxdepth\maxdepth
\@texttop\dimen128=\dp\@leftcolumn
\unvbox\@leftcolumn\vfill%
\vskip-\dimen128\@textbottom}
\global\maxdepth\@maxdepth\fi
\if@rightfootnoteexist\setbox\@outputbox\vbox to\@colht{\boxmaxdepth\maxdepth
\@texttop\dimen128=\dp\@outputbox
\unvbox\@outputbox\vskip\skip\footins\footnoterule%
%\message{in outputdblcol rightcolumnfootnote=ht:\the\ht\@rightcolumnfootnote, dp:\the\dp\@rightcolumnfootnote, wd:\the\wd\@rightcolumnfootnote}
\unvbox\@rightcolumnfootnote
\vskip-\dimen128\@textbottom}
\global\maxdepth\@maxdepth\else
\setbox\@outputbox\vbox to\@colht{\boxmaxdepth\maxdepth
\@texttop\dimen128=\dp\@outputbox
\unvbox\@outputbox\vfill%
\vskip-\dimen128\@textbottom}
\global\maxdepth\@maxdepth\fi
   \setbox\@outputbox\vbox{\hbox to\textwidth{\hbox to\columnwidth
     {\box\@leftcolumn \hss}\hfil \vrule width\columnseprule\hfil
      \hbox to\columnwidth{\box\@outputbox \hss}}}\@combinedblfloats
      \@outputpage \begingroup \@dblfloatplacement \@startdblcolumn
      \@whilesw\if@fcolmade \fi{\@outputpage\@startdblcolumn}\endgroup
   \fi}

\def\balancing@outputdblcol{\if@firstcolumn %\global\@firstcolumnfalse
 \else \global\@firstcolumntrue
  \setbox\@outputbox\vbox{\unvbox\@leftcolumn\vskip\@interlinegap%
   \penalty\interlinepenalty\unvbox\@outputbox}
 \fi
 \balancecolumns \@combinedblfloats
\if@leftfootnoteexist\setbox\@outputbox\vbox to\@colht{\boxmaxdepth\maxdepth
\@texttop\dimen128=\dp\@outputbox
\unvbox\@outputbox\vfill\vskip\skip\footins\footnoterule%
%\message{in balancing@outputdblcol leftcolumnfootnote=ht:\the\ht\@leftcolumnfootnote, dp:\the\dp\@leftcolumnfootnote, wd:\the\wd\@leftcolumnfootnote}
\unvbox\@leftcolumnfootnote
\if@rightfootnoteexist
%\message{in balancing@outputdblcol rightcolumnfootnote=ht:\the\ht\@rightcolumnfootnote, dp:\the\dp\@rightcolumnfootnote, wd:\the\wd\@rightcolumnfootnote}
\unvbox\@rightcolumnfootnote\fi
\vskip-\dimen128
\@textbottom}
\global\maxdepth\@maxdepth\fi
 \@outputpage \begingroup \@dblfloatplacement \@startdblcolumn
 \@whilesw\if@fcolmade \fi{\@outputpage\@startdblcolumn}\endgroup}

\newbox\@leftcolumnfootnote \newbox\@rightcolumnfootnote
\newif\if@leftfootnoteexist \newif\if@rightfootnoteexist

% Extra \@texttop somehow found its way into \@makecol.  Deleted
% 1 Dec 86.  (Found by Mike Harrison)

\def\@makecol{\if@firstcolumn
\ifvoid\footins \setbox\@outputbox\box\@cclv
%\message{in makecol @outputbox=(\the\ht\@outputbox, \the\dp\@outputbox, \the\wd\@outputbox)}
\global\@leftfootnoteexistfalse
   \else\global\@leftfootnoteexisttrue
\setbox\@outputbox\box\@cclv
%\message{in makecol with foot @outputbox=(\the\ht\@outputbox, \the\dp\@outputbox, \the\wd\@outputbox)}
\global\setbox\@leftcolumnfootnote\box\footins
%\message{in makecol leftcolumnfootnote=ht:\the\ht\@leftcolumnfootnote, dp:\the\dp\@leftcolumnfootnote, wd:\the\wd\@leftcolumnfootnote}
\fi
     \xdef\@freelist{\@freelist\@midlist}\gdef\@midlist{}\@combinefloats
\else
\ifvoid\footins \setbox\@outputbox\box\@cclv
%\message{in makecol right @outputbox=(\the\ht\@outputbox, \the\dp\@outputbox, \the\wd\@outputbox)}
\global\@rightfootnoteexistfalse
   \else\global\@rightfootnoteexisttrue
\setbox\@outputbox\box\@cclv
%\message{in makecol right with foot @outputbox=(\the\ht\@outputbox, \the\dp\@outputbox, \the\wd\@outputbox)}
\global\setbox\@rightcolumnfootnote\box\footins
%\message{in makecol rightcolumnfootnote=ht:\the\ht\@rightcolumnfootnote, dp:\the\dp\@rightcolumnfootnote, wd:\the\wd\@rightcolumnfootnote}
\fi
     \xdef\@freelist{\@freelist\@midlist}\gdef\@midlist{}\@combinefloats
\fi
  }
}

\ExecuteOptions{a4paper}
\ProcessOptions\relax

\allowdisplaybreaks[3]

\def\@typeset{Japanese Journal of Applied Physics}
\def\@apexmark{Appl.{\space}Phys.{\space}Express}
\def\@APEXmark{Applied Physics Express}
\def\@jjapmark{PR Seminar, FY2015}
\def\@JJAPmark{PR Seminar, FY2015} 
\def\@rpmark{Rei Kawashima}
\def\@ipmark{INVITED REVIEW PAPER}
\def\@hrmark{COMPREHENSIVE REVIEW}
\def\@stmark{SELECTED TOPICS IN APPLIED PHYSICS}
\def\@rcmark{RAPID COMMUNICATION}
\def\@bnmark{BRIEF NOTE}
\def\@ermark{ERRATA}


\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}

\let\makefigurecaptions\relax

\frenchspacing
\@twosidetrue
\AtBeginDocument\sloppy
\AtEndDocument{\label{lastpage}}

\newcommand{\degC}{\mbox{$^{\circ}$C}}
\let\Cdot\textbullet
\newcommand{\mbi}[1]{\boldsymbol{#1}}
\newcommand{\mathendash}{\mbox{--}}
\newcommand{\address}{\relax}
\newcommand{\name}{\relax}
\newcommand{\surname}{\relax}
\newcommand{\note}{\relax}
\newcommand{\tbreak}{\relax}
\newcommand{\conline}{(Color{\space}online)}
\newcommand{\MEXT}{Ministry of Education, Culture, Sports, Science and Technology}
\newcommand{\JSPS}{Japan Society for the Promotion of Science}

\DeclareMathOperator{\textIm}{Im}
\DeclareMathOperator{\textRe}{Re}
\DeclareMathOperator{\tr}{tr}
\DeclareMathOperator{\Tr}{Tr}

\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{1}
\setlength\parskip{\z@}

\if@twocolumn
	\setlength{\textwidth}{504\p@}
	\setlength{\textheight}{720\p@}
	\setlength{\columnsep}{14\p@}
	\hoffset-23.5\p@
	%\setlength\oddsidemargin{-1cm}
	%\setlength\evensidemargin{-1cm}
	\setlength\topmargin{-1.2cm}
	\setlength\headheight{8\p@}
	\setlength\headsep{10\p@} 
	\setlength\topskip{24\p@}
	\setlength\parindent{10\p@}
	\setlength\footskip{20\p@}
\else
	\setlength{\textwidth}{444\p@}
	\setlength{\textheight}{660\p@}
	\hoffset5\p@
	\if@lp
		\voffset-1.2cm
	\else
		\voffset-.5cm
	\fi
	\setlength\headheight{10\p@}
	\setlength\headsep{12\p@} 
	\setlength\topskip{20\p@}
	\setlength\parindent{18\p@}
	\setlength\footskip{36\p@}
\fi
\addtolength\textheight{\topskip}

\setlength{\skip\footins}{9\p@ \@plus 4\p@ \@minus 2\p@}
\setlength\floatsep    {24\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\textfloatsep{24\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\intextsep   {24\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\dblfloatsep    {24\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\dbltextfloatsep{24\p@ \@plus 4\p@ \@minus 4\p@}
\setlength\@fptop{0\p@ \@plus 1fil}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\setlength\partopsep{\z@}

\setlength\arraycolsep {6\p@} 
\setlength\tabcolsep {6\p@}
\setlength\arrayrulewidth {.25\p@}
\setlength\doublerulesep {2\p@}
\setlength\tabbingsep {\labelsep}
\setlength\fboxsep{3\p@} 
\setlength\fboxrule {.2\p@}
\setlength\marginparwidth {.75\p@} 
\setlength\marginparsep {12\p@}

\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301

\setcounter{topnumber}{5}
\renewcommand\topfraction{.9}
\setcounter{bottomnumber}{5}
\renewcommand\bottomfraction{.9}
\setcounter{totalnumber}{10}
\renewcommand\textfraction{.1}
\renewcommand\floatpagefraction{.9}
\setcounter{dbltopnumber}{3}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.1}

\newcommand\abstsize{%
   \if@twocolumn
     \@setfontsize\abstsize\@ixpt{11}
   \else
     \@setfontsize\abstsize\@xpt{16}
   \fi
}

\renewcommand\normalsize{%
   \if@twocolumn
   \@setfontsize\normalsize\@xpt{12}
   \else
   \@setfontsize\normalsize\@xiipt{20}
   \fi
   \abovedisplayskip 6\p@
   \abovedisplayshortskip 6\p@
   \belowdisplayshortskip 6\p@
   \belowdisplayskip \abovedisplayskip
   \setlength\jot{6\p@}
   \let\@listi\@listI}
\normalsize
\newcommand\small{%
   \if@twocolumn
   \@setfontsize\small\@ixpt{11}
   \else
   \@setfontsize\small\@xpt{16}%
   \fi
   \abovedisplayskip 6\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \def\@listi{\leftmargin\leftmargini
               \topsep \z@ %%%%% 4\p@ \@plus2\p@ \@minus2\p@
               \parsep \z@ %%%%% 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\footnotesize{%
   \if@twocolumn
   \@setfontsize\footnotesize\@viiipt{10}
   \else
   \@setfontsize\footnotesize\@ixpt{14}%
   \fi
   \abovedisplayskip 5\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \def\@listi{\leftmargin\leftmargini
               \topsep \z@ %%%%% 3\p@ \@plus\p@ \@minus\p@
               \parsep \z@ %%%%% 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\scriptsize{%
    \if@twocolumn
    \@setfontsize\scriptsize\@viipt\@ixpt
      \else
    \@setfontsize\scriptsize\@viiipt\@ixpt
    \fi
}
\newcommand\tiny{%
    \if@twocolumn
     \@setfontsize\tiny\@vipt\@viiipt
       \else
     \@setfontsize\tiny\@viipt\@viiipt
     \fi
}
\newcommand\large{%
    \if@twocolumn
      \@setfontsize\large\@xiipt{16}
        \else
      \@setfontsize\large\@xivpt{20}
      \fi
}
\newcommand\Large{\@setfontsize\Large\@xviipt{22}}
\newcommand\LARGE{\@setfontsize\LARGE\@xxpt{25}}
\newcommand\huge{\@setfontsize\huge\@xxvpt{30}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\newcommand\titlesize{\@setfontsize\titlesize\@xiipt{16}}

\@mparswitchfalse
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

\renewcommand\thesection {\@arabic\c@section}
\renewcommand\thesubsection {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection {\thesubsection.\@arabic\c@subsubsection}
\renewcommand\theparagraph {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}

\if@seceq
\@addtoreset{equation}{section}
\renewcommand{\theequation}{%
 \ifnum\c@section>\z@\thesection.\fi \@arabic\c@equation}
\fi

\newcommand\@secapp{}

\newcommand\section{\@startsection{section}{1}{\z@}%
                   {\if@twocolumn -8\p@ \@plus\p@ \@minus\p@ \else -.5\baselineskip \fi}%
                   {\if@twocolumn 4\p@ \@plus\p@ \@minus\p@ \else \p@ \fi}%
                   {\normalfont\normalsize\sffamily\bfseries\boldmath\raggedright}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                   {-12\p@ \@plus \p@ \@minus \p@}{\p@ \@plus\p@ \@minus\p@}%
                   {\normalfont\normalsize\sffamily\raggedright}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                   {-12\p@ \@plus \p@ \@minus \p@}{\p@ \@plus\p@ \@minus\p@}%
                   {\normalfont\normalsize\sffamily}} 
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                   {\z@}{-10\p@}%
                   {\normalfont\normalsize\sffamily\bfseries\boldmath}}
\newcommand\subparagraph{\@startsection{subparagraph}{4}{\parindent}%
                   {\z@}{-10\p@}%
                   {\normalfont\normalsize\sffamily}}

\def\@nullstring{}
\def\@appendixstring{Appendix}

\def\@sect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \def\@svsec{}\else
       \ifnum #2=1		% TRUE if \section{}
         \refstepcounter{#1}
         \def\@appendixtitle{#8}
         \ifx\@secapp\@appendixstring
         \ifx\@appendixtitle\@nullstring
             \edef\@svsec{\@secapp \csname the#1\endcsname\hskip10\p@}
         \else
             \edef\@svsec{\@secapp \csname the#1\endcsname:\hskip10\p@}
         \fi
         \else
           \ifx\thesection\@nullstring
             \edef\@svsec{\@secapp \csname the#1\endcsname\hskip10\p@}
           \else
             \edef\@svsec{\@secapp \csname the#1\endcsname.\hskip10\p@}
           \fi
         \fi
       \else
         \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname\hskip10\p@}
       \fi
     \fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
        \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
        \endgroup
       \csname #1mark\endcsname{#7}
      \else
        \def\@svsechd{#6\hskip #3\@svsec #8\csname #1mark\endcsname{#7}}
      \fi
     \@xsect{#5}}

\setcounter{secnumdepth}{3}

\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{(\theenumi)}
\newcommand\labelenumii{\theenumii.}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}

\setlength\leftmargini {20\p@}
\setlength\leftmarginii {10\p@} 
\setlength\leftmarginiii {10\p@} 
\setlength\leftmarginiv {10\p@} 
\setlength\leftmarginv {10\p@}
\setlength\leftmarginvi {10\p@}
\setlength\leftmargin{\leftmargini}
\setlength\labelwidth{\leftmargini\advance\labelwidth-\labelsep}
\setlength\labelsep {5\p@} 
\setlength\parsep {\z@}
 
\def\@listi{\leftmargin\leftmargini
            \parsep \z@ %%%4\p@ \@plus2\p@ \@minus2\p@
            \topsep \z@ %%%4\p@ \@plus2\p@ \@minus2\p@
            \itemsep \z@ %%%\@plus2\p@ \@minus2\p@
            }
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    \z@ \@plus2\p@ \@minus\p@
              \parsep    \z@ \@plus2\p@ \@minus\p@
              \itemsep   \parsep}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    \z@ \@plus2\p@ \@minus\p@
              \parsep    \z@
              \partopsep \z@ \@plus2\p@ \@minus\p@
              \itemsep   \topsep}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
 
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -15\p@%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 15\p@}%
                \item\relax}
               {\endlist}

\newenvironment{quotation}
               {\list{}{\listparindent 15\p@%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}

\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}

\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

%%%% bibliography %%%%
\def\thebibliography#1{%
\raggedright
\if@twocolumn
\vskip1.25\baselineskip
%\hbox to\columnwidth{\hrulefill}
\footnotesize
\noindent\hrulefill\vskip4\p@
\else
\clearpage
%\section*{\refname}
\noindent\hskip-.5em\normalsize{\sffamily\bfseries\refname}\par\vskip4\p@
\fi
\list{\arabic{enumi})}{\settowidth\labelwidth{[#1]}
 \leftmargin\labelwidth
 \advance\leftmargin\labelsep
 \usecounter{enumi}}
% \sloppy
 \sfcode`\.=500\relax}
\let\endthebibliography=\endlist

\renewcommand\footnoterule{%
	\noindent\vskip1\baselineskip \kern-3\p@
	\hrule\@width\columnwidth
	\kern 2.6\p@
}

\newcommand\refname{References}
\newcommand\figurename{Fig.}
\newcommand\tablename{Table}
\newcommand\appendixname{Appendix}

\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{1em}
\setlength\belowcaptionskip{\z@} %{1em}

\newcount\halfflag
\halfflag=0

%%%%%%%%%%captions%%%%%%%%%%%%%%%%%%%%%%
\long\def\@makecaption#1#2{\raggedright%
\if@twocolumn
\footnotesize
	\ifnum \halfflag>0
	\par
	\vskip\abovecaptionskip
	\hsize=86mm
	\setbox\@tempboxa\hbox{#1\hskip1em#2}
		\ifdim \wd\@tempboxa<\hsize
		\global \@minipagefalse
		\hb@xt@\hsize{\hfil\box\@tempboxa\hfil}
		%\noindent\hbox to \hsize{\hfill\box\@tempboxa\hfill}\par
		\else
		#1\hskip1em#2\par
		%\par\unhbox\raggedright\@tempboxa\par
		\fi
	\else
	\par
	\vskipp\abovecaptionskip
	\hsize=177mm
	\setbox\@tempboxa\hbox{#1\hskip1em#2}
		\ifdim \wd\@tempboxa<\hsize
		\global \@minipagefalse
		\hb@xt@\hsize{\hfil\box\@tempboxa\hfil}
		%\noindent\hbox to \hsize{\hfill\box\@tempboxa\hfill}\par
		\else
        #1\hskip1em#2\par
		%% \centerline{%
		%% \parbox{155mm}
		%% \par
		%%\unhbox\raggedright\@tempboxa %}
		%%\par
		\fi
    \fi
\else
\small
\vskip\abovecaptionskip
\sbox\@tempboxa{#1\hskip1em#2}%
%%%	\ifdim \wd\@tempboxa > \hsize
%%%	\par
	#1\hskip1em#2\par
%%%	\else
%%%	\global \@minipagefalse
%%%	\hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%%%	\fi
\vskip\belowcaptionskip
\fi
}

\long\def\@caption#1[#2]#3{\par
     \begingroup
     \@parboxrestore
     \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}
     \par
     \endgroup
}

\newcounter{figure}
\def\thecaption{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\fnum@figure{%
	\if@twocolumn
	\footnotesize
	\else
	\small
	\fi
	\textsf{\bfseries\figurename\hskip.5em\thefigure.
}}
\newcounter{table}
\def\thetable{\@Roman\c@table}
\def\fps@table{tbp}
\def\ftype@table{2} 
\def\fnum@table{%
	\if@twocolumn
	\footnotesize
	\else
	\small
	\fi
	\textsf{\bfseries\tablename\hskip.5em\thetable.}
}

\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height 1\p@ \futurelet
   \reserved@a\@xhline}
\def\@xhline{\ifx\reserved@a\hline
               \vskip\doublerulesep
               \vskip-\arrayrulewidth
             \fi
      \ifnum0=`{\fi}}

\newenvironment{figure}
               {\halfflag=1\@float{figure}}
               {\end@float}
               
\newenvironment{figure*}
               {\halfflag=0\@dblfloat{figure}}
               {\end@dblfloat}

\newenvironment{table}
               {\halfflag=1\@float{table}}
               {\end@float}
               
\newenvironment{table*}
               {\halfflag=0\@dblfloat{table}}
               {\end@dblfloat}

\newenvironment{halffigure}
               {\halfflag=1\@float{figure}}
               {\end@float}
               
\newenvironment{fullfigure}
               {\halfflag=0\@dblfloat{figure}}
               {\end@dblfloat}

\newenvironment{halftable}
               {\halfflag=1\@float{table}}
               {\end@float}
               
\newenvironment{fulltable}
               {\halfflag=0\@dblfloat{table}}
               {\end@dblfloat}

\def\tabular{%
	\if@twocolumn
	\begin{center}
		\else
		\relax
	\fi
	\vskip.5\baselineskip
	\let\@halignto\@empty%
	\if@twocolumn\footnotesize\else\small\relax\fi\@tabular
}
\def\endtabular{\crcr\egroup\egroup $\egroup
	\if@twocolumn
	\end{center}
		\else
		\relax
	\fi
}

\let\halftabular\tabular\relax
\let\fulltabular\tabular\relax
\let\endhalftabular\endtabular\relax
\let\endfulltabular\endtabular\relax

%%%%%%%%%%%%%
% maketitle %
%%%%%%%%%%%%%

\def\maketitle{\par
\begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\hskip1\p@%\rlap{
	\@textsuperscript{\normalfont\@thefnmark}} %}%
    \long\def\@makefntext##1{%\vskip2\p@ 
	        \hangindent8\p@ \hangafter1 \noindent
            \hb@xt@1em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
 \if@twocolumn
   \twocolumn[\@maketitle]
 \else
   \global\@topnum\z@ 
 \@maketitle 
 \fi
\thispagestyle{firstpage}\@thanks
\endgroup
 \setcounter{footnote}{0}
 \global\let\thanks\relax
 \global\let\maketitle\relax
 \global\let\@maketitle\relax
 \global\let\@thanks\@empty
 \global\let\@author\@empty
 \global\let\@title\@empty
 \global\let\title\relax
 \global\let\author\relax
}

\def\@maketitle{%
\null
\relax
{\raggedright\large\sffamily\bfseries\boldmath\raggedright%
	\if@twocolumn\vskip24\p@
	\else
	\relax
	\fi
	\@title \par
}

\if@twocolumn
\vskip8\p@
  \else
\vskip12\p@
  \fi
{\raggedright\normalsize\sffamily\raggedright\@author \par}

\if@twocolumn
\vskip4\p@
  \else
\vskip12\p@
  \fi
{\raggedright\small\sffamily\itshape\raggedright\@inst \par}

\if@twocolumn
\vskip10\p@
\@ifundefined{@recdate}{\def\@recdate{}}{}
\@ifundefined{@revdate}{\def\@revdate{}}{}
\@ifundefined{@accdate}{\def\@accdate{}}{}
\else
\vskip6\p@
  \@ifundefined{@recdate}{}{%
  \noindent\footnotesize (Received\space\@recdate)
  \par}
\fi

\gdef\@recdate{}
\gdef\@revdate{}
\gdef\@accdate{}

\if@twocolumn
  \vskip4\p@
  \@ifundefined{@abst}{}{%
  \fboxsep6\p@%
  \abstsize\noindent\framebox{\parbox{.975\textwidth}{\sffamily\@abst}}
  \par
}
\else
  \vskip6\p@
  \noindent\hbox to\textwidth{\hrulefill}\par
  \@ifundefined{@abst}{}{%
  \abstsize\sffamily\noindent\@abst \par
}
\fi

%\if@twocolumn
% \vskip6\p@
% \@ifundefined{@kword}{}{%
%  \noindent\small\sffamily KEYWORDS:\hskip.5em\parbox[t]{144mm}{\raggedright\@kword}
%  \par
%}
%\else
% \vskip12\p@
 \@ifundefined{@kword}{}{%
   \relax %%%\noindent\small\sffamily KEYWORDS:\hskip.5em\parbox[t]{124mm}{\raggedright\@kword}
  \par
}
%\fi

\if@twocolumn
  \vskip2\baselineskip
\else
  \noindent\vskip-4pt\hbox to\textwidth{\hrulefill}\par
  \vskip2\baselineskip
\fi
}

\long\def\abst#1{\long\gdef\@abst{#1}}
\long\def\inst#1{\long\gdef\@inst{#1}}
\long\def\kword#1{\long\gdef\@kword{#1}}
%
\def\recdate#1{\gdef\@recdate{#1}}
\def\revdate#1{\gdef\@revdate{#1}}
\def\accdate#1{\gdef\@accdate{#1}}

%%%%%%Definition of headings%%%%%%
\mark{{}{}}
\def\ps@firstpage{%
  \let\@mkboth\@gobbletwo%
  \def\@evenfoot{\small\sffamily%
  \if@twocolumn
  \hfill\thepage\hfill
  \else
  \hfill\thepage/\pageref{lastpage}\hfill
  \fi
}
  \let\@oddfoot\@evenfoot
  \def\@evenhead{\small\sffamily%
  \underline{\hbox to\textwidth{\if@ap\@APEXmark\else\@JJAPmark\fi%
  	\hfill%
  	  \if@ap\relax
  	    \else
      \if@ip\bfseries\@ipmark
        \else
      \if@hr\bfseries\@hrmark
        \else
      \if@st\bfseries\@stmark
        \else
      \if@rc\bfseries\@rcmark
        \else
      \if@bn\bfseries\@bnmark
        \else
      \if@cr\relax
        \else
      \if@er\relax
        \else
      \@rpmark
       \fi
         \fi
           \fi
             \fi
               \fi
                 \fi
                   \fi
                     \fi
     		 }
     	 }
      }
	\let\@oddhead\@evenhead
}

\def\ps@jjap{%
  \let\@mkboth\@gobbletwo%
  \def\@evenfoot{\small\sffamily%
  \if@twocolumn
  \hfill\thepage\hfill
  \else
  \hfill\thepage/\pageref{lastpage}\hfill
  \fi
  }
  \let\@oddfoot\@evenfoot
  \def\@evenhead{\small\sffamily%
  \underline{\hbox to\textwidth{\if@ap\@apexmark\else\@jjapmark\fi\hfill%
      \if@ap\relax
        \else
      \if@ip\bfseries\@ipmark
        \else
      \if@hr\bfseries\@hrmark
        \else
      \if@st\bfseries\@stmark
        \else
      \if@rc\bfseries\@rcmark
        \else
      \if@bn\bfseries\@bnmark
        \else
      \if@cr\relax
        \else
      \if@er\relax
        \else
      \@rpmark
       \fi
         \fi
           \fi
             \fi
               \fi
                 \fi
                   \fi
                     \fi
     		 }
     	 }
      }
	\let\@oddhead\@evenhead
}

%\renewcommand{\month}[1]{\ifcase#1\relax ???\or
% January\or February\or March\or April\or May\or June\or
% July\or August\or September\or October\or November\or December
% \else ???\fi}

\def\today{\ifcase\month\or
 January\or February\or March\or April\or May\or June\or
 July\or August\or September\or October\or November\or December\fi
 \space\number\day, \number\year}

\pagestyle{jjap}
\pagenumbering{arabic}

%
% Acknowledgments
%

\newcommand{\acknowledgment}{\par%
  \if@rc\par\vskip2\p@\footnotesize\raggedright\textsf{\bfseries Acknowledgment}\hskip1em
    \else
  \if@bn\par\vskip2\p@\footnotesize\raggedright\textsf{\bfseries Acknowledgment}\hskip1em
    \else
      \section*{Acknowledgment}
    \fi
  \fi
}
\newcommand{\acknowledgments}{\par%
  \if@rc\raggedright\paragraph{\small Acknowledgments}\small
    \else
  \if@bn\raggedright\paragraph{\small Acknowledgments}\small
    \else
      \section*{Acknowledgments}
    \fi
  \fi
}
\newcommand{\acknowledgement}{\par%
  \if@rc\raggedright\paragraph{\small Acknowledgement}\small
    \else
  \if@bn\raggedright\paragraph{\small Acknowledgement}\small
    \else
      \section*{Acknowledgement}
    \fi
  \fi
}
\newcommand{\acknowledgements}{\par%
  \if@rc\raggedright\paragraph{\small Acknowledgements}\small
    \else
  \if@bn\raggedright\paragraph{\small Acknowledgements}\small
    \else
    \section*{Acknowledgements}
    \fi
  \fi
}


%%%%%%%%% Appendix %%%%%%%%%%%
\newcount\@appendixflag
\@appendixflag=0
\newcount\@appendixindex
\@appendixindex=0

\def\getappendix#1{\global\@namedef{appnedix@#1}{1}}
\def\appendix{\if@twocolumn\par \else\clearpage \fi
\renewcommand{\section}[1]{\sectiona{##1}}
\normalsize}

\def\sectiona{%
 \c@table=0
 \c@figure=0
 \c@equation=0
 \setcounter{equation}{0}
 \ifnum \@appendixflag=0
   \setcounter{section}{0} \@appendixflag=1
   \def\@secapp{Appendix}
   \def\thetable{\Alph{section}\textperiodcentered\@arabic\c@table}%
   \def\thefigure{\Alph{section}\textperiodcentered\@arabic\c@figure}
   \def\theequation{\Alph{section}\textperiodcentered\arabic{equation}}
 \fi

 \@ifundefined {appnedix@1}{\def\thesection{}}
   {\def\thesection{\ \Alph{section}}}

 \def\thesubsection {\Alph{section}.\arabic{subsection}}
 \if@filesw \immediate\write\@auxout
       {\string\getappendix{\the\@appendixindex}}
 \advance \@appendixindex by 1 \fi
\@startsection{section}{1}{\z@}%
  {\if@twocolumn -8\p@ \@plus\p@ \@minus\p@ \else -.5\baselineskip \fi}%
  {\if@twocolumn 4\p@ \@plus\p@ \@minus\p@ \else \p@ \fi}%
  {\normalfont\normalsize\sffamily\bfseries\boldmath}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\endinput
%% End of file `jjap3.cls'.
