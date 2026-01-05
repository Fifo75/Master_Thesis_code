# Master_Thesis_code
\documentclass[11pt, a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=1in]{geometry}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{enumitem}

% Setup hyperlinks
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=cyan,
}

% Custom styling for file names to look like code
\newcommand{\code}[1]{\texttt{\textbf{\textcolor{black!80}{#1}}}}

\title{\textbf{Project Documentation: Temporal Networks and Centrality}}
\author{}
\date{}

\begin{document}

\maketitle

\section*{Overview}
This repository contains the source code for the Master's thesis titled \textbf{\textit{``Temporal network, centrality measures, and ordinary differential equations''}}.

The core implementation is written in \textbf{MATLAB}. All source files can be found in the directory \code{Thesis\_code}.

\section*{Prerequisites}
For the optimal execution of this code, the following MATLAB library is required:
\begin{itemize}
    \item \code{chebfun} (Open-source software for numerical computing with functions)
\end{itemize}

\section*{File Descriptions}
The \code{Thesis\_code} directory contains the following files and functions:

\begin{description}[style=nextline, leftmargin=1.5em, font=\ttfamily\bfseries]

    \item[temporal\_node\_evolution]
    This function calculates receiver centralities for both Katz and matrix exponential methods. It generates summary plots representing the temporal node evolution for each of the centralities.

    \item[load\_FBsoc\_dataset]
    This function loads the \code{FBsoc} dataset and constructs the adjacency matrices according to the specified time granularity (month, week, or day).

    \item[lagrange\_cardinal\_basis]
    Calculates the function handles for the Lagrange Cardinal Basis using Barycentric Interpolation for equidistant points.

    \item[interpolation\_solver]
    This function computes the matrix version of the Time-Ordered Exponential based on Section 3.5, utilizing the theoretical framework established in Chapter 3 of the thesis.

    \item[init]
    This is the base script for all experiments presented in the thesis. The individual experiments are divided into sections within this file.

    \item[genCoeffMatrix\_SP\_interval]
    Computes the coefficient matrix of $f(t)\Theta(t-s)$, where $f$ is a smooth function and $\Theta$ is the Heaviside theta function, within a basis of orthonormal Legendre polynomials.

    \item[FBsoc.txt]
    The raw text file containing all data for the FBsoc dataset.

    \item[exact\_sol]
    This function computes the exact solution of the Time-Ordered Exponential using the matrix exponential for the piecewise constant case, or \code{ode45} for the continuous case.

\end{description}

\end{document}
