\name{NEWS}
\alias{NEWS}
\title{News for Package \pkg{prefmod}}
\encoding{UTF-8}
\section{Changes in Version 0.9-0}{
  \itemize{
    \item \code{bin2int()} now throws and error if \code{bin} contains more than 31 elements (causing and overflow).
    \item Lots of improvements suggested by \CRANpkg{lintr}.
    \item \CRANpkg{prefmod}'s DLL now unloads upon unloading the package.
    \item Implemented \CRANpkg{testthat} tests.
  }
}
\section{Changes in Version 0.8-38}{
  \itemize{
    \item Minor fixes.
  }
}
\section{Changes in Version 0.8-37}{
  \itemize{
    \item Minor fixes.
  }
}
\section{Changes in Version 0.8-36}{
  \itemize{
    \item Minor fixes.
  }
}
\section{Changes in Version 0.8-35}{
  \itemize{
    \item Minor fixes and updated documentation.
  }
}
\section{Changes in Version 0.8-34}{
  \itemize{
    \item Corrected the coding of the rating-scale items of the ISSP 2000.
      \code{CAR}, \code{IND}, \code{FARM}, \code{WATER}, \code{TEMP}, and \code{GENE} were coded 0 (dangerous) to 4 (not dangerous) instead of 1 (dangerous) to 5 (not dangerous).
    \item \code{plotworth()} is defunct -- use \code{plot()} (or \code{plot.wmat}) instead.
    \item Minor fixes and improvements.
  }
}
