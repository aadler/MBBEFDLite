% Copyright (c) 2024, Avraham Adler All rights reserved
% SPDX-License-Identifier: MPL-2.0+

\name{NEWS}
\alias{NEWS}
\title{NEWS file for the MBBEFDLite package}

\section{Developmental Note: Version 0.0.0+}{
  \itemize{
    \item So long as the package remains in its experimental development
    state---noted by a 0 major version---the API may change without incrementing
    the major version. Please read these development notes carefully. Breaking
    changes will be prefaced by \dQuote{\strong{Breaking:}}.
  }
}

\section{Version 0.1.0 (2026-01-07)}{
  \subsection{Added}{
    \itemize{
      \item \strong{Breaking:} The \code{mommb} algorithm now calculates
      \eqn{E[x^2]} as \eqn{(E[x]) ^ 2 + Var[x]} instead of directly as
      \eqn{E[x ^ 2]}. Using the unbiased estimator allows for a slightly larger
      initial estimate. Given that the algorithm fails earlier, this allows for
      a slighltly larger legal starting location. The results are functionally
      identical, but not always equal within machine tolerances, so the change
      is breaking.
      \item \strong{Breaking:} The \code{mommb} algorithm now allows for passing
      a couplet consisting of the first two raw moments through \code{x} when
      new parameter \code{m} is set to \code{TRUE}. When \code{m} is
      \code{FALSE}, the default, \code{x} is treated as a vector of observations
      regardless of its length. This is breaking as the positions of the
      parameters has now changed.
      \item A new parameter \code{trace} has been added that when \code{TRUE}
      will cause \code{mommb} to print the values of \code{g} and \code{b} at
      each iteration \code{i}.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item \strong{Breaking:} The \code{mommb} algorithm now fails immediately
      when the integration would imply a negative \code{p}, instead of
      restarting with a different \code{g}. Experimentation indicates that the
      restart \emph{never} helped.
      \item Use \R \acronym{API} \code{Memzero} instead of \code{memset}.
      \item Changed Github-related files and workflows.
      \item Updated unit tests.
    }
  }
}

\section{Version 0.0.5 (2025-01-31)}{
  \subsection{Changed}{
    \itemize{
      \item Replaced \code{\()} shorthand with full \code{function} call in
      \code{mommb} to remove dependency on \R \eqn{\ge}{>=} 4.1.0.
      \item Updated package metadata test.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Use \code{log1p} instead of \code{log 1.0 + ...} where appropriate.
    }
  }
}

\section{Version 0.0.4 (2024-07-24)}{
  \subsection{Fixed}{
    \itemize{
      \item Corrected documentation in \code{ecmb}.
      \item Fixed one unit test.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Use \code{log1p} instead of \code{log 1.0 + ...} where appropriate.
    }
  }
}

\section{Version 0.0.3 (2024-06-20)}{
  \subsection{Changed}{
    \itemize{
      \item Adjust \R code for \code{mommb} to be more efficient.
      \item Update documentation.
      \item Use native \acronym{CRAN} \acronym{DOI}.
      \item Expand package metadata tests.
    }
  }
}

\section{Version 0.0.2 (2024-05-16)}{
  \subsection{Fixed}{
    \itemize{
      \item Corrected package link in \file{CITATION}.
      \item Removed extra parentheses in description of \code{g} in docs.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Expanded documentation of \code{ecmb}.
      \item Added \acronym{DOI}.
    }
  }
}

\section{Version 0.0.1 (2024-05-14)}{
    \itemize{
      \item Initial release.
  }
}
