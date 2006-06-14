\name{savejpg}
\alias{savejpg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ save device to jpg image file }
\description{
  save plot device to jpg image file
}
\usage{
savejpg(x, width = 1024, height = 768)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ file name to be saved to }
  \item{width}{ The width for the figure in pixal }
  \item{height}{ The height for the figure }

}

\value{
  For windows version, it produce bmp formatted image, otherwise,
  produce jpg images.
}

\author{ Yongming Sun }

\examples{
%- savejpg(file)
}
\keyword{ device }% at least one, from doc/KEYWORDS
