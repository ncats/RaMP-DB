% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RaMP.R
\name{show,RaMP-method}
\alias{show,RaMP-method}
\alias{RaMP}
\alias{show}
\alias{listRaMPVersions}
\title{Connection to a RaMP database}
\usage{
\S4method{show}{RaMP}(object)

RaMP(version = character())

listRaMPVersions(local = FALSE)
}
\arguments{
\item{version}{`character(1)` specifying the RaMP version to load. By
default (`version = character()`), the most recent release will be
used.}

\item{local}{`logical(1)` for `listRaMPVersion`: whether remote
(`local = FALSE`, default) or locally (`local = TRUE`) available RaMP
versions should be listed.}
}
\description{
Connections to a *RaMP* database can be established and managed with the
`RaMP` function. The returned `RaMP` object provides the reference to the
database and it can be passed to the various functions to query that
specific database. RaMP databases are provided as self-contained SQLite
databases that are automatically downloaded and locally cached with the
`RaMP` function. The caching mechanism prevents repeated downloads of the
same database version.

- `RaMP`: eventually download and connect to a RaMP database. Parameter
  `version` allows to specify the RaMP release version to which a connection
  should be established. If the specified version is not available locally,
  it will be downloaded and cached. Use `listRaMPVersions()` to list
  available local or remote databases. Alternatively, the connection to a
  RaMP database can be directly provided through parameter `dbcon`.

- `listRaMPVersions`: list available local or remote RaMP database releases.
}
