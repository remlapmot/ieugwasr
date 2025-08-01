% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{api_query}
\alias{api_query}
\title{Wrapper for sending queries and payloads to API}
\usage{
api_query(
  path,
  query = NULL,
  opengwas_jwt = get_opengwas_jwt(),
  method = "GET",
  silent = TRUE,
  encode = "json",
  timeout = 300,
  override_429 = FALSE,
  x_api_source = paste0("ieugwasr/", utils::packageVersion("ieugwasr"))
)
}
\arguments{
\item{path}{Either a full query path (e.g. for get) or an endpoint (e.g. for post) queries}

\item{query}{If post query, provide a list of arguments as the payload. \code{NULL} by default}

\item{opengwas_jwt}{Used to authenticate protected endpoints. Login to \url{https://api.opengwas.io} to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.}

\item{method}{\code{"GET"} (default) or \code{"POST"}, \code{"DELETE"} etc}

\item{silent}{\code{TRUE}/\code{FALSE} to be passed to httr call. \code{TRUE} by default}

\item{encode}{Default = \code{"json"}, see \code{\link[httr:POST]{httr::POST}} for options}

\item{timeout}{Default = \code{300}, avoid increasing this, preferentially
simplify the query first.}

\item{override_429}{Default=\code{FALSE}. If allowance is exceeded then the query will error before submitting a request to avoid getting blocked. If you are sure you want to submit the request then set this to TRUE.}

\item{x_api_source}{Default = \code{paste0("ieugwasr/", utils::packageVersion("ieugwasr"))}.}
}
\value{
httr response object
}
\description{
There are a number of different GET and POST endpoints in the GWAS database API.
This is a generic way to access them.
}
