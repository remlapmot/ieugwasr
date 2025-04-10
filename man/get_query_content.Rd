% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{get_query_content}
\alias{get_query_content}
\title{Parse out json response from httr object}
\usage{
get_query_content(response)
}
\arguments{
\item{response}{Output from \code{\link[=api_query]{api_query()}}}
}
\value{
Parsed json output from query, often in form of data frame.
If status code is not successful then return the actual response
}
\description{
Parse out json response from httr object
}
