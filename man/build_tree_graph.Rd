% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{build_tree_graph}
\alias{build_tree_graph}
\title{Build tree graph}
\usage{
build_tree_graph(
  clusterings,
  prefix,
  count_filter,
  prop_filter,
  metadata,
  node_aes_list
)
}
\arguments{
\item{clusterings}{numeric matrix containing clustering information, each
column contains clustering at a separate resolution}

\item{prefix}{string indicating columns containing clustering information}

\item{count_filter}{count threshold for filtering edges in the clustering
graph}

\item{prop_filter}{in proportion threshold for filtering edges in the
clustering graph}

\item{metadata}{data.frame containing metadata on each sample that can be
used as node aesthetics}

\item{node_aes_list}{nested list containing node aesthetics}
}
\value{
\link[tidygraph:tbl_graph]{tidygraph::tbl_graph} object containing the tree graph
}
\description{
Build a tree graph from a set of clusterings, metadata and associated
aesthetics
}
\keyword{internal}
