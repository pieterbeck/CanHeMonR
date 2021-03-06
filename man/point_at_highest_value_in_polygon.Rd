% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PWN_functions.R
\name{point_at_highest_value_in_polygon}
\alias{point_at_highest_value_in_polygon}
\title{Place Point At Maximum Location}
\usage{
point_at_highest_value_in_polygon(r, pol)
}
\arguments{
\item{r}{A raster layer}

\item{pol}{A SpatialPolygon}
}
\value{
A SpatialPoint
}
\description{
Convert a SpatialPolygon to a SpatialPoint, assigning the point to the location inside the polygong that
has the highest value in a matching raster
}

