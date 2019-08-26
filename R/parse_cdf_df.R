#' Read ESM netCDF data for arbitrary variables
#'
#' The variable data should be a netCDF file of annual grid values.  The data
#' is assumed to have its dimensions of time, lat, lon, where lon is the most
#' rapidly varying index, and time is the least.
#'
#' The output will be a list with six fields:
#' \describe{
#'   \item{\strong{vardata}}{Matrix of variable data.}
#'   \item{\strong{globalop}}{Vector operator for global mean value of the variable.}
#'   \item{\strong{lat}}{Vector of latitude values.  These are only replicated
#' once, so their length is nlat, not ngrid.  You have to use \code{rep(lat,
#' nlon)} to get latitudes for each grid cell.}
#'   \item{\strong{lon}}{Vector of longitude values.  These are only replicated
#' once, so their lingth is nlon, not ngrid.  Getting longitudes for each grid
#' cell is tricksy.  Try \code{as.vector(matrix(rep(lon, nlat), nrow=nlat,
#' byrow=TRUE))}.  Fortunately, you don't need them very often.}
#'   \item{\strong{time}}{Vector of time values, given as years since the base
#' year of the dataset.}
#'   \item{\strong{tags}}{A list of datasets that were concatenated to get this
#' structure.  The names of the list are the tags given to the dataset, and the
#' values are a vector of start-row, end-row pairs.}
#' }
#'
#' The data at each time is represented as a flattened vector of grid cells.
#' The flattening is performed by transposing to lat, lon ordering, so that lat
#' will be the most rapidly varying index (because R uses
#' Fortran-style indexing).  Then the spatial dimensions are discarded,
#' resulting in a 1D vector.  The time dimension is kept, resulting in a matrix
#' with years in rows and grid cells in columns.  The dimensions of the matrix
#' will be Nyear x Ngrid, where Nyear is the number of years in the data set,
#' and Ngrid is the number of grid cells.
#'
#' The \code{globalop} vector is a vector whose dot product with a flattened grid is
#' the area-weighted global mean for the grid.  It is stored as a column vector,
#' so \code{vardata \%*\% globalop} is the time series of global mean values.
#'
#' The lat and lon dimension variables from the input file  are also
#' stored in the structure.  These are primarily useful for writing out
#' generated grids as netCDF files.  The time dimension variable is converted to
#' integer years, starting at 0.
#'
#' Conventionally, we refer to the output list as \code{griddata}.  Notably, any
#' other function with a \code{griddata} argument is expecting one of these
#' structures.
#'
#' @param filename Name of the input netCDF file
#' @param len Maximum length of the time series to read.  If the data read is
#' longer, it will be trimmed.  (Default: read entire time series, regardless of
#' length.)
#' @param varname Name of the variable to read from the netcdf file.
#' @param tag A string identifying the name of the scenario.  If omitted, the
#' tag will default to the filename.
#' @param latvar Name of the latitude dimension variable in the input file.
#' @param lonvar Name of the longitude dimension variable in the input file.
#' @param timevar Name of the time dimension variable in the input file.
#' @return A \code{griddata} list (see details).
#' @importFrom assertthat assert_that
#' @export
#'

read.general <- function(filename, len=NULL, tag=basename(filename), varname='tas',
                         latvar='lat', lonvar='lon', timevar='time')
{
    vann <- ncdf4::nc_open(filename)

    ## v3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    v3d <- ncdf4::ncvar_get(vann, var=varname)
    lat <- ncdf4::ncvar_get(vann, var=latvar)
    nlat <- length(lat)
    lon <- ncdf4::ncvar_get(vann, var=lonvar)
    nlon <- length(lon)
    time <- ncdf4::ncvar_get(vann, var=timevar)
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(vann)

    assertthat::assert_that(all(dim(v3d) == c(nlon, nlat, ntime)))

    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            time <- time[1:ntime]
            timeout <- timeout[1:ntime]
            v3d <- v3d[ , , 1:ntime]
        }
    }


    ## reorganize and flatten the 3-D array into a 2d array of ntime x ngrid
    ## As we do this, we will also make latitude the most rapidly varying index
    ## for the individual time slices.

    variable <- aperm(v3d, c(3,2,1))
    dim(variable) <- c(ntime, nlat*nlon)

    ## create a grid cell area factor that has the same shape as the flattened
    ## grid.  Since latitude is the most rapidly varying dimension in the grid,
    ## we just need to replicate cos(latitude), nlon times
    areafac <- rep(cos(lat * pi/180.0), times=nlon)
    ## normalize areafac so that sum(area*grid) is the global weighted average
    areafac <- areafac/sum(areafac)

    ## The tag element of the output is a named list of 2-element vectors,
    ## giving the start and end of the data in the output.  When datasets are
    ## concatenated, this allows us to recover the originals
    tagout <- list(c(1,ntime))
    names(tagout) <- tag

    gd <- list(vardata=as.matrix(variable), globalop=as.matrix(areafac), lat=lat, lon=lon,
               time=timeout, tags=tagout)
    class(gd) <- 'griddata'
    gd
}

read.geo.aggregate <- function(filename, len=NULL, tag=basename(filename), varname='tas',
                               timevar='time')
{
    vann <- ncdf4::nc_open(filename)

    ## v3d should have dimensions of time x lat x lon in the netcdf file.
    ## Because R uses Fortran array layout, this gets reversed to lon x lat x time.
    v3d <- ncdf4::ncvar_get(vann, var=varname)
    time <- ncdf4::ncvar_get(vann, var=timevar)
    ntime <- length(time)
    timeout <- seq_along(time) - 1
    ncdf4::nc_close(vann)

    assertthat::assert_that(all(dim(v3d) == c(ntime)))

    if(!is.null(len)) {
        ## Trim the input to the requested length.
        if(ntime < len) {
            warning('Input time series shorter than requested length.  Requested: ', len, ' Read: ', ntime)
        }
        else {
            ntime <- len
            time <- time[1:ntime]
            timeout <- timeout[1:ntime]
            v3d <- v3d[ , , 1:ntime]
        }
    }


    ## reorganize and flatten the 3-D array into a 2d array of ntime x ngrid
    ## As we do this, we will also make latitude the most rapidly varying index
    ## for the individual time slices.

    variable <- aperm(v3d, c(1))
    dim(variable) <- c(ntime)

    ## The tag element of the output is a named list of 2-element vectors,
    ## giving the start and end of the data in the output.  When datasets are
    ## concatenated, this allows us to recover the originals
    tagout <- list(c(1,ntime))
    names(tagout) <- tag

    gd <- list(vardata=as.matrix(variable),
               time=timeout, tags=tagout)
    class(gd) <- 'griddata'
    gd
}
