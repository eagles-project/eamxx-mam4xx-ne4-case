#!/usr/bin/env python3

# This script adds aerosol tracers and prescribes wind velocities for a
# Cosine Bell problem that demonstrates aerosol nucleation on a global
# cubed-sphere grid within the given file.

#
# TODO: add git repo information to output file metadata
#

import argparse
import netCDF4 as nc
from math import pi, sqrt, cos, sin, atan2, nan
from shutil import copyfile
import pathlib


# ------------
# parameters
# ------------
# cross product in cartesian coordinates
def cross(x, y):
    """
    cross product in R3
    """
    return (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )


# magnitude of a vector in cartesian coordinates
def magnitude(x):
    """
    magnitude of a vector in R3
    """
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])


# dot product in cartesian coordinates
def dot(x, y):
    """
    dot product in R3
    """
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]


def cosine_bell(lat, lon, hmax=1, R=0.5):
    """
    scalar function usually used for tracer mixing ratio initial data

    lat = latitude [in]
    lon = longitude [in]
    hmax = maximum value of function [in]
    R = characteristic radius of cosine bell [in]
    returns cosine bell height
    """
    lat0, lon0 = 0, 5 * pi / 6  # lat, lon of signal center

    # find the distance r from the center of the signal to our point
    # (via cartesian coordinates)
    center = (cos(lon0) * cos(lat0), sin(lon0) * cos(lat0), sin(lat0))
    x = (cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat))
    r = atan2(magnitude(cross(x, center)), dot(x, center))

    return 0.5 * hmax * (1 + cos(pi * r / R))


def rigid_rotation_wind(lat, lon, alpha=0):
    """
    horizontal wind vectors corresponding to a rigid rotation
    about an axis inclined at an angle of alpha to the positive z-axis
    with a period of 5 units of time.

    lat [in] latitude
    lon [in] longitude
    alpha [in] angle between solid body rotation and polar axis
    retrun (u, v) zonal and meridional wind components
    """
    u0 = 2 * pi / 5  # wind velocity magnitude
    alpha = 0
    u = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
    v = -u0 * sin(lon) * sin(alpha)
    return u, v


def copy_atts(src, dst, verbose=False):
    """ """
    if verbose:
        print("copying netcdf attributes.")
    for att in src.ncattrs():
        if verbose:
            print("\tcopying {} attribute".format(att))
        dst.setncattr(att, src.getncattr(att))


def copy_dims(src, dst, verbose=False):
    """
    copies netcdf dimensions from src to dst
    """
    if verbose:
        print("copying netCDF dimensions.")
    for dim in src.dimensions:
        if verbose:
            print("\tcopying {} dimension".format(dim))
        dst.createDimension(src.dimensions[dim].name, size=src.dimensions[dim].size)


def copy_vars(src, dst, verbose=False):
    """
    copies netcdf variables from src to dst
    """
    if verbose:
        print("copying netCDF variables.")
    for varname in src.variables:
        if verbose:
            print("\tcopying {} variable".format(varname))
        src_var = src.variables[varname]
        dst_var = dst.createVariable(src_var.name, src_var.dtype, src_var.dimensions)
        atts = {}
        for a in src_var.ncattrs():
            atts[a] = src_var.getncattr(a)
        dst_var.setncatts(atts)
        dst_var[:] = src_var[:]


def get_parser():
    """
    creates a command line parser for this test case
    """
    parser = argparse.ArgumentParser(
        prog="add_aerosol_ics",
        description="adds initial conditions for aerosol tests to an existing eamxx initial condition file (does not change the input file --- output is copied to a new file).",
        epilog="Warning: if your output file exists, it will be overwritten.",
    )
    # first positional argument : input file
    parser.add_argument(
        "input_file",
        help="path to an eamxx initial condition file",
        metavar="SOURCE.nc",
    )
    # second positional argument : output file
    parser.add_argument(
        "output_file",
        help="path and name of new initial condition file with aerosol data",
        metavar="DEST.nc",
    )
    # optional argument
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=0,
        help="angle of rotation relative to polar axis",
    )
    # optional argument
    parser.add_argument(
        "--h2so4",
        type=float,
        default=2e-9,
        help="initial maximum mixing ratio (kg h2so4 / kg total air)",
    )
    # optional argument
    parser.add_argument(
        "--n_aitken",
        type=float,
        default=0,
        help="initial Aitken mode number mixing ratio (1 / kg total air)",
    )
    # optional argument
    parser.add_argument(
        "--q_aitken_so4",
        type=float,
        default=0,
        help="initial sulfate mixing ratio for Aitken mode (kg Aitken so4 / kg total air)",
    )
    # optional argument
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="write extra info to console while processing",
    )
    return parser


def check_dims(required_dims, src_ncdata, verbose=False):
    """
    checks that the required dimensions for aerosols are present in the source Dataset
    """
    if verbose:
        print("checking for required dimensions {}".format(required_dims))
    for d in required_dims:
        if d not in src_ncdata.dimensions:
            raise KeyError("dimension {} not found in source file".format(d))
        elif verbose:
            print("found {} dimension".format(d))


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    #
    # verify file extensions
    #
    if pathlib.Path(args.input_file).suffix != ".nc":
        raise ValueError("expected a .nc input file")
    if pathlib.Path(args.output_file).suffix != ".nc":
        raise ValueError("eamxx will need a .nc output file")

    #
    # copy existing data to new file
    #
    if args.verbose:
        print(
            "copying input file {} to output file {}".format(
                args.input_file, args.output_file
            )
        )
    #     copyfile(args.input_file, args.output_file)
    if args.verbose:
        print("adding additional aerosol data to output file.")
    #     ds = nc.Dataset(args.output_file, 'as')
    src_ncdata = nc.Dataset(args.input_file, "r")
    dst_ncdata = nc.Dataset(args.output_file, "w", format="NETCDF3_CLASSIC")
    copy_atts(src_ncdata, dst_ncdata)

    copy_dims(src_ncdata, dst_ncdata, args.verbose)
    copy_vars(src_ncdata, dst_ncdata, args.verbose)
    #
    # update file metadata
    #
    dst_ncdata.setncatts(
        {
            "source_file": args.input_file,
            "additional_history": "aerosol data added by eamxx-mam4xx add_aerosol_ics.py",
        }
    )

    #
    # add aerosol initial data
    #
    datatype = "float64"
    required_dims = ("time", "ncol", "lev")
    check_dims(required_dims, src_ncdata, args.verbose)

    if args.verbose:
        print("creating aerosol variables")
    q_h2so4 = dst_ncdata.createVariable(
        "q_h2so4", datatype, required_dims, fill_value=nan
    )
    q_aitken_so4 = dst_ncdata.createVariable(
        "q_aitken_so4", datatype, required_dims, fill_value=nan
    )
    n_aitken = dst_ncdata.createVariable(
        "n_aitken", datatype, required_dims, fill_value=nan
    )
    q_h2so4.setncatts(
        {
            "mdims": 1,
            "units": "kg/kg",
            "long_name": "H2SO4 gas total air mass mixing ratio",
        }
    )
    q_aitken_so4.setncatts(
        {
            "mdims": 1,
            "units": "kg/kg",
            "long_name": "SO4 aerosol total air mass mixing ratio (Aitken mode)",
        }
    )
    n_aitken.setncatts(
        {
            "mdims": 1,
            "units": "1/kg",
            "long_name": "aerosol Aitken mode total air number mixing ratio",
        }
    )

    #
    # define winds for test case from cosine_bell_winds function
    #
    horiz_winds = dst_ncdata["horiz_winds"]
    # some existing nc initial condition files have the wrong long_name for horizontal wind; we'll fix it here for the new output file.
    horiz_winds.setncattr("long_name", "horizontal wind (zonal u, meridional v)")
    lat, lon = dst_ncdata["lat"], dst_ncdata["lon"]

    if args.verbose:
        print("\tcreating aerosol data values")
    q_aitken_so4[:] = args.q_aitken_so4
    n_aitken[:] = args.n_aitken
    ncol = q_h2so4.shape[1]
    for icol in range(ncol):
        q_h2so4[0, icol, :] = cosine_bell(lat[icol], lon[icol], args.h2so4)
        u, v = rigid_rotation_wind(lat[icol], lon[icol], args.alpha)
        horiz_winds[0, icol, 0, :] = u
        horiz_winds[0, icol, 1, :] = v
    dst_ncdata.close()
    if args.verbose:
        print("output file {} is ready.".format(args.output_file))
