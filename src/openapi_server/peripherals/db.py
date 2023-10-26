from datetime import datetime

from typing import Optional, Tuple, Union
from astropy.units import Quantity

from astropy.coordinates import SkyCoord
import astropy.units as u

import pandas as pd
import numpy as np
import sqlite3

## code adapted from https://github.com/lgrcia/twirl and https://github.com/ppp-one/gaia-tmass-sqlite

def db_query(db: str, min_dec: float, max_dec: float, min_ra: float, max_ra: float) -> pd.DataFrame:
    """
    Queries a federated database for astronomical data within a specified range of declination and right ascension.

    Args:
        db (str): The path to the SQLite database file.
        min_dec (float): The minimum declination value to query.
        max_dec (float): The maximum declination value to query.
        min_ra (float): The minimum right ascension value to query.
        max_ra (float): The maximum right ascension value to query.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the queried astronomical data.
    """

    conn = sqlite3.connect(db)

    # Determine the relevant shard(s) based on the query parameters.
    arr = np.arange(np.floor(min_dec), np.ceil(max_dec) + 1, 1)
    relevant_shard_ids = set()
    for i in range(len(arr) - 1):
        shard_id = f"{arr[i]:.0f}_{arr[i+1]:.0f}"
        relevant_shard_ids.add(shard_id)

    # Execute the federated query across the relevant shard(s).
    df_total = pd.DataFrame()
    for shard_id in relevant_shard_ids:
        shard_table_name = f'{shard_id}'
        if min_ra > max_ra:
            q = f'SELECT * FROM `{shard_table_name}` WHERE dec BETWEEN {min_dec} AND {max_dec} AND (ra BETWEEN {min_ra} AND 360 OR ra BETWEEN 0 AND {max_ra})'
        else:
            q = f'SELECT * FROM `{shard_table_name}` WHERE dec BETWEEN {min_dec} AND {max_dec} AND ra BETWEEN {min_ra} AND {max_ra}'
        df = pd.read_sql_query(q, conn)
        df_total = pd.concat([df, df_total], axis=0)

    # Close the conn and return the results.
    conn.close()
    return df_total

def gaia_db_query(
    center: Union[Tuple[float, float], SkyCoord],
    fov: Union[float, Quantity],
    limit: int = 10000,
    tmass: bool = False,
    dateobs: Optional[datetime] = None,
    db_path: str = '',
) -> np.ndarray:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query. By default, it is set to True.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from twirl import gaia_radecs
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> fov = 0.1
    >>> radecs = gaia_radecs(center, fov)
    """

    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    min_dec = dec - dec_fov/2
    max_dec = dec + dec_fov/2
    min_ra = ra - ra_fov/2
    max_ra = ra + ra_fov/2

    if min_ra < 0:
        min_ra += 360

    if max_ra > 360:
        max_ra -= 360

    if min_dec < -90:
        min_dec = -90

    if max_dec > 90:
        max_dec = 90

    table = db_query(db_path, min_dec, max_dec, min_ra, max_ra)
    if tmass:
        table = table.sort_values(by=['j_m']).reset_index(drop=True)
    else:
        table = table.sort_values(by=['phot_g_mean_mag']).reset_index(drop=True)

    table.replace('', np.nan, inplace=True)
    table.dropna(inplace=True)
    
    # limit number of stars
    table = table[0:limit]

    # add proper motion to ra and dec
    if dateobs is not None:
        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25 # type: ignore
        
        years = dateobs - 2015.5 # type: ignore
        table["ra"] += years * table["pmra"] / 1000 / 3600
        table["dec"] += years * table["pmdec"] / 1000 / 3600    
        
    return np.array([table["ra"].values, table["dec"].values]).T, table['j_m'].values

def gaia_radecs(
    center: Union[Tuple[float, float], SkyCoord],
    fov: Union[float, Quantity],
    limit: int = 10000,
    circular: bool = True,
    tmass: bool = False,
    dateobs: Optional[datetime] = None,
) -> np.ndarray:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query. By default, it is set to True.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from twirl import gaia_radecs
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> fov = 0.1
    >>> radecs = gaia_radecs(center, fov)
    """
    from astroquery.gaia import Gaia

    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    radius = np.max([ra_fov, dec_fov]) / 2

    if circular and not tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY gaia.phot_rp_mean_flux DESC
            """
        )
    elif circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux, tmass.j_m
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY tmass.j_m
            """
        )
    elif not circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux, tmass.j_m
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY tmass.j_m
            """
        )
    else:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux
            FROM gaiadr2.gaia_source AS gaia
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY gaia.phot_rp_mean_flux DESC
            """
        )

    table = job.get_results()
    
    # add proper motion to ra and dec
    if dateobs is not None:
        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25 # type: ignore
        
        years = dateobs - 2015.5 # type: ignore
        table["ra"] += years * table["pmra"] / 1000 / 3600
        table["dec"] += years * table["pmdec"] / 1000 / 3600

    
    if tmass:
        table.remove_rows(np.isnan(table['j_m']))
        return np.array([table["ra"].value.data, table["dec"].value.data]).T, table['j_m'].value.data
    else:
        table.remove_rows(np.isnan(table['phot_rp_mean_flux']))
        return np.array([table["ra"].value.data, table["dec"].value.data]).T, table['phot_rp_mean_flux'].value.data