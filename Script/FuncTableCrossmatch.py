"""
File: FuncTableCrossmatch.py
Name: Chia-Lin Ko
Create Date: Jun 09, 2021
Last Modified Date: Jun 11, 2021
------------------------
This program aims to crossmatch tables.

Modified from Zhen-Kai Gao's script
Modified from astroML.crossmatch.crossmatch
"""
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

def crossmatch(X1, X2, max_distance=np.inf):
    """
    Purpose
        Cross-match the values between dataset X1 and X2.
        By default, this uses a KD Tree for speed.
        This function is modified from astroML.crossmatch.crossmatch
    ---------------------
    Input Parameter
        X1             [array]: first dataset, shape(N1, D1).
        X2             [array]: second dataset, shape(N2, D2), D1 should be the same as D2.
        max_distance   [float]: optional, maximum radius of search.
                                If no point is within the given radius, then inf will be returned.
    ---------------------
    Return
        ind         [ndarrays]: The index of the closest point in X2 to each point in X1
                                Array is length N1.
                                Locations with no match are indicated by ind[i] = []
    """
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    N1, D1 = X1.shape
    N2, D2 = X2.shape
    if D1 != D2:
        raise ValueError('Arrays must have the same second dimension')

    kdt = cKDTree(X2)
    ind = kdt.query_ball_point(X1, r=max_distance)
    
    return ind


def crossmatch_angular(X1, X2, max_distance=np.inf, return_dist=False):
    """
    Purpose
        Cross-match angular values between dataset X1 and X2.
        By default, this uses a KD Tree for speed.
        Because the KD Tree only handles cartesian distances, the angles
        are projected onto a 3D sphere.
        This function is modified from astroML.crossmatch.crossmatch_angular
    ---------------------
    Input Parameter
        X1             [array]: first dataset, shape(N1, 2).
                                X1[:, 0] is the RA in degrees.
                                X1[:, 1] is the DEC in degrees.
        X2             [array]: second dataset, shape(N2, 2).
                                X2[:, 0] is the RA in degrees.
                                X2[:, 1] is the DEC in degrees.
        max_distance   [float]: optional, maximum radius of search, measured in degrees.
                                If no point is within the given radius, then inf will be returned.
        return_dist     [bool]: optional, If True, return the distance array,
                                otherwise only return the index array (default)/
    ---------------------
    Return
        ind, (dist) [ndarrays]: the index and angular distance of each pair.
                                Locations with no match are indicated by dist[i] = inf, ind[i] = []
    """
    # Convert the RA, Dec, and max_distance from deg to rad
    X1 = np.deg2rad(X1)
    X2 = np.deg2rad(X2)
    max_distance = np.deg2rad(max_distance)

    # Convert 2D RA/DEC to 3D cartesian coordinates
    Y1 = np.transpose(np.vstack([np.cos(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 1])]))
    Y2 = np.transpose(np.vstack([np.cos(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 1])]))

    # Law of cosines to compute 3D distance
    max_y = np.sqrt(2 - 2 * np.cos(max_distance))
    ind = crossmatch(Y1, Y2, max_y)
    
    # Calculate distances
    if return_dist:
        def get_dist(y1, y2):
            if len(y2) == 0:
                dist = np.array([np.inf])
            else:
                y1 = np.broadcast_to(y1, y2.shape)
                dist = np.diag(cdist(y1, y2))
            return dist
        dist = [get_dist(v, Y2[ind[i]]) for i, v in enumerate(Y1)]
        dist = np.array(dist, dtype='object')

        # Convert distances back to angles using the law of tangents
        not_inf = np.vectorize(lambda x: ~np.isinf(x).any())(dist)
        def convert_dist(d):
            x = 0.5 * d
            return (180. / np.pi * 2 * np.arctan2(x, np.sqrt(np.maximum(0, 1 - x ** 2))))
        dist[not_inf] = np.array([convert_dist(i) for i in dist[not_inf]], dtype='object')        
        return ind, dist    
    else:
        return ind


def ravel_list(lst):
    import functools
    import operator
    return functools.reduce(operator.iconcat, lst, [])


def set_low_values_to_zero(x, tol=1e-9):
    x[x < tol] = 0.
    return x


def merge_df(df1, df2, ind, dist=None):
    """
    Purpose
        Merge the dataframe df1 and df2.
        Join type is All from df1. 
        Match selection is All matches.
    ---------------------
    Input Parameter
        df1    [DataFrame]: first dataset.
        df2    [DataFrame]: second dataset.
        ind     [ndarrays]: the index of each pair.
        dist    [ndarrays]: (optional), the angular distance (in degree) of each pair.
                            If given, then the new column 'Separation' (in arcsec) will be 
                            added to the returned merged dataset.
    ---------------------
    Return
        mfd    [DataFrame]: the merge dataset. 
                            If the row number of merged dataset is different from the row number of df1, 
                            then thes new column 'Count' will be added to the returned merged dataset.
    """
    
    # Left ind for merge
    rep = [max(len(i), 1) for i in ind]
    left_ind = np.repeat(np.arange(len(df1)), rep)    
    # Right ind for merge
    df2 = df2.copy()
    df2.loc[df2.index.max() + 1] = None
    right_ind = ravel_list([i if len(i) > 0 else [df2.index.max()] for i in ind])    
    # Merge
    mdf = pd.merge(left=df1.iloc[left_ind].reset_index(drop=True), 
                   right=df2.iloc[right_ind].reset_index(drop=True),
                   left_index=True, right_index=True, suffixes=['_1', '_2'])   
    if 'Unnamed: 0' in mdf.columns:
        mdf.drop(mdf.columns[mdf.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    # Column: Separation 
    if dist is not None:    
        dist = dist*3600  # deg to arcsec
        sep = ravel_list([set_low_values_to_zero(i, tol=1e-9) for i in dist])
        # if 'Separation' in mdf.columns:
        #     sep_id      = mdf.columns.str.contains('Separation', case = False).sum()
        #     column_sep  = 'Separation' + '_%s'%(sep_id)
        # else:
        column_sep  = 'Separation'
        mdf[column_sep] = sep
        mdf[column_sep].replace(np.inf, np.nan, inplace=True)
        mdf.sort_values(by=[mdf.columns[0], column_sep])
        
    # Column: Count 
    if df1.shape[0] != mdf.shape[0]:
        rep0 = [max(len(i), 0) for i in ind]
        count = np.repeat(rep0, rep)
        # if 'Count' in mdf.columns:
        #     count_id      = mdf.columns.str.contains('Count', case = False).sum()
        #     column_count  = 'Count' + '_%s'%(count_id)
        # else:
        column_count  = 'Count'
        mdf[column_count] = count

    return mdf


def best_match(mdf, column):
    mdf = mdf.drop_duplicates(subset=[column], keep='first').reset_index(drop=True)
    return mdf


def inner_join(mdf, column='Separation'):
    mdf = mdf.dropna(subset=[column]).reset_index(drop=True)
    return mdf