"""
Converts NaviEdit exported .xyz and .nav files into std_data
"""
import logging
import os
import sys

import numpy as np
import pandas as pd

from auvlib.data_tools import all_data, std_data


def process_xyz_pings_file(filepath):
    """Read a .xyz file containing info about each beam hits and return a processed dataframe
    where each row contains info for 1 MBES swath"""
    logging.info(f'\tProcessing xyz file: {filepath}...')
    naviedit_pings = pd.read_csv(filepath, delimiter='\t', header=0)

    # Aggregate xyz hits of each beam into a hits array
    naviedit_pings['hit'] = naviedit_pings[['X', 'Y',
                                            'Z']].apply(lambda row: list(row),
                                                        axis=1)

    # Parse string datetime values and produce a datetime object
    naviedit_pings['time_string'] = naviedit_pings[[
        'yyyy', 'mmddhhmm', 'ss.ss'
    ]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    naviedit_pings['time_string'] = pd.to_datetime(
        naviedit_pings['time_string'], format='%Y.%m%d%H%M.0%S.%f')
    naviedit_pings['time_stamp'] = naviedit_pings['time_string'].map(
        pd.Timestamp.timestamp)

    # Group the beams by Ping No -> aggregate hits from the same ping into a list
    beams = naviedit_pings.groupby('Ping No').agg({'hit': list})

    # Get unique pings
    ping_info = naviedit_pings[[
        'time_string', 'time_stamp', 'Ping No', 'Tide', 'Heading', 'Heave',
        'Pitch', 'Roll'
    ]].drop_duplicates().set_index('Ping No')

    # Final results with 1 row = 1 MBES swath
    mbes_ping_result = ping_info.join(beams)
    return mbes_ping_result


def process_nav_file(filepath):
    """Read a .nav file and return a processed dataframe where each row = 1 positional info with timestamps."""
    logging.info(f'\tProcessing nav file: {filepath}...')
    naviedit_nav = pd.read_csv(filepath, header=0)

    # Parse string datetime values and produce a datetime object
    naviedit_nav['time_string'] = naviedit_nav[['DATE', 'TIME']].apply(
        lambda row: ''.join(row.values.astype(str)), axis=1)
    naviedit_nav['time_string'] = pd.to_datetime(naviedit_nav['time_string'],
                                                 format='%Y%m%d%H:%M:%S.%f')
    naviedit_nav['time_stamp'] = naviedit_nav['time_string'].map(
        pd.Timestamp.timestamp)
    naviedit_nav.rename(columns={
        'EASTING': 'easting',
        'NORTHING': 'northing',
        'DEPTH': 'depth'
    },
                        inplace=True)

    return naviedit_nav


def merge_pings_and_nav_df(pings_df, nav_df):
    """Merge the processed pings and navigation dataframes based on their timestamps. Return a merged dataframe
    where each row contains info about 1 MBES swath and the corresponding vehicle position in ENU coordinates."""
    logging.info('\tMerging pings and nav df...')
    merged = pings_df.copy()
    merged[['easting', 'northing']] = 0
    for row_idx, row in merged.iterrows():
        ping_time = row['time_stamp']
        nav_row_idx = nav_df['time_stamp'].le(ping_time).idxmin()
        nav_row = nav_df.iloc[nav_row_idx]
        nav_row_prev = nav_df.iloc[nav_row_idx - 1]

        ratio = (ping_time - nav_row_prev['time_stamp']) / (
            nav_row['time_stamp'] - nav_row_prev['time_stamp'])
        east = nav_row_prev['easting'] + ratio * (nav_row['easting'] -
                                                  nav_row_prev['easting'])
        north = nav_row_prev['northing'] + ratio * (nav_row['northing'] -
                                                    nav_row_prev['northing'])
        depth = nav_row_prev['depth'] + ratio * (nav_row['depth'] -
                                                 nav_row_prev['depth'])

        merged.loc[row_idx, ['easting', 'northing', 'depth']] = nav_row[[
            'easting', 'northing', 'depth'
        ]]
    return merged


def construct_std_data_from_merged_df(merged_df):
    """Translate the merged_df dataframe into a list of std_data objects"""
    logging.info('\tConstructing std_data from merged df...')
    merged_data = []
    for row_idx, row in merged_df.iterrows():
        row_data = std_data.mbes_ping()

        if row_idx == merged_df.index.min():
            row_data.first_in_file_ = True

        row_data.id_ = row_idx
        row_data.beams = [np.array(x) for x in row['hit']]
        row_data.heading_ = row['Heading']
        row_data.heave_ = row['Heave']
        row_data.pitch_ = row['Pitch']
        row_data.pos_ = np.array(
            [row['easting'], row['northing'], row['depth']])
        row_data.roll_ = row['Roll']
        row_data.time_stamp_ = int(row['time_stamp'])
        row_data.time_string_ = row['time_string'].strftime(
            '%Y-%m-%d %H:%M:%S.%f')

        merged_data.append(row_data)
    return merged_data


def construct_std_data_from_naviedit_export(folder, store=True):
    """Given a folder with .xyz and .nav exports from NaviEdit, extract relevant fields
    and construct std_data from these files, optionally store the std_data in .cereal
    format to disk."""
    filenames = sorted(
        set(
            x.split('.')[0] for x in os.listdir(folder)
            if x.split('.')[-1] in ['xyz', 'nav']))
    for filename in filenames:
        if not os.path.exists(os.path.join(
                folder, f'{filename}.xyz')) or not os.path.exists(
                    os.path.join(folder, f'{filename}.nav')):
            raise IOError(
                f'{filename}.xyz or {filename}.nav missing. Aborting...')

    std_data_dict = {}
    for i, filename in enumerate(filenames):
        logging.info(f'Processing {i}/{len(filenames)} file {filename}...')
        pings_df = process_xyz_pings_file(
            os.path.join(folder, f'{filename}.xyz'))
        nav_df = process_nav_file(os.path.join(folder, f'{filename}.nav'))
        merged_df = merge_pings_and_nav_df(pings_df, nav_df)
        merged_std_data = construct_std_data_from_merged_df(merged_df)

        std_data_dict[filename] = merged_std_data
        outpath = os.path.join(folder, f'{filename}.cereal')
        if store:
            std_data.write_data(merged_std_data, outpath)
    return std_data_dict


def main():
    folder = sys.argv[1]

    # Set up logger
    logger = logging.getLogger(name='data processor logger')
    fhandler = logging.FileHandler(filename=os.path.join(
        folder, 'data_processing.log'),
                                   mode='a')
    fhandler.setLevel(logging.INFO)
    logger.addHandler(fhandler)

    # Construct std_data from the .xyz and .nav files from the given folder
    data_dict = construct_std_data_from_naviedit_export(folder=folder)

    melt_data_combined = []
    for section in data_dict.values():
        melt_data_combined.extend(section)
    std_data.write_data(melt_data_combined,
                        os.path.join(folder, 'merged.cereal'))


if __name__ == '__main__':
    main()
