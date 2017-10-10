import os
import multiprocessing as mp
import re
import itertools

#-------------------------------------------------------------------------------

def find_all_datasets(top_dir, processes=2):
    """Find all of the files in the filesystem
    
    Parameters
    ----------
    top_dir : str
        top level directory of filesystem you are using
    processes : int, optional
        number of processes to be performed in parallel
    
    Returns
    -------
    results_as_list : list
        All of the resulting files as a list object
    """
    top_levels = []

    for item in os.listdir(top_dir):
        full_path = os.path.join(top_dir, item)
        pattern = re.compile('(\d{5})|CCI')
        if pattern.match(item) is not None:
            top_levels.append(full_path)


    
    pool = mp.Pool(processes)
    results = pool.map_async(find_datasets, top_levels)
    results_as_list = list(itertools.chain.from_iterable(results.get()))

    return results_as_list

#-------------------------------------------------------------------------------

def find_datasets(data_dir):
    """Iterator to yield all datasets recursively from the base.

    Only files that include the .fits extension will be returned.

    Parameters
    ----------
    data_dir : str
        directory to search over for files

    Yields
    ----------
    root, filename : tuple
        root path and filename of each found .fits file

    """

    datasets = []

    for root, dirs, files in os.walk(data_dir):
        for filename in files:
            if not '.fits.gz' in filename:
                continue
            datasets.append((root, filename))

    return datasets

#-------------------------------------------------------------------------------
