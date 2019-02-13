from glob import glob
from tqdm import tqdm
import glob
import shutil
import os
from os import path
import datetime as dt
import sys



output_base_directory = "D:\\singleRAW"

def age_of_file(filename):
    """Returns age of file creation, in hours."""
    st = os.stat(filename)
    mtime = dt.datetime.fromtimestamp(st.st_mtime)
    age = dt.datetime.now() - mtime
    return age.total_seconds() / 60. / 60.

def get_valid_maxquant_filenames(directory=os.path.join(r'Z:\maintenance', '2018'), minsize=10000000, creation_age=120,
                                 age_min=0.1):
    """
    Returns a list of filenames in "directory" that are younger than "creation_age" hours ago and bigger than "minsize" bytes.
    """
    Filesall = []
    for path in tqdm(glob(os.path.join(directory, '**/*.raw'), recursive=True)):
        if age_of_file(path) < creation_age and os.path.getsize(path) > minsize and age_of_file(path) > age_min:
            Filesall.append(path)

    return Filesall


# Get list of valid files for analysis
Filestodo = []
for filename in get_valid_maxquant_filenames('Z:/maintenance', minsize=10000000, creation_age=120, age_min=0.1):

    # completed files
    output_directory = path.join(output_base_directory, path.splitext(path.basename(filename))[0])
    if path.exists(path.join(output_directory, "combined", "txt", "summary.txt")):
        continue

    # crashed past maxquant analyses
    # todo: handle repeated-crashing files
    elif path.exists(output_directory):
        if age_of_file(output_directory) < 8.:
            continue
        else:
            print ("!!! Detected Incomplete Analysis:",  (filename), "!!!")
            continue
    else:
        Filestodo.append((output_directory, filename))

print('Files to go: {}'.format(len(Filestodo)))

if len(Filestodo) == 0:
    sys.exit()

# Select the file to analyze
analysis_directory, raw_filename = Filestodo[0]
print('Starting {}'.format(raw_filename))

# Copy raw into the new directory.
os.mkdir(analysis_directory)
shutil.copy(raw_filename, analysis_directory)
