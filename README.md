> Collect here scripts, snippets and ideas  
> Later also the implementations for merging into JKCS

# Readme/manual

### TODO
- [ ] `n_clusters_out` away and replaced by just specifying `n_structures_out`
- [ ] WHat to do when someone has renamed their files?
  - [ ] Specify them in the input4Clustering.csv
- [ ] Add input/output to the function docstrings

## Main: `StructureSelection.py`

### Variables
| Variable          | Datatype        | Meaning        |
| :---------------- | :-------------- | :------------- |
| `n_jobs`          | int             | Number of jobs for parallelization |
| `n_clusters_init` | int             | k of k-means |
| `n_clusters_out`  | int             | n of clusters selected from initial k of k-means |
| `n_structures_out`| int             | n of structures outputted from the selection |
| `level`           | str             | parameter defining whether XTB of DFT is used |
| `normEd`          | boolean         | Option to normalise energies to zero and convert to kcal/mol |
| `own_cmap`        | [20*str]        | Define colornames to be used in plotting |
| `wrkdir`          | str             | Path to current $WRKDIR |
| `data_df`         | DataFrame       | Structure info from resultsXXX.dat |
| `structure_index` | Index           | Indeces from `data_df.index`|
| `xyz_df`          | DataFrame       | Structure xyz from movieXXX.dat |


| `datafile`        | str             | Path of resultsXXX.dat |
| `xyzfile`         | str             | Path of movieXXX.xyz |




| `path_output`     | str             | Path to output folder inside $WRKDIR |
| `dirname`         | str             | Input the name of scanned dir to `scan_xtb` |
| `filetypes`       | [str]           | Input desired filetypes (other than .log or .xyz) to `scan_xtb` |
| `subfolders`      | [str]           | Output subfolders from `scan_xtb`, used by `scan_xtb` for recursion |
| `files`           | [str]           | Output list of files (other than .log or .xyz) from `scan_xtb` |
| `zyxs`            | [str]           | Output list of .xyz names from `scan_xtb` |
| `zyx_paths`       | [str]           | Output list of .xyz filepaths from `scan_xtb` |
| `logs`            | [str]           | Output list of .log names from `scan_xtb` |
| `log_paths`       | [str]           | Output list of .log filepaths from `scan_xtb` |
| `path_XTB`        | str             | Path to XTB-folder |
| `zyxs_df`         | DataFrame       | Sorted DataFrame of .xyz paths |
| `logs_df`         | DataFrame       | Sorted DataFrame of .log paths |
| `structures`      | DataFrame       | ASE-data of all structures indexed |
| `zero_structure`  | ASE molecule    | For aquiring the elements and the size of the cluster |
| `zero_name`       | str             | The name (index) of the first structure |
| `chemsyms`        |                 | List of all atomsnames in the cluster |
| `chemsyms_uniques`| [str]           | List of which elements the clusters have |
| `n_atoms`         | int             | Number of atoms in the clusters |
| `logtext`         | [str]           | input of `findEdip`, a logfile to be read |
| `logdata`         | DataFrame       | Energy and Dipole (and GyRadius) values from XTB-logs |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |
| ``        |  |  |

### Functions
| Function       | in             | out            | Purpose        |
| :------------- | :------------- | :------------- | :------------- |
| `makedir` | `dirname` |  | Make folders if they don't exist |
| `scan_xtb` | `dirname`, `filetypes` | `subfolders`, `files`, `zyxs`, `zyx_paths`, `logs,log_paths` | Scan XTB folder recursively and return paths and names of .log and .xyz -files |
| `view` | ASE structure | visualization | Visualize the molecule/cluster |
| `findEdip` | XTB logfile as array of lines read | Values of E and dip | Find Energy and dipoles from XTB logfile |
| `calcgr` | `xyzfile` | GyRadius | Calculate Radius of Gyration from given xyzfile |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |
| `` | `` | `` |  |




### Argparse
> Not needed when input4Clustering.csv is used
```python3
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    # Add own arguments here
    parser.add_argument("-n", "--n_jobs",type=int,required=True,help="Provide the number of parallel jobs")
    parser.add_argument("-k", "--n_clusters_init",type=int,required=True,help="Provide the number of initial clusters for k-means")
    parser.add_argument("-c", "--n_clusters_out",type=int,required=True,help="Provide the number of clusters selected from k-means")
    parser.add_argument("-s", "--n_structures_out",type=int,required=True,help="Provide the number of structures = local minima outputted")
    parser.add_argument("-l", "--level", type=str,required=True,help="[XTB]/[DFT] Whether to perform selection on XTB or DFT structures.")
    parser.add_argument("--normEd", action='store_true',help="Normalise energies and convert to kcal/mol")
    parser.add_argument("--verbose", action='store_true',help="Whether to print out some progress.")
    args = parser.parse_args()
else:
    print("I don't know what you are doing but some parameters may not work..")
```
