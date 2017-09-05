# extract_info_from_ITERDB

## Description
This project contains code to extract information from [ITERDB](https://inis.iaea.org/search/search.aspx?orig_q=RN:46090459) files.

`read_iterdb_x.py' constructs a dictionary based on the data from ITERDB file.
* te (electron temperature), ne (electron density), ti, ni, nz, vrot
* for each quantity, there's a radial coordinate that goes with it

`write_iterdb.py` contains methods that could take the profiles and write into the ITERDB file format.

`diffFromMTM.py` computes diffusivity estimates based on magnetic fluctuation: `D = pi * v_th_e * Ls * B_tilda^2`, where `Ls = q * R / shat` 

`D_chi_from_source.py` computes particle diffusivity (`D`) and heat diffusivity (`chi`) using total stored particles and energy divided by heating rate and particle loss.

`doubleCheckITERDB.py` compares and plots profiles from two ITERDB files.

`rhoiL.py` computes the ratio between ion gyroradius `rhoi` and pedestal scale length `L`.

## Usage

Usage of `read_iterdb_x.py` and `write_iterdb` is in `TestDrive_iterdb.py`. To run `TestDrive_iterdb.py`, make sure that there's an ITERDB file and an EFIT file in the working directory:
```
python TestDrive_efit.py CMOD1120815027May12.iterdb g1120907032.01012
```
`diffFromMTM.py` and `D_chi_from_source.py` follow the same format, taking two arguments: ITERDB file, EFIT file.

To compare two ITERDB files, use `doubleCheckITERDB.py`:
```
python doubleCheckITERDB.py CMOD1120815027May12.iterdb CMOD1120815027May08.iterdb
```
## Note
