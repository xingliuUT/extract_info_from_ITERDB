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

## Note
