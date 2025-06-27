# CM SAF Metadata checker

This tool is used to check data files compliance with [CM SAF agreed metadata conventions](https://github.com/cmsaf/metadata-conventions).

## Usage

```
cmsaf-checker [-h] [-s CMSAF_METADATA_STANDARD] [-v VERSION] [-r REFERENCE] [-i IGNORE_ATTR] [-c] [-m [MISSING]] [-l] [-d DIRECTORY] files [files ...]

positional arguments:
  files

options:
  -h, --help            show this help message and exit
  -s, --cmsaf_metadata_standard CMSAF_METADATA_STANDARD
                        location of the CM SAF Metadata Standards
  -v, --version VERSION
                        CM SAF standards version to apply
  -r, --reference REFERENCE
                        Test files using this as the reference
  -i, --ignore_attr IGNORE_ATTR
                        Ignore list of attributes
  -c, --coordinates     Test coordinates time, latitude, longitude
  -m, --missing [MISSING]
                        Test for missing files
  -l, --lazy            Turn some errors to warnings
  -d, --directory DIRECTORY
                        Search for files with pattern in this directory.
```

## Examples

The CM SAF metadata-conventions project provides [sample files](https://github.com/cmsaf/metadata-conventions?tab=readme-ov-file#sample-files). These can be tested as follows
```
cmsaf-checker -c TSTin20200101000000120IMPGS01GL.nc TSTdm20200101000000120IMPGS01GL.nc
```

To check all netCDF files in the directory `foo`
```
cmsaf-checker -c -d foo "*.nc"
```
