# CM SAF Metadata checker

This tool is used to check data files compliance with [CM SAF agreed metadata conventions](https://github.com/cmsaf/metadata-conventions).

## Usage

```
cmsaf-checker [-h] [-s PATH] [-v VERSION] [-r REFERENCE] [-i IGNORE_ATTR] [-c] [-m [MISSING]] [-l] [-d DIRECTORY] files [files ...]

positional arguments:
  files

options:
  -h, --help            show this help message and exit
  -s, --cmsaf_metadata_standard PATH
                        additional search path for CM SAF standard and GCMD
                        keyword files. Files are resolved in order: this path,
                        $CMSAF_CHECKER_PREFIX, the share folder next to the
                        installed script
  -v, --version VERSION
                        CM SAF standards version to apply
  -f, --standard_file STANDARD_FILE
                        Test files using this as the CM SAF Metadata Standard
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

## Standard and keyword file search paths

The checker resolves each reference file (CM SAF metadata standard XML, included
XML, and GCMD keyword CSV files) independently, searching in the following order:

| Priority | Source | How to set |
|----------|--------|------------|
| highest | Explicit path | `-s PATH` command-line option |
|        | Prefix directory | `CMSAF_CHECKER_PREFIX` environment variable |
| lowest | Script share folder | Installed automatically alongside the package |

Because each file is looked up individually, it is possible to override only a
single keyword CSV or the included XML with a newer version by placing it in a
higher-priority path, while the remaining files are still taken from the default
share folder.

## Examples

The CM SAF metadata-conventions project provides [sample files](https://github.com/cmsaf/metadata-conventions?tab=readme-ov-file#sample-files). These can be tested as follows
```
cmsaf-checker -c TSTin20200101000000120IMPGS01GL.nc TSTdm20200101000000120IMPGS01GL.nc
```

To check all netCDF files in the directory `foo`
```
cmsaf-checker -c -d foo "*.nc"
```
