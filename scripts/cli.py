#!/usr/bin/env python3
from netCDF4 import Dataset, num2date
from astropy.time import Time
from dateutil.relativedelta import relativedelta
import calendar as cal
import csv
import datetime
import os
import re
from types import SimpleNamespace
from typing import NamedTuple
import xml.etree.ElementTree as ET

import numpy as np
import pytz

__version__ = "3.3.0"
__prefix__  = ""
STANDARD = ''

CMSAF_NAMING_STANDARD = r'^([A-Z]{3})([cdhimpswa])([ncfhmds])((?:19|20)\d\d)(0[1-9]|1[012])(0[1-9]|[12]\d|3[01])([01]\d|2[0-3])([0-5]\d)(\d{3})(\d{2})([0-9A-Z]{5})([A-Z0-9]{2})([A-Z0-9]{2})((?:\.hdf|\.hdf\.gz|\.gz|\.nc|\.nc\.gz){0,1})$'

RC_ERR  = "## ERROR ##"
RC_WARN = "## WARNING ##"
RC_OK   = "## OK ##"
RC_FAIL = "## FAILED ##"
RC_INFO = "## INFORMATION ##"



def normalize_whitespace(text):
    "Remove redundant whitespace from a string."
    return ' '.join(text.split())


def significant_digits(value):
    indx = np.nonzero(value)
    res  = np.multiply(value, 0)
    res[indx[0]] = np.log10(np.absolute(value[indx[0]]))
    res = np.power(np.float64(10), np.floor(res)+1)

    return res


def float_spacing (x, y):
    ax  = np.abs(np.spacing(x))
    ay  = np.abs(np.spacing(y))
    res = np.where(ax > ay, ax, ay)

    return res



# Grid code → resolution in degrees, as defined in the CM SAF naming convention (PDF table, page 8).
# Lat/lon grids only; projection-based codes are intentionally absent because they have no
# single scalar resolution and are handled as non-regular grids (return None).
_GRID_RESOLUTION: dict[str, np.float64] = {
    "17": np.float64("0.03"),
    "23": np.float64("0.05"),
    "26": np.float64("0.1"),
    "19": np.float64("0.25"),
    "13": np.float64("0.5"),
    "22": np.float64("0.625"),
    "20": np.float64("1.0"),
}

# Known native-projection grid codes (no scalar resolution).
# Used to suppress spurious warnings for recognised non-regular grids.
# "05" = native SEVIRI grid.
_PROJECTION_GRIDS: frozenset[str] = frozenset({"05"})

# Known polar grids.
# "15" = Lambert azimuthal equal-area projection (15x15 km²).
# "21" = Lambert azimuthal equal-area projection (25x25 km²).
_POLAR_GRIDS: frozenset[str] = frozenset({"15", "21"})

# Mapping from dtype token (used in XML 'dtype' attributes) to the set of numpy
# scalar types that token accepts.  Type-class shorthands ('int', 'float', 'str')
# group entire families; specific tokens map to a singleton.
# Used by _matches_dtype() for both variable-level and attribute-level checks.
_DTYPE_CLASSES: dict[str, frozenset] = {
    'int':     frozenset({np.int8,   np.int16,  np.int32,  np.int64,
                          np.uint8,  np.uint16, np.uint32, np.uint64}),
    'float':   frozenset({np.float32, np.float64}),
    'str':     frozenset({str, np.bytes_, np.str_}),
    'int8':    frozenset({np.int8}),
    'int16':   frozenset({np.int16}),
    'int32':   frozenset({np.int32}),
    'int64':   frozenset({np.int64}),
    'uint8':   frozenset({np.uint8}),
    'uint16':  frozenset({np.uint16}),
    'uint32':  frozenset({np.uint32}),
    'uint64':  frozenset({np.uint64}),
    'float32': frozenset({np.float32}),
    'float64': frozenset({np.float64}),
}


def _matches_dtype(value, dtype_spec: str) -> bool:
    """Return True if *value* matches any token in the comma-separated *dtype_spec*.

    *value* may be:
      - np.dtype  — the variable's own .dtype (e.g. np.dtype('int16'))
      - str       — a plain Python string attribute value
      - numpy scalar or ndarray — a numeric attribute value

    An empty *dtype_spec* always returns True (no constraint declared).
    Unknown tokens are tried as literal numpy dtype names; invalid names are
    silently ignored (no match).
    """
    if not dtype_spec:
        return True
    if isinstance(value, np.dtype):
        actual = value.type
    elif isinstance(value, str):
        actual = str
    else:
        actual = np.asarray(value).dtype.type
    for token in (t.strip() for t in dtype_spec.split(',')):
        allowed = _DTYPE_CLASSES.get(token)
        if allowed is not None:
            if actual in allowed:
                return True
        else:
            try:
                if actual == np.dtype(token).type:
                    return True
            except TypeError:
                pass
    return False


def cmsaf_decode_grid(filename):
    """
    Decode the grid resolution (degrees) from a CM SAF standard filename.

    Returns the resolution as np.float64 for known lat/lon grid codes,
    or None if the filename does not match the naming convention or the
    grid code has no scalar resolution (e.g. satellite projection grids).
    """
    decode = re.match(CMSAF_NAMING_STANDARD, filename)
    if decode is None:
        return None
    code = decode.group(10)
    resolution = _GRID_RESOLUTION.get(code)
    if resolution is None and code not in _PROJECTION_GRIDS and code not in _POLAR_GRIDS:
        print(f"{RC_WARN} Unknown grid code '{code}' in filename '{filename}'")
    return resolution


class TimeDuration(NamedTuple):
    """Decoded ISO 8601 duration with each component as a plain integer."""
    year:   int = 0
    month:  int = 0
    week:   int = 0
    day:    int = 0
    hour:   int = 0
    minute: int = 0
    second: int = 0


def decode_timeDuration(duration):
    """
    Decode time duration strings

    Parses ISO 8601 duration format (e.g. 'P1M', 'PT1H', 'P0000-01-00T00:00:00')
    and returns a TimeDuration with each component as an integer.
    Returns None if the string does not match either supported format.
    """
    regDuration1 = "^P(?:(?P<year>[0-9]+)Y)?(?:(?P<month>[0-9]+)M)?(?:(?P<week>[0-9]+)W)?(?:(?P<day>[0-9]+)D)?(?:T(?:(?P<hour>[0-9]+)?H)?(?:(?P<minute>[0-9]+)?M)?(?:(?P<second>[0-9]+)?S)?)?$"
    regDuration2 = "^P(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})T(?P<hour>[0-9]{2}):(?P<minute>[0-9]{2}):(?P<second>[0-9]+)$"

    timeDuration = re.match(regDuration1, duration)
    if timeDuration is None:
        timeDuration = re.match(regDuration2, duration)

    if timeDuration is None:
        return None

    def _int(group):
        v = timeDuration.group(group)
        return int(v) if v is not None else 0

    return TimeDuration(
        year   = _int('year'),
        month  = _int('month'),
        week   = _int('week') if 'week' in timeDuration.groupdict() else 0,
        day    = _int('day'),
        hour   = _int('hour'),
        minute = _int('minute'),
        second = _int('second'),
    )


def _parse_attribute_rule(attr_el) -> dict:
    """Parse a single <attribute> element into a rule dict.

    Shares the same content/regex sub-element structure as <entry> so the
    validation logic in _checkVariableAttributes can reuse the same patterns
    as _checkGlobalAttributes.

    The optional 'dtype' attribute (comma-separated token list) is validated
    by _matches_dtype(); it is independent of the existing 'type' attribute
    which preserves its original s/f32/f64 semantics for global-attribute rules.
    """
    return {
        'id':       normalize_whitespace(attr_el.get('id', '')),
        'type':     normalize_whitespace(attr_el.get('type', 's')),
        'dtype':    normalize_whitespace(attr_el.get('dtype', '')),
        'required': normalize_whitespace(attr_el.get('required', 'no')),
        'join':     attr_el.get('join', 'or').strip(),
        'content':  [
            {'value': normalize_whitespace(c.text or ''),
             'type':  c.get('type', 'string').strip()}
            for c in attr_el.findall('content')
        ],
        'regex':    [
            {'value': (r.text or '').strip(),
             'warn':  r.get('warn', '').strip(),
             'type':  r.get('type', 'string').strip()}
            for r in attr_el.findall('regex')
        ],
    }


def _parse_standard(filename: str) -> SimpleNamespace:
    """
    Parse a CM SAF metadata standard XML file.

    Returns a SimpleNamespace with attributes:
      .dict           — {id: entry_dict} for all <entry> elements
      .variables      — list of variable rule dicts for all <variable> elements
      .version_number — str
      .last_modified  — str
      .include        — str | None  (name of included base standard, if present)

    Each entry_dict has keys: type, required, evaluate, list, join,
    content (list of dicts), regex (list of dicts), keywords (list of str).

    Each variable rule dict has keys:
      match_by    — 'standard_name' or 'name'
      match_value — str
      required    — 'yes' | 'no'
      scope       — str (comma-separated list: 'all', 'grid', 'swath',
                    'regular_grid', 'native_grid', 'seviri_grid')
      coordinate  — bool: if True, only CF coordinate variables are matched
                    (1-D, name == dimension name per CF Conventions §1.2)
      attributes  — list of attribute rule dicts (see _parse_attribute_rule)
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    include_el = root.find('include')

    entry_dict = {}
    for entry in root.iter('entry'):
        eid = normalize_whitespace(entry.get('id', ''))
        entry_dict[eid] = {
            'type':     normalize_whitespace(entry.get('type',     's')),
            'required': normalize_whitespace(entry.get('required', 'no')),
            'evaluate': normalize_whitespace(entry.get('evaluate', 'no')),
            'list':     normalize_whitespace(entry.get('list',     '')),
            'join':     entry.get('join', 'or').strip(),
            'content':  [
                {'value': normalize_whitespace(c.text or ''),
                 'type':  c.get('type', 'string').strip()}
                for c in entry.findall('content')
            ],
            'regex':    [
                {'value': (r.text or '').strip(),
                 'warn':  r.get('warn', '').strip(),
                 'type':  r.get('type', 'string').strip()}
                for r in entry.findall('regex')
            ],
            'keywords': [
                normalize_whitespace(k.text or '')
                for k in entry.findall('keywords')
            ],
        }

    variable_list = []
    for var_el in root.findall('variable'):
        if var_el.get('standard_name') is not None:
            match_by    = 'standard_name'
            match_value = normalize_whitespace(var_el.get('standard_name'))
        elif var_el.get('name') is not None:
            match_by    = 'name'
            match_value = normalize_whitespace(var_el.get('name'))
        else:
            continue  # malformed element — skip silently
        variable_list.append({
            'match_by':    match_by,
            'match_value': match_value,
            'required':    normalize_whitespace(var_el.get('required',   'no')),
            'scope':       normalize_whitespace(var_el.get('scope',      'all')),
            'coordinate':  normalize_whitespace(var_el.get('coordinate', 'no')) == 'yes',
            'dtype':       normalize_whitespace(var_el.get('dtype',      '')),
            'attributes':  [_parse_attribute_rule(a) for a in var_el.findall('attribute')],
        })

    return SimpleNamespace(
        version_number = normalize_whitespace(root.findtext('version_number', '')),
        last_modified  = normalize_whitespace(root.findtext('last_modified',  '')),
        include        = normalize_whitespace(include_el.text) if include_el is not None else None,
        dict           = entry_dict,
        variables      = variable_list,
    )


def _find_file(filename: str, search_paths: list) -> str | None:
    """Return the first match for *filename* across *search_paths*, or None.

    Each entry in *search_paths* is checked in order. As a final fallback
    the filename is treated as a literal path and tested directly.
    """
    for path in search_paths:
        candidate = os.path.join(path, filename)
        if os.path.isfile(candidate):
            return candidate
    if os.path.isfile(filename):
        return filename
    return None


class Keywords:
    """
    class to read GCMD keywords
    """

    def __init__(self, filename, search_paths=None):
        self.filename = None
        self.groups = None
        self.keywordList = {}

        # resolve keyword file; fall back to literal path so readFile() can
        # report a meaningful IOError if it is genuinely missing
        self.filename = _find_file(filename, search_paths or []) or filename

    def readFile(self):
        try:
            fh = open(self.filename, "r")
        except IOError as detail:
            print(detail)
            return 2

        with fh:
            # csv reader
            reader = csv.reader(fh, delimiter=',', quotechar='"')

            # loop lines
            for row in reader:
                if not row:
                    continue
                if 'Keyword Version' in row[0]:
                    meta = re.match(r'.*Keyword Version:\s*([\d.]+)', row[0])
                    if meta:
                        self.keywordList['Version'] = meta.group(1)
                elif not self.groups and 'UUID' in row:
                    self.groups = row
                    self.groups.pop()
                elif self.groups:
                    uuid = row.pop().strip('"')
                    if len(row) == len(self.groups):
                        kw = {}
                        for index, item in enumerate(self.groups):
                            kw[item] = row[index].strip()
                        self.keywordList[uuid] = kw

        return 0

    def findKeyword(self, keyword="Version"):
        result = []

        keyword = keyword.upper()
        for key, kw in self.keywordList.items():
            if key == keyword:
                result.append(kw)
                continue

            # Search all columns derived from the CSV header, with one
            # special case: 'Term' is only matched when Variable_Level_1
            # is empty (i.e. the keyword lives at the Term level, not deeper).
            for col in self.groups or []:
                if col not in kw:
                    continue
                if col == 'Term' and kw.get('Variable_Level_1', '') != '':
                    continue
                if kw[col].upper() == keyword:
                    result.append(kw)
                    break

        return result

    def findKeywordList(self,keyword="Version"):
        result = self.findKeyword(keyword)
        if len(result) == 0:
            return

        # make list
        kwList = []
        for entry in result:
            list = []
            for item in self.groups:
                if len(entry[item]) > 0: list.append(entry[item])
            kwList.append(" > ".join(list))

        # return as list of strings
        return kwList


def _attributes_equal(ar, ac) -> bool:
    """Return True if two NetCDF attribute values are equal.

    Handles the three cases that arise in practice:
      - str           — plain equality
      - np.ndarray    — element-wise via np.array_equal, NaN-aware for floats
      - scalar numeric (np.float32/64, np.int*, etc.) — NaN-aware where applicable
    """
    if isinstance(ar, str):
        return ar == ac
    if isinstance(ar, np.ndarray):
        if np.issubdtype(ar.dtype, np.floating):
            nan_r, nan_c = np.isnan(ar), np.isnan(ac)
            if not np.array_equal(nan_r, nan_c):
                return False
            return np.array_equal(ar[~nan_r], ac[~nan_c])
        return np.array_equal(ar, ac)
    # scalar numeric — guard np.isnan with try/except since it raises
    # TypeError for integer types
    try:
        if np.isnan(ar) and np.isnan(ac):
            return True
    except (TypeError, ValueError):
        pass
    return bool(ar == ac)


def _is_coordinate_variable(var) -> bool:
    """Return True if *var* is a CF coordinate variable.

    A coordinate variable is one-dimensional and its name must equal its
    sole dimension name (CF Conventions §1.2).
    """
    return len(var.dimensions) == 1 and var.name == var.dimensions[0]


class DatasetX():
    """
    Expand standard python Dataset netcdf4 class
    """

    def __init__(self, *args, **kwargs):
        self._ds = Dataset(*args, **kwargs);


    def __getattribute__(self, item):
        if item in ["_ds", "getCoordinates", "getVariableByStandardName", "getVariableList", "isSwathData",
                "matchCoordinate", "matchCoordinateTime", "getVariableByName", "getvar", "getgrp"]:
            return object.__getattribute__(self, item)
        elif item == "ds":
            return object.__getattribute__(self, "_ds")
        else:
            return Dataset.__getattribute__(self._ds, item)


    def getvar(self, name):
        return self._ds[name]


    def getgrp(self, name):
        return self._ds[name]


    def close(self):
        if self._ds.isopen():
            self._ds.close()


    def isSwathData(self):
        """
        """
        rc = False

        if hasattr(self,"cdm_data_type"):
            if self.cdm_data_type == "swath":
                rc = True

        return rc


    def getVariableList(self):
        """
        Compile list of all variables
        """

        vList = []
        for v in self.variables:
            vList.append(os.path.join("/",v))
        for g in self.groups:
            for v in self.groups[g].variables:
                vList.append(os.path.join("/",g,v))

        return vList


    def getVariableByName(self, name):
        """
        Find variables by name in root and all data groups
        """

        result = {}

        for item in self.getVariableList():
            if os.path.basename(item) == name:
                result[item] = self._ds[item]

        return result


    def getVariableByStandardName(self, name):
        """
        Find variables by standard name in root and all data groups
        """

        result = {}

        varList = self.get_variables_by_attributes(standard_name=name)
        for item in varList:
            result[os.path.join(self.path,item.name)] = item

        for gName in self.groups:
            grp = self.groups[gName]
            varList = grp.get_variables_by_attributes(standard_name=name)
            for item in varList:
                result[os.path.join(grp.path,item.name)] = item

        return result


    def getCoordinates(self, standardName, shortName=[]):
        """
        Find coordinate variables.

        A candidate must be a CF coordinate variable: one-dimensional with its
        variable name equal to its dimension name (_is_coordinate_variable).
        """
        axis = self.getVariableByStandardName(standardName)
        if len(axis) == 0:
            for vName in shortName:
                axis = self.getVariableByName(vName)
                for item in list(axis.keys()):
                    if not _is_coordinate_variable(self._ds[item]):
                        axis.pop(item)
                        continue
                    print(f"## WARNING ##: missing standard name attribute for '{item}'")
        else:
            for item in list(axis.keys()):
                if not _is_coordinate_variable(self._ds[item]):
                    axis.pop(item)
                    continue
                var = self._ds[item]
                if hasattr(var, "long_name"):
                    # filter sub satellite geolocation
                    if re.match(r'^.*sub.satellite.*$', var.long_name):
                        axis.pop(item)

        return(axis)


    def matchCoordinate(self, dim, axis):
        """
        Find first matching coordinate
        """

        dPath = os.path.join(dim.group().path,dim.name)
        if dPath in axis.keys():
            return dPath
        else:
            return None


    def matchCoordinateTime(self, var, axisTime):
        """
        Find first matching time coordinate
        """

        keyTime = None
        iPath   = var.group().path
        while True:
            for _keyTime in axisTime.keys():
                _path = axisTime[_keyTime].group().path
                if _path == iPath:
                    keyTime = _keyTime
                    break
            if keyTime is not None or iPath == "/":
                break
            else:
                iPath = os.path.dirname(iPath)

        return(keyTime)


class CMSAFChecker:
    """
    CM SAF Checking class
    """

    def __init__(self, search_paths=None, version=None, referenceFile=None,
        coordinates=False, ignore=None, lazy=False, standard_file=None):

        self.search_paths  = search_paths or []
        self.standard_file = standard_file
        self.Dataset       = None
        self.testFn        = None
        self.refFile       = referenceFile
        self.refDataset    = None
        self.coordinates   = coordinates
        self.lazy          = lazy
        self.std_name_dh   = None
        self.gIgnoreAtt    = []
        self.vIgnoreAtt    = []
        self.err = 0
        self.errAttr = []
        self.warn = 0
        self.warnAttr = []
        self.info = 0
        self.infoAttr = []
        if version:
            self.version = "_v"+version.replace(".","-")
        else:
            self.version = ""

        # open reference NetCDF file
        if self.refFile is not None:
            try:
                print(f"Reference File: '{self.refFile}'\n")
                self.refDataset = DatasetX(self.refFile, mode='r')
            except Exception:
                print("\nCould not open reference file, please check that NetCDF is formatted correctly.\n")
                raise

        # attributes to ignore
        if self.refFile is not None:
            self.gIgnoreAtt = ["date_created", "time_coverage_start", "time_coverage_end", "filename", "date_modified", "history"]
        if ignore is not None:
            attList = ignore.split(",")
            for att in attList:
                if att.find('@') >= 0:
                    self.vIgnoreAtt.append(att)
                else:
                    self.gIgnoreAtt.append(att)

    def _loadStandard(self):
        """
        Load the XML metadata standard, searching each path.
        Each referenced file (main XML, included XML, GCMD keyword CSVs) is
        resolved independently through the same search_paths list so a newer
        version of any individual file in an earlier path takes precedence.
        """
        # locate the standard XML
        if self.standard_file is not None:
            fn = _find_file(self.standard_file, self.search_paths)
            if fn is None:
                print(f"Could not find '{self.standard_file}' in any search path: {self.search_paths}")
                exit(1)
        else:
            default_name = "cmsaf_metadata_standard" + self.version + ".xml"
            fn = _find_file(default_name, self.search_paths)
            if fn is None:
                print(f"Could not find '{default_name}' in any search path: {self.search_paths}")
                exit(1)

        print(f"Using standard file: '{fn}'")
        try:
            self.std_name_dh = _parse_standard(fn)
        except IOError as detail:
            print(detail)
            raise

        # locate and load any included XML, searched independently
        if self.std_name_dh.include:
            inc_fn = _find_file(self.std_name_dh.include, self.search_paths)
            if inc_fn is None:
                print(f"No such file '{self.std_name_dh.include}'.")
                exit(1)
            print(f"Including '{inc_fn}'.")
            try:
                include_ = _parse_standard(inc_fn)
            except IOError as detail:
                print(detail)
                raise

            # included file is the base; main file entries override
            merged = include_.dict
            merged.update(self.std_name_dh.dict)
            self.std_name_dh.dict = merged

            # variable rules: base first, child appends (no replacement)
            self.std_name_dh.variables = include_.variables + self.std_name_dh.variables


    def __del__(self):
        if self.refDataset:
            self.refDataset.close();


    def _reset(self):
        self.err = 0
        self.errAttr = []
        self.warn = 0
        self.warnAttr = []
        self.info = 0
        self.infoAttr = []


    def checker(self, file):
        """
        check wrapping procedure
        """

        # Load standard file for this file
        if self.refDataset is None:
            self._loadStandard()
            print(f"Using CM SAF Metadata Standard Version {self.std_name_dh.version_number} ({self.std_name_dh.last_modified})")

        # Check for valid filename
        fileSuffix = re.compile(r'^\S+\.nc$')
        if not fileSuffix.match(file):
            print(f"{RC_ERR} Filename must have '.nc' suffix")
            exit(1)

        # Read in netCDF file
        try:
            self.Dataset = DatasetX(file, mode='r')
            self.File    = os.path.basename(os.path.realpath(self.Dataset.filepath()))
        except RuntimeError as detail:
            print(f"{RC_ERR} ", detail)
            return 1
        except Exception:
            print(f"{RC_ERR} Could not open file, please check that NetCDF is formatted correctly.\n".upper())
            return 1

        # test compression
        print(f"\n{'':=^80}\n>>> checking compression\n{'':=^80}")
        rcCompress = self._checkCompression()
        if rcCompress == 0:
            print(f"\n{RC_OK} <<< compression")
        else:
            print(f"\n{RC_FAIL} <<< compression")

        # test variables
        print(f"\n{'':=^80}\n>>> checking variables\n{'':=^80}")
        rcVariables = self._checkVariables()
        if rcVariables == 0:
            print(f"\n{RC_OK} <<< variables")
        else:
            print(f"\n{RC_FAIL} <<< variables")

        # test against reference
        try:
            if self.refDataset is not None:
                print(f"\n{'':=^80}\n>>> checking metadata reference file\n{'':=^80}")
                rc = self._checkReferenceFile()
                if rc == 0:
                    print(f"\n{RC_OK} <<< metadata reference file")
                else:
                    print(f"\n{RC_FAIL} <<< metadata reference file")
            else:
                print(f"\n{'':=^80}\n>>> checking metadata standard\n{'':=^80}")
                rc = self._checkStandard()

                print(f"\n{'':=^80}\n>>> checking variable attributes\n{'':=^80}")
                rcVarAttr = self._checkVariableAttributes()
                if rcVarAttr == 0:
                    print(f"\n{RC_OK} <<< variable attributes")
                else:
                    print(f"\n{RC_FAIL} <<< variable attributes")
                rc += rcVarAttr

            if self.coordinates:
                print(f"\n{'':=^80}\n>>> Checking coordinates\n{'':=^80}")
                rcCoord = self._checkCoordinates()
                if rcCoord == 0:
                    print(f"\n{RC_OK} <<< coordinates")
                else:
                    print(f"\n{RC_FAIL} <<< coordinates")
                rc += rcCoord
        finally:
            self.Dataset.close()

        return rc+rcCompress+rcVariables


    def _matches_scope(self, scope_str: str) -> bool:
        """Return True if the current file matches any scope token in scope_str.

        Scope tokens (comma-separated):
          all          — always matches
          swath        — isSwathData() is True
          grid         — isSwathData() is False (regular or native projection)
          regular_grid — grid and grid code in _GRID_RESOLUTION
          native_grid  — grid and grid code in _PROJECTION_GRIDS
          polar_grid   — grid and grid code in _POLAR_GRIDS
          seviri_grid  — grid and grid code == '05'
        """
        scopes = [s.strip() for s in scope_str.split(',')]
        if 'all' in scopes:
            return True
        is_swath  = self.Dataset.isSwathData()
        decode    = re.match(CMSAF_NAMING_STANDARD, self.File)
        grid_code = decode.group(10) if decode else None
        for scope in scopes:
            if scope == 'swath'        and is_swath:                                         return True
            if scope == 'grid'         and not is_swath:                                     return True
            if scope == 'regular_grid' and not is_swath and grid_code in _GRID_RESOLUTION:   return True
            if scope == 'native_grid'  and not is_swath and grid_code in _PROJECTION_GRIDS:  return True
            if scope == 'polar_grid'   and not is_swath and grid_code in _POLAR_GRIDS:       return True
            if scope == 'seviri_grid'  and not is_swath and grid_code == '05':               return True
        return False


    def _checkOneVariableAttribute(self, var, vpath: str, rule: dict) -> int:
        """Check a single attribute rule against a variable.

        Returns 0 (pass) or 1 (fail). Mirrors the content/regex matching logic
        of _checkGlobalAttributes but without keyword or list support, which are
        not needed for variable attribute rules.
        """
        attr_id  = rule['id']
        required = rule['required']

        if not hasattr(var, attr_id):
            if required == 'yes':
                print(f"    {RC_ERR} missing required attribute '{attr_id}'")
                self.err += 1
                self.errAttr.append(f"{vpath}@{attr_id}")
                return 1
            return 0

        value = getattr(var, attr_id)
        print(f"    {attr_id}: {value!r}")

        # dtype check: verify the attribute's data type before content/regex checks
        if rule['dtype']:
            if not _matches_dtype(value, rule['dtype']):
                actual = 'str' if isinstance(value, str) else str(np.asarray(value).dtype)
                print(f"    {RC_ERR} incorrect data type for attribute '{attr_id}': '{actual}'")
                print(f"    Expecting dtype: {rule['dtype']}")
                self.err += 1
                self.errAttr.append(f"{vpath}@{attr_id}")
                return 1

        # content check: any matching entry is a pass (join=or default)
        if rule['content']:
            if not any(item['value'] == str(value) for item in rule['content']):
                print(f"    {RC_ERR} incorrect value for '{attr_id}': '{value}'")
                if len(rule['content']) == 1:
                    print(f"    Expecting: '{rule['content'][0]['value']}'")
                else:
                    print(f"    Expecting one of: {[i['value'] for i in rule['content']]}")
                self.err += 1
                self.errAttr.append(f"{vpath}@{attr_id}")
                return 1

        # regex check: pass if at least one non-warn pattern matches
        if rule['regex']:
            matched = False
            for item in rule['regex']:
                if re.search(item['value'], str(value)):
                    if item['warn']:
                        print(f"    {RC_WARN} {item['warn']}")
                        self.warn += 1
                        self.warnAttr.append(f"{vpath}@{attr_id}")
                    else:
                        matched = True
            all_warn = all(item['warn'] for item in rule['regex'])
            if not matched and not all_warn:
                print(f"    {RC_ERR} '{attr_id}' value '{value}' does not match expected pattern")
                self.err += 1
                self.errAttr.append(f"{vpath}@{attr_id}")
                return 1

        return 0


    def _checkVariableAttributes(self) -> int:
        """Check variable attribute rules defined in the XML standard.

        Iterates over std_name_dh.variables, resolves each rule to actual
        variables in the dataset, and validates their attributes.
        Returns 0 if all checks pass, 1 otherwise.
        """
        rc = 0
        ds = self.Dataset

        for var_rule in self.std_name_dh.variables:
            match_by    = var_rule['match_by']
            match_value = var_rule['match_value']
            required    = var_rule['required']
            scope_str   = var_rule['scope']
            is_coord    = var_rule['coordinate']

            # skip if file type is outside the declared scope
            if not self._matches_scope(scope_str):
                continue

            # resolve matching variables from the dataset
            if match_by == 'standard_name':
                found_vars = ds.getVariableByStandardName(match_value)
            else:
                found_vars = ds.getVariableByName(match_value)

            # when coordinate="yes", restrict to CF coordinate variables only
            # (1-D, variable name == its sole dimension name, CF Conventions §1.2)
            if is_coord:
                found_vars = {k: v for k, v in found_vars.items()
                              if _is_coordinate_variable(v)}

            if not found_vars:
                if required == 'yes':
                    print(f"  {RC_ERR} missing required variable with {match_by}='{match_value}'")
                    self.err += 1
                    self.errAttr.append(match_value)
                    rc = 1
                continue

            # check attribute rules for each matched variable
            for vpath, var in found_vars.items():
                print(f"\n  {vpath} ({match_by}='{match_value}')")

                # variable-level dtype check
                dtype_spec = var_rule['dtype']
                if dtype_spec and not _matches_dtype(var.dtype, dtype_spec):
                    print(f"  {RC_ERR} incorrect data type '{var.dtype}' for '{vpath}'")
                    print(f"  Expecting dtype: {dtype_spec}")
                    self.err += 1
                    self.errAttr.append(vpath)
                    rc = 1

                for attr_rule in var_rule['attributes']:
                    if self._checkOneVariableAttribute(var, vpath, attr_rule) != 0:
                        rc = 1

        return rc


    def _checkStandard(self):
        """
        check global metadata against CM SAF standard
        """
        # check global attributes
        rc = self._checkGlobalAttributes()

        print(f"\n{'':=^80}\nMetadata Summary\n{'':=^80}")
        print(f"ERRORs given: {self.err}")
        if self.err > 0:
            print("  in attributes: ", [s for s in self.errAttr])
        print(f"WARNINGS given: {self.warn}")
        if self.warn > 0:
            print("  in attributes: ", [s for s in self.warnAttr])
        print(f"INFORMATION messages: {self.info}")
        if self.info > 0:
            print("  in attributes: ", [s for s in self.infoAttr])

        return rc


    def _evaluatePlaceholders(self, attr: dict, ds) -> None:
        """Expand ${references} and ${year} placeholders in attr['content'] values in place.

        Called for every standard entry before any value checks so that the
        expanded values are visible in both the required-attribute pass and the
        content-matching pass.  The mutation is intentional and matches the
        original behaviour.
        """
        if attr['evaluate'] != 'yes':
            return
        if hasattr(ds, 'references'):
            ref_val = getattr(ds, 'references')
            for item in attr['content']:
                item['value'] = item['value'].replace('${references}', ref_val)
        now = os.environ.get('CMSAF_RELEASE_YEAR') or datetime.datetime.now().strftime('%Y')
        for item in attr['content']:
            item['value'] = item['value'].replace('${year}', now)

    def _resolveKeywordFilename(self, keywordsFn: str, ds) -> str:
        """Resolve a ${..._version} placeholder in a keyword CSV filename.

        Returns the filename with the version string substituted when the
        corresponding dataset attribute is present; otherwise returns the
        original (template) filename unchanged.
        """
        if '${' not in keywordsFn:
            return keywordsFn
        vocabulary_name = None
        vocabulary_version = None
        decode = re.match(r'^.*\$\{([a-z_]*)_version\}.*$', keywordsFn)
        if decode is not None:
            vocabulary_name = decode.group(1)
        if vocabulary_name is not None and hasattr(ds, vocabulary_name):
            decode = re.match(r'^.*Version +([0-9\.]*)$', getattr(ds, vocabulary_name))
            if decode is not None:
                vocabulary_version = decode.group(1)
        if vocabulary_name is not None and vocabulary_version is not None:
            keywordsFn = keywordsFn.replace('${' + vocabulary_name + '_version}', vocabulary_version)
        return keywordsFn

    def _getKeywords(self, raw_filename: str, key: str, kwList: dict):
        """Load (or retrieve from cache) the Keywords object for *raw_filename*.

        *raw_filename* may contain a ${..._version} placeholder; it is resolved
        against the current dataset before lookup.  The resolved filename is used
        as the cache key so the same CSV is never read twice.

        Returns the Keywords object on success, or None on failure (error already
        printed and self.err already incremented).
        """
        resolved = self._resolveKeywordFilename(raw_filename, self.Dataset)
        if resolved in kwList:
            return kwList[resolved]
        kw = Keywords(filename=resolved, search_paths=self.search_paths)
        if kw.readFile() != 0:
            print(f"{RC_ERR} Test incomplete")
            self.err += 1
            self.errAttr.append(key)
            return None
        kwList[resolved] = kw
        return kw

    def _checkKeywordValue(self, a: str, kw, key: str) -> bool:
        """Validate a single keyword path *a* against vocabulary *kw*.

        Returns True if the path is unambiguously found (caller should set
        or_passed).  On failure the appropriate error is printed and self.err
        is incremented.
        """
        entry_list = re.split(r' *> *', a)
        entry_item = ' > '.join(entry_list)

        kw_item = kw.findKeywordList(entry_list[-1])
        if not kw_item:
            print(f"{RC_ERR} '{entry_list[-1]}' not found as a keyword leaf in the vocabulary")
            self.err += 1
            if key not in self.errAttr:
                self.errAttr.append(key)
            return False

        # Build a match pattern anchored appropriately for the vocabulary type.
        # Short_Name-structured vocabularies (providers, instruments): a single
        # term may be the Short_Name at the start of the path, so match any
        # complete path segment.  Science keyword hierarchies: a single term
        # must be the leaf (rightmost element), so anchor to the end only.
        if len(entry_list) == 1:
            if 'Short_Name' in (kw.groups or []):
                entry_pattern = '(^| > )' + re.escape(entry_item) + '( > |$)'
            else:
                entry_pattern = '(^| > )' + re.escape(entry_item) + '$'
        else:
            entry_pattern = '.*' + re.escape(entry_item) + '$'

        matching_paths = [p for p in kw_item if re.search(entry_pattern, p) is not None]
        n_hits = len(matching_paths)

        if n_hits == 1:
            print(f"decoded as '{matching_paths[0]}'")
            return True
        if n_hits == 0:
            print(f"{RC_ERR} keyword path '{entry_item}' not found; "
                  f"valid paths containing '{entry_list[-1]}':")
            for p in kw_item:
                print(f"  {p}")
        else:
            print(f"{RC_ERR} keyword '{entry_item}' is ambiguous "
                  f"({n_hits} matches); be more specific, e.g.:")
            for p in matching_paths:
                print(f"  {p}")
        self.err += 1
        if key not in self.errAttr:
            self.errAttr.append(key)
        return False

    def _checkOneGlobalAttribute(self, key: str, attr, std: dict, fn: str, kwList: dict) -> None:
        """Check a single global attribute value against its standard rule.

        Handles type checking, empty-string reporting, list splitting,
        per-element content/regex/keyword validation, and the final
        join=or / join=and evaluation.  Updates self.err/warn/info counters
        and their Attr lists directly; does not return a value.
        """
        # --- type check ---
        attr_type = type(attr)
        if attr_type == type(np.array([])):
            if attr.size == 1 and type(attr[0]) == type(np.float32(1.0)):
                attr_type = 'f32'
            elif attr.size == 1 and type(attr[0]) == type(np.float64(1.0)):
                attr_type = 'f64'
        else:
            if attr_type == type(np.float64(1.0)):
                attr_type = 'f64'
            elif attr_type == type(np.float32(1.0)):
                attr_type = 'f32'
            elif attr_type == type(str('s')):
                attr_type = 's'

        if str.find(attr_type, std['type']) == -1:
            print(f"{RC_ERR} Incorrect attribute data type")
            print(f"Expecting: {std['type']}, found: {attr_type}")
            self.err += 1
            self.errAttr.append(key)

        # --- empty-string guard (string attributes only) ---
        if attr_type == 's':
            if len(attr) == 0:
                if std['required'] == 'yes':
                    if key in self.gIgnoreAtt:
                        print(f"{RC_INFO} Ignoring empty required attribute")
                        self.info += 1
                        if key not in self.infoAttr:
                            self.infoAttr.append(key)
                    else:
                        print(f"{RC_ERR} empty required attribute")
                        self.err += 1
                        if key not in self.errAttr:
                            self.errAttr.append(key)
                else:
                    print(f"{RC_INFO} empty attribute")
                    self.info += 1
                    if key not in self.infoAttr:
                        self.infoAttr.append(key)
                return

        # non-string attributes: print value and skip content/regex/keyword checks
        if attr_type != 's':
            print(attr)
            return

        # --- split attribute into a list if a delimiter is defined ---
        attr_list = []
        try:
            if std['list']:
                for row in csv.reader([attr], delimiter=std['list']):
                    attr_list = row
            else:
                attr_list = [attr]
        except UnicodeEncodeError as detail:
            print(f"{RC_ERR} {detail}")
            self.err += 1
            self.errAttr.append(key)

        # join=or:  did any validator pass for any list element?
        or_passed = False
        # join=and: which required content values have been seen?
        and_seen     = set()
        and_required = {item['value'] for item in std['content']}
        # tracks whether a content/regex check already reported an error
        # (used to suppress the redundant "missing a correct value" or message)
        key_rc = 0

        for a in attr_list:
            attr_match = []

            # strip surrounding whitespace and quotes
            a_stripped = a.strip().strip('"')
            if a != a_stripped:
                print(f"{RC_WARN} white spaces or quotes detected")
                self.warn += 1
                if key not in self.warnAttr:
                    self.warnAttr.append(key)
            a = a_stripped
            print(a)

            # --- content check ---
            if std['content']:
                matched_item = next(
                    (item for item in std['content'] if item['value'] == a),
                    None,
                )
                if matched_item is not None:
                    or_passed = True
                    and_seen.add(matched_item['value'])
                    attr_match.append(matched_item)
                else:
                    if key in self.gIgnoreAtt:
                        print(f"{RC_INFO} Ignoring incorrect attribute content '{key}'")
                        self.info += 1
                        if key not in self.infoAttr:
                            self.infoAttr.append(key)
                    else:
                        print(f"{RC_ERR} incorrect attribute content :: '{a}'")
                        if len(std['content']) == 1:
                            print(f"Expecting: '{std['content'][0]['value']}'")
                        key_rc = 1
                        self.err += 1
                        if key not in self.errAttr:
                            self.errAttr.append(key)

            # --- regex check ---
            if std['regex']:
                regex_matched = False
                for item in std['regex']:
                    if re.search(item['value'], a):
                        if item['warn']:
                            print(f"{RC_WARN} {item['warn']}")
                            self.warn += 1
                            if key not in self.warnAttr:
                                self.warnAttr.append(key)
                        else:
                            regex_matched = True
                            attr_match.append(item)
                if regex_matched:
                    or_passed = True
                else:
                    if key in self.gIgnoreAtt:
                        print(f"{RC_INFO} Ignoring incorrect attribute content '{key}'")
                        self.info += 1
                        if key not in self.infoAttr:
                            self.infoAttr.append(key)
                    else:
                        print(f"{RC_ERR} incorrect attribute content :: '{a}'")
                        key_rc = 1
                        self.err += 1
                        if key not in self.errAttr:
                            self.errAttr.append(key)

            # --- keyword check ---
            if std['keywords']:
                kw = self._getKeywords(std['keywords'][0], key, kwList)
                if kw is None:
                    continue
                for m_item in attr_match:
                    if m_item['type'] == 'keyword':
                        if self._checkKeywordValue(a, kw, key):
                            or_passed = True

        # --- evaluate hits ---
        if attr_list:
            if std['join'].lower() == 'or':
                if not or_passed and key_rc == 0:
                    print(f"{RC_ERR} missing a correct value for attribute '{key}'")
                    self.err += 1
                    if key not in self.errAttr:
                        self.errAttr.append(key)
            elif std['join'].lower() == 'and':
                for value in sorted(and_required - and_seen):
                    print(f"{RC_ERR} missing required specific attribute content :: '{value}'")
                    self.err += 1
                    if key not in self.errAttr:
                        self.errAttr.append(key)

    def _checkGlobalAttributes(self):
        """
        Check global attributes
        """

        print("\ncheck validity of global attributes.")
        ds = self.Dataset
        rc = 0

        fn = os.path.basename(os.path.realpath(ds.filepath()))
        kwList = {}

        # pass 1: expand placeholders and check required/none presence
        for key, attr in self.std_name_dh.dict.items():
            self._evaluatePlaceholders(attr, ds)
            if attr['required'] == 'yes' and not hasattr(ds, key):
                if key in self.gIgnoreAtt:
                    print(f"{RC_INFO} Ignoring missing required attribute '{key}'")
                    self.info += 1
                    self.infoAttr.append(key)
                else:
                    print(f"{RC_ERR} Missing required attribute '{key}'")
                    self.err += 1
                    self.errAttr.append(key)
            elif attr['required'] == 'none' and hasattr(ds, key):
                print(f"{RC_ERR} Found improper attribute '{key}'")
                self.err += 1
                self.errAttr.append(key)

        # pass 2: check values of attributes that are present
        for key in ds.ncattrs():
            attr = getattr(ds, key)

            if key == 'filename':
                print(f"\n{key}:\n{attr}")
                if fn != ds.filename:
                    print(f"{RC_ERR} incorrect file name :: '{ds.filename}'")
                    self.err += 1
                    self.errAttr.append(key)

            if key in self.std_name_dh.dict:
                print(f"\n{key}:")
                self._checkOneGlobalAttribute(key, attr, self.std_name_dh.dict[key], fn, kwList)

        if self.err > 0:
            rc = 1

        return rc


    def _checkReferenceAttributes(self, new, ref, parent, ignore=None):
        """
        check all attributes in group against refence file
        """
        rc = 0
        attCheck = {}

        if ignore is None:
            ignore=self.gIgnoreAtt

        # compare attributes from reference file
        for attName in ref.ncattrs():
            if parent != "/":
                attNameFull = parent+"@"+attName
            else:
                attNameFull = attName

            # just test if attribute is there
            if not hasattr(new, attName):
                if attName in ignore:
                    print(f"{RC_INFO} Missing attribute :: '{attNameFull}'")
                else:
                    print(f"{RC_ERR} Missing attribute :: '{attNameFull}'")
                    rc = 1
            else:
                # mark as read
                attCheck[attName] = 1

                # file name
                if attName == 'filename':
                    if self.File != new.filename:
                        print(f"{RC_ERR} incorrect file name :: '{new.filename}'")
                        rc = 1
                    continue

                # skip attributes that are allowed to change
                elif attName in ignore:
                    print(f"{RC_INFO} changing attribute {attNameFull} :: '{getattr(new,attName)}'")
                    continue

                # test if attributes are identical
                else:
                    ar = getattr(ref, attName)
                    ac = getattr(new, attName)
                    if not _attributes_equal(ar, ac):
                        print(f"{RC_ERR} attribute '{attNameFull}' differ")
                        print(f"{'':<4}{attNameFull} :: expecting '{ar}', found '{ac}'")
                        rc = 1

        # check for new attributes
        for attName in new.ncattrs():
            if parent != "/":
                attNameFull = parent+"@"+attName
            else:
                attNameFull = attName

            if not attName in attCheck:
                # skip attributes that are allowed to change
                if attName in ignore:
                    print(f"{RC_INFO} changing new attribute {attNameFull} :: '{getattr(new,attName)}'")
                else:
                    print(f"{RC_ERR} New attribute {attNameFull} :: '{getattr(new,attName)}'")
                    rc = 1

        return rc


    def _checkReferenceVariables(self, new, ref):
        """
        check all variables in group against refence file
        """
        rc = 0
        varCheck = {}

        # compare variables from reference file
        for varName in ref.variables:
            # test if variable exists
            if not varName in new.variables:
                print(f"{RC_ERR} missing variable :: '{varName}'")
                rc = 1
            else:
                # mark as read
                varCheck[varName] = 1

                # define variables
                vr = ref.variables[varName]
                vc = new.variables[varName]

                # list of attributes to ignore
                ignore = []
                for item in self.vIgnoreAtt:
                    if item.find('@') >= 0:
                        itemList = item.split("@")
                        if itemList[0] == varName or len(itemList[0]) == 0 :
                            ignore.append(itemList[1])

                # check data type
                if vr.dtype != vc.dtype:
                    print(f"{RC_ERR} type of variable '{varName}' differs from reference.")
                    print(f"## ref='{vr.dtype}', file='{vc.dtype}'")
                    rc = 1

                # check array shape
                if not self.refDataset.isSwathData() and vr.shape != vc.shape:
                    print(f"{RC_ERR} shape of variable '{varName}' differs from reference.")
                    print(f"## ref='{vr.shape}', file='{vc.shape}'")
                    rc = 1

                # check variable attributes
                if self._checkReferenceAttributes(vc, vr, varName, ignore=ignore) > 0:
                    rc = 1

        # check for new variables
        for varName in new.variables:
            if not varName in varCheck:
                print(f"{RC_ERR} new variable '{varName}'")
                rc = 1

        return rc


    def _checkReferenceFile(self):
        """
        Check metadata from a file against a reference file.
        """

        rc = 0

        # get file name
        fn = os.path.realpath(self.Dataset.filepath())
        fn = os.path.basename (fn)

        # process global attributes
        print("global attributes")

        # track global attributes
        grpCheck = {}

        # check root group attributes
        rc = self._checkReferenceAttributes(self.Dataset, self.refDataset, "/")

        # check group attributes
        print("\ngroup attributes")
        for grpName in self.refDataset.groups:
            if grpName in self.Dataset.groups:
                rc = self._checkReferenceAttributes(self.Dataset.groups[grpName], self.refDataset.groups[grpName], grpName)
                grpCheck[grpName] = 1
            else:
                print(f"{RC_ERR} missing group '{grpName}'")

        # check for new groups
        for grpName in self.Dataset.groups:
            if not grpName in grpCheck:
                print(f"{RC_ERR} new group '{grpName}'")
                rc = 1

        # print global attribute check result
        if rc == 0:
            print(f"\n{RC_OK} <<< reference attributes")
        else:
            print(f"\n{RC_FAIL} <<< reference attributes")
        rcAll = rc

        # process variables
        print("\nvariables")
        rc = self._checkReferenceVariables(self.Dataset, self.refDataset)

        # print result
        if rc == 0:
            print(f"\n{RC_OK} <<< variables")
        else:
            print(f"\n{RC_FAIL} <<< variables")

        # return to caller
        return rcAll + rc


    def _checkCoordinates(self):
        """
        Check coordinates of a netcdf file.
        """

        rc = 0
        ds = self.Dataset
        tests = {}

        # decode duration
        timeDuration   = None
        timeResolution = None
        if hasattr(ds,"time_coverage_duration"):
            timeDuration = decode_timeDuration(ds.time_coverage_duration)
        if hasattr(ds, "time_coverage_resolution"):
            timeResolution = decode_timeDuration(ds.time_coverage_resolution)

        # decode stats and cycle
        decode = re.match(CMSAF_NAMING_STANDARD, self.File)
        expRecords = None
        expClimate = False
        expTimeDuration   = None
        expTimeResolution = None

        if decode is None:
            print(f"{RC_WARN} filename '{self.File}' does not match CM SAF naming convention, skipping filename-derived checks")
        else:
            # test for diurnal cycle
            if decode.group(3) == 'd':
                expRecords = 24
                if timeResolution is None:
                    timeResolution = decode_timeDuration("PT1H")
                elif timeResolution != decode_timeDuration("PT1H"):
                    print(f"{RC_ERR} ## expecting 'PT1H' as time_coverage_resolution for diurnal cycle")
                    tests['time'] = 1
                    timeResolution = decode_timeDuration("PT1H")
            # monthly climatology
            elif decode.group(2) == 'c' and decode.group(3) == 'm':
                expRecords = 12
            elif decode.group(2) != 'i':
                expRecords = 1

            # get expected number of records from time_coverage_resolution and time_coverage_duration
            if decode.group(2) == 'i':
                if (timeDuration is not None) and (timeResolution is not None):
                    td = datetime.timedelta(days=timeDuration.day,   hours=timeDuration.hour,   minutes=timeDuration.minute,   seconds=timeDuration.second)
                    tr = datetime.timedelta(days=timeResolution.day, hours=timeResolution.hour, minutes=timeResolution.minute, seconds=timeResolution.second)
                    expRecords = divmod(td.total_seconds(), tr.total_seconds())
                    expRecords = int(expRecords[0])

            # expect climate attribute
            if decode.group(2) == 'c' or decode.group(3) == 'd':
                expClimate = True

            # expecting duration and resolution
            if decode.group(2) == 'm':
                expTimeDuration   = ('P1M','P0000-01-00T00:00:00')
                expTimeResolution = ('P1M','P0000-01-00T00:00:00')
            if decode.group(2) == 'd':
                expTimeDuration   = ('P1D', 'P0000-00-01T00:00:00')
                expTimeResolution = ('P1D', 'P0000-00-01T00:00:00')
            if decode.group(3) == 'd':
                expTimeResolution = ('PT1H', 'P0000-00-00T01:00:00')
            if decode.group(2) == 'h':
                expTimeResolution = ('PT1H', 'P0000-00-00T01:00:00')
                expTimeDuration   = ('P1D', 'P0000-00-01T00:00:00', 'PT1H', 'P0000-00-00T01:00:00')
                if (timeDuration is not None) and (timeResolution is not None):
                    td = datetime.timedelta(days=timeDuration.day,   hours=timeDuration.hour,   minutes=timeDuration.minute,   seconds=timeDuration.second)
                    tr = datetime.timedelta(days=timeResolution.day, hours=timeResolution.hour, minutes=timeResolution.minute, seconds=timeResolution.second)
                    expRecords = divmod(td.total_seconds(), tr.total_seconds())
                    expRecords = int(expRecords[0])

        # find possible time axis
        axisTime = ds.getCoordinates("time", shortName=["time"])

        # find record_status variables
        tmp = ds.getVariableByName("record_status")
        recordStatus = {}
        print("\n>>> record_status")

        if len(tmp) == 0:
            if ds.isSwathData():
                print(f"{RC_INFO} record_status variable not required for 'swath' data")
                tests['record_status'] = 0
            else:
                print(f"{RC_ERR} missing record_status variable")
                tests['record_status'] = 1
        else:
            tests['record_status'] = 0

        # test all record status variables
        for key in tmp.keys():
            print(f"{'':<4}{key}")
            item = tmp[key]
            recordStatus[key] = { "var": item, "dict": None, "val": None, "mask": None, "time": None, }

            if hasattr(item,'flag_meanings') and hasattr(item,'flag_values'):
                recordStatus[key]["dict"] = dict(zip(item.flag_values,item.flag_meanings.split(" ")))
                recordStatus[key]["val"]  = item[:]
                recordStatus[key]["mask"] = np.ma.getmaskarray(item[:])
            else:
                print(f"{'':<4}{RC_ERR} missing valid variable '{key}'")
                tests['record_status'] = 1
                rc = 1

            # select valid time variable
            keyTime = ds.matchCoordinateTime(item, axisTime)
            if keyTime is not None:
                recordStatus[key]["time"] = keyTime
                if os.path.basename(keyTime) not in item.dimensions:
                    print(f"{'':<4}{RC_ERR} missing time dimension for '{key}'")
                    tests['record_status'] = 1
                    rc = 1
            else:
                print(f"{'':<4}{RC_ERR} missing valid time for '{key}'")
                tests['record_status'] = 1
                rc = 1

        # print result
        if tests['record_status'] == 0:
            print(f"\n{RC_OK} <<< record_status")
        else:
            print(f"\n{RC_FAIL} <<< record_status")
            rc = 1

        print("\n>>> time")

        # test time coordinates
        tests['time'] = 0
        if len(axisTime) == 0:
            print(f"{RC_ERR} missing time variable")
            tests['time'] = 1

        for vTime in axisTime.keys():
            print(f"\n{'':<4}{vTime}")

            itRc = self._checkCoordinatesTime(axisTime[vTime], expClimate=expClimate, expRecords=expRecords, \
                    expResolution=timeResolution, recordStatus=recordStatus)

            # print result
            if itRc == 0:
                print(f"\n{'':<4}{RC_OK} <<< coordinate '{vTime}'")
            else:
                print(f"\n{'':<4}{RC_FAIL} <<< coordinate '{vTime}'")
                tests['time'] = 1

        # time related global attributes
        if expTimeDuration is not None and hasattr(ds,"time_coverage_duration"):
            if ds.time_coverage_duration not in expTimeDuration:
                print(f"\n{'':<4}{RC_ERR} Unexpected time_coverage_duration '{ds.time_coverage_duration}'")
                tests['time'] = 1
            else:
                print(f"\n{'':<4}time coverage duration: '{ds.time_coverage_duration}'")

        if expTimeResolution is not None and hasattr(ds,"time_coverage_resolution"):
            if ds.time_coverage_resolution not in expTimeResolution:
                print(f"{'':<4}{RC_ERR} Unexpected time_coverage_resolution '{ds.time_coverage_resolution}'")
                tests['time'] = 1
            else:
                print(f"{'':<4}time coverage resolution: '{ds.time_coverage_resolution}'")

        # print final time check result
        if tests['time'] == 0:
            print(f"\n{'':<4}{RC_OK} <<< time")
        else:
            print(f"\n{'':<4}{RC_FAIL} <<< time")
            rc = 1

        resFile = cmsaf_decode_grid (self.File)

        # test latitude
        tests['lat'] = 0
        axisLat = ds.getCoordinates("latitude", shortName=["lat","latitude"])

        print("\n>>> latitude")

        if len(axisLat) == 0:
            if resFile is None:
                print(f"{'':<4}{RC_INFO} no regular grid expected from filename, skipping latitude coordinate check")
                tests['lat'] = 0
            else:
                print(f"{'':<4}{RC_ERR} missing latitude coordinate")
                tests['lat'] = 1
                rc = 1

        for vLat in axisLat.keys():
            print(f"\n{'':<4}{vLat}")

            keyTime = ds.matchCoordinateTime(axisLat[vLat], axisTime)
            itRc = self._checkCoordinatesGeo(axisLat[vLat], axisTime[keyTime], shortName="lat", longName="latitude", expAxis="Y")

            # print result
            if itRc == 0:
                print(f"\n{'':<4}{RC_OK} <<< coordinate '{vLat}'")
            else:
                print(f"\n{'':<4}{RC_FAIL} <<< coordinate '{vLat}'")
                tests['lat'] = 1

        # print final latitude check result
        if tests['lat'] == 0:
            print(f"\n{'':<4}{RC_OK} <<< latitude")
        else:
            print(f"\n{'':<4}{RC_FAIL} <<< latitude")
            rc = 1

        ## test longitude ##
        tests['lon'] = 0
        axisLon = ds.getCoordinates("longitude", shortName=["lon","longitude"])

        print("\n>>> longitude")

        if len(axisLon) == 0:
            if resFile is None:
                print(f"{'':<4}{RC_INFO} no regular grid expected from filename, skipping longitude coordinate check")
                tests['lon'] = 0
            else:
                print(f"{'':<4}{RC_ERR} missing longitude coordinate")
                tests['lon'] = 1
                rc = 1

        for vLon in axisLon.keys():
            print(f"\n{'':<4}{vLon}")

            keyTime = ds.matchCoordinateTime(axisLon[vLon], axisTime)
            itRc = self._checkCoordinatesGeo(axisLon[vLon], axisTime[keyTime], shortName="lon", longName="longitude", expAxis="X")

            # print result
            if itRc == 0:
                print(f"\n{'':<4}{RC_OK} <<< coordinate '{vLon}'")
            else:
                print(f"\n{'':<4}{RC_FAIL} <<< coordinate '{vLon}'")
                tests['lon'] = 1

        # print final longitude check result
        if tests['lon'] == 0:
            print(f"\n{'':<4}{RC_OK} <<< longitude")
        else:
            print(f"\n{'':<4}{RC_FAIL} <<< longitude")
            rc = 1

        print("\nsummary coordinates")
        tmp=["OK", "FAILED"]
        for key in tests:
            print(f"{'':<4}{key} :: {tmp[tests[key]]}", end='')
        print('')

        return rc


    def _checkTimeAxis(self, timeC, expClimate: bool) -> tuple:
        """Validate the axis, climatology, and units attributes on a time coordinate variable.

        Returns (rc, tUnits, sinceYr, calendar):
          rc       — 0 OK, -1 axis/range error (time steps may still be decoded),
                     1 units missing or unparseable (skip time step decoding)
          tUnits   — str value of the 'units' attribute, or None if absent/invalid
          sinceYr  — int year extracted from units, or None if units invalid
          calendar — str calendar name (default 'standard'); only set when rc <= 0
        """
        rc = 0
        tUnits   = None
        sinceYr  = None
        calendar = "standard"

        # test axis attribute
        if not hasattr(timeC, "axis"):
            print(f"{'':<8}{RC_ERR} missing mandatory attribute 'axis'")
            rc = -1
        elif timeC.axis != "T":
            print(f"{'':<8}{RC_ERR} invalid value attribute 'axis={timeC.axis}'")
            rc = -1

        # test for climate bounds attribute
        if expClimate:
            if not hasattr(timeC, 'climatology'):
                print(f"{'':<8}{RC_WARN} Expecting attribute 'climatology' as attribute for time bounds time")
        else:
            if hasattr(timeC, 'climatology'):
                print(f"{'':<8}{RC_WARN} Found unexpected attribute 'climatology' as attribute for time bounds time")

        # decode time units
        if not hasattr(timeC, 'units'):
            print(f"{'':<8}{RC_ERR} missing mandatory attribute 'units' for time axis")
            rc = 1
        else:
            tUnits = timeC.units
            m = re.match('(.*since )((-)?[-0-9]{4})(.*)', tUnits)
            if m is None:
                print(f"{'':<8}{RC_ERR} invalid time axis: '{tUnits}'")
                rc = 1
            else:
                sinceYr = int(m.group(2))
                if sinceYr < 1958:
                    print(f"{'':<8}{RC_ERR} invalid time axis: '{tUnits}'")
                    rc = -1

        # calendar is only extracted when units are valid
        if rc <= 0 and hasattr(timeC, 'calendar'):
            calendar = timeC.calendar

        return rc, tUnits, sinceYr, calendar


    def _decodeTimeSteps(self, timeC, sinceYr: int, calendar: str) -> tuple:
        """Decode the numeric time axis values to a numpy array of datetime objects.

        For swath data only the first and last record are decoded.
        Returns (tSteps, rc): tSteps is None on unrecoverable error.
        """
        rc = 0
        ds = self.Dataset

        axisTmp = timeC
        if ds.isSwathData():
            axisTmp = np.empty(2, dtype=timeC.dtype)
            axisTmp[0] = timeC[0]
            axisTmp[1] = timeC[-1]
            print(f"{'':<8}{RC_INFO} checking only first an last record for swath data.")

        tSteps = None
        try:
            tSteps = np.empty(axisTmp.shape, dtype=datetime.datetime)
            it = np.nditer(axisTmp, flags=['c_index'])
            while not it.finished:
                if sinceYr < 1:
                    t = Time(axisTmp[it.index], format='jd', scale='utc', precision=4)
                    t = datetime.datetime.strptime(t.isot, "%Y-%m-%dT%H:%M:%S.%f")
                    if t is not None:
                        tSteps[it.index] = t.replace(tzinfo=pytz.utc)
                else:
                    try:
                        t = num2date(axisTmp[it.index], timeC.units, calendar=calendar)
                        if isinstance(t, datetime.datetime):
                            t = t.replace(tzinfo=pytz.utc)
                        tmp = np.around(t.microsecond * np.float64(0.01)).astype(np.int64) * 100
                        if t.microsecond > 0:
                            print(f"{'':<8}{RC_WARN} time record not exact (mus={t.microsecond})")
                        t = t + datetime.timedelta(microseconds=int(tmp - t.microsecond))
                    except ValueError:
                        print(f"{'':<8}{RC_ERR} invalid time record")
                        rc = 1
                        t = None

                    if t is not None:
                        tSteps[it.index] = t
                it.iternext()
        except Exception:
            rc = 1
            tSteps = None
            print(f"{'':<8}{RC_ERR} invalid time axis.")

        return tSteps, rc


    def _decodeTimeBounds(self, timeC, tUnits, sinceYr, calendar: str) -> tuple:
        """Locate, validate, and decode the bounds or climatology variable for a time axis.

        Also determines whether missing bounds (due to a shape error) constitute
        an error, since that decision is tightly coupled to the variable lookup.

        Returns (timeBoundsKey, timeBounds, rc):
          timeBoundsKey — 'bounds' or 'climatology'
          timeBounds    — decoded datetime array, or None
          rc            — 0 OK, 1 on any error
        """
        rc = 0
        ds = self.Dataset

        timeBoundsKey = "bounds"
        timeBoundsVar = None
        timeBounds    = None

        if tUnits is None or sinceYr is None:
            print(f"{'':<8}{RC_ERR} cannot decode time bounds: time units are invalid")
            return timeBoundsKey, None, 1

        if hasattr(timeC, "climatology"):
            timeBoundsKey = 'climatology'

        if hasattr(timeC, timeBoundsKey):
            tmp = getattr(timeC, timeBoundsKey)
            if timeC.group().name == "/":
                if tmp in ds.variables:
                    timeBoundsVar = ds.getvar(tmp)
            else:
                if tmp in ds.getgrp(timeC.group().path).variables:
                    timeBoundsVar = ds.getvar(os.path.join(timeC.group().path, tmp))
            if timeBoundsVar is None:
                print(f"{'':<8}{RC_ERR} Missing configured bounds variable '{tmp}'.")
                rc = 1

        if timeBoundsVar is not None:
            if hasattr(timeBoundsVar, 'units'):
                if timeBoundsVar.units != tUnits:
                    print(f"{'':<8}{RC_ERR} time bounds must have same axis as time, but found: {timeBoundsVar.units}")
                    rc = 1
            if timeBoundsVar.shape != (timeC[:].size, 2):
                print(f"{'':<8}{RC_ERR} time bounds must have shape (size_of_time,2), but found: {timeBoundsVar.shape}")
                rc = 1
            else:
                timeBounds = np.empty(timeBoundsVar.shape, dtype=datetime.datetime)
                it = np.nditer(timeBoundsVar, flags=['multi_index'])
                while not it.finished:
                    if sinceYr < 1:
                        t = Time(it[0], format='jd', scale='utc', precision=4)
                        t = datetime.datetime.strptime(t.isot, "%Y-%m-%dT%H:%M:%S.%f")
                        timeBounds[it.multi_index] = t.replace(tzinfo=pytz.utc)
                    else:
                        try:
                            t = num2date(it[0], tUnits, calendar=calendar)
                            if isinstance(t, datetime.datetime):
                                t = t.replace(tzinfo=pytz.utc)
                            timeBounds[it.multi_index] = t
                        except OverflowError:
                            print(f"{'':<8}{RC_ERR} Invalid time bound record value: {it[0]}")
                            timeBounds = None
                            rc = 1
                            break
                    it.iternext()

            # determine whether missing bounds (shape error) constitute an error
            if timeBounds is None:
                expTimeBounds = True

                # no bounds required for swath data
                if ds.isSwathData():
                    expTimeBounds = False

                # no bounds required if all variable_id variables use cell_methods: point
                elif hasattr(ds, "variable_id"):
                    haveCellMethodsPoint = True
                    for item in ds.variable_id.split(","):
                        hasCellMethodsPoint = False
                        if item in ds.variables:
                            if hasattr(ds.variables[item], "cell_methods"):
                                if re.match('^.*time: *point.*$', ds.variables[item].cell_methods):
                                    hasCellMethodsPoint = True
                        if not hasCellMethodsPoint:
                            haveCellMethodsPoint = False
                            break
                    if haveCellMethodsPoint:
                        expTimeBounds = False

                if expTimeBounds:
                    print(f"{'':<8}{RC_ERR} Missing time bounds")
                    if not self.lazy:
                        rc = 1
                    else:
                        print(f"{'':<8}Ignoring while beeing lazy")
                else:
                    print(f"{'':<8}{RC_INFO} No time bounds required ")

        return timeBoundsKey, timeBounds, rc

    def _checkTimeRecords(self, tSteps, timeBounds, timeBoundsKey, timeDuration,
                          expResolution, recordStatus, timeC) -> int:
        """Validate each time record against bounds, expected resolution, and record_status.

        Checks per record:
          - expected record position (from expResolution arithmetic, including leap-year pentad)
          - record_status flag value
          - time bounds containment, gap, and overlap (when timeBounds is available)
          - coverage continuity between successive records

        Returns rc: 0 OK, 1 on any failure.
        """
        rc = 0
        prevTimeEnd    = None
        timeBoundRight = None
        itRecord       = tSteps[0]
        nRecord        = None

        it = np.nditer(tSteps, flags=['c_index', 'refs_ok'])
        while not it.finished:
            itRc = 'OK'
            itRecStatus = None

            # update current expected record
            if nRecord is not None:
                itRecord = nRecord

            # compute next expected record from resolution
            if expResolution is not None:
                nRecord = itRecord
                if expResolution.year > 0:
                    nRecord = nRecord.replace(year=nRecord.year + expResolution.year)
                if expResolution.month > 0:
                    nm = nRecord.month + expResolution.month
                    ny = nRecord.year
                    if nm > 12:
                        ny += 1
                        nm -= 12
                    nRecord = nRecord.replace(year=ny, month=nm)
                tr = datetime.timedelta(days=expResolution.day, hours=expResolution.hour,
                                        minutes=expResolution.minute, seconds=expResolution.second)
                nRecord += tr

                # add one day during a leap year for the last February pentad period
                if expResolution.day == 5:
                    if nRecord.timetuple().tm_yday == 61 and cal.isleap(nRecord.year):
                        nRecord += datetime.timedelta(days=1)

            # decode matching record_status flag
            rsKey = os.path.join(timeC.group().path, 'record_status')
            if rsKey in recordStatus:
                rsItem = recordStatus[rsKey]
                if len(rsItem["var"].shape) == 0:
                    recordStatusValC  = rsItem["val"][0]
                    recordStatusMaskC = rsItem["mask"][0]
                else:
                    recordStatusValC  = rsItem["val"][it.index]
                    recordStatusMaskC = rsItem["mask"][it.index]

                if not recordStatusMaskC and recordStatusValC in rsItem["dict"]:
                    itRecStatus = rsItem["dict"][recordStatusValC]
                else:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} invalid record_status value [{recordStatusValC}]")
                    itRecStatus = None

            # compute effective right bound, adjusting for climatology duration
            if timeBounds is not None:
                timeBoundRight = timeBounds[it.index, 1]
                if timeBoundsKey == 'climatology' and timeDuration is not None:
                    lb = timeBounds[it.index, 0]
                    lb = datetime.datetime(year=lb.year, month=lb.month, day=lb.day,
                                           hour=lb.hour, minute=lb.minute, second=lb.second)
                    rb = timeBounds[it.index, 1]
                    rb = datetime.datetime(year=rb.year, month=rb.month, day=rb.day,
                                           hour=rb.hour, minute=rb.minute, second=rb.second)
                    num_days = 0

                    # years
                    if timeDuration.year > 1:
                        et = lb + relativedelta(years=timeDuration.year - 1)
                        num_days = (et - lb).days
                        if cal.isleap(et.year) and (rb - et).days == 29:
                            num_days += 1

                    # months
                    if timeDuration.month > 0:
                        et = lb + relativedelta(months=timeDuration.month - 1)
                        num_days = (et - lb).days + cal.monthrange(et.year, et.month)[1] - 1

                    # days
                    num_days += max(timeDuration.day - 1, 0)

                    # sub-day remainder
                    rest = (relativedelta(hours=max(timeDuration.hour - 1, 0))
                            + relativedelta(minute=max(timeDuration.minute - 1, 0))
                            + relativedelta(second=max(timeDuration.second - 1, 0)))

                    timeBoundRight = rb - (relativedelta(days=num_days) + rest)

            # test record against time bounds
            if timeBounds is not None:
                if (timeBounds[it.index, 0] > tSteps[it.index]) or \
                   (timeBoundRight          < tSteps[it.index]):
                    itRc = 'FAILED (record not in bounds)'
                    rc = 1
                if nRecord is not None:
                    if timeBoundRight < nRecord:
                        itRc = 'FAILED (gap in right bound)'
                        rc = 1
                    if timeBoundRight > nRecord:
                        itRc = 'FAILED (overlap in right bound)'
                        rc = 1
                if itRecStatus is not None:
                    itRc = itRc + ' [status=' + itRecStatus + ']'
                    itRecStatus = None
                print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} "
                      f"[{timeBounds[it.index,0].isoformat()}, {timeBounds[it.index,1].isoformat()}]"
                      f" -> {itRc}")

                # test coverage continuity
                if prevTimeEnd is not None:
                    timeDiff = (timeBounds[it.index, 0] - prevTimeEnd).total_seconds()
                    if timeDiff > 0.:
                        rc = 1
                        print(f"{'':<8}{RC_ERR} gap in time coverage {timeDiff} seconds")
                    elif timeDiff < 0.:
                        rc = 1
                        print(f"{'':<8}{RC_ERR} overlap in time coverage {timeDiff} seconds")

            # test record against expected resolution (no bounds available)
            elif expResolution is not None:
                timeDiff = itRecord - tSteps[it.index]
                if timeDiff.total_seconds() > 0.:
                    itRc = 'FAILED'
                    rc = 1
                if itRecStatus is not None:
                    itRc = itRc + ' [status=' + itRecStatus + ']'
                    itRecStatus = None
                if timeDiff.total_seconds() == 0.:
                    print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} -> {itRc}")
                else:
                    print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} expected {itRecord} -> {itRc}")

            # bare record listing (neither bounds nor resolution available)
            else:
                if tSteps[it.index] is None:
                    itRc = 'FAILED'
                if itRecStatus is not None:
                    itRc = itRc + ' [status=' + itRecStatus + ']'
                print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} -> {itRc}")

            if timeBounds is not None and timeBoundRight is not None:
                prevTimeEnd = timeBoundRight

            it.iternext()

        return rc


    def _checkTimeCoverage(self, tSteps, timeBounds, timeBoundsKey) -> tuple:
        """Check time_coverage_start and time_coverage_end global attributes.

        Validates both against the first/last decoded time records and, where
        applicable, against the first/last time bounds.

        Returns (rc, timeCoverStart, timeCoverEnd).
        """
        rc = 0
        ds = self.Dataset
        timeCoverStart = None
        timeCoverEnd   = None

        if hasattr(ds, "time_coverage_start"):
            try:
                t = num2date(0, "days since {}".format(ds.time_coverage_start))
                if isinstance(t, datetime.datetime):
                    t = t.replace(tzinfo=pytz.utc, microsecond=0)
                timeCoverStart = t
            except ValueError:
                print(f"{'':<8}{RC_ERR} Unexpected time format: '{ds.time_coverage_start}'")
            else:
                if timeCoverStart > tSteps[0]:
                    print(f"{'':<8}{RC_ERR} first time record '{tSteps[0].isoformat()}' "
                          f"not within time_coverage_start attribute: '{timeCoverStart.isoformat()}'")
                    rc = 1
                if timeBounds is not None and timeCoverStart != timeBounds[0, 0]:
                    print(f"{'':<8}{RC_ERR} time bound [0,0] is not matching "
                          f"time_coverage_start attribute: '{timeCoverStart.isoformat()}'")
                    rc = 1

        if hasattr(ds, "time_coverage_end"):
            try:
                t = num2date(0, "days since {}".format(ds.time_coverage_end))
                if isinstance(t, datetime.datetime):
                    t = t.replace(tzinfo=pytz.utc, microsecond=0)
                timeCoverEnd = t
            except ValueError:
                print(f"{'':<8}{RC_ERR} Unexpected time format: {ds.time_coverage_end}")
            else:
                if timeCoverEnd < tSteps[-1]:
                    print(f"{'':<8}{RC_ERR} last time record '{tSteps[-1].isoformat()}' "
                          f"not within time_coverage_end attribute: '{timeCoverEnd.isoformat()}'")
                    rc = 1
                if (timeBounds is not None and timeBoundsKey != 'climatology'
                        and timeCoverEnd != timeBounds[-1, 1]):
                    print(f"{'':<8}{RC_ERR} time bound [-1,1] is not matching "
                          f"time_coverage_end attribute: '{timeCoverEnd.isoformat()}'")
                    rc = 1

        return rc, timeCoverStart, timeCoverEnd


    def _checkCoordinatesTime(self, timeC, expClimate=False, expRecords=None, expResolution=None, recordStatus=None):
        """
        Check time coordinates of a netcdf file.
        """

        rc = 0
        ds = self.Dataset

        if recordStatus is None:
            recordStatus = {}

        timeDuration = (decode_timeDuration(ds.time_coverage_duration)
                        if hasattr(ds, "time_coverage_duration") else None)

        # validate axis, climatology, and units attributes
        axis_rc, tUnits, sinceYr, calendar = self._checkTimeAxis(timeC, expClimate)
        rc = axis_rc

        # decode time steps (only when units are valid)
        tSteps = None
        if axis_rc <= 0:
            tSteps, step_rc = self._decodeTimeSteps(timeC, sinceYr, calendar)
            if step_rc:
                rc = 1

        # check filename-derived time and record count
        decode = re.match(CMSAF_NAMING_STANDARD, self.File)
        try:
            timeStepFn = datetime.datetime.strptime(self.File[5:17], "%Y%m%d%H%M")
            timeStepFn = timeStepFn.replace(tzinfo=pytz.utc)
        except ValueError:
            timeStepFn = None

        if tSteps is not None:
            if decode is not None and not ds.isSwathData():
                if timeStepFn is not None and decode.group(9) == "002" and decode.group(1) == "UTH":
                    timeStepFn = timeStepFn - datetime.timedelta(minutes=30)
                if timeStepFn is not None and timeStepFn != tSteps[0]:
                    print(f"{'':<8}{RC_ERR} time record mismatch, expecting {timeStepFn.isoformat()} as first record")
                    rc = 1
            else:
                if timeStepFn is not None and tSteps[0] is not None and timeStepFn > tSteps[0]:
                    print(f"{'':<8}{RC_ERR} time record mismatch, expecting {timeStepFn.isoformat()} before first record")
                    rc = 1

            if expRecords is not None and expRecords != tSteps.size:
                print(f"{'':<8}{RC_ERR} Expecting {expRecords} records but found {tSteps.size}.")
                rc = 1

        # decode and validate time bounds (only when units are valid)
        timeBoundsKey, timeBounds, bounds_rc = 'bounds', None, 0
        if axis_rc <= 0:
            timeBoundsKey, timeBounds, bounds_rc = self._decodeTimeBounds(
                timeC, tUnits, sinceYr, calendar)
        else:
            print(f"{'':<8}{RC_INFO} skipping time bounds check: time units are invalid")
        if bounds_rc:
            rc = 1

        # per-record checks and time coverage attributes
        timeCoverStart = None
        timeCoverEnd   = None
        if tSteps is not None:
            print("")
            rec_rc = self._checkTimeRecords(tSteps, timeBounds, timeBoundsKey,
                                            timeDuration, expResolution, recordStatus, timeC)
            if rec_rc:
                rc = 1

            cov_rc, timeCoverStart, timeCoverEnd = self._checkTimeCoverage(
                tSteps, timeBounds, timeBoundsKey)
            if cov_rc:
                rc = 1

            print(f"\n{'':<8}first time record: {tSteps[0].isoformat()}")
            print(f"{'':<8}last  time record: {tSteps[-1].isoformat()}")

        if timeCoverStart and timeCoverEnd:
            print(f"{'':<8}time coverage: [{timeCoverStart.isoformat()}, {timeCoverEnd.isoformat()}]")

        return rc


    def _checkCoordinatesGeo(self, coordVar, axisTime, shortName, longName, expAxis=None):
        """
        Check geo coordinates of a netcdf file.
        """

        rc = 0
        ds = self.Dataset
        resFile = cmsaf_decode_grid(self.File)

        # test axis attribute
        if not hasattr(coordVar, "axis"):
            # exclude from checks if not fixed
            if axisTime.name in coordVar.dimensions:
                print(f"{'':<4}{RC_INFO} coordinate is not fixed in time")
                rc = 10
            else:
                print(f"{'':<8}{RC_ERR} missing mandatory attribute 'axis'")
                rc = -1
        elif expAxis is not None:
            if coordVar.axis != expAxis:
                print(f"{'':<8}{RC_ERR} invalid value attribute 'axis={coordVar.axis}'")
                rc = -1

        # test axis values
        if rc <= 0:
            coord = coordVar[:]
            coordOrder = 1
            if len(coord.shape) == 1:
                if coord[0] > coord[-1]:
                    coordOrder = -1
                    coord[:] = coord[::-1]
            precision = np.finfo(coord.dtype).precision
            np.set_printoptions(precision=precision, suppress=False, threshold=10,linewidth=80)

            coordMin = np.min(coord)
            coordMax = np.max(coord)

            # test global geospatial bounds
            geoMinAttrName = f"geospatial_{shortName}_min"
            geoMinAttr = None
            if hasattr(ds,geoMinAttrName):
                try:
                    tmp_ = getattr(ds,geoMinAttrName)
                    if tmp_ > coordMin:
                        print(f"{'':<8}{RC_ERR} {geoMinAttrName} mismatch: {tmp_} > {coordMin}")
                        rc = 1
                except TypeError:
                    print(f"{'':<8}{RC_ERR} {geoMinAttrName}: unexpected data format")
                    rc = 1
                else:
                    geoMinAttr = tmp_

            geoMaxAttrName = f"geospatial_{shortName}_max"
            geoMaxAttr = None
            if hasattr(ds,geoMaxAttrName):
                try:
                    tmp_ = getattr(ds,geoMaxAttrName)
                    if tmp_ < coordMax:
                        print(f"{'':<8}{RC_ERR} {geoMaxAttrName} mismatch: {tmp_} < {coordMax}")
                        rc = 1
                except TypeError:
                    print(f"{'':<8}{RC_ERR} {geoMaxAttrName}: unexpected data format")
                    rc = 1
                else:
                    geoMaxAttr = tmp_

            # get grid definition
            resAttr = None
            geoResAttrName = f"geospatial_{shortName}_resolution"
            if hasattr(ds,geoResAttrName):
                tmp_ = getattr(ds,geoResAttrName)
                if type(tmp_) == type(str('s')):
                    resAttr = re.match('^([0-9.]+) +(degree).*$', tmp_)
                    if resAttr is not None:
                        resAttr = np.float64(resAttr.group(1))
                else:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} {geoResAttrName} :: must be a text type")

            # select grid
            if (resAttr is None) and (resFile is None):
                coordRes = None
            elif resAttr is not None:
                coordRes = resAttr
            elif resFile is not None:
                coordRes = resFile
            else:
                coordRes = resAttr
                if resFile != resAttr:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} grid definition from file name '{resFile}' <--> and attributes '{resAttr}'")

            # check grid
            if coordRes is not None:
                if np.ma.count_masked(coord) > 0:
                    print(f"{'':<8}{RC_ERR} {longName} contains missing data")
                    rc = 1
                else:
                    finfo    = np.finfo(coord.dtype)
                    eps      = significant_digits(coord) * finfo.resolution
                    indx     = np.arange(np.int64(0),len(coord),dtype=np.int64)
                    meshExp  = np.multiply(indx, np.multiply(np.int64(10000),coordRes)).astype(np.int64)
                    tmp      = np.rint(1./eps[0]).astype(np.float64)
                    tmp      = np.divide(np.rint(np.divide(coord[0], eps[0])), tmp)
                    tmp      = np.multiply(tmp, np.int64(10000)).astype(np.int64)
                    meshExp  = np.multiply(np.add(tmp,meshExp),np.float64(0.0001))
                    tmp      = np.absolute(np.subtract(meshExp, coord))
                    eps      = float_spacing(meshExp, coord)
                    indx     = np.where(tmp > eps)[0]
                    if len(indx) > 0:
                        print(f"{'':<8}{RC_ERR} at {len(indx)} locations:")
                        print(f"{'':<10}found:    {coord[indx]}")
                        print(f"{'':<10}expecting:{meshExp[indx]}")
                        if not self.lazy:
                            rc = 1
                        else:
                            print(f"{'':<8}{RC_INFO} Ignoring while beeing lazy")

                    # test 0,0 [must not be in center]
                    tmp  = np.absolute(coord)
                    indx = np.where(tmp < eps)[0]
                    if len(indx) > 0:
                        rc = 1
                        print(f"{'':<8}{RC_ERR} {shortName}=0 is not allowed as {longName} center value.")

            # test coordinate bounds
            boundsKey = 'bounds'
            boundsVar = None
            bounds    = None
            if hasattr(coordVar,boundsKey):
                tmp = getattr(coordVar,boundsKey)
                if coordVar.group().name == "/":
                    if tmp in ds.variables:
                        boundsVar = ds.getvar(tmp)
                else:
                    if tmp in ds.getgrp(coordVar.group().path).variables:
                        boundsVar = ds.getgrp(coordVar.group().path).variables[tmp]
                if boundsVar is not None:
                    if coordOrder == 1:
                        bounds = boundsVar[:]
                    else:
                        bounds = np.flip(boundsVar[:])
                else:
                    print(f"{'':<8}{RC_ERR} Missing configured bounds variable '{tmp}'.")
                    rc = 1

            if bounds is None:
                rc = 1
                print(f"{'':<8}{RC_ERR} missing bounds for {longName} coordinate")
            else:
                if shortName == 'lat':
                    leftName = "lower"
                    rightName = "upper"
                elif shortName == 'lon':
                    leftName = "left"
                    rightName = "right"
                leftMaxName = leftName+"most"
                rightMaxName = rightName+"most"

                print(f"{'':<8}[{bounds[0,0]!s} -> {bounds[-1,0]!s}] {leftName} bounds")
                print(f"{'':<8}[{bounds[0,1]!s} -> {bounds[-1,1]!s}] {rightName} bounds")

                # test coordinate within bounds
                indx1 = np.where(bounds[:,0] > coord)[0]
                indx2 = np.where(bounds[:,1] < coord)[0]
                if len(indx1) > 0 or len(indx2) > 0:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} {longName} values not within bounds")

                # test for gaps in coordinate bounds
                tmp = np.subtract(boundsVar[0:-2,1],boundsVar[1:-1,0])
                indx = np.where(tmp > 0)[0]
                if len(indx) > 0:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} gaps in {longName} bounds")

                # test for ovarlap in coordinate bounds
                indx = np.where(tmp < 0)[0]
                if len(indx) > 0:
                    rc = 1
                    print(f"{'':<8}{RC_ERR} {longName} bounds overlap")

                # test bounds against global attribute
                if geoMinAttr is not None:
                    if geoMinAttr != bounds[0,0]:
                        rc = 1
                        print(f"{'':<8}{RC_ERR} mismatch between {leftMaxName} {longName} bound '{boundsVar[0,0]}' and {geoMinAttrName} '{geoMinAttr}'")
                if geoMaxAttr is not None:
                    if geoMaxAttr != bounds[-1,1]:
                        rc = 1
                        print(f"{'':<8}{RC_ERR} mismatch between {rightMaxName} {longName} bound '{boundsVar[-1,1]}' and {geoMaxAttrName} '{geoMaxAttr}'")

            # print result
            if coordRes is not None:
                print(f"{'':<8}[{coordMin!s} -> {coordMax!s} by {coordRes!s}]")
            else:
                print(f"{'':<8}[{coordMin!s} -> {coordMax!s}]")

        rc = 0 if rc>=10 else rc
        return(rc)


    def _checkCompression(self):
        """
        Check compression settings for all 3-dimensional variables.
        """

        rc = 0
        ds = self.Dataset

        # find time coordinate
        axisTime = ds.getCoordinates("time", shortName=["time"])

        if (ds.data_model == "NETCDF4") or (ds.data_model == "NETCDF4_CLASSIC"):
            vList = ds.getVariableList()
            vFilter = []

            # filter variables to test
            if ds.isSwathData():
                for vName in vList:
                    aTime = None
                    for dim in ds.getvar(vName).get_dims():
                        aTime = ds.matchCoordinate(dim, axisTime)
                        if aTime is not None: break
                    if aTime is not None: vFilter.append(vName)
            else:
                for vName in vList:
                    if len(ds.getvar(vName).shape) >= 3: vFilter.append(vName)

            for vName in vFilter:
                xRc = 1
                filters = ds.getvar(vName).filters()
                if filters is None:
                    xRc = 1
                else:
                    if 'zlib' in filters:
                        if filters['zlib']:
                            print(f"{vName:<15} level={filters['complevel']}")
                            if filters['complevel'] > 0: xRc = 0
                if xRc == 1:
                    rc = 1
                    print(f"{RC_ERR} Variable {vName} is not compressed.")

        else:
            rc = 1
            print(f"{RC_ERR} file data dype is not netcdf4: '{ds.data_model}'")

        return rc


    def _checkVariables(self):
        """
        Check variables.
        """

        rc = 0
        ds = self.Dataset
        vList = ds.getVariableList()

        # recommended attributes
        listRec = ['units', 'standard_name', 'grid_mapping']
        listRecAxis = ['grid_mapping']
        listRecSwath = ['grid_mapping']

        # mandatory attributes
        listMan = ['long_name']

        # items to ignore
        listManSkip = {'cell_index': ['long_name']}
        listRecSkip = {'cell_index': ['grid_mapping']}

        # skip variables
        listSkip = ['record_status']

        # find coordinates
        axisLon  = ds.getCoordinates("longitude", shortName=["lon","longitude"])
        axisLat  = ds.getCoordinates("latitude", shortName=["lat","latitude"])
        axisTime = ds.getCoordinates("time", shortName=["time"])

        # loop variables and test attributes
        for vName in vList:
            vNameB = os.path.basename(vName)
            if vNameB in listSkip:
                continue

            print(vName)
            var = ds.getvar(vName)
            if var.group().name != "/":
                grp = ds.groups[var.group().name]
            else:
                grp = ds.ds

            # test grid mapping variable
            if hasattr(var, 'grid_mapping_name'):

                if var.grid_mapping_name == "latitude_longitude":
                    if len(axisLat) == 0:
                        print(f"{RC_ERR} missing required latitude coordinate")
                        rc = 1
                    if len(axisLon) == 0:
                        print(f"{RC_ERR} missing required longitude coordinate")
                        rc = 1

                varGridMapping = grp.get_variables_by_attributes(grid_mapping_name=var.grid_mapping_name)
                if len(varGridMapping) == 0:
                    varGridMapping = ds.get_variables_by_attributes(grid_mapping_name=var.grid_mapping_name)
                if len(varGridMapping) > 1:
                    print(f"{RC_ERR} grid_mapping_name='{var.grid_mapping_name}' ambiguous")
                    rc = 1

                continue

            # find dimensions
            aLon  = None
            aLat  = None
            aTime = None
            for dim in var.get_dims():
                if aLon is None:
                    aLon = ds.matchCoordinate(dim, axisLon)
                if aLat is None:
                    aLat = ds.matchCoordinate(dim, axisLat)
                if aTime is None:
                    aTime = ds.matchCoordinate(dim, axisTime)

            # test dimensions, all variables with lat and lon must also have a time coordinate
            if aLon is not None and aLat is not None:
                if aTime is None:
                    rc = 1
                    print(f"{'':<4}{RC_ERR} {vName} :: missing time dimension")

            # test attributes
            for item in [*listMan, *listRec]:
                itRc = 0
                if item in listRecAxis and vName in axisTime:
                    itRc = 0
                elif item in listRecAxis and vName in axisLon:
                    itRc = 0
                elif item in listRecAxis and vName in axisLat:
                    itRc = 0
                elif item in listRecSwath and ds.isSwathData():
                    itRc = 0
                elif not hasattr(var, item):
                    itRc = 1
                    for gvName in grp.variables:
                        gv = grp.variables[gvName]
                        if hasattr(gv, 'bounds') and gv.bounds == vNameB:
                            itRc = 0
                            break
                        if hasattr(gv, 'climatology') and gv.climatology == vNameB:
                            itRc = 0
                            break
                if itRc == 1:
                    xList = listManSkip[vNameB] if vNameB in listManSkip else []
                    if item in listMan and item not in xList:
                        rc = 1
                        print(f"{'':<4}{RC_ERR} {vName} :: missing mandatory attribute '{item}'")
                    xList = listRecSkip[vNameB] if vNameB in listRecSkip else []
                    if item in listRec and item not in xList:
                        print(f"{'':<4}{RC_WARN} {vName} :: missing recommended attribute '{item}'")

            # test flags
            if hasattr(var,'flag_values') and hasattr(var,'flag_meanings'):
                flagV = np.atleast_1d(var.flag_values)
                flagM = np.atleast_1d(var.flag_meanings.split(" "))
                if len(flagV) != len(flagM):
                    rc = 1
                    print(f"{'':<4}{RC_ERR} {vName} :: mismatch between flag_values and flag_value")

            # test grid mapping
            if hasattr(var, 'grid_mapping'):
                varGrid = ds.getVariableByName(var.grid_mapping)
                if not ((os.path.join(grp.path,var.grid_mapping) in varGrid.keys()) or
                       (os.path.join("/",var.grid_mapping) in varGrid.keys())):
                    print(f"{RC_ERR} missing defined 'grid_mapping' variable '{var.grid_mapping}'")
                    rc = 1

        # test variables defined in variable_id
        if hasattr(ds, 'variable_id'):
            print(f"\nvariable_id :: '{ds.variable_id}'")
            for item in ds.variable_id.split(","):
                try:
                    var = ds.getvar(item)
                except (KeyError, IndexError):
                    rc = 1
                    print(f"{'':<4}{RC_ERR} missing variable '{item}'")
                else:
                    print(f"{'':<4}{RC_OK} '{item}'")

        return rc


def main():
    """
    Main program
    """

    import argparse
    import fnmatch
    from os  import getenv
    from sys import argv,exit

    # get standards search paths in order: (c) -s option, (b) prefix, (a) script share
    prefix = getenv('CMSAF_CHECKER_PREFIX')
    path_a = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "share")
    path_b = prefix if prefix is not None else None
    default_search_paths = [p for p in [path_b, path_a] if p is not None]

    print(f"CMSAF Checker Version {__version__}")

    # argument parser
    parser = argparse.ArgumentParser(prog='cmsaf-checker')
    parser.add_argument('-s', '--cmsaf_metadata_standard',
        default=None, metavar='PATH',
        help='additional search path for CM SAF standard and GCMD keyword files. '
             'Files are resolved in order: this path, $CMSAF_CHECKER_PREFIX, '
             f'{path_a}')
    parser.add_argument('-v', '--version',
        help='CM SAF standards version to apply')
    parser.add_argument('-f', '--standard_file',
        help='Test files using this as the reference')
    parser.add_argument('-r', '--reference',
        help='Test files using this as the reference')
    parser.add_argument('-i', '--ignore_attr',
        help='Ignore list of attributes')
    parser.add_argument('-c', '--coordinates', action='store_true',
        help='Test coordinates time, latitude, longitude')
    parser.add_argument('-m', '--missing', nargs='?', const='filename',
        help='Test for missing files')
    parser.add_argument('-l', '--lazy', action='store_true',
        help='Turn some errors to warnings')
    parser.add_argument('-d', '--directory',
        help='Search for files with pattern in this directory.')
    parser.add_argument('files', nargs='+')

    args = parser.parse_args()

    # build final search path list: prepend explicit -s if given
    search_paths = default_search_paths
    if args.cmsaf_metadata_standard is not None:
        search_paths = [args.cmsaf_metadata_standard] + search_paths
    for i, p in enumerate(search_paths, start=ord('a')):
        print(f"CMSAF Standards path ({chr(i)}) '{p}'")

    # get a new checker object
    if args.reference == None:
        inst = CMSAFChecker(search_paths=search_paths, version=args.version,
            coordinates=args.coordinates, lazy=args.lazy, ignore=args.ignore_attr,
            standard_file=args.standard_file)
    else:
        inst = CMSAFChecker(referenceFile=args.reference, ignore=args.ignore_attr,
            coordinates=args.coordinates, lazy=args.lazy,
            standard_file=args.standard_file)

    # file pattern expansion
    files = []
    if args.directory is not None:
        for root, dirs, names in os.walk(args.directory, followlinks=True):
            tmp = fnmatch.filter(names, args.files[0])
            for fn in tmp:
                if args.directory != root:
                    continue
                fnp=os.path.join(root,fn)
                if os.access(fnp, os.R_OK):
                    files.append(fnp)
                else:
                    print(f"Skipping file '{fnp}'")
    else:
        files.extend(args.files)
    files.sort(key=lambda s: os.path.basename(s))

    # result dict
    res = {'OK': 0, 'FAILED': 0, 'MISSING' : 0}

    # decode first file name
    fileDelta = None
    if len(files) > 1 and args.missing is not None:
        if args.missing == "filename":
            fileName = os.path.realpath(files[0])
            fileName = os.path.basename (fileName)
            fileAttr = re.match(CMSAF_NAMING_STANDARD, fileName)
            fileStep = fileAttr.group(2)
        else:
            fileStep = args.missing
        if fileStep == 'd':
            fileDelta = datetime.timedelta(days=1, hours=0, minutes=0, seconds=0)
        elif fileStep == 'm':
            fileDelta = datetime.timedelta(days=31, hours=0, minutes=0, seconds=0)
        elif fileStep == 'h':
            fileDelta = datetime.timedelta(days=0, hours=1, minutes=0, seconds=0)
        elif fileStep == 'M15':
            fileDelta = datetime.timedelta(days=0, hours=0, minutes=15, seconds=0)

    # loop files
    lastTime = None
    for index, file in enumerate(files):
        # test missing
        if fileDelta is not None:
            currentFile = os.path.realpath(file)
            currentFile = os.path.basename (currentFile)
            currentTime = datetime.datetime.strptime(currentFile[5:17], "%Y%m%d%H%M")
            currentTime = currentTime.replace(tzinfo=pytz.utc)
            if lastTime is not None:
                while 1:
                    nextTime = lastTime + fileDelta
                    if fileStep == 'm':
                        nextTime = nextTime.replace(day=1)
                    fileDiff = currentTime-nextTime
                    if fileDiff.total_seconds() > 0:
                        print(f"\n{'':=^80}\nMissing File for {nextTime.isoformat('T')}\n{'':=^80}")
                        res['MISSING'] += 1
                        lastTime = nextTime
                    else:
                        break
            lastTime = currentTime

        # print info
        print(f"\n{'':=^80}\nChecking File {index+1}/{len(files)}\n{'':=^80}\n'{file}'")

        # check current file
        inst._reset()
        rc = inst.checker(file)
        if rc == 0:
            res['OK'] += 1
            rcMsg = RC_OK
        else:
            res['FAILED'] += 1
            rcMsg = RC_FAIL
        print(f"\n{'':-^80}\n{rcMsg} <<< result for {file}\n{'':-^80}")

    # close reference file
    if inst.refDataset is not None:
        inst.refDataset.close();
        inst.refDataset = None

    # final result
    print(f"\n{'':=^80}\nOverall Summary\n{'':=^80}")
    print(f"Out of {len(files)}, {res['FAILED']} FAILED")
    if (args.missing) and (fileDelta is not None):
        print(f"{res['MISSING']} files MISSING")

    if res['FAILED'] == 0:
        exit(0)
    else:
        exit(1)

if __name__ == '__main__':
    main()

# ex: set foldmethod=indent:
