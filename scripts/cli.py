#!/usr/bin/env python3
from netCDF4 import Dataset, num2date
from astropy.time import Time
import re, string, types, csv, os.path
import glob
import numpy as np
import pytz, datetime

__version__ = "3.0.0"
__prefix__  = ""
STANDARD = ''

CMSAF_NAMING_STANDARD = r'^([A-Z]{3})([dhimpswa])([ncfhmds])((?:19|20)\d\d)(0[1-9]|1[012])(0[1-9]|[12]\d|3[01])([01]\d|2[0-3])([0-5]\d)(\d{3})(\d{2})([0-9A-Z]{5})([A-Z0-9]{2})([A-Z0-9]{2})((?:\.hdf|\.hdf\.gz|\.gz|\.nc|\.nc\.gz){0,1})$'

RC_ERR  = "## ERROR ##"
RC_WARN = "## WARNING ##"
RC_OK   = "## OK ##"
RC_FAIL = "## FAILED ##"
RC_INFO = "## INFORMATION ##"

from xml.sax import ContentHandler
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces


def normalize_whitespace(text):
    "Remove redundant whitespace from a string."
    return ' '.join(text.split())


def significant_digits(value, digits):
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


def leap_year(year):
    """
    Check if year is a leap year
    """

    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
                return True
            else:
                return False
        else:
             return True
    else:
        return False


def cmsaf_decode_grid(filename):
    """
    Decode the grid definition from CM SAF standarda file name
    """
    decode = re.match(CMSAF_NAMING_STANDARD, filename)
    if decode == None:
        return None

    if (decode.group(10)) == "19":
        return np.float64("0.25")
    elif (decode.group(10)) == "23":
        return np.float64("0.05")
    elif (decode.group(10)) == "26":
        return np.float64("0.1")
    elif (decode.group(10)) == "20":
        return np.float64("1.0")

    return None


def decode_timeDuration(duration):
    """
    Decode time duration strings
    """
    regDuration1 = "^P(?:(?P<year>[0-9]+)Y)?(?:(?P<month>[0-9]+)M)?(?:(?P<week>[0-9]+)W)?(?:(?P<day>[0-9]+)D)?(?:T(?:(?P<hour>[0-9]+)?H)?(?:(?P<minute>[0-9]+)?M)?(?:(?P<second>[0-9]+)?S)?)?$"
    regDuration2 = "^P(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})T(?P<hour>[0-9]{2}):(?P<minute>[0-9]{2}):(?P<second>[0-9]+)$"

    timeDuration = re.match(regDuration1, duration)
    if (timeDuration is None):
        timeDuration = re.match(regDuration2, duration)

    if (timeDuration is None):
      return None

    result = np.zeros([7], dtype=np.float64)
    if (timeDuration.group('year') is not None):
        result[0] = np.float64(timeDuration.group('year'))
    if (timeDuration.group('month') is not None):
        result[1] = np.float64(timeDuration.group('month'))
    if "week" in timeDuration.groupdict():
        if (timeDuration.group('week') is not None):
            result[2] = np.float64(timeDuration.group('week'))
    if (timeDuration.group('day') is not None):
        result[3] = np.float64(timeDuration.group('day'))
    if (timeDuration.group('hour') is not None):
        result[4] = np.float64(timeDuration.group('hour'))
    if (timeDuration.group('minute') is not None):
        result[5] = np.float64(timeDuration.group('minute'))
    if (timeDuration.group('second') is not None):
        result[6] = np.float64(timeDuration.group('second'))

    return result


class CMSAFStandard(ContentHandler):
    """
    XML Parser class for CM SAF attributes

    Parse the xml standard table, reading all entries into a dictionary
    """

    def __init__(self):
        self.inContent = 0
        self.inRegex = 0
        self.inKeywords = 0
        self.inVersionNoContent = 0
        self.inLastModifiedContent = 0
        self.inInclude = 0
        self.dict = {}
        self.entry = {}
        self.content = ""
        self.regex = {}
        self.keywords = ""

    def startElement(self, name, attrs):
      # attributes
      if name == 'entry':
        self.entry = {
          'type':     normalize_whitespace(attrs.get('type', "s")),
          'required': normalize_whitespace(attrs.get('required', "no")),
          'evaluate': normalize_whitespace(attrs.get('evaluate', "no")),
          'list' :    normalize_whitespace(attrs.get('list', "")),
          'join':     attrs.get('join', "or").strip(' '),
          'content':  [],
          'regex':    [],
          'keywords':  []
        }
        self.this_id = normalize_whitespace(attrs.get('id', ""))

      # possible content
      elif name == 'content':
        self.inContent = 1
        self.content = {}
        self.content['value'] = ""
        self.content['type']  = attrs.get('type', "string").strip(' ');

      # regular expression
      elif name == 'regex':
        self.inRegex = 1
        self.regex = {}
        self.regex['value'] = ""
        self.regex['warn']  = attrs.get('warn', "").strip(' ')
        self.regex['type']  = attrs.get('type', "string").strip(' ');

      # keyword list
      elif name == 'keywords':
        self.inKeywords = 1
        self.keywords = ""

      elif name == 'version_number':
        self.inVersionNoContent = 1
        self.version_number = ""

      elif name == 'last_modified':
        self.inLastModifiedContent = 1
        self.last_modified = ""

      elif name == 'include':
        self.inInclude = 1
        self.include = ""

    def characters(self, ch):

        if self.inContent:
          self.content['value'] = self.content ['value']+ ch

        elif self.inRegex:
          self.regex['value'] = self.regex['value'] + ch

        elif self.inKeywords:
          self.keywords = self.keywords + ch

        elif self.inVersionNoContent:
            self.version_number = self.version_number + ch

        elif self.inLastModifiedContent:
            self.last_modified = self.last_modified + ch

        elif self.inInclude:
            self.include = self.include + ch

    def endElement(self, name):
        if name == 'entry':
            self.dict[self.this_id] = self.entry.copy()
            self.entry.clear()

        elif name == 'content':
            self.inContent = 0
            self.entry['content'].append(self.content)

        elif name == 'regex':
            self.inRegex = 0
            self.entry['regex'].append(self.regex)

        elif name == 'keywords':
          self.inKeywords = 0
          self.entry['keywords'].append(self.keywords)

        # If it's the end of the version_number element, save it
        elif name == 'version_number':
            self.inVersionNoContent = 0
            self.version_number = normalize_whitespace(self.version_number)

        # If it's the end of the last_modified element, save the last modified date
        elif name == 'last_modified':
            self.inLastModifiedContent = 0
            self.last_modified = normalize_whitespace(self.last_modified)

        # If it's the end of the include element, save it
        elif name == 'include':
            self.inInclude = 0
            self.include = normalize_whitespace(self.include)


class Keywords:
    """
    class to read GCMD keywords
    """

    def __init__(self, filename, standardsPath=None):
        self.filename = filename
        self.groups = None
        self.keywordList = {}

        if standardsPath:
          self.filename = standardsPath + "/" + self.filename

        # try to open data file
        try:
            fh = open(self.filename, "r")
            fh.close()

        # fallback to global STANDARD
        except IOError as detail:
            self.filename = STANDARD+"/"+filename

    def readFile(self):
        try:
            fh = open(self.filename, "r")
        except IOError as detail:
            print (detail)
            return 2

        # cvs reader
        reader = csv.reader(fh, delimiter=',', quotechar='"')

        # loop lines
        for row in reader:
            if 'Keyword Version' and 'Revision' in row:
                meta = re.match(r'.*Keyword Version: (\d\.\d).*', row[0])
                self.keywordList['Version'] = meta.group(1)
            elif (not self.groups) and ('UUID' in row):
                self.groups = row
                self.groups.pop()
            elif self.groups:
                uuid = row.pop().strip('"')
                if (len(row) == len(self.groups)):
                    kw = {}
                    for index, item in enumerate(self.groups):
                        kw[item] = row[index].strip()
                        self.keywordList[uuid] = kw
        fh.close()

        return 0

    def findKeyword(self,keyword="Version"):
        result = []

        keyword = keyword.upper()
        for key, kw in self.keywordList.items():
            if key == keyword:
                result.append(kw)
            if 'Short_Name' in kw:
                if kw['Long_Name'].upper() == keyword:
                    result.append(kw)
                elif kw['Short_Name'].upper() == keyword:
                    result.append(kw)
                elif 'Sub_Category' in kw:
                    if kw['Sub_Category'].upper() == keyword:
                        result.append(kw)
            if 'Series_Entity' in kw:
                if kw['Series_Entity'].upper() == keyword:
                    result.append(kw)
            if 'Variable_Level_1' in kw:
                if kw['Variable_Level_1'].upper() == keyword:
                    result.append(kw)
            if 'Variable_Level_2' in kw:
                if kw['Variable_Level_2'].upper() == keyword:
                    result.append(kw)
            if 'Variable_Level_3' in kw:
                if kw['Variable_Level_3'].upper() == keyword:
                    result.append(kw)
            if 'Term' in kw and kw['Variable_Level_1'] == "":
                if kw['Term'].upper() == keyword:
                    result.append(kw)

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
                if (len(entry[item]) > 0): list.append(entry[item])
            kwList.append(" > ".join(list))

        # return as list of strings
        return kwList


class DatasetX(Dataset):
    """
    Expand standard python Dataset netcdf4 class
    """

    def __init__(self, *args, **kwargs):
       Dataset.__init__(self, *args, **kwargs)


    def __getattribute__(self, item):
        if item == "__dict__":
            return Dataset.__getattr__(self, item)
        else:
            return Dataset.__getattribute__(self, item)


    def isSwathData(self):
        """
        """
        rc = False

        if hasattr(self,"cdm_data_type"):
            if (self.cdm_data_type == "swath"):
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
            if (os.path.basename(item) == name):
                result[item] = self[item]

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
        Find coordinate variables
        """
        axis = self.getVariableByStandardName(standardName)
        if len(axis) == 0:
            for vName in shortName:
                axis = self.getVariableByName(vName)
                for item in axis.keys():
                    print (f"## WARNING ##: missing standard name attribute for '{item}'")

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
                if (_path == iPath):
                    keyTime = _keyTime
                    break
            if (keyTime is not None or iPath == "/"):
                break
            else:
                iPath = os.path.dirname(iPath)

        return(keyTime)


class CMSAFChecker:
    """
    CM SAF Checking class
    """

    def __init__(self, metadataStandards=None, version=None, referenceFile=None,
        coordinates=False, ignore=None, lazy=False):

        self.metadataStandards = metadataStandards
        self.Dataset           = None
        self.testFn            = None
        self.refFile           = referenceFile
        self.refDataset        = None
        self.coordinates       = coordinates
        self.lazy              = lazy
        self.gIgnoreAtt        = []
        self.vIgnoreAtt        = []
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

        # open file
        if self.refFile is not None:
            try:
                print(f"Reference File: '{self.refFile}'\n")
                self.refDataset = DatasetX(self.refFile, mode='r')
            except:
                print("\nCould not open reference file, please check that NetCDF is formatted correctly.\n")
                raise

        # Set up dictionary of defined standards
        if (self.refDataset == None):
            parser = make_parser()
            parser.setFeature(feature_namespaces, 0)
            self.std_name_dh = CMSAFStandard()
            parser.setContentHandler(self.std_name_dh)
            fn=self.metadataStandards
            if (os.path.isdir(fn)):
                fn = fn+"/cmsaf_metadata_standard"+self.version+".xml"
            else:
                self.metadataStandards = STANDARD
            try:
                parser.parse(fn)
            except IOError as detail:
                print (detail)
                raise

            # read included file
            if (hasattr(self.std_name_dh, 'include')):
                print (f"Including {self.std_name_dh.include}.")
                include_ = CMSAFStandard()
                parser.setContentHandler(include_)
                try:
                    parser.parse(self.std_name_dh.include)
                except IOError as detail:
                    print (detail)
                    raise

                # update dict
                tmp = self.std_name_dh.dict
                self.std_name_dh.dict = include_.dict
                self.std_name_dh.dict.update(tmp)

        # attributes to ignore
        if (self.refFile is not None):
            self.gIgnoreAtt = ["date_created", "time_coverage_start", "time_coverage_end", "filename", "date_modified", "history"]
        if (ignore is not None):
            attList = ignore.split(",")
            for att in attList:
                if (att.find('@') >= 0):
                    self.vIgnoreAtt.append(att)
                else:
                    self.gIgnoreAtt.append(att)


    def __del__(self):
        if (self.refDataset):
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
        except:
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

        return (rc+rcCompress+rcVariables)


    def _checkStandard(self):
        """
        check global metadata against CM SAF standard
        """
        rc = 1

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


    def _checkGlobalAttributes(self):
        """
        Check global attributes
        """

        print("\ncheck validity of global attributes.")
        ds = self.Dataset
        rc = 0

        # get file name
        fn = os.path.realpath(ds.filepath())
        fn = os.path.basename (fn)

        kwList = {}

        # loop and find required attributes
        for key in self.std_name_dh.dict:
            attr = self.std_name_dh.dict[key]

            # evaluate placeholders
            if attr['evaluate'] == "yes":

                # reference attribute
                if hasattr(ds, 'references'):
                    tmp = getattr(ds,'references')
                    for index, item in enumerate(attr['content']):
                        item['value'] = item['value'].replace("${references}", tmp)

                # year
                now = os.environ.get('CMSAF_RELEASE_YEAR')
                if now is None:
                    now = datetime.datetime.now().strftime("%Y")
                for index, item in enumerate(attr['content']):
                    item['value'] = item['value'].replace("${year}", now)

            # test if required is defined
            if attr['required'] == "yes":
                if not hasattr(ds, key):
                    if key in self.gIgnoreAtt:
                        print(f"{RC_INFO} Ignoring missing required attribute '{key}'")
                        self.info += 1
                        self.infoAttr.append(key)
                    else:
                        print(f"{RC_ERR} Missing required attribute '{key}'")
                        self.err += 1
                        self.errAttr.append(key)

            # test improper attributes
            if attr['required'] == "none":
                if hasattr(ds,key):
                    print(f"{RC_ERR} Found improper attribute '{key}'")
                    self.err += 1
                    self.errAttr.append(key)

        # loop global attributes and check
        for key in ds.ncattrs():

            attr  = getattr(ds,key)
            keyRc = 0

            # file name
            if key == 'filename':
                print(f"\n{key}:\n{attr}")
                if fn != ds.filename:
                    print(f"{RC_ERR} incorrect file name :: '{ds.filename}'")
                    self.err += 1
                    self.errAttr.append(key)

            if key in self.std_name_dh.dict:
                print(f"\n{key}:")
                std = self.std_name_dh.dict[key]

                # check attributes type
                attrType = type(attr)
                if (attrType == type(np.array([]))):
                    if (attr.size == 1) and (type(attr[0]) == type(np.float32(1.0))):
                        attrType = "f32"
                    elif (attr.size == 1) and (type(attr[0]) == type(np.float64(1.0))):
                        attrType = "f64"
                else:
                    if attrType == type(np.float64(1.0)):
                        attrType = "f64"
                    elif attrType == type(np.float32(1.0)):
                        attrType = "f32"
                    elif attrType == type(str('s')):
                        attrType = 's'

                if str.find(attrType, std['type']) == -1:
                    print(f"{RC_ERR} Incorrect attribute data type")
                    print(f"Expecting: {std['type']}, found: {attrType}")
                    keyRc = 1
                    self.err += 1
                    self.errAttr.append(key)

                # report empty attribute
                if attrType == 's':
                    if len(attr) == 0:
                        if (std['required'] == "yes"):
                            print(f"{RC_ERR} empty required attribute")
                            self.err += 1
                            self.errAttr.append(key)
                        else:
                            print(f"{RC_INFO} empty attribute")
                            self.info += 1
                            self.infoAttr.append(key)
                        continue

                # skip non string from further checks
                if attrType != 's':
                    print(attr)
                    continue

                # start with empty list
                attrList = []
                attrHits = []

                # empty hits lists
                for i in std['content']:
                    attrHits.append(0)
                for i in std['regex']:
                    attrHits.append(0)
                for i in std['keywords']:
                    attrHits.append(0)

                # make a list if attribute is defined as a list of values
                try:
                    if len(std['list']) > 0:
                        for row in csv.reader([attr], delimiter=std['list']):
                            attrList = row
                    else:
                        attrList = [attr]
                except UnicodeEncodeError as detail:
                    print(f"{RC_ERR} {detail}")
                    self.err += 1
                    self.errAttr.append(key)

                # loop content list
                for a in attrList:
                    attrMatch = []

                    # remove white spaces, and quotes
                    a_ = a.strip()
                    a_ = a_.strip('"')
                    if (a != a_):
                        print(f"{RC_WARN} white spaces or quotes detected")
                        self.warn += 1
                        if key not in self.warnAttr:
                            self.warnAttr.append(key)
                    a = a_
                    print(a)

                    # check attribute content
                    hitsIndex = 0
                    if len(std['content']) > 0:

                        # loop possible entries
                        validIndex = None
                        for index, item in enumerate(std['content']):
                            if item['value'] == a:
                                validIndex = hitsIndex+index
                                attrMatch.append(item)

                        # test if valid attribute
                        if not validIndex == None:
                            attrHits[validIndex] += 1
                        else:
                            print(f"{RC_ERR} incorrect attribute content :: '{a}'")
                            if len(std['content']) == 1:
                                print(f"Expecting: '{std['content'][0]['value']}'")
                            keyRc = 1
                            self.err += 1
                            if key not in self.errAttr:
                                self.errAttr.append(key)

                    # check attribute content with regular expression
                    hitsIndex += len(std['content'])
                    if len(std['regex']) > 0:
                        validIndex = None
                        for index, item in enumerate(std['regex']):
                          if re.search(item['value'], a):
                            if (item['warn'] != ""):
                                print(f"{RC_WARN} {item['warn']}")
                                self.warn += 1
                                if key not in self.warnAttr:
                                    self.warnAttr.append(key)
                            else:
                                validIndex = hitsIndex+index
                                attrMatch.append(item)

                        if not validIndex == None:
                            attrHits[validIndex] += 1
                        else:
                            print(f"{RC_ERR} incorrect attribute content :: '{a}'")
                            keyRc = 1
                            self.err += 1
                            if key not in self.errAttr:
                                self.errAttr.append(key)

                    # check keyword list
                    hitsIndex += len(std['regex'])
                    if len(std['keywords']) > 0:
                        # read new keyword list from file
                        keywordsFn = std['keywords'][0]
                        if not keywordsFn in kwList:
                            # evaluate keyword version number
                            if keywordsFn.find('${') >= 0:
                                vocabulary_version = None
                                vocabulary_name    = None
                                decode = re.match(r'^.*\$\{([a-z_]*)_version\}.*$', keywordsFn)
                                if decode is not None:
                                    vocabulary_name = decode.group(1)
                                if vocabulary_name is not None and hasattr(ds,vocabulary_name):
                                    decode = re.match(r'^.*Version +([0-9\.]*)$', getattr(ds,vocabulary_name))
                                    if (decode is not None):
                                        vocabulary_version = decode.group(1)
                                if (vocabulary_version is not None and vocabulary_name is not None):
                                    keywordsFn = keywordsFn.replace('${'+vocabulary_name+'_version}',vocabulary_version)
                            kw = Keywords(filename=keywordsFn, standardsPath=self.metadataStandards)
                            keyRc = kw.readFile()
                            if (keyRc != 0):
                                print(f"{RC_ERR} Test incomplete")
                                self.err += 1
                                self.errAttr.append(key)
                                continue
                            kwList[keywordsFn] = kw;
                        else:
                            kw = kwList[keywordsFn]

                        # loop matches
                        for mIndex, mItem in enumerate(attrMatch):
                            if mItem['type'] == "keyword":
                                # split keywords
                                entryList = re.split(" *> *", a)
                                # find elements using last entry in list
                                kwItem = kw.findKeywordList(entryList[-1])
                                if not kwItem:
                                    print(f"{RC_ERR} Keyword not found in List :: '{entryList[-1]}'")
                                    self.err += 1
                                    if key not in self.errAttr:
                                        self.errAttr.append(key)
                                else:
                                    # check all entries
                                    entryItem = " > ".join(entryList)
                                    entryP = ""
                                    if len(entryList) == 1:
                                        entryP = ".*" + re.escape(entryItem) + ".*"
                                    else:
                                        entryP = ".*" + re.escape(entryItem) + "$"
                                    entryHits = 0
                                    for item in kwItem:
                                        if re.search(entryP, item) is not None:
                                            entryHits += 1
                                            print(f"decoded as '{item}'")
                                    if entryHits != 1:
                                      self.err += 1
                                      print(f"{RC_ERR}keyword classification incorrect, expecting one of:\n##")
                                      print("\n## ".join(kwItem))
                                      if key not in self.errAttr:
                                          self.errAttr.append(key)
                                    else:
                                        attrHits[hitsIndex] += 1

                # evaluate hits
                if len(attrList) > 1:
                    if (std['join'].lower() == "or"):
                        total = 0
                        for i in attrHits:
                            total += i
                        if (total == 0) and (keyRc == 0):
                          print(f"{RC_ERR}missing a correct value for attribute '{key}'")
                          keyRc = 1
                          self.err += 1
                          if key not in self.errAttr:
                              self.errAttr.append(key)
                    elif (std['join'].lower() == "and"):
                        if attrHits.count(0) > 0:
                            hitsIndex = 0
                            for index, item in enumerate(std['content']):
                                if attrHits[hitsIndex+index] == 0:
                                    print(f"{RC_ERR}missing required specific attribute content :: '{item['value']}'")
                                    keyRc = 1
                                    self.err += 1
                                    if key not in self.errAttr:
                                        self.errAttr.append(key)

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
                print(f"{RC_ERR} Missing attribute :: '{attNameFull}'")
                rc = 1
            else:
                # mark as read
                attCheck[attName] = 1

                # file name
                if attName == 'filename':
                    if (self.File != new.filename):
                        print(f"{RC_ERR} incorrect file name :: '{new.filename}'")
                        rc = 1
                    continue

                # skip attributes that are allowed to change
                elif attName in ignore:
                    print(f"{RC_INFO} changing attribute {attNameFull} :: '{getattr(new,attName)}'")
                    continue

                    # test if attributes are identical
                    ar = getattr(ref, vAttName)
                    ac = getattr(new. vAttName)
                    if type(ar) == np.ndarray or type(ar) == np.float32 or type(ar) == np.float64:
                        if (np.isnan(ar).all() and np.isnan(ac).all()):
                            tmp = np.array([])
                        elif (np.isnan(ar).any() or np.isnan(ac).any()):
                            tmp = np.array([1])
                        else:
                            tmp = np.where(np.absolute(ar-ac) > 0)[0]
                        if (len(tmp) > 0):
                            print(f"{RC_ERR} attribute '{attNameFull}' differ")
                            print(f"{'':<4}{attNameFull} :: expecting '{ar}', found '{ac}'")
                            rc = 1
                    else:
                        if (ar != ac):
                            print(f"{RC_ERR} attribute '{attNameFull}' differ")
                            print(f"{'':<4}{attNameFull} ::  expecting '{ar}', found '{ac}'")
                            rc = 1

        # check for new attributes
        for attName in new.ncattrs():
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
                    if (item.find('@') >= 0):
                        itemList = item.split("@")
                        if (itemList[0] == varName or len(itemList[0]) == 0) :
                            ignore.append(itemList[1])

                # check data type
                if (vr.dtype != vc.dtype):
                    print(f"{RC_ERR} type of variable '{varName}' differs from reference.")
                    print(f"## ref='{vr.dtype}', file='{vc.dtype}'")
                    rc = 1

                # check array shape
                if (not self.refDataset.isSwathData() and vr.shape != vc.shape):
                    print(f"{RC_ERR} shape of variable '{varName}' differs from reference.")
                    print(f"## ref='{vr.shape}', file='{vc.shape}'")
                    rc = 1

                # track global attributes
                vAttCheck = {}

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
        if (rc == 0):
            print(f"\n{RC_OK} <<< reference attributes")
        else:
            print(f"\n{RC_FAIL} <<< reference attributes")
        rcAll = rc

        # process variables
        print("\nvariables")
        rc = self._checkReferenceVariables(self.Dataset, self.refDataset)

        # print result
        if (rc == 0):
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
        if decode is not None:
            # test for diurnal cycle
            if decode.group(3) == 'd':
                expRecords = 24
                if (timeResolution is None):
                    timeResolution = decode_timeDuration("PT1H")
                elif (not np.array_equal(timeResolution,decode_timeDuration("PT1H"))):
                    print(f"{RC_ERR} ## expecting 'PT1H' as time_coverage_resolution for diurnal cycle")
                    tests['time'] = 1
                    timeResolution = decode_timeDuration("PT1H")
            elif decode.group(2) != 'i':
                expRecords = 1

            # get expected number of records from time_coverage_resolution and time_coverage_duration
            if (decode.group(2) == 'i'):
                if (timeDuration is not None) and (timeResolution is not None):
                    td = datetime.timedelta(days=timeDuration[3],   hours=timeDuration[4],   minutes=timeDuration[5],   seconds=timeDuration[6])
                    tr = datetime.timedelta(days=timeResolution[3], hours=timeResolution[4], minutes=timeResolution[5], seconds=timeResolution[6])
                    expRecords = divmod(td.total_seconds(), tr.total_seconds())
                    expRecords = int(expRecords[0])

            # expect climate attribute
            if decode.group(2) == 'm' and decode.group(3) == 'd':
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
                    td = datetime.timedelta(days=timeDuration[3],   hours=timeDuration[4],   minutes=timeDuration[5],   seconds=timeDuration[6])
                    tr = datetime.timedelta(days=timeResolution[3], hours=timeResolution[4], minutes=timeResolution[5], seconds=timeResolution[6])
                    expRecords = divmod(td.total_seconds(), tr.total_seconds())
                    expRecords = int(expRecords[0])

        # find possible time axis
        axisTime = ds.getCoordinates("time", shortName=["time"])

        # find record_status variables
        tmp = ds.getVariableByName("record_status")
        recordStatus = {}
        print("\nrecord_status")

        if len(tmp) == 0:
            print(f"{RC_ERR} missing record_status variable")
            tests['record_status'] = 1
        else:
            tests['record_status'] = 0

        # test all record status variables
        for key in tmp.keys():
            print(f"{'':<4}{key}")
            item = tmp[key]
            recordStatus[key] = { "var": item, "dict": None, "val": None, "mask": None, "time": None, }

            itemPath = os.path.dirname(key)
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

        # test time coordinates
        tests['time'] = 0
        if len(axisTime) == 0:
            print(f"{RC_ERR} missing time variable")
            tests['time'] = 1
        else:
            print("\ntime")

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

        # test latitude
        tests['lat'] = 0
        resFile = cmsaf_decode_grid (self.File)
        axisLat = ds.getCoordinates("latitude", shortName=["lat","latitude"])

        if len(axisLat) == 0:
            print(f"{RC_ERR} missing latitude coordinate")
            tests['lat'] = 1
            rc = 1
        else:
            print("\nlatitude")

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

        if len(axisLon) == 0:
            print(f"{RC_ERR} missing longitude coordinate")
            tests['lon'] = 1
            rc = 1
        else:
            print("\nlongitude")

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


    def _checkCoordinatesTime(self, timeC, expClimate=False, expRecords=None, expResolution=None, recordStatus={}):
        """
        Check time coordinates of a netcdf file.
        """

        rc = 0
        ds = self.Dataset
        tUnits = None
        tSteps = None
        timeCoverStart = None
        timeCoverEnd = None
        calendar = "standard"

        # test axis attribute
        if not hasattr(timeC,"axis"):
            print(f"{'':<8}{RC_ERR} missing mandatory attribute 'axis'")
            rc = 1
        elif timeC.axis != "T":
            print(f"{'':<8}{RC_ERR} invalid value attribute 'axis={timeC.axis}'")
            rc = 1

        # test for climate bounds
        if expClimate:
            if not hasattr(timcC,'climatology'):
                print(f"{'':<8}{RC_WARN} Expecting attribute 'climatology' as attribute for time bounds time")

        # decode time record in file
        if not hasattr(timeC,'units'):
            print(f"{'':<8}{RC_ERR} missing mandatory attribute 'units' for time axis")
            rc = 1
        else:
            tUnits = timeC.units
            sinceYr = re.match('(.*since )((-)?[-0-9]{4})(.*)', tUnits)
            if sinceYr == None:
                print(f"{'':<8}{RC_ERR} invalid time axis: '{tUnits}'")
                rc = 1
            else:
                sinceYr = int(sinceYr.group(2))
                if sinceYr < 1958:
                    print(f"{'':<8}{RC_ERR} invalid time axis: '{tUnits}'")
                    rc = -1

        # found valid time unit
        if rc <= 0:
            if hasattr(timeC,'calendar'):
                calendar = timeC.calendar

            # decode all time steps if not swath files, otherwise just first and last step
            axisTmp = timeC
            if ds.isSwathData():
                axisTmp = np.empty(2, dtype=timeC.dtype)
                axisTmp[0] = timeC[0]
                axisTmp[1] = timeC[-1]

            try:
                tSteps = np.empty(axisTmp.shape, dtype=datetime.datetime)
                it = np.nditer(axisTmp, flags=['c_index'])
                while not it.finished:
                    if sinceYr < 1:
                        t = Time(axisTmp[it.index], format='jd', scale='utc', precision=4)
                        t = datetime.datetime.strptime(t.isot, "%Y-%m-%dT%H:%M:%S.%f")
                        if (t is not None):
                            tSteps[it.index] = t.replace(tzinfo=pytz.utc);
                    else:
                        try:
                            t = num2date (axisTmp[it.index], timeC.units, calendar=calendar)
                            if (isinstance(t, datetime.datetime)):
                                t = t.replace(tzinfo=pytz.utc)
                            tmp = np.around(t.microsecond * np.float64(0.01)).astype(np.int64)*100
                            if (t.microsecond > 0):
                                print(f"{'':<8}{RC_WARN} time record not exact (mus={t.microsecond})")
                            t  = t + datetime.timedelta(microseconds=int(tmp-t.microsecond))
                        except ValueError:
                            print(f"{'':<8}{RC_ERR} invalid time record")
                            rc = 1
                            t = None

                        if (t is not None):
                            tSteps[it.index] = t
                    it.iternext()
            except:
                rc = 1
                print(f"{'':<8}{RC_ERR} invalid time axis.")

        # decode time record in file name
        decode = re.match(CMSAF_NAMING_STANDARD, self.File)
        try:
            timeStepFn = datetime.datetime.strptime(self.File[5:17], "%Y%m%d%H%M")
            timeStepFn = timeStepFn.replace(tzinfo=pytz.utc);
        except ValueError:
            print(f"{'':<8}{RC_INFO} Non standard file name '{self.File}'")
            timeStepFn = None

        # checks on time steps
        if (tSteps is not None):
            # compare dates
            if (decode is not None and ds.isSwathData() == False):
                if (timeStepFn is not None and decode.group(9) == "002" and decode.group(1) == "UTH"):
                    timeStepFn = timeStepFn - datetime.timedelta(days=0, hours=0, minutes=30, seconds=0)
                if (timeStepFn is not None) and (timeStepFn != tSteps[0]):
                    print(f"{'':<8}{RC_ERR} time record mismatch, expecting {timeStepFn.isoformat()} as first record")
                    rc = 1
            else:
                if (timeStepFn is not None) and (tSteps[0] is not None) and (timeStepFn > tSteps[0]):
                    print(f"{'':<8}{RC_ERR} time record mismatch, expecting {timeStepFn.isoformat()} before first record")
                    rc = 1

        # test records
        if expRecords is not None:
            if expRecords != tSteps.size:
                print(f"{'':<8}{RC_ERR} Expecting {expRecords} records but found {tSteps.size}.")
                rc = 1

        # test time bounds
        timeBoundsKey = "bounds"
        timeBoundsVar = None
        timeBounds    = None
        if hasattr(timeC,"climatology"):
            timeBoundsKey  = 'climatology'
        if hasattr(timeC,timeBoundsKey):
            tmp = getattr(timeC,timeBoundsKey)
            if timeC.group().name == "/":
                if tmp in ds.variables:
                    timeBoundsVar = ds[tmp]
            else:
                if tmp in ds[timeC.group().path].variables:
                    timeBoundsVar = ds[os.path.join(timeC.group().path,tmp)]
            if timeBoundsVar is None:
                print(f"{'':<8}{RC_ERR} Missing configured bounds variable '{tmp}'.")
                rc = 1

        # check bound units
        if timeBoundsVar is not None:
            if hasattr(timeBoundsVar,'units'):
                if timeBoundsVar.units != tUnits:
                    print(f"{'':<8}{RC_ERR} time bounds must have same axis as time, but found: {timeBoundsVar.units}")
                    rc = 1
            if timeBoundsVar.shape != (timeC[:].size,2):
                print(f"{'':<8}{RC_ERR} time bounds must have shape (size_of_time,2), but found: {timeBoundsVar.shape}")
                rc = 1
            else:
                timeBounds = np.empty(timeBoundsVar.shape, dtype=type(tSteps))
                it = np.nditer(timeBoundsVar, flags=['multi_index'])
                while not it.finished:
                    if sinceYr < 1:
                        t = Time(it[0], format='jd', scale='utc', precision=4)
                        t = datetime.datetime.strptime(t.isot, "%Y-%m-%dT%H:%M:%S.%f")
                        timeBounds[it.multi_index] = t.replace(tzinfo=pytz.utc)
                    else:
                        t = num2date(it[0], tUnits, calendar=calendar)
                        if (isinstance(t, datetime.datetime)):
                            t = t.replace(tzinfo=pytz.utc)
                        timeBounds[it.multi_index] = t
                    it.iternext()

            # check if time bounds are required
            if timeBounds is None:
                expTimeBounds = True

                # no bounds for swath data
                if ds.isSwathData():
                    expTimeBounds = False

                # test for cell_methods in main variables
                elif hasattr(ds,"variable_id"):
                    haveCellMethodsPoint = True
                    for item in ds.variable_id.split(","):
                        hasCellMethodsPoint = False
                        if item in ds.variables:
                            if hasattr(ds.variables[item],"cell_methods"):
                                cellMethods = ds.variables[item].cell_methods
                                cellMethods = re.match('^.*time: *point.*$', cellMethods)
                                if cellMethods is not None:
                                    hasCellMethodsPoint = True
                        if not hasCellMethodsPoint:
                            haveCellMethodsPoint = False
                            break
                    if haveCellMethodsPoint:
                        expTimeBounds = False

                if expTimeBounds:
                    print("f{'':<8}{RC_ERR} Missing time bounds")
                    if not self.lazy:
                        rc = 1
                    else:
                        print(f"{'':<8}Ignoring while beeing lazy")
                else:
                    print(f"{'':<8}{RC_INFO} No time bounds required ")

        # loop records
        if tSteps is not None:
            print("")
            prevTimeEnd = None
            itRecord    = tSteps[0]
            nRecord     = None
            it          = np.nditer(tSteps, flags=['c_index','refs_ok'])
            while not it.finished:
                itRc = 'OK'
                itRecStatus = None

                # update current expected record
                if nRecord is not None:
                    itRecord = nRecord

                # next expected record
                if expResolution is not None:
                    nRecord = itRecord
                    if (expResolution[0] > 0):
                        nRecord = nRecord.replace(year=nRecord.year+int(expResolution[0]))
                    if (expResolution[1] > 0):
                        nm = nRecord.month + int(expResolution[1])
                        ny = nRecord.year
                        if (nm > 12):
                            ny += 1
                            nm -= 12
                        nRecord = nRecord.replace(year=ny,month=nm)
                    tr = datetime.timedelta(days=expResolution[3], hours=expResolution[4], minutes=expResolution[5], seconds=expResolution[6])
                    nRecord += tr

                    # add one day during a leap year for last feb pentad period
                    if expResolution[3] == 5:
                        if (nRecord.timetuple().tm_yday == 61 and leap_year(nRecord.year)):
                            nRecord += datetime.timedelta(days=1)

                # decode matching record status
                rsKey = os.path.join(timeC.group().path, 'record_status')
                if rsKey in recordStatus.keys():
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

                # test time records against time bounds
                if timeBounds is not None:
                    if (timeBounds[it.index,0] > tSteps[it.index]) or \
                       (timeBounds[it.index,1] < tSteps[it.index]):
                        itRc = 'FAILED (record not in bounds)'
                        rc = 1
                    if nRecord is not None:
                        if timeBounds[it.index,1] < nRecord:
                            itRc = 'FAILED (gap in right bound)'
                            rc = 1
                        if timeBounds[it.index,1] > nRecord:
                            itRc = 'FAILED (overlap in right bound)'
                            rc = 1
                    if itRecStatus is not None:
                        itRc = itRc + ' [status='+itRecStatus+']'
                        itRecStatus = None
                    print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} [{timeBounds[it.index,0].isoformat()}, {timeBounds[it.index,1].isoformat()}] -> {itRc}")

                    # test coverage
                    if prevTimeEnd is not None:
                        timeDiff = timeBounds[it.index,0] - prevTimeEnd
                        timeDiff = timeDiff.total_seconds()
                        if timeDiff > 0.:
                            itRc = 'FAILED'
                            rc = 1
                            print(f"{'':<8}{RC_ERR} gap in time coverage {timeDiff} seconds")
                        elif timeDiff < 0.:
                            itRc = 'FAILED'
                            rc = 1
                            print(f"{'':<8}{RC_ERR} overlap in time coverage {timeDiff} seconds")

                # test time records against file name time resolution
                elif expResolution is not None:
                    itRc = 'OK'
                    timeDiff = itRecord - tSteps[it.index]
                    if timeDiff.total_seconds() > 0.:
                        itRc = 'FAILED'
                        rc = 1
                    if itRecStatus is not None:
                        itRc = itRc + ' [status='+itRecStatus+']'
                        itRecStatus = None
                    if timeDiff.total_seconds() == 0.:
                        print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} -> {itRc}")
                    else:
                        print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} expected {itRecord} -> {itRc}")

                else:
                    itRc = 'OK'
                    if tSteps[it.index] is None:
                        itRc = 'FAILED'
                    if itRecStatus is not None:
                        itRc = itRc + ' [status='+itRecStatus+']'
                    print(f"{'':<8}{it.index+1: >3} {tSteps[it.index]} -> {itRc}")

                if timeBounds is not None:
                    prevTimeEnd = timeBounds[it.index,1]

                it.iternext()

            # test time coverage range start
            if hasattr(ds,"time_coverage_start"):
                try:
                    t = num2date(0, "days since {}".format(ds.time_coverage_start))
                    if (isinstance(t, datetime.datetime)):
                        t = t.replace(tzinfo=pytz.utc, microsecond=0)
                    timeCoverStart = t
                except:
                    print("f{'':<8}{RC_ERR} Unexpected time format: '{ds.time_coverage_start}'")
                else:
                    if timeCoverStart > tSteps[0]:
                        print(f"{'':<8}{RC_ERR} first time record '{tSteps[0].isoformat()}' not within time_coverage_start attribute: '{timeCoverStart.isoformat()}'")
                        rc = 1
                    if (timeBounds is not None) and (timeCoverStart != timeBounds[0,0]):
                        print(f"{'':<8}{RC_ERR} time bound [0,0] is not matching time_coverage_start attribute: '{timeCoverStart.isoformat()}'")
                        rc = 1

            # test time coverage range
            if hasattr(ds,"time_coverage_end"):
                try:
                    t = num2date(0, "days since {}".format(ds.time_coverage_end))
                    if (isinstance(t, datetime.datetime)):
                        t = t.replace(tzinfo=pytz.utc, microsecond=0)
                    timeCoverEnd = t
                except:
                    print(f"{'':<8}{RC_ERR} Unexpected time format: {ds.time_coverage_end}")
                else:
                    if timeCoverEnd < tSteps[-1]:
                        print(f"{'':<8}{RC_ERR} last time record '{tSteps[-1].isoformat()}' not within time_coverage_end attribute: '{timeCoverEnd.isoformat()}'")
                        rc = 1
                    if (timeBounds is not None) and (timeBoundsKey != 'climatology') and timeCoverEnd != timeBounds[-1,1]:
                        print(f"{'':<8}{RC_ERR} time bound [-1,1] is not matching time_coverage_end attribute: '{timeCoverEnd.isoformat()}'")
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
            print(f"{'':<8}{RC_ERR} missing mandatory attribute 'axis'")
            rc = 1
        elif (expAxis is not None):
            if coordVar.axis != expAxis:
                print(f"{'':<8}{RC_ERR} invalid value attribute 'axis={coordVar.axis}'")
                rc = 1

        # exclude from checks if not fixed
        if axisTime.name in coordVar.dimensions:
            print(f"{'':<4} not fixed")
            rc = 1

        # test axis values
        if rc == 0:
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
                except:
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
                except:
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
            if (coordRes is not None):
                if (np.ma.count_masked(coord) > 0):
                    print(f"{'':<8}{RC_ERR} {longName} contains missing data")
                    rc = 1
                else:
                    finfo    = np.finfo(coord.dtype)
                    eps      = significant_digits(coord,finfo.precision)*finfo.resolution
                    indx     = np.arange(np.int64(0),len(coord),dtype=np.int64)
                    meshExp  = np.multiply(indx, np.multiply(np.int64(10000),coordRes)).astype(np.int64)
                    tmp      = np.rint(1./eps[0]).astype(np.float64)
                    tmp      = np.divide(np.rint(np.divide(coord[0], eps[0])), tmp)
                    tmp      = np.multiply(tmp, np.int64(10000)).astype(np.int64)
                    meshExp  = np.multiply(np.add(tmp,meshExp),np.float64(0.0001))
                    tmp      = np.absolute(np.subtract(meshExp, coord))
                    eps      = float_spacing(meshExp, coord)
                    indx     = np.where(tmp > eps)[0]
                    if (len(indx) > 0):
                        print("f{'':<8}{RC_ERR} at {len(indx)} locations:")
                        print("f{'':<10}found:    {coord[indx]}")
                        print("f{'':<10}expecting:{meshExp[indx]}")
                        if (not self.lazy):
                            rc = 1
                        else:
                            print(f"{'':<8}{RC_INFO} Ignoring while beeing lazy")

                    # test 0,0 [must not be in center]
                    tmp  = np.absolute(coord)
                    indx = np.where(tmp < eps)[0]
                    if (len(indx) > 0):
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
                        boundsVar = ds[tmp]
                else:
                    if tmp in ds[coordVar.group().path].variables:
                        boundsVar = ds[coordVar.group().path].variables[tmp]
                if boundsVar is not None:
                    if (coordOrder == 1):
                        bounds = boundsVar[:]
                    else:
                        bounds = np.flip(boundsVar[:])
                else:
                    print(f"{'':<8}{RC_ERR} Missing configured bounds variable '{tmp}'.")
                    rc = 1

            if (bounds is None):
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

                print(f"{'':<8}[{bounds[0,0]} -> {bounds[-1,0]}] {leftName} bounds")
                print(f"{'':<8}[{bounds[0,1]} -> {bounds[-1,1]}] {rightName} bounds")

                # test coordinate within bounds
                indx1 = np.where(bounds[:,0] > coord)[0]
                indx2 = np.where(bounds[:,1] < coord)[0]
                if (len(indx1) > 0 or len(indx2) > 0):
                    rc = 1
                    print(f"{'':<8}{RC_ERR} {longName} values not within bounds")

                # test for gaps in coordinate bounds
                tmp = np.subtract(boundsVar[0:-2,1],boundsVar[1:-1,0])
                indx = np.where(tmp > 0)[0]
                if (len(indx) > 0):
                    rc = 1
                    print(f"{'':<8}{RC_ERR} gaps in {longName} bounds")

                # test for ovarlap in coordinate bounds
                indx = np.where(tmp < 0)[0]
                if (len(indx) > 0):
                    rc = 1
                    print("{:<8}{RC_ERR} {longName} bounds overlap")

                # test bounds against global attribute
                if (geoMinAttr is not None):
                    if (geoMinAttr != boundsVar[0,0]):
                        rc = 1
                        print(f"{'':<8}{RC_ERR} mismatch between {leftMaxName} {longName} bound '{boundsVar[0,0]}' and {geoMinAttrName} '{geoMinAttr}'")
                if (geoMaxAttr is not None):
                    if (geoMaxAttr != boundsVar[-1,1]):
                        rc = 1
                        print("{:<8}{RC_ERR} mismatch between {rightMaxName} {longName} bound '{boundsVar[-1,1]}' and {geoMaxAttrName} '{geoMaxAttr}'")

            # print result
            if (coordRes is not None):
                print(f"{'':<8}[{coordMin} -> {coordMax} by {coordRes}]")
            else:
                print(f"{'':<8}[{coordMin} -> {coordMax}]")

        return(rc)


    def _checkCompression(self):
        """
        Check compression settings for all 3-dimensional variables.
        """

        rc = 0
        ds = self.Dataset

        if (ds.data_model == "NETCDF4") or (ds.data_model == "NETCDF4_CLASSIC"):
            vList = ds.getVariableList()

            for vName in vList:
                if (len(ds[vName].shape) >= 3):
                    xRc = 1
                    filters = ds[vName].filters()
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

        # record_status must be defined
        recordStatus = ds.getVariableByName("record_status")
        if len(recordStatus) > 0:
            for key in recordStatus:
                item = recordStatus[key]
                if (not (item.dtype == np.dtype(np.uint8) or
                         item.dtype == np.dtype(np.int8)  or
                         item.dtype == np.dtype(np.int16) or
                         item.dtype == np.dtype(np.int32)    ) ):
                    rc = 1
                    print(f"{RC_ERR} incorrect data type of '{key}'\n")

        # no record status for swath data required
        elif (not ds.isSwathData()):
            rc = 1
            print(f"{RC_ERR} missing mandatory variable 'record_status'\n")

        # recommended attributes
        listRec = ['units', 'standard_name', 'grid_mapping']
        listRecAxis = ['grid_mapping']

        # mandatory attributes
        listMan = ['long_name']

        # skip variables
        listSkip = ['record_status', 'latlon_grid']

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
            var = ds[vName]
            if var.group().name != "/":
                grp = ds.groups[var.group().name]
            else:
                grp = ds

            if hasattr(var, 'grid_mapping_name'):
                continue

            # find dimensions
            aLon  = None
            aLat  = None
            aTime = None
            for dim in var.get_dims():
                if (aLon is None):
                    aLon = ds.matchCoordinate(dim, axisLon)
                if (aLat is None):
                    aLat = ds.matchCoordinate(dim, axisLat)
                if (aTime is None):
                    aTime = ds.matchCoordinate(dim, axisTime)

            # test dimensions, all variables with lat and lon must also have a time coordinate
            if (aLon is not None and aLat is not None):
                if (aTime is None):
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
                elif not hasattr(var, item):
                    itRc = 1
                    for gvName in grp.variables:
                        gv = grp.variables[gvName]
                        if (hasattr(gv, 'bounds') and gv.bounds == vNameB):
                            itRc = 0
                            break
                if (itRc == 1):
                    if item in listMan:
                        print(f"{'':<4}{RC_ERR} {vName} :: missing mandatory attribute '{item}'")
                    elif item in listRec:
                        print(f"{'':<4}{RC_WARN} {vName} :: missing recommended attribute '{item}'")

            # test flags
            if hasattr(var,'flag_values') and hasattr(var,'flag_meanings'):
                flagV = var.flag_values
                flagM = var.flag_meanings.split(" ")
                if (len(flagV) != len(flagM)):
                    rc = 1
                    print(f"{'':<4}{RC_ERR} {vName} :: mismatch between flag_values and flag_value")

        # test variables defined in variable_id
        if hasattr(ds, 'variable_id'):
            print(f"\nvariable_id :: '{ds.variable_id}'")
            for item in ds.variable_id.split(","):
                try:
                    var = ds[item]
                except:
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

    # get standards
    prefix = getenv('CMSAF_CHECKER_PREFIX')
    if (prefix is not None):
        STANDARD = prefix + '/share/'
    else:
        STANDARD = os.path.dirname(__file__)
        STANDARD = os.path.join(os.path.dirname(STANDARD),"share")

    print(f"CMSAF Checker Version {__version__}")
    print(f"CMSAF Standards path '{STANDARD}'")

    # argument parser
    parser = argparse.ArgumentParser(prog='cmsaf-checker')
    parser.add_argument('-s', '--cmsaf_metadata_standard',
        default=STANDARD, help='location of the CM SAF Metadata Standards')
    parser.add_argument('-v', '--version',
        help='CM SAF standards version to apply')
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

    # get a new checker object
    if (args.reference == None):
        inst = CMSAFChecker(metadataStandards=args.cmsaf_metadata_standard, version=args.version,
            coordinates=args.coordinates, lazy=args.lazy, ignore=args.ignore_attr)
    else:
        inst = CMSAFChecker(referenceFile=args.reference, ignore=args.ignore_attr,
            coordinates=args.coordinates, lazy=args.lazy)

    if (inst.refDataset == None):
        print(f"Using CM SAF Metadata Standard Version {inst.std_name_dh.version_number} ({inst.std_name_dh.last_modified})")

    # file pattern expansion
    files = []
    if args.directory is not None:
        for root, dirs, names in os.walk(args.directory, followlinks=True):
            tmp = fnmatch.filter(names, args.files[0])
            for fn in tmp:
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
    if (len(files) > 1 and args.missing is not None):
        if (args.missing == "filename"):
            fileName = os.path.realpath(files[0])
            fileName = os.path.basename (fileName)
            fileAttr = re.match(CMSAF_NAMING_STANDARD, fileName)
            fileStep = fileAttr.group(2)
        else:
            fileStep = args.missing
        if (fileStep == 'd'):
            fileDelta = datetime.timedelta(days=1, hours=0, minutes=0, seconds=0)
        elif (fileStep == 'm'):
            fileDelta = datetime.timedelta(days=31, hours=0, minutes=0, seconds=0)
        elif (fileStep == 'h'):
            fileDelta = datetime.timedelta(days=0, hours=1, minutes=0, seconds=0)
        elif (fileStep == 'M15'):
            fileDelta = datetime.timedelta(days=0, hours=0, minutes=15, seconds=0)

    # loop files
    lastTime = None
    for index, file in enumerate(files):
        # test missing
        if (fileDelta is not None):
            currentFile = os.path.realpath(file)
            currentFile = os.path.basename (currentFile)
            currentTime = datetime.datetime.strptime(currentFile[5:17], "%Y%m%d%H%M")
            currentTime = currentTime.replace(tzinfo=pytz.utc)
            if (lastTime is not None):
                while 1:
                    nextTime = lastTime + fileDelta
                    if (fileStep == 'm'):
                        nextTime = nextTime.replace(day=1)
                    fileDiff = currentTime-nextTime
                    if (fileDiff.total_seconds() > 0):
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

    # final result
    print(f"\n{'':=^80}\nOverall Summary\n{'':=^80}")
    print(f"Out of {len(files)}, {res['FAILED']} FAILED")
    if (args.missing) and (fileDelta is not None):
        print(f"{res['MISSING']} files MISSING")

    if (res['FAILED'] == 0):
        exit(0)
    else:
        exit(1)

if __name__ == '__main__':
    main()

# ex: set foldmethod=indent:
