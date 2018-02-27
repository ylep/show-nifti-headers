#! /usr/bin/env python2
# -*- coding: utf-8 -*-

"""Print the header of a NIFTI-1 file in human-readable form.

Gzip-compressed files are recognized if their name ends in '.gz'.
Standard input is used if the file name is a single hyphen '-'.
"""

# This program is based on the NIfTI-1.1 specification as described on
# http://nifti.nimh.nih.gov/nifti-1

from __future__ import print_function, unicode_literals
from future_builtins import zip

import collections
import exceptions
import logging
import math
import os.path
import struct
import sys


logging.basicConfig(format='%(message)s', level=logging.INFO)


class NiftiFormatError(exceptions.Exception):
    """Exception signalling errors encountered during NIfTI file parsing"""
    pass


class NotNiftiError(NiftiFormatError):
    """Exception signalling that a file is not in NIfTI format"""
    pass


class InconsistentNiftiError(NiftiFormatError):
    """Exception signalling an inconsistency in a NIfTI file"""
    pass


# These tables represent static data from the official definition of
# NIfTI-1 <http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h>
NIFTI1_FIELDS_DESCRIPTION = [
    ('i', 'sizeof_hdr', "MUST be 348"),
    ('10s', 'data_type', "++UNUSED++"),
    ('18s', 'db_name', "++UNUSED++"),
    ('i', 'extents', "++UNUSED++"),
    ('h', 'session_error', "++UNUSED++"),
    ('b', 'regular', "++UNUSED++"),
    ('b', 'dim_info', "MRI slice ordering"),
    ('8h', 'dim', "Data array dimensions"),
    ('f', 'intent_p1', "1st intent parameter"),
    ('f', 'intent_p2', "2nd intent parameter"),
    ('f', 'intent_p3', "3rd intent parameter"),
    ('h', 'intent_code', "NIFTIINTENT code"),
    ('h', 'datatype', "Defines data type"),
    ('h', 'bitpix', "Number bits/voxel"),
    ('h', 'slice_start', "First slice index"),
    ('8f', 'pixdim', "Grid spacings"),
    ('f', 'vox_offset', "Offset into .nii file"),
    ('f', 'scl_slope', "Data scaling: slope"),
    ('f', 'scl_inter', "Data scaling: offset"),
    ('h', 'slice_end', "Last slice index"),
    ('b', 'slice_code', "Slice timing order"),
    ('b', 'xyzt_units', "Units of pixdim[1..4]"),
    ('f', 'cal_max', "Max display intensity"),
    ('f', 'cal_min', "Min display intensity"),
    ('f', 'slice_duration', "Time for 1 slice"),
    ('f', 'toffset', "Time axis shift"),
    ('i', 'glmax', "++UNUSED++"),
    ('i', 'glmin', "++UNUSED++"),
    ('80s', 'descrip', "any text you like"),
    ('24s', 'aux_file', "auxiliary filename"),
    ('h', 'qform_code', "NIFTIXFORM code"),
    ('h', 'sform_code', "NIFTIXFORM code"),
    ('f', 'quatern_b', "Quaternion b param"),
    ('f', 'quatern_c', "Quaternion c param"),
    ('f', 'quatern_d', "Quaternion d param"),
    ('f', 'qoffset_x', "Quaternion x shift"),
    ('f', 'qoffset_y', "Quaternion y shift"),
    ('f', 'qoffset_z', "Quaternion z shift"),
    ('4f', 'srow_x', "1st row affine transform"),
    ('4f', 'srow_y', "2nd row affine transform"),
    ('4f', 'srow_z', "3rd row affine transform"),
    ('16s', 'intent_name', "name or meaning of data"),
    ('4s', 'magic', 'MUST be "ni1\\0" or "n+1\\0"'),
]

# Concatenate all field descriptions in struct module syntax
NIFTI1_STRUCT_FORMAT = "".join(next(zip(*NIFTI1_FIELDS_DESCRIPTION)))
# The call to str() is necessary in Python 2.6 because struct does not
# accept unicode strings
assert struct.calcsize(str('=' + NIFTI1_STRUCT_FORMAT)) == 348


# These tables represent static data from the official definition of
# NIfTI-2 <https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h>
NIFTI2_FIELDS_DESCRIPTION = [
    ('i', 'sizeof_hdr', "MUST be 540"),
    ('8s', 'magic', 'MUST be valid signature'),
    ('h', 'datatype', "Defines data type"),
    ('h', 'bitpix', "Number bits/voxel"),
    ('8q', 'dim', "Data array dimensions"),
    ('d', 'intent_p1', "1st intent parameter"),
    ('d', 'intent_p2', "2nd intent parameter"),
    ('d', 'intent_p3', "3rd intent parameter"),
    ('8d', 'pixdim', "Grid spacings"),
    ('q', 'vox_offset', "Offset into .nii file"),
    ('d', 'scl_slope', "Data scaling: slope"),
    ('d', 'scl_inter', "Data scaling: offset"),
    ('d', 'cal_max', "Max display intensity"),
    ('d', 'cal_min', "Min display intensity"),
    ('d', 'slice_duration', "Time for 1 slice"),
    ('d', 'toffset', "Time axis shift"),
    ('q', 'slice_start', "First slice index"),
    ('q', 'slice_end', "Last slice index"),
    ('80s', 'descrip', "any text you like"),
    ('24s', 'aux_file', "auxiliary filename"),
    ('i', 'qform_code', "NIFTIXFORM code"),
    ('i', 'sform_code', "NIFTIXFORM code"),
    ('d', 'quatern_b', "Quaternion b param"),
    ('d', 'quatern_c', "Quaternion c param"),
    ('d', 'quatern_d', "Quaternion d param"),
    ('d', 'qoffset_x', "Quaternion x shift"),
    ('d', 'qoffset_y', "Quaternion y shift"),
    ('d', 'qoffset_z', "Quaternion z shift"),
    ('4d', 'srow_x', "1st row affine transform"),
    ('4d', 'srow_y', "2nd row affine transform"),
    ('4d', 'srow_z', "3rd row affine transform"),
    ('i', 'slice_code', "Slice timing order"),
    ('i', 'xyzt_units', "Units of pixdim[1..4]"),
    ('i', 'intent_code', "NIFTI_INTENT_* code"),
    ('16s', 'intent_name', "'name' or meaning of data"),
    ('b', 'dim_info', "MRI slice ordering"),
    ('15s', 'unused_str', "unused, filled with \\0"),
]

# Concatenate all field descriptions in struct module syntax
NIFTI2_STRUCT_FORMAT = "".join(next(zip(*NIFTI2_FIELDS_DESCRIPTION)))
# The call to str() is necessary in Python 2.6 because struct does not
# accept unicode strings
assert struct.calcsize(str('=' + NIFTI2_STRUCT_FORMAT)) == 540


INTENT_DICT = {
    0: 'NIFTI_INTENT_NONE',
    2: 'NIFTI_INTENT_CORREL',
    3: 'NIFTI_INTENT_TTEST',
    4: 'NIFTI_INTENT_FTEST',
    5: 'NIFTI_INTENT_ZSCORE',
    6: 'NIFTI_INTENT_CHISQ',
    7: 'NIFTI_INTENT_BETA',
    8: 'NIFTI_INTENT_BINOM',
    9: 'NIFTI_INTENT_GAMMA',
    10: 'NIFTI_INTENT_POISSON',
    11: 'NIFTI_INTENT_NORMAL',
    12: 'NIFTI_INTENT_FTEST_NONC',
    13: 'NIFTI_INTENT_CHISQ_NONC',
    14: 'NIFTI_INTENT_LOGISTIC',
    15: 'NIFTI_INTENT_LAPLACE',
    16: 'NIFTI_INTENT_UNIFORM',
    17: 'NIFTI_INTENT_TTEST_NONC',
    18: 'NIFTI_INTENT_WEIBULL',
    19: 'NIFTI_INTENT_CHI',
    20: 'NIFTI_INTENT_INVGAUSS',
    21: 'NIFTI_INTENT_EXTVAL',
    22: 'NIFTI_INTENT_PVAL',
    23: 'NIFTI_INTENT_LOGPVAL',
    24: 'NIFTI_INTENT_LOG10PVAL',
    1001: 'NIFTI_INTENT_ESTIMATE',
    1002: 'NIFTI_INTENT_LABEL',
    1003: 'NIFTI_INTENT_NEURONAME',
    1004: 'NIFTI_INTENT_GENMATRIX',
    1005: 'NIFTI_INTENT_SYMMATRIX',
    1006: 'NIFTI_INTENT_DISPVECT',
    1007: 'NIFTI_INTENT_VECTOR',
    1008: 'NIFTI_INTENT_POINTSET',
    1009: 'NIFTI_INTENT_TRIANGLE',
    1010: 'NIFTI_INTENT_QUATERNION',
    1011: 'NIFTI_INTENT_DIMLESS',
    2001: 'NIFTI_INTENT_TIME_SERIES',
    2002: 'NIFTI_INTENT_NODE_INDEX',
    2003: 'NIFTI_INTENT_RGB_VECTOR',
    2004: 'NIFTI_INTENT_RGBA_VECTOR',
    2005: 'NIFTI_INTENT_SHAPE',
}

DATATYPE_INFO_DICT = {
    0: ('', None),
    1: ('DT_BINARY', 1),
    2: ('NIFTI_TYPE_UINT8', 8),
    4: ('NIFTI_TYPE_INT16', 16),
    8: ('NIFTI_TYPE_INT32', 32),
    16: ('NIFTI_TYPE_FLOAT32', 32),
    32: ('NIFTI_TYPE_COMPLEX64', 64),
    64: ('NIFTI_TYPE_FLOAT64', 64),
    128: ('NIFTI_TYPE_RGB24', 24),
    255: ('DT_ALL', None),
    256: ('NIFTI_TYPE_INT8', 8),
    512: ('NIFTI_TYPE_UINT16', 16),
    768: ('NIFTI_TYPE_UINT32', 32),
    1024: ('NIFTI_TYPE_INT64', 64),
    1280: ('NIFTI_TYPE_UINT64', 64),
    1536: ('NIFTI_TYPE_FLOAT128', 128),
    1792: ('NIFTI_TYPE_COMPLEX128', 128),
    2048: ('NIFTI_TYPE_COMPLEX256', 256),
    2304: ('NIFTI_TYPE_RGBA32', 32),
}

XYZT_UNITS_DICT = {
    0: 'NIFTI_UNITS_UNKNOWN',
    1: 'NIFTI_UNITS_METER',
    2: 'NIFTI_UNITS_MM',
    3: 'NIFTI_UNITS_MICRON',
    8: 'NIFTI_UNITS_SEC',
    16: 'NIFTI_UNITS_MSEC',
    24: 'NIFTI_UNITS_USEC',
    32: 'NIFTI_UNITS_HZ',
    40: 'NIFTI_UNITS_PPM',
    48: 'NIFTI_UNITS_RADS',
}

SPACE_UNIT_TEXT_DICT = {
    1: " m",
    2: " mm",
    3: " µm",
}

TIME_UNIT_TEXT_DICT = {
    8: " s",
    16: " ms",
    24: " µs",
    32: " Hz",
    40: " ppm",
    48: " rad",
}

SLICE_CODE_DICT = {
    0: 'NIFTI_SLICE_UNKNOWN',
    1: 'NIFTI_SLICE_SEQ_INC',
    2: 'NIFTI_SLICE_SEQ_DEC',
    3: 'NIFTI_SLICE_ALT_INC',
    4: 'NIFTI_SLICE_ALT_DEC',
}

XFORM_CODE_DICT = {
    0: 'NIFTI_XFORM_UNKNOWN',
    1: 'NIFTI_XFORM_SCANNER_ANAT',
    2: 'NIFTI_XFORM_ALIGNED_ANAT',
    3: 'NIFTI_XFORM_TALAIRACH',
    4: 'NIFTI_XFORM_MNI_152',
}


class NiftiHeader(object):
    @classmethod
    def _guess_byte_order(cls, binary_header):
        """Guess the byte order of a raw binary NIfTI header.

        Return '<' for little-endian, '>' for big-endian (suitable for
        configuring the struct module).
        """
        native = '>' if sys.byteorder == 'big' else '<'
        if cls._test_byte_order(binary_header, native):
            return '<'
        else:
            other = '<' if sys.byteorder == 'big' else '>'
            if cls._test_byte_order(binary_header, other):
                return '>'
            else:
                raise InconsistentNiftiError("dim[0] must be in range 1..7")

    @staticmethod
    def _test_magic_string_version(magic_field):
        """Test NIfTI magic string, returning NIfTI version or None."""
        if magic_field[0:1] != b'n' or magic_field[3:4] != b'\0':
            return None
        if magic_field[1:2] not in (b'i', b'+'):
            return None
        if b'1' <= magic_field[2:3] <= b'9':
            return int(str(magic_field[2:3]))
        # else return None is implied

    def __init__(self, source):
        """Initialize header data.

        The data source can be a raw binary header contained in a
        Python byte string, or a mapping containing raw header fields.
        """
        # Once Python 2.6 support can be dropped, raw could be
        # of type collections.OrderedDict
        if isinstance(source, collections.Mapping):
            self.raw = dict(source)
        else:  # TODO check if this is bytes type (bytes/bytearray)
            self.raw = {}
            self.from_binary(source)
        self.check_format()

    def read_binary_extender(self, binary_extender):
        if binary_extender:
            return (binary_extender[0:1] != b'\0')
        else:
            return False

    def list_inconsistencies(self):
        """List the inconsistencies of the header's raw data."""
        self.check_format()

        if not 1 <= self.dim[0] <= 7:
            yield InconsistentNiftiError("dim[0] must be in range 1..7")

        if self.onefile and not self.vox_offset >= self.TOTAL_HEADER_SIZE:
            yield InconsistentNiftiError(
                "vox_offset should be {0} or more (is {1})"
                .format(self.TOTAL_HEADER_SIZE, self.vox_offset))

        datatype = self.datatype
        try:
            datatype_info = DATATYPE_INFO_DICT[datatype]
        except KeyError:
            yield InconsistentNiftiError(
                "unknown datatype {0}".format(datatype))

        if self.bitpix != datatype_info[1]:
            yield InconsistentNiftiError(
                "bitpix ({0}) does not match datatype {1}"
                .format(self.bitpix, datatype_info[0]))

        raw_vox_offset = self.raw['vox_offset']
        if raw_vox_offset != round(raw_vox_offset):
            yield InconsistentNiftiError("vox_offset must be an integer")

        # It is guaranteed by the byte order check that 1 <= dim[0] <= 7
        for i, dim_i in enumerate(self.dim[1:(self.dim[0] + 1)], start=1):
            if not dim_i > 0:
                yield InconsistentNiftiError(
                    "dim[{0}] must be positive (is {1})".format(i, dim_i))

        if self.vox_offset % 16 != 0:
            yield InconsistentNiftiError(
                "vox_offset should be an integral multiple of 16")

        if self.intent_code not in INTENT_DICT:
            yield InconsistentNiftiError(
                "unknown intent_code value {0}".format(self.intent_code))
        # Based on the intent code, other checks could be done
        # (e.g. dimensionality)

    @property
    def nifti_version(self):
        """Version of the NIfTI format used."""
        return self._test_magic_string_version(self.raw['magic'])

    @property
    def vox_offset(self):
        """The vox_offset field as an integer."""
        return int(self.raw['vox_offset'])

    @property
    def descrip(self):
        """The descrip field as a byte string."""
        return self.raw['descrip'].rstrip(b'\0')

    @property
    def aux_file(self):
        """The aux_file field as a byte string."""
        return self.raw['aux_file'].rstrip(b'\0')

    @property
    def intent_name(self):
        """The intent_name field as a byte string."""
        return self.raw['intent_name'].rstrip(b'\0')

    def __getattr__(self, name):
        """Access raw header fields."""
        try:
            return self.raw[name]
        except KeyError:
            raise AttributeError("no such attribute {0!r}".format(name))

    def __dir__(self):
        type_dir = set(dir(type(self)))
        dict_attributes = set(self.__dict__)
        raw_attributes = set(self.raw.iterkeys())
        return list(type_dir | dict_attributes | raw_attributes)

    def print_raw(self, file=sys.stdout, describe_fields=False):
        """Print raw header fields, optionally with description."""
        for _, name, description in self.FIELDS_DESCRIPTION:
            value = self.raw[name]
            if describe_fields:
                print("{0:14} {1!r} [{2}]".format(name, value, description),
                      file=file)
            else:
                print("{0:14} {1!r}".format(name, value), file=file)

    def print_interpreted(self, file=sys.stdout):
        """Print interpreted header data on the supplied stream."""
        print("General file information", file=file)
        print("===========================", file=file)
        print("NIfTI version: {0}".format(self.nifti_version), file=file)
        if self.onefile:
            print("Data is stored in this file", file=file)
            if self.vox_offset > self.TOTAL_HEADER_SIZE:
                print("{0} bytes of padding and/or extensions exist after "
                      "the header".format(self.vox_offset
                                          - self.TOTAL_HEADER_SIZE))
        else:
            print("Data is stored in an associated .img file", file=file)
            if self.vox_offset != 0:
                print("Data is stored at {0} bytes offset"
                      .format(self.vox_offset), file=file)

        print(file=file)
        print("General dataset information", file=file)
        print("===========================", file=file)
        print("datatype = {0}".format(self.readable_datatype), file=file)
        print("dim: {0}".format(self.meaningful_dim), file=file)

        readable_xyzt_units = " | ".join(self.readable_xyzt_units)
        print("xyzt_units = {0}".format(readable_xyzt_units), file=file)

        print(file=file)
        print("3D image (volume) orientation and location in space", file=file)
        print("===================================================", file=file)
        print(file=file)
        print("Method 1", file=file)
        print("--------", file=file)

        readable_pixdim = "(" + ", ".join(self.readable_pixdim) + ")"
        print("pixdim: {0}".format(readable_pixdim), file=file)

        if self.qform_code > 0:
            print(file=file)
            print("Method 2", file=file)
            print("--------", file=file)
            print("qform_code = {0}"
                  .format(self.readable_qform_code), file=file)

            print("pixdim: {0}".format(readable_pixdim), file=file)
            print("qfac = {0}".format(self.readable_qfac), file=file)

            print("quatern_abcd: {0}".format(self.quatern_abcd), file=file)
            print("qoffset_xyz: {0}".format(self.qoffset_xyz), file=file)

            if not self.qfac_is_valid:
                print("WARNING: the following matrix could be wrong because "
                      "the value of qfac is invalid, assuming qfac = 1",
                      file=file)
                logging.warn("the transformation matrix")
            print("Equivalent affine matrix:", file=file)
            print_affine_matrix(self.affine_matrix_method2, file=file)

        if self.sform_code > 0:
            print(file=file)
            print("Method 3", file=file)
            print("--------", file=file)
            print("sform_code = {0}".format(self.readable_sform_code),
                  file=file)

            print("Affine matrix:", file=file)
            print_affine_matrix(self.affine_matrix_method3, file=file)

        if self.intent_code or self.intent_name:
            print(file=file)
            print("Interpretation of voxel data", file=file)
            print("============================", file=file)

            print("intent_code = {0}".format(self.readable_intent_code),
                  file=file)

            if self.intent_code:
                intent_p = (self.intent_p1, self.intent_p2,
                            self.intent_p3)
                print("intent_p = {0}".format(intent_p), file=file)

            if self.intent_name:
                print("intent_name = {0}".format(self.intent_name), file=file)

        if self.dim_info or (self.slice_dim and self.slice_duration > 0):
            print(file=file)
            print("MRI-specific spatial and temporal information", file=file)
            print("=============================================", file=file)
            print("dim_info = {0}".format(self.readable_dim_info), file=file)

            if self.slice_dim and self.slice_duration > 0:
                print("slice_duration = {0}"
                      .format(self.slice_duration), file=file)

            if self.slice_code and self.slice_dim and self.slice_duration > 0:
                print("slice_code = {0}"
                      .format(self.readable_slice_code), file=file)

            if self.slice_code:
                print("slice_start = {0}".format(self.slice_start), file=file)
                print("slice_end = {0}".format(self.slice_end), file=file)

        misc_title = {"printed": False}
        def print_misc_title():
            """Print the title as needed, at most once."""
            # Work around the absence of a "nonlocal" keyword in Python2 by
            # using a dictionary
            if not misc_title["printed"]:
                print(file=file)
                print("Miscellaneous", file=file)
                print("=============", file=file)
                misc_title["printed"] = True

        if self.descrip:
            print_misc_title()
            print("descrip = {0!r}".format(self.descrip), file=file)

        if self.aux_file:
            print_misc_title()
            print("aux_file = {0!r}".format(self.aux_file), file=file)

        if self.scl_slope and (self.scl_slope != 1.0 or self.scl_inter != 0):
            print_misc_title()
            print("scl_slope = {0}".format(self.scl_slope), file=file)
            print("scl_inter = {0}".format(self.scl_inter), file=file)

        if self.toffset:
            print_misc_title()
            print("toffset = {0}".format(self.toffset), file=file)

        if self.cal_min or self.cal_max:
            print_misc_title()
            print("cal_min = {0}".format(self.cal_min), file=file)
            print("cal_max = {0}".format(self.cal_max), file=file)

        if self.extensions_present:
            print("WARNING: unsupported NIfTI extensions were detected.\n"
                  "The data or metadata contained in these extensions"
                  " cannot be interpreted.\n", file=file)
            logging.warn("unsupported NIfTI extensions are present")

    @property
    def readable_qform_code(self):
        """Readable value of the qform_code field."""
        try:
            return XFORM_CODE_DICT[self.qform_code]
        except KeyError:
            return "{0} /* invalid qform_code */".format(self.qform_code)

    @property
    def readable_sform_code(self):
        """Readable value of the sform_code field."""
        try:
            return XFORM_CODE_DICT[self.sform_code]
        except KeyError:
            return "{0} /* invalid sform_code */".format(self.sform_code)

    @property
    def qfac(self):
        """Value of qfac (used in method 2)."""
        raw_qfac = self.pixdim[0]
        if raw_qfac == -1.0:
            return -1
        else:
            return 1

    @property
    def qfac_is_valid(self):
        """Assert if the value of qfac is stored in a valid way."""
        return self.pixdim[0] in (1.0, -1.0)

    @property
    def readable_qfac(self):
        """Readable value of qfac."""
        if self.qfac_is_valid:
            return str(self.qfac)
        else:
            return ("1 /* assumed, invalid pixdim[0] = {0} */"
                    .format(self.pixdim[0]))

    @property
    def quatern_abcd(self):
        """Quaternion values (a, b, c, d)."""
        b = self.quatern_b
        c = self.quatern_c
        d = self.quatern_d
        a2 = 1.0 - (b * b + c * c + d * d)
        if a2 >= 0:
            a = math.sqrt(a2)
        else:
            a = float("NaN")
        return (a, b, c, d)

    @property
    def qoffset_xyz(self):
        """Quaternion offsets (x, y, z)."""
        return (self.qoffset_x, self.qoffset_y, self.qoffset_z)

    @property
    def affine_matrix_method2(self):
        """Equivalent affine matrix for method 2 quaternion data.

        The calculation is as described in the reference NIfTI header,
        under the name "Method 2".
        """
        (a, b, c, d) = self.quatern_abcd

        pi = self.pixdim[1]
        pj = self.pixdim[2]
        pk = self.pixdim[3] * self.qfac

        row_x = (pi * (a * a + b * b - c * c - d * d),
                 pj * (2 * b * c - 2 * a * d),
                 pk * (2 * b * d + 2 * a * c),
                 self.qoffset_x)
        row_y = (pi * (2 * b * c + 2 * a * d),
                 pj * (a * a + c * c - b * b - d * d),
                 pk * (2 * c * d - 2 * a * b),
                 self.qoffset_y)
        row_z = (pi * (2 * b * d - 2 * a * c),
                 pj * (2 * c * d + 2 * a * b),
                 pk * (a * a + d * d - c * c - b * b),
                 self.qoffset_z)
        return (row_x, row_y, row_z)

    @property
    def affine_matrix_method3(self):
        """Affine matrix of method 3."""
        return (self.srow_x, self.srow_y, self.srow_z)

    @property
    def slice_dim(self):
        """MRI slice selection direction."""
        return (self.dim_info >> 4) & 0x03

    @property
    def phase_dim(self):
        """MRI phase encoding direction."""
        return (self.dim_info >> 2) & 0x03

    @property
    def freq_dim(self):
        """MRI frequency encoding direction."""
        return self.dim_info & 0x03

    @property
    def readable_slice_code(self):
        """Readable flag for the slice_code field."""
        try:
            return SLICE_CODE_DICT[self.slice_code]
        except KeyError:
            return "{0} /* unknown slice code */".format(self.slice_code)

    @property
    def readable_datatype(self):
        """Readable flag for the datatype field."""
        try:
            return DATATYPE_INFO_DICT[self.datatype][0]
        except KeyError:
            return "{0} /* unknown datatype value */".format(self.datatype)

    @property
    def meaningful_dim(self):
        """Meaningful list of sizes along each dimension of the dataset."""
        return self.dim[1:(self.dim[0] + 1)]

    @property
    def readable_pixdim(self):
        """Yield readable values for the pixdim field, with units."""
        for i in xrange(1, self.dim[0] + 1):
            if 1 <= i <= 3:  # space dimensions
                yield "{0:.3}{1}".format(self.pixdim[i], self.space_unit_text)
            elif i == 4:  # time dimension
                yield "{0:.3}{1}".format(self.pixdim[i], self.time_unit_text)
            else:  # supplementary dimensions
                yield "{0}".format(self.pixdim[i])

    @property
    def readable_xyzt_units(self):
        """Readable value for the xyzt_units field."""
        if self.xyzt_units == 0x00:
            yield XYZT_UNITS_DICT[0]
        else:
            space_unit = self.xyzt_units & 0x07
            time_unit = self.xyzt_units & 0x38
            if space_unit:
                try:
                    yield XYZT_UNITS_DICT[space_unit]
                except KeyError:
                    yield ("0x{0:02X} /* invalid space unit */"
                           .format(space_unit))
            if time_unit:
                try:
                    yield XYZT_UNITS_DICT[time_unit]
                except KeyError:
                    yield "0x{0:02X} /* invalid time unit */".format(time_unit)
            if self.xyzt_units & ~0x3F:
                yield ("0x{0:X} /* invalid xyzt_unit bits */"
                       .format(self.xyzt_units))

    @property
    def space_unit_text(self):
        """Unit along space dimensions as a SI abbreviation."""
        return SPACE_UNIT_TEXT_DICT.get(self.xyzt_units & 0x07, "")

    @property
    def time_unit_text(self):
        """Unit along time dimension as a SI abbreviation."""
        return TIME_UNIT_TEXT_DICT.get(self.xyzt_units & 0x38, "")

    @property
    def separate_dim_info(self):
        """Yield readable sub-values of the dim_info field."""
        if self.freq_dim:
            yield "{0} /* frequency */".format(self.freq_dim)
        if self.phase_dim:
            yield "({0} << 2) /* phase */".format(self.phase_dim)
        if self.slice_dim:
            yield "({0} << 4) /* slice */".format(self.slice_dim)
        if self.dim_info & ~0x7F:
            yield "0x{0:X} /* invalid bits */"

    @property
    def readable_dim_info(self):
        """Readable value for the dim_info field (in C notation)."""
        if self.dim_info == 0:
            return "0 /* unspecified */"
        else:
            return " | ".join(self.separate_dim_info)

    @property
    def readable_intent_code(self):
        """Readable flag for the intent_code field."""
        try:
            return INTENT_DICT[self.intent_code]
        except KeyError:
            return "{0} /* unknown intent_code */".format(self.intent_code)


class Nifti1Header(NiftiHeader):
    """In-memory representation of a NIfTI-1 header.

    Can be constructed from a byte string containing a raw header. The
    fields can then be accessed as attributes.  Specialized attributes
    exist that produce readable, symbolic values for certain fields
    (as opposed to raw, bit-packed values).

    The header can be output to a stream using the print_raw() and
    print_interpreted() members.
    """

    FIELDS_DESCRIPTION = NIFTI1_FIELDS_DESCRIPTION
    STRUCT_FORMAT = NIFTI1_STRUCT_FORMAT
    TOTAL_HEADER_SIZE = 352

    @staticmethod
    def _test_byte_order(binary_header, byte_order):
        """Test if a binary header uses the given byte order ('<' or '>')."""
        dim0_binary = binary_header[40:42]
        dim0 = struct.unpack(str(byte_order + 'h'), dim0_binary)[0]
        return 1 <= dim0 <= 7

    def from_binary(self, binary_header):
        """Load from a raw binary header."""
        if len(binary_header) < 348:
            raise NotNiftiError(
                "file too short ({0} bytes, header is at least 348 bytes)"
                .format(len(binary_header)))

        nifti_version = self._test_magic_string_version(binary_header[344:348])
        if nifti_version is None:
            raise NotNiftiError("missing NIfTI-1 magic string")
        elif nifti_version != 1:
            raise NotNiftiError("unsupported NIfTI version {0}"
                                .format(nifti_version))

        self.byte_order = self._guess_byte_order(binary_header)

        unpack = struct.unpack(str(self.byte_order + self.STRUCT_FORMAT),
                               binary_header[:348])
        assert len(unpack) == 66

        self.raw['sizeof_hdr'] = unpack[0]
        self.raw['data_type'] = unpack[1].rstrip(b'\0')
        self.raw['db_name'] = unpack[2].rstrip(b'\0')
        self.raw['extents'] = unpack[3]
        self.raw['session_error'] = unpack[4]
        self.raw['regular'] = unpack[5]

        self.raw['dim_info'] = unpack[6]
        self.raw['dim'] = unpack[7:15]
        assert 1 <= self.dim[0] <= 7  # Guaranteed by byte order check

        self.raw['intent_p1'] = unpack[15]
        self.raw['intent_p2'] = unpack[16]
        self.raw['intent_p3'] = unpack[17]
        self.raw['intent_code'] = unpack[18]

        self.raw['datatype'] = unpack[19]
        self.raw['bitpix'] = unpack[20]

        self.raw['slice_start'] = unpack[21]
        self.raw['pixdim'] = unpack[22:30]

        self.raw['vox_offset'] = unpack[30]

        self.raw['scl_slope'] = unpack[31]
        self.raw['scl_inter'] = unpack[32]
        self.raw['slice_end'] = unpack[33]
        self.raw['slice_code'] = unpack[34]
        self.raw['xyzt_units'] = unpack[35]
        self.raw['cal_max'] = unpack[36]
        self.raw['cal_min'] = unpack[37]
        self.raw['slice_duration'] = unpack[38]
        self.raw['toffset'] = unpack[39]
        self.raw['glmax'] = unpack[40]
        self.raw['glmin'] = unpack[41]

        self.raw['descrip'] = unpack[42].rstrip('\0')
        self.raw['aux_file'] = unpack[43].rstrip('\0')

        self.raw['qform_code'] = unpack[44]
        self.raw['sform_code'] = unpack[45]

        self.raw['quatern_b'] = unpack[46]
        self.raw['quatern_c'] = unpack[47]
        self.raw['quatern_d'] = unpack[48]
        self.raw['qoffset_x'] = unpack[49]
        self.raw['qoffset_y'] = unpack[50]
        self.raw['qoffset_z'] = unpack[51]

        self.raw['srow_x'] = unpack[52:56]
        self.raw['srow_y'] = unpack[56:60]
        self.raw['srow_z'] = unpack[60:64]

        self.raw['intent_name'] = unpack[64].rstrip('\0')

        self.raw['magic'] = unpack[65]
        self.onefile = (self.magic[1:2] == b'+')

        binary_extender = binary_header[348:352]
        self.extensions_present = self.read_binary_extender(binary_extender)

    def check_format(self):
        """Check that the header satisfies the basic requirements of NIfTI-2"""
        if self.magic not in (b'ni1\0', b'n+1\0'):
            raise NotNiftiError("invalid magic string {0}".format(self.magic))
        if self.sizeof_hdr != 348:
            raise NotNiftiError("sizeof_hdr must be 348")

    def __repr__(self):
        return "Nifti1Header({0!r})".format(self.raw)


class Nifti2Header(NiftiHeader):
    """In-memory representation of a NIfTI-2 header.

    Can be constructed from a byte string containing a raw header. The
    fields can then be accessed as attributes.  Specialized attributes
    exist that produce readable, symbolic values for certain fields
    (as opposed to raw, bit-packed values).

    The header can be output to a stream using the print_raw() and
    print_interpreted() members.
    """

    FIELDS_DESCRIPTION = NIFTI2_FIELDS_DESCRIPTION
    STRUCT_FORMAT = NIFTI2_STRUCT_FORMAT
    TOTAL_HEADER_SIZE = 544

    @staticmethod
    def _test_byte_order(binary_header, byte_order):
        """Test if a binary header uses the given byte order ('<' or '>')."""
        dim0_binary = binary_header[16:24]
        dim0 = struct.unpack(str(byte_order + 'q'), dim0_binary)[0]
        return 1 <= dim0 <= 7

    def from_binary(self, binary_header):
        """Load from a raw binary header."""
        if len(binary_header) < 540:
            raise NotNiftiError(
                "file too short ({0} bytes, header is at least 540 bytes)"
                .format(len(binary_header)))

        nifti_version = self._test_magic_string_version(binary_header[4:12])
        if nifti_version is None:
            raise NotNiftiError("missing NIfTI-2 magic string")
        elif nifti_version != 2:
            raise NotNiftiError("unsupported NIfTI version {0}"
                                .format(nifti_version))

        self.byte_order = self._guess_byte_order(binary_header)

        unpack = struct.unpack(str(self.byte_order + self.STRUCT_FORMAT),
                               binary_header[:540])
        assert len(unpack) == 60

        self.raw['sizeof_hdr'] = unpack[0]
        self.raw['magic'] = unpack[1]

        self.raw['datatype'] = unpack[2]
        self.raw['bitpix'] = unpack[3]
        self.raw['dim'] = unpack[4:12]
        self.raw['intent_p1'] = unpack[12]
        self.raw['intent_p2'] = unpack[13]
        self.raw['intent_p3'] = unpack[14]
        self.raw['pixdim'] = unpack[15:23]
        self.raw['vox_offset'] = unpack[23]
        self.raw['scl_slope'] = unpack[24]
        self.raw['scl_inter'] = unpack[25]
        self.raw['cal_max'] = unpack[26]
        self.raw['cal_min'] = unpack[27]
        self.raw['slice_duration'] = unpack[28]
        self.raw['toffset'] = unpack[29]
        self.raw['slice_start'] = unpack[30]
        self.raw['slice_end'] = unpack[31]
        self.raw['descrip'] = unpack[32].rstrip(b'\0')
        self.raw['aux_file'] = unpack[33].rstrip(b'\0')
        self.raw['qform_code'] = unpack[34]
        self.raw['sform_code'] = unpack[35]
        self.raw['quatern_b'] = unpack[36]
        self.raw['quatern_c'] = unpack[37]
        self.raw['quatern_d'] = unpack[38]
        self.raw['qoffset_x'] = unpack[39]
        self.raw['qoffset_y'] = unpack[40]
        self.raw['qoffset_z'] = unpack[41]
        self.raw['srow_x'] = unpack[42:46]
        self.raw['srow_y'] = unpack[46:50]
        self.raw['srow_z'] = unpack[50:54]
        self.raw['slice_code'] = unpack[54]
        self.raw['xyzt_units'] = unpack[55]
        self.raw['intent_code'] = unpack[56]
        self.raw['intent_name'] = unpack[57].rstrip(b'\0')
        self.raw['dim_info'] = unpack[58]
        self.raw['unused_str'] = unpack[59].rstrip(b'\0')

        self.onefile = (self.magic[1:2] == b'+')

        binary_extender = binary_header[540:544]
        self.extensions_present = self.read_binary_extender(binary_extender)

    def check_format(self):
        """Check that the header satisfies the basic requirements of NIfTI-2"""
        if self.magic not in (b'ni2\0\r\n\032\n', b'n+2\0\r\n\032\n'):
            raise NotNiftiError("invalid magic string {0}".format(self.magic))
        if self.sizeof_hdr != 540:
            raise NotNiftiError("sizeof_hdr must be 540")

    def list_inconsistencies(self):
        """List the inconsistencies of the header's raw data."""
        self.check_format()
        for i in super(Nifti2Header, self).list_inconsistencies():
            yield i

        if self.unused_str != '':
            yield InconsistentNiftiError(
                "unused_str should be set to all zero bytes")

    def __repr__(self):
        return "Nifti2Header({0!r})".format(self.raw)


def print_affine_matrix(R, file=sys.stdout):
    """Print an affine matrix to the given file."""
    # TODO handle double precision of NIfTI-2 matrices
    print("{0:8.5f} {1:8.5f} {2:8.5f} {3:8.3f}".format(*R[0]), file=file)
    print("{0:8.5f} {1:8.5f} {2:8.5f} {3:8.3f}".format(*R[1]), file=file)
    print("{0:8.5f} {1:8.5f} {2:8.5f} {3:8.3f}".format(*R[2]), file=file)


def read_nifti_header(binary_header):
    sizeof_hdr_le = struct.unpack(str("<i"), binary_header[:4])[0]
    sizeof_hdr_be = struct.unpack(str(">i"), binary_header[:4])[0]
    if sizeof_hdr_le == 348 or sizeof_hdr_be == 348:
        if binary_header[344:348] in (b'ni1\0', b'n+1\0'):
            return Nifti1Header(binary_header)
    elif sizeof_hdr_le == 540 or sizeof_hdr_be == 540:
        if binary_header[4:12] in (b'ni2\0\r\n\032\n', b'n+2\0\r\n\032\n'):
            return Nifti2Header(binary_header)


def main():
    """The script's entry point."""
    # Process command line
    from optparse import OptionParser
    parser = OptionParser(description=__doc__,
                          usage="%prog [-nrv] [--] nifti_file",
                          version="%prog unreleased")
    parser.add_option("-n", "--no-interpreted", action="store_false",
                      dest="interpreted",
                      help="do not print the interpreted header")
    parser.add_option("-r", "--raw", action="store_true",
                      dest="raw", help="print raw header fields")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose",
                      help="print the description of raw header fields")

    parser.set_defaults(interpreted=True, raw=False, verbose=False)

    (options, args) = parser.parse_args()

    if len(args) == 1:
        filename = args[0]
    elif len(args) > 1:
        parser.error("one file at a time, please")
    else:
        parser.error("no input specified")

    # Read NIfTI file header
    if os.path.splitext(filename)[1] == '.gz':
        import gzip
        import contextlib
        try:
            with contextlib.closing(gzip.open(filename, 'rb')) as f:
                binary_header = f.read(544)
        except IOError as exc:
            sys.exit("error reading {0!r} as a gzip-compressed file: {1}"
                     .format(filename, exc))
    elif filename == "-":
        try:
            binary_header = sys.stdin.read(544)
            sys.stdin.close()
        except IOError as exc:
            sys.exit("error reading standard input: {0}".format(exc))
    else:
        import io
        try:
            with io.open(filename, 'rb', buffering=0) as f:
                binary_header = f.read(544)
        except IOError as exc:
            sys.exit("error reading {0!r}: {1}".format(filename, exc))

    # Parse header
    header = read_nifti_header(binary_header)
    if header is None:
        sys.exit("error: {0} is not a NIfTI-1 or NIfTI-2 file"
                 .format(filename))

    # Print output to stdout, ensuring proper failure if printing fails
    try:
        if options.raw:
            header.print_raw(describe_fields=options.verbose)
        if options.interpreted:
            header.print_interpreted()

        inconsistencies = list(header.list_inconsistencies())
        if inconsistencies:
            print()
            print("INCONSISTENCIES")
            print("===============")
            for i in inconsistencies:
                print("- {0}".format(i))

        sys.stdout.close()
    except IOError as exc:
        sys.exit("error printing output: {0}".format(exc))


if __name__ == '__main__':
    main()
