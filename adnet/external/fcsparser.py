"""
Modified from https://github.com/jfrelinger/fcm
"""

"""Provide access to data stored in FCS files"""

from warnings import warn
from operator import and_
from math import log
from struct import calcsize, unpack
from io import IOBase
import numpy
import os
from functools import reduce
import six


class FCSreader(object):

    """
    class to hold object to read and parse fcs files.  main usage is to get
    a FCMdata object out of a fcs file
    """

    def __init__(self, filename, transform=None, sidx=None, spill=None):
        # self.filename = filename
        # if type(filename) == str:
        #    self.filename = filename
        #    self._fh = open(filename, 'rb')
        # else: # we should have gotten a filehandle then
        #    self.filename = filename.name
        #    self._fh = filename
        # self._fh = cStringIO.StringIO(open(filename, 'rb').read())
        # self.filename = filename

        if isinstance(filename, six.string_types):
            try:
                self._fh = open(filename, 'rb')
            except IOError:
                raise
            except:
                raise Exception('Unknown exception occurred reading file')
        elif isinstance(filename, (file, IOBase)):
            self._fh = filename
        else:
            raise TypeError(
                "Filename must be a file path or a file handle (either 'file' type or io.IOBase")
        if transform not in ['log', 'logicle', None]:
            raise TypeError(
                'transform argument must be one of [\'log\', \'logicle\', None]')
        self.transform = transform
        self.cur_offset = 0
        self.spill = spill
        self.sidx = sidx

    def get_FCMdata(self, auto_comp=False, **kwargs):
        """Return the next FCM data set stored in a FCS file"""

        # parse headers
        header = self.parse_header(self.cur_offset)

        # parse text
        text = self.parse_text(
            self.cur_offset,
            header['text_start'],
            header['text_stop'])

        # parse annalysis
        try:
            astart = text['beginanalysis']
        except KeyError:
            astart = header['analysis_start']
        try:
            astop = text['endanalysis']
        except KeyError:
            astop = header['analysis_end']
        analysis = self.parse_analysis(self.cur_offset, astart, astop)
        # parse data
        try:
            dstart = int(text['begindata'])
        except KeyError:
            dstart = header['data_start']
        try:
            dstop = int(text['enddata'])
        except KeyError:
            dstop = header['data_end']

        # account for LMD reporting the wrong values for the size of the data
        # segment
        lmd = self.fix_lmd(
            self.cur_offset,
            header['text_start'],
            header['text_stop'])
        dstop = dstop + lmd
        data = self.parse_data(self.cur_offset, dstart, dstop, text)

        # build fcmdata object
        channels = []
        scchannels = []
        scchannel_indexes = []
        to_transform = []
        base_chan_name = []
        for i in range(1, int(text['par']) + 1):
            base_chan_name.append(text['p%dn' % i])
            try:
                if text['p%ds' % i] not in ['', ' ']:
                    name = text['p%ds' % i]
                else:
                    name = text['p%dn' % i]
            except KeyError:
                name = text['p%dn' % i]
            channels.append(name)
            # if not name.lower().startswith('fl'):
            if not is_fl_channel(name):
                scchannels.append(name)
                if name != 'Time':
                    scchannel_indexes.append(i - 1)
            else:  # we're a FL channel
                try:
                    if text['p%dr' % i] == '262144':
                        to_transform.append(i - 1)
                except KeyError:
                    pass

        try:
            unused_path, name = os.path.split(self._fh.name)
        except AttributeError:
            name = 'InMemoryFile'
        name, unused_ext = os.path.splitext(name)
        return name, channels, data


    def read_bytes(self, offset, start, stop):
        """Read in bytes from start to stop inclusive."""

        self._fh.seek(offset + start)

        return self._fh.read(stop - start + 1)

    def parse_header(self, offset):
        """
        Parse the FCM data in fcs file at the offset (supporting multiple
        data segments in a file
        """

        header = {}
        header['version'] = float(self.read_bytes(offset, 3, 5))
        header['text_start'] = int(self.read_bytes(offset, 10, 17))
        header['text_stop'] = int(self.read_bytes(offset, 18, 25))
        header['data_start'] = int(self.read_bytes(offset, 26, 33))
        header['data_end'] = int(self.read_bytes(offset, 34, 41))
        try:
            header['analysis_start'] = int(self.read_bytes(offset, 42, 49))
        except ValueError:
            header['analysis_start'] = -1
        try:
            header['analysis_end'] = int(self.read_bytes(offset, 50, 57))
        except ValueError:
            header['analysis_end'] = -1

        return header

    def fix_lmd(self, offset, start, stop):
        """function to handle lmd counting differently then most other FCS data"""
        text = self.read_bytes(offset, start, stop)
        if text[0] == text[-1]:
            return 0
        else:
            return -1

    def parse_text(self, offset, start, stop):
        """return parsed text segment of fcs file"""

        text = self.read_bytes(offset, start, stop)
        # TODO: add support for suplement text segment
        return parse_pairs(text)

    def parse_analysis(self, offset, start, stop):
        """return parsed analysis segment of fcs file"""

        if start == stop:
            return {}
        else:
            text = self.read_bytes(offset, start, stop)
            return parse_pairs(text)

    def parse_data(self, offset, start, stop, text):
        """return numpy.array of data segment of fcs file"""

        dtype = text['datatype']
        mode = text['mode']
        tot = int(text['tot'])
        if mode == 'c' or mode == 'u':
            raise('not supported')

        if text['byteord'] == '1,2,3,4' or text['byteord'] == '1,2':
            order = '<'
        elif text['byteord'] == '4,3,2,1' or text['byteord'] == '2,1':
            order = '>'
        else:
            warn(
                "unsupported byte order %s , using default @" %
                text['byteord'])
            order = '@'
        # from here on out we assume mode l (list)

        bitwidth = []
        drange = []
        for i in range(1, int(text['par']) + 1):
            bitwidth.append(int(text['p%db' % i]))
            drange.append(int(float(text['p%dr' % i])))

        if dtype.lower() == 'i':
            data = self.parse_int_data(
                offset,
                start,
                stop,
                bitwidth,
                drange,
                tot,
                order)
        elif dtype.lower() == 'f' or dtype.lower() == 'd':
            data = self.parse_float_data(
                offset,
                start,
                stop,
                dtype.lower(),
                tot,
                order)
        else:  # ascii
            data = self.parse_ascii_data(
                offset,
                start,
                stop,
                bitwidth,
                dtype,
                tot,
                order)
        return data

    def parse_int_data(
            self,
            offset,
            start,
            stop,
            bitwidth,
            drange,
            tot,
            order):
        """Parse out and return integer list data from fcs file"""

        if reduce(and_, [item in [8, 16, 32] for item in bitwidth]):
            if len(set(bitwidth)) == 1:  # uniform size for all parameters
                # calculate how much data to read in.
                num_items = (stop - start + 1) / \
                    calcsize(fmt_integer(bitwidth[0]))
                try:
                    fmt = numpy.dtype(
                        '%s%d%s' %
                        (order,
                         num_items,
                         fmt_integer(
                             bitwidth[0])))
                    tmp = numpy.memmap(
                        self._fh,
                        fmt,
                        mode='c',
                        offset=start +
                        offset,
                        shape=1)
                except IOError:
                    # unpack into a list
                    tmp = unpack(
                        '%s%d%s' %
                        (order, num_items, fmt_integer(
                            bitwidth[0])), self.read_bytes(
                            offset, start, stop))
                    tmp = numpy.array(tmp)

            # parameter sizes are different e.g. 8, 8, 16,8, 32 ... do one at a
            # time
            else:
                unused_bitwidths = map(int, map(log2, drange))
                tmp = []
                cur = start
                while cur < stop:
                    for i, curwidth in enumerate(bitwidth):
                        bitmask = mask_integer(curwidth, unused_bitwidths[i])
                        nbytes = curwidth / 8
                        bin_string = self.read_bytes(
                            offset,
                            cur,
                            cur +
                            nbytes -
                            1)
                        cur += nbytes
                        val = bitmask & unpack(
                            '%s%s' %
                            (order, fmt_integer(curwidth)), bin_string)[0]
                        tmp.append(val)
        else:  # non starndard bitwiths...  Does this happen?
            warn('Non-standard bitwidths for data segments')
            return None
        return numpy.array(tmp).reshape((tot, len(bitwidth)))

    def parse_float_data(self, offset, start, stop, dtype, tot, order):
        """Parse out and return float list data from fcs file"""

        # count up how many to read in
        num_items = (stop - start + 1) / calcsize(dtype)
        try:
            fmt = numpy.dtype('%s%d%s' % (order, num_items, dtype))
            tmp = numpy.memmap(
                self._fh,
                fmt,
                mode='c',
                offset=start +
                offset,
                shape=1)
        except IOError:
            tmp = unpack(
                '%s%d%s' %
                (order, num_items, dtype), self.read_bytes(
                    offset, start, stop))
            tmp = numpy.array(tmp)
        return tmp.reshape((tot, num_items / tot))

    def parse_ascii_data(
            self,
            offset,
            start,
            stop,
            bitwidth,
            dtype,
            tot,
            order):
        """Parse out ascii encoded data from fcs file"""

        num_items = (stop - start + 1) / calcsize(dtype)
        # tmp = unpack('%s%d%s' % (order, num_items, dtype), self.read_bytes(offset, start, stop))
        try:
            fmt = numpy.dtype('%s%d%s' % (order, num_items, dtype))
            tmp = numpy.memmap(
                self._fh,
                fmt,
                mode='c',
                offset=start +
                offset,
                shape=1)
        except IOError:
            tmp = unpack(
                '%s%d%s' %
                (order, num_items, dtype), self.read_bytes(
                    offset, start, stop))
            tmp = numpy.array(tmp)
        return tmp.reshape((tot, num_items / tot))


def parse_pairs(text):
    """return key/value pairs from a delimited string"""
    delim = text[0]

    if delim != text[-1]:
        warn("text in segment does not start and end with delimiter")

    # no longer needed as we're not using regex
    # if delim == r'|':
    #    delim = '\|'
    # elif delim == r'\a'[0]:  # test for delimiter being \
    #    delim = '\\\\'  # regex will require it to be \\

    tmp = text[1:-1]
    # match the delimited character unless it's doubled
    # regex = re.compile('(?<=[^%s])%s(?!%s)' % (delim, delim, delim))
    # tmp = regex.split(tmp)
    # return dict(zip([x.lower().replace(delim +
    #                                   delim, delim) for x in tmp[::2]], [x.replace(delim +
    # delim, delim) for x in tmp[1::2]]))
    sp = tmp.split(delim)
    rslts = {}
    for i in range(len(sp) / 2):
        rslts[sp[(i * 2)].lower().replace('$', '')] = sp[((i * 2) + 1)]
    return rslts


def fmt_integer(b):
    """return binary format of an integer"""

    if b == 8:
        return 'B'
    elif b == 16:
        return 'H'
    elif b == 32:
        return 'I'
    else:
        print("Cannot handle integers of bit size %d" % b)
        return None


def mask_integer(b, ub):
    """return bitmask of an integer and a bitwitdh"""

    if b == 8:
        return (0xFF >> (b - ub))
    elif b == 16:
        return (0xFFFF >> (b - ub))
    elif b == 32:
        return (0xFFFFFFFF >> (b - ub))
    else:
        print("Cannot handle integers of bit size %d" % b)
        return None


def log_factory(base):
    """constructor of various log based functions"""

    def f(x):
        return log(x, base)
    return f

log2 = log_factory(2)


def loadFCS(
        filename,
        transform=None,
        auto_comp=False,
        spill=None,
        sidx=None,
        file_index=0,
        **kwargs):
    """Load and return a FCM data object from an FCS file"""

    tmp = FCSreader(filename, transform, spill=spill, sidx=sidx)
    for _ in range(file_index + 1):
        data = tmp.get_FCMdata(auto_comp, **kwargs)
    tmp._fh.close()
    del tmp
    return data


def loadMultipleFCS(
        files,
        transform=None,
        auto_comp=False,
        spill=None,
        sidx=None,
        **kwargs):
    for filename in files:
        tmp = loadFCS(filename, transform, auto_comp, spill, sidx, **kwargs)
        try:
            if 'r' not in kwargs.keys():
                kwargs['r'] = tmp._r
        except AttributeError:
            pass
        yield tmp


def is_fl_channel(name):
    """
    Try and decide if a channel is a flourescent channel or if it's some other type
    returns a boolean
    """
    name = name.lower()
    if name.startswith('cs'):
        return False
    elif name.startswith('fs'):
        return False
    elif name.startswith('ss'):
        return False
    elif name.startswith('ae'):
        return False
    elif name.startswith('cv'):
        return False
    elif name.startswith('time'):
        return False
    else:
        return True

if __name__ == '__main__':
    fn = '/home/vitoz/mnt/imls-bod/XiaokangL/EGF transfection/20150422/gated/BaseFileName_B2_cell_SRC-GFP EGF 0min.fcs'
    channels, data = loadFCS(fn)
    print(data.shape)
    print(len(channels))
    print(data.mean)


