#!/usr/bin/env python2
#
# Intel HEX to ST microelectronics DfuSe file converter
# Copyright (C)2015 Thomas Kindler <mail@t-kindler.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import os
import struct
import argparse
from zlib import crc32
from intelhex import IntelHex


def load_hex():
    """
    Load hex or binary file.
    :return:    intelhex object
    """
    if args.verbose:
        print "Loading %s..." % args.source

    try:
        ih = IntelHex()

        if args.format == "hex":
            ih.loadhex(args.source)
        else:
            ih.loadbin(args.source, args.start_addr)

    except Exception, e:
        print e
        exit(1)

    if args.verbose:
        print "  Start: 0x%08x" % ih.minaddr()
        print "  End  : 0x%08x" % ih.maxaddr()

    return ih


def save_dfu(ih):
    """
    Save as STMicroelectronics DfuSe file.
    see UM0391 - DfuSe File Format Specification

    :param ih:  intelhex object
    """
    if args.verbose:
        print "Saving %s..." % args.target
        print "  Device ID: 0x%04x:0x%04x" % (args.vid, args.pid)
        print "  Target name: %s" % args.target_name

    # Image element
    #
    image_data = ih.tobinstr()

    data = struct.pack(
        "<II",
        ih.minaddr(),       # dwElementAddress
        len(image_data)     # dwElementSize
    ) + image_data          # Data

    # Target prefix
    #
    szTargetName = args.target_name.encode("ascii")

    data = struct.pack(
        "<6sBI255sII",
        b"Target",          # szSignature
        0,                  # bAlternateSetting
        1,                  # bTargetNamed
        szTargetName,       # szTargetName
        len(data),          # dwTargetSize
        1                   # dwNbElements
    ) + data

    # Prefix
    #
    data = struct.pack(
        "<5sBIB",
        b"DfuSe",           # szSignature
        0x01,               # bVersion,
        len(data)+11,       # DFUImageSize,
        1                   # bTargets
    ) + data

    # Suffix
    #
    data += struct.pack(
        "<HHHH3sB",
        0xFFFF,         # bcdDevice
        args.pid,       # idProduct
        args.vid,       # idVendor
        0x011a,         # bdcDFU
        b"UFD",         # ucDfuSignature
        16              # bLength
    )

    dwCRC = ~crc32(data) & 0xFFFFFFFF

    data += struct.pack(
        "<I",
        dwCRC           # dwCRC
    )

    try:
        open(args.target, "wb").write(data)

    except Exception, e:
        print e
        exit(1)


# Parse arguments
#
DEFAULT_DEVICE = "0x0483:0xdf11"

parser = argparse.ArgumentParser(
    description="Convert hex files to STMicroelectronics DfuSe format"
)

parser.add_argument(
    "--version",
    action="version",
    version="%(prog)s 1.1.0"
)

parser.add_argument(
    "-q", "--quiet", dest="verbose",
    default=True,
    action="store_false",
    help="do not print status messages"
)

parser.add_argument(
    "source", help="source file"
)

parser.add_argument(
    "target", nargs = "?",
    help="target file"
)

parser.add_argument(
    "-f", "--format",
    default="hex", choices=["hex", "bin"],
    help="source file format (default: %(default)s)"
)

parser.add_argument(
    "-s", "--start", dest="start_addr",
    help="set start address (for bin files)"
)

parser.add_argument(
    "-d", "--device",
    default=DEFAULT_DEVICE,
    help="device VID:PID (default: %(default)s)"
)

parser.add_argument(
    "-n", "--name", dest="target_name",
    default="application",
    help="target name (default: %(default)s)"
)

args=parser.parse_args()

# Check arguments
#
if args.target == None:
    (root, ext) = os.path.splitext(args.source)
    args.target = root + ".dfu"


if args.format == "bin":
    if args.start_addr == None:
        print "option --start required for binary files"
        exit(1)

    args.start_addr = int(args.start_addr, 0)
else:
    if args.start_addr != None:
        print "option --start not allowed for hex files"
        exit(1)


args.vid = int(args.device.split(':', 1)[0], 0)
args.pid = int(args.device.split(':', 1)[1], 0)


# Convert file
#
save_dfu( load_hex() )
