#!/usr/bin/env python
#
#   Horus Binary (and oldschool) Telemetry Uploader
#
#   Mark Jessop 2015-12-31
#   <vk5qi@rfhead.net>
#
#   This script takes either a hex representation of the binary payload, or 
#   a 'classic' ASCII sentence, and uploads it to Habitat.
#
#   Currently this script tells the two apart by looking for 'HORUS' at the start
#   of the argument to determine if it's an ASCII sentence.
#   
#   It's designed to be called from fsk_horus_stream.m, and is tailored for it's output.
#
#   Dependencies:
#   - Python 2.7 (Will probably break in Python 3)
#   - crcmod (pip install crcmod)
#

import time, struct, json, socket, httplib, crcmod, argparse, sys
from base64 import b64encode
from hashlib import sha256
from datetime import datetime

def crc16_ccitt(data):
    """
    Calculate the CRC16 CCITT checksum of *data*.
    
    (CRC16 CCITT: start 0xFFFF, poly 0x1021)
    """
    crc16 = crcmod.predefined.mkCrcFun('crc-ccitt-false')
    return crc16(data)

# Binary packet format, from https://github.com/darksidelemm/PicoHorusBinary/tree/master/PicoPayloadGPS
# struct TBinaryPacket
# {
# uint8_t       PayloadID;
# uint16_t  Counter;
# uint8_t       Hours;
# uint8_t       Minutes;
# uint8_t       Seconds;
# float     Latitude;
# float     Longitude;
# uint16_t      Altitude;
# uint8_t   Speed; // Speed in Knots (1-255 knots)
# uint8_t   Sats;
# int8_t   Temp; // Twos Complement Temp value.
# uint8_t   BattVoltage; // 0 = 0.5v, 255 = 2.0V, linear steps in-between.
# uint16_t Checksum; // CRC16-CCITT Checksum.
# };  //  __attribute__ ((packed));

def decode_horus_binary_telemetry(payload):

    horus_format_struct = "<BHBBBffHBBbBH"
    try:
        unpacked = struct.unpack(horus_format_struct, payload)
    except:
        print "Wrong string length. Packet contents:"
        print ":".join("{:02x}".format(ord(c)) for c in payload)
        sys.exit(1)

    telemetry = {}
    telemetry['payload_id'] = unpacked[0]
    telemetry['counter'] = unpacked[1]
    telemetry['time'] = "%02d:%02d:%02d" % (unpacked[2],unpacked[3],unpacked[4])
    telemetry['latitude'] = unpacked[5]
    telemetry['longitude'] = unpacked[6]
    telemetry['altitude'] = unpacked[7]
    telemetry['speed'] = unpacked[8]
    telemetry['sats'] = unpacked[9]
    telemetry['temp'] = unpacked[10]
    telemetry['batt_voltage_raw'] = unpacked[11]
    telemetry['checksum'] = unpacked[12]

    # Convert some of the fields into more useful units.
    telemetry['batt_voltage'] = 0.5 + 1.5*telemetry['batt_voltage_raw']/255.0

    return telemetry

# Compatible with Habitat Payload ID 55f39c02e518bdd2885c8aac7d1cdd7c
def telemetry_to_sentence(telemetry):
    sentence = "$$PICOHORUSBINARY,%d,%s,%.5f,%.5f,%d,%d,%d,%d,%.2f" % (telemetry['counter'],telemetry['time'],telemetry['latitude'],
        telemetry['longitude'],telemetry['altitude'],telemetry['speed'],telemetry['sats'],telemetry['temp'],telemetry['batt_voltage'])

    checksum = hex(crc16_ccitt(sentence[2:]))[2:].upper().zfill(4)
    output = sentence + "*" + checksum + "\n"
    return output

# Habitat Upload Functions
def habitat_upload_sentence(sentence, callsign="N0CALL"):

    sentence_b64 = b64encode(sentence)

    date = datetime.utcnow().isoformat("T") + "Z"

    data = {
        "type": "payload_telemetry",
        "data": {
            "_raw": sentence_b64
            },
        "receivers": {
            callsign: {
                "time_created": date,
                "time_uploaded": date,
                },
            },
    }
    try:
        c = httplib.HTTPConnection("habitat.habhub.org",timeout=4)
        c.request(
            "PUT",
            "/habitat/_design/payload_telemetry/_update/add_listener/%s" % sha256(sentence_b64).hexdigest(),
            json.dumps(data),  # BODY
            {"Content-Type": "application/json"}  # HEADERS
            )

        response = c.getresponse()
        sys.exit(0)
    except Exception as e:
        print("Failed to upload to Habitat.")
        sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("raw_data", help="Raw Data, either hex binary data, or ASCII payload string.")
parser.add_argument("-c","--callsign",default="N0CALL",help="Habitat Upload Callsign")
args = parser.parse_args()


uploader_callsign = args.callsign
raw_data = args.raw_data
print(raw_data)

if raw_data.startswith("AX5ARG"):
    # Assume the data is a standard telemetry string and just upload it.
    # Append a Newline and checksum (if not alread there) to the end, and "$$"'s to the start.
    if not '*' in raw_data:
        # Assume there is no checksum on the end of the string, and add one.
        checksum = hex(crc16_ccitt(raw_data))[2:].upper().zfill(4)
        raw_data = raw_data + "*" + checksum
    if raw_data[:-1] != '\n':
        raw_data = "$$" + raw_data + '\n'

    habitat_upload_sentence(raw_data, callsign = uploader_callsign)
else:
    # Attempt to decode some hex data.
    data = raw_data.decode("hex")
    telem = decode_horus_binary_telemetry(data)
    # Only convert and upload if checksum passes.
    if(crc16_ccitt(data[:-2]) != telem['checksum']):
        print("Checksum Failed!")
        sys.exit(1)

    sentence = telemetry_to_sentence(telem)
    print("Uploading: %s"%(sentence))
    habitat_upload_sentence(sentence, callsign = uploader_callsign)

