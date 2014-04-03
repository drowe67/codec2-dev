#Asterisk 11 Codec 2 support
===========================

##Description
These patches add Codec 2 2400 support to Asterisk 11.
The following patches are provided:

* asterisk-11.8.1-codec2.patch: plain Asterisk 11.      
* asterisk-11.8.1-opus-codec2.patch: Asterisk patched with Meetecho's Opus codec support.

##Building
Building and installing are integrated within Asterisk building environment. libcodec2 must be installed beforehand.

##Credits
I've followed the example of [asterisk-opus](https://github.com/meetecho/asterisk-opus), by [@meetecho](https://github.com/meetecho), to adapt Codec2 Asterisk 1.8 patch to version 11.

Many thanks to the [Codec2](http://www.rowetel.com/blog/codec2.html) team for developing such great codec!

Developed by [Antonio Eugenio Burriel](https://github.com/aeburriel)
