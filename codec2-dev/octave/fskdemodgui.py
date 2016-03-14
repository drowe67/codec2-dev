#!/usr/bin/env python
#
#	fsk_demod Statistics GUI
#	Accepts the stats output from fsk_demod on stdin, and plots it.
#
#	Mark Jessop 2016-03-13 <vk5qi@rfhead.net>
#
#	NOTE: This is intended to be run on a 'live' stream of samples, and hence expects
#	updates at about 10Hz. Anything faster will fill up the input queue and be discarded.
#
#	Call using: 
#	<producer>| ./fsk_demod 2X 8 923096 115387 - - S 2> >(python ~/Dev/codec2-dev/python/fskdemodgui.py) | <consumer>
#
#
import sys, time, json, Queue
from threading import Thread
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

# Some settings...
update_rate = 10 # Hz
history_size = 10*10 # 10 seconds at 10Hz...
history_scale = np.linspace((-1*history_size+1)/10.0,0,history_size)

# Input queue
in_queue = Queue.Queue(history_size) # 1-element FIFO... 

win = pg.GraphicsWindow()
win.setWindowTitle('FSK Demodulator Modem Statistics')

# Plot objects
ebno_plot = win.addPlot(title="Eb/No")
ebno_plot.setLabel('left','Eb/No (dB)')
ebno_plot.setLabel('bottom','Time (seconds)')
ppm_plot = win.addPlot(title="Symbol Rate Offset")
ppm_plot.setLabel('left','Offset (ppm)')
ppm_plot.setLabel('bottom','Time (seconds)')
win.nextRow()
fest_plot = win.addPlot(title="Tone Frequency Estimation")
fest_plot.setLabel('left','Frequency (Hz)')
fest_plot.setLabel('bottom','Time (seconds)')
eye_plot = win.addPlot(title="Eye Diagram")

# Data arrays...
ebno_data = np.zeros(history_size)
ppm_data = np.zeros(history_size)
fest_data = np.zeros((4,history_size))

# Curve objects, so we can update them...
ebno_curve = ebno_plot.plot(x=history_scale,y=ebno_data)
ppm_curve = ppm_plot.plot(x=history_scale,y=ppm_data)
fest1_curve = fest_plot.plot(x=history_scale,y=fest_data[0,:],pen=(255,0,0)) # f1 = Red
fest2_curve = fest_plot.plot(x=history_scale,y=fest_data[1,:],pen=(0,255,0)) # f2 = Blue

# Plot update function. Reads from queue, processes and updates plots.
def update_plots():
	global timeout,timeout_counter,eye_plot,ebno_curve, ppm_curve, fest1_curve, fest2_curve, ebno_data, ppm_data, fest_data, in_queue

	try:
		if in_queue.empty():
			return
		in_data = in_queue.get_nowait()
		in_data = json.loads(in_data)
	except Exception as e:

		sys.stderr.write(str(e))
		return

	try:
		new_ebno = in_data['EbNodB']
		new_ppm = in_data['ppm']
		new_fest1 = in_data['f1_est']
		new_fest2 = in_data['f2_est']
	except Exception as e:
		print("ERROR reading dict: %s" % e)

	try:
		# Roll data arrays
		ebno_data[:-1] = ebno_data[1:]
		ppm_data[:-1] = ppm_data[1:]
		fest_data = np.roll(fest_data,-1,axis=1)

		# Add in new data points
		ebno_data[-1] = new_ebno
		ppm_data[-1] = new_ppm
		fest_data[0,-1] = new_fest1
		fest_data[1,-1] = new_fest2

		# Update plots
		ebno_curve.setData(x=history_scale,y=ebno_data)
		ppm_curve.setData(x=history_scale,y=ppm_data)
		fest1_curve.setData(x=history_scale,y=fest_data[0,:],pen=(255,0,0)) # f1 = Red
		fest2_curve.setData(x=history_scale,y=fest_data[1,:],pen=(0,255,0)) # f2 = Blue

	except Exception as e:
		print(e.strerror)

	# Now try reading in the eye diagram and other tone estimates.
	try:
		eye_data = np.array(in_data['eye_diagram'])
		eye_plot.clear()
		for line in eye_data:
			eye_plot.plot(line)

	except Exception as e:
		pass

	

timer = pg.QtCore.QTimer()
timer.timeout.connect(update_plots)
timer.start(30)


# Thread to read from stdin and push into a queue to be processed.
def read_input():
	global in_queue

	while True:
		in_line = sys.stdin.readline()

		# Only push actual data into the queue...
		# This stops sending heaps of empty strings into the queue when fsk_demod closes.
		if in_line == "":
			time.sleep(0.1)
			continue

		if not in_queue.full():
			in_queue.put_nowait(in_line)

read_thread = Thread(target=read_input)
read_thread.daemon = True # Set as daemon, so when all other threads die, this one gets killed too.
read_thread.start()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
	import sys
	if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
		try:
			QtGui.QApplication.instance().exec_()
		except KeyboardInterrupt:
			sys.exit(0)
