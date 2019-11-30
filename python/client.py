#!/usr/bin/env python2

from __future__ import print_function
from time import sleep
import sys
import zmq

# port = 5051
port = 5052



# Socket to talk to server
context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.setsockopt(zmq.SUBSCRIBE, "")


print ("connect...")
socket.connect ("tcp://192.168.0.235:%s" % port)
# socket.connect ("tcp://127.0.0.1:%s" % port)

while True:
    m = socket.recv()
    print (m)