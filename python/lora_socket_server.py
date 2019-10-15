#!/usr/bin/env python2

from __future__ import print_function
import argparse
import time
import zmq
import socket

parser = argparse.ArgumentParser(description='Connect to udp port for receiving decoded LoRa signals,'
' if an ACK is required publish ACK iq samples via zmq socket for Remote Radio Head to receive.')

parser.add_argument('-p', '--port', default=5052, help='zmq port to publish downstream (i.e ACK) iq samples')
args = parser.parse_args()

context = zmq.Context()
zmq_socket = context.socket(zmq.PUB)
port = args.port
zmq_socket.bind("tcp://*:%s" % port)
print ('zmq publishes on all interfaces on port ' + str(port))


with open ("/home/sili/Documents/signal_recordings/gateway_downlink_trimmed.raw") as f:
    ack = f.read()

upd_port = 40868
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.bind(("127.0.0.1", upd_port))
print ('listen for decoded lora packages on udp port ' + str(upd_port))
while True:
    data, addr = s.recvfrom(1024)

    print ("decoded data is: ", data)
    print ("as hex: ", data.encode("hex"))

    # TODO Better to check if ACK is last 3 bytes of actual payload, attention lora_receiver also passes header and crc to us
    if ("ACK" in data):
        print ("received package request ACK, send ACK to remote radio head")
        #TODO Acutally genearate appropriate down packet instead of sending recorded packet
        #zmq_socket.send("hii") does not work with gnuradio
        zmq_socket.send(ack) # works! gnuradio probably only works with passing complex types
        print ('Done')
    else:
        print ("received package requests no ACK")