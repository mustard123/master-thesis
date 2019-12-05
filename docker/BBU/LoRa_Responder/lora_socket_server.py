#!/usr/bin/env python2

from __future__ import print_function
import argparse
import time
import zmq
import socket
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='Connect to udp port for receiving decoded LoRa signals,'
' if an ACK is required publish ACK iq samples via zmq socket for Remote Radio Head to receive and send out (TX).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-o', '--out-port', default=5052, help='zmq port to publish downstream (i.e ACK) iq samples')
parser.add_argument('-i', '--input-port', default=40868, help='UDP port to connect for receiving decoded lora messages')
args = parser.parse_args()

context = zmq.Context()
zmq_socket = context.socket(zmq.PUB)
port = args.out_port
udp_port = args.input_port
zmq_socket.bind("tcp://*:%s" % port)
print ('zmq publishes on all interfaces on port ' + str(port))


with open (dir_path + "/gateway_downlink_trimmed_23_bytes_SF12_CR4.raw") as f:
    ack = f.read()

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.bind(("127.0.0.1", udp_port))
print ('listen for decoded lora packages on udp port ' + str(udp_port))
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
