#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Zero Mq Split B
# Generated: Thu Dec  5 18:43:23 2019
##################################################


from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import lora


class zero_mq_split_b(gr.top_block):

    def __init__(self, bandwidth=125000, capture_freq=868500000, decoded_out_port=40868, samp_rate=1000000, spreading_factor=12, zmq_address_iq_in='tcp://127.0.0.1:5051'):
        gr.top_block.__init__(self, "Zero Mq Split B")

        ##################################################
        # Parameters
        ##################################################
        self.bandwidth = bandwidth
        self.capture_freq = capture_freq
        self.decoded_out_port = decoded_out_port
        self.samp_rate = samp_rate
        self.spreading_factor = spreading_factor
        self.zmq_address_iq_in = zmq_address_iq_in

        ##################################################
        # Variables
        ##################################################
        self.offset = offset = 0

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_sub_source_0 = zeromq.sub_source(gr.sizeof_gr_complex, 1, zmq_address_iq_in, -1, False, -1)
        self.lora_message_socket_sink_0 = lora.message_socket_sink('127.0.0.1', decoded_out_port, 0)
        self.lora_lora_receiver_0_0 = lora.lora_receiver(samp_rate, capture_freq, (capture_freq, ), bandwidth, spreading_factor, False, 4, True, False, False, 1, False, False)

        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.lora_lora_receiver_0_0, 'frames'), (self.lora_message_socket_sink_0, 'in'))
        self.connect((self.zeromq_sub_source_0, 0), (self.lora_lora_receiver_0_0, 0))

    def get_bandwidth(self):
        return self.bandwidth

    def set_bandwidth(self, bandwidth):
        self.bandwidth = bandwidth

    def get_capture_freq(self):
        return self.capture_freq

    def set_capture_freq(self, capture_freq):
        self.capture_freq = capture_freq

    def get_decoded_out_port(self):
        return self.decoded_out_port

    def set_decoded_out_port(self, decoded_out_port):
        self.decoded_out_port = decoded_out_port

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_spreading_factor(self):
        return self.spreading_factor

    def set_spreading_factor(self, spreading_factor):
        self.spreading_factor = spreading_factor
        self.lora_lora_receiver_0_0.set_sf(self.spreading_factor)

    def get_zmq_address_iq_in(self):
        return self.zmq_address_iq_in

    def set_zmq_address_iq_in(self, zmq_address_iq_in):
        self.zmq_address_iq_in = zmq_address_iq_in

    def get_offset(self):
        return self.offset

    def set_offset(self, offset):
        self.offset = offset


def argument_parser():
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option)
    parser.add_option(
        "", "--bandwidth", dest="bandwidth", type="intx", default=125000,
        help="Set bandwidth [default=%default]")
    parser.add_option(
        "", "--capture-freq", dest="capture_freq", type="eng_float", default=eng_notation.num_to_str(868500000),
        help="Set capture_freq [default=%default]")
    parser.add_option(
        "", "--decoded-out-port", dest="decoded_out_port", type="intx", default=40868,
        help="Set decoded_out_port [default=%default]")
    parser.add_option(
        "", "--samp-rate", dest="samp_rate", type="eng_float", default=eng_notation.num_to_str(1000000),
        help="Set samp_rate [default=%default]")
    parser.add_option(
        "", "--spreading-factor", dest="spreading_factor", type="intx", default=12,
        help="Set spreading_factor [default=%default]")
    parser.add_option(
        "", "--zmq-address-iq-in", dest="zmq_address_iq_in", type="string", default='tcp://127.0.0.1:5051',
        help="Set zmq_address_iq_in [default=%default]")
    return parser


def main(top_block_cls=zero_mq_split_b, options=None):
    if options is None:
        options, _ = argument_parser().parse_args()

    tb = top_block_cls(bandwidth=options.bandwidth, capture_freq=options.capture_freq, decoded_out_port=options.decoded_out_port, samp_rate=options.samp_rate, spreading_factor=options.spreading_factor, zmq_address_iq_in=options.zmq_address_iq_in)
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
