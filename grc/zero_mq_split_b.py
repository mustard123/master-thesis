#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Zero Mq Split B
# Generated: Tue Oct 15 15:40:11 2019
##################################################


from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import lora


class zero_mq_split_b(gr.top_block):

    def __init__(self, decoded_out_port=40868, spreading_factor=12, zmq_address_iq_in='tcp://127.0.0.1:5051'):
        gr.top_block.__init__(self, "Zero Mq Split B")

        ##################################################
        # Parameters
        ##################################################
        self.decoded_out_port = decoded_out_port
        self.spreading_factor = spreading_factor
        self.zmq_address_iq_in = zmq_address_iq_in

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 1e6
        self.offset = offset = 0
        self.capture_freq = capture_freq = 868.5e6

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_sub_source_0 = zeromq.sub_source(gr.sizeof_gr_complex, 1, zmq_address_iq_in, -1, False, -1)
        self.lora_message_socket_sink_0 = lora.message_socket_sink('127.0.0.1', decoded_out_port, 0)
        self.lora_lora_receiver_0_0 = lora.lora_receiver(1e6, capture_freq, (capture_freq, ), 125000, spreading_factor, False, 4, True, False, False, 1, False, False)

        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.lora_lora_receiver_0_0, 'frames'), (self.lora_message_socket_sink_0, 'in'))
        self.connect((self.zeromq_sub_source_0, 0), (self.lora_lora_receiver_0_0, 0))

    def get_decoded_out_port(self):
        return self.decoded_out_port

    def set_decoded_out_port(self, decoded_out_port):
        self.decoded_out_port = decoded_out_port

    def get_spreading_factor(self):
        return self.spreading_factor

    def set_spreading_factor(self, spreading_factor):
        self.spreading_factor = spreading_factor
        self.lora_lora_receiver_0_0.set_sf(self.spreading_factor)

    def get_zmq_address_iq_in(self):
        return self.zmq_address_iq_in

    def set_zmq_address_iq_in(self, zmq_address_iq_in):
        self.zmq_address_iq_in = zmq_address_iq_in

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_offset(self):
        return self.offset

    def set_offset(self, offset):
        self.offset = offset

    def get_capture_freq(self):
        return self.capture_freq

    def set_capture_freq(self, capture_freq):
        self.capture_freq = capture_freq


def argument_parser():
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option)
    parser.add_option(
        "", "--decoded-out-port", dest="decoded_out_port", type="intx", default=40868,
        help="Set decoded_out_port [default=%default]")
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

    tb = top_block_cls(decoded_out_port=options.decoded_out_port, spreading_factor=options.spreading_factor, zmq_address_iq_in=options.zmq_address_iq_in)
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
