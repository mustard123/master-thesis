#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Zero Mq Split A
# Generated: Thu Jan 23 22:06:59 2020
##################################################


from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import limesdr


class zero_mq_split_a(gr.top_block):

    def __init__(self, zmq_address_iq_in='tcp://127.0.0.1:5052', zmq_address_iq_out='tcp://*:5051'):
        gr.top_block.__init__(self, "Zero Mq Split A")

        ##################################################
        # Parameters
        ##################################################
        self.zmq_address_iq_in = zmq_address_iq_in
        self.zmq_address_iq_out = zmq_address_iq_out

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
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, zmq_address_iq_out, -1, False, -1)
        self.limesdr_source_0 = limesdr.source('0009072C0287211A', 0, '')
        self.limesdr_source_0.set_sample_rate(samp_rate)
        self.limesdr_source_0.set_center_freq(capture_freq, 0)
        self.limesdr_source_0.set_bandwidth(5e6,0)
        self.limesdr_source_0.set_gain(30,0)
        self.limesdr_source_0.set_antenna(2,0)
        self.limesdr_source_0.calibrate(5e6, 0)

        self.limesdr_sink_0 = limesdr.sink('', 0, '', '')
        self.limesdr_sink_0.set_sample_rate(samp_rate)
        self.limesdr_sink_0.set_center_freq(capture_freq, 0)
        self.limesdr_sink_0.set_bandwidth(5e6,0)
        self.limesdr_sink_0.set_gain(50,0)
        self.limesdr_sink_0.set_antenna(255,0)
        self.limesdr_sink_0.calibrate(5e6, 0)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.limesdr_source_0, 0), (self.zeromq_pub_sink_0, 0))
        self.connect((self.zeromq_sub_source_0, 0), (self.limesdr_sink_0, 0))

    def get_zmq_address_iq_in(self):
        return self.zmq_address_iq_in

    def set_zmq_address_iq_in(self, zmq_address_iq_in):
        self.zmq_address_iq_in = zmq_address_iq_in

    def get_zmq_address_iq_out(self):
        return self.zmq_address_iq_out

    def set_zmq_address_iq_out(self, zmq_address_iq_out):
        self.zmq_address_iq_out = zmq_address_iq_out

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
        self.limesdr_source_0.set_center_freq(self.capture_freq, 0)
        self.limesdr_sink_0.set_center_freq(self.capture_freq, 0)


def argument_parser():
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option)
    parser.add_option(
        "", "--zmq-address-iq-in", dest="zmq_address_iq_in", type="string", default='tcp://127.0.0.1:5052',
        help="Set zmq_address_iq_in [default=%default]")
    parser.add_option(
        "", "--zmq-address-iq-out", dest="zmq_address_iq_out", type="string", default='tcp://*:5051',
        help="Set zmq_address_iq_out [default=%default]")
    return parser


def main(top_block_cls=zero_mq_split_a, options=None):
    if options is None:
        options, _ = argument_parser().parse_args()

    tb = top_block_cls(zmq_address_iq_in=options.zmq_address_iq_in, zmq_address_iq_out=options.zmq_address_iq_out)
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
