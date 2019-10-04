#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Zero Mq Split A
# Generated: Thu Sep 12 19:54:52 2019
##################################################


from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser


class zero_mq_split_a(gr.top_block):

    def __init__(self, zmq_address='tcp://127.0.0.1:5051'):
        gr.top_block.__init__(self, "Zero Mq Split A")

        ##################################################
        # Parameters
        ##################################################
        self.zmq_address = zmq_address

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 1e6
        self.offset = offset = 0
        self.capture_freq = capture_freq = 868.5e6

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, zmq_address, -1, False, -1)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, '/home/sili/Documents/signal_recordings/trimmed', True)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0, 0), (self.zeromq_pub_sink_0, 0))

    def get_zmq_address(self):
        return self.zmq_address

    def set_zmq_address(self, zmq_address):
        self.zmq_address = zmq_address

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
        "", "--zmq-address", dest="zmq_address", type="string", default='tcp://127.0.0.1:5051',
        help="Set zmq_address [default=%default]")
    return parser


def main(top_block_cls=zero_mq_split_a, options=None):
    if options is None:
        options, _ = argument_parser().parse_args()

    tb = top_block_cls(zmq_address=options.zmq_address)
    tb.start()
    try:
        raw_input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()
