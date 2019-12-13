#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Channelizer
# Generated: Thu Dec 12 17:02:04 2019
##################################################


from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser


class channelizer(gr.top_block):

    def __init__(self, bandwidth=125e3, input_file='', output_file='channelized.raw', samp_rate=1e6):
        gr.top_block.__init__(self, "Channelizer")

        ##################################################
        # Parameters
        ##################################################
        self.bandwidth = bandwidth
        self.input_file = input_file
        self.output_file = output_file
        self.samp_rate = samp_rate

        ##################################################
        # Variables
        ##################################################

        self.variable_low_pass_filter_taps_0 = variable_low_pass_filter_taps_0 = firdes.low_pass(1.0, samp_rate, (bandwidth/2)+15000, 10000, firdes.WIN_HAMMING, 6.76)


        ##################################################
        # Blocks
        ##################################################
        self.freq_xlating_fir_filter_xxx_0 = filter.freq_xlating_fir_filter_ccc(1, (variable_low_pass_filter_taps_0), 0, samp_rate)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, input_file, False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, output_file, False)
        self.blocks_file_sink_0.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0, 0), (self.freq_xlating_fir_filter_xxx_0, 0))
        self.connect((self.freq_xlating_fir_filter_xxx_0, 0), (self.blocks_file_sink_0, 0))

    def get_bandwidth(self):
        return self.bandwidth

    def set_bandwidth(self, bandwidth):
        self.bandwidth = bandwidth

    def get_input_file(self):
        return self.input_file

    def set_input_file(self, input_file):
        self.input_file = input_file
        self.blocks_file_source_0.open(self.input_file, False)

    def get_output_file(self):
        return self.output_file

    def set_output_file(self, output_file):
        self.output_file = output_file
        self.blocks_file_sink_0.open(self.output_file)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_variable_low_pass_filter_taps_0(self):
        return self.variable_low_pass_filter_taps_0

    def set_variable_low_pass_filter_taps_0(self, variable_low_pass_filter_taps_0):
        self.variable_low_pass_filter_taps_0 = variable_low_pass_filter_taps_0
        self.freq_xlating_fir_filter_xxx_0.set_taps((self.variable_low_pass_filter_taps_0))


def argument_parser():
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option)
    parser.add_option(
        "", "--bandwidth", dest="bandwidth", type="eng_float", default=eng_notation.num_to_str(125e3),
        help="Set bandwidth [default=%default]")
    parser.add_option(
        "", "--input-file", dest="input_file", type="string", default='',
        help="Set input_file [default=%default]")
    parser.add_option(
        "", "--output-file", dest="output_file", type="string", default='channelized.raw',
        help="Set output_file [default=%default]")
    parser.add_option(
        "", "--samp-rate", dest="samp_rate", type="eng_float", default=eng_notation.num_to_str(1e6),
        help="Set samp_rate [default=%default]")
    return parser


def main(top_block_cls=channelizer, options=None):
    if options is None:
        options, _ = argument_parser().parse_args()

    tb = top_block_cls(bandwidth=options.bandwidth, input_file=options.input_file, output_file=options.output_file, samp_rate=options.samp_rate)
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
