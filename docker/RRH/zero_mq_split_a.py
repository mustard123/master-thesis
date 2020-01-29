#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Zero Mq Split A
# Generated: Wed Jan 29 18:30:19 2020
##################################################


from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import limesdr


class zero_mq_split_a(gr.top_block):

    def __init__(self, RX_device_serial='', TX_device_serial='', capture_freq=868500000, samp_rate=1000000, zmq_address_iq_in='tcp://127.0.0.1:5052', zmq_address_iq_out='tcp://*:5051'):
        gr.top_block.__init__(self, "Zero Mq Split A")

        ##################################################
        # Parameters
        ##################################################
        self.RX_device_serial = RX_device_serial
        self.TX_device_serial = TX_device_serial
        self.capture_freq = capture_freq
        self.samp_rate = samp_rate
        self.zmq_address_iq_in = zmq_address_iq_in
        self.zmq_address_iq_out = zmq_address_iq_out

        ##################################################
        # Variables
        ##################################################
        self.offset = offset = 0

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_sub_source_0 = zeromq.sub_source(gr.sizeof_gr_complex, 1, zmq_address_iq_in, -1, False, -1)
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, zmq_address_iq_out, -1, False, -1)
        self.limesdr_source_0 = limesdr.source(RX_device_serial, 0, '')
        self.limesdr_source_0.set_sample_rate(samp_rate)
        self.limesdr_source_0.set_center_freq(capture_freq, 0)
        self.limesdr_source_0.set_bandwidth(5e6,0)
        self.limesdr_source_0.set_gain(30,0)
        self.limesdr_source_0.set_antenna(2,0)
        self.limesdr_source_0.calibrate(5e6, 0)

        self.limesdr_sink_0 = limesdr.sink(TX_device_serial, 0, '', '')
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

    def get_RX_device_serial(self):
        return self.RX_device_serial

    def set_RX_device_serial(self, RX_device_serial):
        self.RX_device_serial = RX_device_serial

    def get_TX_device_serial(self):
        return self.TX_device_serial

    def set_TX_device_serial(self, TX_device_serial):
        self.TX_device_serial = TX_device_serial

    def get_capture_freq(self):
        return self.capture_freq

    def set_capture_freq(self, capture_freq):
        self.capture_freq = capture_freq
        self.limesdr_source_0.set_center_freq(self.capture_freq, 0)
        self.limesdr_sink_0.set_center_freq(self.capture_freq, 0)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_zmq_address_iq_in(self):
        return self.zmq_address_iq_in

    def set_zmq_address_iq_in(self, zmq_address_iq_in):
        self.zmq_address_iq_in = zmq_address_iq_in

    def get_zmq_address_iq_out(self):
        return self.zmq_address_iq_out

    def set_zmq_address_iq_out(self, zmq_address_iq_out):
        self.zmq_address_iq_out = zmq_address_iq_out

    def get_offset(self):
        return self.offset

    def set_offset(self, offset):
        self.offset = offset


def argument_parser():
    parser = OptionParser(usage="%prog: [options]", option_class=eng_option)
    parser.add_option(
        "", "--RX-device-serial", dest="RX_device_serial", type="string", default='',
        help="Set RX_device_serial [default=%default]")
    parser.add_option(
        "", "--TX-device-serial", dest="TX_device_serial", type="string", default='',
        help="Set TX_device_serial [default=%default]")
    parser.add_option(
        "", "--capture-freq", dest="capture_freq", type="eng_float", default=eng_notation.num_to_str(868500000),
        help="Set capture_freq [default=%default]")
    parser.add_option(
        "", "--samp-rate", dest="samp_rate", type="eng_float", default=eng_notation.num_to_str(1000000),
        help="Set samp_rate [default=%default]")
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

    tb = top_block_cls(RX_device_serial=options.RX_device_serial, TX_device_serial=options.TX_device_serial, capture_freq=options.capture_freq, samp_rate=options.samp_rate, zmq_address_iq_in=options.zmq_address_iq_in, zmq_address_iq_out=options.zmq_address_iq_out)
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
