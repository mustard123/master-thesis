# C-RAN for LoRa

## Run with Docker

1. Clone the repo
2. Go to the docker directory

### ***Info***

* The container run in priviledged mode to easily access plugged in USB devices
* The container run in network mode host (No NAT or Bridge has to be considered). This means the containers have the ip address of the host machine. If RRH and BBU run on different machines, find out their respective IP with *ifconfig* and pass the address as arguments in the docker-compose.yml, see below.

---

## RRH

In the RRH directory run:

``` 
docker-compose up 
```

This starts the Remote Radio Head. The RRH looks for a LimeSDR, it prints errors if it cannot find one. You can plug one in after the container has started and it should get detectet. By default it uses the first LimeSDR it can find.

### Parameters

There are various parameters which you can specify in the *docker-compose.yml* file

run this to see what the possible params are:
``` 
./zero_mq_split_a.py -h
```

Output:
```
Usage: zero_mq_split_a.py: [options]

Options:
  -h, --help            show this help message and exit
  --RX-device-serial=RX_DEVICE_SERIAL
                        Set RX_device_serial [default=]
  --TX-device-serial=TX_DEVICE_SERIAL
                        Set TX_device_serial [default=]
  --capture-freq=CAPTURE_FREQ
                        Set capture_freq [default=868.5M]
  --samp-rate=SAMP_RATE
                        Set samp_rate [default=1.0M]
  --zmq-address-iq-in=ZMQ_ADDRESS_IQ_IN
                        Set zmq_address_iq_in [default=tcp://127.0.0.1:5052]
  --zmq-address-iq-out=ZMQ_ADDRESS_IQ_OUT
                        Set zmq_address_iq_out [default=tcp://*:5051]

```

| Param               | Explanation        |
| -------------       | ------------- |
| RX-device-serial  | By default, the program will use the first LimeSDR it can find for receiving and transmitting signal. If you have two devices you can specify which should receive by passing the device Serial (See section **Help** for more info)  |
| TX-device-serial        | By default, the program will use the first LimeSDR it can find for receiving and transmitting signal. If you have two devices you can specify which should transmit by passing the device Serial (See section **Help** for more info)  |
|capture-freq|The frequency in Hz at which the RRH listens for signals. Default value  is 86850000|
|samp-rate|How many samples per second. Default value is 1000000. Must be at least double the bandwidth of the expected signal see *Nyquist-Shannon principle*|
|zmq-address-iq-in| ZMQ address to which the RRH subscribes to receive an IQ samples stream (from the BBU) to then send out (TX). Default value is tcp://127.0.0.1:5052 meaning the IQ samples are expected to come from localhost on port 5052. Normally RRH and BBU are on different devices but on  the same network |
|--zmq-address-iq-out|ZMQ address on which the RRH streams out the IQ samples (to the BBU) it receives (RX). Default is tcp://*:5051 meaning it publishes the stream on all interface on port 5051|

---
To pass the parameters you have to specify them in the docker-compose.yml

Example: 

To pass a capture frequencey of 915M and a sample rate of 250k enter the params in the following way in the command field:

*docker-compose.yml*
``` 
version: '3'
services:
    rrh:
        build: .
        privileged: true
        network_mode: host
        volumes:
                - /dev/bus/usb:/dev/bus/usb
        command: ["--capture-freq", "915000000", "--samp_rate", "250000"]
```

---
## BBU

The BBU has two components: 
* LoRa_Decoder: receives a stream of IQ samples from the RRH, decodes the LoRa signal and sends the decoded message out on a UDP socket
* LoRa_Responder: receives the messages from that UDP socket and, depending on message content, streams response IQ samples to the RRH or does not give a response

In the BBU directory run:

``` 
docker-compose up 
```

This starts both components of the BBU

### Params

The LoRa_Decoder has the following params:

```
Usage: zero_mq_split_b.py: [options]

Options:
  -h, --help            show this help message and exit
  --bandwidth=BANDWIDTH
                        Set bandwidth [default=125000]
  --capture-freq=CAPTURE_FREQ
                        Set capture_freq [default=868.5M]
  --decoded-out-port=DECODED_OUT_PORT
                        Set decoded_out_port [default=40868]
  --samp-rate=SAMP_RATE
                        Set samp_rate [default=1.0M]
  --spreading-factor=SPREADING_FACTOR
                        Set spreading_factor [default=12]
  --zmq-address-iq-in=ZMQ_ADDRESS_IQ_IN
                        Set zmq_address_iq_in [default=tcp://127.0.0.1:5051]

```
| Param               | Explanation        |
| -------------       | ------------- |
| bandwith  | The bandwidth in Hz of the LoRa signal. Default is 125000.   |
| capture-freq        | The frequency in Hz of the LoRa signal. The RRH of course must also listen on this frequeny. Default is 868500000. |
|decoded-out-port|On which port the decoded messages will be sent out. Localhost only. The LoRa_Responder needs to be configured to listen on this port. Default is 40868.|
|samp-rate| How many samples per second to expect from the RRH. Default is 1000000|
|spreading-factor| The spreading factor of the incoming LoRa signal. From [7-12] inclusive. Default is 12 |
|--zmq-address-iq-in|ZMQ address to which the BBU subscribes to receive an IQ samples stream (from the RRH) to decode. Default value is tcp://127.0.0.1:5051 meaning the IQ samples are expected to come from localhost on port 5051. Normally RRH and BBU are on different devices but on the same network|

---

The LoRa_Responder has the following params:

```
usage: lora_socket_server.py [-h] [-o OUT_PORT] [-i INPUT_PORT]

Connect to udp port for receiving decoded LoRa signals, if an ACK is required
publish ACK iq samples via zmq socket for Remote Radio Head to receive and
send out (TX).

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PORT, --out-port OUT_PORT
                        zmq port to publish downstream (i.e ACK) iq samples
                        (default: 5052)
  -i INPUT_PORT, --input-port INPUT_PORT
                        UDP port to connect for receiving decoded lora
                        messages (default: 40868)


```

| Param               | Explanation        |
| -------------       | ------------- |
| out-port  | Publish the response IQ samples on all interface on this port. Default is  5052. (The response is 23 bytes long and SF 12. This is hardcoded for now)|
|input-port| UDP port to receive the decoded messages sent by the LoRa_Decoder. Default is 40868|
---

To pass the parameters you have to specify them in the docker-compose.yml

Example: 

To have the LoRa_Decoder send the decoded messages out on port 30300 and the Lora_Responder to listen on port 30300 accordingly pass the arguments like below to the respective command field: 

*docker-compose.yml*
``` 
version: '3'
services:
        lora_decoder:
                build: ./LoRa_Decoder
                network_mode: host
                tty: true
                command: ["--decoded-out-port", "30300"]
        lora_responder:
                build: ./LoRa_Responder
                network_mode: host
                tty: true
                command: ["--input-port", "30300"] 
```

## LimeSDR

* Plug in the antennas on the LimeSDR board on *RX1_L* and *TX1_1*

## Help

* LimeSDR calibration error: 
  * [Download LimeSuite Toolkit](https://wiki.myriadrf.org/Lime_Suite) to calibrate the LimeSDR

* LimeSDR find device serial:
  * With LimeSuite installed run *LimeUtil --find*
  * Or run *lsusb -v* and look for the LimeSDR device

