FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y  python python-pip \
&& pip install pyzmq==18.1.0



WORKDIR /

COPY ./* /BBU_SEND_RESPONSE/

ENTRYPOINT ["BBU_SEND_RESPONSE/forever.py", "BBU_SEND_RESPONSE/lora_socket_server.py"]
