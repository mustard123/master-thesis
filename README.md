# master-thesis

## Docker installation
The docker image contains gnuaradio with the compiled LimeSuite SDR modules and the lora modules from https://github.com/rpp0/gr-lora. This allows you to get quickly setup without tedious manual installation

* clone the repo
* cd into docker directory
* run  ``` docker build -t <myimagename> . ```
* then run ```./run_image.sh <myimagename> ```
* inside the container run ``` ./gnuradio.sh ```
* In the gnuradio GUI click on open file an select a .grc file in the /home directory