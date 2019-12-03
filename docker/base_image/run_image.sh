image_name=$1
if [ -z $image_name ]
then 
	echo "please pass the name of the image you want to run"
       	exit
fi
xhost +
DOCKER_XAUTH=/tmp/.docker.xauth
touch /tmp/.docker.xauth
xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $DOCKER_XAUTH nmerge -

docker run --net=host -i -t --entrypoint="bin/bash" --rm --privileged -e DISPLAY=$DISPLAY -e XAUTHORITY=$DOCKER_XAUTH -v /tmp/.docker.xauth:/tmp/.docker.xauth:ro -v /tmp/.X11-unix:/tmp/.X11-unix:ro -v /dev/bus/usb:/dev/bus/usb $image_name 
