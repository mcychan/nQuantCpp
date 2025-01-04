FROM ubuntu:22.04
WORKDIR /tmp
RUN apt update -y
RUN apt install -y build-essential cmake gcc-12 g++-12 libomp-dev wget
RUN apt install -y libgdiplus locales fontconfig
RUN wget -O - https://dl.winehq.org/wine-builds/winehq.key | apt-key add -
RUN add-apt-repository 'deb https://dl.winehq.org/wine-builds/ubuntu/ jammy main'
RUN apt install -y --install-recommends winehq-stable-dev
ADD . /tmp/nQuantCpp
WORKDIR /tmp/nQuantCpp
RUN cmake -S . -B ../build
# RUN cmake --build ../build
# RUN cp *.gif /tmp/build/nQuantCpp/
# WORKDIR /tmp/build/nQuantCpp
# docker system prune -a
# docker build -t nquantcpp .
# docker run -it nquantcpp bash
# docker cp <containerId>:/file/path/within/container /host/path/target
# docker cp foo.txt <containerId>:/foo.txt
