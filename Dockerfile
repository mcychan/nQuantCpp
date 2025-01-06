FROM ubuntu:22.04
WORKDIR /tmp
RUN apt update -y
RUN apt install -y build-essential cmake gcc-12 g++-12 libomp-dev
RUN DEBIAN_FRONTEND="noninteractive" apt install -y libgdiplus mingw-w64 locales fontconfig
ADD . /tmp/nQuantCpp
WORKDIR /tmp/nQuantCpp
RUN cmake -D CMAKE_CXX_COMPILER=g++-12 /usr/bin/g++-12 -S . -B ../build
# RUN cmake --build ../build
# RUN cp *.gif /tmp/build/nQuantCpp/
# WORKDIR /tmp/build/nQuantCpp
# docker system prune -a
# docker build -t nquantcpp .
# docker run -it nquantcpp bash
# docker cp <containerId>:/file/path/within/container /host/path/target
# docker cp foo.txt <containerId>:/foo.txt
