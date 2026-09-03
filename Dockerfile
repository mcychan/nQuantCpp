FROM ubuntu:24.04

ENV WINEARCH=win64
ENV WINEPREFIX=/root/.wine64
ENV DISPLAY=""
ENV WINEDEBUG=-all

WORKDIR /tmp

RUN dpkg --add-architecture i386 && apt update -y && apt install -y \
    build-essential \
    cmake \
    mingw-w64 \
    wget \
    fontconfig \
    locales \
    wine \
    wine64 \
    wine32:i386 \
    libgdiplus \
    libgif7 \
    libjpeg-turbo8 \
    libpng16-16 \
    && rm -rf /var/lib/apt/lists/*

RUN wineboot --init

RUN echo "set(CMAKE_SYSTEM_NAME Windows)" > /tmp/mingw-toolchain.cmake && \
    echo "set(CMAKE_C_COMPILER /usr/bin/x86_64-w64-mingw32-gcc)" >> /tmp/mingw-toolchain.cmake && \
    echo "set(CMAKE_CXX_COMPILER /usr/bin/x86_64-w64-mingw32-g++)" >> /tmp/mingw-toolchain.cmake && \
    echo "set(CMAKE_RC_COMPILER /usr/bin/x86_64-w64-mingw32-windres)" >> /tmp/mingw-toolchain.cmake

ADD . /tmp/nQuantCpp
WORKDIR /tmp/nQuantCpp

RUN cmake -S . -B ../build \
    -DCMAKE_TOOLCHAIN_FILE=/tmp/mingw-toolchain.cmake \
    -DCMAKE_CXX_FLAGS="-fopenmp" \
    -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++ -Wl,-Bstatic -lstdc++ -lgomp -lwinpthread -lole32 -Wl,-Bdynamic" \
    -DCMAKE_CXX_STANDARD_LIBRARIES="-Wl,-Bstatic -lstdc++ -lgomp -lwinpthread -lole32 -Wl,-Bdynamic"
RUN cmake --build ../build
RUN cp *.gif /tmp/build/nQuantCpp/
WORKDIR /tmp/build/nQuantCpp

# docker system prune -a
# docker build -t nquantcpp .
# docker run -it nquantcpp bash
# wine ./nQuantCpp.exe ImgV64.gif /a pnnlab
# docker cp <containerId>:/file/path/within/container /host/path/target
# docker cp foo.txt <containerId>:/foo.txt
