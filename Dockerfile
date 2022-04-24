FROM gcc:latest

RUN apt update && apt install -y libgsl-dev

COPY . .

RUN make

ENTRYPOINT [ "./main" ]
