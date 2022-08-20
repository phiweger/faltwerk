FROM continuumio/miniconda3:4.10.3
ARG DEBIAN_FRONTEND=noninteractive

# https://3dmol.csb.pitt.edu/doc/tutorial-learning_environment.html

RUN \
    apt update --allow-releaseinfo-change && \
    apt install -y procps git npm python3-pip wget gzip pigz bc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install flask flask_socketio eventlet

RUN \
    git clone https://github.com/3dmol/3Dmol.js.git && \
    cd 3Dmol.js && \
    npm install

RUN \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

# Default: Launch a server on port 5000, e.g. http://HOSTNAME:5000/
# Command in exec mode: Ctrl+C works, use "" not '' bc/ parsed into json
# https://docs.docker.com/engine/reference/builder/#cmd
CMD ["python3", "/3Dmol.js/packages/learning-environment/server.py", "-p", "5000"]
