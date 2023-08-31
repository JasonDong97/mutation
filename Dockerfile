FROM ubuntu:22.04
LABEL maintainer="lingyu zeng <pylyzeng@gmail.com>"

WORKDIR /root
# COPY ../../share_data/opus_mut.zip ./
COPY 4i24.pdb test.list _foldxLinux64.tar_.gz EvoEF2-master.zip install_scwrl4.0.2_64bit_2020_linux mutation.py pyrosetta-2023.31+release.1799523-py311_0.tar.bz2 noarch/repodata.json ./
ENV PATH="/root/bin:/root/micromamba/bin:${PATH}"
ENV CONDA_PREFIX="/root/micromamba/envs/pyrosetta"
ENV TZ="Asia/Shanghai"
ENV DEBIAN_FRONTEND="noninteractive"

# Update and install necessary packages
RUN apt-get update && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    apt-get install tzdata -y && \
    apt-get install git zip wget bzip2 libgl1-mesa-glx g++ -y && \
    unzip EvoEF2-master.zip && \
    chmod +x EvoEF2-master/build.sh && \
    cd EvoEF2-master && \
    ./build.sh && \
    cd .. && \
    tar zxvf _foldxLinux64.tar_.gz && \
    chmod +x install_scwrl4.0.2_64bit_2020_linux && \
    echo -e "Y\nLicense Holder Name" | ./install_scwrl4.0.2_64bit_2020_linux ./ && \
    # Install Micromamba
    wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    ./bin/micromamba shell init --shell bash --root-prefix=~/micromamba && \
    echo -e "show_channel_urls: true\n\
channels:\n\
    - conda-forge\n\
    - r\n\
    - defaults\n\
    - bioconda\n\
    - https://levinthal:paradox@conda.graylab.jhu.edu\n\
custom_channels:\n\
    conda-forge: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
    msys2: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
    bioconda: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
    menpo: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
    pytorch: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
    simpleitk: https://mirrors.bfsu.edu.cn/anaconda/cloud\n\
default_channels:\n\
    - https://mirrors.bfsu.edu.cn/anaconda/pkgs/main\n\
    - https://mirrors.bfsu.edu.cn/anaconda/pkgs/r\n\
    - https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2" >> ~/.condarc && \
    mkdir -p /root/noarch && \
	mv repodata.json /root/noarch && \
    mv pyrosetta-2023.31+release.1799523-py311_0.tar.bz2 /root/noarch && \
    ./bin/micromamba create -n pyrosetta -c conda-forge -c bioconda -c defaults python=3.11 click loguru biopython pymol-open-source pyrosetta-2023.31+release.1799523-py311_0.tar.bz2 -y && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* EvoEF2-master.zip install_scwrl4.0.2_64bit_2020_linux _foldxLinux64.tar_.gz /root/noarch

WORKDIR /work
VOLUME ["/work"]
ENTRYPOINT ["/root/micromamba/envs/pyrosetta/bin/python", "/root/mutation.py"]
CMD ["--help"]
