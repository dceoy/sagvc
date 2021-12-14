FROM ubuntu:20.04 AS builder

ENV DEBIAN_FRONTEND noninteractive

COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bcftools:latest /usr/local/src/bcftools /usr/local/src/bcftools
COPY --from=dceoy/bedtools:latest /usr/local/src/bedtools2 /usr/local/src/bedtools2
COPY --from=dceoy/gatk:latest /opt/conda /opt/conda
COPY --from=dceoy/gatk:latest /opt/gatk /opt/gatk
COPY --from=dceoy/delly:latest /usr/local/bin/delly /usr/local/bin/delly
COPY --from=dceoy/manta:latest /opt/manta /opt/manta
COPY --from=dceoy/strelka:latest /opt/strelka /opt/strelka
COPY --from=dceoy/msisensor:latest /usr/local/bin/msisensor /usr/local/bin/msisensor
COPY --from=dceoy/msisensor-pro:latest /usr/local/bin/msisensor-pro /usr/local/bin/msisensor-pro
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD . /tmp/sagvc

RUN set -e \
      && ln -sf bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        g++ gcc libbz2-dev libcurl4-gnutls-dev libgsl-dev liblzma-dev \
        libncurses5-dev libperl-dev libssl-dev libz-dev make perl pkg-config \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/gatk/bin:/opt/conda/envs/gatk/bin:/opt/conda/bin:${PATH}

RUN set -e \
      && /opt/conda/bin/python3 /tmp/get-pip.py \
      && /opt/conda/bin/python3 -m pip install -U --no-cache-dir \
        cnvkit cutadapt https://github.com/dceoy/ftarc/archive/main.tar.gz \
        /tmp/sagvc \
      && /opt/conda/bin/conda update -n base -c defaults conda \
      && echo >> /opt/gatk/gatkcondaenv.yml \
      && echo -e '# CNVkit' >> /opt/gatk/gatkcondaenv.yml \
      && echo -e '- bioconductor-dnacopy' >> /opt/gatk/gatkcondaenv.yml \
      && /opt/conda/bin/conda env update -n gatk -f /opt/gatk/gatkcondaenv.yml \
      && /opt/conda/bin/conda clean -yaf \
      && find /opt/conda -follow -type f -name '*.a' -delete \
      && find /opt/conda -follow -type f -name '*.pyc' -delete \
      && rm -rf /root/.cache/pip /tmp/get-pip.py

RUN set -e \
      && cd /usr/local/src/samtools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/samtools \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/bcftools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && cd /usr/local/src/bcftools \
      && make clean \
      && ./configure --enable-libgsl --enable-perl-filters \
      && make \
      && make install \
      && cd /usr/local/src/bedtools2 \
      && make clean \
      && make \
      && make install


FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local /usr/local
COPY --from=builder /opt /opt

RUN set -e \
      && ln -sf bash /bin/sh \
      && echo '. /opt/conda/etc/profile.d/conda.sh' >> /etc/profile \
      && echo 'source activate gatk' >> /etc/profile \
      && echo 'source /opt/gatk/gatk-completion.sh' >> /etc/profile

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates curl gnupg gnuplot \
        libcurl3-gnutls libgsl23 libgkl-jni libncurses5 openjdk-8-jre pbzip2 \
        perl pigz python texlive-fonts-extra texlive-fonts-recommended \
        texlive-latex-base texlive-latex-extra wget

RUN set -eo pipefail \
      && echo 'deb http://packages.cloud.google.com/apt cloud-sdk-bionic main' \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
      && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | apt-key add - \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        google-cloud-sdk \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && unlink /usr/lib/ssl/openssl.cnf \
      && echo -e 'openssl_conf = default_conf' > /usr/lib/ssl/openssl.cnf \
      && echo >> /usr/lib/ssl/openssl.cnf \
      && cat /etc/ssl/openssl.cnf >> /usr/lib/ssl/openssl.cnf \
      && echo >> /usr/lib/ssl/openssl.cnf \
      && echo -e '[default_conf]' >> /usr/lib/ssl/openssl.cnf \
      && echo -e 'ssl_conf = ssl_sect' >> /usr/lib/ssl/openssl.cnf \
      && echo >> /usr/lib/ssl/openssl.cnf \
      && echo -e '[ssl_sect]' >> /usr/lib/ssl/openssl.cnf \
      && echo -e 'system_default = system_default_sect' >> /usr/lib/ssl/openssl.cnf \
      && echo >> /usr/lib/ssl/openssl.cnf \
      && echo -e '[system_default_sect]' >> /usr/lib/ssl/openssl.cnf \
      && echo -e 'MinProtocol = TLSv1.2' >> /usr/lib/ssl/openssl.cnf \
      && echo -e 'CipherString = DEFAULT:@SECLEVEL=1' >> /usr/lib/ssl/openssl.cnf

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV CLASSPATH /opt/gatk/gatk.jar:${CLASSPATH}
ENV BCFTOOLS_PLUGINS /usr/local/src/bcftools/plugins
ENV PYTHONPATH /opt/manta/lib/python:/opt/strelka/lib/python:${PYTHONPATH}
ENV PATH /opt/gatk/bin:/opt/conda/envs/gatk/bin:/opt/conda/bin:/opt/manta/bin:/opt/strelka/bin:${PATH}
ENV MPLCONFIGDIR /tmp/mpl

ENTRYPOINT ["/opt/conda/bin/sagvc"]
