FROM ubuntu:rolling

RUN apt-get update -yq \
    && apt-get install locales perl curl wget build-essential tar gzip pkg-config python3 python3-dev python3-pip git -yq \
    && localedef -i fr_FR -c -f UTF-8 -A /usr/share/locale/locale.alias fr_FR.UTF-8 \
    && pip3 install Biopython numpy pandas matplotlib joblib pyvcf pysam seaborn pytest xlsxwriter dna_features_viewer \
    && mkdir -p /root/dev && mkdir -p /root/bin

ENV LANG fr_FR.utf8

ADD . /root/dev/GenomicPackage

WORKDIR /root

CMD ["/bin/bash", "/root/dev/GenomicPackage/docker/execute.sh"]

