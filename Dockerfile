FROM ppc64le/ubuntu:18.04
RUN apt-get update && \
    apt-get install -y python && \
    apt-get install -y wget \
    curl \
    bc \
    vim \
    time \
    unzip \
    less \
    bedtools \
    samtools \
    openjdk-8-jdk \
    tabix \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove
ENV HOME /root
ENV export JAVA_HOME /usr/lib/jvm/java-8-openjdk-ppc64el
ENV JAVA_LIBRARY_PATH /usr/lib/powerpc64le-linux-gnu/jni
CMD ["bash"]
RUN mkdir -p /apps/bin
RUN mkdir -p /work
RUN mkdir -p /data
WORKDIR  /work
ADD gatk-4.1.4.1 /apps/gatk-4.1.4.1
COPY bin /apps/bin
ADD benchmarks /apps/benchmarks
ENV GATK_HOME /apps
ENV PATH $GATK_HOME/bin:$PATH
ENV GATK_LOCAL_JAR $GATK_HOME/gatk-4.1.4.1/libs/gatk.jar
ENV GATK_SPARK_JAR $GATK_HOME/gatk-4.1.4.1/libs/gatk-spark.jar
ENV LD_LIBRARY_PATH $GATK_HOME/gatk-4.1.4.1/libs:$LD_LIBRARY_PATH
