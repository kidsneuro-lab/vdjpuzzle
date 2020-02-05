FROM syin2018/ccia:2018workshop2
ENV LANG=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
RUN yum install -y wget bzip2 git ca-certificates curl which
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda && \
rm ~/miniconda.sh && \
mkdir -p /opt/vdjpuzzle && \
mkdir -p /data && \
mkdir -p /usr/share/Modules/modulefiles/vdjpuzzle && \
/opt/conda/bin/conda clean -tipsy && \
ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
echo "conda activate base" >> ~/.bashrc && \
echo 'export PATH=/opt/vdjpuzzle/bin:$PATH' >> ~/.bashrc && \
echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/pbsuser/.bashrc && \
echo "conda activate vdjpuzzle" >> /home/pbsuser/.bashrc && \
echo 'export PATH=/opt/vdjpuzzle/bin:$PATH' >> /home/pbsuser/.bashrc
ADD ./environment.yml /opt/vdjpuzzle/environment.yml
ADD ./bin /opt/vdjpuzzle/bin
ADD ./azcopy/azcopy /opt/vdjpuzzle/bin/azcopy
ADD ./azcopy/aget /opt/vdjpuzzle/bin/aget
ADD ./azcopy/als /opt/vdjpuzzle/bin/als
ADD ./azcopy/aput /opt/vdjpuzzle/bin/aput
ADD ./azcopy/arm /opt/vdjpuzzle/bin/arm
ADD ./azcopy/armdir /opt/vdjpuzzle/bin/armdir
ADD ./azcopy/asas /opt/vdjpuzzle/bin/asas
ADD ./azcopy/asy /opt/vdjpuzzle/bin/asy
ADD ./db /opt/vdjpuzzle/db
ADD ./README.md /opt/vdjpuzzle/README.md
ADD ./share /opt/vdjpuzzle/share
ADD ./Example /opt/vdjpuzzle/Example
ADD ./LICENSE.txt /opt/vdjpuzzle/LICENSE.txt
ADD ./scripts /opt/vdjpuzzle/scripts
RUN sed -i 's/#PBS -l wd/cd "${PBS_O_WORKDIR}"/g' /opt/vdjpuzzle/scripts/run_part_1_and_2.sh && \
sed -i 's/#PBS -q normal/#PBS -q workq/g' /opt/vdjpuzzle/scripts/run_part_1_and_2.sh && \
sed -i 's/#PBS -l wd/cd "${PBS_O_WORKDIR}"/g' /opt/vdjpuzzle/scripts/run_part3.sh && \
sed -i 's/#PBS -q normal/#PBS -q workq/g' /opt/vdjpuzzle/scripts/run_part3.sh
ADD ./modulefiles/vdjpuzzle /usr/share/Modules/modulefiles/vdjpuzzle/2019.1
RUN /opt/conda/bin/conda config --add channels defaults && \
/opt/conda/bin/conda config --add channels bioconda && \
/opt/conda/bin/conda config --add channels conda-forge && \
/opt/conda/bin/conda env create -f /opt/vdjpuzzle/environment.yml &&\
chown -R pbsuser:pbsuser /opt/vdjpuzzle &&\
chown -R pbsuser:pbsuser /opt/conda
CMD ["/bin/bash"]
