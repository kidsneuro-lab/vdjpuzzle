FROM continuumio/miniconda3:23.5.2-0

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Copy all the files
COPY . .

# Ensure conda env is activated and PATH variables are set
RUN echo "conda activate vdjpuzzle" >> ~/.bashrc && \
    echo 'export PATH=/app/scripts:$PATH' >> ~/.bashrc && \
    echo 'export PATH=/app/azcopy:$PATH' >> ~/.bashrc

RUN chmod u+x /app/scripts/*.sh && \
    chmod u+x /app/azcopy/*

# Set the default shell. The --login bit will ensure that .bashrc is run
SHELL ["/bin/bash", "--login", "-c"]

ENTRYPOINT ["/bin/bash", "--login", "-c"]
CMD [ "/app/bin/vdjpuzzle" ]