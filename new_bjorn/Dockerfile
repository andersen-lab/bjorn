FROM mambaorg/micromamba:git-a628e26-buster

COPY --chown=$MAMBA_USER:$MAMBA_USER ./ /tmp/bjorn
RUN micromamba install -y -n base -c conda-forge git
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN git clone --depth=1 --branch v4.1.2 https://github.com/cov-lineages/pangolin.git
RUN grep -ve '- defaults$' /tmp/pangolin/environment.yml > /tmp/pangolin.yml
RUN micromamba install -y -n base -f /tmp/pangolin.yml
RUN pip install ./pangolin
RUN pangolin --add-assignment-cache
RUN pangolin --update-data
RUN micromamba install --no-allow-uninstall -y -n base -f /tmp/bjorn/newbjorn.yml
RUN micromamba clean --all -y
WORKDIR /tmp/bjorn
RUN echo "will cite" | parallel --citation || true
