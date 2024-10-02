FROM mambaorg/micromamba:git-a628e26-buster

COPY --chown=$MAMBA_USER:$MAMBA_USER ./ /tmp/bjorn
RUN micromamba install -y -n base -c conda-forge git
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN git clone --depth=1 --branch v4.3.1 https://github.com/cov-lineages/pangolin.git
RUN grep -ve '- defaults$' /tmp/pangolin/environment.yml > /tmp/pangolin.yml
RUN micromamba install -y -n base -f /tmp/pangolin.yml
RUN pip install ./pangolin
ENV GIT_CLONE_PROTECTION_ACTIVE=false
RUN pangolin --add-assignment-cache
RUN pangolin --update-data
RUN sed -iE 's/name,lineage_histogram =/name,_,lineage_histogram =/g' /opt/conda/lib/python3.8/site-packages/pangolin/utils/report_collation.py
RUN micromamba install --no-allow-uninstall -y -n base -f /tmp/bjorn/bjorn.yml
RUN micromamba clean --all -y
WORKDIR /tmp/bjorn
RUN echo "will cite" | parallel --citation || true
ENV PAGER=w3m
