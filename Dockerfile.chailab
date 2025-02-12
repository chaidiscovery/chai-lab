FROM ubuntu:22.04 AS chailab-baseimage

ENV \
  LANG=C.UTF-8 \
  LC_ALL=C.UTF-8 \
  # config for apt
  DEBIAN_FRONTEND=noninteractive \
  # default editor for git cli
  EDITOR=vim \
  # keep (large) mypy cache outside of working tree
  MYPY_CACHE_DIR='/tmp/.chai_lab_mypy_cache' \
  # always flush output from python
  PYTHONUNBUFFERED=TRUE \
  # enable fault handler (print tracebacks even after segfault or NCCL errors).
  PYTHONFAULTHANDLER=1 \
  # keep __pycache__ out of working tree
  PYTHONPYCACHEPREFIX='/tmp/.chai_lab_pycache'


RUN --mount=type=cache,target=/var/cache/apt \
  apt-get -qq update \
  && apt-get -qq install -y \
  # common things
  gnupg ca-certificates wget git curl aria2 lsb-release tzdata \
  rsync sudo tree htop tmux unzip \
  clang \
  # for direct ssh into container
  openssh-server socat \
  # provides `fuser` command
  psmisc \
  # RDMA/InfiniBand
  libibverbs1 librdmacm1 \
  # text editors, needed by git cli
  nano vim \
  build-essential libstdc++6 \
  # python
  python3.10 python3.10-dev \
  # needed for templates logic
  kalign \
  # (run continues)
  # stop git from complaining about dubious ownership.
  && git config --global --add safe.directory "*" \
  #
  # cuda softlinking is needed in podman, but not docker
  && ln -s /lib/x86_64-linux-gnu/libcuda.so.1 /lib/x86_64-linux-gnu/libcuda.so \
  && ldconfig /lib/x86_64-linux-gnu/ \
  # setup timezone, to $TZ, ubuntu-specific
  # && ln -fs /usr/share/zoneinfo/$TZ /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  # change default shell to bash (has no effect during building)
  && chsh -s /bin/bash


ENV \
  # expose CUDA libraries. Now that we don't build anything this is likely redundant
  LD_LIBRARY_PATH="/usr/local/cuda/lib64/stubs/:${LD_LIBRARY_PATH:-}" \
  # Set uv timeout to larger value to account for slow download time of nvidia-cudnn-cu12
  UV_HTTP_TIMEOUT=1000 \
  # where virtual env will be installed
  VIRTUAL_ENV=/opt/venv

# Install dependencies in virtualenv
COPY ./requirements.in /tmp/requirements.in
# from https://pythonspeed.com/articles/activate-virtualenv-dockerfile/
# a trick to have virtualenv "always activated"
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN --mount=type=cache,target=/root/.cache/uv \
  # Install uv
  curl -LsSf https://astral.sh/uv/0.5.4/install.sh | sh \
  && . $HOME/.local/bin/env \
  && uv venv --no-python-downloads $VIRTUAL_ENV \
  # this is sh, not bash, so . not source
  && . $VIRTUAL_ENV/bin/activate \
  && uv pip install uv pip -r /tmp/requirements.in


# making sure envvars are set in all shells
RUN echo "PATH=\"$PATH\"" >> /etc/environment \
  && echo "LANG=\"$LANG\"" >> /etc/environment \
  && echo "LC_ALL=\"$LC_ALL\"" >> /etc/environment \
  && echo "LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH\"" >> /etc/environment \
  && echo "EDITOR=\"$EDITOR\"" >> /etc/environment

# no startup command. 