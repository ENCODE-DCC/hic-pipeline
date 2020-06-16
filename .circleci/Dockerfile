# Defines image for running linting
# You should never need to use this directly.
# Builds with docker build -t encodedcc/hic-pipeline:circleci-lint - < Dockerfile
FROM circleci/python:3.7.3-stretch

# Installs rust with cargo clippy and rustfmt
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Install shellcheck
RUN curl -s https://storage.googleapis.com/shellcheck/shellcheck-stable.linux.x86_64.tar.xz | sudo tar xJ && sudo mv shellcheck-stable/shellcheck /usr/local/bin
