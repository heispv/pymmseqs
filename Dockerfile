# Build stage for Debian-based image
FROM python:3.10-slim-bookworm AS debian-build
WORKDIR /app/src
COPY . .
COPY examples/ ../examples/
COPY docs/ ../docs/

RUN apt-get update \
    && apt-get install -y curl \
    && pip install --upgrade pip \
    && pip install . \
    && pip cache purge \
    && rm -rf /var/lib/apt/lists/* \
    && cd / \
    && rm -rf /app/src

WORKDIR /app

CMD ["bash"]

# Build stage for Alpine-based image 
FROM python:3.10-alpine AS alpine-build
WORKDIR /app/src
COPY . .
COPY examples/ ../examples/
COPY docs/ ../docs/

RUN apk add --no-cache curl gcc python3-dev musl-dev linux-headers \
    && pip install --upgrade pip \
    && pip install . \
    && pip cache purge \
    && cd / \
    && rm -rf /app/src

WORKDIR /app

CMD ["sh"]
