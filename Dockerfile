FROM python:3.12-slim

COPY --from=ghcr.io/astral-sh/uv:0.9.25 /uv /usr/local/bin/uv

WORKDIR /work

COPY pyproject.toml uv.lock README.md ./
RUN uv sync --locked --no-install-project

COPY src/ ./src/
RUN uv sync --locked

ENV PATH="/work/.venv/bin:$PATH"

COPY . .
