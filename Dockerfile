FROM sagemath/sagemath:latest
WORKDIR /app
COPY . /app
RUN sage -python -m pip install --no-cache-dir \
    numpy \
    scipy \
    matplotlib
CMD ["sage", "-python", "src/main.py"]