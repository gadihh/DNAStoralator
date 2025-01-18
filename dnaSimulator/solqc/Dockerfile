# Use an official Python runtime as a parent image
FROM python:3.7

# Set the working directory to /app
WORKDIR /solqc

# Copy the current directory contents into the container at /app
COPY . /solqc

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r ./requirements.txt

# Create an executable from qc_server
RUN pyinstaller --onefile qc_server.py

# Remove source code from image.
RUN rm -r src
RUN rm -r qc_server.py
RUN rm -r main.py
RUN rm -r requirements.txt
RUN mkdir data

# Run qc server when the container launches
CMD ["./dist/qc_server"]