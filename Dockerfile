# Use an official Python slim image
FROM python:3.9-slim

# Install gfortran and build-essential for compiling Fortran code
RUN apt-get update && \
    apt-get install -y gfortran build-essential && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory to /app (our project root)
WORKDIR /app

# Copy the entire project into the container
COPY . /app

# Compile the Fortran code into librefprop.so
RUN gfortran -O2 -fPIC -shared -fno-underscoring \
  -I./FLUIDS \
  -I./FORTRAN \
  ./FORTRAN/DLLFILES/*.FOR \
  ./FORTRAN/CORE_ANC.FOR ./FORTRAN/CORE_FEQ.FOR ./FORTRAN/CORE_PR.FOR ./FORTRAN/FLSH_SUB.FOR \
  ./FORTRAN/MIX_HMX.FOR ./FORTRAN/PROP_SUB.FOR ./FORTRAN/REFPROP.FOR ./FORTRAN/SAT_SUB.FOR ./FORTRAN/SETUP.FOR \
  ./FORTRAN/TRNS_TCX.FOR ./FORTRAN/TRNS_VIS.FOR ./FORTRAN/TRNSP.FOR ./FORTRAN/UTILITY.FOR \
  -o librefprop.so

# Install Python dependencies
RUN pip install --no-cache-dir flask flask_cors numpy future

# Set environment variables that your app.py expects.
ENV REFPROP_ROOT=/app
ENV PYTHONPATH=/app

# Create any necessary directories
RUN mkdir -p API/utils API/endpoints API/models

# Expose the port used by Flask
EXPOSE 5000

# Run the Flask API when the container starts
CMD ["python", "API/app.py"]