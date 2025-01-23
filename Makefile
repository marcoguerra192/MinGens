# Default target: build the Cython module
all: build

# Build the Cython module using setup.py
build:
	@echo "Building Cython module..."
	cd src && python setup.py build_ext --inplace

# Install packages
install:
	@echo "Installing Python dependencies..."
	pip install -r requirements.txt
