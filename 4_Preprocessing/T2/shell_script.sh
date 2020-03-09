#!/bin/bash

# Install dependencies in userspace
pip install --user sitk

# Launch application with arguments passed to the script
python $@