#!/bin/bash

python convert_charmm2xml.py

python mk_vacuumbox.py
python mk_waterbox.py
#
python mk_setup.py
