#!/usr/bin/env python3
#

## Is the molecule/system neutral ?
# If not then add ions
# Else exit with simple no
## Run packmol script to add ions
#First check the charge of molecule/system
# Run in subprocess $PACKMOL <ion_box.inp >ion_box.out


```
    Check the syscharge file in the structre dir. If neutral molecule then don't add charge else check if positive
    or negative and add correct counter ion (CLA/POT) or ask for counter-ion name.
```

# Imports
import os
import subprocess

input_file = "../meso/syscharge.dat"

class System:
    def __init__(
        self,
        syscharge=0, # Assume neutral
    )

    def get_charge()
