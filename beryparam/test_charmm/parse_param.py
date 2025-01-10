import logging
from collections import defaultdict

logger = logging.getLogger(__name__)
######################################
# Read params from sys_param #
######################################

### Parse the charge/params from sys_param.str
def parse_param(sys_param):
    """
    Parse the params from sys_param.str file such as charge, box size, units, counter_ions
    """
    param_dict = defaultdict(str)
    f  = open(sys_param , 'r').readlines()
    for line in f:
        line = line.split()

        if line[0].lower() == "string":
            param_dict[str(line[1]).lower()] = str(line[3])
        elif line[0].lower() == "float":
            param_dict[str(line[1]).lower()] = float(line[3])
        elif line[0].lower() == "int":
            param_dict[str(line[1]).lower()] = int(line[3])
        elif line[0].lower() == "list":
            param_dict[str(line[1]).lower()] = list(line[3:])

    return(param_dict)
