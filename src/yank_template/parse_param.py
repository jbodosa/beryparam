
from collections import defaultdict

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
        param_dict[str(line[1]).lower()] = line[3]

    #    print(param_dict)
    return(param_dict)
