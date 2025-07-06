# Configurations for the TATSP runner

# Max solver times
MAXTIME_LB_LP = 1010
MAXTIME_LB_RLXLAG = 1020
MAXTIME_LB_COLGEN = 1030
MAXTIME_UB_LP = 1040
MAXTIME_UB_RLXLAG = 1050
MAXTIME_UB_COLGEN = 1060
MAXTIME_ILP = 1070

# Other parameters
SEED_NUMBER = 1234
RA_PARAM = 291207
VERBOSE_MODE = --verbose

# Choose between directories to change instances
INSTANCES_DIR = instances/instances_release_1
# INSTANCES_DIR = instances/instances_release_2

# Toy instance for quick testing
TEST_INSTANCE := instances/inst-slide8.txt
