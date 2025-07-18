# Configurations for the TATSP runner

# Max solver times
MAXTIME_LB_LP = 510
MAXTIME_LB_RLXLAG = 520
MAXTIME_LB_COLGEN = 530
MAXTIME_UB_LP = 540
MAXTIME_UB_RLXLAG = 550
MAXTIME_UB_COLGEN = 560
MAXTIME_ILP = 570

# Other parameters
SEED_NUMBER = 1234
RA_PARAM = 291207
VERBOSE_MODE = --verbose

# Choose between directories to change instances
INSTANCES_DIR = instances/instances_release_1
# INSTANCES_DIR = instances/instances_release_2

# Toy instance for quick testing
# TEST_INSTANCE := instances/inst-slide8.txt
# TEST_INSTANCE := instances/instances_release_1/grf1.txt
# TEST_INSTANCE := instances/instances_release_1/grf8.txt
# TEST_INSTANCE := instances/instances_release_1/grf18.txt
TEST_INSTANCE := instances/instances_release_2/grf101.txt
# TEST_INSTANCE := instances/instances_release_2/grf112.txt
# TEST_INSTANCE := instances/instances_release_2/grf129.txt