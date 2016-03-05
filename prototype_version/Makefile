# Makefile for building project program from the CARLsim library

# NOTE: if you are compiling your code in a directory different from
# examples/<example_name> or projects/<project_name> then you need to either
# move the configured user.mk file to this directory or set the path to
# where CARLsim can find the user.mk.
USER_MK_PATH = /CompNeuro/CARLsim-master/
include $(USER_MK_PATH)user.mk

project := mem_mod_gpu
output := *.dot *.dat *.log *.csv

# You should not need to edit the file beyond this point
# ------------------------------------------------------

# we are compiling from lib
CARLSIM_FLAGS += -I$(CARLSIM_LIB_DIR)/include/kernel \
				 -I$(CARLSIM_LIB_DIR)/include/interface \
				 -I$(CARLSIM_LIB_DIR)/include/spike_monitor \
				 -I$(CARLSIM_LIB_DIR)/include/connection_monitor \
				 -I$(CARLSIM_LIB_DIR)/include/spike_generators \
				 -I$(CARLSIM_LIB_DIR)/include/visual_stimulus \
				 -I$(CARLSIM_LIB_DIR)/include/simple_weight_tuner \
				 -I$(CARLSIM_LIB_DIR)/include/stopwatch \
				 -I$(CARLSIM_LIB_DIR)/include/group_monitor
CARLSIM_LIBS  += -L$(CARLSIM_LIB_DIR)/lib -lCARLsim

local_src  := main_$(project).cpp
local_prog := $(project)

# you can add your own local objects
local_objs :=

output_files += $(local_prog) $(local_objs)

.PHONY: clean distclean devtest
# compile from CARLsim lib
$(local_prog): $(local_src) $(local_objs)
	$(NVCC) $(CARLSIM_INCLUDES) $(CARLSIM_FLAGS) $(local_objs) $< -o $@ $(CARLSIM_LFLAGS) $(CARLSIM_LIBS)
	
all: $(local_src) $(local_objs)
	$(NVCC) $(CARLSIM_INCLUDES) $(CARLSIM_FLAGS) $(local_objs) $< -o $@ $(CARLSIM_LFLAGS) $(CARLSIM_LIBS)

clean:
	$(RM) $(output_files)

distclean:
	$(RM) $(output_files) results/*

devtest:
	@echo $(CARLSIM_FLAGS)
