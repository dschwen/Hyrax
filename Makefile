###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Required Environment variables
# LIBMESH_DIR	- location of the libMesh library
#
# Required Make variables
# APPLICATION_NAME  - the name of this application (all lower case)
# MOOSE_DIR	        - location of the MOOSE framework
#
# Optional Environment variables
# CURR_DIR	- current directory (DO NOT MODIFY THIS VARIABLE)
#
# Note: Make sure that there is no whitespace after the word 'yes' if enabling
# an application
###############################################################################
CURR_DIR	?= $(shell pwd)
ROOT_DIR        ?= $(shell dirname `pwd`)
LIBMESH_DIR     ?= $(ROOT_DIR)/libmesh
MOOSE_DIR	?= $(ROOT_DIR)/moose
ELK_DIR         ?= $(ROOT_DIR)/elk
HYRAX_DIR       ?= $(ROOT_DIR)/hyrax

MAKE_LIBRARY := no
APPLICATION_NAME := hyrax

DEP_APPS    ?= $(shell $(MOOSE_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))

include $(MOOSE_DIR)/build.mk
#dependencies
include $(MOOSE_DIR)/moose.mk
include $(ELK_DIR)/elk.mk
include $(HYRAX_DIR)/hyrax.mk

###############################################################################
# Additional special case targets should be added here
