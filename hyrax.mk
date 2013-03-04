hyrax_INC_DIRS := $(shell find $(HYRAX_DIR)/include -type d -not -path "*/.svn*")
hyrax_INCLUDE  := $(foreach i, $(hyrax_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(hyrax_INCLUDE) $(libmesh_INCLUDE)

hyrax_LIB := $(HYRAX_DIR)/libhyrax-$(METHOD).la
LIBS += $(hyrax_LIB)

hyrax_APP := $(HYRAX_DIR)/hyrax-$(METHOD)

# source files
hyrax_srcfiles    := $(shell find $(HYRAX_DIR)/src -name "*.C" -not -name main.C)
hyrax_csrcfiles   := $(shell find $(HYRAX_DIR)/src -name "*.c")
hyrax_fsrcfiles   := $(shell find $(HYRAX_DIR)/src -name "*.f")
hyrax_f90srcfiles := $(shell find $(HYRAX_DIR)/src -name "*.f90")

# object files
hyrax_objects := $(patsubst %.C, %.$(obj-suffix), $(hyrax_srcfiles))
hyrax_objects += $(patsubst %.c, %.$(obj-suffix), $(hyrax_csrcfiles))
hyrax_objects += $(patsubst %.f, %.$(obj-suffix), $(hyrax_fsrcfiles))
hyrax_objects += $(patsubst %.f90, %.$(obj-suffix), $(hyrax_f90srcfiles))

# plugin files
hyrax_plugfiles    := $(shell find $(HYRAX_DIR)/plugins/ -name "*.C" 2>/dev/null)
hyrax_cplugfiles   := $(shell find $(HYRAX_DIR)/plugins/ -name "*.c" 2>/dev/null)
hyrax_fplugfiles   := $(shell find $(HYRAX_DIR)/plugins/ -name "*.f" 2>/dev/null)
hyrax_f90plugfiles := $(shell find $(HYRAX_DIR)/plugins/ -name "*.f90" 2>/dev/null)

# plugins
hyrax_plugins := $(patsubst %.C, %-$(METHOD).plugin, $(hyrax_plugfiles))
hyrax_plugins += $(patsubst %.c, %-$(METHOD).plugin, $(hyrax_cplugfiles))
hyrax_plugins += $(patsubst %.f, %-$(METHOD).plugin, $(hyrax_fplugfiles))
hyrax_plugins += $(patsubst %.f90, %-$(METHOD).plugin, $(hyrax_f90plugfiles))

# hyrax main
hyrax_main_src    := $(HYRAX_DIR)/src/main.C
hyrax_app_objects := $(patsubst %.C, %.$(obj-suffix), $(hyrax_main_src))

# dependency files
hyrax_deps := $(patsubst %.C, %.$(obj-suffix).d, $(hyrax_srcfiles)) \
              $(patsubst %.c, %.$(obj-suffix).d, $(hyrax_csrcfiles)) \
              $(patsubst %.C, %.$(obj-suffix).d, $(hyrax_main_src))

# If building shared libs, make the plugins a dependency, otherwise don't.
ifeq ($(libmesh_shared),yes)
  hyrax_plugin_deps := $(hyrax_plugins)
else
  hyrax_plugin_deps :=
endif

all:: $(hyrax_LIB)

$(hyrax_LIB): $(hyrax_objects) $(hyrax_plugin_deps)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(hyrax_objects) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(HYRAX_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(hyrax_LIB) $(HYRAX_DIR)

# include HYRAX dep files
-include $(hyrax_deps)


# how to build HYRAX application
ifeq ($(APPLICATION_NAME),hyrax)
all:: hyrax

hyrax: $(hyrax_APP)

$(hyrax_APP): $(moose_LIB) $(elk_MODULES) $(hyrax_LIB) $(hyrax_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
          $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(hyrax_app_objects) $(hyrax_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

endif

#
# Maintenance
#
delete_list := $(hyrax_APP) $(hyrax_LIB) $(HYRAX_DIR)/libhyrax-$(METHOD).*

cleanall:: 
	make -C $(HYRAX_DIR) clean 

###############################################################################
# Additional special case targets should be added here
