hyrax_SRC_DIRS := $(HYRAX_DIR)/src/*/*

hyrax_INC_DIRS := $(shell find $(HYRAX_DIR)/include -type d -not -path "*/.svn*")
hyrax_INCLUDE  := $(foreach i, $(hyrax_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(hyrax_INCLUDE) $(libmesh_INCLUDE)

hyrax_LIB := $(HYRAX_DIR)/libhyrax-$(METHOD)$(libext)
LIBS += $(hyrax_LIB)

hyrax_APP := $(HYRAX_DIR)/hyrax-$(METHOD)

# source files
hyrax_srcfiles    := $(shell find $(hyrax_SRC_DIRS) -name *.C)
hyrax_csrcfiles   := $(shell find $(hyrax_SRC_DIRS) -name *.c)
hyrax_fsrcfiles   := $(shell find $(hyrax_SRC_DIRS) -name *.f)
hyrax_f90srcfiles := $(shell find $(hyrax_SRC_DIRS) -name *.f90)
# object files
hyrax_objects	:= $(patsubst %.C, %.$(obj-suffix), $(hyrax_srcfiles))
hyrax_objects	+= $(patsubst %.c, %.$(obj-suffix), $(hyrax_csrcfiles))
hyrax_objects += $(patsubst %.f, %.$(obj-suffix), $(hyrax_fsrcfiles))
hyrax_objects += $(patsubst %.f90, %.$(obj-suffix), $(hyrax_f90srcfiles))

hyrax_app_objects := $(patsubst %.C, %.$(obj-suffix), $(HYRAX_DIR)/src/main.C)

all:: $(hyrax_LIB)

# build rule for lib HYRAX
ifeq ($(enable-shared),yes)
# Build dynamic library
$(hyrax_LIB): $(hyrax_objects)
	@echo "Linking "$@"..."
	@$(libmesh_CC) $(libmesh_CXXSHAREDFLAG) -o $@ $(hyrax_objects) $(libmesh_LDFLAGS)
else
# Build static library
ifeq ($(findstring darwin,$(hostos)),darwin)
$(hyrax_LIB): $(hyrax_objects)
	@echo "Linking "$@"..."
	@libtool -static -o $@ $(hyrax_objects)
else
$(hyrax_LIB): $(hyrax_objects)
	@echo "Linking "$@"..."
	@$(AR) rv $@ $(hyrax_objects)
endif
endif

# include HYRAX dep files
-include $(HYRAX_DIR)/src/*/*.d


# how to build HYRAX application
ifeq ($(APPLICATION_NAME),hyrax)
all:: hyrax

hyrax: $(hyrax_APP)

$(hyrax_APP): $(moose_LIB) $(elk_MODULES) $(hyrax_LIB) $(hyrax_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(hyrax_app_objects) -o $@ $(hyrax_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

-include $(HYRAX_DIR)/src/*.d
endif

#
# Maintenance
#
delete_list := $(hyrax_APP) $(hyrax_LIB)

cleanall:: 
	make -C $(HYRAX_DIR) clean 

###############################################################################
# Additional special case targets should be added here
