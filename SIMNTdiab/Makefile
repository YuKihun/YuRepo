# Wind River Workbench generated Makefile.
# Do not edit!!!
#
# The file ".wrmakefile" is the template used by the Wind River Workbench to
# generate the makefiles of this project. Add user-specific build targets and
# make rules only(!) in this project's ".wrmakefile" file. These will then be
# automatically dumped into the makefiles.

WIND_HOME := $(subst \,/,$(WIND_HOME))
WIND_BASE := $(subst \,/,$(WIND_BASE))
WIND_USR := $(subst \,/,$(WIND_USR))
WRVX_COMPBASE := $(subst \,/,$(WRVX_COMPBASE))

all : pre_build main_all post_build

_clean ::
	@echo "make: removing targets and objects of `pwd`"

TRACE=0
TRACEON=$(TRACE:0=@)
TRACE_FLAG=$(TRACEON:1=)

MAKEFILE := Makefile

FLEXIBLE_BUILD := 1

BUILD_SPEC = SIMNTdiab
DEBUG_MODE = 1
ifeq ($(DEBUG_MODE),1)
MODE_DIR := Debug
else
MODE_DIR := NonDebug
endif
OBJ_DIR := .
WS_ROOT_DIR := D:/WindRiver_workspace
PRJ_ROOT_DIR := $(WS_ROOT_DIR)/Vincentys_Formulae



#Global Build Macros
PROJECT_TYPE = DKM
DEFINES = 
EXPAND_DBG = 0


#BuildSpec specific Build Macros
VX_CPU_FAMILY = simpc
CPU = SIMNT
TOOL_FAMILY = diab
TOOL = diab
TOOL_PATH = 
CC_ARCH_SPEC = -tX86LH:vxworks69
VSB_DIR = $(WIND_BASE)/target/lib
VSB_CONFIG_FILE = $(VSB_DIR)/h/config/vsbConfig.h
LIBPATH = 
LIBS = 

IDE_INCLUDES = -I$(WIND_BASE)/target/h -I$(WIND_BASE)/target/h/wrn/coreip 

IDE_LIBRARIES = 

IDE_DEFINES = -DCPU=_VX_$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL -D_VSB_CONFIG_FILE=\"$(VSB_DIR)/h/config/vsbConfig.h\"   -DIP_PORT_VXWORKS=69 



#BuildTool flags
ifeq ($(DEBUG_MODE),1)
DEBUGFLAGS_C-Compiler = -g
DEBUGFLAGS_C++-Compiler = -g
DEBUGFLAGS_Linker = -g
DEBUGFLAGS_Partial-Image-Linker = 
DEBUGFLAGS_Librarian = 
DEBUGFLAGS_Assembler = -g
else
DEBUGFLAGS_C-Compiler =  -XO -Xsize-opt 
DEBUGFLAGS_C++-Compiler =  -XO -Xsize-opt 
DEBUGFLAGS_Linker =  -XO -Xsize-opt 
DEBUGFLAGS_Partial-Image-Linker = 
DEBUGFLAGS_Librarian = 
DEBUGFLAGS_Assembler =  -XO -Xsize-opt 
endif


#Project Targets
PROJECT_TARGETS = Vincontys_Formulae/$(MODE_DIR)/Vincontys_Formulae.out \
	 Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage.o


#Rules

# Vincontys_Formulae
ifeq ($(DEBUG_MODE),1)
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_C-Compiler = -g
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_C++-Compiler = -g
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Linker = -g
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Partial-Image-Linker = 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Librarian = 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Assembler = -g
else
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_C-Compiler =  -XO -Xsize-opt 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_C++-Compiler =  -XO -Xsize-opt 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Linker =  -XO -Xsize-opt 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Partial-Image-Linker = 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Librarian = 
Vincontys_Formulae/$(MODE_DIR)/% : DEBUGFLAGS_Assembler =  -XO -Xsize-opt 
endif
Vincontys_Formulae/$(MODE_DIR)/% : IDE_INCLUDES = -I$(WIND_BASE)/target/h -I$(WIND_BASE)/target/h/wrn/coreip 
Vincontys_Formulae/$(MODE_DIR)/% : IDE_LIBRARIES = 
Vincontys_Formulae/$(MODE_DIR)/% : IDE_DEFINES = -DCPU=_VX_$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL -D_VSB_CONFIG_FILE=\"$(VSB_DIR)/h/config/vsbConfig.h\"   -DIP_PORT_VXWORKS=69 
Vincontys_Formulae/$(MODE_DIR)/% : PROJECT_TYPE = DKM
Vincontys_Formulae/$(MODE_DIR)/% : DEFINES = 
Vincontys_Formulae/$(MODE_DIR)/% : EXPAND_DBG = 0
Vincontys_Formulae/$(MODE_DIR)/% : VX_CPU_FAMILY = simpc
Vincontys_Formulae/$(MODE_DIR)/% : CPU = SIMNT
Vincontys_Formulae/$(MODE_DIR)/% : TOOL_FAMILY = diab
Vincontys_Formulae/$(MODE_DIR)/% : TOOL = diab
Vincontys_Formulae/$(MODE_DIR)/% : TOOL_PATH = 
Vincontys_Formulae/$(MODE_DIR)/% : CC_ARCH_SPEC = -tX86LH:vxworks69
Vincontys_Formulae/$(MODE_DIR)/% : VSB_DIR = $(WIND_BASE)/target/lib
Vincontys_Formulae/$(MODE_DIR)/% : VSB_CONFIG_FILE = $(VSB_DIR)/h/config/vsbConfig.h
Vincontys_Formulae/$(MODE_DIR)/% : LIBPATH = 
Vincontys_Formulae/$(MODE_DIR)/% : LIBS = 
Vincontys_Formulae/$(MODE_DIR)/% : OBJ_DIR := Vincontys_Formulae/$(MODE_DIR)

OBJECTS_Vincontys_Formulae = Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage.o

Vincontys_Formulae/$(MODE_DIR)/Vincontys_Formulae.out : $(OBJECTS_Vincontys_Formulae)
	$(TRACE_FLAG)if [ ! -d "`dirname "$@"`" ]; then mkdir -p "`dirname "$@"`"; fi;echo "building $@";rm -f "$@";ddump -Ng $(OBJECTS_Vincontys_Formulae) | tclsh $(WIND_BASE)/host/resource/hutils/tcl/munch.tcl -c pentium -tags $(VSB_DIR)/tags/simpc/SIMNT/common/dkm.tags > $(OBJ_DIR)/ctdt.c; $(TOOL_PATH)dcc $(DEBUGFLAGS_Linker) $(CC_ARCH_SPEC) -Xdollar-in-ident -ei1518,4177,4223,4301,4550,5409,1606 -ei4193,4826,4381,4237,1573,4007,4082,4177,4223,4260,4550,5361,5828 -ei1522,4092,4111,4144,4152,4167,4171,4174,4186,4188,4191,4192,4223,4231,4236,4284,4375,4494,4513,5152,5457 -Xforce-declarations  $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES)  $(IDE_DEFINES) $(DEFINES) -o $(OBJ_DIR)/ctdt.o -c $(OBJ_DIR)/ctdt.c; $(TOOL_PATH)dld -tX86LH:vxworks69 -X -r5 -f 0x90,1,1 -r4 -o "$@" $(OBJ_DIR)/ctdt.o $(OBJECTS_Vincontys_Formulae) $(IDE_LIBRARIES) $(LIBPATH) $(LIBS) $(ADDED_LIBPATH) $(ADDED_LIBS) && if [ "$(EXPAND_DBG)" = "1" ]; then plink "$@";fi

Vincontys_Formulae/$(MODE_DIR)/Vincontys_Formulae_compile_file : $(FILE) ;

_clean :: Vincontys_Formulae/$(MODE_DIR)/Vincontys_Formulae_clean

Vincontys_Formulae/$(MODE_DIR)/Vincontys_Formulae_clean : 
	$(TRACE_FLAG)if [ -d "Vincontys_Formulae" ]; then cd "Vincontys_Formulae"; rm -rf $(MODE_DIR); fi


# Vincontys_Formulae_partialImage
ifeq ($(DEBUG_MODE),1)
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_C-Compiler = -g
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_C++-Compiler = -g
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Linker = -g
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Partial-Image-Linker = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Librarian = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Assembler = -g
else
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_C-Compiler =  -XO -Xsize-opt 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_C++-Compiler =  -XO -Xsize-opt 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Linker =  -XO -Xsize-opt 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Partial-Image-Linker = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Librarian = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEBUGFLAGS_Assembler =  -XO -Xsize-opt 
endif
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : IDE_INCLUDES = -I$(WIND_BASE)/target/h -I$(WIND_BASE)/target/h/wrn/coreip 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : IDE_LIBRARIES = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : IDE_DEFINES = -DCPU=_VX_$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL -D_VSB_CONFIG_FILE=\"$(VSB_DIR)/h/config/vsbConfig.h\"   -DIP_PORT_VXWORKS=69 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : PROJECT_TYPE = DKM
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : DEFINES = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : EXPAND_DBG = 0
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : VX_CPU_FAMILY = simpc
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : CPU = SIMNT
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : TOOL_FAMILY = diab
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : TOOL = diab
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : TOOL_PATH = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : CC_ARCH_SPEC = -tX86LH:vxworks69
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : VSB_DIR = $(WIND_BASE)/target/lib
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : VSB_CONFIG_FILE = $(VSB_DIR)/h/config/vsbConfig.h
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : LIBPATH = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : LIBS = 
Vincontys_Formulae_partialImage/$(MODE_DIR)/% : OBJ_DIR := Vincontys_Formulae_partialImage/$(MODE_DIR)

Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CDirectMethod.o : $(PRJ_ROOT_DIR)/Src/CDirectMethod.cpp
	$(TRACE_FLAG)if [ ! -d "`dirname "$@"`" ]; then mkdir -p "`dirname "$@"`"; fi;echo "building $@"; $(TOOL_PATH)dcc $(DEBUGFLAGS_C++-Compiler) $(CC_ARCH_SPEC) -W:c:,-Xclib-optim-off -Xansi -Xlocal-data-area-static-only  -W:c++:.CPP -ei1518,4177,4223,4301,4550,5409,1606 -ei4193,4826,4381,4237,1573,4007,4082,4177,4223,4260,4550,5361,5828 -ei1522,4092,4111,4144,4152,4167,4171,4174,4186,4188,4191,4192,4223,4231,4236,4284,4375,4494,4513,5152,5457 -Xforce-declarations   -Xmake-dependency=0xd $(IDE_DEFINES) $(DEFINES) $(ADDED_C++FLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -o "$@" -c "$<"


Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CInverseMethod.o : $(PRJ_ROOT_DIR)/Src/CInverseMethod.cpp
	$(TRACE_FLAG)if [ ! -d "`dirname "$@"`" ]; then mkdir -p "`dirname "$@"`"; fi;echo "building $@"; $(TOOL_PATH)dcc $(DEBUGFLAGS_C++-Compiler) $(CC_ARCH_SPEC) -W:c:,-Xclib-optim-off -Xansi -Xlocal-data-area-static-only  -W:c++:.CPP -ei1518,4177,4223,4301,4550,5409,1606 -ei4193,4826,4381,4237,1573,4007,4082,4177,4223,4260,4550,5361,5828 -ei1522,4092,4111,4144,4152,4167,4171,4174,4186,4188,4191,4192,4223,4231,4236,4284,4375,4494,4513,5152,5457 -Xforce-declarations   -Xmake-dependency=0xd $(IDE_DEFINES) $(DEFINES) $(ADDED_C++FLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -o "$@" -c "$<"


Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/main.o : $(PRJ_ROOT_DIR)/Src/main.cpp
	$(TRACE_FLAG)if [ ! -d "`dirname "$@"`" ]; then mkdir -p "`dirname "$@"`"; fi;echo "building $@"; $(TOOL_PATH)dcc $(DEBUGFLAGS_C++-Compiler) $(CC_ARCH_SPEC) -W:c:,-Xclib-optim-off -Xansi -Xlocal-data-area-static-only  -W:c++:.CPP -ei1518,4177,4223,4301,4550,5409,1606 -ei4193,4826,4381,4237,1573,4007,4082,4177,4223,4260,4550,5361,5828 -ei1522,4092,4111,4144,4152,4167,4171,4174,4186,4188,4191,4192,4223,4231,4236,4284,4375,4494,4513,5152,5457 -Xforce-declarations   -Xmake-dependency=0xd $(IDE_DEFINES) $(DEFINES) $(ADDED_C++FLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -o "$@" -c "$<"


OBJECTS_Vincontys_Formulae_partialImage = Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CDirectMethod.o \
	 Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CInverseMethod.o \
	 Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/main.o

Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage.o : $(OBJECTS_Vincontys_Formulae_partialImage)
	$(TRACE_FLAG)if [ ! -d "`dirname "$@"`" ]; then mkdir -p "`dirname "$@"`"; fi;echo "building $@"; $(TOOL_PATH)dld -tX86LH:vxworks69 -X -r5 -f 0x90,1,1  -o "$@" $(OBJECTS_Vincontys_Formulae_partialImage) $(ADDED_OBJECTS) $(IDE_LIBRARIES) $(LIBPATH) $(LIBS) $(ADDED_LIBPATH) $(ADDED_LIBS) && if [ "$(EXPAND_DBG)" = "1" ]; then plink "$@";fi

Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage_compile_file : $(FILE) ;

_clean :: Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage_clean

Vincontys_Formulae_partialImage/$(MODE_DIR)/Vincontys_Formulae_partialImage_clean : 
	$(TRACE_FLAG)if [ -d "Vincontys_Formulae_partialImage" ]; then cd "Vincontys_Formulae_partialImage"; rm -rf $(MODE_DIR); fi

DEP_FILES := Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CDirectMethod.d Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/CInverseMethod.d Vincontys_Formulae_partialImage/$(MODE_DIR)/Objects/Vincentys_Formulae/Src/main.d
-include $(DEP_FILES)


WIND_SCOPETOOLS_BASE := $(subst \,/,$(WIND_SCOPETOOLS_BASE))

clean_scopetools :
	$(TRACE_FLAG)rm -rf $(PRJ_ROOT_DIR)/.coveragescope/db

CLEAN_STEP := clean_scopetools


-include $(PRJ_ROOT_DIR)/*.makefile

-include *.makefile

TARGET_JOBS_RULE=$(MAKE) -f $(MAKEFILE) --jobs $(TARGET_JOBS) $(MFLAGS) $* TARGET_JOBS=1
ifeq ($(JOBS),1)
main_all : external_build  $(PROJECT_TARGETS)
	@echo "make: built targets of `pwd`"
else
main_all : external_build 
	@$(MAKE) -f $(MAKEFILE) --jobs $(JOBS) $(MFLAGS) $(PROJECT_TARGETS) TARGET_JOBS=1 &&\
	echo "make: built targets of `pwd`"
endif 

# entry point for extending the build
external_build ::
	@echo ""

# main entry point for pre processing prior to the build
pre_build :: $(PRE_BUILD_STEP) generate_sources
	@echo ""

# entry point for generating sources prior to the build
generate_sources ::
	@echo ""

# main entry point for post processing after the build
post_build :: $(POST_BUILD_STEP) deploy_output
	@echo ""

# entry point for deploying output after the build
deploy_output ::
	@echo ""

clean :: external_clean $(CLEAN_STEP) _clean

# entry point for extending the build clean
external_clean ::
	@echo ""
