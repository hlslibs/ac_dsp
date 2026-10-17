ifeq "$(AC_TYPES_INC)" ""
  $(error Error: The AC_TYPES_INC variable was not set. Cannot locate AC Datatypes installation)
endif

AC_MATH_INC = ../include

ifeq "$(GCC_EXEC)" ""
GCC_EXEC = g++
endif

CXXFLAGS = -g -std=c++11 -I.
LDFLAGS = -s -static-libstdc++
WD = $(shell pwd)
GCOV_ENABLED = false

ifneq "x$(GCOV_HOME)" "x"
GCOV_ENABLED = true
GCOV = $(GCOV_HOME)/gcov
LCOV = $(GCOV_HOME)/lcov
GENHTML = $(GCOV_HOME)/genhtml
CXXFLAGS += -coverage
LDFLAGS += -lgcov
GCOV_OPTS := --gcov-tool $(GCOV) --quiet --rc geninfo_auto_base=1 --rc lcov_branch_coverage=1 --external
GCCVER = $(shell $(GCC_EXEC) -dumpversion)
GCOV_FILT_PAT := '/usr/include/*' '/usr/lib/*' '$(GCCVER)/*' 
endif

DIE        = || exit 1 ;
RM        := rm -f
ECHO      := echo

SOURCES_CPP = \
  rtest_ac_fir_const_coeffs.cpp \
  rtest_ac_fir_load_coeffs.cpp \
  rtest_ac_fir_prog_coeffs.cpp \
  rtest_ac_cic_dec_full.cpp \
  rtest_ac_cic_intr_full.cpp

FFT_R2_OBJS = rtest_ac_fft_r2_dit_inpl.o rtest_ac_fft_r2_dif_inpl.o \
  rtest_ac_fft_r2_dit_sdf.o rtest_ac_fft_r2_dif_sdf.o

OBJS = $(SOURCES_CPP:.cpp=.o) $(FFT_R2_OBJS)

GCOVDAT = $(SOURCES_CPP:.cpp=.gcda) $(SOURCES_CPP:.cpp=.gcno) $(SOURCES_CPP:.cpp=.o.base.info) $(SOURCES_CPP:.cpp=.o.test.info) $(SOURCES_CPP:.cpp=.o.total.info) $(SOURCES_CPP:.cpp=.o.filt.info) 

# Compilation rule
%.o: %.cpp
ifeq "$(GCOV_ENABLED)" "true"
	-@$(RM) *.gcda *.gcno
endif
	-@$(ECHO) "------------------------------ Compile  $< ----------------------------------"
	@$(GCC_EXEC) $(CXXFLAGS) -I$(AC_TYPES_INC) -I$(AC_MATH_INC) $< $(LDFLAGS) $(LINK_LIBNAMES) -o $@ $(DIE)
ifeq "$(GCOV_ENABLED)" "true"
	-@$(ECHO) "------------------------------ Prep     $< ----------------------------------"
	-@$(RM) *.info
	@-$(LCOV) $(GCOV_OPTS) --capture --initial --directory . --output-file $@.base.info
endif
	-@$(ECHO) "------------------------------ Running  $< ----------------------------------"
	@./$@ $(DIE)
ifeq "$(GCOV_ENABLED)" "true"
	-@$(ECHO) "------------------------------ Coverage $< ----------------------------------"
	@-$(LCOV) $(GCOV_OPTS) --capture --directory . --output-file $@.test.info
	@-$(LCOV) $(GCOV_OPTS) --add-tracefile $@.base.info --add-tracefile $@.test.info  --output-file $@.total.info
	@-$(LCOV) $(GCOV_OPTS) --remove $@.total.info $(GCOV_FILT_PAT) -o $@.filt.info
	@-$(GENHTML) --rc geninfo_auto_base=1 --rc lcov_branch_coverage=1 --ignore-errors source $@.filt.info --legend --output-directory=$@-test-lcov
	@-$(ECHO) "HTML coverage report written to: file://$(WD)/$@-test-lcov/index.html"
endif

all: $(OBJS)

# Compile each FFT architecture separately because the twiddle tables are local
# to each core and have include guards.
rtest_ac_fft_r2_dif_inpl.o: FFT_R2_FLAGS = -DTEST_DIF
rtest_ac_fft_r2_dit_sdf.o: FFT_R2_FLAGS = -DTEST_SDF
rtest_ac_fft_r2_dif_sdf.o: FFT_R2_FLAGS = -DTEST_DIF -DTEST_SDF

$(FFT_R2_OBJS): rtest_ac_fft_r2_%.o: rtest_ac_fft_r2.cpp
	-@$(ECHO) "------------------------------ Compile  $@ ----------------------------------"
	@$(GCC_EXEC) $(CXXFLAGS) $(FFT_R2_FLAGS) -I$(AC_TYPES_INC) -I$(AC_MATH_INC) $< $(LDFLAGS) $(LINK_LIBNAMES) -o $@ $(DIE)
	-@$(ECHO) "------------------------------ Running  $@ ----------------------------------"
	@./$@ $(DIE)

.PHONY: clean
clean:
	@-$(RM) -r $(OBJS) $(GCOVDAT)

