# Build the SMBH + binary three-body code against SpaceHub.
#   make            -> bin/smbh_binary_PN and bin/smbh_binary_Newt
#   make pn | newt  -> one of them
#   make SPACEHUB=/path/to/SpaceHub/src   (override the auto-detected location)

# SpaceHub location: $SPACEHUB if set, else what env.sh resolves (env.local.sh or
# <repo>/external/SpaceHub/src).
SPACEHUB ?= $(shell bash -c '. "$(CURDIR)/../env.sh" >/dev/null 2>&1; printf %s "$$SPACEHUB"')
SPACEHUB := $(patsubst %/spaceHub.hpp,%,$(SPACEHUB))

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -pthread
# SpaceHub 190777d6 uses uint64_t without <cstdint>; GCC >= 13 no longer includes it implicitly.
CPPFLAGS += -include cstdint
SRC      := smbh_binary.cpp
BIN      := bin

.PHONY: all pn newt clean check
all: pn newt
pn:   $(BIN)/smbh_binary_PN
newt: $(BIN)/smbh_binary_Newt

check:
	@test -n "$(SPACEHUB)" && test -f "$(SPACEHUB)/spaceHub.hpp" || \
	  { echo "SpaceHub not found: set SPACEHUB=/path/to/SpaceHub/src"; exit 1; }
	@echo "SpaceHub: $(SPACEHUB)"

$(BIN)/smbh_binary_PN: $(SRC) | check $(BIN)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(SPACEHUB) $< -o $@

$(BIN)/smbh_binary_Newt: $(SRC) | check $(BIN)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DNEWTONIAN_ONLY -I$(SPACEHUB) $< -o $@

$(BIN):
	mkdir -p $@

clean:
	rm -f $(BIN)/smbh_binary_PN $(BIN)/smbh_binary_Newt
