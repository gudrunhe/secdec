SECTOR%(sector_index)i_CPP = \
	%(sector_cpp_files)s
SECTOR%(sector_index)i_DISTSRC = \
	%(sector_distsrc_files)s
SECTOR%(sector_index)i_MMA = \
	%(sector_mma_files)s
SECTOR_CPP += $(SECTOR%(sector_index)i_CPP)
SECTOR_MMA += $(SECTOR%(sector_index)i_MMA)

$(SECTOR%(sector_index)i_DISTSRC) $(SECTOR%(sector_index)i_CPP) $(patsubst %%.cpp,%%.hpp,$(SECTOR%(sector_index)i_CPP)) : codegen/sector%(sector_index)i.done ;
$(SECTOR%(sector_index)i_MMA) : codegen/sector%(sector_index)i.mma.done ;
