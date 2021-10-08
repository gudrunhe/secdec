SECTOR%(sector_index)i_CPP = \
	%(sector_cpp_files)s
SECTOR%(sector_index)i_DISTSRC = \
	%(sector_distsrc_files)s
SECTOR_CPP += $(SECTOR%(sector_index)i_CPP)

$(SECTOR%(sector_index)i_DISTSRC) $(SECTOR%(sector_index)i_CPP) $(patsubst %%.cpp,%%.hpp,$(SECTOR%(sector_index)i_CPP)) : codegen/sector%(sector_index)i.done ;
