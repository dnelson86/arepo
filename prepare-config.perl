
# This file processes the configurations options in Config.sh, producing 
# two files:
#
#   arepoconfig.h          to be included in each source file (via allvars.h)
#   compile_time_info.c    code to be compiled in, which will print the configuration 
#
if( @ARGV != 2)
{
    print "usage: perl prepare-config.perl <Config.sh> <build dir>\n";
    exit;
}

open(FILE, @ARGV[0]);
$path = @ARGV[1];


open(OUTFILE, ">${path}/arepoconfig.h");
open(COUTF,   ">${path}/compile_time_info.c");

open(COUTF2,   ">${path}/compile_time_info_hdf5.c");

print COUTF "#include <stdio.h>\n";
print COUTF "void output_compile_time_options(void)\n\{\n";
print COUTF "printf(\n";

print COUTF2 "#include <stdio.h>\n";

print COUTF2 "#include \"arepoconfig.h\"\n";

print COUTF2 "#ifdef HAVE_HDF5\n";
print COUTF2 "#include <hdf5.h>\n\n";

print COUTF2 "hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);\n";
print COUTF2 "hid_t my_H5Screate(H5S_class_t type);\n";
print COUTF2 "herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);\n";
print COUTF2 "herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);\n";
print COUTF2 "herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);\n\n";
print COUTF2 "herr_t my_H5Tclose(hid_t type_id);\n\n";


print COUTF2 "void write_compile_time_options_in_hdf5(hid_t handle)\n\{\n";
print COUTF2 "hid_t hdf5_dataspace, hdf5_attribute;\n";
print COUTF2 "double val;\n";
print COUTF2 "hid_t atype = H5Tcopy(H5T_C_S1);\n";
print COUTF2 "H5Tset_size(atype, 1);\n";


while($line=<FILE>)
{
    chop $line;

    @fields = split ' ' , $line;

    if(substr($fields[0], 0, 1) ne "#")
    {
	if(length($fields[0]) > 0)
	{
	    @subfields = split '=', $fields[0];

	    print OUTFILE "#define $subfields[0] $subfields[1]\n";
            print COUTF   "\"        $fields[0]\\n\"\n";
            
            print COUTF2 "hdf5_dataspace = my_H5Screate(H5S_SCALAR);\n";
            if ($subfields[1] eq "")
            {
                print COUTF2 "hdf5_attribute = my_H5Acreate(handle, \"$subfields[0]\" , atype, hdf5_dataspace, H5P_DEFAULT);\n";
                print COUTF2 "my_H5Awrite(hdf5_attribute, atype, \"\", \"$subfields[0]\");\n";
            }
            else
            {
                print COUTF2 "hdf5_attribute = my_H5Acreate(handle, \"$subfields[0]\" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);\n";
                print COUTF2 "val = $subfields[1];\n";
                print COUTF2 "my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, \"$subfields[0]\");\n";
            }
            print COUTF2 "my_H5Aclose(hdf5_attribute, \"$subfields[0]\");\n";
            print COUTF2 "my_H5Sclose(hdf5_dataspace, H5S_SCALAR);\n\n";
	}
    }
}

print COUTF "\"\\n\");\n";
print COUTF "\}\n";

print COUTF2 "my_H5Tclose(atype);\n";
print COUTF2 "\}\n";
print COUTF2 "#endif\n";
