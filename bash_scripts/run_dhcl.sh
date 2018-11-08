source activate p2.7
rm -f -r ./files/from_dhcl
python ./external_scripts/dhcl/executables/everything.py -d ./files/input_pdb --outdir ./files/from_dhcl
source deactivate
python process_dhcl_output.py

