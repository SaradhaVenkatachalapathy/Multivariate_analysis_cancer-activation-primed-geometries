path_to_programs="E:\\co_culture\\programs\\"
path_to_image_folder="E:\\co_culture\\"
runMacro(path_to_programs+"Binary.ijm",path_to_image_folder)
runMacro(path_to_programs+"3dcrop.ijm",path_to_image_folder+" 4")
runMacro(path_to_programs+"measures.ijm",path_to_image_folder+" 1 3 MKL")
runMacro(path_to_programs+"zproject_2d_measure.ijm",path_to_image_folder)

runMacro(path_to_programs+"compaction_zproject_nucleus.ijm",path_to_image_folder)
runMacro(path_to_programs+"asma_per_nuc.ijm",path_to_image_folder+" 4")
runMacro(path_to_programs+"asma_per_nuc.ijm",asma_binary_positive+" 4 300")

runMacro(path_to_programs+"sorted.ijm",path_to_image_folder)

runMacro(path_to_programs+"hc_node_id.ijm",path_to_image_folder + " indivisual_nuclei_ch2")
runMacro(path_to_programs+"hc_measure.ijm",path_to_image_folder+ " MKL")
runMacro(path_to_programs+"hc_distances.ijm",path_to_image_folder+ " MKL 320")
