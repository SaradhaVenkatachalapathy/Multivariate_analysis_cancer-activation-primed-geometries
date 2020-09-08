# Image segmentation and feature extraction 

The first program to run is the _"move_files_rawimages.ijm"_. It moves all files in the input folder to a new folder it creates in the input directory called raw images. 

In order to then identify nuclei as 3D objects, we use _"binarize_the_nucleus.ijm"_ to obtain binary images and _"3dcrop_nuclei.ijm"_ to segment and to obtain individual nuclear crops.

  1. The program _"measures.ijm"_ quantifies the 3D geometrical and intensity measures for nuclei across specified channels 
  2. The program _"zproject_2d_measure.ijm"_ quantifies the geometrical and intensity measures for a z-projected nuclei across channels 
  3. The program _"compaction_zproject_nucleus.ijm"_ allows us to measure chromatin compaction features.

  
  4. The program _"asma_per_nuc.ijm"_ allows us to obtain cytoplasmic aSMA levels as defined by a 2 µm ring from the nuclear edge.
  5. The program _"asma_binary_positive.ijm"_ calculates the aSMA positive cell fraction

  6. The program _"hc_node_id.ijm"_ identifies heterochrtomatin nodes.
  7. The program _"hc_measure.ijm"_ computes the 3D geometric features of heterochromatin nodes
  8. The program _"hc_distances.ijm"_ computes the spactial distribution metrics for position of heterochromatin nodes inside the nucleus. 
 
  9. The program _"sorted.ijm"_ allows the the user to sor the nuclei as "good" or "bad"  
