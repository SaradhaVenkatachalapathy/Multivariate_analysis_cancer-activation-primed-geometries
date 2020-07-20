dirsa = getArgument();
dirres= dirsa + "2D measures"+ File.separator;
File.makeDirectory(dirres); 
setBatchMode(true);

dir1= dirsa + "indivisual_nuclei_ch1" + File.separator;
list1 = getFileList(dir1);

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i];
	open(path);
	run("8-bit");
//Set Properties based on the microscope objective and the camera specifications//
run("Properties...", "unit=micron");
//Z project//
run("Z Project...", " projection=[Sum Slices]");
// Set the threshold //
	run("8-bit");
setThreshold(1, 255);

// set the measurements to be made//
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");
//Run the analysis//
	run("Measure");
	close();
}
run("Close All");
saveAs("Results",  dirres + "2D_indivisual_nuclei_ch1.csv"); 
run("Clear Results");

dir1= dirsa + "indivisual_nuclei_ch2" + File.separator;
list1 = getFileList(dir1);

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i];
	open(path);
	run("8-bit");
//Set Properties based on the microscope objective and the camera specifications//
run("Properties...", "unit=micron");
//Z project//
run("Z Project...", " projection=[Sum Slices]");
// Set the threshold //
	run("8-bit");
setThreshold(1, 255);
// set the measurements to be made//
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");
//Run the analysis//
	run("Measure");
	close();
}
run("Close All");
saveAs("Results",  dirres + "2D_indivisual_nuclei_ch2.csv"); 
run("Clear Results");

dir1= dirsa + "indivisual_nuclei_ch3" + File.separator;
list1 = getFileList(dir1);

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i];
	open(path);
	run("8-bit");
//Set Properties based on the microscope objective and the camera specifications//
run("Properties...", "unit=micron");
//Z project//
run("Z Project...", " projection=[Sum Slices]");
// Set the threshold //
	run("8-bit");
setThreshold(1, 255);
// set the measurements to be made//
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");
//Run the analysis//
	run("Measure");
	close();
}
run("Close All");
saveAs("Results",  dirres + "2D_indivisual_nuclei_ch3.csv"); 
run("Clear Results");

dir1= dirsa + "indivisual_nuclei_ch4" + File.separator;
list1 = getFileList(dir1);

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i];
	open(path);
	run("8-bit");
//Set Properties based on the microscope objective and the camera specifications//
run("Properties...", "unit=micron");
//Z project//
run("Z Project...", " projection=[Sum Slices]");
// Set the threshold //
	run("8-bit");
setThreshold(1, 255);
// set the measurements to be made//
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");
//Run the analysis//
	run("Measure");
	close();
}
run("Close All");
saveAs("Results",  dirres + "2D_indivisual_nuclei_ch4.csv"); 
run("Clear Results");