/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  MACRO TO IDENTIFY TRANSCRIPTION FACTORS AND POL2 CLUSTERS                											                    					/////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                      /////////////
///////  DATE LAST EDITED: May 1st, 2017                                                                                                                            /////////////
///////  ASSUMPTIONS: The input image is a confocal zstack (slice height=0.5microns) with 3 channels one for nucleus, active pol2, and one transcription factors    /////////////
///////  DESCRIPTION: Two dialog box opens where the user inputs the source (where the images are strored) the results(where results need to be saved) directory	/////////////
///////               The program opens nucleus channel, uses the gaussian filter (sigma=0.15, scaled) to the stack and the thresholds the image in various levels	/////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


dirsa=getArgument();


run("Set Measurements...", "mean standard modal min integrated median skewness area_fraction limit display redirect=None decimal=3");
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
setBatchMode(true);

//set the names of the directories
dir=dirsa + "indivisual_nuclei_ch1"+ File.separator;



//make results folder
dir2 = dirsa + "HC and pol2 hubs" + File.separator;
newDir12 = dir2 + "HC_binary" + File.separator;
newDir5 = dir2 + "datafiles" + File.separator;

File.makeDirectory(dir2); 
File.makeDirectory(newDir12); 

File.makeDirectory(newDir5);


filenames=getFileList(dir);
//filenames=getFileList(dir);
title=newArray(filenames.length);
mean_HC=newArray(filenames.length);
SD_HC=newArray(filenames.length);

run("Set Measurements...", "mean standard modal min display redirect=None decimal=3");
//Identify HC, TF1 and pol2//
for(i=0; i<filenames.length; i++){

//HC
	path=dir+filenames[i];
	open(path);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tif"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	run("8-bit");
	run("Gaussian Blur...", "sigma=0.15 scaled stack");
	title[i]=getTitle(); 
	//get mean and SD
		run("Duplicate...", "duplicate");
		run("Z Project...", "projection=[Max Intensity]");
		setAutoThreshold("Huang dark");
		run("Measure");
		mean_HC[i]=getResult("Mean",0);
		SD_HC[i]=getResult("StdDev",0);
		maxint=getResult("Max",0);
		run("Close");
		run("Clear Results");
	//save mean	+1.5*SD
	selectWindow(title[i]);
	lower_thresh=mean_HC[i]+1.5*SD_HC[i];
	if(lower_thresh<maxint){
	setThreshold(lower_thresh, 255);
	run("Make Binary", "method=Default background=Default");
	saveAs("tiff", newDir12+baseName);		
	}
	run("Close All");

	run("Clear Results");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");

}
print("\\Clear");
//log the mean and SD
	print("label"+","+","+ "mean_HC"+","+"SD_HC");
	for(i=0; i<filenames.length; i++){
	print(title[i]+","+ mean_HC[i]+","+SD_HC[i]);
	}

	selectWindow("Log");  //select Log-window 
	saveAs("Text", newDir5+"_mean_sd.txt"); 
print("\\Clear");
