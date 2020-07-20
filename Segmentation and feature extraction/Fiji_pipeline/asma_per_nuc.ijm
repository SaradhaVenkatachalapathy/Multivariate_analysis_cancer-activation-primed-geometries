args = split(getArgument()," "); 
dirsa=args[0];
asma_ch=args[1];
setBatchMode(true);
//set the names of the directories
dirr=dirsa + "rawimages"+ File.separator;

dir_sma_binary=dirsa +"asma_levels_per_nucleus"+ File.separator;
dir_numb_of_nuc_field=dir_sma_binary +"number_of_nuclei_per_field"+ File.separator;
dir_amsa_per_field=dir_sma_binary +"asma_per_field"+ File.separator;
File.makeDirectory(dir_sma_binary); 
File.makeDirectory(dir_numb_of_nuc_field); 
File.makeDirectory(dir_amsa_per_field); 

//nucleus
filenames=getFileList(dirr);
numb_nuc=newArray(filenames.length);
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");
for(i=0; i<filenames.length; i++){
	path=dirr+filenames[i];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	run("Z Project...", "projection=[Max Intensity]");
	run("Gaussian Blur...", "sigma=4");
	setAutoThreshold("Otsu dark");
	run("Analyze Particles...", "display include");
	numb_nuc[i]=nResults;
	saveAs("Results",  dir_numb_of_nuc_field + baseName+"nuclear_dat.csv"); 
	run("Close All");
	run("Clear Results");
}
//asma
filenames=getFileList(dirr);
asma_int=newArray(filenames.length);
asma_mean=newArray(filenames.length);
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis display redirect=None decimal=4");
for(i=0; i<filenames.length; i++){
	path=dirr+filenames[i];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=asma_ch c_end=asma_ch c_step=1");
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	run("Z Project...", "projection=[Sum Slices]");
	run("Measure");
	asma_int[i]=getResult("RawIntDen", i);
	asma_mean[i]=getResult("Mean", i);
	run("Close All");
}
saveAs("Results",  dir_amsa_per_field + "asma_dat.csv"); 

run("Clear Results");

for(i=0; i<filenames.length; i++){
	setResult("numb_nuc_per_field", i, numb_nuc[i]);
	setResult("tot_asma_per_field", i, asma_int[i]);
	setResult("tot_asma_per_nuc", i, asma_int[i]/numb_nuc[i]);
	setResult("mean_asma_per_field", i, asma_mean[i]);
	setResult("mean_asma_per_nuc", i, asma_mean[i]/numb_nuc[i]);
}

saveAs("Results",  dir_amsa_per_field + "asma_per_nuc.csv"); 
run("Clear Results");