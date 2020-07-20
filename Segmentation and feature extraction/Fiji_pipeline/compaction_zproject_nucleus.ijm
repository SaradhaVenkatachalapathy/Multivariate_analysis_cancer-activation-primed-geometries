dirsa = getArgument();

dir1= dirsa + "indivisual_nuclei_ch1"+ File.separator;
dir_dat= dirsa + "compaction_measure"+ File.separator;
File.makeDirectory(dir_dat); 
setBatchMode(true);
list1 = getFileList(dir1);
n= list1.length;
HCcontent=newArray(n);
ECcontent=newArray(n);
HCvolume=newArray(n);
ECvolume=newArray(n);
volume=newArray(n);
height=newArray(n);
label=newArray(n);
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack limit display redirect=None decimal=3");
for (i=0; i<list1.length; i++) {

		path = dir1+list1[i];
		open(path);
		run("8-bit");
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 255);
		run("Analyze Particles...", "size=10-Infinity circularity=0.00-1.00 show=Nothing display stack");
		num=nResults;
		volume[i]=0;
		height[i]=0;
		nuc_area=newArray(num);
		low_thresh=newArray(num);
		slice_num=newArray(num);
		hcarea=newArray(num);
		tot_intensity=newArray(num);
		hc_intensity=newArray(num);
		
		for(k=0; k<num; k++){
			nuc_area[k]=getResult("Area",k);
			slice_num[k]=getResult("Slice",k);
			low_thresh[k]=getResult("Mean",k) +1.5*getResult("StdDev",k);
			volume[i]=volume[i]+getResult("Area",k)*0.5;
			tot_intensity[k]=getResult("RawIntDen",k);
			if (getResult("Area",k)>0){
				height[i]=height[i]+0.5;
			}
		}
		run("Clear Results");

		for(k=0; k<num; k++){
			setSlice(slice_num[k]);
			//run("Threshold...");
			setThreshold(low_thresh[k], 255);
			run("Measure");
		}
		for(k=0; k<num; k++){
			hcarea[k]=getResult("Area",k);
			hc_intensity[k]=getResult("RawIntDen",k);
		}

		label[i]=File.getName(path);
		HCcontent[i]=0;
		ECcontent[i]=0;
		HCvolume[i]=0;
		ECvolume[i]=0;
		for (j=0; j<nResults; j++){
			HCvolume[i]=HCvolume[i]+hcarea[j];
			ECvolume[i]=ECvolume[i]+(nuc_area[j]-hcarea[j]);
			HCcontent[i]=HCcontent[i]+hc_intensity[j];
			ECcontent[i]=ECcontent[i]+(tot_intensity[j]-hc_intensity[j]);
			
		}
		close();
		HCvolume[i]=HCvolume[i]*0.5;
		ECvolume[i]=ECvolume[i]*0.5;
		run("Clear Results");
}
run("Clear Results");
run("Clear Results");
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
	//run("Threshold...");
	setThreshold(1, 255);
	
// set the measurements to be made//

//Run the analysis//
	run("Measure");
	
	
	close();
run("Close All");
}
run("Close All");

for (k=0; k<list1.length; k++) {
	setResult("Label1",k ,label[k]);
	setResult("HCcontent",k, HCcontent[k]);
	setResult("ECcontent",k, ECcontent[k]);
	setResult("HCvolume",k, HCvolume[k]);
	setResult("ECvolume",k, ECvolume[k]);
	setResult("HC_EC_content",k, HCcontent[k]/ECcontent[k]);
	setResult("HC_EC_volume",k, HCvolume[k]/ECvolume[k]);
	setResult("Volume",k,volume[k]);
	setResult("Height",k,height[k]);
	setResult("DNA density",k,getResult("RawIntDen",k)/volume[k]);
}

saveAs("Results",dir_dat+"compaction.csv" ); 
run("Clear Results");