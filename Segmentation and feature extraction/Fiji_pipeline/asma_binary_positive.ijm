args = split(getArgument()," "); 
dirsa=args[0];
asma_ch=args[1];
threshold_low=args[2];
setBatchMode(true);
//set the names of the directories
dirn=dirsa + "3d objects"+ File.separator;
dirr=dirsa + "rawimages"+ File.separator;

dir_sma_binary=dirsa +"asma_positive_binary"+ File.separator;
File.makeDirectory(dir_sma_binary); 

filenames=getFileList(dirr);

for(i=0; i<filenames.length; i++){

	path=dirr+filenames[i];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=asma_ch c_end=asma_ch c_step=1");
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	
	run("Gaussian Blur...", "sigma=3 stack");
	setThreshold(threshold_low, 4095);
	//waitForUser("ok?");
	run("Make Binary", "method=Default background=Default");
	run("Fill Holes", "stack");
	saveAs("Tiff", dir_sma_binary+baseName+".tiff"); 
	run("Close All");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
}

dirb= dirsa + "3d objects"+ File.separator;
filenames=getFileList(dirb);

Dir_res_asma_positive=dirsa + "Sorting_asma_levels"+ File.separator;
File.makeDirectory(Dir_res_asma_positive); 


for(f=0;f<filenames.length;f++){
	open(dirb+filenames[f]);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tiff"); 
	baseName=substring(obj_image, 0, baseNameEnd); 
	
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	for(p=0;p<nb_obj;p++){
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		name=baseName+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(name);
	}
	Ext.Manager3D_DeselectAll();
	run("Close All");

	path=dir_sma_binary+filenames[f];
	open(path);
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_asma_positive+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");

	Ext.Manager3D_Close();
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
}

