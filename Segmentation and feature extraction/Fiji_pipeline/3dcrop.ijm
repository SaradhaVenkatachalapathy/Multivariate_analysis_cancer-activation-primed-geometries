/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO IDENTIFY SEGMENT NUCLEI IN 3D AND MEASURE THE INTENSITY IN CHANNEL 1 AND 2         											                    /////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   /////////////
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus 
///////  DESCRIPTION: Two dialog box opens where the user inputs the source (where the folder containing images are strored). The program creates 2 subfolders where
///////               cropped 3D nuclei are stored. The program opens nucleus channel, smoothens the stack and the thresholds the image and identifies objects in 3D. 
///////               The program then opens the raw image and crops the image in the all channels individually. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 
dirsa = args[0];
nchannels=args[1];
//dirsa=getDirectory("Please choose the source directory");
dirw= dirsa + "after watershed"+ File.separator;
dirraw= dirsa + "rawimages"+ File.separator;
//dir=getDirectory("choose results");
dir= dirsa + "indivisual_nuclei_ch1"+ File.separator;
dirc2= dirsa + "indivisual_nuclei_ch2"+ File.separator;
dirc3= dirsa + "indivisual_nuclei_ch3"+ File.separator;
dirc4= dirsa + "indivisual_nuclei_ch4"+ File.separator;

dirb= dirsa + "3d objects"+ File.separator;
dirx= dirsa + "binary mask"+ File.separator;
dir1= dirsa + "data"+ File.separator;
File.makeDirectory(dir); 
File.makeDirectory(dirc2); 
File.makeDirectory(dirc3); 
File.makeDirectory(dirc4); 
File.makeDirectory(dirb); 
File.makeDirectory(dirx); 
File.makeDirectory(dir1); 
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");


setBatchMode(true);
filenames=getFileList(dirw);
filenamesraw=getFileList(dirraw);
baseName=newArray(filenames.length);

for(f=0;f<filenames.length;f++){

	path=dirraw+filenamesraw[f];
	open(dirw+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 
	
	getVoxelSize(width, height, depth, unit);
	a=2000/(width*height*depth);// maximum volume is 1500 cu.microns
	a_1=200/(width*height*depth);// minimum volume is 200 cu.microns
	
	run("3D Objects Counter", "threshold=128 slice=12 min.=a_1 max.=a objects"); //Size filter
	saveAs("Tiff", dirb+baseName[f]+".tiff"); 
	obj_image=getTitle();
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	run("Close All");

	
	zlow=newArray(nb_obj);
	zhigh=newArray(nb_obj);
	wid=newArray(nb_obj);
	heig=newArray(nb_obj);
	xv=newArray(nb_obj);
	yv=newArray(nb_obj);
	name=newArray(nb_obj);
	print("\\Clear"); 
	
	for(p=0;p<nb_obj;p++){
		t=p+1;
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		zlows=z0+1;
		zhighs=z1+1;
		wids=x1-x0;
		heigs=y1-y0;
		
		open(dirb+obj_image);
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate range=zlows-zhighs");
		run("3D Convex Hull");//convex hull
		setThreshold(t, t);		
		run("Make Binary", "method=Default background=Default");
		names=baseName[f]+"_nucleus_mask_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirx+names);
		img1=getTitle();
		Ext.Manager3D_CloseResult("M");
		
		//Channel1
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=1 c_end=1 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dir+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");

		if(nchannels>1){
		//Channel2
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=2 c_end=2 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc2+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");
		}

		if(nchannels>2){
		//Channel3
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=3 c_end=3 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc3+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");
		}
		
		if(nchannels>3){
		//Channel4
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=4 c_end=4 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc4+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		
		selectImage(img1); 
		close("\\Others");
		}
		
		
		
		xv[p]=x0;
		yv[p]=y0;
		wid[p]=wids;
		heig[p]=heigs;
		name[p]=names;
		zlow[p]=z0+1;
		zhigh[p]=z1+1;
		run("Select None");
     	run("Close All");


     	
	}


	saveAs("Results",dir1+baseName[f]+"_objects_statistics.csv" ); 
	run("Clear Results");
	for(ol=0; ol<nb_obj; ol++){
		setResult("Name",ol,name[ol]);
		setResult("BX0",ol,xv[ol]);
		setResult("BY0",ol,yv[ol]);
		setResult("Width",ol,wid[ol]);
		setResult("Height",ol,heig[ol]);
		setResult("zstart",ol,zlow[ol]);
		setResult("zend",ol,zhigh[ol]);
	}
	saveAs("Results",dir1+baseName[f]+"_log.csv" ); 
	run("Close All");
    Ext.Manager3D_Close();
   	run("Clear Results");
}



