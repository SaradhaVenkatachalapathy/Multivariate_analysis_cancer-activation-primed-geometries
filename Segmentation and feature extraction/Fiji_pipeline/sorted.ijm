////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO SORT CELLS THAT HAVE BEEN SEGMENTED       											                    /////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   /////////////
///////  ASSUMPTIONS: The nuclei have been segmented using the 3d crop 
///////  DESCRIPTION: The user can input the sample type and then based on the nuclei sort the nuclei into NIH3T3, MCF7 or improper folders
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//dirsa=getDirectory("Please choose the directory from the 3d crop program");
dirsa = getArgument()
setBatchMode(false);
filenamesa=getFileList(dirsa);
dir= dirsa + "indivisual_nuclei_ch1"+ File.separator;
dir_sorted=dirsa + "data"+ File.separator;

		
		//open individial nuclei and use user input to sort them and the corresponding images of the other channels
		filenames=getFileList(dir);
		for(f=0;f<filenames.length;f++){
			//make the dialog box
			Dialog.create("Cell type");
  			items = newArray("NIH3T3", "Improper");
  			Dialog.addRadioButtonGroup("Cell type", items, 2, 1, "NIH3T3");
  		
			path1=dir+filenames[f];
			open(path1);
			t=getTitle();
			run("Make Montage...", "columns=6 rows=4 scale=1");
			//run("Brightness/Contrast...");
			run("Enhance Contrast", "saturated=0.35");

			Dialog.show;
  			a=Dialog.getRadioButton;
  			
			print(f+","+t+","+a);
			run("Close All");	
		}
		
selectWindow("Log");  //select Log-window 
saveAs("Text", dir_sorted+"Log.txt"); 	
print("\\Clear");