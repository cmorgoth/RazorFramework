
#!/usr/bin/env python

import subprocess, time, sys, os, glob

if __name__ == "__main__":
	directory = sys.argv[1]
	original_file_names = glob.glob(directory+"/*.pdf")
	for index in range(len(original_file_names)):		
		new_name =  original_file_names[index].replace("_0.","_0")
		os.system("mv "+original_file_names[index]+" "+new_name)

		
