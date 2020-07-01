import os
from distutils.dir_util import copy_tree
from shutil import copyfile


def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]
#import tail


#alpha = [0,  -10**(-4), -5*10**(-4), -10**(-3)];
#phi = [0, 0.0001, 0.0005,0.001,0.00005,0.0001]
# Re = [10**(-5), 10**(-4), 10**(-3)]
#phi = [0.00005, 0.0001]
#Main folder (do not use because of error in absolute path)
casesFolder = 'newCases_highPhi'

os.chdir(casesFolder)

removeDatFolder = False

if removeDatFolder:
        os.system('rm -r dat')
        
os.system('mkdir dat_phi')
os.system('mkdir vels_phi')
for geoId in range(4):

	for phi_id in range(4):

#		phi_id += 4
#		for re in range(len(Re)):

		# Pattern for creating folders
		subFolder = 'geo' + str(geoId) + 'phi' + str(phi_id) + 'nu0' 
		#os.system('mkdir ' + subFolder)    
		os.chdir(subFolder);
		
		f = open("./postProcessing/averageSample.dat" ,"r")
		#f = tail(f,1000)
		newfile = open("../dat_phi/" + subFolder +"_up.dat","w")
		newList = tail(f,100)
		for item in newList:
			newfile.write("%s\n"%item )

		#fromFile = "./postProcessing/averageSample.dat"
		#toFile = "../dat/" + subFolder +"_up.dat"
		#copy_tree(fromFile, toFile)
		#copyfile(fromFile, toFile)

                fromFile = "./permeationVel.csv"
                toFile = "../vels_phi/" + subFolder +"_vw.csv"
                #copy_tree(fromFile, toFile)
                copyfile(fromFile, toFile)


                timeList = []
                for time in next(os.walk('./postProcessing/sample'))[1]:
#                        if time not in ['system','constant','dynamicCode','postProcessing','VTK' ]:
                	timeList.append(int(time))
                #timeList = [int(time) for time in next(os.walk('./'))[1]]                        
                latest_time = str(max(timeList))
                        
                fromFile = "./postProcessing/sample/" + latest_time + "/data_U.xy"
		toFile = "../vels_phi/U_alongZ" + subFolder + ".xy"
		copyfile(fromFile, toFile)

		os.chdir('../')

    

