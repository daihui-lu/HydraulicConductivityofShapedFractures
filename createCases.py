import os
from distutils.dir_util import copy_tree
import subprocess

alpha = [0,  -10**(-4), -5*10**(-4), -10**(-3)];
phi = [0, 0.0005,0.001,0.005]
#Re = [1, 10, 100]
nu = [1]
removeOld = True;

if removeOld:
	os.system('rm -r newCases_highPhi')


 #Main folder (do not use because of error in absolute path)
casesFolder = 'newCases_highPhi'

os.system('mkdir ' + casesFolder)
os.chdir(casesFolder)

for geoId in range(4):

	for phi_id in range(4):

		for nu_id in range(1):
#		phi_id += 4

#		for re in range(len(Re)):

			# Pattern for creating folders
			subFolder = 'geo' + str(geoId) + 'phi' + str(phi_id) +'nu' + str(nu_id) 
			os.system('mkdir ' + subFolder)    
			#os.system('cd ' + subFolder);
                	os.chdir(subFolder)
			# Case creation			
			fromDirectory = "../../baseCase"
			toDirectory = "./"
			copy_tree(fromDirectory, toDirectory)

			file = open("parameters.H","w") 
			 
			file.write("alpha   " + str(alpha[geoId])+  "; phi    " + str(phi[phi_id]) + "; nu	"+ str(nu[nu_id])+";")

			file.close() 
	
			os.system('./Allclean &')
			os.system('./Allrun &')
           

		
			os.chdir('../')

    

