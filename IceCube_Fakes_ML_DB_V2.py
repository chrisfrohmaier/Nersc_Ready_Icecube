#2014-01-28 !!!!LATEST VERSION!!!!

#Before you run this make sure you have the input images (both science and mask) in this directory along with the following files:
#default.conv
#default.nnw
#Pipe_sexfile_CMD.sex
#PTF_Transform_Param.param

#!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!--- YOU NEED THE MODULES BELOW

import numpy, os, random, glob, shutil, time, subprocess, math
from multiprocessing import Pool
from astropy.io import fits
from astropy import wcs
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True) #Ignores UserWarnings otherwise Astropy spits out loads when it overwrites files
warnings.filterwarnings('ignore', category=Warning, append=True)
from astropy.io.fits import getheader 


def file_structure(): #This definition creates the file structure for results to be added into

	if not os.path.exists('Results_V2'):
		os.makedirs('Results_V2')
	if not os.path.exists('Results_V2/Catalog'):
		os.makedirs('Results_V2/Catalog')
	if not os.path.exists('Results_V2/Fake_Star_Catalog'):
		os.makedirs('Results_V2/Fake_Star_Catalog')
	if not os.path.exists('Results_V2/Fakes_added'):
		os.makedirs('Results_V2/Fakes_added')	
	if not os.path.exists('Results_V2/Galaxies'):
		os.makedirs('Results_V2/Galaxies')

def Sextract(science_image,zeropoint,seeing,saturation,gain): #Runs Sextractor and creates a catalog of all the stars
	
	subprocess.call('sex -c Pipe_sexfile_CMD.sex '+science_image[0]+science_image[1]+'.fits -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv -CATALOG_NAME Results_V2/Catalog/'+science_image[1]+'_Catalog_V2.cat -WEIGHT_IMAGE '+science_image[0]+science_image[1]+'.weight.fits -MAG_ZEROPOINT'+'	'+str(zeropoint)+' -SEEING_FWHM '+str(seeing)+' -SATUR_LEVEL '+str(saturation)+' -GAIN '+str(gain)+' -VERBOSE_TYPE QUIET',shell=True)

def Enough_Objects(science_image): #Checks that sextractor has found at least 200 objects
	enough=True
	test=os.popen('wc -l Results_V2/Catalog/'+science_image[1]+'_Catalog_V2.cat').read()
	rows=test.split()
	if float(rows[0])<217: #217 because the catalog have commented lines which are counted
			return False
	
def Selecting_Bright(science_image): #Selected the top 20 brightest stars in the catalog
	f=open('Results_V2/Catalog/'+science_image[1]+'_Catalog_V2.cat') #CHANGE THIS SO THAT IT CAN TAKE ANY CATALOG INPUT
	fin=f.readline()
	
	xcord=[]
	ycord=[]
	flux_array=[]
	mag_array=[]
	background_array=[]
	mag_best_array=[]
	Alpha_Sky_array=[]
	Delta_Sky_array=[]

	while fin: #DEFINING BINS BASED ON MAGNITUDE OF STAR
		
		if fin.startswith('#'):
			#h.write(fin)
			#print fin
			fin=f.readline()
			continue
		ln=fin.split()
		#print ln[15]
		mag=float(ln[2]) #magnitude of star	
		x=float(ln[3])
		y=float(ln[4])
		flux=float(ln[1])	
		background=float(ln[5])
		mag_best=float(ln[17])
		alpha=float(ln[18])
		delta=float(ln[19])

		
		if float(ln[7])<0.3: #Not Elliptical
			if float(ln[9])>0.5: #Considered a good star
				if x>100.0 and x<1948.0 and y>100.0 and y<3996.0: #No Edge Stars
					if int(ln[8])==0:
						if mag>12:
							
							xcord.append(x); ycord.append(y); mag_array.append(mag); flux_array.append(flux); background_array.append(background); mag_best_array.append(mag_best); Alpha_Sky_array.append(alpha); Delta_Sky_array.append(delta)
		fin=f.readline()					
	f.close()
	mag_array, xcord, ycord,flux_array, background_array, mag_best_array, Alpha_Sky_array, Delta_Sky_array= (list(x) for x in zip(*sorted(zip(mag_array ,xcord, ycord, flux_array, background_array, mag_best_array, Alpha_Sky_array, Delta_Sky_array))))
	mag_array=mag_array[:20]
	xcord=xcord[:20]
	ycord=ycord[:20]
	flux_array=flux_array[:20]
	background_array=background_array[:20]
	mag_best_array=mag_best_array[:20]
	Alpha_Sky_array=Alpha_Sky_array[:20]
	Delta_Sky_array=Delta_Sky_array[:20]
		
	return xcord, ycord, mag_array, flux_array, background_array, mag_best_array, Alpha_Sky_array, Delta_Sky_array

def selecting_galaxies(science_image,): #Finds and creates a catalog of Galxies
	science_image=science_image
	f=open('Results_V2/Catalog/'+science_image[1]+'_Catalog_V2.cat')
	g=open('Results_V2/Galaxies/'+science_image[1]+'_Galaxy_Catalog_V2.cat','w')
	l=open('Results_V2/Galaxies/'+science_image[1]+'_Galaxy_regions_V2.reg','w')
	m=open('Results_V2/Galaxies/'+science_image[1]+'_Star_regions_V2.reg','w')
	fin=f.readline()
	X2_array=[]
	Y2_array=[]
	FWHM_array=[]
	Ellipticity_array=[]
	class_s_array=[]
	counts=0
	while fin:
		if fin.startswith('#'):
			#h.write(fin)
			#print fin
			fin=f.readline()
			continue
		ln=fin.split()
		class_s=float(ln[9])
		gal_mag=float(ln[2])
		xcord=float(ln[3])
		ycord=float(ln[4])
		X2=float(ln[13])
		Y2=float(ln[14])
		X2_array.append(X2)
		Y2_array.append(Y2)
		
		FWHM=float(ln[16])
		FWHM_array.append(FWHM)
		elips=float(ln[7])
		Ellipticity_array.append(elips)
		class_s_array.append(class_s)
		if class_s<0.5 and gal_mag>14:
			if FWHM<15:
				if xcord>50.0 and xcord<1998.0 and ycord>50.0 and ycord<4046.0: #No Edge Galaxies
					g.write(fin)
					#g.write(str((ln[0]))+' '+str((ln[1]))+' '+str((ln[2]))+' '+str((ln[3]))+' '+str((ln[4]))+' '+str(ln[5])+' '+str((ln[6]))+' '+str((ln[7]))+' '+str((ln[8]))+' '+str((ln[9]))+' '+str((ln[10]))+' '+str((ln[11]))+' '+str(ln[12])+' '+str((ln[13]))+' '+str((ln[14]))+' '+str((ln[15]))+' '+str((ln[16]))+'\n')
					l.write(str(xcord)+' '+str(ycord)+'\n')
					counts+=1
		else:
			m.write(str(xcord)+' '+str(ycord)+'\n')		

				
		fin=f.readline()


def Scaling(science_image ,xcord, ycord, mag_array, flux_array, background_array, zpt, fake_stars, CCD_Num, magnitude_best,alpha_sky, delta_sky):

	ranmagarray=[]
	xcord_star=[]
	ycord_star=[]
	newx_star=[]
	newy_star=[]
	mag_array_star=[]
	flux_array_star=[]
	ran_mag_star=[]
	ran_flux_star=[]
	background_array_star=[]
	scaling_factor_star=[]
	CCD_Num_star=[]
	faint_fake=max(mag_array)+1.0
	best_mag_array=[]
	alpha_array=[]
	delta_array=[]
	#print 'faint_fake', faint_fake
	for i in range(0,fake_stars):
		ran_mag=random.uniform(faint_fake, 23.0) #The fake stars will be in this range of magnitudes
		ran_flux=10.0**((ran_mag-zpt)/(-2.5))
		ranmagarray.append(ran_mag)
		star=int(random.uniform(0,len(xcord)-1))
		
		scaling_factor=((ran_flux)/flux_array[star])
		
		newX=random.uniform(100.0,1948.0) #This lines don't actually do anything anymore! The new x and y co-ordinates are based on galaxy locations.
		newY=random.uniform(100.0, 3996.0)

		xcord_star.append(xcord[star]); ycord_star.append(ycord[star]); newx_star.append(newX); newy_star.append(newY); mag_array_star.append(mag_array[star]); flux_array_star.append(flux_array[star]); ran_mag_star.append(ran_mag); ran_flux_star.append(ran_flux); background_array_star.append(background_array[star]); scaling_factor_star.append(scaling_factor); CCD_Num_star.append(CCD_Num); best_mag_array.append(magnitude_best[star]); alpha_array.append(alpha_sky[star]); delta_array.append(delta_sky[star])
		i+=1
	return xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, faint_fake, best_mag_array, alpha_array, delta_array


def add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star,  mag_best_star, alpha_array, delta_array, zeropoint, seeing, saturation, gain, readnoise, MOONILLF, AIRMASS, ELLIP, MEDSKY, SKYSIG,LMT_MG, MJD,  MoonRA, MoonDec, PTFFIELD):
	#This step finds the galaxies and adds fake stars to them
	h=open('Results_V2/Galaxies/'+science_image[1]+'_Galaxy_Catalog_V2.cat') #Opens the Galaxy catalog
	f=open('Results_V2/Fake_Star_Catalog/'+science_image[1]+'_Fake_Star_Catalog_V2.dat','w') #Opens the fake star catalog
	reg=open('Results_V2/Fake_Star_Catalog/'+science_image[1]+'_Fakes_Star_Regions_V2.reg','w') #creates region file
	regworld=open('Results_V2/Fake_Star_Catalog/'+science_image[1]+'_Fakes_Star_Regions_V2_World.reg','w')
	Galworld=open('Results_V2/Fake_Star_Catalog/'+science_image[1]+'_Galaxy_Regions_V2_Sex.reg','w')
	hin=h.readline() #reads first line
	
	fake_star_array=[] #prepares and array for fake stars
	
	for stars in range(0,len(xcord_star)):
		fake_star_array.append(stars) #creates an array [0,1,2,3, etc]

	gal_line_array=[] #prepares and array for galaxies
	while hin: #Adds all galxies into an array
		gal_line_array.append(hin)
		hin=h.readline()
	#while fin: #adds all stars into an array
	#	fake_star_array.append(fin)
	#	fin=f.readline()
	h.close()
	#f.close()

	hdulist_sci= fits.open(science_image[0]+science_image[1]+'.fits',ignore_missing_end=True) #THE SCIENCE DATA THAT WILL OPENED AND MODIFIED
	science_data= hdulist_sci[0].data
	

	resy=science_data.shape[0]
	resx=science_data.shape[1]




	#print len(fake_star_array), ' Fake Stars have been added to ', len(fake_star_array), ' Galaxies'
	#print gal_line_array[1]
	#print 'Number of Fakes to Be added: ', len(xcord_star)
	num_of_gal_fakes=0
	j=open('Results_V2/Fakes_added/'+science_image[1]+'_Flux_Boxes_V2.dat','w')
	galaxy_mask=numpy.ones((resy,resx),dtype=bool)
	for i in range(0,len(xcord_star)): #Will only add n fake stars to n Galaxies
		#host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice
		
		source_star=fake_star_array.pop(random.randrange(0,len(fake_star_array))) #selecting a random source star. Used .pop() so that the same star isnt chosen twice
			
		
		#print y 
		#print 'len: ',len(gal_line_array)
		#ln=host_galaxy.split()
		#x=float(ln[3])
		#y=float(ln[4])
		

		
		#print 'Lenth of Possible Galaxies: ', len(gal_line_array)
		while len(gal_line_array)>0: #and num_of_gal_fakes<len(xcord_star):
			#host_galaxy=random.choice(gal_line_array)
			#print 'Host Galaxy: ', host_galaxy
			host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice
			ln=host_galaxy.split()
			x=float(ln[3]) #Galaxy X-cord
			y=float(ln[4]) #Galaxy Y-cord
			galaxy_mag_auto=float(ln[2])
			galaxy_mag_best=float(ln[17])
			galaxy_flux=float(ln[1])
			galaxy_background=float(ln[5])

			
			galaxy_mask[0:50,0:2048]=False
			galaxy_mask[4056:4096,0:2048]=False
			galaxy_mask[0:4096,0:50]=False
			galaxy_mask[0:4096,1998:2048]=False

			if galaxy_mask[y,x]==False:
				#print 'Cant Go there'
				continue
			else:

				r=40 #Radius of Mask
				ym,xm = numpy.ogrid[-y:resy-y, -x:resx-x] #Some clever numpy stuff
				mask = xm*xm + ym*ym <= r*r
				galaxy_mask[mask]=False
				
				'Doing Usual Business'
				#
				CXX=float(ln[10])
				CYY=float(ln[11])
				CXY=float(ln[12])
				R=3.0
				#----------------THIS IS A GRID USED FOR CALCULATING STAR POSITIONS !!NOT!! FOR SCALING STARS
				#Draw a large grid around the galaxy of say 20,20 pixels. Run through that grid for every x and y and if it satisfies the equation on page 32 of the sextractor manual then append it to an
				#array. Then randomly choose a coordinate and insert a fake star there.
				grid_size=20
				grid_startx=int(x-grid_size/2.0)
				grid_starty=int(y-grid_size/2.0)		
				grid_finx=int(x+grid_size/2.0)
				grid_finy=int(y+grid_size/2.0)
				#----------------
				#good_x=[]
				#good_y=[]
				good_coords=[]
				#This loops runs through the grid and finds the pixel co-ordinates that satisfy the inequality below
				for q in range(grid_starty,grid_finy):
					for w in range(grid_startx,grid_finx):
						ell_p=(CXX*((w-x)*(w-x)))+(CYY*((q-y)*(q-y)))+(CXY*((w-x)*(q-y))) #IF THIS INEQUALITY IS SATISFIED THEN A STAR CAN GO AT THIS LOCATION 
						if ell_p<=3.0:
							good_coords.append([w,q])
							#good_y.append(q)
							#good_x.append(w)
				#print i, good_coords
				fake_star_positon=random.choice(good_coords) #Chooes a random co-ordinate
				#kn=source_star.split()
				sourcex=xcord_star[source_star] #stars current x location
				sourcey=ycord_star[source_star] #stars current y location
				newx=fake_star_positon[0] #where the star will go x
				newy=fake_star_positon[1] #where the star will go y

				#Saving Flux info Pre Fakes Addition

				pixel_fluxes=[]
				
				#print science_image
				#print 'newx newy sci:	', newx, newy, (science_data[newy,newx])
				pixel_fluxes.append((science_data[newy,newx]))
				for s in boxsize:
					start_box_x=newx-int(s/2)
					fin_box_x=newx+int(s/2)+1
					start_box_y=newy-int(s/2)
					fin_box_y=newy+int(s/2)+1
					#print 'Box X Dimensions: ', fin_box_x - start_box_x
					#print 'Box Y Dimensions: ', fin_box_y - start_box_y
					Box_Matrix=[[0 for l in xrange(s)] for l in xrange(s)]
					#print Box_Matrix
					countg=0
					
					for g in range (start_box_y,fin_box_y):
						counth=0
						for h in range(start_box_x,fin_box_x):
							Box_Matrix[countg][counth]=(science_data[g,h])
							#print Box_Matrix[countg][counth]
							counth+=1
							#print Box_Matrix
						countg+=1
					pixel_fluxes.append(Box_Matrix)
							
					#print 's newx newy sci:	', s, newx, newy, (science_data[newy,newx])
				j.write(str(pixel_fluxes)+'\n')
				
				reg.write(str(newx)+' '+str(newy)+'\n') #fake star region file
				
				scale_fac=scaling_factor_star[source_star] #scale factor
				back=background_array_star[source_star] #background
				
				#---Old area to be scaled---
				startx=int(sourcex-10.0)
				starty=int(sourcey-10.0)		
				finx=int(sourcex+10.0)
				finy=int(sourcey+10.0)
				
				#---New area to have flux added---
				Nstartx=newx-10.0
				Nstarty=newy-10.0
				Nfinx=newx+10.0
				Nfiny=newy+10.0


				newdata=numpy.ones((20,20)) #Preparing a blank gird for scaled objects

				newdata[0:20,0:20]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object

				science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata #Modifying the science image
				

				f.write(str(xcord_star[source_star])+' '+str(ycord_star[source_star])+' '+str(alpha_array[source_star])+' '+str(delta_array[source_star])+' '+str(newx)+' '+str(newy)+' '+str(mag_array_star[source_star])+' '+str(mag_best_star[source_star])+' '+str(flux_array_star[source_star])+' '+str(ran_mag_star[source_star])+' '+str(ran_flux_star[source_star])+' '+str(background_array_star[source_star])+' '+str(scaling_factor_star[source_star])+' '+str(int(PTFFIELD))+' '+str(CCD_Num_star[source_star])+' '+str(x)+' '+str(y)+' '+str(galaxy_mag_auto)+' '+str(galaxy_mag_best)+' '+str(galaxy_flux)+' '+str(galaxy_background)+' '+str(gain)+' '+str(readnoise)+' '+str(MOONILLF)+' '+str(MoonRA)+' '+str(MoonDec)+' '+str(AIRMASS)+' '+str(seeing)+' '+str(ELLIP)+' '+str(MEDSKY)+' '+str(SKYSIG)+' '+str(zeropoint)+' '+str(LMT_MG)+' '+str(MJD)+'\n')
				
				num_of_gal_fakes+=1
				break
	
	hdulist_sci.writeto(science_image[0]+science_image[1]+'_fakesV2.fits', output_verify='ignore', clobber=True) #Saving image after loop of 200 Stars is complete
	j.close()
	reg.close()	
	f.close()

	#print num_of_gal_fakes, 'fake Stars Added to Galaxies in the Image: ', science_image[1]

	#Creating a Galaxy Mask Fits file
	
	
	galaxy_mask_float=galaxy_mask.astype(int)
	hdu=fits.PrimaryHDU(galaxy_mask_float)
	hdu.scale(type='int16')
	hdulist=fits.HDUList([hdu])
	
	#print hdulist.info()
	hdulist.writeto('Results_V2/Fake_Star_Catalog/'+science_image[1]+'_GMask_V2.fits', clobber=True, output_verify='ignore')

def Sub_ML_DB(science_image):
	#print 'Sub Stuff for:', science_image[1]
	ref=science_image[1]
	new_image=science_image[1]+'_fakesV2.fits'

	subprocess.call('cd '+str(science_image[0])+'; /project/projectdirs/deepsky/rates/icecube/scripts/diffem $'+str(ref)+'.fits $'+str(new_image)+'', shell=True)

def Execute(run):
	#print '!!!!!!', run
	science_image=run

	
	#print '@@@@@', run[0], run[1]
	#print '!!!!!',science_image[0], science_image[1]

	sci_fil=science_image[0]+science_image[1]+'.fits'
	#maskfile=science_image[0]+science_image[1]+'.weight'
	#print  '######', science_image[1],'.weight'
	#print sci_fil
	#print maskfile

	
	#print 'Name: ', science_image

	try:
		hdulist_multi_sci=fits.open(science_image[0]+science_image[1]+'.fits')
		#print '++++ multi_mask assign ', science_image
		
	except IOError or Warning or UnboundLocalError:
		bad_images=open('Results_V2/Bad_Images_V2.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Astropy Could not Open the .fits file')+'\n')
		bad_images.close()
		print 'Cant open Science'
		
		

		return
	
		
	zeropoint=float(hdulist_multi_sci[0].header['UB1_ZP'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=55000.0 #float(hdulist_multi_sci[0].header['SATURATE'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])
	CCD_Num=float(hdulist_multi_sci[0].header['CCDID'])
	PTFFIELD=int(hdulist_multi_sci[0].header['PTFFIELD'])
	readnoise=float(hdulist_multi_sci[0].header['READNOI'])
	MOONILLF=float(hdulist_multi_sci[0].header['MOONILLF'])
	AIRMASS=float(hdulist_multi_sci[0].header['AIRMASS'])
	ELLIP=float(hdulist_multi_sci[0].header['ELLIP'])
	MEDSKY=float(hdulist_multi_sci[0].header['MEDSKY'])
	SKYSIG=float(hdulist_multi_sci[0].header['SKYSIG'])
	LMT_MG=float(hdulist_multi_sci[0].header['LMT_MG'])
	MJD=float(hdulist_multi_sci[0].header['OBSMJD'])
	MoonRA=float(hdulist_multi_sci[0].header['MOONRA'])
	MoonDec=float(hdulist_multi_sci[0].header['MOONDEC'])




	fake_stars= 50 #number of fake stars per image (integer please!)

	hdulist_multi_sci.close()


	Sextract(science_image,zeropoint,seeing,saturation,gain)
	
	catsize=Enough_Objects(science_image)
	if catsize==False:
			print science_image, 'didn\'t have enough objects detected so it was moved to Results_V2/Bad_Images/ and the newly created weight map, sex file and catalog have been deleted'
			bad_images=open('Results_V2/Bad_Images_V2.dat','a')
			bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Sextractor did not detect enough objects (<200)')+'\n')
			os.remove('Results_V2/Catalog/'+science_image[1]+'_Catalog_V2.cat')
			return

	x, y, mag, flux, back, magnitude_best, alpha_sky, delta_sky = Selecting_Bright(science_image)
	#print 'Selecting Bright Done'
	xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, faint_fake, mag_best_star, alpha_array, delta_array =Scaling(science_image, x, y, mag, flux, back, zeropoint, fake_stars, CCD_Num, magnitude_best, alpha_sky, delta_sky)
	#print 'Scaling Done'
	mag_log=open('Results_V2/Magnitude_Log_File.dat','a')
	mag_log.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str(mag[0])+' '+str(mag[-1])+' '+str(faint_fake)+' '+str('23.0')+'\n')
	
	selecting_galaxies(science_image)
	#print 'Selected Galaxies'
	
	boxsize=[3,5,7]
	add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, mag_best_star, alpha_array, delta_array, zeropoint, seeing, saturation, gain, readnoise, MOONILLF, AIRMASS, ELLIP, MEDSKY, SKYSIG, LMT_MG, MJD, MoonRA, MoonDec, PTFFIELD)
	good_images=open('Results_V2/Good_Images_V2.dat','a')
	good_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+'\n')

	Sub_ML_DB(science_image)
#-----------------------------------RUN PIPELINE------------------------------------------

def Run_All_V2():
	file_structure()
	all_fits=[] #Establishing an array to find the files
	#path=[]
	#fnames=[]
	for dirpath,dirname,filenames in os.walk(os.path.abspath('../../fakes')): #Traverses through a directory tree

		for file in filenames:
			fileex=os.path.splitext(file)[-1] #Splits the file name, [-1] means it will look at the extension
			if fileex== '.fits': #wanted all .fits files 

				all_fits.append([dirpath, file])
	#print all_fits
	test_files=0
	testlim=200


	science_fits=[]
	for i in range(len(all_fits)):
		#fname=all_fits[1]
		ln=all_fits[i]
		fname=ln[1].split('.')
		#print fname
		
		
		if fname[-2]=='w'and test_files<testlim:
		
			science_fits.append([ln[0]+str('/'), (os.path.splitext(ln[1])[0])])
			test_files+=1

		elif test_files>testlim:
			break

	bad_images=open('Results_V2/Bad_Images_V2.dat','w')
	bad_images.close()
	good_images=open('Results_V2/Good_Images_V2.dat','w')
	good_images.close()
	mag_log=open('Results_V2/Magnitude_Log_File.dat','w')
	mag_log.close()
	mag_log_col=open('Results_V2/Magnitude_Log_File.columns','w')
	mag_log_col.write(str('1. Path to Image')+'\n'+str('2. Brightest Source')+'\n'+str('3. Faintest Source Mag')+'\n'+str('4. Brightest Fake')+'\n'+str('5. Faintest Fake'))
	mag_log_col.close()

	
	t0=time.time()
	processors=8
	pool=Pool(processors)
	pool.map(Execute,science_fits)
	pool.close()

	print 'V2 took: ', time.time()-t0, 'seconds'


