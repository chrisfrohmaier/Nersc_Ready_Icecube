#2014_07_09

#Before you run this make sure you have the input images (both science and mask) in this directory along with the following files:
#default.conv
#default.nnw
#Pipe_sexfile_CMD.sex
#PTF_Transform_Param.param

#!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!--- YOU NEED THE MODULES BELOW
#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy, os, random, glob, shutil, time, subprocess, math
import multiprocessing
from multiprocessing import Pool
from astropy.io import fits
from astropy import wcs
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True) #Ignores UserWarnings otherwise Astropy spits out loads when it overwrites files
warnings.filterwarnings('ignore', category=Warning, append=True)
from astropy.io.fits import getheader
import sys

global vnum
vnum=int(sys.argv[1])

print 'THIS VNUM IS: ', vnum
#!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!



def file_structure(): #This definition creates the file structure for results to be added into

	if not os.path.exists('Results_V'+str(vnum)+''):
		os.makedirs('Results_V'+str(vnum)+'')
	if not os.path.exists('Results_V'+str(vnum)+'/Catalog'):
		os.makedirs('Results_V'+str(vnum)+'/Catalog')
	if not os.path.exists('Results_V'+str(vnum)+'/Fake_Star_Catalog'):
		os.makedirs('Results_V'+str(vnum)+'/Fake_Star_Catalog')
	if not os.path.exists('Results_V'+str(vnum)+'/Fakes_added'):
		os.makedirs('Results_V'+str(vnum)+'/Fakes_added')
	if not os.path.exists('Results_V'+str(vnum)+'/Galaxies'):
		os.makedirs('Results_V'+str(vnum)+'/Galaxies')

def Sextract(science_image,zeropoint,seeing,saturation,gain): #Runs Sextractor and creates a catalog of all the stars

	subprocess.call('sex -c Pipe_sexfile_CMD.sex '+science_image[0]+science_image[1]+'.fits -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv -CATALOG_NAME Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat -WEIGHT_IMAGE '+science_image[0]+science_image[1]+'.weight.fits -MAG_ZEROPOINT'+'	'+str(zeropoint)+' -SEEING_FWHM '+str(seeing)+' -SATUR_LEVEL '+str(saturation)+' -GAIN '+str(gain)+' -PHOT_FLUXFRAC 0.2,0.5,0.9 -VERBOSE_TYPE QUIET',shell=True)

def Enough_Objects(science_image): #Checks that sextractor has found at least 200 objects
	enough=True
	test=os.popen('wc -l Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat').read()
	rows=test.split()
	if float(rows[0])<217: #217 because the catalog have commented lines which are counted
			return False

def Selecting_Bright(science_image): #Selected the top 20 brightest stars in the catalog
	f=open('Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat') #CHANGE THIS SO THAT IT CAN TAKE ANY CATALOG INPUT
	fin=f.readline()

	xcord=[]
	ycord=[]
	flux_array=[]
	mag_array=[]
	background_array=[]
	mag_best_array=[]
	Alpha_Sky_array=[]
	Delta_Sky_array=[]

	while fin:

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

	#print 'Doing Galaxy Stuff for: ', science_image
	f=open('Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat')
	g=open('Results_V'+str(vnum)+'/Galaxies/'+science_image[1]+'_Galaxy_Catalog_V'+str(vnum)+'.cat','w')
	l=open('Results_V'+str(vnum)+'/Galaxies/'+science_image[1]+'_Galaxy_regions_V'+str(vnum)+'.reg','w')

	fin=f.readline()


	counts=0

	#Creating the Numpy Grid of Galaxies
	hdulist_sci= fits.open(science_image[0]+science_image[1]+'.fits',ignore_missing_end=True) #THE SCIENCE DATA THAT WILL OPENED AND MODIFIED
	science_data= hdulist_sci[0].data
	header_data= hdulist_sci[0].header
	resy=science_data.shape[0]
	resx=science_data.shape[1]
	galaxy_area=numpy.ones((resy,resx),dtype=bool)
	#galtest=numpy.ones((resy,resx),dtype=bool)

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
		CXX=float(ln[10])
		CYY=float(ln[11])
		CXY=float(ln[12])


		FWHM=float(ln[16])


		if class_s<0.5 and gal_mag>14 and gal_mag<30:
			#if FWHM<15:
			if xcord>21.0 and xcord<2027.0 and ycord>21.0 and ycord<4075.0: #No Edge Galaxies
				g.write(fin)
				#g.write(str((ln[0]))+' '+str((ln[1]))+' '+str((ln[2]))+' '+str((ln[3]))+' '+str((ln[4]))+' '+str(ln[5])+' '+str((ln[6]))+' '+str((ln[7]))+' '+str((ln[8]))+' '+str((ln[9]))+' '+str((ln[10]))+' '+str((ln[11]))+' '+str(ln[12])+' '+str((ln[13]))+' '+str((ln[14]))+' '+str((ln[15]))+' '+str((ln[16]))+'\n')
				l.write(str(xcord)+' '+str(ycord)+'\n')
				counts+=1
				gyo,gxo= numpy.indices((40,40))

				gx=gxo+(xcord-20)
				gy=gyo+(ycord-20)

				a_galaxy=numpy.where(((CXX*((gx-xcord)*(gx-xcord)))+(CYY*((gy-ycord)*(gy-ycord)))+(CXY*((gx-xcord)*(gy-ycord))) <= 3))



				galaxy_area[a_galaxy[0]+int(ycord),a_galaxy[1]+int(xcord)]=False

				#galtest[a_galaxy[0]+int(ycord-36),a_galaxy[1]+int(xcord-36)]=0.0





		fin=f.readline()

	#galaxies_int=galtest.astype(float)
	#hdu_gals=fits.PrimaryHDU(data=galaxies_int,header=header_data)
	#hdu_gals.scale(type='int16')
	#hdulist_gals=fits.HDUList([hdu_gals])

	#print hdulist.info()
	#hdulist_gals.writeto('Results_V'+str(vnum)+'/Galaxies/'+science_image[1]+'_GALAXIES_V'+str(vnum)+'.fits', clobber=True, output_verify='ignore')
	#hdulist_gals.close()		
	#numpy.savetxt('Results_V'+str(vnum)+'/Galaxies/'+science_image[1]+'_Galaxy_Grids.dat', galaxy_area, delimiter=' ',fmt='%d')
	#print 'Finished Doing Galaxy Stuff: ', science_image
	f.close()
	g.close()
	l.close()
	return galaxy_area

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
		ran_mag=random.uniform(faint_fake, 22.5) #The fake stars will be in this range of magnitudes
		ran_flux=10.0**((ran_mag-zpt)/(-2.5))
		ranmagarray.append(ran_mag)
		star=int(random.uniform(0,len(xcord)-1))

		scaling_factor=((ran_flux)/flux_array[star])

		newX=random.uniform(100.0,1948.0) #This lines don't actually do anything anymore! The new x and y co-ordinates are based on galaxy locations and the hostless parameters later on.
		newY=random.uniform(100.0, 3996.0)

		xcord_star.append(xcord[star]); ycord_star.append(ycord[star]); newx_star.append(newX); newy_star.append(newY); mag_array_star.append(mag_array[star]); flux_array_star.append(flux_array[star]); ran_mag_star.append(ran_mag); ran_flux_star.append(ran_flux); background_array_star.append(background_array[star]); scaling_factor_star.append(scaling_factor); CCD_Num_star.append(CCD_Num); best_mag_array.append(magnitude_best[star]); alpha_array.append(alpha_sky[star]); delta_array.append(delta_sky[star])
		i+=1
	return xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, faint_fake, best_mag_array, alpha_array, delta_array


def add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star,  mag_best_star, alpha_array, delta_array, zeropoint, seeing, saturation, gain, readnoise, MOONILLF, AIRMASS, ELLIP, MEDSKY, SKYSIG,LMT_MG, MJD,  MoonRA, MoonDec, PTFFIELD, galareas):
	#This step finds the galaxies and adds fake stars to them
	h=open('Results_V'+str(vnum)+'/Galaxies/'+science_image[1]+'_Galaxy_Catalog_V'+str(vnum)+'.cat') #Opens the Galaxy catalog
	f=open('Results_V'+str(vnum)+'/Fake_Star_Catalog/'+science_image[1]+'_Fake_Star_Catalog_V'+str(vnum)+'.dat','w') #Opens the fake star catalog
	reg=open('Results_V'+str(vnum)+'/Fake_Star_Catalog/'+science_image[1]+'_Fakes_Star_Regions_V'+str(vnum)+'.reg','w') #creates region file
	ch=open('Results_V'+str(vnum)+'/Fake_Star_Catalog/'+science_image[1]+'_add_F2G_progress_V'+str(vnum)+'.dat','a') #Debugging file


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
	ch.write(str('Resolution:')+' '+str(resy)+' '+str(resx)+'\n')



	#print len(fake_star_array), ' Fake Stars have been added to ', len(fake_star_array), ' Galaxies'
	#print gal_line_array[1]
	#print 'Number of Fakes to Be added to hosts: ', int(len(xcord_star)*0.9)
	num_of_fakes_all=0
	#j=open('Results_V'+str(vnum)+'/Fakes_added/'+science_image[1]+'_Flux_Boxes_V'+str(vnum)+'.dat','w')
	galaxy_mask=numpy.ones((resy,resx),dtype=bool)
	for i in range(0,int(len(xcord_star)*0.9)): #Will only add n*0.9 fake stars to n Galaxies
		#host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice

		source_star=fake_star_array.pop(random.randrange(0,len(fake_star_array))) #selecting a random source star. Used .pop() so that the same star isnt chosen twice

		ch.write(str('!!!!!')+' '+str(num_of_fakes_all)+' '+str('!!!!!')+'\n')
		#print y
		#print 'len: ',len(gal_line_array)
		#ln=host_galaxy.split()
		#x=float(ln[3])
		#y=float(ln[4])



		#print 'Lenth of Possible Galaxies: ', len(gal_line_array)
		while len(gal_line_array)>0: #and num_of_fakes_all<len(xcord_star):
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


			galaxy_mask[0:40,0:2048]=False
			galaxy_mask[4056:4096,0:2048]=False
			galaxy_mask[0:4096,0:40]=False
			galaxy_mask[0:4096,2008:2048]=False
			ch.write(str('Galaxy Mask part 1 done')+'\n')
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
				FR_02=float(ln[20])
				FR_05=float(ln[21])
				FR_09=float(ln[22])
				Host_Elongation=float(ln[6])
				Host_Ellipticity=float(ln[7])
				Host_Alpha=float(ln[18])
				Host_Dec=float(ln[19])


				R=3.0
				#----------------THIS IS A GRID USED FOR CALCULATING STAR POSITIONS !!NOT!! FOR SCALING STARS
				#Draw a large grid around the galaxy of say 20,20 pixels. Run through that grid for every x and y and if it satisfies the equation on page 32 of the sextractor manual then append it to an
				#array. Then randomly choose a coordinate and insert a fake star there.


				gyo,gxo= numpy.indices((40,40))

				gx=gxo+(x-20)
				gy=gyo+(y-20)

				#print 'GX GY: ', gx, gy
				#print  'Done gx gy'
				goodcords=numpy.where(((CXX*((gx-x)*(gx-x)))+(CYY*((gy-y)*(gy-y)))+(CXY*((gx-x)*(gy-y))) <= 3)==True)
				ch.write(str('Good Cords section CXX CYY: Done')+'\n')
				if len(goodcords[0])-1<1:
					continue
				#print 'Done Good Cords'
				#print 'Length of goodcords xy: ', len(goodcords[0]), len(goodcords[1])
				#print 'Good Cords: ', goodcords
				gc=random.randint(0,len(goodcords[0])-1)
				ch.write(str('Choosing gc: Done')+'\n')
				#print 'Done gc'
				newy=(goodcords[0][gc])+(y-20)
				newx=(goodcords[1][gc])+(x-20)
				#print 'DOne newx newy'

				sourcex=xcord_star[source_star] #stars current x location
				sourcey=ycord_star[source_star] #stars current y location
				ch.write(str('Newy and Newx: Done')+'\n')

				##Creating the fboxes

				fbox1=numpy.sum(science_data[newy,newx])
				fbox2=numpy.sum(science_data[newy,newx]) + numpy.sum(science_data[newy-1.0,newx]) + numpy.sum(science_data[newy+1.0,newx]) + numpy.sum(science_data[newy,newx-1.0]) + numpy.sum(science_data[newy,newx+1.0])
				fbox3=numpy.sum(science_data[newy-1.0:newy+2.0, newx-1.0:newx+2.0])
				fbox4=numpy.sum(science_data[newy-1.0:newy+2.0, newx-1.0:newx+2.0]) + numpy.sum(science_data[newy-2.0,newx]) + numpy.sum(science_data[newy+2.0,newx]) + numpy.sum(science_data[newy, newx-2.0]) + numpy.sum(science_data[newy, newx+2.0])
				fbox5=numpy.sum(science_data[newy-2.0:newy+3.0, newx-2.0:newx+3.0])
				fbox6=numpy.sum(science_data[newy-5.0:newy+6.0, newx-5.0:newx+6.0])
				ch.write(str('Fboxes: Done')+'\n')
				reg.write(str(newx)+' '+str(newy)+'\n') #fake star region file

				scale_fac=scaling_factor_star[source_star] #scale factor
				back=background_array_star[source_star] #background

				#---Old area to be scaled---
				startx=int(sourcex-5.0)
				starty=int(sourcey-5.0)
				finx=int(sourcex+5.0)
				finy=int(sourcey+5.0)

				#---New area to have flux added---
				Nstartx=newx-5.0
				Nstarty=newy-5.0
				Nfinx=newx+5.0
				Nfiny=newy+5.0


				newdata=numpy.ones((10,10)) #Preparing a blank gird for scaled objects

				newdata[0:10,0:10]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object
				ch.write(str('New scaled Data: Added')+'\n')
				#print x,y
				#print 'New Data Shape: ', newdata.shape
				#print 'Science Shape: ', science_data[starty:finy,startx:finx].shape
				science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata #Modifying the science image


				f.write(str(xcord_star[source_star])+' '+str(ycord_star[source_star])+' '+str(alpha_array[source_star])+' '+str(delta_array[source_star])+' '+str(newx)+' '+str(newy)+' '+str(mag_array_star[source_star])+' '+str(mag_best_star[source_star])+' '+str(flux_array_star[source_star])+' '+str(ran_mag_star[source_star])+' '+str(ran_flux_star[source_star])+' '+str(background_array_star[source_star])+' '+str(scaling_factor_star[source_star])+' '+str(int(PTFFIELD))+' '+str(CCD_Num_star[source_star])+' '+str(x)+' '+str(y)+' '+str(Host_Alpha)+' '+str(Host_Dec)+' '+str(galaxy_mag_auto)+' '+str(galaxy_mag_best)+' '+str(galaxy_flux)+' '+str(galaxy_background)+' '+str(CXX)+' '+str(CYY)+' '+str(CXY)+' '+str(Host_Elongation)+' '+str(Host_Ellipticity)+' '+str(FR_02)+' '+str(FR_05)+' '+str(FR_09)+' '+str(fbox1)+' '+str(fbox2)+' '+str(fbox3)+' '+str(fbox4)+' '+str(fbox5)+' '+str(fbox6)+' '+str(gain)+' '+str(readnoise)+' '+str(MOONILLF)+' '+str(MoonRA)+' '+str(MoonDec)+' '+str(AIRMASS)+' '+str(seeing)+' '+str(ELLIP)+' '+str(MEDSKY)+' '+str(SKYSIG)+' '+str(zeropoint)+' '+str(LMT_MG)+' '+str(MJD)+'\n')

				num_of_fakes_all+=1
				ch.write(str('Host Galaxy: Done')+'\n')
				break

	for g in range(0,int(len(xcord_star)-int(len(xcord_star)*0.9))):
		Star_Location=True
		#print 'How Many Hostless: ',
		#print len(fake_star_array)
		source_star=fake_star_array.pop(random.randrange(0,len(fake_star_array)))
		ch.write(str('Hostless Source Star: Chosen')+'\n')
		while Star_Location==True:
			hostlessx=int(random.uniform(40.0,2008.0))
			hostlessy=int(random.uniform(40.0,4056.0))
			sourcex=xcord_star[source_star] #stars current x location
			sourcey=ycord_star[source_star] #stars current y location
			reg.write(str(hostlessx)+' '+str(hostlessy)+'\n') #fake star region file
			scale_fac=scaling_factor_star[source_star] #scale factor
			back=background_array_star[source_star] #background
			ch.write(str('Hostless Location: Chosen')+'\n')
			if galaxy_mask[hostlessy,hostlessx]==False and galareas[hostlessy,hostlessx]==False:
				#print 'Cant Go there<-- Hostless'
				continue
			else:

				r=40 #Radius of Mask
				ym,xm = numpy.ogrid[-hostlessy:resy-hostlessy, -hostlessx:resx-hostlessx] #Some clever numpy stuff
				mask = xm*xm + ym*ym <= r*r
				galaxy_mask[mask]=False
				ch.write(str('Hostless r and Mask: Done')+'\n')
				#---Old area to be scaled---
				startx=int(sourcex-10.0)
				starty=int(sourcey-10.0)
				finx=int(sourcex+10.0)
				finy=int(sourcey+10.0)

				#---New area to have flux added---
				Nstartx=hostlessx-10.0
				Nstarty=hostlessy-10.0
				Nfinx=hostlessx+10.0
				Nfiny=hostlessy+10.0


				fbox1=numpy.sum(science_data[hostlessy,hostlessx])
				fbox2=numpy.sum(science_data[hostlessy,hostlessx]) + numpy.sum(science_data[hostlessy-1.0,hostlessx]) + numpy.sum(science_data[hostlessy+1.0,hostlessx]) + numpy.sum(science_data[hostlessy,hostlessx-1.0]) + numpy.sum(science_data[hostlessy,hostlessx+1.0])
				fbox3=numpy.sum(science_data[hostlessy-1.0:hostlessy+2.0, hostlessx-1.0:hostlessx+2.0])
				fbox4=numpy.sum(science_data[hostlessy-1.0:hostlessy+2.0, hostlessx-1.0:hostlessx+2.0]) + numpy.sum(science_data[hostlessy-2.0,hostlessx]) + numpy.sum(science_data[hostlessy+2.0,hostlessx]) + numpy.sum(science_data[hostlessy, hostlessx-2.0]) + numpy.sum(science_data[hostlessy, hostlessx+2.0])
				fbox5=numpy.sum(science_data[hostlessy-2.0:hostlessy+3.0, hostlessx-2.0:hostlessx+3.0])
				fbox6=numpy.sum(science_data[hostlessy-5.0:hostlessy+6.0, hostlessx-5.0:hostlessx+6.0])
				ch.write(str('Hostless Fbox: Done')+'\n')

				newdata=numpy.ones((20,20)) #Preparing a blank gird for scaled objects

				newdata[0:20,0:20]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object
				#print x,y
				#print 'New Data Shape: ', newdata.shape
				#print 'Science Shape: ', science_data[starty:finy,startx:finx].shape
				science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata #Modifying the science image


				f.write(str(xcord_star[source_star])+' '+str(ycord_star[source_star])+' '+str(alpha_array[source_star])+' '+str(delta_array[source_star])+' '+str(hostlessx)+' '+str(hostlessy)+' '+str(mag_array_star[source_star])+' '+str(mag_best_star[source_star])+' '+str(flux_array_star[source_star])+' '+str(ran_mag_star[source_star])+' '+str(ran_flux_star[source_star])+' '+str(background_array_star[source_star])+' '+str(scaling_factor_star[source_star])+' '+str(int(PTFFIELD))+' '+str(CCD_Num_star[source_star])+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(99999.99)+' '+str(fbox1)+' '+str(fbox2)+' '+str(fbox3)+' '+str(fbox4)+' '+str(fbox5)+' '+str(fbox6)+' '+str(gain)+' '+str(readnoise)+' '+str(MOONILLF)+' '+str(MoonRA)+' '+str(MoonDec)+' '+str(AIRMASS)+' '+str(seeing)+' '+str(ELLIP)+' '+str(MEDSKY)+' '+str(SKYSIG)+' '+str(zeropoint)+' '+str(LMT_MG)+' '+str(MJD)+'\n')

				num_of_fakes_all+=1
				Star_Location=False
				ch.write(str('All Hostless Done')+'\n')
	hdulist_sci.writeto(science_image[0]+science_image[1]+'_fakesV'+str(vnum)+'.fits', output_verify='ignore', clobber=True) #Saving image after loop of 200 Stars is complete
	#j.close()
	reg.close()
	f.close()
	hdulist_sci.close()


	#print num_of_fakes_all, 'fake Stars Added to Galaxies and hostless in the Image: ', science_image[1]
	ch.write(str('Num of Fakes Added:')+' '+str(num_of_fakes_all)+'\n')
	ch.close()
	#Creating a Galaxy Mask Fits file

	'''
	galaxy_mask_float=galaxy_mask.astype(int)
	hdu=fits.PrimaryHDU(galaxy_mask_float)
	hdu.scale(type='int16')
	hdulist=fits.HDUList([hdu])

	#print hdulist.info()
	hdulist.writeto('Results_V'+str(vnum)+'/Fake_Star_Catalog/'+science_image[1]+'_GMask_V'+str(vnum)+'.fits', clobber=True, output_verify='ignore')
	hdulist.close()
	'''
def Execute(run):
	#print '!!!!!!', run
	science_image=run

	tstart=time.time()
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
		bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Astropy Could not Open the .fits file')+'\n')
		bad_images.close()
		#print 'Cant open Science'



		return

	hdulist_multi_sci.verify('fix')
	zeropoint=float(hdulist_multi_sci[0].header['UB1_ZP'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=55000.0 #float(hdulist_multi_sci[0].header['SATURATE'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])
	CCD_Num=float(hdulist_multi_sci[0].header['CCDID'])
	PTFFIELD=int(hdulist_multi_sci[0].header['PTFFIELD'])
	readnoise=float(hdulist_multi_sci[0].header['READNOI'])
	MOONILLF=float(hdulist_multi_sci[0].header['MOONILLF'])
	AIRMASS=float(hdulist_multi_sci[0].header['AIRMASS'])
	ELLIP=(hdulist_multi_sci[0].header['ELLIP'])
	if ELLIP=='NAN.0':
		bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: ELLIP has a NAN')+'\n')
		bad_images.close()
		#print science_image[0]+science_image[1], ' Has a NAN'
		return
	else:
		ELLIP=float(hdulist_multi_sci[0].header['ELLIP'])
	MEDSKY=float(hdulist_multi_sci[0].header['MEDSKY'])
	SKYSIG=float(hdulist_multi_sci[0].header['SKYSIG'])
	LMT_MG=float(hdulist_multi_sci[0].header['LMT_MG'])
	MJD=float(hdulist_multi_sci[0].header['OBSMJD'])
	MoonRA=float(hdulist_multi_sci[0].header['MOONRA'])
	MoonDec=float(hdulist_multi_sci[0].header['MOONDEC'])




	fake_stars= 60 #number of fake stars per image (integer please!)

	hdulist_multi_sci.close()


	Sextract(science_image,zeropoint,seeing,saturation,gain)

	catsize=Enough_Objects(science_image)
	if catsize==False:
			#print science_image, 'didn\'t have enough objects detected so it was moved to Results_V'+str(vnum)+'/Bad_Images/ and the newly created weight map, sex file and catalog have been deleted'
			bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
			bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Sextractor did not detect enough objects (<200)')+'\n')
			os.remove('Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat')
			return

	x, y, mag, flux, back, magnitude_best, alpha_sky, delta_sky = Selecting_Bright(science_image)
	#print 'Selecting Bright Done'
	xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, faint_fake, mag_best_star, alpha_array, delta_array =Scaling(science_image, x, y, mag, flux, back, zeropoint, fake_stars, CCD_Num, magnitude_best, alpha_sky, delta_sky)
	#print 'Scaling Done'
	mag_log=open('Results_V'+str(vnum)+'/Magnitude_Log_File.dat','a')
	#print 'Maglog Open'
	mag_log.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str(mag[0])+' '+str(mag[-1])+' '+str(faint_fake)+' '+str('22.5')+'\n')
	#print 'Maglogwrite'
	galareas=selecting_galaxies(science_image)
	#print 'Selected Galaxies'

	boxsize=[3,5,7] #This is redundant now, please do not consider this useful.
	add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star, mag_best_star, alpha_array, delta_array, zeropoint, seeing, saturation, gain, readnoise, MOONILLF, AIRMASS, ELLIP, MEDSKY, SKYSIG, LMT_MG, MJD, MoonRA, MoonDec, PTFFIELD, galareas)

	t_total=time.time()-tstart

	good_images=open('Results_V'+str(vnum)+'/Good_Images_V'+str(vnum)+'.dat','a')
	good_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str(t_total)+'\n')

	#Sub_ML_DB(science_image)
#-----------------------------------RUN PIPELINE------------------------------------------

def Run_All(masterlist):
	file_structure()
	'''
	all_fits=[] #Establishing an array to find the files
	#path=[]
	#fnames=[]
	for dirpath,dirname,filenames in os.walk(os.path.abspath('../../fakes')): #Traverses through a directory tree

		for file in filenames:
			fileex=os.path.splitext(file)[-1] #Splits the file name, [-1] means it will look at the extension
			if fileex== '.fits': #wanted all .fits files

				all_fits.append([dirpath, file])
	#print all_fits



	science_fits=[]
	for i in range(len(all_fits)):
		#fname=all_fits[1]
		ln=all_fits[i]
		fname=ln[1].split('.')
		#print fname


		if fname[-2]=='w':

			science_fits.append([ln[0]+str('/'), (os.path.splitext(ln[1])[0])])
	'''
	science_fits=[]
	k=open(masterlist)
	for line in k:
		koo=line.strip('.fits\n')
		kn=koo.split(' ')
		science_fits.append([kn[0]+str('/'),kn[1]])

	#print 'Science_Fits', science_fits
	bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
	bad_images.close()
	good_images=open('Results_V'+str(vnum)+'/Good_Images_V'+str(vnum)+'.dat','a')
	good_images.close()
	mag_log=open('Results_V'+str(vnum)+'/Magnitude_Log_File.dat','a')
	mag_log.close()
	mag_log_col=open('Results_V'+str(vnum)+'/Magnitude_Log_File.columns','w')
	mag_log_col.write(str('1. Path to Image')+'\n'+str('2. Brightest Source')+'\n'+str('3. Faintest Source Mag')+'\n'+str('4. Brightest Fake')+'\n'+str('5. Faintest Fake'))
	mag_log_col.close()


	t0=time.time()
	processors=multiprocessing.cpu_count()
	pool=Pool(processors)
	pool.map(Execute,science_fits)
	pool.close()

	print 'V'+str(vnum)+' took: ', time.time()-t0, 'seconds'
Run_All('Nam_List.dat')