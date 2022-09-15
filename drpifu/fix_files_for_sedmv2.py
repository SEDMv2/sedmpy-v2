import numpy as np 
from astropy.io import fits
import glob2,os,subprocess,argparse
from astropy.time import Time

SEDMREDUXPATH = '/media/yashvi/Data/Work/SEDMv2/sedmv2/raw'

def change_filenames(date):
	files = glob2.glob(f'{SEDMREDUXPATH}/{date}/*.fits')
	for fl in files:
		if('ifu' in fl):
			continue
		dt = Time(fits.open(fl)[0].header['DATE'],format='isot').isot[0:]
		newfname = 'ifu'+dt.replace('-','').replace('T','_').replace(':','_')[:-4]+'.fits'
		subprocess.call(f'mv {fl} {SEDMREDUXPATH}/{date}/{newfname}',shell=True)

def add_header_keywords(filename):
	hdr = fits.open(filename)[0].header
	if 'MJD-OBS' in hdr:
		return
	date = Time(hdr['DATE'],format='isot')
	fits.setval(filename,'DATE-OBS',value=date.isot)
	fits.setval(filename,'JD',value=date.jd)
	fits.setval(filename,'MJD-OBS',value=date.mjd)
	fits.setval(filename,'ADCSPEED',value=2.0)
	fits.setval(filename,'FILENAME',value=filename)
	if 'bias' in hdr['IMGTYPE']:
		fits.setval(filename,'EXPTIME',value=0)
		fits.setval(filename,'DOMEST',value='closed')
		fits.setval(filename,'LAMPCUR',value=0.0)
		fits.setval(filename,'OBJECT',value='Calib: bias of')
		fits.setval(filename,'NAME',value='Calib: bias of')
		if hdr['MODE_NUM']!=0:
			fits.setval(filename,'ADCSPEED',value=0.1)
	if 'flat' in hdr['IMGTYPE']:
		fits.setval(filename,'EXPTIME',value=120)
		fits.setval(filename,'DOMEST',value='closed')
		fits.setval(filename,'LAMPCUR',value=1.0)
		fits.setval(filename,'OBJECT',value='Calib: dome of')
		fits.setval(filename,'NAME',value='Calib: dome of')
	if ('test' in hdr['IMGTYPE'] or 'Test' in hdr['IMGTYPE']) and ('xenon' in hdr['ORIGNAME'] or 'Xe' in hdr['ORIGNAME']):
		fits.setval(filename,'EXPTIME',value=60)
		fits.setval(filename,'DOMEST',value='closed')
		fits.setval(filename,'LAMPCUR',value=0.0)
		fits.setval(filename,'OBJECT',value='Calib: Xe of')
		fits.setval(filename,'NAME',value='Calib: Xe of')
	if ('test' in hdr['IMGTYPE'] or 'Test' in hdr['IMGTYPE']) and ('Hg' in hdr['ORIGNAME']):
		fits.setval(filename,'EXPTIME',value=60)
		fits.setval(filename,'DOMEST',value='closed')
		fits.setval(filename,'LAMPCUR',value=0.0)
		fits.setval(filename,'OBJECT',value='Calib: Hg of')
		fits.setval(filename,'NAME',value='Calib: Hg of')
	if ('test' in hdr['IMGTYPE'] or 'Test' in hdr['IMGTYPE']) and ('Cd' in hdr['ORIGNAME']):
		fits.setval(filename,'EXPTIME',value=60)
		fits.setval(filename,'DOMEST',value='closed')
		fits.setval(filename,'LAMPCUR',value=0.0)
		fits.setval(filename,'OBJECT',value='Calib: Cd of')
		fits.setval(filename,'NAME',value='Calib: Cd of')
	if ('object' in hdr['IMGTYPE']):
		fits.setval(filename,'DOMEST',value='open')
		fits.setval(filename,'LAMPCUR',value=0.0)
	print('Added relevant keywords to ',filename)

def rotate_image(filename):
	hdu = fits.open(filename)[0]
	hdr = hdu.header 
	if 'ROTFLIP' in hdr:
		# fits.setval(filename,'ROTFLIP',value=1)
		return
	data = np.fliplr(np.rot90(hdu.data))
	fits.HDUList([fits.PrimaryHDU(data,header=hdr)]).writeto(filename,overwrite=True)
	fits.setval(filename,'ROTFLIP',value=1)
	print('Rotated and flipped ',filename)


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--date',default=Time.now().iso[0:10].replace('-',''),help='Date to process YYMMDD')
	args = parser.parse_args()

	change_filenames(args.date)
	files = sorted(glob2.glob(f'{SEDMREDUXPATH}/{args.date}/ifu*.fits'))
	for fl in files:
		add_header_keywords(fl)
		rotate_image(fl)