import numpy as np
from astropy.io import fits
import glob2, os, subprocess, argparse
from astropy.time import Time

SEDMRAWPATH = os.getenv('SEDMRAWPATH')


def change_filenames(date):
    files = glob2.glob(f'{SEDMRAWPATH}/{date}/*.fits')
    for fl in files:
        hdr = fits.open(fl)[0].header
        dt = Time(hdr['DATE'], format='isot').isot[0:]
        if 'iXon' in hdr['CCDNAME']:
            prefix = 'rc'
        else:
            prefix = 'ifu'
        newfname = prefix + dt.replace('-', '').replace('T', '_').replace(':', '_')[:-4] + '.fits'
        if 'ifu' in newfname:
            subprocess.call(f'mv {fl} {SEDMRAWPATH}/{date}/{newfname}', shell=True)
        else:
            subprocess.call(f'mv {fl} {SEDMRAWPATH}/{date}/{newfname}', shell=True)


def add_header_keywords(filename,sciexptime=60):
    if '.fz' in filename:
        hdrnum = 1
    else:
        hdrnum = 0
    hdr = fits.open(filename)[hdrnum].header
    date = Time(hdr['DATE'], format='isot')
    fits.setval(filename, 'DATE-OBS', value=date.isot)
    fits.setval(filename, 'JD', value=date.jd)
    fits.setval(filename, 'MJD-OBS', value=date.mjd)
    # fits.setval(filename, 'ADCSPEED', value=2.0)
    # fits.setval(filename, 'FILENAME', value=filename)
    if 'bias' in hdr['IMGTYPE']:
        fits.setval(filename, 'EXPTIME', value=0)
        fits.setval(filename, 'DOMEST', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='Calib: bias of')
        fits.setval(filename, 'NAME', value='Calib: bias of')
        # if hdr['MODE_NUM'] != 0:
        #     fits.setval(filename, 'ADCSPEED', value=0.1)
    if 'flat' in hdr['IMGTYPE']:
        fits.setval(filename, 'EXPTIME', value=5)
        fits.setval(filename, 'DOMEST', value='closed')
        fits.setval(filename, 'LAMPCUR', value=1.0)
        fits.setval(filename, 'OBJECT', value='Calib: dome of')
        fits.setval(filename, 'NAME', value='Calib: dome of')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_xe' in hdr['FILENAME']):
        fits.setval(filename, 'EXPTIME', value=120)
        fits.setval(filename, 'DOMEST', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='Calib: Xe of')
        fits.setval(filename, 'NAME', value='Calib: Xe of')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_hg' in hdr['FILENAME']):
        fits.setval(filename, 'EXPTIME', value=120)
        fits.setval(filename, 'DOMEST', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='Calib: Hg of')
        fits.setval(filename, 'NAME', value='Calib: Hg of')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_cd' in hdr['FILENAME']):
        fits.setval(filename, 'EXPTIME', value=120)
        fits.setval(filename, 'DOMEST', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='Calib: Cd of')
        fits.setval(filename, 'NAME', value='Calib: Cd of')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('object' in hdr['IMGTYPE']):
        fits.setval(filename, 'DOMEST', value='open')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'EXPTIME', value=sciexptime)
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    print('Added relevant keywords to ', filename)


def rotate_image(filename):
    if '.fz' in filename:
        hdrnum = 1
    else:
        hdrnum = 0
    hdu = fits.open(filename)[hdrnum]
    hdr = hdu.header
    if 'ROTFLIP' in hdr:
        # fits.setval(filename,'ROTFLIP',value=1)
        return
    data = np.fliplr(np.rot90(hdu.data))
    fits.HDUList([fits.PrimaryHDU(data, header=hdr)]).writeto(filename, overwrite=True)
    fits.setval(filename, 'ROTFLIP', value=1)
    print('Rotated and flipped ', filename)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--date', default=Time.now().iso[0:10].replace('-', ''), help='Date to process YYMMDD')
    args = parser.parse_args()

    # Move darks to unused directory
    if glob2.glob(f'{SEDMRAWPATH}/{args.date}/unused') == []:
        os.mkdir(f'{SEDMRAWPATH}/{args.date}/unused')
    darkfiles = glob2.glob(f'{SEDMRAWPATH}/{args.date}/speccal_dk*')
    for df in darkfiles:
        subprocess.call(f'mv {df} {SEDMRAWPATH}/{args.date}/unused/', shell=True)
    # change_filenames(args.date)
    files = sorted(glob2.glob(f'{SEDMRAWPATH}/{args.date}/speccal*.fits.fz'))
    for fl in files:
        subprocess.call(f"funpack -F {fl}", shell=True)
        add_header_keywords(fl[0:-3])
        rotate_image(fl[0:-3])

    # change_filenames(args.date) for sci
    scifiles = sorted(glob2.glob(f'{SEDMRAWPATH}/{args.date}/spec_*.fits'))
    for fl in scifiles:
        print(fl)
        sciexptime = float(input('Input exptime for this file: '))
        add_header_keywords(fl,sciexptime)
        rotate_image(fl)
