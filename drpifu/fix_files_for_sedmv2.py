import numpy as np
from astropy.io import fits
import glob2, os, subprocess, argparse
from astropy.time import Time
from astroplan import Observer
from astropy.coordinates import SkyCoord
import astropy.units as u

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


def add_header_keywords(filename):
    if '.fz' in filename:
        hdrnum = 1
    else:
        hdrnum = 0
    hdr = fits.open(filename)[hdrnum].header
    if 'JD' in hdr.keys():
        return
    date = Time(hdr['DATE'], format='isot')
    fits.setval(filename, 'DATE-OBS', value=date.isot)
    fits.setval(filename, 'JD', value=date.jd)
    fits.setval(filename, 'MJD-OBS', value=date.mjd)
    if 'bias' in hdr['IMGTYPE']:
        fits.setval(filename, 'EXPTIME', value=0)
        fits.setval(filename, 'DOMESTAT', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='bias')
        fits.setval(filename, 'NAME', value='bias')
    if 'flat' in hdr['IMGTYPE']:
        fits.setval(filename, 'EXPTIME', value=5)
        fits.setval(filename, 'DOMESTAT', value='closed')
        fits.setval(filename, 'LAMPCUR', value=1.0)
        fits.setval(filename, 'OBJECT', value='dome')
        fits.setval(filename, 'NAME', value='dome')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_xe' in filename):
        fits.setval(filename, 'EXPTIME', value=30)
        fits.setval(filename, 'DOMESTAT', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='arc')
        fits.setval(filename, 'NAME', value='Xe')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_hg' in filename):
        fits.setval(filename, 'EXPTIME', value=30)
        fits.setval(filename, 'DOMESTAT', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='arc')
        fits.setval(filename, 'NAME', value='Hg')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('speccal_cd' in filename):
        fits.setval(filename, 'EXPTIME', value=30)
        fits.setval(filename, 'DOMESTAT', value='closed')
        fits.setval(filename, 'LAMPCUR', value=0.0)
        fits.setval(filename, 'OBJECT', value='arc')
        fits.setval(filename, 'NAME', value='Cd')
        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    if ('object' in hdr['IMGTYPE']):
        if 'EXPTIME' not in hdr.keys():
            if float(hdr['EXPOSURE']) < 1e-2:
                sciexptime = float(input('Input exptime for this file: '))
            else:
                sciexptime = float(hdr['EXPOSURE'])
            fits.setval(filename, 'EXPTIME', value=sciexptime)
        fits.setval(filename, 'DOMESTAT', value='open')
        fits.setval(filename, 'LAMPCUR', value=0.0)

        objname = hdr['QCOMMENT'].split()[0]
        if 'STD-' in objname:
            imtype = 'standard'
        else:
            imtype = 'science'
        fits.setval(filename, 'OBJECT', value=objname)
        fits.setval(filename, 'NAME', value=objname)
        fits.setval(filename, 'IMGTYPE', value=imtype)

        loc = Observer.at_site('Kitt Peak')
        airmass = loc.altaz(time=date, target=SkyCoord(ra=hdr['RAD'],dec=hdr['DECD'],unit=(u.deg))).secz
        fits.setval(filename, 'AIRMASS', value=airmass.value)

        fits.setval(filename, 'CRVAL1', value=155.4749009102372)
        fits.setval(filename, 'CRVAL2', value=88.6437849176057)
    print('Added relevant keywords to ', filename)
    return


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
        subprocess.call(f"cp {fl} {SEDMRAWPATH}/{args.date}/unused/", shell=True)
        subprocess.call(f"funpack -F {fl}", shell=True)

    files = sorted(glob2.glob(f'{SEDMRAWPATH}/{args.date}/speccal*.fits'))
    for fl in files:
        add_header_keywords(fl)
        rotate_image(fl)

    # change_filenames(args.date) for sci
    scifiles = sorted(glob2.glob(f'{SEDMRAWPATH}/{args.date}/sedm2_*.fits*'))
    for fl in scifiles:
        print(fl)
        if ".fz" in fl:
            subprocess.call(f"cp {fl} {SEDMRAWPATH}/{args.date}/unused/", shell=True)
            subprocess.call(f"funpack -F {fl}", shell=True)
            fl = fl[0:-3]
        add_header_keywords(fl)
        rotate_image(fl)
