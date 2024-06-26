# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""
import subprocess
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

import numpy as np
import aplpy
try:
    import coordinates_conversor
except ImportError:
    import drprc.coordinates_conversor as coordinates_conversor
# try:
#     import fitsutils
# except ImportError:
#     import drprc.fitsutils as fitsutils
import datetime
import os
import sys
import glob
import argparse
import json
import sedmpy_version

from matplotlib import pylab as plt
plt.switch_backend('Agg')

try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, 'sedmconfig.json')
with open(configfile) as config_file:
    sedm_cfg = json.load(config_file)

# Default paths
_rcpath = sedm_cfg['paths']['photpath']
_reduxpath = sedm_cfg['paths']['reduxpath']
_srcpath = sedm_cfg['paths']['srcpath']
_altrcpath = sedm_cfg['paths']['rawpath']


def finder(acqfile, imhdr, findername, searchrad=0.2/60.):
    acqhdr = fits.open(acqfile)[0].header
    kra = imhdr['RA']
    kdec = imhdr['DEC']
    ocoord = SkyCoord(ra=kra,dec=kdec,unit=(u.hourangle,u.deg))
    ora = ocoord.ra.deg
    odec = ocoord.dec.deg
    utc = acqhdr['DATE']
    hdulist = fits.open(acqfile)[0]
    img = hdulist.data * 1.            
    img = img.T

    wcs = WCS(acqhdr)
    print(wcs)

    target_pix = wcs.wcs_world2pix([(np.array([ora, odec], np.float_))], 1)[0]
    print(target_pix)
    # check bounds
    good_coords = False
    if (target_pix[0] < 0 or target_pix[0] > img.shape[0] or
       target_pix[1] < 0 or target_pix[1] > img.shape[1]):
        print("ERROR - target outside acquisition image: x,y = ", target_pix)
        target_pix = np.array([[150., 150.]], dtype=np.float)
        ra, dec = wcs.wcs_pix2world(target_pix, 1)[0]
        target_pix = target_pix[0]
        findername = findername.replace(".png", "_failed.png")
    else:
        good_coords = True
        ra, dec = ora, odec

    corner_pix = wcs.wcs_world2pix([(np.array([ra, dec + searchrad],
                                              np.float_))], 1)[0]
    dx = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))

    # imgslice = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx,
    #                int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]
    imgslice_target = img[int(target_pix[0])-dx:int(target_pix[0])+dx,
                          int(target_pix[1])-dx:int(target_pix[1])+dx]

    # Maybe the object has moved out of this frame.
    # In this case, make the finder larger.
    x, y = imgslice_target.shape

    print(img.shape, target_pix, corner_pix, dx, int(target_pix[0])-dx,
          int(target_pix[0])+dx, int(target_pix[1])-dx, int(target_pix[1])+dx)

    if (x < 2*dx-1) or (y < 2*dx-1):
        imgslice_target = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx,
                              int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]

    # zmin, zmax = zscale.zscale()
    zmin = np.percentile(imgslice_target.flatten(), 5)
    zmax = np.percentile(imgslice_target.flatten(), 98.5)
   
    print("Min: %.1f, max: %.1f" % (zmin, zmax))
    gc = aplpy.FITSFigure(acqfile, figsize=(10, 9), north=True)
    gc.show_grayscale(vmin=zmin, vmax=zmax, smooth=1, kernel="gauss")
    gc.add_scalebar(0.1/60.)
    gc.scalebar.set_label('10 arcsec')
    gc.scalebar.set_color('white')
    gc.recenter(ra, dec, searchrad)

    ras = np.array([ra, ra])
    decs = np.array([dec, dec])
    dxs = np.array([0, searchrad/10 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/10, 0])

    if good_coords:
        gc.show_arrows(ras, decs, dxs, dys, edgecolor="red", facecolor="red",
                       head_width=0)

    ras = np.array([ra+searchrad*0.7 / np.cos(np.deg2rad(dec)),
                    ra+searchrad*0.7 / np.cos(np.deg2rad(dec))])
    decs = np.array([dec-searchrad*0.9, dec-searchrad*0.9])
    dxs = np.array([0, searchrad/5 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/5, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="white", facecolor="white")
    gc.add_label(ras[0]+dxs[0]*1.1, decs[0]+dys[0]*1.1, 'N', relative=False,
                 color="white", horizontalalignment="center")
    gc.add_label(ras[1]+dxs[1]*1.1, decs[1]+dys[1]*1.1, 'E', relative=False,
                 color="white", horizontalalignment="center")

    img_name, img_filter = os.path.basename(finderpath).split('.')[0].split('_')[-2:]
    if img_filter=='cl':
        img_filter = 'clear'
    gc.add_label(0.05, 0.95, 'Object: %s' % img_name, relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.9, 'Filter: %s' % img_filter, relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.85, 'Coordinates: RA=%s DEC=%s' %
                 (kra, kdec), relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.8, "UTC: %s" % utc, relative=True,
                 color="white", horizontalalignment="left")
    if not good_coords:
        c = SkyCoord(ra=ra,dec=dec,unit=(u.deg)).to_string('hmsdms').split()
        gc.add_label(0.05, 0.75, 'Acquired coords: RA=%s DEC=%s' %
                     (c[0], c[1]), relative=True,
                     color="red", horizontalalignment="left")
        gc.add_label(0.05, 0.70, 'FAILED ACQUISITION', relative=True,
                     color="red", horizontalalignment="left")
    
    gc.save(findername)
    print("Created %s" % findername)
    

def simple_finder(acqfile, findername):

    hdulist = fits.open(acqfile)[0]
    img = hdulist.data * 1.            
    img = img.T
    img = img[0:500, 0:500]
    newimg = img
    # ndimage.filters.gaussian_filter(img, 1, order=0, mode='constant',
    # cval=0.0, truncate=20.0)

    # name = fitsutils.get_par(myfile, "NAME")
    # filter = fitsutils.get_par(myfile, "FILTER")

    # zmin, zmax = zscale.zscale()
    zmin = np.percentile(newimg.flatten(), 10)
    zmax = np.percentile(newimg.flatten(), 99)
    plt.figure(figsize=(10, 9))
    plt.imshow(newimg, origin="lower", cmap=plt.get_cmap('gray'),
               vmin=zmin, vmax=zmax)

    plt.savefig(findername)

    print("Created ", findername)


def simple_finder_astro(acqfile, imhdr, findername, searchrad=28./3600):

    hdulist = fits.open(acqfile)[0]
    img = hdulist.data * 1.            

    # name = fitsutils.get_par(myfile, "NAME")
    # filter = fitsutils.get_par(myfile, "FILTER")
    ocoord = SkyCoord(ra=imhdr['RA'], dec=imhdr['DEC'], unit=(u.hourangle,u.deg))
    ra = ocoord.ra.deg
    dec = ocoord.dec.deg
    
    wcs = WCS(hdulist.header)

    target_pix = wcs.wcs_world2pix([(np.array([ra, dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_world2pix([(np.array([ra+searchrad, dec+searchrad],
                                              np.float_))], 1)[0]
    targ_x = int(target_pix[0])
    targ_y = int(target_pix[1])
    # Size of the finder in pixels
    
    dx = int(np.abs(np.ceil(corner_pix[0] - target_pix[0])))
    dy = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    # size = int( (searchrad/0.394)/2)
    
    # zmin, zmax = zscale.zscale()
    newimg = img[targ_x-dx: targ_x+dx, targ_y-dy: targ_y+dy]

    zmin = np.percentile(newimg.flatten(), 5)
    zmax = np.percentile(newimg.flatten(), 98.5)
    
    print("X %d Y %d Size %d, %d zmin=%.2f zmax=%.2f. Size = %s" %
          (targ_x, targ_y, dx, dy, zmin, zmax, newimg.shape))

    from astropy.visualization.wcsaxes import SphericalCircle

    plt.figure(figsize=(10, 9))
    ax = plt.subplot(projection=wcs)
    ax.imshow(np.flip(newimg, axis=0),
              origin="lower", cmap=plt.get_cmap('gray'), vmin=zmin, vmax=zmax)
    r = SphericalCircle((ra * u.deg, dec * u.deg), 5./3600 * u.degree,
                        edgecolor='red', facecolor='none',
                        transform=ax.get_transform('fk5'))
    ax.add_patch(r)

    # ax = plt.gca()
    # ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=20,
    #       edgecolor='red', facecolor='none')

    # plt.plot(dy, dx, "+", color="r", ms=20, mfc=None, mew=2)
    # plt.plot(Y, X, "+", color="r", ms=20, mfc=None, mew=2)
    # plt.xlim(Y-dy, Y+dy, X-dx, X+dx)

    plt.savefig(findername)

    print("Created ", findername)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    
    Creates a finder chart for every acquisition image in the folder
    specified as a parameter.
    As a final step, it copies the acquisition image to the "agn" machine
    to visualize it.
        
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-d', '--rcdir', type=str, dest="rcdir",
                        help='Directory with rc images from tonight.',
                        default=None)
    parser.add_argument('-r', '--reduxdir', type=str, dest="reduxdir",
                        help='Directory with reduced ifu images from tonight.',
                        default=None)
    parser.add_argument('-i', '--imfile', type=str, dest="imfile",
                        help='IFU image that requires a finder',
                        default=None)
    
    args = parser.parse_args()

    imfile = args.imfile

    if imfile is None:
        print('No imfile given, exiting')
        sys.exit(1)

    timestamp = imfile.split('/')[-2]
    rcdir = os.path.join(_rcpath, timestamp)
    if not os.path.isdir(rcdir):
        os.mkdir(rcdir)
    reduxdir = '/'.join(imfile.split('/')[0:-1])
    imhdr = fits.open(imfile)[0].header
    objnam = imhdr['OBJECT']
    imtime = Time(imhdr['JD'],format='jd')
    imreadtime = Time(imhdr['DATE'],format='isot')
        # objnam = fitsutils.get_par(imfile, "OBJECT").split()[0]
        # if 'STD' in objnam:
        #     objnam = objnam.split('STD-')[-1].split()[0]
    # else:
    #     rcdir = args.rcdir
    #     reduxdir = args.reduxdir
    #     objnam = None
    #     if rcdir is None:
    #         timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
    #         timestamp = timestamp.split("T")[0].replace("-", "")
    #         rcdir = os.path.join(_rcpath, timestamp)
    #         if not os.path.isdir(rcdir):
    #             rcdir = os.path.join(_altrcpath, timestamp)
    #     else:
    #         timestamp = os.path.basename(os.path.abspath(rcdir))
    #
    #     if reduxdir is None:
    #         reduxdir = os.path.join(_reduxpath, timestamp)

    os.chdir(reduxdir)
    print("Making finder for object: %s" % objnam)
    print("Changed to directory where the reduced data is: %s" % reduxdir)
    print("Getting acquisition images from directory: %s" % rcdir)

    if not (os.path.isdir("finders")):
        os.makedirs("finders")

    # We gather all point images to locate the last acquisition.
    files = sorted(glob.glob(os.path.join(rcdir, f"point*{objnam}*.new")))
    # If no files in rcdir for the given object, download them from robo machine
    if len(files)==0:
        # First try downloading from /Data
        subprocess.call(
            f'rsync -av -e "ssh -p22222" sedm@140.252.53.120:/Data/{timestamp}/point*{objnam}*.new {rcdir}/',
            shell=True)
        # Check if files downloaded, if not download from /Backup/Data
        if len(glob.glob(os.path.join(rcdir, f"point*{objnam}*.new"))) == 0:
            subprocess.call(
                f'rsync -av -e "ssh -p22222" sedm@140.252.53.120:/Backup/Data/{timestamp}/point*{objnam}*.new {rcdir}/',
                shell=True)
        # Check again, if still no files downloaded then print error
        if len(glob.glob(os.path.join(rcdir, f"point*{objnam}*.new"))) == 0:
            print(f'Could not find pointing files for {objnam}')
    filesacq = []
    for fl in files:
        dt,tm = os.path.basename(fl).split('_')[1:3]
        filetime = Time(f'{dt[0:4]}-{dt[4:6]}-{dt[6:8]}T{tm[0:2]}:{tm[2:4]}:{tm[4:]}',format='isot')
        if (filetime > imtime) & (filetime < imreadtime):
            filesacq.append(fl)

    # for f in files:
    #     try:
    #         ff = pf.open(f)
    #     except OSError:
    #         print("WARNING - corrupt fits file: %s" % f)
    #         continue
    #     if "IMGTYPE" in ff[0].header:
    #         imgtype = ff[0].header["IMGTYPE"]
    #     else:
    #         imgtype = ''
    #     if "OBJECT" in ff[0].header:
    #         obj = ff[0].header["OBJECT"]
    #     else:
    #         obj = ''
    #     if (imgtype.upper() == "ACQUISITION" or "ACQ" in imgtype.upper() or
    #        "ACQ" in obj) and ("TEST" not in imgtype.upper()):
    #         if objnam:
    #             if objnam in obj:
    #                 filesacq.append(f)
    #         else:
    #             filesacq.append(f)
    #     else:
    #         if 'OBJECT' in ff[0].header:
    #             if 'finding' in ff[0].header['OBJECT'] and \
    #                     objnam in ff[0].header['OBJECT']:
    #                 filesacq.append(f)

    n_acq = len(filesacq)
    print("Found %d files for finders:\n%s" % (n_acq, filesacq))

    n_find = 0
    for f in filesacq:
        n_find += 1
        print("Trying finder %d of %d: %s" % (n_find, n_acq, f))
        # try:
        #     objnam = fitsutils.get_par(f, "OBJECT")
        # except:
        #     print('There is no object in this file %s. Skipping'
        #           ' and moving to the next file.' % f)
        #     continue

        # We generate only one finder for each object.
        name = f.split('/')[-1].split(".")[0]
        # if 'finding' in objnam:
        #     filt = objnam.split()[1].split('[')[-1].split(']')[0]
        #     objnam = objnam.split()[0].split('STD-')[-1]
        # else:
        #     objnam = objnam.split()[0]
        #     filt = fitsutils.get_par(f, "FILTER")
        filt = 'cl'
        finderplotf = 'finder_%s_%s_%s.png' % (name, objnam, filt)
        finderpath = os.path.join(reduxdir, os.path.join("finders/",
                                                         finderplotf))
        # Check if it was already done
        if not os.path.isfile(finderpath):
            # # Check for existing astrometry file
            # astrof = os.path.join(rcdir, "a_%s" % f.split('/')[-1])
            # if not os.path.exists(astrof):
            #     # link rc image into reduxdir
            #     dest = os.path.join(reduxdir, f.split('/')[-1])
            #     if os.path.isfile(dest):
            #         print("RC ACQ image already exists: %s" % dest)
            #     else:
            #         os.symlink(f, dest)
            #
            #     print("Solving astrometry", dest)
            #     # Solving for astrometry
            #     astrof = dest.replace(".fits", "_astrom.fits").replace(".gz",
            #                                                            "")
            #     if not os.path.exists(astrof):
            #         returncode = subprocess.call([_srcpath+'bin/do_astrom',
            #                                       dest])
            #         if returncode != 0:
            #             print("Astrometry failed, perform median subtraction")
            #             returncode = subprocess.call(
            #                 [_srcpath+'bin/spy',
            #                  _srcpath+'drpifu/med_sub.py', '-i',
            #                  f.split('/')[-1]])
            #             if returncode != 0:
            #                 print("Astrometry failed for %s, "
            #                       "skipping finder %s" % (dest, finderpath))
            #                 continue
            #             returncode = subprocess.call(
            #                 [_srcpath+'bin/do_astrom', dest])
            #         if returncode != 0:
            #             print("Astrometry failed for %s, skipping finder %s" %
            #                   (dest, finderpath))
            #             continue
            #     else:
            #         print("Astrometry file already exists: %s" % astrof)
            #     # Check results
            #     if not os.path.isfile(astrof):
            #         print("Astrometry results not found %s" % astrof)
            #         continue

            print("Using astrometry in %s" % f)
            # try:
            finder(f, imhdr, finderpath)
            # except ValueError:
            #     print("Bad astrometry for this file: %s" % f)
            #     continue
            # except AttributeError:
            #     print("Error when generating the finder for file %s" % f)
            #     print(sys.exc_info()[0])
            #     simple_finder_astro(f, imhdr, finderpath)
            #
            # except:
            #     print("Error when generating the finder for file %s. "
            #           "Probably montage is broken." % f)
            #     print(sys.exc_info()[0])
            #     simple_finder_astro(f, imhdr, finderpath)
        else:
            print("Finder already exists: %s" % finderpath)
