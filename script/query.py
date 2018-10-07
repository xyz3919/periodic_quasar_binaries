import os, sys, string,urllib,copy
from math import log10, radians, pi,cos,sin
import time, traceback, datetime
import threading
import StringIO
import numpy

global limit_per_min, recent_query_list, do_limit
limit_per_min = 19.99
recent_query_list = []
do_limit = True


class Stripe82:
	"""
        Written by J. Bloom (jsbloom@astro.berkeley.edu), Aug. 2010.
        Version 1.0
        
        Edited by Yu-Ching(Tony) Chen (ycchen@illinois.ede), Oct 2018.
	"""
    
	dr_url="http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp"
	formats = ['csv','xml','html']
	def_fmt = "csv"
	verbose = True
	outdir =  "data/SDSS/"
	
	def __init__(self):
		if not os.path.exists(self.outdir):
			os.mkdir(self.outdir)
		
	def _filtercomment(self,sql):
		"Get rid of comments starting with --"
		fsql = ''
		for line in sql.split('\n'):
			fsql += line.split('--')[0] + ' ' + os.linesep;
		return fsql

	def q(self,ra,dec,name,dist=5.0,clobber=False,idn=0):
		outname = self.outdir + name +".csv"
		if not clobber and os.path.exists(outname):
			print "file %s already exists" % outname
			return
			 
		sss = """SELECT p.objid,p.ra,p.dec,dbo.fdistancearcmineq(p.ra,p.dec,%f,%f)*60 as darcsec,
						p.psfmag_g,p.psfMagErr_g,
						p.psfmag_r,p.psfMagErr_r,
						p.psfmag_i,p.psfMagErr_i,
		                                p.psfmag_z,p.psfMagErr_z,
		                case
						  when p.parentID != 0 then "True"
						  else "False"
						end as "deblended",
		                f.mjd_r,f.mjd_g,f.mjd_i,f.mjd_z,f.nobjects,f.nstars
			     FROM PhotoObjAll p
				 JOIN field as f on p.fieldid = f.fieldid
		         WHERE ra between %f and %f
		           and dec between %f and %f
		           order by f.mjd_r
		      """ % (ra, dec, ra - cos(radians(dec))*0.5*dist/3600.0, ra + cos(radians(dec))*0.5*dist/3600.0,\
               dec - 0.5*dist/3600.0, dec + 0.5*dist/3600.0)		
		fff = self._query(sss)
		line = fff.readline()
		if line.startswith("ERROR") or line.startswith("No objects"):
			if self.verbose:
				print line
				print fff.readlines()
			return
		## get the keys
		kkk = line.split(",")
		print kkk
		f = open(outname,"w")
		f.write("# stripe82 phot for (ra,dec) = (%f,%f)\n" % (ra,dec))
		f.write("# time run %s (utc)\n" % str(datetime.datetime.utcnow()))
		f.write(",".join(kkk))
		ooo = fff.readlines()
		print " .... %i data points" % len(ooo)
		f.writelines(ooo)
		return

        def q_orig(self,ra,dec,dist=5.0,clobber=False,idn=0):
                outname = self.outdir + "sdss_%i_%f%s%f.phot" % (idn,ra,"-" if dec < 0 else "+",abs(dec))
                if not clobber and os.path.exists(outname):
                        print "file %s already exists" % outname
                        return

                sss = """SELECT p.objid,p.ra,p.dec,dbo.fdistancearcmineq(p.ra,p.dec,%f,%f)*60 as darcsec,
                                p.dered_r,
                                                p.psfmag_u,p.psfMagErr_u,p.dered_u,p.u,p.err_u,
                                                p.psfmag_g,p.psfMagErr_g,p.dered_g,p.g,p.err_g,
                                                p.psfmag_r,p.psfMagErr_r,p.dered_r,p.r,p.err_r,
                                                p.psfmag_i,p.psfMagErr_i,p.dered_i,p.i,p.err_i,
                                p.psfmag_z,p.psfMagErr_z,p.dered_z,p.z,p.err_z,
                                                p.petroR90_r,p.petroR90Err_r,
                                case
                                                  when p.parentID != 0 then "True"
                                                  else "False"
                                                end as "deblended",
                                f.mjd_r,f.mjd_g,f.mjd_i,f.nobjects,f.nstars
                             FROM PhotoObjAll p
                                 JOIN field as f on p.fieldid = f.fieldid
                         WHERE ra between %f and %f
                           and dec between %f and %f
                           order by f.mjd_r
                      """ % (ra, dec, ra - cos(radians(dec))*0.5*dist/3600.0, ra + cos(radians(dec))*0.5*dist/3600.0,\
               dec - 0.5*dist/3600.0, dec + 0.5*dist/3600.0)
                fff = self._query(sss)
                line = fff.readline()
                if line.startswith("ERROR") or line.startswith("No objects"):
                        if self.verbose:
                                print line
                                print fff.readlines()
                        return
                ## get the keys
                kkk = line.split(",")
                print kkk
                f = open(outname,"w")
                f.write("# stripe82 phot for (ra,dec) = (%f,%f)\n" % (ra,dec))
                f.write("# time run %s (utc)\n" % str(datetime.datetime.utcnow()))
                f.write("#" + ",".join(kkk))
                ooo = fff.readlines()
                print " .... %i data points" % len(ooo)
                f.writelines(ooo)
                return
	
	def runnat(self,f="qso_not_in_sesar.txt"):
		a= numpy.loadtxt(f)
		for i in range(len(a[:,0])):
			ra,dec = a[i,0],a[i,1]
			print i, ra, dec
			self.q(ra,dec)
			
	def _query(self,sql,url=dr_url,fmt=def_fmt,wait_period=62):
		"Run query and return file object"
		global limit_per_min, recent_query_list, do_limit
		if do_limit:
			recent_query_list.append((time.time(),sql))
			## cull out all the old calls
			recent_query_list = [x for x in recent_query_list if time.time() - x[0] < wait_period]
			if len(recent_query_list) > limit_per_min:
				## ug, we've got to wait
				tmp = [time.time() - x[0] for x in recent_query_list if time.time() - x[0] < wait_period]
				wait = wait_period - max(tmp) + 1
				if self.verbose:
					print "Date: %s Query length is %i in the last %f sec" % (str(datetime.datetime.now()),len(recent_query_list) - 1 , wait_period)  
					print "waiting %f sec, so as not to block the SDSS query %s" % (wait,sql)
				time.sleep(wait)

		fsql = self._filtercomment(sql)
		params = urllib.urlencode({'cmd': fsql, 'format': fmt})
		try:
			return urllib.urlopen(url+'?%s' % params) 
		except:
			print "TRIED: " + url+'?%s' % params
			print "EXCEPT: sdss.py._query()"
			return StringIO.StringIO() # This is an empty filehandler
