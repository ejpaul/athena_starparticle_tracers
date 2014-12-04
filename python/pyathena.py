import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import glob
import os
import re

def parse_filename(filename):
	"""
	#   PARSE_FILENAME    Break up a full-path filename into its component
	#   parts to check the extension, make it more readable, and extract the step
	#   number.  
	#
	#   PATH,BASENAME,STEP,EXT = PARSE_FILENAME(FILENAME)
	#
	#   E.g. If FILENAME='/home/Blast.0000.bin', then PATH='/home',
	#   BASENAME='Blast', STEP='0000', and EXT='bin'.
	#
	"""


	path=os.path.dirname(filename)
	if path[-3:] == 'id0': 
		path=path[:-3]
		mpi_mode=True
	else:
		path=path+os.path.sep
		mpi_mode=False

	base=os.path.basename(filename)
	id=base[:-9]
	step, ext=base[-8:].split('.')

	return path,id,step,ext,mpi_mode

def parse_line(line, grid):
	sp = line.strip().split()

	if "vtk" in sp:
		grid['vtk_version'] = sp[-1]
	elif "time=" in sp:
		time_index = sp.index("time=")
        	grid['time'] = float(sp[time_index+1].rstrip(','))
        	if 'level' in sp: grid['level'] = int(sp[time_index+3].rstrip(','))
        	if 'domain' in sp: grid['domain'] = int(sp[time_index+5].rstrip(','))  
		if sp[0] == "PRIMITIVE": 
			grid['prim_var_type']=True
	elif "DIMENSIONS" in sp:
		grid['Nx'] = np.array(sp[-3:]).astype('int')
	elif "ORIGIN" in sp:
		grid['left_edge'] = np.array(sp[-3:]).astype('float64')
	elif "SPACING" in sp:
		grid['dx'] = np.array(sp[-3:]).astype('float64')
	elif "CELL_DATA" in sp:
		grid['ncells'] = int(sp[-1])
	elif "SCALARS" in sp:
		grid['read_field'] = sp[1]
		grid['read_type'] = 'scalar'
	elif "VECTORS" in sp:
		grid['read_field'] = sp[1]
		grid['read_type'] = 'vector'
	elif "NSTARS" in sp:
		grid['nstar'] = eval(sp[1])


class AthenaDomain(object):
	def __init__(self,filename,grids=None):
		self.flist = glob.glob(filename)
		if len(self.flist) == 0: 
			print 'no such file: %s' % filename
		dir, id, step, ext, mpi = parse_filename(filename)
		self.dir = dir
		self.id = id
		self.step = step
		self.ext = ext
		self.starfile = os.path.join(dir+'id0/','%s.%s.%s.%s' % (id,step,'starpar',ext))
		if mpi: self.flist += glob.glob(os.path.join(dir,'id*/%s-id*.%s.%s' % (id, step, ext)))
		self.flist.sort()
		self.ngrids = len(self.flist)
		if grids==None: self.grids=self.setup_grid()
		else: 
			if grids[0]['filename'] != self.flist[0]:
				for g,f in zip(grids,self.flist): g['filename']=f
			self.grids=grids
		self.domain=self.setup_domain(self.grids)
		self.setup()
	
	def setup(self):
		self.domain['data']={}

	def setup_domain(self,grids):
		domain = {}
		left_edges = np.empty((self.ngrids,3), dtype='float64')
		dxs = np.empty((self.ngrids,3), dtype='float64')
		Nxs = np.ones_like(dxs)
		for nproc,g in enumerate(grids):
			left_edges[nproc,:] = g['left_edge']
			Nxs[nproc,:] = g['Nx']
			dxs[nproc,:] = g['dx']

		right_edges = left_edges + Nxs*dxs

		left_edge = left_edges.min(0)
		right_edge = right_edges.max(0)

		gis = np.round((left_edges - left_edge)/dxs)
		
		for nproc,g in enumerate(grids):
			g['is']=gis[nproc,:]

		domain['left_edge'] = left_edge
		domain['right_edge'] = right_edge
		domain['dx'] = dxs[0,:]
		domain['Lx'] = right_edge - left_edge
		domain['center'] = 0.5*(right_edge + left_edge)
		domain['Nx'] = np.round(domain['Lx']/domain['dx']).astype('int')
		domain['ndim'] = 3 # should be revised
		file = open(self.flist[0],'rb')
		tmpgrid = {}
		tmpgrid['time']=None
		while tmpgrid['time'] is None:
			line = file.readline()
			parse_line(line,tmpgrid)
		file.close()
		domain['time'] = tmpgrid['time']
		domain['data'] = {}
		domain['field_map']=None

		return domain


	def setup_grid(self):
		grids=[]
		for nproc in range(self.ngrids):
			file = open(self.flist[nproc],'rb')
			grid = {}
			grid['filename']=self.flist[nproc]
			grid['read_field'] = None
			grid['read_type'] = None
			while grid['read_field'] is None:
				grid['data_offset']=file.tell()
				line = file.readline()
				parse_line(line, grid)
			file.close()
			grid['Nx'] -= 1
			grid['Nx'][grid['Nx'] == 0] = 1
			grid['dx'][grid['Nx'] == 1] = 1.
			grid['right_edge'] = grid['left_edge'] + grid['Nx']*grid['dx']
			#grid['field_map']=None

			grids.append(grid)

		return grids

class AthenaDataSet(AthenaDomain):
	def setup(self):
		for i,g in enumerate(self.grids):
#			if g['field_map']==None: 
#				print "Setting %d-th grid" % (i)
#				self.set_field_map(g)
#			else: 
			g['data']={}
		sd=self.domain
		fm=sd['field_map']
		if fm==None:
			self.domain['field_map'] = self.set_field_map(self.grids[0])
		
		
	def set_field_map(self,grid):
		file=open(grid['filename'],'rb')
		file.seek(0,2)
		eof = file.tell()
		offset = grid['data_offset']
		file.seek(offset)

		field_map={}

		while offset < eof:

			line=file.readline()
			sp = line.strip().split()
			field_map[sp[1]] = {}
			field_map[sp[1]]['read_table']=False

			if "SCALARS" in line:
				tmp=file.readline()
				field_map[sp[1]]['read_table']=True
				field_map[sp[1]]['nvar'] = 1
			elif "VECTORS" in line:
				field_map[sp[1]]['nvar'] = 3
			else:
				print 'Error: '+sp[0] + ' is unknown type'
				raise TypeError

			field_map[sp[1]]['offset']=offset
			field_map[sp[1]]['ndata']=field_map[sp[1]]['nvar']*grid['ncells']
			if sp[2]=='int': dtype='i'
			elif sp[2]=='float': dtype='f'
			field_map[sp[1]]['dtype']=dtype
			field_map[sp[1]]['dsize']=field_map[sp[1]]['ndata']*struct.calcsize(dtype)
			file.seek(field_map[sp[1]]['dsize'],1)
			file.readline()
			offset = file.tell()

		#grid['field_map'] = field_map
		#grid['data']={}
		return field_map


	def _read_field(self,file_pointer,field_map):
		ndata=field_map['ndata']
		dtype=field_map['dtype']
		file_pointer.seek(field_map['offset'])
		file_pointer.readline() # HEADER
		if field_map['read_table']: file_pointer.readline()
		data = file_pointer.read(field_map['dsize'])
		var = np.asarray(struct.unpack('>'+ndata*dtype,data))

		return var

	def _read_grid_data(self,grid,field):
		file=open(grid['filename'],'rb')
		#fm=grid['field_map']
		fm=self.domain['field_map']
		nx1=grid['Nx'][0]
		nx2=grid['Nx'][1]
		nx3=grid['Nx'][2]

		nvar=fm[field]['nvar']
		var = self._read_field(file,fm[field])
		if nvar == 1: 
			var.shape = (nx3, nx2, nx1)
			grid['data'][field]=var
		else: 
			var.shape = (nx3, nx2, nx1, nvar)
			for i in range(nvar):
				var3d=var[:,:,:,i]
				grid['data'][field+"%s" % (i+1)]=var3d
		file.close()


	def set_grid_data(self,grid,field):
		gd=grid['data']
		if gd.has_key(field):
 			return gd[field]

		#fm = grid['field_map']
		fm = self.domain['field_map']
		if field in fm.keys():
			nvar=fm[field]['nvar']
			self._read_grid_data(grid,field)
		elif field == 'magnetic_pressure' or field == 'plasma_beta':
			self._read_grid_data(grid,'cell_centered_B')
			gd['magnetic_pressure']=0.5*(gd['cell_centered_B1']**2+gd['cell_centered_B2']**2+gd['cell_centered_B3']**2)
			if not gd.has_key('pressure'): self._read_grid_data(grid,'pressure')
			gd['plasma_beta']=gd['pressure']/gd['magnetic_pressure']
			nvar=1
			
		if nvar == 1: return gd[field]
		elif nvar == 3: return gd[field+'1'],gd[field+'2'],gd[field+'3']
		
	def read_all_data(self,field):
		#fm=self.grids[0]['field_map']
		fm=self.domain['field_map']
		nvar=fm[field]['nvar']
		dnx=self.domain['Nx']
		if nvar==1:
			data=np.empty((dnx[2],dnx[1],dnx[0]),dtype='float64')
		elif nvar==3:
			data1=np.empty((dnx[2],dnx[1],dnx[0]),dtype='float64')
			data2=np.empty((dnx[2],dnx[1],dnx[0]),dtype='float64')
			data3=np.empty((dnx[2],dnx[1],dnx[0]),dtype='float64')

		for g in self.grids:
			gis=g['is']
			gnx=g['Nx']
			gie=gis+gnx
			if nvar==1:
				gd=self.set_grid_data(g,field)
				data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd
			elif nvar==3:
				gd1,gd2,gd3=self.set_grid_data(g,field)
				data1[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd1
				data2[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd2
				data3[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd3
			
		if nvar==1:
			return data
		elif nvar==3:
			return data1,data2,data3

	def read_starvtk(self):
		file=open(self.starfile,'rb')
		star = {}
		star['filename']=self.starfile
		star['read_field'] = None
		star['read_type'] = None
		while star['read_field'] is None:
			star['data_offset']=file.tell()
			line = file.readline()
			parse_line(line, star)

		nstar=star['nstar']
		fm=self.set_field_map(star)
		id=self._read_field(file,fm['star_particle_id'])
		mass=self._read_field(file,fm['star_particle_mass'])
		pos=self._read_field(file,fm['star_particle_position']).reshape(nstar,3)
		vel=self._read_field(file,fm['star_particle_velocity']).reshape(nstar,3)
		file.close()
		star=[]
		for i in id:
			star.append({})

		for i in range(nstar):
			star_dict = star[id[i]]
			star_dict['id']=id[i]
			star_dict['mass']=mass[i]
			star_dict['velocity']=vel[i]
			star_dict['position']=pos[i]

		return star

class AthenaRegion(object):
	def __init__(self,ds,*args,**kwargs):
		self.ds=ds
		self._region_selector(*args,**kwargs)

	def _plane_grid_selector(self):

		grid_list=[]
		for g in self.ds.grids:
			axis=self.region['axis']
			c=self.region['center'][axis]
			le=self.region['left_edge']
			re=self.region['right_edge']
			aidx=self.region['axis_idx']

			zmin=g['left_edge'][axis]
			zmax=g['right_edge'][axis]
			gle=g['left_edge'][aidx]
			gre=g['right_edge'][aidx]
			
			if (zmin-c)*(zmax-c)<=0:
				dle = gle-le
				dre = gre-re
				if dle[0]*dre[0]<=0 and dle[1]*dre[1]<=0:
					#print g['left_edge'],g['right_edge'],c
					grid_list.append(g)

		print len(grid_list)
		self.region['grid_list']=grid_list


class AthenaSlice(AthenaRegion):
	def _region_selector(self,*args,**kwargs):
		self.slice(*args,**kwargs)

	def slice(self,*args,**kwargs):
		if kwargs.has_key('axis'): axis=kwargs['axis']
		else: print 'need to specify symmetric axis'

		if kwargs.has_key('center'): c=kwargs['center']
		else: print 'need to specify center of the plane'

		if kwargs.has_key('field'): field=kwargs['field']
		else: print 'need to specify field'

		if axis == 'x': axis = 0
		elif axis == 'y': axis = 1
		elif axis == 'z': axis = 2

		if c=='center': c=self.ds.domain['center']

		aidx=[0,1,2]
		aidx.remove(axis)
		self.axis_labels=['x','y','z']
		self.axis_labels.remove(self.axis_labels[axis])
		le=self.ds.domain['left_edge'][aidx]
		re=self.ds.domain['right_edge'][aidx]
		nx=self.ds.domain['Nx'][aidx]

		self.region={'axis':axis,'center':c,'left_edge':le,'right_edge':re,'axis_idx':aidx}
		self.data=np.empty((nx[1],nx[0]),dtype='float64')
		self.bound=(le[0],re[0],le[1],re[1])

		self._plane_grid_selector()

		#ng=0
		for g in self.region['grid_list']:
			#print "Reading:",g['filename']
			gd=g['data']
			cidx=cc_idx(self.ds.domain,c)-g['is']

			gis=g['is'][aidx]
			gnx=g['Nx'][aidx]
			gie=gis+gnx
			if not gd.has_key(field): 
				#ng = ng+1
				#print "Reading:",g['filename'],ng
				self.ds.set_grid_data(g,field)
			if axis==0:
				data=gd[field][:,:,cidx[0]]
			if axis==1:
				data=gd[field][:,cidx[1],:]
			if axis==2:
				data=gd[field][cidx[2],:,:]
			self.data[gis[1]:gie[1],gis[0]:gie[0]]=data


def cc_pos(domain,idx):
	le=domain['left_edge']
	dx=domain['dx']
	return le+0.5*dx+dx*np.array(idx)

def cc_idx(domain,pos):
	le=domain['left_edge']
	dx=domain['dx']
	return (pos-le-0.5*dx)/dx
