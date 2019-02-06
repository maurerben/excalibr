#!/usr/bin/python

from lxml import etree
import sys
import numpy as np
from math import *
import os
import matplotlib.pyplot as plt

class bandstr:
    """
    Object to store the data contained in the bandstructure.xml file
    of an exciting calculation.
    Args:
    name (string): filename of the bandstructure xml information
    .. attribute :: name
        associated filename
    """
    def __init__(self,name):
        if os.path.isfile(name):
            self._setparameter(name)
        else:
            print("\nScript needs xml file\"",name,"\"to run\n")
            sys.exit()
    def _setparameter(self,name):
        self.tree = etree.parse(name)
        self.filename=name
        try:
            if  self.tree.xpath('/bandstructure/@character')[0]=='true':
                self.character=True
            else:
                self.character=False
        except IndexError:
            self.character=False
        if self.character:
            self.pts=self.tree.xpath('/bandstructure/species[1]/atom[1]/band[1]/point/@distance')
            self.bands=len(self.tree.xpath('/bandstructure/species[1]/atom[1]/band'))
        else:
            self.pts = self.tree.xpath('/bandstructure/band[1]/point/@distance')
            self.bands =len(self.tree.xpath('/bandstructure/band'))

    def _getband(self, bandnr):
        if not self.character:
            y = [float(xe)*27.211 for xe in self.tree.xpath("/bandstructure/band[%i]/point/@eval"%(bandnr))]
        else:
            y = [float(xe)*27.211 for xe in self.tree.xpath("/bandstructure/species[1]/atom[1]/band[%i]/point/@eval"%(bandnr))]

        return self.pts, y

    def getbands(self, bandmin, bandmax):
        """
        returns the energy bands for the path in reciprocal space when the band character of the
        bands is not calculated. As often high-energy or core states are not of interest,
        the function allows to specify a window for the bands.
        Args:
            param1 (int):: bandmin
                lowest band to be included
            param2 (int):: bandmax
                highest band to be included
        Returns:
            list :: list of points along the high-symmetry path
            list :: list of bands of size (bandmax-bandmin+1), where every item in the list
                is a list of the energy values for every entry in x
            float :: lowest energy value for lowest band in output
            float :: highest energy value for highest band in output
        Raises:
            ValueError: Subroutine can not be used when the character of the bands is
            included in the class. 
        """
        
        y=[]
        for i in range(bandmin,bandmax):
            if (i==bandmin):
                x=self._getband(i)[0]
            y.append(self._getband(i)[1])
            
        windowlow=min(y[0])
        windowhigh=max(y[-1])
        return x,y, windowlow, windowhigh

    def getlabels(self):
        """
        Function to get the labels of high-symmetry points along the path in reciprocal
        space as defined in bandstructure.xml
        Returns:
            float:: list of the position for the labels along the x-axis
            str :: list that contains the label for each point along the axis. For the GAMMA point,
                the label is written as $\Gamma to allow the direct usage in matplotlib
        """
        nrlabel=len(self.tree.xpath('/bandstructure/vertex'))
        x=[]
        label=[]
        for i in range(0,nrlabel):
            x.append(float(self.tree.xpath('/bandstructure/vertex[%i]/@distance'%(i+1))[0]))
            label.append(str(self.tree.xpath('/bandstructure/vertex[%i]/@label'%(i+1))[0]))
            if (label[i]=='GAMMA'): 
                label[i]='$\Gamma$'
            else:
                label[i]='$'+label[i]+'$'
        return x,label

    def _get_character_bands(self):
        if self.character:
            out=[]
            bands=self.tree.xpath('/bandstructure/species[1]/atom[1]/band')
            for band in bands:
                for xe in band.xpath("./point/@eval"):
                    out.append(float(xe)*27.211)
        else:
            raise ValueError('Subroutine not defined if character is false!')
        return self.pts, out

    def _get_bandcharacter(self,species,atom,l):
        if self.character:
            out_bands=[]
            bands=self.tree.xpath('/bandstructure/species[%s]/atom[%s]/band'%(species, atom)) 
            for band in bands:
                for pts in band.xpath('./point'):
                    out_bands.append(pts.xpath('./bc[%s]/@character'%(int(l)+1))[0])
                    
        else:
            raise ValueError('Subroutine not defined if character is false!')

        return out_bands

    def atoml(self,species, atom, l):
        """
        Function to get the bandstructure and the character of the bands, where the
        projection is on a specific atom and l-channel. Note that the output of this function is 
        different from that of getbands, to be consistent with 
        matplotlib.pyplot.scatter.
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: atom
                atom number of the species as defined in the exciting input file. Note 
                that the first atom of each species has atomnumber=1
            int :: l
                l-channel to be projected on. Only projection on l=0, l=1, l=2, or 
                l=3 has been implemented in exciting yet.

        Returns:
            list :: energy points for each band in one list of size 
                len(self.pts)*self.bands
            list :: energy value for each entry in the first output
            list :: character for each entry in the first output
        """
       
        pts, bands=self._get_character_bands()
        characters=self._get_bandcharacter(species,atom, l)
        pts_out=[]
        for i in range(self.bands):
            for j in range(len(self.pts)):
                pts_out.append(pts[j])
        return pts_out, bands, characters
    
    def atom(self, species, atom):
        """
        Function to get the bandstructure and the character of the bands, where the
        projection is on a specific atom. Note that the output of this function is 
        different from that of getbands, to be consistent with 
        matplotlib.pyplot.scatter.
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: atom
                atom number of the species as defined in the exciting input file. Note 
                that the first atom of each species has atomnumber=1
        Returns:
            list :: energy points for each band in one list of size 
                len(self.pts)*self.bands
            list :: energy value for each entry in the first output
            list :: character for each entry in the first output
        """
        pts,bands=self._get_character_bands()
        out_bands=[]

        pts_out=[]
        for i in range(self.bands):
            for j in range(len(self.pts)):
                pts_out.append(pts[j])       
        __bands=self.tree.xpath('/bandstructure/species[%s]/atom[%s]/band'%(species, atom))
        for band in __bands:
            for pts in band.xpath('./point'):
                out_bands.append(float(pts.xpath('./@sum')[0]))
        return pts_out, bands, out_bands  

    def speciesl(self, species, l):
        """
        Function to get the bandstructure and the character of the bands, where the
        projection is on a specific species and l-channel. Note that the output of this function is 
        different from that of getbands, to be consistent with 
        matplotlib.pyplot.scatter.
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: l
                l-channel to be projected on. Only projection on l=0, l=1, l=2, or 
                l=3 has been implemented in exciting yet.
        Returns:
            list :: energy points for each band in one list of size 
                len(self.pts)*self.bands
            list :: energy value for each entry in the first output
            list :: character for each entry in the first output
        """

        pts, bands=self._get_character_bands()

        pts_out=[]
        for i in range(self.bands):
            for j in range(len(self.pts)):
                pts_out.append(pts[j])

        atoms=self.tree.xpath('/bandstructure/species[%s]'%species)
        i=0
        for atom in atoms:
            if i==0:
                character=np.array(self.atoml(species,i+1,l)[2]).astype(np.float)
                i+=1
            else:
                character=character+np.array(self.atoml(species,i+1,l)[2]).astype(np.float)
                i+=1

        return pts_out, bands, character

    def species(self, species):
        """
        Function to get the bandstructure and the character of the bands, where the
        projection is on a specific species. Note that the output of this function is 
        different from that of getbands, to be consistent with 
        matplotlib.pyplot.scatter.
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
        Returns:
            list :: energy points for each band in one list of size 
                len(self.pts)*self.bands
            list :: energy value for each entry in the first output
            list :: character for each entry in the first output
        """
       
        pts, bands=self._get_character_bands()

        pts_out=[]
        for i in range(self.bands):
            for j in range(len(self.pts)):
                pts_out.append(pts[j])

        atoms=self.tree.xpath('/bandstructure/species[%s]/atom'%species)
        character = np.array([])
        i = 0
        for atom in atoms:
            if i==0:
                character=np.array(self.atom(species,i+1)[2]).astype(np.float)
                i+=1
            else:
                character=character+np.array(self.atom(species,i+1)[2]).astype(np.float)
                i+=1

        return pts_out, bands, character.tolist()

    """ 
    Functions for plotting bandstructures and visulizing band characters:
    """

    def plot_band(self,**kwargs):
        """ 
        Creates a matplotlib.pyplot.plot object wich can be showed by plt.show().
        Args (optional):
            string :: label
                defines label in the plot
            string :: color
                defines the color of the bands
            list[float] :: ylim
                defines the limits of the y axis in the plot
            float :: offset
                allows to shift the bands by an offset.
        """
        label = kwargs.pop("label","")
        color = kwargs.pop("color","tab:blue")
        ylim  = kwargs.pop("ylim",[-10,10])
        offset = kwargs.pop("offset",0)

        x0,y0,miny, maxy = self.getbands(1,self.bands+1)
        labelx, labely = self.getlabels()
        
        y=y0
        if (offset != 0):
            y = [[val - offset for val in band] for band in y0]
        
        x = [float(x) for x in x0]

        plt.ylabel("Energy [eV]")
        plt.ylim(ylim)
        plt.xticks(labelx,labely)
        for i in range(0,len(y)):  
            if (i==0):
                plt.plot(x,y[i], zorder = 1, color=color,label=label)
            else :
                plt.plot(x,y[i],  zorder = 1, color=color)
        plt.tick_params(axis='both', which='major', labelsize=20)
        for i in range(0,len(labely)):
            plt.plot((labelx[i],labelx[i]),(ylim[0],ylim[1]),'-k')
        plt.plot((labelx[0],labelx[-1]),(0,0),'--k')

    def scatter_atoml(self,species, atom, l,**kwargs):
        """ 
        Creates matplotlib.pyplot.scatter object wich can be showed by plt.show().
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: atom
                atom number of the species as defined in the exciting input file. Note 
                that the first atom of each species has atomnumber=1
            int :: l
                l-channel to be projected on. Only projection on l=0, l=1, l=2, or 
                l=3 has been implemented in exciting yet.
        Args (optional):
            string :: label
                defines label in the plot. 
                Default: empty
            string :: color
                defines the color of the bands
                Default: tab:blue
            list[float] :: ylim
                defines the limits of the y axis in the plot
                Default: [-10,10]
            float :: offset
                allows to shift the bands by an offset.
                Default: 0
            float :: scale
                multiplies with the character to make it visible.
                Default: 1
        """
        label = kwargs.pop("label","")
        color = kwargs.pop("color","tab:blue")
        ylim  = kwargs.pop("ylim",[-10,10])
        offset = kwargs.pop("offset",0)
        scale = kwargs.pop("scale", 1)

        x,y,char = self.atoml(species, atom, l)
        labelx, labely = self.getlabels()

        x = [float(pt) for pt in x]
        y = [e-offset  for e in y]
        char = [float(ch) for ch in char]
        char = np.array(char) * scale 

        plt.scatter(x, y, s=char, marker='o', facecolors='none', zorder = 2, edgecolors=color, label=label)
        plt.xticks(labelx,labely)
        for i in range(0,len(labely)):
            plt.plot((labelx[i],labelx[i]),(ylim[0],ylim[1]),'-k')
        plt.plot((labelx[0],labelx[-1]),(0,0),'--k')
    
    def scatter_atom(self, species, atom, **kwargs):
        """ 
        Creates matplotlib.pyplot.scatter object wich can be showed by plt.show().
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: atom
                atom number of the species as defined in the exciting input file. Note 
                that the first atom of each species has atomnumber=1
        Args (optional):
            string :: label
                defines label in the plot. 
                Default: empty
            string :: color
                defines the color of the bands
                Default: tab:blue
            list[float] :: ylim
                defines the limits of the y axis in the plot
                Default: [-10,10]
            float :: offset
                allows to shift the bands by an offset.
                Default: 0
            float :: scale
                multiplies with the character to make it visible.
                Default: 1
        """
        label = kwargs.pop("label","")
        color = kwargs.pop("color","tab:blue")
        ylim  = kwargs.pop("ylim",[-10,10])
        offset = kwargs.pop("offset",0)
        scale = kwargs.pop("scale", 1)

        x,y,char = self.atom(species,atom)
        labelx, labely = self.getlabels()

        x = [float(pt) for pt in x]
        y = [e-offset  for e in y]
        char = np.array(char) * scale

        plt.scatter(x, y, s=char, marker='o', facecolors='none', zorder = 2, edgecolors=color, label=label)
        plt.xticks(labelx,labely)
        for i in range(0,len(labely)):
            plt.plot((labelx[i],labelx[i]),(ylim[0],ylim[1]),'-k')
        plt.plot((labelx[0],labelx[-1]),(0,0),'--k')

    def scatter_speciesl(self, species, l, **kwargs):
        """ 
        Creates matplotlib.pyplot.scatter object wich can be showed by plt.show().
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
            int :: l
                l-channel to be projected on. Only projection on l=0, l=1, l=2, or 
                l=3 has been implemented in exciting yet.
        Args (optional):
            string :: label
                defines label in the plot. 
                Default: empty
            string :: color
                defines the color of the bands
                Default: tab:blue
            list[float] :: ylim
                defines the limits of the y axis in the plot
                Default: [-10,10]
            float :: offset
                allows to shift the bands by an offset.
                Default: 0
            float :: scale
                multiplies with the character to make it visible.
                Default: 1
        """
        label = kwargs.pop("label","")
        color = kwargs.pop("color","tab:blue")
        ylim  = kwargs.pop("ylim",[-10,10])
        offset = kwargs.pop("offset",0)
        scale = kwargs.pop("scale", 1)

        x,y,char = self.speciesl(species,l)
        labelx, labely = self.getlabels()

        x = [float(pt) for pt in x]
        y = [e-offset  for e in y]
        char = np.array(char) * scale

        

        plt.scatter(x, y, s=char, marker='o', facecolors='none', zorder = 2, edgecolors=color, label=label)
        plt.xticks(labelx,labely)
        for i in range(0,len(labely)):
            plt.plot((labelx[i],labelx[i]),(ylim[0],ylim[1]),'-k')
        plt.plot((labelx[0],labelx[-1]),(0,0),'--k')

    def scatter_species(self,species,**kwargs):
        """ 
        Creates matplotlib.pyplot.scatter object wich can be showed by plt.show().
        Args:
            int :: species
                species number as defined in the exciting input file. Note that in this
                definition, the first species has speciesnumber=1
        Args (optional):
            string :: label
                defines label in the plot. 
                Default: empty
            string :: color
                defines the color of the bands
                Default: tab:blue
            list[float] :: ylim
                defines the limits of the y axis in the plot
                Default: [-10,10]
            float :: offset
                allows to shift the bands by an offset.
                Default: 0
            float :: scale
                multiplies with the character to make it visible.
                Default: 1
        """
        label = kwargs.pop("label","")
        color = kwargs.pop("color","tab:blue")
        ylim  = kwargs.pop("ylim",[-10,10])
        offset = kwargs.pop("offset",0)
        scale = kwargs.pop("scale", 1)

        x,y,char = self.species(species)
        labelx, labely = self.getlabels()

        x = [float(pt) for pt in x]
        y = [e-offset  for e in y]
        char = np.array(char) * scale
        
        plt.scatter(x, y, s=char, marker='o', facecolors='none', zorder = 2, edgecolors=color, label=label)
        plt.xticks(labelx,labely)
        for i in range(0,len(labely)):
            plt.plot((labelx[i],labelx[i]),(ylim[0],ylim[1]),'-k')
        plt.plot((labelx[0],labelx[-1]),(0,0),'--k')


        






    


