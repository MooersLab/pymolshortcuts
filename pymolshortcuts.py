from __future__ import division
# -*- coding: utf-8 -*-
"""
    DESCRIPTION

        On the startup of PyMOL, this script defines a number of aliases.
        The aliases are listed here instead of in the pymolrc file
        to avoid clutter of the command history window. Source this 
        file from your .pymolrc file on the mac or linux or 
        from your pymolrc.pml file on Windows by adding the command:
                           
        run ~/mg18OU/startupAliases.py  
        
        Requires quat.py from the PyMOL Wiki 
        (http://www.pymolwiki.org/index.php/BiologicalUnit/Quat) 
        
        Store quat.py in ~/mg18OU/.
        
        Tested on PyMOL versions 1.5.0.5, 1.8.0.5, 1.8.1.0, 1.8.2.0, 1.8.2.2, 2.1.0
        
        Alternately, launch PyMOL from commandline with alias specified. For example:
        
        pymol -d 'fetch 1lw9, async=0; AO'

        No guarantee is given that this script will work with older
        or newer versions of PyMOL.


  Copyright Notice
  ================
  
     Copyright (C) 2016  Blaine Mooers

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details:
    http://www.gnu.org/licenses/.

  The source code in this file is copyrighted, but you can
  freely use and copy it as long as you don't change or remove any of
  the copyright notices.
  
  Blaine Mooers , PhD 
  blaine@ouhsc.edu
  975 NE 10th St, BRC 466
  University of Oklahoma Health Sciences Center, 
  Oklahoma City, OK, USA 73104

""" 
import webbrowser, datetime, numpy, subprocess, bs4, requests, time
from pymol import cmd, stored, cgo, xray
from math import cos, sin, radians, sqrt
# from subprocess import *


__author__ = "Blaine Mooers"
__copyright__ = "Blaine Mooers, University of Oklahoma Health Sciences Center, Oklahoma City, OK, USA 73104"
__license__ = "GPL-3"
__version__ = "0.1"
__credits__ = [""] 
# people who reported bug fixes, made suggestions, etc. 
__date__ = "22 June 2018"
__maintainer__ = "Blaine Mooers"
__email__ = "blaine@ouhsc.edu"
__status__ = "Testing" 

cmd.set('ray_opaque_background','on')

#################################################################################

#category: Print shortcuts and their descriptions

def SC():
    '''
    DESCRIPTION
    
    Print to screen list of the shortcuts that are available in the script pymolshortcuts.py.  
    Copyright (C) 2018  Blaine Mooers. 
    This script comes with ABSOLUTELY NO WARRANTY; for details, 
    please see the license file. 

    Copy to ~/pymolshortcuts.py where ~ is the path to your home directory.
    Add the following run command on one line in your .pymolrc (hidden) or pymolrc file 
    (pymolrc.pml on Windows) in your home directory to enable loading of the functions in
    pymolshortcuts.py everytime that you start PyMOL:
    
    run ~/pymolshortcuts.py

    Format of list of shortcuts below:
    shortcut name, description of shortcut: PDB code (where applicable).


    Information outside of PyMOL (internet connection need and some customization is required)

    ACA, American Crystallographic Association
    AX, arXiv search
    ALS, Advanced Light Source
    APS, Advanced Photon Source
    BC, BIOCAT BioSAXS beamline at APS
    BD, Biological SAXS Database
    BX, bioarvix search
    CH, UCSF Chimera
    CHESS, CHESS
    EMDB, Electron Microscopy Data Bank
    EP, EasyPyMOL github site
    GM, gmail
    GO, google search
    GS, google scholar
    GSM, google scholar multiple search terms
    JM, JMol
    LBSF, Laboratory of Molecular Biology at OUHSC
    MCL, link to Macromolecular Crystallography Laboratory at OU
    MG, links for OUHSC molecular graphics course 
    NDB, Nucleic Acid Database 
    NSLSII, 
    notPyMOL, list of other molecular graphics programs
    PM, pubmed search
    PML, PyMOL Users Mail Service
    PS, Protein Soceity
    PW, PyMOL wiki
    PPC, Protien Production Facility at OU
    PDB, Protein Data Bank 
    RS, RNA Society
    SAXS, SAXS links at OUHSC
    SSRL, SSRL Structural Molecular Biology
    SB, BioSAXS at SSRL
    SSURF, Soceity for Science at User Research Facilities
    WM, webmail
    WS, National Weather Service (time to head for the storm shelter?)


    Molecules in standard orientations (need internet connection): 

    GGT, gamma glutamyl transpeptidase as cartoon (1.67 ang): 4gdx.
    GU, 10-mer RNA with eight GU base pairs (1.32 ang): 4pco.
    N9, neuraminidase as cartoon, biological unit (1.55 ang): 4dgr.
    T4L, WT T4 lysozyme (1.09 ang) as a ribbon diagram: 3fa0. 
    U8, 16-mer dsRNA with 8 contiguous Us. U-helix RaNA (1.37 ang): 3nd3.
    WC8, 16-mer RNA with all Watson-Crick base pairs (1.52 ang): 3nd4.


    Complex figures to serve as templates: 
    
    BST, Base-stacking figure, (1.32 ang): 4pco. 
    LG, Electron density map of nine sugar glycan (1.55 ang): 4dgr. 
    NA, Sodium cation in major groove of 16-mer RNA (1.52 ang): 3nd4.

    Shortcuts defined with L append in the front of the above shortcuts 
    load files from the directory ~/mgOU18. (e.g., LT4L for T4L). There are
    for leading workshops where the internet connection is not reliable.

    LBST, Base-stacking figure, (1.32 ang): 4pco. 
    LGGT, gamma glutamyl transpeptidase as cartoon (1.67 ang): 4gdx.
    LGU, 10-mer RNA with eight GU base pairs (1.32 ang): 4pco.
    LLG, Electron density map of nine sugar glycan (1.55 ang): 4dgr. 
    LNA, Sodium cation in major groove of 16-mer RNA (1.52 ang): 3nd4.
    LN9, neuraminidase as cartoon, biological unit (1.55 ang): 4dgr.
    LT4L, WT T4 lysozyme (1.09 ang) as a ribbon diagram: 3fa0. 
    LU8, 16-mer dsRNA with 8 contiguous Us. U-helix RaNA (1.37 ang): 3nd3.
    LWC8, 16-mer RNA with all Watson-Crick base pairs (1.52 ang): 3nd4.

    Molecular representations that can be applied to any visible molecular object:
    
    AO, Make ambient occlusion image. Requires global view of protein.
    AOD, Make ambient occlusion image with dark carbons.
    BU, Display biological unit. 
    BW, Make black and white ribbon cartoon on white background.
    CB, Define color blind compatible coloring scheme. 
    CBSS, Color ribbon and cartoons with colorblind friendly colors. 
    CR, Make colored filled-ring cartoon of nucleic acids..
    CSS, Color ribbon and cartoons by secondary structure: red, green and yellow. 
    FR, Make filled-ring cartoon of nucleic acids.
    SE, Make SAXS envelope from a bead model.


    File management and scripting functions:

    dsc, delete supercell and symmetry mates
    gitinit, print the steps for making a git repository
    ms2pdf, make 3D pdfs from current scene of cartoons or molecular surfaces (requires IDTFconverter and LaTeX)
    nmr, show all of the models in a nmr structure
    nmroff, hide all but the average model in a nmr structure
    rv, return the viewport settings on one line
    sc2, generate supercell and symmetry mates with 2 cells in each direction
    sc3, generate supercell and symmetry mates with 3 cells in each direction
    spng, save current scene to png file with time stamp in filename
    spse, save current scene to session file with time stamp in filename
    timcolor, color by atomic models by biophysical properties with Tim Mather's color scheme

    Type 'help <ShortCutName>' (e.g., help spng) for a description of the
    function and for two sets of commands. The first set of commands has
    line breaks for easy selection single lines of commands. The second
    set of commands is one one line for easy copying and pasting of the
    entire horizontal script. The commands can be copied from the
    command history window and pasted onto the command line for code
    reuse. Some aliases require additional scripts. 
    
    Type 'SC' to refresh the list of aliases.
    Type 'help rline' to see commands for moving cursor on the command line.

    '''
    print(SC.__doc__)
cmd.extend('SC',SC)


############################### Show many models (NMR and crystal packing) ##############

#category: Show many models (NMR and crystal packing)

def nmr():
    """ 
    Description
    
    Show all of the models in nmr structure. 
    I can never remmber this command.
    """
    cmd.do('set all_states, on')
cmd.extend("nmr", nmr)


def nmroff():
    """ 
    Description
    
    Hide all but first model in a nmr structure. 
    """
    cmd.do('set all_states, off')
cmd.extend("nmr", nmr)


def rmsc():
    """
    Description
    
    Remove supercell and the symmetry mates.
    """
    cmd.do('delete supercell*;delete m*_*')
cmd.extend("rmsc", rmsc)


def sc111():
    """
    Description
    
    Make a lattice of 1 x 1 x 1 unit cells. 
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.
    """
    cmd.do('run $HOME/mg18OU/supercell.py')
    cmd.do('supercell 1, 1, 1, , orange, supercell111, 1')
cmd.extend("sc111", sc111)


def sc221():
    """
    Description
    
    Make a lattice of 2 x 2 x 1 unit cells. 
    Use 'rmsc' to remove supercell objects. 
    Requires Thomas Holder's supercell.py script.
    """
    cmd.do('run $HOME/mg18OU/supercell.py')
    cmd.do('supercell 2, 2, 1, , orange, supercell221, 1')
cmd.extend("sc221", sc221)


def sc112():
    """
    Description
    
    Make a lattice of 1 x 1 x 2 unit cells. 
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.
    """
    cmd.do('run $HOME/mg18OU/supercell.py')
    cmd.do('supercell 1, 1, 2, , orange, supercell112, 1')
cmd.extend("sc221", sc112)


def sc222():
    """
    Description
    
    Make a lattice of 2 x 2 x 2 unit cells. 
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.
    """
    cmd.do('run $HOME/mg18OU/supercell.py')
    cmd.do('supercell 2, 2, 2, , orange, supercell222, 1')
cmd.extend("sc222", sc222)


def sc333():
    """
    Description
    
    Make a lattice of 3 x 3 x 3 unit cells. 
    Use 'rmsc' to remove supercell objects. 
    Requires Thomas Holder's supercell.py script.
    
    """
    cmd.do('run $HOME/mg18OU/supercell.py')
    cmd.do('supercell 3, 3, 3, , green, supercell333, 1')
cmd.extend("sc333", sc333)


#################### Save files with date and time in filename #########################spng

#category: Save files with date and time in filename

def saln(stemName="saved"):
    """
    DESCRIPTION

    Save a aln file (alignment file) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    saln currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".aln") 
cmd.extend('saln',saln)


def scif(stemName="saved"):
    """
    DESCRIPTION

    Save a cif file (Crystallographic Information File) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    scif currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".cif") 
cmd.extend('scif',scif)


def sccp4(stemName="saved"):
    """
    DESCRIPTION

    Save a ccp4 file (CCP4 electron density map file) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    sccp4 currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".ccp4") 
cmd.extend('sccp4',sccp4)


def sdae(stemName="saved"):
    """
    DESCRIPTION

    Save a dae file (Collada File) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    sdae currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".dae") 
cmd.extend('sdae',sdae)


def sdat(stemName="saved"):
    """
    DESCRIPTION

    Save dat file (output data file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".dat") 
cmd.extend('sdat',sdat)


def sfasta(stemName="saved"):
    """
    DESCRIPTION

    Save a fasta file (sequence file) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    sfasta currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".fasta") 
cmd.extend('sfasta',sfasta)


def sidtf(stemName="saved"):
    """
    DESCRIPTION

    Save a idtf file (Intermediate Data Text Format) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    sidtf currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".idtf") 
cmd.extend('sidtf',sidtf)


def smae(stemName="saved"):
    """
    DESCRIPTION

    Save mae file (Maestro file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smoe currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mae") 
cmd.extend('smae',smae)


def smmd(stemName="saved"):
    """
    DESCRIPTION

    Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smmd currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mmd") 
cmd.extend('smmd',smmd)


def smmod(stemName="saved"):
    """
    DESCRIPTION

    Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smmd currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mmod") 
cmd.extend('smmod',smmod)


def spmo(stemName="saved"):
    """
    DESCRIPTION

    Save pmo file (XYZ, binary format file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spmo currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pmo") 
cmd.extend('spmo',spmo)


def smoe(stemName="saved"):
    """
    DESCRIPTION

    Save moe file (Molecular Operating Environment) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smoe currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".moe") 
cmd.extend('smoe',smoe)


def smol(stemName="saved"):
    """
    DESCRIPTION

    Save mol file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mol") 
cmd.extend('smol',smol)


def smol2(stemName="saved"):
    """
    DESCRIPTION

    Save mol2 (Sybyl file format) file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol2 currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mol2") 
cmd.extend('smol2',smol2)


def smtl(stemName="saved"):
    """
    DESCRIPTION

    Save mtl (Wavefront Material file format) file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smtl currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".mtl") 
cmd.extend('smtl',smtl)


def sobj(stemName="saved"):
    """
    DESCRIPTION

    Save obj file (Wavefront mesh file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".obj") 
cmd.extend('sobj',sobj)


def sout(stemName="saved"):
    """
    DESCRIPTION

    Save out file (output data file) with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".out") 
cmd.extend('sout',sout)


def spdb(stemName="saved"):
    """
    DESCRIPTION

    Save pdb file with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spng currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pdb") 
cmd.extend('spdb',spdb)


def spkl(stemName="saved"):
    """
    DESCRIPTION

    Save a pkl file (Python pickle file) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spkl currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pkl") 
cmd.extend('spkl',spkl)


def spkla(stemName="saved"):
    """
    DESCRIPTION

    Save a pkla file (Python pickle file) with a time stamp included in the filename to avoid overwriting work..
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spkl currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pkla") 
cmd.extend('spkla',spkla)



def spng(stemName="saved"):
    """
    DESCRIPTION

    Save png file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spng currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".png") 
cmd.extend('spng',spng)


def spov(stemName="saved"):
    """
    DESCRIPTION

    Save pov (POV-ray tracing file format) file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spov currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pov") 
cmd.extend('spov',spov)


def spqr(stemName="saved"):
    """
    DESCRIPTION

    Save pqr file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spng currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pqr") 
cmd.extend('spqr',spqr)


def spse(stemName="saved"):
    """ DESCRIPTION

    Save session file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    spse currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".pse") 
cmd.extend('spse',spse)


def ssdf(stemName="saved"):
    """
    DESCRIPTION

    Save sdf file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    smol currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".sdf") 
cmd.extend('ssdf',ssdf)


def swrl(stemName="saved"):
    """
    DESCRIPTION

    Save wrl (VRML 2 file format) file with a time stamp included in the filename to avoid overwriting work.
    Read as a commandline argument, a string as the filename stem or 
    use the default filename stem "saved".

    USAGE:
    
    swrl currentScene

    """
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT) 
    cmd.save(stemName+s+".wrl") 
cmd.extend('swrl',swrl)

############ Make molecular representations that are not available in PyMOL#############

#category: Make molecular representations that are not available in PyMOL


# ##### Commands applicable to displayed molecules. ##########
def AO():
    '''
    DESCRIPTION
    
    Commands to make ambient occlusion image like those in Qutemole. 

  
    USAGE

    Type 'AO' to execute. Type 'help AO' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    From 
    
    September 16, 2016: modified colors of H and C atoms following
    recommendation at http://cupnet.net/ambient-occlusion-pymol/.
    Changed these two commands to address warnings from PyMOL:
    Changed 'set depth_cue, off' to 'set depth_cue,0.'  
    Changed 'set light_count,10' to 'set light_count,8';  
 
    
    The commands with linebreaks:
    
    set_color oxygen, [1.0,0.4,0.4];
    set_color nitrogen, [0.5,0.5,1.0];
    remove solvent;
    as spheres;
    util.cbaw;    
    bg white;
    set light_count,10;
    set spec_count,1;
    set shininess, 10;
    set specular,0.25;
    set ambient,0;
    set direct,0;
    set reflect,1.5;
    set ray_shadow_decay_factor, 0.1;
    set ray_shadow_decay_range, 2;
    set depth_cue, 0;
    ray
    
    The commands without linebreaks:
    
    set_color oxygen, [1.0,0.4,0.4];set_color nitrogen, [0.5,0.5,1.0];remove solvent;as spheres;util.cbaw;bg white;set light_count,8;set spec_count,1;set shininess, 10;set specular,0.25;set ambient,0;set direct,0;set reflect,1.5;set ray_shadow_decay_factor, 0.1;set ray_shadow_decay_range, 2;set depth_cue,0;color gray20, symbol c;ray 
    
    '''
    cmd.set_color('oxygen', '[1.0,0.4,0.4]')
    cmd.set_color('nitrogen', '[0.5,0.5,1.0]')
    cmd.remove('solvent')
    cmd.show_as('spheres')
    cmd.util.cbaw()
    cmd.bg_color('white')
    cmd.set('light_count', '8')
    cmd.set('spec_count', '1')
    cmd.set('shininess', '10')
    cmd.set('specular', '0.25')
    cmd.set('ambient', '0')
    cmd.set('direct', '0')
    cmd.set('reflect', '1.5')
    cmd.set('ray_shadow_decay_factor', '0.1')
    cmd.set('ray_shadow_decay_range', '2')
    cmd.set('depth_cue','0')
    cmd.set('ray_opaque_background','on')
    cmd.ray()
cmd.extend('AO',AO)


# ##### Commands applicable to displayed molecules. ##########
def AOD():
    '''
    DESCRIPTION
    
    Make ambient occlusion image of any with dark carbon atoms. 

  
    USAGE

    Type 'AO' to execute. Type 'help AO' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    
    September 16, 2016: 
    Changed these two commands to address warnings from PyMOL:
    Changed 'set depth_cue, off' to 'set depth_cue,0.'  
    Changed 'set light_count,10' to 'set light_count,8';  
 
    
    The commands with linebreaks:
    
    set_color oxygen, [1.0,0.4,0.4];
    set_color nitrogen, [0.5,0.5,1.0];
    remove solvent;
    as spheres;
    util.cbaw;    
    bg white;
    set light_count,8;
    set spec_count,1;
    set shininess, 10;
    set specular,0.25;
    set ambient,0;
    set direct,0;
    set reflect,1.5;
    set ray_shadow_decay_factor, 0.1;
    set ray_shadow_decay_range, 2;
    set depth_cue, 0;
    color gray20, symbol c
    color gray90, symbol h
    ray
    
    The commands without linebreaks:
    
    set_color oxygen, [1.0,0.4,0.4];set_color nitrogen, [0.5,0.5,1.0];remove solvent;as spheres;util.cbaw;bg white;set light_count,8;set spec_count,1;set shininess, 10;set specular,0.25;set ambient,0;set direct,0;set reflect,1.5;set ray_shadow_decay_factor, 0.1;set ray_shadow_decay_range, 2;set depth_cue,0;color gray20, symbol c;color gray70, symbol h;ray 
    
    '''
    cmd.set_color('oxygen', '[1.0,0.4,0.4]')
    cmd.set_color('nitrogen', '[0.5,0.5,1.0]')
    cmd.remove('solvent')
    cmd.show_as('spheres')
    cmd.util.cbaw()
    cmd.bg_color('white')
    cmd.set('light_count', '8')
    cmd.set('spec_count', '1')
    cmd.set('shininess', '10')
    cmd.set('specular', '0.25')
    cmd.set('ambient', '0')
    cmd.set('direct', '0')
    cmd.set('reflect', '1.5')
    cmd.set('ray_shadow_decay_factor', '0.1')
    cmd.set('ray_shadow_decay_range', '2')
    cmd.set('depth_cue','0')
    cmd.color('gray20', 'symbol c')
    cmd.color('gray90', 'symbol h')
    cmd.set('ray_opaque_background','on')
    cmd.ray()
cmd.extend('AOD',AOD)


def BW():
    '''
    DESCRIPTION
    
    Commands to make black-and white-ribbon cartoon on a white background.
    Good for avoiding color figure charges. Requires a pdb file. Can cause PyMOL
    to crash if applied to a cartoon representation. Best applied while object is 
    still shown as lines.
    
    
    USAGE
    
    Orient struture as desired. Then type 'BW' to execute the function. Type
    'help BW' to see this documentation printed to the command history window.
    Select from the command history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    show cartoon; 
    hide lines; 
    hide nonbonded; 
    # black and white cartoon;
    # note how the dcomment is on a separate line and not to the right of a command; 
    set ray_trace_mode, 2; 
    bg_color white; 
    set antialias, 2; 
    ray 1600,1600; 
    png test.png
    
    Commands without linebreaks:
    
    show cartoon; hide lines; hide nonbonded; set ray_trace_mode, 2; # black and white cartoon; bg_color white; set antialias, 2; ray 1600,1600; png test.png
 
    '''
    cmd.show_as("cartoon", "all"); 
    cmd.hide('lines'); 
    cmd.hide('nonbonded'); 
    # black and white cartoon; 
    cmd.set('ray_trace_mode', '2'); 
    cmd.bg_color('white'); 
    cmd.set('antialias', '2'); 
    cmd.set('ray_opaque_background','on')
    cmd.ray('600','600'); 
    cmd.png('test.png')
cmd.extend('BW',BW) 


def BU():
    '''
    DESCRIPTION
    
    Commands to make biological unit. Requires a pdb file. There are
    other ways of displaying the biological unit in PyMOL. Depends on
    the quat.py script by Thomas Holder.
    
    USAGE
    
    Type 'BU' to execute. Type 'help BU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    run ~/mg18OU/quat.py; 
    quat 

    The commands without linebreaks:
    
    run ~/mg18OU/quat.py; quat 

    '''
#    cmd.alias('aBU', 'run ~/mg18OU/quat.py; quat') 
    cmd.do('run $HOME/mg18OU/quat.py')
#    cmd.run('$HOME/mg18OU/quat.py')
    cmd.do('quat')
cmd.extend('BU',BU)



def CB():
    '''
    DESCRIPTION
    
    Loads Jared Sampson's script "colorblindfriendly.py" from the
    ~/Pymol-script-repo directory. The colorblind-friendly color
    names are printed to the command history window and are available
    for use like standard colors in PyMOL.

    Read the header of colorblindfriendly.py for more information. Although
    this script is also in my PyMOL plugin collection, it is out of sight
    and out of mind. The listing of the alias to this script is a
    reminder to use it. 
    
    USAGE

    Type 'CB' to execute. Type 'help CB' to see this documentation
    printed to the command history window. 
    
    Command:
    
    run ~/Pymol-script-repo/colorblindfriendly.py 
    
    '''
    cmd.run('run $HOME/mg18OU/colorBlindFriendly.py')
cmd.extend('CB',CB)


def CR():
    '''
    DESCRIPTION
    
    Commands to make colored filled-ring cartoon of nucleic acids. May
    need to 'hide everything' first. If asymmetric unit has one strand
    of a dsRNA or dsDNA, remember to apply the BU alias to display the
    second strand. ddd
    
    Adapted from KP Wu's blog post:
    https://kpwu.wordpress.com/2012/05/24/pymol-different-colors-of-
    nucleic-acid-rings/
    
    USAGE

    Type 'CR' to execute. Type 'help CR' to see this documentation
    printed to the command history window. Select from the command 
    history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    hide everything;
    bg_color white;
    cartoon oval;
    set cartoon_ring_mode,3;
    set cartoon_nucleic_acid_color, blue;
    select rna_A, resn A;
    select rna_C, resn C;
    select rna_G, resn G;
    select rna_U, resn U;
    color yellow, rna_A;
    color red, rna_C;
    color gray40, rna_G;
    color palecyan, rna_U;
    as cartoon 
    
    The commands without linebreaks:
    
    hide everything;bg_color white;cartoon oval;set cartoon_ring_mode,3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C, resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A;color red, rna_C;color gray40, rna_G;color palecyan, rna_U;as cartoon 
    
    '''
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.cartoon('oval')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_nucleic_acid_color', 'blue')
    cmd.select('rna_A', 'resn A')
    cmd.select('rna_C', 'resn C')
    cmd.select('rna_G', 'resn G')
    cmd.select('rna_U', 'resn U')
    cmd.color('yellow', 'rna_A')
    cmd.color('red', 'rna_C')
    cmd.color('gray40', 'rna_G')
    cmd.color('palecyan', 'rna_U')
    cmd.show_as('cartoon')
    cmd.disable('rna_U')    
cmd.extend('CR',CR)


def CSS():
    '''
    DESCRIPTION
    
    Commands to color ribbon or cartoon representations of proteins by
    secondary structures. 
    
    USAGE

    Type 'CSS' to activate execute. Type 'help CSS' to see this documentation
    printed to the command history window. Select from the command 
    history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:

    as cartoon;
    color red, ss H;
    color yellow,ss S;
    color green, ss L+;
    
    The commands without linebreaks:
    
    as cartoon; color red, ss H; color yellow,ss S; color green, ss L+;
    
    '''
    cmd.show_as('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')
cmd.extend('CSS',CSS)


def CBSS():
    '''
    DESCRIPTION
    
    Apply colorblind-friendly coloring to ribbon or cartoon representations.
    Depends on colorblindfriendly.py. 
    Script is assumed to be stored in $HOME/Pymol-script-repo/.
    Adjust the path as needed.
    Script can be attained via the PyMOL wiki or the plugin manager or the git repo.
    
    USAGE

    Type 'CBSS' to execute. Type 'help CBSS' to print the
    documentation in the command history window. Select from the command 
    history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    run ~/Pymol-script-repo/colorblindfriendly.py;
    as cartoon;
    color cb_red, ss H;
    color cb_yellow,ss S;
    color cb_green, ss L+; 
    
    The commands without linebreaks:
    
    run $HOME/Pymol-script-repo/colorBlindFriendly.py;as cartoon;color cb_red, ss H;color cb_yellow,ss S;color cb_green, ss L+; 
 
    '''
    cmd.run('$HOME/mg18OU/colorBlindFriendly.py')
    cmd.show_as('cartoon')
    cmd.color('cb_red', 'ss H')
    cmd.color('cb_yellow', 'ss S')
    cmd.color('cb_green', 'ss L+')
cmd.extend('CBSS',CBSS)


def DU():
    '''
    DESCRIPTION
    
    Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins. 
    Set 50% transparent cartoon so it can be combined with lines, sticks, and ball-and-sticks (try BS shortcut).
    
    USAGE

    Type 'DU' to execute. Type 'help DU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    Commands with linebreaks:
    
    cartoon dumbbell;
    set cartoon_dumbbell_width, 0.2;
    set cartoon_dumbbell_radius, 0.4;
    show cartoon; 
    
    Commands without linebreaks:
    
    cartoon dumbbell;set cartoon_dumbbell_width, 0.2;set cartoon_dumbbell_radius, 0.4;show cartoon; 

    '''
    cmd.cartoon('dumbbell')
    cmd.set('cartoon_dumbbell_width', '0.2')
    cmd.set('cartoon_dumbbell_radius', '0.4')
    cmd.show('cartoon')
cmd.extend('DU',DU)    
    
    
def FR():
    '''
    DESCRIPTION
    
    Make filled-ring cartoon of nucleic acids. May need to enter 'hide everything' first. 
    Adapted from the script on http://www-cryst.bioc.cam.ac.uk/members/zbyszek/figures_pymol. 
    
    USAGE

    Type 'FR' to execute. Type 'help FR' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the PyMOL gui.
    
    Commands with linebreaks:
    
    show sticks;
    set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;
    set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;
    as cartoon;
    
    Commands without linebreaks:
    
    show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon; 
 
    '''
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')
cmd.extend('FR',FR)    


def HH():
    '''
    DESCRIPTION
    
    Hide hydrogen atoms of currently visible molecular objects. 
    
    USAGE

    Type 'HH' to execute. Type 'help HH' to see this documentation
    printed to the command history window. 
    
    The command:
    
    hide everything, name H* 
    
    '''
    cmd.hide('everything', 'name H*') 
cmd.extend('HH',HH)    


def PU():
    '''
    DESCRIPTION
    
    Make putty cartoon of main chain of nucleic acids and proteins.
    The radius of the cartoon is inversely proportional to
    the B-factors. Set 50% transparent cartoon so it can be combined with
    lines, sticks, and ball-and-sticks (try BS alias).
    
    USAGE

    Type 'PU' to execute. Type 'help PU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    The commands with linebreaks: 
    
    cartoon putty;
    set cartoon_ladder_mode, 0;
    set cartoon_transparency,0.5;
    set cartoon_ring_finder, 0;
    show cartoon 
    
    The commands without linebreaks: 
    
    cartoon putty;set cartoon_ladder_mode, 0;set cartoon_transparency,0.5;set cartoon_ring_finder, 0;show cartoon 

    '''
    cmd.cartoon('putty')
    cmd.set('cartoon_ladder_mode', '0')
    cmd.set('cartoon_transparency', '0.5')
    cmd.set('cartoon_ring_finder', '0')
    cmd.show('cartoon')
cmd.extend('PU',PU)


def SE():
    '''
    DESCRIPTION
    
    Commands to make SAXS envelope from a bead model.
    
    USAGE

    Enter 'SE' to execute. Type 'help SE' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
    
    Commands with linebreaks:
    
    alter all,vdw=4.5;
    set solvent_radius, 4.3;
    show surface; 
    
    Commands without linebreaks:
    alter all,vdw=4.5;set solvent_radius, 4.3;show surface; 

    '''
    cmd.alter('all,vdw=4.5')
    cmd.set('solvent_radius', '4.3')
    cmd.show('surface')
cmd.extend('SE',SE)    


def getchem():
    """
    DESCRIPTION
    
    Create selections based on the biophysical properties of each residue.
    This is mainly used by timcolor but could be run alone. 
    Developed by Tim Mathers's graduate student Phillipe.
    
    """
    cmd.select ('calcium','resn ca or resn cal')
    cmd.select ('acid','resn asp or resn glu or resn cgu')
    cmd.select ('basic','resn arg or resn lys or resn his')
    cmd.select ('nonpolar','resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
    cmd.select ('polar','resn ser or resn thr or resn asn or resn gln or resn tyr')
    cmd.select ('cys','resn cys or resn cyx')
    cmd.select ('backbone','name ca or name n or name c or name o')
    cmd.select ('none')
cmd.extend ("getchem",getchem)


def timcolor(selection='all'):
    '''
    DESCRIPTION
    
    Use Tim Mather's coloring scheme applied to the selections defined in getchem(). 
    The scheme is as follows: 

          'acid'    :  'red' 
          'basic'   :  'blue' 
          'nonpolar':  'orange'
          'polar'   :  'green' 
          'cys'     :  'yellow'
    
    USAGE: 
    
    timcolor <selection>
    '''
    print selection
    code={'acid'    :  'red'    ,
          'basic'   :  'blue'   ,
          'nonpolar':  'orange' ,
          'polar'   :  'green'  ,
          'cys'     :  'yellow'}
    cmd.do ('getchem')
    cmd.select ('none')
    for elem in code:
        line='color '+code[elem]+','+elem+'&'+selection
        print line
        cmd.do (line)
    word='color white,backbone &'+selection
    print word
    cmd.do (word)                  #Used to be in code, but looks like
                                   #dictionnaries are accessed at random
    cmd.hide ('everything','resn HOH')
cmd.extend ('timcolor',timcolor)


##################### Launch a full-featured text editor from PyMOL ############################

#category: Launch a full-featured text editor from PyMOL

def atom(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with the text editor Atom from within PyMOL. 
    Adjust the path as needed for your system.

    USAGE
    atom test.pml

    Note that Atom is slow to start. 
    You might want to keep it open all of the time!
    Consider Sublime Text 3 or Visual Code Studio which startup are faster. 
    '''
    arg = ("/usr/local/bin/atom " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('atom',atom)


def bbedit(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with the text editor bbedit from within PyMOL. 
    Adjust the path as needed for your system.

    USAGE
    bbedit test.pml

    Note that Atom is slow to start. 
    You might want to keep it open all of the time!
    Consider Sublime Text 3 or Visual Code Studio which startup are faster. 
    '''
    arg = ("/usr/local/bin/bbedit " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('bbedit',bbedit)


def code(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with Visual Studio Code from within PyMOL. 
    Install the bioSyntax extension (free) from the Visual Studio Marketplace
    to get color syntax highlighting of pml files along with Fasta and other
    sequence files. 
    https://marketplace.visualstudio.com/items?itemName=reageyao.biosyntax
    

    USAGE
    
    code test.pml

    Note:
    
    See https://code.visualstudio.com/docs/editor/command-line
    '''
    arg = ("/usr/local/bin/code " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('code',code)


def emacs(fileName="testme.pml"):
    '''
    DESCRIPTION

    Open file with emacs from within PyMOL. 
    Adjust path to emacs on your computer as needed.

    USAGE
    nv test.pml
    
    Currently opening emacs in new shell and testme.pml in buffer 2. 
    In Normal mode, enter :.bnext to switch to second buffer.
    '''
    arg = ("/opt/local/bin/emacs " + fileName)
    arg2 = ('--file ' + fileName)
    subprocess.call(['open', '-W', '-a', 'iTerm.app', '/opt/local/bin/emacs', '--args', arg2])
    #Popen(shlex.split("""x-terminal-emulator -e 'nv -c "testme.pml"'"""), stdout=PIPE)
    #process.wait()
    ########subprocess.call(arg,shell=True)
    #handle = Popen(arg, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
    #print handle.stdout.read()
    #handle.flush()
    return
cmd.extend('emacs',emacs)


def gedit(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with gedit from within PyMOL. 
    Adjust url for your location.
    Can be installed via macports on the mac.

    USAGE
    mate test.pml
    '''
    arg = ("/opt/local/bin/gedit -w " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('gedit',gedit)


def jedit(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with jedit from within PyMOL. 
    Adjust url for your location.
    Can be installed via macports on the mac.

    USAGE
    mate test.pml
    '''
    arg = ("open -a jedit " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('jedit',jedit)


def mate(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with Textmate (Mac OS only) from within PyMOL. 
    Adjust path to Textmate on your computer as needed.

    USAGE
    
    mate test.pml
    '''
    arg = ("/usr/local/bin/mate -w " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('mate',mate)


def notepadpp(fileName="test.pml"):
    '''
    DESCRIPTION

    Open file with notepadpp (Mac OS only) from within PyMOL. 
    Adjust path to notepadpp on your computer as needed.

    USAGE
    
    notepadpp test.pml
    '''
    arg = ("/Applications/Notepad++7.app " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('notepadpp',notepadpp)


def nv(fileName="testme.pml"):
    '''
    DESCRIPTION

    Open file with neovim from within PyMOL. 
    Adjust path to neovim on your computer as needed.

    USAGE
    nv test.pml
    
    Currently opening neovim in new shell and testme.pml in buffer 2. 
    In Normal mode, enter :bnext to switch to second buffer.
    '''
    arg = ("/opt/local/bin/nvim " + fileName)
    subprocess.call(['open', '-W', '-a', 'iTerm.app', '/opt/local/bin/nvim', '--args', fileName])
    #Popen(shlex.split("""x-terminal-emulator -e 'nv -c "testme.pml"'"""), stdout=PIPE)
    #process.wait()
    ########subprocess.call(arg,shell=True)
    #handle = Popen(arg, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
    #print handle.stdout.read()
    #handle.flush()
    return
cmd.extend('nv',nv)


def oni(fileName="test.pml"):
    '''
    DESCRIPTION

    Open the editor Oni from within PyMOL. 
    The is an editor based on neovim.

    USAGE
        oni test.pml
    '''
    arg = ("/Applications/Oni.app/Contents/MacOS/Oni " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('oni',oni)


def pdbed(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open PDBEditor.jar from within PyMOL. 
    Adjust url for your location.
    https://sourceforge.net/projects/pdbeditorjl/

    USAGE
    pdbed test.pdb

    Notes:
    Needs exception for not opening pdb file. 
    '''
    print("Please wait. Editor is slow to start.")
    arg = ("java -jar /Applications/jars/PDB_Editor_FIX090203.jar " + fileName)
    subprocess.call(arg,shell=True)

    return
cmd.extend('pdbed',pdbed)


def st3(fileName="test.pml"):
    '''
    DESCRIPTION

    Open sublime text 3 from within PyMOL. 
    Adjust url for your location.

    USAGE
    st3 test.pml
    '''
    arg = ("/usr/local/bin/subl -w " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('st3',st3)


def vim(fileName="test.pml"):
    '''
    DESCRIPTION

    Open vim from within PyMOL. 
    Adjust file path to vim on your computer.

    USAGE
    vim testme.pml
    '''
    arg = ("/opt/local/bin/vim " + fileName)
    subprocess.call(['open', '-W', '-a', 'Terminal.app', '/opt/local/bin/vim', '--args', fileName])
    return
cmd.extend('vim',vim)


############################### Open word processor ########################################

#category: Open word processord


def word(fileName="Script.docx"):
    '''
    DESCRIPTION

    Open word from within PyMOL. 
    Adjust file path to MS Word on your computer.

    USAGE
    word testme.docx

    The subprocess command may need more work. 
    '''
    arg = ("open " + fileName)
    subprocess.call(['open', '-W', '-a', 'Microsoft Word.app', '--args', fileName])
    return
cmd.extend('word',word)


############################### Open data analysis programs ########################################

#category: Open data analysis programs

def cranR():
    '''
    DESCRIPTION

    Open the Cran R from within PyMOL. 

    USAGE
        cranR
    '''
    arg = ("/usr/local/bin/R")
    subprocess.call(arg,shell=True)
    return
cmd.extend('cranR',cranR)


def ddb():
    '''
    DESCRIPTION

    Open the DBBrowserSQLite. 

    USAGE
        ddb
    '''
    arg = ("open -a DBBrowserSQLite")
    subprocess.call(arg,shell=True)
    return
cmd.extend('ddb',ddb)


def excel():
    '''
    DESCRIPTION

    Open the excel from within PyMOL. 

    USAGE
        cranR
    '''
    arg = ("excel")
    subprocess.call(arg,shell=True)
    return
cmd.extend('excel',excel)


def JASP():
    '''
    DESCRIPTION

    Open the JASP from within PyMOL. 
    The is an editor based on JASP.

    USAGE
        JASP
    '''
    arg = ("/Applications/Jasp.app/Contents/MacOS/JASP " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('cranR',cranR)


def JMP():
    '''
    DESCRIPTION

    Open the JMP from within PyMOL. 

    USAGE
        JMP
    '''
    arg = ("/Applications/JMP\ Pro\ 14.app/Contents/MacOS/JMP")
    subprocess.call(arg,shell=True)
    return
cmd.extend('JMP',JMP)



def jabref():
    '''
    DESCRIPTION

    Open the jabref from within PyMOL. 

    USAGE

        Jabref
    '''
    arg = ('open -a "/Applications/JabRef/JabRef.app"')
    subprocess.call(arg,shell=True)
    return
cmd.extend('jabref',jabref)


def julia():
    '''
    DESCRIPTION

    Open the jabref from within PyMOL. 

    USAGE

        julia
    '''
    arg = ('/Applications/Julia-0.6.app/Contents/MacOS/applet')
    subprocess.call(arg,shell=True)
    return
cmd.extend('julia',julia)


def oc():
    '''
    DESCRIPTION

    Open the jabref from within PyMOL. 

    USAGE

        julia
    '''
    arg = ('octave --no-gui-lib')
    subprocess.call(arg,shell=True)
    return
cmd.extend('oc',oc)


def ppt():
    '''
    DESCRIPTION

    Open the powerpoint from within PyMOL. 

    USAGE

        julia
    '''
    arg = ('open -a /Applications/Microsoft\ PowerPoint.app')
    subprocess.call(arg,shell=True)
    return
cmd.extend('ptt',ppt)




############################### Terminal windows ########################################

#category: Open terminal windows

def iterm():
    '''
    DESCRIPTION

    Open iTerm2 window on MacOS. 

    USAGE
    iterm
    '''
    subprocess.call(['open', '-a', 'iTerm'])
    return
cmd.extend('iterm',iterm)

def term():
    '''
    DESCRIPTION

    Open a Terminal window on MacOS. 

    USAGE
    iterm
    '''
    subprocess.call(['open', '-a', 'Terminal'])
    return
cmd.extend('term',term)



############################### Molecular Graphics  programs ########################################

#category: Open other molecular graphics programs


def ccp4mg(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open ccp4mg from within PyMOL. 
    Adjust url for your location.

    USAGE
    ccp4mg test.pbd
    '''
    arg = ("ccp4mg " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('ccp4mg',ccp4mg)


def chimera(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open Chimera from within PyMOL. 
    Adjust url for your location.

    USAGE
    chimera test.pbd
    '''
    arg = ("/Applications/Chimera.app/Contents/MacOS/chimera " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('chimera',chimera)


def coot(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open coot from within PyMOL. 
    Adjust url for your location.

    USAGE
    coot test.pbd
    '''
    arg = ("/Applications/ccp4-7.0/bin/coot " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('coot',coot)


def jmol(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open Jmol from within PyMOL. 
    Adjust file path for your location of Jmol.

    USAGE
    jmol test.pdb
    '''
    arg = ("sh /Applications/jars/jmol-14.29.16/jmol.sh " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('jmol',jmol)


def vmd(fileName="test.pdb"):
    '''
    DESCRIPTION

    Open vmd from within PyMOL. 
    Adjust url for your location.

    USAGE
    vmd test.pml
    '''
    arg = ("/Applications/VMD194.app/Contents/MacOS/startup.command " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('vmd',vmd)


def yasara(fileName="test.pml"):
    '''
    DESCRIPTION

    Open the molecular graphics prograom YASASRA from within PyMOL. 

    USAGE
        yasara test.pml
    '''
    arg = ("/Applications/YASARA.app/Contents/MacOS/yasara.app " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('yasara',yasara)

############################### Graphics programs ########################################

#category: Image manipulation programs


def gimp(fileName="test.png"):
    '''
    DESCRIPTION

    Open the molecular graphics program with gimp from within PyMOL. 

    USAGE
        gimp test.png
    '''
    arg = ("/opt/local/bin/gimp " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('gimp',gimp)


def inkscape(fileName="test.svg"):
    '''
    DESCRIPTION

    Open the molecular graphics program with gimp from within PyMOL. 

    USAGE
        inkscape test.png
    '''
    arg = ("/opt/local/bin/inkscape " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('inkscape',inkscape)


############################### Webapps ########################################

#category: Open certain webapps

def gcal():
    '''
    DESCRIPTION

    Open Google Calendar. 

    USAGE

    GM
    '''
    webbrowser.open('https://calendar.google.com/calendar/r')
cmd.extend('gcal',gcal)


def GM():
    '''
    DESCRIPTION

    Open gmail. 

    USAGE

    GM
    '''
    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')
cmd.extend('GM',GM)


def WM():
    '''
    DESCRIPTION
    
    Open Web Mail in defualt browser. Adjust url for your institution.

    USAGE
    WM
    '''
    webbrowser.open('https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f')
cmd.extend('WM',WM)


def WS():
    '''
    DESCRIPTION

    Open National Weather Service website for locale. 
    Adjust url for your location.

    USAGE
    WS
    '''
    webbrowser.open('https://www.weather.gov/oun/')
cmd.extend('WS',WS)


############################### Samples ########################################

#category: Samples

def GGT():
    '''
    DESCRIPTION

    WT human gamma glutamyl transpeptidase at 1.67 Angstrom
    resolution as cartoon. PDB Code 4gdx.
    
    USAGE

    Type 'GGT' to activate. Type 'help GGT' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    fetch 4gdx, type=pdb, async=0;
    remove name H*;
    as cartoon;
    bg_color white; 
    hide (name c+o+n);
    set cartoon_side_chain_helper, on;color red, 4gdx and ss H; 
    color yellow,4gdx and ss S;color green,4gdx and ss L+; 
    select ASNNAG, resn NAG or resi 95 or i. 120  or i. 230 or
    i. 266 or i. 344 or i. 511 or i. 381; 
    color red, elem o and ASNNAG; 
    color blue, elem n and ASNNAG;
    color yellow, elem c and ASNNAG;
    show sticks,ASNNAG;
    disable ASNNAG;
    set_view (0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); 
    draw 
    
    The commands without linebreaks:
    
    delete all;fetch 4gdx, type=pdb, async=0;remove name H*;as cartoon;bg_color white; hide (name c+o+n);set cartoon_side_chain_helper,  on;color red, 4gdx and ss H; color yellow,4gdx and ss S;color green,4gdx and ss L+; select ASNNAG,resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 ori. 511 or i. 381; color red, elem o and ASNNAG; color blue, elem n and ASNNAG;color yellow, elem c  and ASNNAG;show sticks,ASNNAG;disable ASNNAG; set_view(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); draw 
    
    '''
    cmd.reinitialize()
    cmd.fetch('4gdx', type='pdb', async='0')
    cmd.remove('name H*')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.hide('(name c+o+n)')
    cmd.set('cartoon_side_chain_helper', 'on')
    cmd.color('red', '4gdx and ss H')
    cmd.color('yellow', '4gdx and ss S')
    cmd.color('green', '4gdx and ss L+')
    cmd.select('ASNNAG', 'resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 or i. 511 or i. 381')
    cmd.color('red', 'elem o and ASNNAG')
    cmd.color('blue', 'elem n and ASNNAG')
    cmd.color('yellow', 'elem c and ASNNAG')
    cmd.show('sticks', 'ASNNAG')
    cmd.disable('ASNNAG')
    cmd.set_view('(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0)')
    cmd.draw()
cmd.extend('GGT',GGT)


def GU():
    '''
    DESCRIPTION

    10-mer dsRNA with 8 contiguous Us. U-helix RNA. 
    1.32 Angstrom resolution: 4PCO. Has five strands in 
    the asymmetric unit. Deleted chain E and cobalt 
    hexammine 102. Cartoon with filled rings and
    bases cartoon.
    
    
    USAGE

    Type 'GU' to activate. 
    
    Type 'help GU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    fetch 4PCO,type=pdb,async=0;
    hide everything; 
    bg_color white; 
    cartoon oval; 
    set cartoon_ring_mode, 3; 
    set cartoon_nucleic_acid_color, blue;
    select rna_A, resn A; 
    select rna_C, resn C;
    select rna_G, resn G; 
    select rna_U, resn U;
    color yellow, rna_A; 
    color red, rna_C;
    color gray40, rna_G; 
    color palecyan, rna_U;
    as cartoon; 
    disable rna_U; 
    set stick_radius, 0.12;
    set nb_spheres_size, 0.3; 
    show nb_spheres; 
    set stick_ball, on;
    set stick_ball_ratio, 1.8; 
    show sticks, resn NCO; 
    show spheres, name Cl; 
    set_view (0.34,-0.81,0.48,0.89,0.11,
    -0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,
    233.71,-20.0);
    draw 
    
    The commands without linebreaks: 
    
    delete all;fetch 4PCO,type=pdb,async=0;hide everything;bg_color white; cartoon oval;set cartoon_ring_mode, 3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C,resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A; color red, rna_C;color gray40, rna_G; color palecyan, rna_U;as cartoon;disable rna_U; set stick_radius, 0.12;set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;set stick_ball_ratio, 1.8; show sticks, resn NCO;show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0);draw 
    
    '''
    cmd.reinitialize();
    cmd.fetch('4PCO', type='pdb', async='0')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.cartoon('oval')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_nucleic_acid_color', 'blue')
    cmd.select('rna_A', 'resn A')
    cmd.select('rna_C', 'resn C')
    cmd.select('rna_G', 'resn G')
    cmd.select('rna_U', 'resn U')
    cmd.color('yellow', 'rna_A')
    cmd.color('red', 'rna_C')
    cmd.color('gray40', 'rna_G')
    cmd.color('palecyan', 'rna_U')
    cmd.show_as('cartoon')
    cmd.disable('rna_U')
    cmd.set('stick_radius', '0.12')
    cmd.set('nb_spheres_size', '0.3')
    cmd.show('nb_spheres')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.8')
    cmd.show('sticks', 'resn NCO')
    cmd.show('spheres', 'name Cl')
    cmd.set_view('(0.34,-0.81, 0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0)')
    cmd.draw()
cmd.extend('GU',GU)


def N9():
    '''
    DESCRIPTION
    
    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr.
    The biological unit has four copies of the asymmetric unit.
    View is down the four-fold axis. Requires the quat.py script by
    Thomas Holder and available at the PyMOL Wiki page. Store quat.py
    in ~/mg18OU.

    USAGE

    Type 'N9' to activate. Type 'help N9' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:

    delete all;
    fetch 4dgr, type=pdb, async=0;
    run $HOME/mg18OU/quat.py;
    quat 4dgr;
    as cartoon; 
    bg_color white;
    color red, 4dgr_1 and ss H;
    color yellow,4dgr_1 and ss S;
    color green, 4dgr_1 and ss L+;
    color cyan, (not 4dgr_1 and ss H);
    color magenta, (not 4dgr_1 and ss S);
    color orange, (not 4dgr_1 and ss L+);
    set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0);
    draw 

    The commands without linebreaks:

    delete all;fetch 4dgr, type=pdb, async=0;run $HOME/mg18OU/quat.py; quat 4dgr;as cartoon; bg_color white;color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;color green, 4dgr_1 and ss L+;color cyan, (not 4dgr_1 and ss H);color magenta, (not 4dgr_1 and ss S);color orange, (not 4dgr_1 and ss L+);set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0); draw 

    '''
    cmd.reinitialize()
    cmd.fetch('4dgr', type='pdb', async='0')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 4dgr')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.color('red', '4dgr_1 and ss H')
    cmd.color('yellow', '4dgr_1 and ss S')
    cmd.color('green', '4dgr_1 and ss L+')
    cmd.color('cyan', '(not 4dgr_1 and ss H)')
    cmd.color('magenta', '(not 4dgr_1 and ss S)')
    cmd.color('orange', '(not 4dgr_1 and ss L+)')
    cmd.set_view('(0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0)')
    cmd.draw()
cmd.extend('N9',N9)


def T4L():
    '''
    DESCRIPTION
    
    WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0. 
    
    USAGE

    Type 'T4L' to activate. Type 'help T4L' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
        
    delete all;
    fetch 3fa0,type=pdb,async=0;
    orient;
    turn z,-90;
    turn y,-5;
    turn x,10; 
    hide everything; 
    bg_color white; 
    show cartoon;
    color red, ss H;
    color yellow, ss S;
    color green, ss L+;
    set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); 
    ray 1500,1600;

    The commands without linebreaks:
    
    delete all;fetch 3fa0,type=pdb,async=0;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white;show cartoon;color red, ss H;color yellow, ss S;color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); ray 1500,1600; 
    
    '''
    cmd.reinitialize()
    cmd.fetch('3fa0', type='pdb', async='0')
    cmd.orient()
    cmd.turn('z', '-90')
    cmd.turn('y', '-5')
    cmd.turn('x', '10')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')
    cmd.set_view('(-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0)')
    cmd.ray('1500', '1600')
    cmd.png("T4L.png")
cmd.extend('T4L',T4L)


def U8():
    '''
    DESCRIPTION

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.
    Has one strand in the asymmetric unit. Uses quat.py to generate
    the second strand. Cartoon with filled rings and bases cartoon.
    
    USAGE

    Type 'U8' to activate. Type 'help U8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    fetch 3nd3,type=pdb,async=0;
    run $HOME/mg18OU/quat.py;
    quat 3nd3;
    hide everything;
    bg_color white; 
    show sticks;
    set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;
    set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;
    as cartoon;
    set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 
    
    The commands without linebreaks:
    
    delete all;fetch 3nd3,type=pdb,async=0;run $HOME/mg18OU/quat.py;quat 3nd3;hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon;set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 

    '''
    
    cmd.reinitialize()
    cmd.fetch('3nd3', type='pdb', async='0')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd3')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')
    cmd.set_view('(-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0)')
    cmd.draw()
cmd.extend('U8',U8)


def WC8():
    '''
    DESCRIPTION

    16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom 
    resolution: 3nd4.  Has one strand in the asymmetric unit. 
    Needs quat.py to generate the second strand. Use the 
    BU alias. Cartoon with filled rings and bases cartoon.
    
    
    USAGE

    Type 'WC8' to activate. Type 'help WC8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all; 
    fetch 3nd4,type=pdb,async=0;
    hide everything;
    run $HOME/mg18OU/quat.py;
    quat 3nd4;
    bg_color white; 
    show sticks; 
    set stick_radius, 0.12; 
    set nb_spheres_size, 0.25;
    show nb_spheres;
    set stick_ball, on; 
    set stick_ball_ratio, 1.8;
    set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);
    hide everything,name H*;
    rock

    The commands without linebreaks:
    
    delete all; fetch 3nd4,type=pdb,async=0;hide everything; run $HOME/mg18OU/quat.py; quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; set nb_spheres_size, 0.25; show nb_spheres; set stick_ball, on; set stick_ball_ratio, 1.8;set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);hide everything, name H*;rock 

    '''
    cmd.reinitialize()
    cmd.fetch('3nd4', type='pdb', async='0')
    cmd.remove('name H*')
    cmd.hide('everything')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd4')
    cmd.bg_color('white')
    cmd.do('show stick')
    cmd.do('set stick_radius, 0.12') 
    cmd.do('set nb_spheres_size, 0.25')
    cmd.do('show nb_spheres')
    cmd.do('set stick_ball, on')
    cmd.do('set stick_ball_ratio, 1.8')
    cmd.set_view('(-0.96,-0.03,0.3,-0.31,0.02,-0.95,0.03,-1.0,-0.03,0.0,0.0,-231.24,8.16,15.68,-1.66,200.47,262.01,-20.0)')
    cmd.rock()
cmd.extend('WC8',WC8)


####### Commands to display complex scenes. #############

#category: Commands to display complex scenes.

def BST():
    '''
    DESCRIPTION
    
    G2G3/U9U8 base step , PDB code 4PCO. 
    From the 1.32 Angstrom resolution structure 
    of the RNA decamer with 8 GU base pairs.
    
    USAGE

    Type 'BST' to execute. Type 'help BST' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    fetch 4PCO, type=pdb, async=0;
    select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));
    remove not G2G3;
    bg_color white;
    show sticks;
    set stick_radius=0.14;
    set stick_ball, on; 
    set stick_ball_ratio,1.9;
    set_view 
    (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0);
    remove name H*;
    select carbon1, element C and (resi 3 or resi 8) 
    # select lower base pair;
    select carbon2, element C and (resi 2 or resi 9) 
    #select upper base pair;
    color gray70, carbon1;
    color gray10, carbon2;
    show sticks;
    space cmyk;
    distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;
    distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;
    distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;
    distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;
    color black, hbond1;
    color black, hbond2;
    color gray70, hbond3;
    color gray70, hbond4;
    show nb_spheres;
    set nb_spheres_size, 0.35;
    hide labels;
    ray 1600,1000;
    png 4PCO.png
    
    Commands without linebreaks: 
    
    delete all;fetch 4PCO, type=pdb, async=0;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png

    '''
    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb', async='0')
    cmd.select('G2G3', '( ((resi 2 or resi 3) and chain A)or ((resi 8 or resi 9) and chain B) )')
    cmd.remove('not G2G3')
    cmd.bg_color('white')
    cmd.set('stick_radius', '0.14')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.9')
    cmd.set_view('(-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0)')
    cmd.remove('name H*')
    cmd.select('carbon1', 'element C and (resi 3 or resi 8)')
    cmd.select('carbon2', 'element C and (resi 2 or resi 9)')
    cmd.color('gray70', 'carbon1')
    cmd.color('gray10', 'carbon2')
    cmd.show('sticks')
    cmd.space('cmyk')
    cmd.distance('hbond1', '/4PCO//B/U`9/N3', '/4PCO//A/G`2/O6')
    cmd.distance('hbond2', '/4PCO//B/U`9/O2', '/4PCO//A/G`2/N1')
    cmd.distance('hbond3', '/4PCO//A/U`3/N3', '/4PCO//B/G`8/O6')
    cmd.distance('hbond4', '/4PCO//A/U`3/O2', '/4PCO//B/G`8/N1')
    cmd.color('black', 'hbond1')
    cmd.color('black', 'hbond2')
    cmd.color('gray70', 'hbond3')
    cmd.color('gray70', 'hbond4')
    cmd.show('nb_spheres')
    cmd.set('nb_spheres_size', '0.35')
    cmd.hide('labels')
    cmd.ray('1600', '1000')
    cmd.png('4PCO.png')
cmd.extend('BST',BST)


def LG():
    '''
    DESCRIPTION
    
    Nine sugar glycan in influenza N9 neuraminidase at 
    1.55 Angstrom  resolution, PDB code 4dgr. 
    The electron density map is contoured at 1.0 sigma. 
    39 commands were used to make this figure.  
    
    USAGE

    Type 'LG' to execute. Type 'help LG' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    fetch 4dgr, async=0;
    fetch 4dgr, type=2fofc,async=0;
    select LongGlycan, resi 469:477;
    orient LongGlycan;
    remove not LongGlycan;
    remove name H*;
    isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;
    color density, 2fofcmap; 
    show sticks;
    show spheres;
    set stick_radius, .07;
    set sphere_scale, .19;
    set sphere_scale, .13, elem H;
    set bg_rgb=[1, 1, 1];
    set stick_quality, 50;
    set sphere_quality, 4;
    color gray85, elem C;
    color red, elem O;
    color slate, elem N;
    color gray98, elem H;
    set stick_color, gray50;
    set ray_trace_mode, 1;
    set ray_texture, 2;
    set antialias, 3;
    set ambient, 0.5;
    set spec_count, 5;
    set shininess, 50;
    set specular, 1;
    set reflect, .1;
    set dash_gap, 0;
    set dash_color, black;
    set dash_gap, .15;
    set dash_length, .05;
    set dash_round_ends, 0;
    set dash_radius, .05;
    set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);
    preset.ball_and_stick("all",mode=1);
    draw 
    
    Commands without linebreaks:
    
    delete all;fetch 4dgr, async=0;fetch 4dgr, type=2fofc, async=0;select LongGlycan, resi 469:477;orient LongGlycan;remove not LongGlycan;remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;color density, 2fofcmap; show sticks;show spheres;set stick_radius, .07;set sphere_scale, .19;set sphere_scale, .13, elem H;set bg_rgb=[1, 1, 1];set stick_quality, 50;set sphere_quality, 4;color gray85, elem C;color red, elem O;color slate, elem N;color gray98, elem H;set stick_color, gray50;set ray_trace_mode, 1;set ray_texture, 2;set antialias, 3;set ambient, 0.5;set spec_count, 5;set shininess, 50;set specular, 1;set reflect, .1;set dash_gap, 0;set dash_color, black;set dash_gap, .15;set dash_length, .05;set dash_round_ends, 0;set dash_radius, .05;set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);preset.ball_and_stick("all",mode=1);draw 
 
    '''
    cmd.reinitialize()
    cmd.fetch('4dgr', async='0')
    cmd.fetch('4dgr', type='2fofc', async='0')
    cmd.select('LongGlycan', 'resi 469:477')
    cmd.orient('LongGlycan')
    cmd.remove('not LongGlycan')
    cmd.remove('name H*')
    cmd.isomesh('2fofcmap', '4dgr_2fofc', '1', 'LongGlycan', carve ='1.8')
    cmd.color('density', '2fofcmap')
    cmd.show('sticks')
    cmd.show('spheres')
    cmd.set('stick_radius', '.07')
    cmd.set('sphere_scale', '.19')
    cmd.set('sphere_scale', '.13', 'elem H')
    cmd.set('bg_rgb', '[1, 1, 1]')
    cmd.set('stick_quality', '50')
    cmd.set('sphere_quality', '4')
    cmd.color('gray85', 'elem C')
    cmd.color('red', 'elem O')
    cmd.color('slate', 'elem N')
    cmd.color('gray98', 'elem H')
    cmd.set('stick_color', 'gray50')
    cmd.set('ray_trace_mode', '1')
    cmd.set('ray_texture', '2')
    cmd.set('antialias', '3')
    cmd.set('ambient', '0.5')
    cmd.set('spec_count', '5')
    cmd.set('shininess', '50')
    cmd.set('specular', '1')
    cmd.set('reflect', '.1')
    cmd.set('dash_gap', '0')
    cmd.set('dash_color', 'black')
    cmd.set('dash_gap', '.15')
    cmd.set('dash_length', '.05')
    cmd.set('dash_round_ends', '0')
    cmd.set('dash_radius', '.05')
    cmd.set_view('(0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0)')
    preset.ball_and_stick("all",mode=1);
    cmd.draw()
cmd.extend('LG',LG)


def NA():
    '''
    DESCRIPTION

    Hydrated sodium cation bound in major groove of a 
    16-mer RNA of Watson-Crick base pairs.
    The sodium is bound to the N7 nitrogen atom of 
    Adenine 3 at 1.55 Angstrom resolution, PDB code 3nd4. 
    57 commands were used to make this figure. 

    More than one label in a horizontal script is not 
    allowed. This one label has to be at the end of the line.
    Labels can be imported from a Labels.pml file.
    Store the label commands one per row in this file.
    Import the file with the @Labels.pml command. 
    Include the path to the file if the labels file is not 
    in the current working directory of PyMOL. 
    
    USAGE

    Type 'NA' to execute. Type 'help NA' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    viewport 900,600;
    fetch 3nd4, type=pdb, async=0;
    run ~/mg18OU/quat.py; 
    quat 3nd4; 
    show sticks;
    set stick_radius=0.125;
    hide everything, name H*;
    bg_color white;
    create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));
    bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 519);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 319);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);
    show nb_spheres; 
    set nb_spheres_size, .35;
    distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;
    distance hbond2,/3nd4_1/1/A/HOH`319/O,/3nd4_1/1/A/A`3/OP2;
    distance hbond3,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/HOH`119/O;
    distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;
    distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;
    distance hbond6,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    distance hbond7,/3nd4_1/1/A/HOH`319/O,/3nd4_1/1/A/G`2/OP2;
    distance hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;
    distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;
    distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;
    distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;
    distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;
    distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;
    distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;
    distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;
    distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;
    distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;
    distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    set depth_cue=0;
    set ray_trace_fog=0;
    set dash_color, black;
    set label_font_id, 5;
    set label_size, 36;
    set label_position, (0.5, 1.0, 2.0);
    set label_color, black;
    set dash_gap, 0.2;
    set dash_width, 2.0;
    set dash_length, 0.2;
    set label_color, black;
    set dash_gap, 0.2;
    set dash_width, 2.0;
    set dash_length, 0.2;
    select carbon, element C; 
    color yellow, carbon;
    disable carbon;
    set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0)
    
    The commands without linebreaks:
    
    delete all;viewport 900,600;fetch 3nd4, type=pdb,async=0;run ~/mg18OU/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 

    '''
    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.fetch('3nd4', type='pdb', async='0');
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd4');
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'name H*');
    cmd.bg_color('white');
    cmd.create('coorCov', '(3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)))');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/A`3/N7)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`119/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`219/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`319/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`419/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`519/O)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 519)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 419)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 119)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 319)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 219)');
    cmd.show('nb_spheres');
    cmd.set('nb_spheres_size', '.35');
    cmd.distance('hbond1', '/3nd4_1/1/A/HOH`119/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond2', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond3', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/HOH`119/O');
    cmd.distance('hbond4', '/3nd4_1/1/A/G`4/N7', '/3nd4_1/1/A/HOH`91/O');
    cmd.distance('hbond5', '/3nd4_1/1/A/G`4/O6', '/3nd4_1/1/A/HOH`419/O');
    cmd.distance('hbond6', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.distance('hbond7', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/G`2/OP2');
    cmd.distance('hbond9', '/3nd4_1/1/A/HOH`419/O', '/3nd4_2/2/A/HOH`74/O');
    cmd.distance('hbond10', '/3nd4_2/2/A/C`15/O2', '/3nd4_1/1/A/G`2/N2');
    cmd.distance('hbond11', '/3nd4_2/2/A/C`15/N3', '/3nd4_1/1/A/G`2/N1');
    cmd.distance('hbond12', '/3nd4_2/2/A/C`15/N4', '/3nd4_1/1/A/G`2/O6');
    cmd.distance('hbond13', '/3nd4_2/2/A/U`14/N3', '/3nd4_1/1/A/A`3/N1');
    cmd.distance('hbond14', '/3nd4_2/2/A/U`14/O4', '/3nd4_1/1/A/A`3/N6');
    cmd.distance('hbond15', '/3nd4_2/2/A/C`13/N4', '/3nd4_1/1/A/G`4/O6');
    cmd.distance('hbond16', '/3nd4_2/2/A/C`13/N3', '/3nd4_1/1/A/G`4/N1');
    cmd.distance('hbond17', '/3nd4_1/1/A/G`4/N2', '/3nd4_2/2/A/C`13/O2');
    cmd.distance('hbond18', '/3nd4_1/1/A/G`2/N2', '/3nd4_2/2/A/C`15/O2');
    cmd.distance('hbond19', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.set('depth_cue', '0');
    cmd.set('ray_trace_fog', '0');
    cmd.set('dash_color', 'black');
    cmd.set('label_font_id', '5');
    cmd.set('label_size', '36')
    cmd.set('label_position', '(0.5, 1.0,2.0)');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.select('carbon', 'element C');
    cmd.color('yellow', 'carbon');
    cmd.disable('carbon');
    cmd.set_view('-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0');
cmd.extend('NA',NA)


################# VARIANTS OF The ABOVE that are NOT internet dependent  #################

#category: Commands to display complex scenes with pdb files on computer.

def LGGT():
    '''
    DESCRIPTION

    WT human gamma glutamyl transpeptidase at 1.67 Angstrom
    resolution as cartoon. PDB code 4gdx. 
    4gdx.pdb must be in the current working directory. 
    
    USAGE

    Type 'LGGT' to activate. Type 'help LGGT' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    load 4gdx.pdb;
    remove name H*;
    as cartoon;
    bg_color white; 
    hide (name c+o+n);
    set cartoon_side_chain_helper, on;color red, 4gdx and ss H; 
    color yellow,4gdx and ss S;color green,4gdx and ss L+; 
    select ASNNAG, resn NAG or resi 95 or i. 120  or i. 230 or
    i. 266 or i. 344 or i. 511 or i. 381; 
    color red, elem o and ASNNAG; 
    color blue, elem n and ASNNAG;
    color yellow, elem c and ASNNAG;
    show sticks,ASNNAG;
    disable ASNNAG;
    set_view (0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); 
    draw 
    
    The commands without linebreaks:
    
    delete all;load 4gdx.pdb;remove name H*;as cartoon;bg_color white; hide (name c+o+n);set cartoon_side_chain_helper,  on;color red, 4gdx and ss H; color yellow,4gdx and ss S;color green,4gdx and ss L+; select ASNNAG,resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 ori. 511 or i. 381; color red, elem o and ASNNAG; color blue, elem n and ASNNAG;color yellow, elem c  and ASNNAG;show sticks,ASNNAG;disable ASNNAG; set_view(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); draw 
    
    '''
    cmd.reinitialize()
    cmd.load('4gdx.pdb')
    cmd.remove('name H*')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.hide('(name c+o+n)')
    cmd.set('cartoon_side_chain_helper', 'on')
    cmd.color('red', '4gdx and ss H')
    cmd.color('yellow', '4gdx and ss S')
    cmd.color('green', '4gdx and ss L+')
    cmd.select('ASNNAG', 'resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 or i. 511 or i. 381')
    cmd.color('red', 'elem o and ASNNAG')
    cmd.color('blue', 'elem n and ASNNAG')
    cmd.color('yellow', 'elem c and ASNNAG')
    cmd.show('sticks', 'ASNNAG')
    cmd.disable('ASNNAG')
    cmd.set_view('(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0)')
    cmd.draw()
cmd.extend('LGGT',LGGT)


def LGU():
    '''
    DESCRIPTION

    10-mer dsRNA. 
    1.32 Angstrom resolution: 4PCO. Has five strands in 
    the asymmetric unit. Deleted chain E and cobalt 
    hexammine 102. Cartoon with filled rings and
    bases cartoon.
    
    
    USAGE

    Type 'GU' to activate. Type 'help GU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    load 4PCO.pdb;
    hide everything; 
    bg_color white; 
    cartoon oval; 
    set cartoon_ring_mode, 3; 
    set cartoon_nucleic_acid_color, blue;
    select rna_A, resn A; 
    select rna_C, resn C;
    select rna_G, resn G; 
    select rna_U, resn U;
    color yellow, rna_A; 
    color red, rna_C;
    color gray40, rna_G; 
    color palecyan, rna_U;
    as cartoon; 
    disable rna_U; 
    set stick_radius, 0.12;
    set nb_spheres_size, 0.3; 
    show nb_spheres; 
    set stick_ball, on;
    set stick_ball_ratio, 1.8; 
    show sticks, resn NCO; 
    show spheres, name Cl; 
    set_view (0.34,-0.81,0.48,0.89,0.11,
    -0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,
    233.71,-20.0);
    draw 
    
    The commands without linebreaks: 
    
    delete all;load 4PCO.pdb,async=0;hide everything;bg_color white; cartoon oval;set cartoon_ring_mode, 3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C,resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A; color red, rna_C;color gray40, rna_G; color palecyan, rna_U;as cartoon;disable rna_U; set stick_radius, 0.12;set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;set stick_ball_ratio, 1.8; show sticks, resn NCO;show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0);draw 
    
    '''
    cmd.reinitialize();
    cmd.load('4PCO.pdb')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.cartoon('oval')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_nucleic_acid_color', 'blue')
    cmd.select('rna_A', 'resn A')
    cmd.select('rna_C', 'resn C')
    cmd.select('rna_G', 'resn G')
    cmd.select('rna_U', 'resn U')
    cmd.color('yellow', 'rna_A')
    cmd.color('red', 'rna_C')
    cmd.color('gray40', 'rna_G')
    cmd.color('palecyan', 'rna_U')
    cmd.show_as('cartoon')
    cmd.disable('rna_U')
    cmd.set('stick_radius', '0.12')
    cmd.set('nb_spheres_size', '0.3')
    cmd.show('nb_spheres')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.8')
    cmd.show('sticks', 'resn NCO')
    cmd.show('spheres', 'name Cl')
    cmd.set_view('(0.34,-0.81, 0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0)')
    cmd.draw()
cmd.extend('LGU',LGU)


def LN9():
    '''
    DESCRIPTION
    
    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code
    4dgr. The biological unit has four copies of the asymmetric unit.
    View is down the four-fold axis. Requires the quat.py script by
    Thomas Holder and available at the PyMOL Wiki page. Store quat.py
    in ~/mg18OU.

    USAGE

    Type 'LN9' to activate. Type 'help LN9' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:

    delete all;
    load 4dgr.pdb;
    run $HOME/mg18OU/quat.py;
    quat 4dgr;
    as cartoon; 
    bg_color white;
    color red, 4dgr_1 and ss H;
    color yellow,4dgr_1 and ss S;
    color green, 4dgr_1 and ss L+;
    color cyan, (not 4dgr_1 and ss H);
    color magenta, (not 4dgr_1 and ss S);
    color orange, (not 4dgr_1 and ss L+);
    set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0);
    draw 

    The commands without linebreaks:

    delete all;load 4dgr.pdb;run $HOME/mg18OU/quat.py; quat 4dgr;as cartoon; bg_color white;color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;color green, 4dgr_1 and ss L+;color cyan, (not 4dgr_1 and ss H);color magenta, (not 4dgr_1 and ss S);color orange, (not 4dgr_1 and ss L+);set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0); draw 

    '''
    cmd.reinitialize()
    cmd.load('4dgr.pdb')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 4dgr')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.color('red', '4dgr_1 and ss H')
    cmd.color('yellow', '4dgr_1 and ss S')
    cmd.color('green', '4dgr_1 and ss L+')
    cmd.color('cyan', '(not 4dgr_1 and ss H)')
    cmd.color('magenta', '(not 4dgr_1 and ss S)')
    cmd.color('orange', '(not 4dgr_1 and ss L+)')
    cmd.set_view('(0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0)')
    cmd.draw()
cmd.extend('LN9',LN9)


def LT4L():
    '''
    DESCRIPTION
    
    Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0. 
    The file 3FA0 must be in the current working directory. 
    
    USAGE

    Type 'LT4L' to activate. Type 'help LT4L' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
        
    delete all;
    load 3fa0.pdb;
    orient;
    turn z,-90;
    turn y,-5;
    turn x,10; 
    hide everything; 
    bg_color white; 
    show cartoon;
    color red, ss H;
    color yellow, ss S;
    color green, ss L+;
    set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); 
    ray 1500,1600;

    The commands without linebreaks:
    
    delete all;load 3fa0.pdb;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white;show cartoon;color red, ss H;color yellow, ss S;color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); ray 1500,1600; 
    
    '''
    cmd.reinitialize()
    cmd.fetch('3fa0', type='pdb', async='0')
    cmd.orient()
    cmd.turn('z', '-90')
    cmd.turn('y', '-5')
    cmd.turn('x', '10')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')
    cmd.set_view('(-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0)')
    cmd.ray('1500', '1600')
    cmd.png("T4L.png")
cmd.extend('LT4L',LT4L)


def LU8():
    '''
    DESCRIPTION

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.
    Has one strand in the asymmetric unit. Uses quat.py to generate
    the second strand. Cartoon with filled rings and bases cartoon.
    The file 3nd3.pdb needs to be in the current working directory.

    USAGE

    Type 'LU8' to activate. Type 'help LU8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    load 3nd3.pdb;
    run $HOME/mg18OU/quat.py;
    quat 3nd3;
    hide everything;
    bg_color white; 
    show sticks;
    set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;
    set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;
    as cartoon;
    set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 
    
    The commands without linebreaks:
    
    delete all;load 3nd3.pdb;run $HOME/mg18OU/quat.py;quat 3nd3;hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon;set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 

    '''
    
    cmd.reinitialize()
    cmd.load('3nd3.pdb')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd3')
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')
    cmd.set_view('(-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0)')
    cmd.draw()
cmd.extend('LU8',LU8)


def LWC8():
    '''
    DESCRIPTION

    16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom 
    resolution: 3nd4.  Has one strand in the asymmetric unit. 
    Needs quat.py to generate the second strand. 
    Cartoon with filled rings and bases cartoon.
    The file 3nd4.pdb must be in the current working directory
    
    USAGE

    Type 'LWC8' to activate. Type 'help LWC8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all; 
    load 3nd4.pdb;
    hide everything;
    run $HOME/mg18OU/quat.py;
    quat 3nd4;
    bg_color white; 
    show sticks; 
    set stick_radius, 0.12; 
    set nb_spheres_size, 0.25;
    show nb_spheres;
    set stick_ball, on; 
    set stick_ball_ratio, 1.8;
    set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);
    hide everything,name H*;
    rock

    The commands without linebreaks:
    
    delete all;load 3nd4.pdb;hide everything;run $HOME/mg18OU/quat.py; quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; set nb_spheres_size, 0.25; show nb_spheres; set stick_ball, on; set stick_ball_ratio, 1.8;set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);hide everything, name H*;rock 

    '''
    cmd.reinitialize()
    cmd.load('3nd4.pdb')
    cmd.remove('name H*')
    cmd.hide('everything')
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd4')
    cmd.bg_color('white')
    cmd.do('show stick')
    cmd.do('set stick_radius, 0.12') 
    cmd.do('set nb_spheres_size, 0.25')
    cmd.do('show nb_spheres')
    cmd.do('set stick_ball, on')
    cmd.do('set stick_ball_ratio, 1.8')
    cmd.set_view('(-0.96,-0.03,0.3,-0.31,0.02,-0.95,0.03,-1.0,-0.03,0.0,0.0,-231.24,8.16,15.68,-1.66,200.47,262.01,-20.0)')
    cmd.rock()
cmd.extend('LWC8',LWC8)


def LBST():
    '''
    DESCRIPTION
    
    G2G3/U9U8 base step , PDB code 4PCO. 
    From the 1.32 Angstrom resolution structure 
    of an RNA decamer with 8 GU base pairs.
    
    USAGE

    Type 'BST' to execute. Type 'help BST' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    load 4pco.pdb;
    select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));
    remove not G2G3;
    bg_color white;
    show sticks;
    set stick_radius=0.14;
    set stick_ball, on; 
    set stick_ball_ratio,1.9;
    set_view 
    (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0);
    remove name H*;
    select carbon1, element C and (resi 3 or resi 8) 
    # select lower base pair;
    select carbon2, element C and (resi 2 or resi 9) 
    #select upper base pair;
    color gray70, carbon1;
    color gray10, carbon2;
    show sticks;
    space cmyk;
    distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;
    distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;
    distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;
    distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;
    color black, hbond1;
    color black, hbond2;
    color gray70, hbond3;
    color gray70, hbond4;
    show nb_spheres;
    set nb_spheres_size, 0.35;
    hide labels;
    ray 1600,1000;
    png 4pco.png
    
    Commands without linebreaks: 
    
    delete all;load 4PCO.pdb;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png

    '''
    cmd.reinitialize()
    cmd.load('4pco.pdb')
    cmd.select('G2G3', '( ((resi 2 or resi 3) and chain A)or ((resi 8 or resi 9) and chain B) )')
    cmd.remove('not G2G3')
    cmd.bg_color('white')
    cmd.set('stick_radius', '0.14')
    cmd.set('stick_ball', 'on')
    cmd.set('stick_ball_ratio', '1.9')
    cmd.set_view('(-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7.24,9.55,11.78,29.46,57.91,-20.0)')
    cmd.remove('name H*')
    cmd.select('carbon1', 'element C and (resi 3 or resi 8)')
    cmd.select('carbon2', 'element C and (resi 2 or resi 9)')
    cmd.color('gray70', 'carbon1')
    cmd.color('gray10', 'carbon2')
    cmd.show('sticks')
    cmd.space('cmyk')
    cmd.distance('hbond1', '/4PCO//B/U`9/N3', '/4PCO//A/G`2/O6')
    cmd.distance('hbond2', '/4PCO//B/U`9/O2', '/4PCO//A/G`2/N1')
    cmd.distance('hbond3', '/4PCO//A/U`3/N3', '/4PCO//B/G`8/O6')
    cmd.distance('hbond4', '/4PCO//A/U`3/O2', '/4PCO//B/G`8/N1')
    cmd.color('black', 'hbond1')
    cmd.color('black', 'hbond2')
    cmd.color('gray70', 'hbond3')
    cmd.color('gray70', 'hbond4')
    cmd.show('nb_spheres')
    cmd.set('nb_spheres_size', '0.35')
    cmd.hide('labels')
    cmd.ray('1600', '1000')
    cmd.png('4PCO.png')
cmd.extend('LBST',LBST)


def LLG():
    '''
    DESCRIPTION
    
    Nine sugar glycan in influenza N9 neuraminidase at 
    1.55 Angstrom  resolution, PDB code 4dgr. 
    The electron density map is contoured at 1.0 sigma. 
    39 commands were used to make this figure.  
    
    USAGE

    Type 'LLG' to execute. Type 'help LLG' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    load 4dgr.pdb;
    load 4dgr2FoFc.ccp4;
    select LongGlycan, resi 469:477;
    orient LongGlycan;
    remove not LongGlycan;
    remove name H*;
    isomesh 2fofcmap, 4dgr2FoFc, 1, LongGlycan, carve = 1.8;
    color density, 2fofcmap; 
    show sticks;
    show spheres;
    set stick_radius, .07;
    set sphere_scale, .19;
    set sphere_scale, .13, elem H;
    set bg_rgb=[1, 1, 1];
    set stick_quality, 50;
    set sphere_quality, 4;
    color gray85, elem C;
    color red, elem O;
    color slate, elem N;
    color gray98, elem H;
    set stick_color, gray50;
    set ray_trace_mode, 1;
    set ray_texture, 2;
    set antialias, 3;
    set ambient, 0.5;
    set spec_count, 5;
    set shininess, 50;
    set specular, 1;
    set reflect, .1;
    set dash_gap, 0;
    set dash_color, black;
    set dash_gap, .15;
    set dash_length, .05;
    set dash_round_ends, 0;
    set dash_radius, .05;
    set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);
    preset.ball_and_stick("all",mode=1);
    draw 
    
    Commands without linebreaks:
    
    delete all;load 4dgr.pdb;fetch 4dgr2FoFc.mtz;select LongGlycan, resi 469:477;orient LongGlycan;remove not LongGlycan;remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;color density, 2fofcmap; show sticks;show spheres;set stick_radius, .07;set sphere_scale, .19;set sphere_scale, .13, elem H;set bg_rgb=[1, 1, 1];set stick_quality, 50;set sphere_quality, 4;color gray85, elem C;color red, elem O;color slate, elem N;color gray98, elem H;set stick_color, gray50;set ray_trace_mode, 1;set ray_texture, 2;set antialias, 3;set ambient, 0.5;set spec_count, 5;set shininess, 50;set specular, 1;set reflect, .1;set dash_gap, 0;set dash_color, black;set dash_gap, .15;set dash_length, .05;set dash_round_ends, 0;set dash_radius, .05;set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);preset.ball_and_stick("all",mode=1);draw 
 
    '''
    cmd.reinitialize()
    cmd.load('4dgr.pdb')
    cmd.load('4dgr2FoFc.ccp4')
    cmd.select('LongGlycan', 'resi 469:477')
    cmd.orient('LongGlycan')
    cmd.remove('not LongGlycan')
    cmd.remove('name H*')
    cmd.isomesh('2fofcmap','4dgr2FoFc', '1', 'LongGlycan', carve ='1.8')
    cmd.color('density','2fofcmap')
    cmd.show('sticks')
    cmd.show('spheres')
    cmd.set('stick_radius', '.07')
    cmd.set('sphere_scale', '.19')
    cmd.set('sphere_scale', '.13', 'elem H')
    cmd.set('bg_rgb', '[1, 1, 1]')
    cmd.set('stick_quality', '50')
    cmd.set('sphere_quality', '4')
    cmd.color('gray85', 'elem C')
    cmd.color('red', 'elem O')
    cmd.color('slate', 'elem N')
    cmd.color('gray98', 'elem H')
    cmd.set('stick_color', 'gray50')
    cmd.set('ray_trace_mode', '1')
    cmd.set('ray_texture', '2')
    cmd.set('antialias', '3')
    cmd.set('ambient', '0.5')
    cmd.set('spec_count', '5')
    cmd.set('shininess', '50')
    cmd.set('specular', '1')
    cmd.set('reflect', '.1')
    cmd.set('dash_gap', '0')
    cmd.set('dash_color', 'black')
    cmd.set('dash_gap', '.15')
    cmd.set('dash_length', '.05')
    cmd.set('dash_round_ends', '0')
    cmd.set('dash_radius', '.05')
    cmd.set_view('(0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0)')
    preset.ball_and_stick("all",mode=1);
    cmd.draw()
cmd.extend('LLG',LLG)


def LNA():
    '''
    DESCRIPTION

    Hydrated sodium cation bound in major groove of a 
    16-mer RNA of Watson-Crick base pairs.
    The sodium is bound to the N7 nitrogen atom of 
    Adenine 3 at 1.55 Angstrom resolution, PDB code 3nd4. 
    57 commands were used to make this figure. 

    More than one label in a horizontal script is not 
    allowed. This one label has to be at the end of the line.
    Labels can be imported from a Labels.pml file.
    Store the label commands one per row in this file.
    Import the file with the @Labels.pml command. 
    Include the path to the file if the labels file is not 
    in the current working directory of PyMOL. 
    
    USAGE

    Type 'NA' to execute. Type 'help NA' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.
        
    The commands with linebreaks:
    
    delete all;
    viewport 900,600;
    load 3nd4.pdb;
    run ~/mg18OU/quat.py; 
    quat 3nd4; 
    show sticks;
    set stick_radius=0.125;
    hide everything, name H*;
    bg_color white;
    create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));
    bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); 
    bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 519);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 319);
    distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);
    show nb_spheres; 
    set nb_spheres_size, .35;
    distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;
    distance hbond2,/3nd4_1/1/A/HOH`319/O,/3nd4_1/1/A/A`3/OP2;
    distance hbond3,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/HOH`119/O;
    distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;
    distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;
    distance hbond6,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    distance hbond7,/3nd4_1/1/A/HOH`319/O,/3nd4_1/1/A/G`2/OP2;
    distance hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;
    distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;
    distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;
    distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;
    distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;
    distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;
    distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;
    distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;
    distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;
    distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;
    distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;
    set depth_cue=0;
    set ray_trace_fog=0;
    set dash_color, black;
    set label_font_id, 5;
    set label_size, 36;
    set label_position, (0.5, 1.0, 2.0);
    set label_color, black;
    set dash_gap, 0.2;
    set dash_width, 2.0;
    set dash_length, 0.2;
    set label_color, black;
    set dash_gap, 0.2;
    set dash_width, 2.0;
    set dash_length, 0.2;
    select carbon, element C; 
    color yellow, carbon;
    disable carbon;
    set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0);
    rock

    The commands without linebreaks:
    
    delete all;viewport 900,600;load 3nd4.pdb;hide cartoon;run ~/mg18OU/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;rock;AOset_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 

    '''
    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.load('3nd4.pdb');
    cmd.hide('cartoon');
    cmd.do('run $HOME/mg18OU/quat.py')
    cmd.do('quat 3nd4');
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.set('sphere_scale', '0.225');
    cmd.hide('everything', 'name H*');
    cmd.bg_color('white');
    cmd.create('coorCov', '(3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)))');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/A`3/N7)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`119/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`219/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`319/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`419/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`519/O)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 519)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 419)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 119)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 319)');
    cmd.distance('(3nd4_1 and chain A and resi 19 and name NA)','(3nd4_1 and chain A and resi 219)');
    cmd.show('nb_spheres');
    cmd.set('nb_spheres_size', '.35');
    cmd.distance('hbond1', '/3nd4_1/1/A/HOH`119/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond2', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/A`3/OP2');
    cmd.distance('hbond3', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/HOH`119/O');
    cmd.distance('hbond4', '/3nd4_1/1/A/G`4/N7', '/3nd4_1/1/A/HOH`91/O');
    cmd.distance('hbond5', '/3nd4_1/1/A/G`4/O6', '/3nd4_1/1/A/HOH`419/O');
    cmd.distance('hbond6', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.distance('hbond7', '/3nd4_1/1/A/HOH`319/O', '/3nd4_1/1/A/G`2/OP2');
    cmd.distance('hbond9', '/3nd4_1/1/A/HOH`419/O', '/3nd4_2/2/A/HOH`74/O');
    cmd.distance('hbond10', '/3nd4_2/2/A/C`15/O2', '/3nd4_1/1/A/G`2/N2');
    cmd.distance('hbond11', '/3nd4_2/2/A/C`15/N3', '/3nd4_1/1/A/G`2/N1');
    cmd.distance('hbond12', '/3nd4_2/2/A/C`15/N4', '/3nd4_1/1/A/G`2/O6');
    cmd.distance('hbond13', '/3nd4_2/2/A/U`14/N3', '/3nd4_1/1/A/A`3/N1');
    cmd.distance('hbond14', '/3nd4_2/2/A/U`14/O4', '/3nd4_1/1/A/A`3/N6');
    cmd.distance('hbond15', '/3nd4_2/2/A/C`13/N4', '/3nd4_1/1/A/G`4/O6');
    cmd.distance('hbond16', '/3nd4_2/2/A/C`13/N3', '/3nd4_1/1/A/G`4/N1');
    cmd.distance('hbond17', '/3nd4_1/1/A/G`4/N2', '/3nd4_2/2/A/C`13/O2');
    cmd.distance('hbond18', '/3nd4_1/1/A/G`2/N2', '/3nd4_2/2/A/C`15/O2');
    cmd.distance('hbond19', '/3nd4_1/1/A/HOH`91/O', '/3nd4_1/1/A/G`4/OP2');
    cmd.set('depth_cue', '0');
    cmd.set('ray_trace_fog', '0');
    cmd.set('dash_color', 'black');
    cmd.set('label_font_id', '5');
    cmd.set('label_size', '36')
    cmd.set('label_position', '(0.5, 1.0,2.0)');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.set('label_color', 'black');
    cmd.set('dash_gap', '0.2');
    cmd.set('dash_width', '2.0');
    cmd.set('dash_length', '0.2');
    cmd.select('carbon', 'element C');
    cmd.color('yellow', 'carbon');
    cmd.disable('carbon');
    cmd.set_view('-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0');
    cmd.rock()
cmd.extend('LNA',LNA)


######################## Re-orient molecule ##############################

#category: Re-orient molecule

def oy():
    '''
    DESCRIPTION

    Align long axis of molecule along z-axis. 
    '''
    cmd.orient(); 
    cmd.turn('y',90) 
cmd.extend("oy",oy)


def omxy():
    '''
    DESCRIPTION

    Align long axis of molecule along minus x-y axis. 
    '''
    cmd.orient(); 
    cmd.turn('z',315) 
cmd.extend("omxy",omxy)


def oxy():
    '''
    DESCRIPTION

    Align long axis of molecule along x-y axis. 
    '''
    cmd.orient(); 
    cmd.turn('z',225) 
cmd.extend("oxy",oxy)


def oz():
    '''
    DESCRIPTION

    Align long axis of molecule along y-axis. 
    '''
    cmd.orient(); 
    cmd.turn('z',270) 
cmd.extend("oz",oz)


######################## Horizontal Scripting ##############################

#category: Horizontal scripting

def cntfiles():
    '''
    DESCRIPTION

    Count number of files in current directory. 

    USAGE
        
        cntpdb
    '''
    arg = (echo "Count the files in the directory.";echo "Usage: cntfiles.";find . -type f | wc -l)
    subprocess.call(arg,shell=True)
    return
cmd.extend('cntfiles',cntfiles)


def cntpdb():
    '''
    DESCRIPTION

    Count number of pdb files in current directory. 

    USAGE
        
        cntpdb
    '''
    arg = ("echo "Count the number of pdb files in subfolders."; echo "Usage: cntpdb"; \
        find ./ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e "{}/*_([1-9999]|****).pdb" ';' -print | wc -l")
    subprocess.call(arg,shell=True)
    return
cmd.extend('cntpdb',cntpdb)


def rline():
       '''
    DESCRIPTION

    Enter "help(rline)" to refresh memory of the readline commands.
      
    These commands are sufficient for most editing tasks:  
      To edit code, positon cursor on command line with left mouse button.  
      Control-e moves the cursor to the end of the line, even when it is out of view.
      Control-a moves the cursor to the beginning of the line, even when it is out of view.    
      Up arrow key recalls last line of commands for editing.

    These commands may not be available on all systems:
      Shift-control-a selects everything from the right of the cursor to the end of the line.
      Shift-control-e selects everything to the left of the cursor to the end of the line.
      Command-f moves the cursor to the end of the current word.  
      Command-b moves the cursor to the begining of the current word.
      Control-f moves the cursor to the right by one character.   
      Control-b moves the cursor to the left by one character.
        '''
cmd.extend("rline",rline)


def rv(StoredView=0, decimal_places=2, outname="roundedview.txt"):
    """
    DESCRIPTION
    
    Get the view settings in a compact format on one line.
            The default StoredView is "0", the current view.
            You can get a stored view assigned to some
            other digit with the view command) and rounds to two decimal
            places (two digits to the right of the decimal point) the
            18-element viewing matrix and rewrites the matrix elements
            on a single line with no whitespaces and a semicolon at the
            end. The saved space eases the making of a single line of
            PyMOL commands separated by semicolons. A semicolon 
            with nothing to the right of it at the end of a line of grouped commands
            is harmless.

        ARGUMENTS
            Usage: rv [view, decimal_places, outname] 
                Note that the values in the [] are optional.
            The default values  for the arguments of the function
            are "0,2, roundedview.txt". 

            Simple one-line example with roundview.py script in current working
            directory--check by typing 'pwd' and 'ls *.py' on the command line. PyMOL
            should return 'roundview.py' in the lisf of files in the external (top) gui.
            Next, paste the following command on the external (top) commandline, hit
            return, and wait 5-10 seconds:

        EXAMPLE
            fetch 1lw9, async=0; rv 0,2

        The following view setting will be returned without the blackslashes:

        set_view (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,-155.16,\
        35.13,11.48,9.72,122.33,187.99,-20.0);

        fetch 1lw9, async=0; rv 0,1

        set_view (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,-155.2,\
        35.1,11.5,9.7,122.3,188.0,-20.0);


            18 elements of the view settings (0-17)

            0 - 8 = column-major 3x3 matrix. Rotates the model axes
            to the camera axes. 

            9 - 11 = origin of rotation relative to the camera
            in camera space

            12 - 14 = origin of rotation in model space

            15 = front plane distance from the camera

            16 = rear plane distance from the camera

            17 = orthoscopic flag 
            (not implemented in older versions)

    """
    
    #convert the commandline arguments from strings to integers

    StoredView = int(StoredView)
    decimal_places = int(decimal_places)

    #call the get_view function

    m = cmd.get_view(StoredView)


    #Make a list of the elements in the orientation matrix.

    myList = [m[0], m[1], m[2], m[3], m[4], m[5], m[6],m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14],m[15], m[16], m[17]]

    #Round off the matrix elements to two decimal places (two fractional places)
    #This rounding approach solved the problem of unwanted
    #whitespaces when I tried using a string format statement

    myRoundedList = [round(elem, decimal_places) for elem in myList]
    
    #x is the string template for the output. The whitespace is required
    #between the "set_view" and "("

    x = 'set_view ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},\
{10},{11},{12},{13},{14},{15},{16},{17});'

    #print to the external gui.
    print x.format(*myRoundedList)

    #print 'set_view ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},\
    #{10},{11},{12},{13},{14},{15},{16},{17})'.format(*myRoundedList)

    #Write to a text file.
    myFile = open("roundedview.txt", "a")
    myFile.write(x.format(*myRoundedList) + "\n")
    myFile.close()
    return
    #The extend command makes roundview into a PyMOL command.
cmd.extend("rv", rv)

#################### Print commands for using git #################### 

#category: Print commands for using git.

def gitAdd():
       '''
    DESCRIPTION

    Enter help(gitAdd) to print steps for adding updates to a file under version control.

      Step 1: add new file to an existing repository. 

        git add fileToBeAdded 

    This command updates the index using the current content found in the working tree, 
    to prepare the content staged for the next commit.
    '''
cmd.extend("gitAdd",gitAdd)


def gitCommit():
       '''
    DESCRIPTION

    Enter help(gitInit) to print steps for saving updates to a file under version control.

      Step 1: commit changes to files 
        git -m commit "Message" file to be updated.  
    '''
cmd.extend("gitCommit",gitCommit)


def gitInit():
       '''
    DESCRIPTION

    Enter help(gitInit) to print steps for creating a git repository.
        Plain text files can be put under version control.
        Binary and bitmap files should not be put under version control.

      Step 1: initialize repository in current directory
        git init 

      Step 2: make list of binary file types to be ignored (e.g. *pse, *ccp4 *mtz, *png,
        *.dmg, *.pdf, *.tiff, *.jpeg, *.jpg, *.zip, *.tar, *.rar, *jar, *.iso, *.7z, *.o, 
        *.so, *.pyc, *.doc, *.docx, *.idtf, *.u3d, *.aux ...)


        touch .gitignore
        git add .gitignore
        git commit -m "message" .gitignore

      Step 3: add all other files to repository
        git add .
      Later add new files one at a time
        git add new.pml

      Step 4: commit new files or changes with message
        git commit -m "Edited new.pml." new.pml

        or to commit changes to a groupd of changed files
        git commit -a -m "Changed all *.pml files"

        To automate the process, add the following function to a .bashAliases file in
        your home directory and the source the .bashAliases file from your .bashrc file.
        Then in another terminal window change to the directory with the current pml
        file that you are editing and enter "gcpml new" to save your changes to "new.pml".
        Use the up arrow key to rerun this command as you make changes.  

        # Function to git commit changes to one pml file.
        # Takes the basename as a command line argument. 
        gcpml()
        {
        if [ $# -lt 1 ]; then
          echo 1>&2 "$0: not enough arguments"
          echo "Usage: gct baseOfTexFileName"
          echo "Note absence of file extension .pml"
          exit 2
        elif [ $# -gt 1 ]; then
          echo 1>&2 "$0: too many arguments"
          echo "Usage: gct baseOfTexFileName"
          echo "Note absence of file extension .pml"
          exit 2
        fi

        git add "$1".pml
        git commit -m "Added new text in $1.tex" "$1".pml
        }
        '''
cmd.extend("gitInit",gitInit)


def gitPull():
       '''
    DESCRIPTION

    Enter help(gitPush) to print steps to send to updates to a repository on github.com. 

      Step 1: pull file from an existing repository. 

        git pull 
    '''
cmd.extend("gitPull",gitPull)


def gitPush():
       '''
    DESCRIPTION

    Enter help(gitPush) to print steps to send to updates to a repository on github.com. 

      Step 1: push updated file to an existing repository. 

        git push 
    '''
cmd.extend("gitPush",gitPush)


############################## Send search tems to websites with search boxes #########################################

#category: Send search term(s) to websites with search boxes.

def AB(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to Amazon.com Books in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 

    USAGE

    AB search term(s), number of hits to returned

    EXAMPLE

    AB pymol
    '''
    url = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords='
    webbrowser.open(url+searchTerm)
cmd.extend('AB',AB)


# def AN(searchTerm="pymol"):
#     '''
#     DESCRIPTION

#     Send search term or phrase to anaconda.com books in default browser.
#     The search phrase does not need to be enclosed in quotes. 
#     The second argument is the number of hits to return. 
#     The default web browser is used. 

#     USAGE

#     AN search term(s), number of hits to returned

#     EXAMPLE

#     AB pymol
#     '''
#     url = ''
#     webbrowser.open(url+searchTerm)
# cmd.extend('AN',AN)


def GB(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to Amazon.com in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 

    USAGE

    GB search term(s), number of hits to returned

    EXAMPLE

    GB pymol
    '''
    url = 'https://www.google.com/search?tbm=bks&q='
    webbrowser.open(url+searchTerm)
cmd.extend('GB',GB)


def GH(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to GitHub in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 

    USAGE

    GH search term(s), number of hits to returned

    EXAMPLE

    GH pymol
    '''
    url = 'https://www.github.com/search?q='
    webbrowser.open(url+searchTerm)
cmd.extend('GH',GH)


def GHN(searchTerm="pymol", numHits=5):
    '''
    DESCRIPTION

    Send search term or phrase to GitHub in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 

    USAGE

    GHN search term(s), number of hits to returned

    EXAMPLE

    GHN pymol, 10
    '''
    print 'Searching Github...'  # display text while downloading the Google page
    url = 'https://www.github.com/search?q='
    res = requests.get(url + ' '.join(searchTerm))
    res.raise_for_status()
    soup = bs4.BeautifulSoup(res.text)
    linkElems = soup.select('.r a')
    numOpen = min(numHits, len(linkElems))
    for i in range(numOpen):
        t0 = time.time()
        webbrowser.open(url + linkElems[i].get('href'))
        response_delay = time.time() - t0
        time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
    print 'Finished searching Github.'  # display text while downloading the Google page
cmd.extend('GHN',GHN)


def GO(searchTerm="pymol",numHits="200"):
    '''
    DESCRIPTION

    Send search term or phrase Google in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 

    USAGE

    GO search term(s), number of hits to returned

    EXAMPLE

    GO Nobel Prize in Chemistry, 30
    '''
    webbrowser.open('https://www.google.com/search?q='+searchTerm+'&num='+str(numHits))
cmd.extend('GO',GO)


def GON(searchTerm="pymol",numHits="5"):
    '''
    DESCRIPTION

    Send search term or phrase Google in default browser and opens the top N results in N new tabs.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The second parameter is optional; its defult value is 5. 
    Each hit will be opened in a separate tab thereby saving a time consuming step.
    If the number of results is fewer than the number requested,
    all of the results will be shown.

    The default web browser is used. 

    Requires the Python modules requests and beautifulsoup4 (bs4).
    They may already be available to open source PyMOL, but they
    must be installed for the proprietary PyMOL. Use the following command
    from a terminal window outside of PyMOL:

    conda install requests beautifulsoup4

    You can launch this command from the commandline in PyMOL, but 
    the execution of the install can be slow and will tie up your
    PyMOL session for 10-20 minutes. 

    USAGE

    GON search term(s), number of hits to returned

    If the second argument is not given, the default value is used.

    EXAMPLE

    GON Nobel Prize in Chemistry, 7

    Prints message when the last page has finished loading.

    '''
    print 'Googling', searchTerm, 'and displaying the top', numHits, 'in separate tabs of the default brower.' 
    res = requests.get('http://google.com/search?q=' + ' '.join(searchTerm))
    res.raise_for_status()
    soup = bs4.BeautifulSoup(res.text)
    linkElems = soup.select('.r a')
    numOpen = min(numHits, len(linkElems))
    for i in range(numOpen):
        t0 = time.time()
        webbrowser.open('https://www.google.com' + linkElems[i].get('href'))
        response_delay = time.time() - t0
        time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
    print 'Finished googling ', searchTerm, '.' # display text while downloading the Google page
cmd.extend('GON',GON)


def GS(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to Google Scholar in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 
    The default search term is pymol.

    USAGE

    GS search term(s)

    EXAMPLES
    Single search phrase:

    GS Linus Pauling

    Multiple search terms:

    GS Linus Pauling; GS Francis Crick; GS Alexander Rich
    '''
    url = 'https://scholar.google.se/scholar?hl=en&q='
    webbrowser.open(url+searchTerm)
cmd.extend('GS',GS)


def GV(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to Google Videos in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 
    The default search term is pymol.

    USAGE

    GV search term(s)

    EXAMPLES
    Single search phrase:

    GV Linus Pauling

    Multiple search terms:

    GV Linus Pauling; GS Francis Crick; GS Alexander Rich
    '''
    url = 'https://www.google.com/search?q=video+'
    webbrowser.open(url+searchTerm)
cmd.extend('GV',GV)


def MA(searchTerm='pymol'):
    '''
    DESCRIPTION

    Send search term to all searchable websites in pymolshortcuts:

    arXiv
    bioRxiv
    GitHub
    Google
    Google Books
    Google Scholar
    Google Video
    PBD
    PubMed
    Pymol Mailing List
    Pymol Wiki
    Research Gate
    Science Direct
    Springer
    Source Forge
    Stackoverflow

    Example

    MA pymol plugin
    '''
    AX(searchTerm)
    BX(searchTerm)
    GB(searchTerm)
    GH(searchTerm)
    GO(searchTerm)
    GS(searchTerm)
    GV(searchTerm)
    PDB(searchTerm)
    PM(searchTerm)
    PML(searchTerm)
    PW(searchTerm)
    RG(searchTerm)
    SD(searchTerm)
    SF(searchTerm)
    SO(searchTerm)
    SP(searchTerm)
cmd.extend('MA',MA)


def MB(searchTerm='pymol'):
    '''
    DESCRIPTION

    Send search term to search multiple sites for term in books:

    Google Books
    Science Direct
    Springer

    Example

    MB pymol plugin
    '''
    GB(searchTerm)
    SD(searchTerm)
    SP(searchTerm)
cmd.extend('MB',MB)


def MC(searchTerm='pymol'):
    '''
    DESCRIPTION

    Send search term to search ten core websites in pymolshortcuts:

    bioRxiv
    GitHub
    Google
    Google Scholar
    PubMed
    Pymol Mailing List
    Pymol Wiki
    Research Gate
    Science Direct
    Stackoverflow

    Example

    MA pymol plugin
    '''
    BX(searchTerm)
    GH(searchTerm)
    GO(searchTerm)
    GS(searchTerm)
    PM(searchTerm)
    PML(searchTerm)
    PW(searchTerm)
    RG(searchTerm)
    SD(searchTerm)
    SF(searchTerm)
cmd.extend('MC',MC)


def MM(searchTerm='pymol'):
    '''
    DESCRIPTION

    Send search term to search for manuscripts in pymolshortcuts:

    arXiv
    bioRxiv
    Google Scholar
    PubMed
    Research Gate
    Science Direct
    Springer

    Example

    MM pymol plugin
    '''
    AX(searchTerm)
    BX(searchTerm)
    GS(searchTerm)
    PM(searchTerm)
    RG(searchTerm)
    SD(searchTerm)
    SP(searchTerm)
cmd.extend('MM',MM)


def PDB(searchTerm="3fa0"):
    '''
    DESCRIPTION
    
    Submit a search term to the Protein Data Bank.

    USAGE:
    PBB 3fa0
    '''
    webbrowser.open('https://www.rcsb.org/structure/'+searchTerm)
cmd.extend('PDB',PDB)


def PML(searchTerm="3d_pdf"):
    '''
    DESCRIPTION
    
    Submit a search term to the PyMOL Users Mail Service.

    USAGE:

    Single term search (multi word searches do NOT have to be inside quotes):
    PML session file

    Multiple term search: 
    PML text editor; PML 3d pdf; PML black and white cartoon;
    '''
    webbrowser.open('https://sourceforge.net/p/pymol/mailman/search/?q='+searchTerm)
cmd.extend('PML',PML)


def PM(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to PubMed.
    The default web browser is used.
    The multi word search terms do not need to be enclosed in quotes. 
    Takes one search term but multiple commands can be submitted at once (see below).

    USAGE
    PM search term

    EXAMPLES
    single
    PM molecular graphics

    Multiple search:
    PM molecular graphics;  PM molecular representation; PM ambient occlusion
    '''
    webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+searchTerm)
cmd.extend('PM',PM)


def IPM(searchTerms = [], *args):
    '''
    DESCRIPTION

    Read list of search terms and submit each term to PubMed in a separate browser tab.
    There a time delay based on the response time of the site to which the request is made.
    The default web browser is used.
    Must enclose each search term (can be of multiple words) in single or double quotes.
    Has a time delay to avoid overwhelming the webserver.

    USAGE

    search=[string,string]; IPM(search)

    Note that the name of the list is arbitrary.

    EXAMPLES

    search=["pymol","vmd","jmol"]; IPM(search)
    '''
    #termList = searchTerms.split(",")
    print 'Sending', searchTerms, 'to Pubmed and display list of search results in separate tabs of the default brower.'
    for term in searchTerms:
        t0 = time.time()
        sterm = str(term)
        webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+sterm)
        response_delay = time.time() - t0
        time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
        print 'Finished searching PubMed for', sterm, '.'
    print 'Finished searching PubMed for', searchTerms, '.' 
cmd.extend('IPM',IPM)


def IPMN(searchTerms = [], *args):
    '''
    DESCRIPTION

    Read list of search terms and submit each term to PubMed in a separate browser tab.
    There a time delay based on the response time of the site to which the request is made.
    The default web browser is used.
    Must enclose each search term (can be of multiple words) in single or double quotes.
    Has a time delay to avoid overwhelming the webserver.

    USAGE

    search=[string,string]; IPM(int,search)

    Note that the name of the list is arbitrary.

    EXAMPLES

    search=["pymol","vmd","jmol"]; IPMN(10,search)
    '''
    #termList = searchTerms.split(",")
    numHits="5"
    print 'Sending', searchTerms, 'to PubMed and display top N search results for each term in separate tabs of the default browser.' 
    for term in searchTerms:
        sterm = str(term)
        res = requests.get('https://www.ncbi.nlm.nih.gov/pubmed/?term=' + ' '.join(sterm))
        res.raise_for_status()
        soup = bs4.BeautifulSoup(res.text)
        linkElems = soup.select('.r a')
        numOpen = min(numHits, len(linkElems))
        for i in range(numOpen):
            t0 = time.time()
            webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term=' + linkElems[i].get('href'))
            response_delay = time.time() - t0
            time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
        print 'Finished searching PubMed for', sterm, '.'
    print 'Finished searching PubMed for', searchTerms, '.' 
cmd.extend('IPMN',IPMN)


def RG(searchTerm='best molecular graphics program'):
    '''
    DESCRIPTION
    
    Submit a search query of Research Gate. 

    Usage:

    RG best molecular graphics program
    '''
    webbrowser.open('https://www.researchgate.net/search.Search.html?type=researcher&query='+searchTerm)
cmd.extend('RG',RG)


def SD(searchTerm="pymol"):
    '''
    DESCRIPTION
    
    Submit a search term to Science Direct.

    USAGE:

    Single term search (multi word searches do NOT have to be inside quotes):
    SD session file

    Multiple term search: 
    SD text editor; SD 3d pdf; SD black and white cartoon;
    '''
    url1 = 'https://www.sciencedirect.com/search/advanced?qs='
    url2 = '&show=100&sortBy=relevance'
    webbrowser.open(url1+searchTerm+url2)
cmd.extend('SD',SD)


def SF(searchTerm='pymol'):
    '''
    DESCRIPTION
    
    Send search term to sourceforge.

    USAGE

    Single search:
    
    SF pymol

    Multiple search: 

    SF pymol; SF jmol; 
    '''
    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)
cmd.extend('SF',SF)


def SP(searchTerm="pymol"):
    '''
    DESCRIPTION
    
    Submit a search term to Springer Books

    USAGE:

    Single term search (multi word searches do NOT have to be inside quotes):
    SP session file

    Multiple term search: 
    SP text editor; SP 3d pdf; SP black and white cartoon;
    '''
    url1 = 'https://www.springer.com/gp/search?query='
    url2 = '&submit=Submit+Query'
    webbrowser.open(url1+searchTerm+url2)
cmd.extend('SP',SP)



################################ Open static web sites  ####################################

#category: Open static web sites

def ACA():
    '''
    DESCRIPTION

    Open the American Crystallographic Association Annual Meeting webpage.
    '''
    webbrowser.open_new_tab('http://www.amercrystalassn.org/2018-meeting-homepage')
cmd.extend('ACA',ACA)


def ALS():
    '''
    DESCRIPTION
    
    Open website of the Advanced Light Source.
    '''
    webbrowser.open('https://als.lbl.gov/')
cmd.extend('ALS',ALS)


def APS():
    '''
    DESCRIPTION
    
    Open website of the Advanced Photon Source.
    '''
    webbrowser.open('https://www.aps.anl.gov/')
cmd.extend('APS',APS)


def AX(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to arXiv.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 

    USAGE
    
    AX search term(s)

    EXAMPLE

    AX molecular graphics
    '''
    webbrowser.open('https://arxiv.org/search/?query='+searchTerm+'&searchtype=all&order=-announced_date_first&size=50')
cmd.extend('AX',AX)


def BC():
    '''
    DESCRIPTION
    
    Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.
    '''
    webbrowser.open('http://www.bio.aps.anl.gov/')
cmd.extend('BC',BC)


def BD():
    '''
    DESCRIPTION
    
    Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB). 
    '''
    webbrowser.open('https://www.sasbdb.org/')
cmd.extend('BD',BD)


def BX(searchTerm="pymol"):
    '''
    DESCRIPTION

    Send search term or phrase to bioRxiv 
    which is maintained by Cold Spring Harbor Laboratory.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 

    USAGE
    
    BX search term(s)

    EXAMPLES

    Single search:

    BX molecular graphics

    Multiple search:

    BX molecular graphics; BX pymol
    '''
    url = 'https://www.biorxiv.org/search/'
    webbrowser.open(url+searchTerm)
cmd.extend('BX',BX)


def CH():
    '''
    DESCRIPTION
    
    Open the webste of UCSF Chimera.
    '''
    webbrowser.open('https://www.cgl.ucsf.edu/chimera/')
cmd.extend('CH',CH)


def CHESS():
    '''
    DESCRIPTION
    
    Open the website of CHESS. 
    '''
    webbrowser.open('https://www.chess.cornell.edu/')
cmd.extend('CHESS',CHESS)


def EMDB():
    '''
    DESCRIPTION
    
    Open the website of the Electron Microscopy Data Bank.
    '''
    webbrowser.open('https://www.ebi.ac.uk/pdbe/emdb/')
cmd.extend('EMDB',EMDB)


def EP():
    '''
    DESCRIPTION
    
    EasyPyMOL github site.
    '''
    webbrowser.open('https://github.com/MooersLab/EasyPyMOL')


def JM():
    '''
    DESCRIPTION
    
    Open the Jmol wiki.
    '''
    webbrowser.open('http://wiki.jmol.org/index.php/Main_Page')
cmd.extend('JM',JM)


def IUCR(searchTerm="pymol"):
    '''
    DESCRIPTION

    Open website of the IUCr Journals.

    USAGE
    IUCR

    '''
    webbrowser.open('https://journals.iucr.org/')
cmd.extend('IUCR',IUCR)


def LBSF():
    '''
    DESCRIPTION
    
    Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.
    '''
    webbrowser.open('https://research.ouhsc.edu/CoreFacilities/LaboratoryofBiomolecularStructureandFunction.aspx')
cmd.extend('LBSF',LBSF)


def MCL():
    '''
    DESCRIPTION
    
    Open website of Macromolecular Crystallography Laboratory at the University of Oklahoma. 
    '''
    webbrowser.open('http://structuralbiology.ou.edu/mcl')
cmd.extend('MCL',MCL)


def MG():
    '''
    DESCRIPTION
    
    Open website of the OUHSC molecular graphics course.
    '''
    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/moleculargraphicslinks.html')
cmd.extend('MG',MG)


def NDB():
    '''
    DESCRIPTION
    
    Open website of the Nucleic Acid Database.
    '''
    webbrowser.open('http://ndbserver.rutgers.edu/')
cmd.extend('NDB',NDB)


def notPyMOL():
    '''
    DESCRIPTION
    
    Open website with list of other molecular graphics programs.
    '''
    webbrowser.open('https://en.wikipedia.org/wiki/List_of_molecular_graphics_systems')
cmd.extend('notPyMOL',notPyMOL)


def NSLSII():
    '''
    DESCRIPTION
    
    Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory.
    '''
    webbrowser.open('https://www.bnl.gov/ps/')
cmd.extend('NSLSII',NSLSII)


def PPC():
    '''
    DESCRIPTION
    
    Open the website of the Protein Production Facility at the University of Oklahoma in Norman.
    '''
    webbrowser.open('http://www.ou.edu/cas/chemistry/research/research-support-services/protein-production-core')
cmd.extend('PPC',PPC)


def PS():
    '''
    DESCRIPTION
    
    Open the home page of the Protein Soceity. 
    '''
    webbrowser.open('https://www.proteinsociety.org/')
cmd.extend('PS',PS)


def PW(searchTerm="3d_pdf"):
    '''
    DESCRIPTION
    
    Submit search of the PyMOL Wiki. 
    
    Usage:
    
    PW 3d_pdf
    '''
    webbrowser.open('https://pymolwiki.org/index.php/'+searchTerm)
cmd.extend('PW',PW)


def RS():
    '''
    DESCRIPTION
    
    Open the homepage of the RNA Society.
    '''
    webbrowser.open('https://www.rnasociety.org/')
cmd.extend('RS',RS)


def SAXS():
    '''
    DESCRIPTION
    
    Open the webpage of SAXS links at OUHSC. 
    
    '''
    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/small-angle-scattering-links-27aug2014.html?sfvrsn=0')
cmd.extend('SAXS',SAXS)


def SB():
    '''
    DESCRIPTION
    
    Open the webpage of SSRL Biological SAXS at BL 4-2.
    
    '''
    webbrowser.open('https://www-ssrl.slac.stanford.edu/~saxs/')
cmd.extend('SB',SB)


def SBGRID():
    '''
    DESCRIPTION
    
    Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel.
    
    '''
    webbrowser.open('https://www.youtube.com/user/SBGridTV/videos')
cmd.extend('SBGRID',SBGRID)


def SciPy18():
    '''
    DESCRIPTION
    
    Open the SciPy 2018 YouTube Channel.
    
    '''
    webbrowser.open('https://www.youtube.com/playlist?list=PLYx7XA2nY5Gd-tNhm79CNMe_qvi35PgUR')
cmd.extend('SciPy18',SciPy18)


def SSRL():
    '''
    DESCRIPTION
    
    Open the webpage of SSRL Structural Molecular Biology.
    '''
    webbrowser.open('http://ssrl.slac.stanford.edu/smb/index.html')
cmd.extend('SSRL',SSRL)

def SSURF():
    '''
    DESCRIPTION
    
    Open the webpage of the Society for Science at User Research Facilities (SSURF).
    SSURF is nonprofit organization in the US that serves as the
    nexus for users and user executive committees at
    national laboratories. SUURF is not a lobbying organization, but
    it help organize visits to Congress to educate legislators about
    the importance of national laboratories in science. Membership
    is free of students. The annual fee is nominal for PIs. 
    '''
    webbrowser.open('http://www.ssurf.org/')
cmd.extend('SSURF',SSURF)


def SO(searchTerm="3d_pdf"):
    '''
    DESCRIPTION
    
    Submit a search term to Stackoverflow.

    USAGE:

    Single term search (multi word searches do NOT have to be inside quotes):
    SO session file

    Multiple term search: 
    SO text editor; SO 3d pdf; SO black and white cartoon;
    '''
    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)
cmd.extend('SO',SO)




########################## 3D-PDFs #################################

#category: 3D-PDFs

def ms2pdf(InFile="pymol", caption="Wild-type T4 lysozyme, 3fa0, 1.09 Ang, Click in the center of the page to activate"):
    """
    DESCRIPTION
    
    Send molecular surface or ribbon cartoon from PyMOL to 3dpdf. 
    
    Does not ramp colors well. 
    Does not work on all other representations.
    Can instead export pse file, import the pse file into Jmol, 
    and enter in the console 'write fileNameStem.idtf'.
    Then enter the commands below in arg = on the command line outside of PyMOL.
    This works with stick models. 

    Usage:
        ms2pdf [,<fileNameStem> (do not include any file extensions)], [,caption]]

    Reads in an optional filename stem. The default is "pymol".
    Also reads in an optional caption.
    Saves an idtf file.
    Passes the idtf file to IDTFConverter, a u3d file to pdflatex,
            and opens the pdf with Adobe Acrobat Reader DC. 

    Prerequisites:
        Install u3d-1.4.5; may have to install zlib, libjpeg, and cairo (libpng)
        Install tex (on unix, texlive and texlive extra with the movie15 and media9 packages)
        Install Abode Acrobat Reader DC version 9+.
        See Jason Vertrees webpage for details (https://pymolwiki.org/index.php/3d_pdf)

        You can run in pymol "save inFile.idtf" to get new 3Droo and 3Dcoo parameters
        and then edit the corresponding lines below. 

    June 10, 2018
    Blaine Mooers


    """

    x = '''\documentclass[12pt,letter]{article}
\usepackage{hyperref}
\usepackage{media9}
\usepackage{verbatim}
\pagestyle{empty}
\\begin{document}
\\begin{figure}[!htb]
    \\begin{center}
        \\addmediapath{./} % here you can set the path where is been saved the u3d file
        \includemedia[
            label='''+ InFile +''',
            width=0.9\\textwidth,
            height=0.9\\textheight,
            activate=pageopen,
            deactivate=pageclose,
            3Dtoolbar=false,
            3Dnavpane=false,
            3Dmenu,
            3Droo=129.0,
            3Dcoo= 0.0 0.0 -129.0,
            3Dc2c=0.0 0.0 1.0,
            3Daac=20.0,
            3Droll=0.0,
            3Dbg=0 0 0, % to set the background color for 3D vwr; white = 1 1 1; so, you need to do the proportion: '255:1=[RGB]:x'
            transparent=false,
            3Dlights=Headlamp,
            3Drender=Solid,
            3Dpartsattrs=restore,
        ]{}{'''+ InFile +'''.u3d}
\\end{center}
\caption{''' + caption + '''.}
\\end{figure}
\\end{document}'''

    cmd.save("%s.idtf" %(InFile))
    FilePSE = str(InFile+".pse")
    FileIDTF = str(InFile+".idtf")
    FileU3D = str(InFile+".u3d")
    FileTEX = str(InFile+".tex")
    FilePDF = str(InFile+".pdf")
    print "InFile = ", InFile
    print "FileIDTF: ", FileIDTF
    print "FileU3D: ", FileU3D
    print "FileTEX: ", FileTEX
    print "FilePDF: ", FilePDF
    myFile = open(FileTEX, "a")
    myFile.write(x+"\n")
    myFile.close()

# Line continuation in python is done by wrapping the code in parentheses
    arg = ("/usr/local/bin/IDTFConverter -i " + FileIDTF + " -o " +
        FileU3D + "; /opt/local/bin/pdflatex " +
        InFile + "; open -a 'Adobe Acrobat Reader DC.app' " + FilePDF)
    subprocess.call(arg,shell=True)
    print "Argument for call: ", arg
    return
cmd.extend("ms2pdf", ms2pdf)


def topdf(InFile="pymol3"):
    """
    DESCRIPTION
    
    Send stick models as pse file from PyMOL through Jmol to 3DPDF.
    

    Usage:
        topdf [<filename> (do not include any file extensions)]

    Reads in a name for the command. The default is "pymol3".
    Saves the current session.
    Passes the session file to Jmol, IDTFConverter, pdflatex,
            and opens pdf with Adobe Reader. Edit the pymol3.tex

    If using PyMOL's Python, edit paths to the executables.
    Find executable paths with "which <executable>"
    Installtion Notes:
        Add to .bashrc file the path to your JMOL files:


        Copy jmol.sh to jmoldata.sh and edit to replace "Jmol.jar" with 
        "JmolData.jar" on command line at the bottom of the file. 

        JmolData.jar is a windowless Jmol for extracting data. 
        Install u3d-1.4.5; may have to install zlib, libjpeg, and cairo (libpng)
        Install texlive and texlive extra with the movie15 and media9 packages
        Install Abode Acrobat Reader DC version 9+

    Reset the camera parameters in the *.tex by finding the ones for the
    default view by selecting Default View in Adobe Reader and then "Get Current
    View". Wait for another window to open. Scroll down to the camera settings.
    Copying and paste these settings into the *.tex document.

    Worked on October 16, 2014; busted June 10, 2018
    Blaine Mooers

    new feature: Jmol reads PyMOL 1.8 PSE files with "set dump_binary, 1"


    """
    cmd.set("pse_export_version", "1.74")
    cmd.save("%s.pse" %(InFile))
    FilePSE = str(InFile+".pse")
    FileIDTF = str(InFile+".idtf")
    FileU3D = str(InFile+".u3d")
    FileTEX = str(InFile+".tex")
    FilePDF = str(InFile+".pdf")
    print "InFile = ", InFile
    print "FilePSE: ", FilePSE
    print "FileIDTF: ", FileIDTF
    print "FileU3D: ", FileU3D
    print "FileTEX: ", FileTEX
    print "FilePDF: ", FilePDF
# Line continuation in python is done by wrapping the code in parentheses
    arg = ("sh /Applications/jars/jmol-14.29.16/jmoldata.sh -nj" +
        " 'load " + FilePSE + ";write " + FileIDTF +
        "';/usr/local/bin/IDTFConverter -i " + FileIDTF + " -o " +
        FileU3D + "; /opt/local/bin/pdflatex " +
        FileTEX + "; open -a 'Adobe Acrobat Reader DC.app' " + FilePDF)
    subprocess.call(arg,shell=True)
    print "Argument for call: ", arg
    return
cmd.extend("topdf", topdf)


""" Print the shortcuts on startup of PyMOL"""
print(SC.__doc__)
