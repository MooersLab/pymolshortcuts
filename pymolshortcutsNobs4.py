# -*- coding: utf-8 -*-
'''
DESCRIPTION

    Loads functions that enable faster work in PyMOL.

    When sourced in the pymolrc file by including the line 'run pymolshortcuts.py,
    this script defines a number of Python functions upon the startup of PyMOL.
    These functions are called shortcuts becasue they are invoked with
    shortnames in analogy to keyborad shortcuts. Unlike keyboard shortcuts, 
    their names consist only of the home keys for faster entry.
    The shortcuts are designed to save time and address some shortcommings
    of PyMOL that are annoying like no version control for output files.
    
    The variant of this file called pymolshortcutsNobs4.py does not depend 
    on Python module beuafitulsoup4 which does not come with PyMOL.
    


INSTALLATION 
    
    Source this file from your .pymolrc file on the mac or linux or 
    from your pymolrc.pml file on Windows by adding the command 
    below (adjust the file path as appropriate) on a separate line:
                           
    run ~/Scripts/PyMOLscripts/pymolshortcuts.py 

    Enter "SC" to print to the command history window a list of shortcuts. 
    
    You will need to edit the paths to your installation of several external applications.
    See the section "PATHS to Applications" below.
    
    
REQUIREMENTS 
    
    Should work for PyMOL running with Python 2.7 and 3.7 interpreters.
    Incentive PyMOL built with Anaconda Python3.7 is ideal but not essential. 
    
    
    To use the functions that submit searches to multiple websites from within PyMOL.
    install the following packages: 


    bs4>=4.8.0
    requests>=2.22.0

    Note that bs4 == beautifulsoup4
    
    With the incentive PyMOL that uses the Anaconda Python as its Python interpreter, 
    enter the following command on the upper command line in PyMOL:
    
    conda install beautifulsoup4 requests
    
    
    For the Open Source PyMOL, use the software manager for the Python used to 
    build PyMOL. E.g., with macports, use "sudo port install py37-bs4 &&
    sudo port install py37-requests"
    
    Alternatively, use pip with the python used to install PyMOL:
    
    /opt/local/bin/python3.7 -m pip install --upgrade bs4
    /opt/local/bin/python3.7 -m pip install --upgrade requests
    
    
    The above procedure should work for older versions of Open Source PyMOL.
    

    PLAN B -- BS4FREE VERSION

    For older incentrive versions of PyMOL, the installation of external modules is not easy
    and it may simpler to forgo the functions that depend on beautiful soup. 
    For these versions of PyMOL, download pymolshortcutsnobs4.py. 
    This script lacks the bs4 and request dependent functions. 

    Source this file from your .pymolrc file on the mac or linux or 
    from your pymolrc.pml file on Windows by adding the command 
    below (adjust the file path as appropriate) on a separate line:
    
    
    run ~/Scripts/PyMOLscripts/pymolshortcutsNobs4.py 
    
    
HELP
    
    To see the documentation for shortcut, enter at the PyMOL prompt:
    
    
    help <ShortCutName>
    
    e.g.,
    
    help AO

    The documentation for a shortcut includes a 1) description of the shortcut, 
    2) one or more examples of its usage, 3) the corresponding PyMOL Macro 
    Language codeas vertical script, 4) the same as a horizontal script to ease 
    copying and pasting the code on the command line, 
    and 5) the corresponding Python code. 
 

 Copyright Notice
 ================
  
 Copyright (c) 2019 Board of Regents for the University of Oklahoma

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 
   
  Blaine Mooers , PhD 
  blaine@ouhsc.edu
  975 NE 10th St, BRC 466
  Department of Biochemistry and Molecular Biology
  College of Medicine
  Stephenson Cancer Center
  University of Oklahoma Health Sciences Center, 
  Oklahoma City, OK, USA 73104

''' 
from __future__ import division, print_function

import subprocess
from math import cos, sin, radians, sqrt
import datetime, time, webbrowser, random

from pymol import cmd, stored, cgo, xray
import numpy


__author__ = "Blaine Mooers"
__copyright__ = "2019 Board of Regents for the University of Oklahoma"
__license__ = "MIT Licencse"
__version__ = "1.0.0" 
# Versioning follows follows MAJOR.MINOR[.PATCH] where major releases are not backward compatable. 
__credits__ = ["Miriam Shakir", "Safra Shakir"] 
# Credits are for people who tested the code, reported bug fixes, made suggestions, etc. 
__date__ = "10 Septermber 2019"
__maintainer__ = "Blaine Mooers"
__email__ = "blaine-mooers@ouhsc.edu"
__status__ = "Production" 

cmd.set('ray_opaque_background','on')

AppPaths='''You may have to edit the file paths to your applications around line XXX in pymolshortcut.py.'''
print(AppPaths)

##################################  PATHS to Applications ###############################################
#
# Set the file paths to the applications on your computer.
# Comment out the paths to the applications that you lack. 
#
# The paths are for my Mac OS 10.14 wiht some applications stored in /usr/local 
# and some in /opt/local/bin thanks to macports.
# Some paths can be found by entering "which <programname>" in a terminal window.
#
# You can replace the file path with a command to open the application.
# For example, you can use the command "open -a bbedit" on the Mac. 
#
##################################  PATHS to Applications ###############################################

# Data analysis
DBBrowserSQLitePath = '/Applications/DBBrowserForSQLite.app/Contents/MacOS/DB\ Browser\ for\ SQLite'
excelPath = '/Applications/Microsoft\ Excel.app/Contents/MacOS/Microsoft\ Excel'
jabrefPath = '/Applications/JabRef.app/Contents/MacOS/JavaApplicationStub'
jaspPath = 'Applications/JASP.app/Contents/MacOS/JASP'
jmpPath = '/Applications/JMP\ Pro\ 14.app/Contents/MacOS/JMP'
juliaPath = '/Applications/Julia-1.2.app/Contents/MacOS/applet'
jupyterPath = '/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter-notebook'
RstudioPath = '/Applications/RStudio.app/Contents/MacOS/Rstudio'


# Local files 
localPDBfilePath = '$HOME/pdbFiles/'
localEMAPfilePath = '$HOME/emapFiles/'
localHKLfilePath = '$HOME/hklFiles/'

# Image manipulation programs
gimp = '/usr/local/bin/gimp'
inkscape = '/opt/local/bin/inkscape'
pptPath =  '/Applications/Microsoft\ Excel.app/Contents/MacOS/Microsoft\ PowerPoint'


# Molecular graphics programs
ccp4mgPath = '/Applications/ccp4-7.0/ccp4i2.app/Contents/MacOS/ccp4mg'
chimeraPath = '/Applications/Chimera.app/Contents/MacOS/chimera'
cootPath = '/usr/local/bin/coot'
jmolPath = 'java -jar /Applications/jars/jmol-14.29.52/Jmol.jar'
yasaraPath = '/Applications/YASARA.app/Contents/MacOS/yasara.app'
vmdPath = '/Applications/VMD194.app/Contents/MacOS/startup.command'

# Text editors
atomPath = '/usr/local/bin/atom'
bbeditPath = '/usr/local/bin/bbedit'
codePath = '/usr/local/bin/code'
emacsPath = '/opt/local/bin/emacs'
geditPath = '/Applications/gedit2.30.2.app/Contents/MacOS/gedit'
jeditPath = '/Applications/jEdit.app/Contents/MacOS/jedit'
notepadppPath = '/Applications/Notepad++.app/Contents/MacOS/startwine'
neovimPath = '/Users/blaine/software/nvim-osx64/bin/nvim'
oniPath = '/Applications/Oni.app/Contents/MacOS/Oni'
pdbedPath = 'java -jar /Applications/jars/PDB_Editor_FIX090203.jar'
sublimeText3Path = '/usr/local/bin/subl'
textmatePath  = '/usr/local/bin/mate'
vimPath = '/opt/local/bin/vim'

# Web sites 
gcalPath = 'open -a Safari.app https://calendar.google.com/calendar/r'
gmailPath = 'open -a Safari.app https://mail.google.com/mail/u/0/#inbox'
webmailPath = 'open -a Safari.app https://webmail.ouhsc.edu'
weatherServicePath = 'open -a Safari.app https://radar.weather.gov/radar.php?rid=TLX'   
 
########################################################################################################
def AB(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to Amazon.com Books in default browser.

    USAGE:

    AB

    Arguments:

    searchTerm

    EXAMPLE:

    AB

    MORE DETAILS:

    Send search term or phrase to Amazon.com Books in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 


    VERTICAL PML SCRIPT:

    url = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords=';
    webbrowser.open(url+searchTerm)


    HORIZONTAL PML SCRIPT:

    url = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords=';webbrowser.open(url+searchTerm)


    PYTHON CODE:

def AB(searchTerm="pymol"):
    url = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords='
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords='
    webbrowser.open(url+searchTerm)

cmd.extend('AB',AB)

def ACA():
    ''' 
    DESCRIPTION:

    Open the American Crystallographic Association Annual Meeting webpage.

    USAGE:

    ACA

    Arguments:

    NA

    EXAMPLE:

    ACA

    MORE DETAILS:

    Open the American Crystallographic Association Annual Meeting webpage.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def ACA():
    webbrowser.open_new_tab('http://www.amercrystalassn.org/2018-meeting-homepage')

    '''

    webbrowser.open_new_tab('http://www.amercrystalassn.org/2018-meeting-homepage')

cmd.extend('ACA',ACA)

def ALS():
    ''' 
    DESCRIPTION:

    Open website of the Advanced Light Source.

    USAGE:

    ALS

    Arguments:

    NA

    EXAMPLE:

    ALS

    MORE DETAILS:

    Open website of the Advanced Light Source.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def ALS():
    webbrowser.open('https://als.lbl.gov/')
    '''

    webbrowser.open('https://als.lbl.gov/')
cmd.extend('ALS',ALS)

def AO():
    ''' 
    DESCRIPTION:

    Commands to make ambient occlusion image like those in Qutemole.

    USAGE:

    AO

    Arguments:

    None

    EXAMPLE:

    AO

    MORE DETAILS:

    Commands to make ambient occlusion image like those in Qutemole.
    Works only with the command line immediately under the command 
    history window at the top of the gui.

    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    set_color oxygen, [1.0,0.4,0.4];set_color nitrogen, [0.5,0.5,1.0];remove solvent;as spheres;util.cbaw;bg white;set light_count,8;set spec_count,1;set shininess, 10;set specular,0.25;set ambient,0;set direct,0;set reflect,1.5;set ray_shadow_decay_factor, 0.1;set ray_shadow_decay_range, 2;set depth_cue,0;color gray20, symbol c;ray

    PYTHON CODE:

def AO():
    cmd.set_color("oxygen", "[1.0,0.4,0.4]")
    cmd.set_color("nitrogen", "[0.5,0.5,1.0]")
    cmd.remove("solvent")
    cmd.show_as("spheres")
    cmd.util.cbaw()
    cmd.bg_color("white")
    cmd.set("light_count", "8")
    cmd.set("spec_count", "1")
    cmd.set("shininess", "10")
    cmd.set("specular", "0.25")
    cmd.set("ambient", "0")
    cmd.set("direct", "0")
    cmd.set("reflect", "1.5")
    cmd.set("ray_shadow_decay_factor", "0.1")
    cmd.set("ray_shadow_decay_range", "2")
    cmd.set("depth_cue","0")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
    '''

    cmd.set_color("oxygen", "[1.0,0.4,0.4]")
    cmd.set_color("nitrogen", "[0.5,0.5,1.0]")
    cmd.remove("solvent")
    cmd.show_as("spheres")
    cmd.util.cbaw()
    cmd.bg_color("white")
    cmd.set("light_count", "8")
    cmd.set("spec_count", "1")
    cmd.set("shininess", "10")
    cmd.set("specular", "0.25")
    cmd.set("ambient", "0")
    cmd.set("direct", "0")
    cmd.set("reflect", "1.5")
    cmd.set("ray_shadow_decay_factor", "0.1")
    cmd.set("ray_shadow_decay_range", "2")
    cmd.set("depth_cue","0")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
cmd.extend('AO',AO)

def AOD():
    ''' 
    DESCRIPTION:

    Make ambient occlusion image of any with dark carbon atoms.

    USAGE:

    AOD

    Arguments:

    None

    EXAMPLE:

    AOD

    MORE DETAILS:

    Type "help AOD" to see this documentation printed to the command history window. 
    Select from the command history individual lines of code to build a new script. 
    Select the hortizontal script at the bottom if retaining most of the commands 
    in your new script. Copy and paste onto the comand line below. Works only 
    with the command line immediately under the command history window at 
    the top of the gui.

    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    set_color oxygen, [1.0,0.4,0.4];set_color nitrogen, [0.5,0.5,1.0];remove solvent;as spheres;util.cbaw;bg white;set light_count,8;set spec_count,1;set shininess, 10;set specular,0.25;set ambient,0;set direct,0;set reflect,1.5;set ray_shadow_decay_factor, 0.1;set ray_shadow_decay_range, 2;set depth_cue,0;color gray20, symbol c;color gray70, symbol h;ray

    PYTHON CODE:

def AOD():
    cmd.set_color("oxygen", "[1.0,0.4,0.4]")
    cmd.set_color("nitrogen", "[0.5,0.5,1.0]")
    cmd.remove("solvent")
    cmd.show_as("spheres")
    cmd.util.cbaw()
    cmd.set_color("carbon", "[0.00 , 0.00 , 0.0]")
    cmd.bg_color("white")
    cmd.set("light_count","8")
    cmd.set("spec_count", "1")
    cmd.set("shininess", "10")
    cmd.set("specular", "0.25")
    cmd.set("ambient", "0")
    cmd.set("direct", "0")
    cmd.set("reflect", "1.5")
    cmd.set("ray_shadow_decay_factor", "0.1")
    cmd.set("ray_shadow_decay_range", "2")
    cmd.set("depth_cue","0")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
    '''

    cmd.set_color("oxygen", "[1.0,0.4,0.4]")
    cmd.set_color("nitrogen", "[0.5,0.5,1.0]")
    cmd.remove("solvent")
    cmd.show_as("spheres")
    cmd.util.cbaw()
    cmd.set_color("carbon", "[0.00 , 0.00 , 0.0]")
    cmd.bg_color("white")
    cmd.set("light_count","8")
    cmd.set("spec_count", "1")
    cmd.set("shininess", "10")
    cmd.set("specular", "0.25")
    cmd.set("ambient", "0")
    cmd.set("direct", "0")
    cmd.set("reflect", "1.5")
    cmd.set("ray_shadow_decay_factor", "0.1")
    cmd.set("ray_shadow_decay_range", "2")
    cmd.set("depth_cue","0")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
cmd.extend("AOD",AOD)

def APS():
    ''' 
    DESCRIPTION:

    Open website of the Advanced Photon Source.

    USAGE:

    APS

    Arguments:

    NA

    EXAMPLE:

    APS

    MORE DETAILS:

    Open website of the Advanced Photon Source.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def APS():
    webbrowser.open('https://www.aps.anl.gov/')
    '''

    webbrowser.open('https://www.aps.anl.gov/')
cmd.extend('APS',APS)

def AX(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to arXiv.

    USAGE:

    AX

    Arguments:

    searchTerm

    EXAMPLE:

    AX

    MORE DETAILS:

    Send search term or phrase to arXiv.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def AX(searchTerm="pymol"):
    webbrowser.open('https://arxiv.org/search/?query='+searchTerm+'&searchtype=all&order=-announced_date_first&size=50')

    '''

    webbrowser.open('https://arxiv.org/search/?query='+searchTerm+'&searchtype=all&order=-announced_date_first&size=50')

cmd.extend('AX',AX)

def BC():
    ''' 
    DESCRIPTION:

    Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.


    USAGE:

    BC

    Arguments:

    NA

    EXAMPLE:

    BC

    MORE DETAILS:

    Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def BC():
    webbrowser.open('http://www.bio.aps.anl.gov/')
    '''

    webbrowser.open('http://www.bio.aps.anl.gov/')
cmd.extend('BC',BC)

def BD():
    ''' 
    DESCRIPTION:

    Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB). 

    USAGE:

    BD

    Arguments:

    NA

    EXAMPLE:

    BD

    MORE DETAILS:

    Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB). 

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def BD():
    webbrowser.open('https://www.sasbdb.org/')
    '''

    webbrowser.open('https://www.sasbdb.org/')
cmd.extend('BD',BD)

def BST():
    ''' 
    DESCRIPTION:

    G2G3/U9U8 base step , PDB code 4PCO. 

    USAGE:

    BST

    Arguments:

    None

    EXAMPLE:

    BST

    MORE DETAILS:

    G2G3/U9U8 base step , PDB code 4PCO. 
    From the 1.32 Angstrom resolution structure 
    of the RNA decamer with 8 GU base pairs.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 4PCO, type=pdb, async=0;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png

    PYTHON CODE:

def BST():
    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb')
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

    '''

    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb')
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

def BU():
    ''' 
    DESCRIPTION:

    Commands to make biological unit.

    USAGE:

    BU

    Arguments:

    None

    EXAMPLE:

    BU

    MORE DETAILS:

    Commands to make biological unit. Requires a pdb file. There are
    other ways of displaying the biological unit in PyMOL. Depends on
    the quat3.py script by Thomas Holder.

   >>>  Edit file path in python code below

    Type 'help BU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    VERTICAL PML SCRIPT:

    run ~/Scripts/PyMOLScripts/quat3.py; 
    quat 


    HORIZONTAL PML SCRIPT:

    run ~/Scripts/PyMOLScripts/quat3.py; quat 

    PYTHON CODE:

def BU():
    cmd.do('run $HOME/Scripts/PyMOLScripts/quat3.py')
    cmd.do('quat')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/quat3.py')
    cmd.do('quat')
cmd.extend('BU',BU)

def BW():
    ''' 
    DESCRIPTION:

    Make black-and white-ribbon cartoon on a white background.

    USAGE:

    Orient struture as desired. Then type 'BW' to execute the function. Type
    'help BW' to see this documentation printed to the command history window.
    Select from the command history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    Arguments:

    NA

    EXAMPLE:

    BW

    MORE DETAILS:

    Commands to make black-and white-ribbon cartoon on a white background.

    VERTICAL PML SCRIPT:

    show cartoon; 
    hide lines; 
    hide nonbonded; 
    # black and white cartoon;
    # note how the dcomment is on a separate line and not to the right of a command; 
    set ray_trace_mode, 2; 
    bg_color white; 
    set opaque_background, on
    # change to off if you are making cover artwork
    set antialias, 2; 
    ray 1600,1600; 
    png test.png


    HORIZONTAL PML SCRIPT:

    show cartoon; hide lines; hide nonbonded; set ray_trace_mode, 2; # black and white cartoon; bg_color white; set antialias, 2; ray 1600,1600; png test.png

    PYTHON CODE:

def BW():
    cmd.show_as("cartoon", "all"); 
    cmd.hide('lines'); 
    cmd.hide('nonbonded'); 
    # black and white cartoon; 
    cms.set("opaque_background", "on")
    cmd.set('ray_trace_mode', '2'); 
    cmd.bg_color('white'); 
    cmd.set('antialias', '2'); 
    cmd.ray('600','600'); 
    cmd.png('test.png')
    '''

    cmd.show_as("cartoon", "all"); 
    cmd.hide('lines'); 
    cmd.hide('nonbonded'); 
    # black and white cartoon; 
    cms.set("opaque_background", "on")
    cmd.set('ray_trace_mode', '2'); 
    cmd.bg_color('white'); 
    cmd.set('antialias', '2'); 
    cmd.ray('600','600'); 
    cmd.png('test.png')
cmd.extend('BW', BW)

def BX(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to bioRxiv 

    USAGE:

    BX

    Arguments:

    searchTerm

    EXAMPLE:

    BX

    MORE DETAILS:

    Send search term or phrase to bioRxiv 
    which is maintained by Cold Spring Harbor Laboratory.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def BX(searchTerm="pymol"):
    url = 'https://www.biorxiv.org/search/'
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://www.biorxiv.org/search/'
    webbrowser.open(url+searchTerm)

cmd.extend('BX',BX)

def CB():
    ''' 
    DESCRIPTION:

    Runs Jared Sampson's script "colorblindfriendly.py". 

    USAGE:

    CB

    Arguments:

    None

    EXAMPLE:

    CB
color cb_sky_blue, 3fa0

    MORE DETAILS:

    Runs Jared Sampson's script "colorblindfriendly.py" which defines colorblind 8 friendly colors.
	
    This shortcut only loads the colors. 
    You still need to apply the colors.
    See the example below. 

    More information and examples can be found at:
    http://www.pymolwiki.org/index.php/color_blind_friendly

    DESCRIPTION

        Certain colors are indistinguishable to people with the various forms of
        color blindness, and therefore are better not used in figures intended for
        public viewing.

        This script generates a palette of named colors for PyMOL that are
        unambiguous both to colorblind and non-colorblind people.

        The colors listed here are defined according to recommendations found at
        http://jfly.iam.u-tokyo.ac.jp/html/color_blind/.  This website is a good
        reference to consult when making all kinds of figures, not just those made
        using PyMOL.

        The colors are:

        * cb_black
        * cb_orange
        * cb_sky_blue (also: cb_skyblue, cb_light_blue, cb_lightblue)
        * cb_bluish_green (also: cb_bluishgreen, cb_green)
        * cb_yellow
        * cb_blue
        * cb_vermillion (also: cb_red, cb_redorange, cb_red_orange)
        * cb_reddish_purple (also: cb_rose, cb_violet, cb_magenta)

    USAGE

        CB
        color myObject, cb_red
        color mySel, cb_yellow

    REQUIREMENTS

        None.

    AUTHOR

        Jared Sampson, NYU Langone Medical Center, 2014

    LICENSE

    Copyright (c) 2014 Jared Sampson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.

    __author__ = 'Jared Sampson'
    __version__ = '0.1'


    VERTICAL PML SCRIPT:

    run ~/Pymol-script-repo/colorblindfriendly.py

    HORIZONTAL PML SCRIPT:

    run ~/Pymol-script-repo/colorblindfriendly.py 

    PYTHON CODE:

def CB():
    # Color blind friendly color list based on information found at:
    # http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
    # The RGB percentage values given on that page are less precise than the 0-255
    # values, so the 0-255 values are converted here (e.g. 230/255 = 0.902).
    cb_colors = (
        ("black", (0.000, 0.000, 0.000),                   # (  0,   0,   0)
                ()),
        ("orange", (0.902, 0.624, 0.000),                   # (230, 159,   0)
         ()),
        ("sky_blue", (0.337, 0.706, 0.914),                   # ( 86, 180, 233)
         ("skyblue", "light_blue", "lightblue")),
        ("bluish_green", (0.000, 0.620, 0.451),                   # (  0, 158, 115)
         ("bluishgreen", "green")),
        ("yellow", (0.941, 0.894, 0.259),                   # (240, 228,  66)
         ()),
        ("blue", (0.000, 0.447, 0.698),                   # (  0, 114, 178)
         ()),
        ("vermillion", (0.835, 0.369, 0.000),                   # (213,  94,   0)
         ("red", "red_orange", "redorange")),
        ("reddish_purple", (0.800, 0.475, 0.655),                   # (204, 121, 167)
         ("reddishpurple", "rose", "violet", "magenta")),  
    )

    for c in cb_colors:
        # main name
        cmd.set_color("cb_%s" % c[0], c[1])
        print("Set color: cb_%s" % c[0])

        # alternate names
        for alt in c[2]:
            cmd.set_color("cb_%s" % alt, c[1])
            print("           cb_%s" % alt)

    '''

    # Color blind friendly color list based on information found at:
    # http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
    # The RGB percentage values given on that page are less precise than the 0-255
    # values, so the 0-255 values are converted here (e.g. 230/255 = 0.902).
    cb_colors = (
        ("black", (0.000, 0.000, 0.000),                   # (  0,   0,   0)
                ()),
        ("orange", (0.902, 0.624, 0.000),                   # (230, 159,   0)
         ()),
        ("sky_blue", (0.337, 0.706, 0.914),                   # ( 86, 180, 233)
         ("skyblue", "light_blue", "lightblue")),
        ("bluish_green", (0.000, 0.620, 0.451),                   # (  0, 158, 115)
         ("bluishgreen", "green")),
        ("yellow", (0.941, 0.894, 0.259),                   # (240, 228,  66)
         ()),
        ("blue", (0.000, 0.447, 0.698),                   # (  0, 114, 178)
         ()),
        ("vermillion", (0.835, 0.369, 0.000),                   # (213,  94,   0)
         ("red", "red_orange", "redorange")),
        ("reddish_purple", (0.800, 0.475, 0.655),                   # (204, 121, 167)
         ("reddishpurple", "rose", "violet", "magenta")),  
    )

    for c in cb_colors:
        # main name
        cmd.set_color("cb_%s" % c[0], c[1])
        print("Set color: cb_%s" % c[0])

        # alternate names
        for alt in c[2]:
            cmd.set_color("cb_%s" % alt, c[1])
            print("           cb_%s" % alt)

cmd.extend('CB',CB)

def CBSS():
    ''' 
    DESCRIPTION:

    Apply colorblind-friendly coloring to ribbon or cartoon representations.

    USAGE:

    CBSS

    Arguments:

    None

    EXAMPLE:

    CBSS

    MORE DETAILS:

    Apply colorblind-friendly coloring to ribbon or cartoon representations.
    Uses Jared Sampson's colorblind friendly colors. 
    See the function CB.
    


    VERTICAL PML SCRIPT:

    CB;
    as cartoon;
    color cb_red, ss H;
    color cb_yellow,ss S;
    color cb_green, ss L+; 


    HORIZONTAL PML SCRIPT:

    CB;as cartoon;color cb_red, ss H;color cb_yellow,ss S;color cb_green, ss L+; 

    PYTHON CODE:

def CBSS():
    cmd.do('CB')
    cmd.show_as('cartoon')
    cmd.color('cb_red', 'ss H')
    cmd.color('cb_yellow', 'ss S')
    cmd.color('cb_green', 'ss L+')

    '''

    cmd.do('CB')
    cmd.show_as('cartoon')
    cmd.color('cb_red', 'ss H')
    cmd.color('cb_yellow', 'ss S')
    cmd.color('cb_green', 'ss L+')

cmd.extend('CBSS',CBSS)

def CH():
    ''' 
    DESCRIPTION:

    Open the webste of UCSF Chimera.

    USAGE:

    CH

    Arguments:

    NA

    EXAMPLE:

    CH

    MORE DETAILS:

    Open the webste of UCSF Chimera.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def CH():
    webbrowser.open('https://www.cgl.ucsf.edu/chimera/')
    '''

    webbrowser.open('https://www.cgl.ucsf.edu/chimera/')
cmd.extend('CH',CH)

def CHESS():
    ''' 
    DESCRIPTION:

    Open the website of CHESS.

    USAGE:

    CHESS

    Arguments:

    NA

    EXAMPLE:

    CHESS

    MORE DETAILS:

    Open the website of CHESS.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def CHESS():
    webbrowser.open('https://www.chess.cornell.edu/')

    '''

    webbrowser.open('https://www.chess.cornell.edu/')

cmd.extend('CHESS',CHESS)

def CR():
    ''' 
    DESCRIPTION:

    Commands to make colored filled-ring cartoon of nucleic acids.

    USAGE:

    CR

    Arguments:

    None

    EXAMPLE:

    CR

    MORE DETAILS:

    Commands to make colored filled-ring cartoon of nucleic acids. May
    need to 'hide everything' first. If asymmetric unit has one strand
    of a dsRNA or dsDNA, remember to apply the BU shortcut to display the
    second strand.
    
    Adapted from KP Wu's blog post:
    https://kpwu.wordpress.com/2012/05/24/pymol-different-colors-of-
    nucleic-acid-rings/


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    hide everything;bg_color white;cartoon oval;set cartoon_ring_mode,3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C, resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A;color red, rna_C;color gray40, rna_G;color palecyan, rna_U;as cartoon 

    PYTHON CODE:

def CR():
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
    DESCRIPTION:

    Commands to color ribbon or cartoon representations of proteins by

    USAGE:

    CSS

    Arguments:

    None

    EXAMPLE:

    CSS

    MORE DETAILS:

    Commands to color ribbon or cartoon representations of proteins by
    secondary structures. Type 'help CSS' to see this documentation
    printed to the command history window. Select from the command 
    history individual lines of code to build a new 
    script. Select the hortizontal script at the bottom if retaining most of 
    the commands in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    VERTICAL PML SCRIPT:

    as cartoon;
    color red, ss H;
    color yellow,ss S;
    color green, ss L+;


    HORIZONTAL PML SCRIPT:

    as cartoon; color red, ss H; color yellow,ss S; color green, ss L+;

    PYTHON CODE:

def CSS():
    cmd.show_as('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')

    '''

    cmd.show_as('cartoon')
    cmd.color('red', 'ss H')
    cmd.color('yellow', 'ss S')
    cmd.color('green', 'ss L+')

cmd.extend('CSS',CSS)

def DU():
    ''' 
    DESCRIPTION:

    Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins. 

    USAGE:

    DU

    Arguments:

    None

    EXAMPLE:

    DU

    MORE DETAILS:

    Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins. 
    Set 50% transparent cartoon so it can be combined with lines, sticks, and ball-and-sticks (try BS shortcut).
    
    Type 'DU' to execute. Type 'help DU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the comand line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

    cartoon dumbbell;
    set cartoon_dumbbell_width, 0.2;
    set cartoon_dumbbell_radius, 0.4;
    show cartoon; 


    HORIZONTAL PML SCRIPT:

    cartoon dumbbell;set cartoon_dumbbell_width, 0.2;set cartoon_dumbbell_radius, 0.4;show cartoon; 


    PYTHON CODE:

def DU():
    cmd.cartoon('dumbbell')
    cmd.set('cartoon_dumbbell_width', '0.2')
    cmd.set('cartoon_dumbbell_radius', '0.4')
    cmd.show('cartoon')

    '''

    cmd.cartoon('dumbbell')
    cmd.set('cartoon_dumbbell_width', '0.2')
    cmd.set('cartoon_dumbbell_radius', '0.4')
    cmd.show('cartoon')

cmd.extend('DU',DU) 

def EMDB():
    ''' 
    DESCRIPTION:

    Open the website of the Electron Microscopy Data Bank.

    USAGE:

    EMDB

    Arguments:

    NA

    EXAMPLE:

    EMDB

    MORE DETAILS:

    Open the website of the Electron Microscopy Data Bank.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def EMDB():
    webbrowser.open('https://www.ebi.ac.uk/pdbe/emdb/')

    '''

    webbrowser.open('https://www.ebi.ac.uk/pdbe/emdb/')

cmd.extend('EMDB',EMDB)

def EP():
    ''' 
    DESCRIPTION:

    EasyPyMOL github site.

    USAGE:

    EP

    Arguments:

    NA

    EXAMPLE:

    EP

    MORE DETAILS:

    EasyPyMOL github site.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def EP():
    webbrowser.open('http://wiki.jmol.org/index.php/Main_Page')
    '''

    webbrowser.open('http://wiki.jmol.org/index.php/Main_Page')
cmd.extend('EP',EP)

def FR():
    ''' 
    DESCRIPTION:

    Make filled-ring cartoon of nucleic acids.

    USAGE:

    FR

    Arguments:

    None

    EXAMPLE:

    FR

    MORE DETAILS:

    Make filled-ring cartoon of nucleic acids. May need to enter 'hide everything' first. 
    Adapted from the script on http://www-cryst.bioc.cam.ac.uk/members/zbyszek/figures_pymol. 
    Type 'FR' to execute. Type 'help FR' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the PyMOL gui.


    VERTICAL PML SCRIPT:

    show sticks;
    set cartoon_ring_mode, 3;
    set cartoon_ring_finder, 1;
    set cartoon_ladder_mode, 1;
    set cartoon_nucleic_acid_mode, 4;
    set cartoon_ring_transparency, 0.5;
    as cartoon;


    HORIZONTAL PML SCRIPT:

    show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon; 

    PYTHON CODE:

def FR():
    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')

    '''

    cmd.show('sticks')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.set('cartoon_ladder_mode', '1')
    cmd.set('cartoon_nucleic_acid_mode', '4')
    cmd.set('cartoon_ring_transparency', '0.5')
    cmd.show_as('cartoon')

cmd.extend('FR',FR)

def GB(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to Google Books in default browser.

    USAGE:

    GB

    Arguments:

    searchTerm

    EXAMPLE:

    GB

    MORE DETAILS:

    Send search term or phrase to Google Books in default browser.

    VERTICAL PML SCRIPT:

    url = 'https://www.google.com/search?tbm=bks&q='
    webbrowser.open(url+searchTerm)


    HORIZONTAL PML SCRIPT:

    url = 'https://www.google.com/search?tbm=bks&q=';webbrowser.open(url+searchTerm)


    PYTHON CODE:

def GB(searchTerm="pymol"):
    url = 'https://www.google.com/search?tbm=bks&q='
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://www.google.com/search?tbm=bks&q='
    webbrowser.open(url+searchTerm)

cmd.extend('GB',GB)

def GGT():
    ''' 
    DESCRIPTION:

    WT human gamma glutamyl transpeptidase at 1.67 Angstrom resolution as cartoon. PDB Code 4gdx.


    USAGE:

    GGT

    Arguments:

    None

    EXAMPLE:

    GGT

    MORE DETAILS:

    WT human gamma glutamyl transpeptidase at 1.67 Angstrom
    resolution as cartoon. PDB Code 4gdx.

    Type 'GGT' to activate. Type 'help GGT' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 4gdx, type=pdb, async=0;remove name H*;as cartoon;bg_color white; hide (name c+o+n);set cartoon_side_chain_helper,  on;color red, 4gdx and ss H; color yellow,4gdx and ss S;color green,4gdx and ss L+; select ASNNAG,resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 ori. 511 or i. 381; color red, elem o and ASNNAG; color blue, elem n and ASNNAG;color yellow, elem c  and ASNNAG;show sticks,ASNNAG;disable ASNNAG; set_view(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); draw 

    PYTHON CODE:

def GGT():
    cmd.reinitialize()
    cmd.fetch('4gdx', type='pdb', async_='0')
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

    '''

    cmd.reinitialize()
    cmd.fetch('4gdx', type='pdb', async_='0')
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

def GH(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to GitHub in default browser.

    USAGE:

    GH

    Arguments:

    searchTerm

    EXAMPLE:

    GH

    MORE DETAILS:

    Send search term or phrase to GitHub in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 


    VERTICAL PML SCRIPT:

    url = 'https://www.github.com/search?q='
    webbrowser.open(url+searchTerm)


    HORIZONTAL PML SCRIPT:

    url = 'https://www.github.com/search?q=';webbrowser.open(url+searchTerm)


    PYTHON CODE:

def GH(searchTerm="pymol"):
    url = 'https://www.github.com/search?q='
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://www.github.com/search?q='
    webbrowser.open(url+searchTerm)

cmd.extend('GH',GH)

def GO(searchTerm="pymol",numHits="200"):
    ''' 
    DESCRIPTION:

    Send search term or phrase Google in default browser.

    USAGE:

    GO

    Arguments:

    searchTerm

    EXAMPLE:

    GO

    MORE DETAILS:

    Send search term or phrase Google in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The second argument is the number of hits to return. 
    The default web browser is used. 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def GO(searchTerm="pymol",numHits="200"):
    webbrowser.open('https://www.google.com/searchq='+searchTerm+'&num='+str(numHits))
    '''

    webbrowser.open('https://www.google.com/searchq='+searchTerm+'&num='+str(numHits))
cmd.extend('GO',GO)

def GS(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to Google Scholar in default browser.

    USAGE:

    GS

    Arguments:

    searchTerm, searchTerm, searchTerm,

    EXAMPLE:

    GS

    MORE DETAILS:

    Send search term or phrase to Google Scholar in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 
    The default search term is pymol.



    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def GS(searchTerm="pymol"):
    url = 'https://scholar.google.se/scholar?hl=en&q='
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://scholar.google.se/scholar?hl=en&q='
    webbrowser.open(url+searchTerm)

cmd.extend('GS',GS)

def GU():
    ''' 
    DESCRIPTION:

    10-mer dsRNA with 8 GU wobble base pairs.


    USAGE:

    GU

    Arguments:

    None

    EXAMPLE:

    GU

    MORE DETAILS:

    10-mer dsRNA with 8 GU wobble base pairs.
    1.32 Angstrom resolution: 4PCO. Has five strands in 
    the asymmetric unit. Deleted chain E and cobalt 
    hexammine 102. Cartoon with filled rings and
    bases cartoon.


Type 'help GU' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 4PCO,type=pdb,async=0;hide everything;bg_color white; cartoon oval;set cartoon_ring_mode, 3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C,resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A; color red, rna_C;color gray40, rna_G; color palecyan, rna_U;as cartoon;disable rna_U; set stick_radius, 0.12;set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;set stick_ball_ratio, 1.8; show sticks, resn NCO;show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0);draw 

    PYTHON CODE:

def GU():
    cmd.reinitialize();
    cmd.fetch('4PCO', type='pdb')
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

    '''

    cmd.reinitialize();
    cmd.fetch('4PCO', type='pdb')
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

def GV(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to Google Videos in default browser.

    USAGE:

    GV

    Arguments:

    GV Linus Pauling; GS Francis Crick; GS Alexander Rich

    EXAMPLE:

    GV

    MORE DETAILS:

    Send search term or phrase to Google Videos in default browser.
    The search phrase does not need to be enclosed in quotes. 
    The default web browser is used. 
    The default search term is pymol.


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def GV(searchTerm="pymol"):
    url = 'https://www.google.com/search?q=video+';
    webbrowser.open(url+searchTerm)

    '''

    url = 'https://www.google.com/search?q=video+';
    webbrowser.open(url+searchTerm)

cmd.extend('GV',GV)

def HH():
    ''' 
    DESCRIPTION:

    Hide hydrogen atoms of currently visible molecular objects.

    USAGE:

    HH

    Arguments:

    None

    EXAMPLE:

    HH

    MORE DETAILS:

    Hide hydrogen atoms of currently visible molecular objects.

    VERTICAL PML SCRIPT:

    hide everything, name H* 

    HORIZONTAL PML SCRIPT:

    hide everything, name H* 

    PYTHON CODE:

def HH():
    cmd.hide('everything', 'name H*')
    '''

    cmd.hide('everything', 'name H*')
cmd.extend('HH',HH)

def IPM(searchTerms = [], *args):
    ''' 
    DESCRIPTION:

    Read list of search terms and submit each term to PubMed in a separate browser tab.

    USAGE:

    IPM

    Arguments:

    search=[string,string]; IPM(search)

    EXAMPLE:

    search=["pymol","vmd","jmol"]; IPM(search)

    MORE DETAILS:

    Read list of search terms and submit each term to PubMed in a separate browser tab.
    There a time delay based on the response time of the site to which the request is made.
    The default web browser is used.
    Must enclose each search term (can be of multiple words) in single or double quotes.
    Has a time delay to avoid overwhelming the webserver.


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def IPM(searchTerms = [], *args):
    print('Sending '+ searchTerms + ' to Pubmed and display list of search results in separate tabs of the default brower.')
    for term in searchTerms:
        t0 = time.time()
        sterm = str(term)
        webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+sterm)
        response_delay = time.time() - t0
        time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
        print('Finished searching PubMed for', sterm, '.')
    print('Finished searching PubMed for ' + searchTerms + '.') 

    '''

    print('Sending '+ searchTerms + ' to Pubmed and display list of search results in separate tabs of the default brower.')
    for term in searchTerms:
        t0 = time.time()
        sterm = str(term)
        webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+sterm)
        response_delay = time.time() - t0
        time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
        print('Finished searching PubMed for', sterm, '.')
    print('Finished searching PubMed for ' + searchTerms + '.') 

cmd.extend('IPM',IPM)

def IUCR(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Open website of the IUCr Journals.

    USAGE:

    IUCR

    Arguments:

    NA

    EXAMPLE:

    IUCR

    MORE DETAILS:

    Open website of the IUCr Journals.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def IUCR(searchTerm="pymol"):
    webbrowser.open('https://journals.iucr.org/')
    '''

    webbrowser.open('https://journals.iucr.org/')
cmd.extend('IUCR',IUCR)

def JASP():
    ''' 
    DESCRIPTION:

    Open JASP from within PyMOL.

    USAGE:

    JASP

    Arguments:

    None

    EXAMPLE:

    JASP

    MORE DETAILS:

    Open JASP from within PyMOL.
    The is a data analysis program that can do Bayesian and frequentist statistics in parallel.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = jaspPath;
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg =jaspPath;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def JASP():
    arg = jaspPath
    subprocess.call(arg,shell=True)
    return

    '''

    arg = jaspPath
    subprocess.call(arg,shell=True)
    return

cmd.extend('JASP',JASP)

def JM():
    ''' 
    DESCRIPTION:

    Open the Jmol wiki.

    USAGE:

    JM

    Arguments:

    NA

    EXAMPLE:

    JM

    MORE DETAILS:

    Open the Jmol wiki.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def JM():
    webbrowser.open('https://github.com/MooersLab/EasyPyMOL')

    '''

    webbrowser.open('https://github.com/MooersLab/EasyPyMOL')

cmd.extend('JM',JM)

def JMP():
    ''' 
    DESCRIPTION:

    Open the JMP from within PyMOL. 

    USAGE:

    JMP

    Arguments:

    None

    EXAMPLE:

    JMP

    MORE DETAILS:

    Open the JMP from within PyMOL. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = jmpPath;
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = jmpPath;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def JMP():
    arg = jmpPath
    subprocess.call(arg,shell=True)
    return

    '''

    arg = jmpPath
    subprocess.call(arg,shell=True)
    return

cmd.extend('JMP',JMP)

def LBSF():
    ''' 
    DESCRIPTION:

    Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.

    USAGE:

    LBSF

    Arguments:

    NA

    EXAMPLE:

    LBSF

    MORE DETAILS:

    Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def LBSF():
    webbrowser.open('https://research.ouhsc.edu/CoreFacilities/LaboratoryofBiomolecularStructureandFunction.aspx')

    '''

    webbrowser.open('https://research.ouhsc.edu/CoreFacilities/LaboratoryofBiomolecularStructureandFunction.aspx')

cmd.extend('LBSF',LBSF)

def LBST():
    ''' 
    DESCRIPTION:

    G2G3/U9U8 base step , PDB code 4PCO. 

    USAGE:

    LBST

    Arguments:

    NA

    EXAMPLE:

    LBST

    MORE DETAILS:

    G2G3/U9U8 base step, PDB code 4PCO. 
    From the 1.32 Angstrom resolution structure 
    of an RNA decamer with 8 GU base pairs.

    Type 'LBST' to execute. Type 'help LBST' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 4PCO.pdb;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png

    PYTHON CODE:

def LBST():
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

def LG():
    ''' 
    DESCRIPTION:

    Nine sugar glycan in influenza N9 neuraminidase, PDB code 4dgr.


    USAGE:

    LG

    Arguments:

    None

    EXAMPLE:

    LG

    MORE DETAILS:

    Nine sugar glycan in influenza N9 neuraminidase at 1.55 Angstrom  resolution, PDB code 4dgr. 
    The electron density map is contoured at 1.0 sigma. 
    39 commands were used to make this figure.  
    Type 'LG' to execute. Type 'help LG' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    delete all;fetch 4dgr, async=0;fetch 4dgr, type=2fofc, async=0;select LongGlycan, resi 469:477;orient LongGlycan;remove not LongGlycan;remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;color density, 2fofcmap; show sticks;show spheres;set stick_radius, .07;set sphere_scale, .19;set sphere_scale, .13, elem H;set bg_rgb=[1, 1, 1];set stick_quality, 50;set sphere_quality, 4;color gray85, elem C;color red, elem O;color slate, elem N;color gray98, elem H;set stick_color, gray50;set ray_trace_mode, 1;set ray_texture, 2;set antialias, 3;set ambient, 0.5;set spec_count, 5;set shininess, 50;set specular, 1;set reflect, .1;set dash_gap, 0;set dash_color, black;set dash_gap, .15;set dash_length, .05;set dash_round_ends, 0;set dash_radius, .05;set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);preset.ball_and_stick("all",mode=1);draw 


    PYTHON CODE:

def LG():
    cmd.reinitialize()
    cmd.fetch('4dgr')
    cmd.fetch('4dgr', type='2fofc')
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

    '''

    cmd.reinitialize()
    cmd.fetch('4dgr')
    cmd.fetch('4dgr', type='2fofc')
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

def LGGT():
    ''' 
    DESCRIPTION:

    WT human gamma glutamyl transpeptidase as cartoon. PDB code 4gdx. 



    USAGE:

    LGGT

    Arguments:

    NA

    EXAMPLE:

    LGGT

    MORE DETAILS:

    WT human gamma glutamyl transpeptidase at 1.67 Angstrom
    resolution as cartoon. PDB code 4gdx. 
    4gdx.pdb must be in the current working directory. 

    Type 'LGGT' to activate. Type 'help LGGT' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    delete all;load 4gdx.pdb;remove name H*;as cartoon;bg_color white; hide (name c+o+n);set cartoon_side_chain_helper,  on;color red, 4gdx and ss H; color yellow,4gdx and ss S;color green,4gdx and ss L+; select ASNNAG,resn NAG or resi 95 or i. 120  or i. 230 or i. 266 or i. 344 ori. 511 or i. 381; color red, elem o and ASNNAG; color blue, elem n and ASNNAG;color yellow, elem c  and ASNNAG;show sticks,ASNNAG;disable ASNNAG; set_view(0.55,-0.83,0.07,0.5,0.26,-0.82,0.66,0.49,0.56,0.0,0.0,-197.16,-22.42,-22.69,-12.01,155.44,238.88,-20.0); draw 


    PYTHON CODE:

def LGGT():
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
    DESCRIPTION:

    10-mer dsRNA. 


    USAGE:

    LGU

    Arguments:

    NA

    EXAMPLE:

    LGU

    MORE DETAILS:

    10-mer dsRNA. 
    1.32 Angstrom resolution: 4PCO. Has five strands in 
    the asymmetric unit. Deleted chain E and cobalt 
    hexammine 102. Cartoon with filled rings and
    bases cartoon.


    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    delete all;load 4PCO.pdb,async=0;hide everything;bg_color white; cartoon oval;set cartoon_ring_mode, 3;set cartoon_nucleic_acid_color, blue;select rna_A, resn A;select rna_C,resn C;select rna_G, resn G;select rna_U, resn U;color yellow, rna_A; color red, rna_C;color gray40, rna_G; color palecyan, rna_U;as cartoon;disable rna_U; set stick_radius, 0.12;set nb_spheres_size, 0.3; show nb_spheres; set stick_ball, on;set stick_ball_ratio, 1.8; show sticks, resn NCO;show spheres, name Cl; set_view (0.34,-0.81,0.48,0.89,0.11,-0.45,0.31,0.58,0.76,-0.0,0.0,-196.36,-9.82,6.76,15.84,159.01,233.71,-20.0);draw 


    PYTHON CODE:

def LGU():
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

def LLG():
    ''' 
    DESCRIPTION:

    Nine sugar glycan in influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr.


    USAGE:

    LLG

    Arguments:

    NA

    EXAMPLE:

    LLG

    MORE DETAILS:

    Nine sugar glycan in influenza N9 neuraminidase at 
    1.55 Angstrom  resolution, PDB code 4dgr. 
    The electron density map is contoured at 1.0 sigma. 
    39 commands were used to make this figure.  

    Type 'LLG' to execute. Type 'help LLG' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 4dgr.pdb;fetch 4dgr2FoFc.mtz;select LongGlycan, resi 469:477;orient LongGlycan;remove not LongGlycan;remove name H*;isomesh 2fofcmap, 4dgr_2fofc, 1, LongGlycan, carve = 1.8;color density, 2fofcmap; show sticks;show spheres;set stick_radius, .07;set sphere_scale, .19;set sphere_scale, .13, elem H;set bg_rgb=[1, 1, 1];set stick_quality, 50;set sphere_quality, 4;color gray85, elem C;color red, elem O;color slate, elem N;color gray98, elem H;set stick_color, gray50;set ray_trace_mode, 1;set ray_texture, 2;set antialias, 3;set ambient, 0.5;set spec_count, 5;set shininess, 50;set specular, 1;set reflect, .1;set dash_gap, 0;set dash_color, black;set dash_gap, .15;set dash_length, .05;set dash_round_ends, 0;set dash_radius, .05;set_view (0.34,-0.72,0.61,0.8,0.56,0.22,-0.51,0.4,0.77,0.0,0.0,-81.31,44.64,-9.02,58.62,65.34,97.28,-20.0);preset.ball_and_stick("all",mode=1);draw 

    PYTHON CODE:

def LLG():
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

def LN9():
    ''' 
    DESCRIPTION:

    Influenza N9 neuraminidase, PDB code 4dgr.

    USAGE:

    LN9

    Arguments:

    NA

    EXAMPLE:

    LN9

    MORE DETAILS:

    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code
    4dgr. The biological unit has four copies of the asymmetric unit.
    View is down the four-fold axis. Requires the quat.py script by
    Thomas Holder and available at the PyMOL Wiki page. Store quat.py
    in ~/mg18OU.

    Type 'LN9' to activate. Type 'help LN9' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 4dgr.pdb;run $HOME/mg18OU/quat.py; quat 4dgr;as cartoon; bg_color white;color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;color green, 4dgr_1 and ss L+;color cyan, (not 4dgr_1 and ss H);color magenta, (not 4dgr_1 and ss S);color orange, (not 4dgr_1 and ss L+);set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0); draw 


    PYTHON CODE:

def LN9():
    cmd.reinitialize()
    cmd.load('4dgr.pdb')
    #cmd.do('run $HOME/mg18OU/quat.py')
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

    '''

    cmd.reinitialize()
    cmd.load('4dgr.pdb')
    #cmd.do('run $HOME/mg18OU/quat.py')
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

def LNA():
    ''' 
    DESCRIPTION:

    Hydrated sodium cation bound in major groove of RNA with 16 Watson-Crick base pairs.


    USAGE:

    LNA

    Arguments:

    NA

    EXAMPLE:

    LNA

    MORE DETAILS:

    Hydrated sodium cation bound in major groove of RNA with 16 Watson-Crick base pairs.
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

    Type 'NA' to execute. Type 'help NA' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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

    HORIZONTAL PML SCRIPT:

    delete all;viewport 900,600;load 3nd4.pdb;hide cartoon;run ~/mg18OU/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;rock;AOset_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 

    PYTHON CODE:

def LNA():
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

def LT4L():
    ''' 
    DESCRIPTION:

    Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0. 


    USAGE:

    LT4L

    Arguments:

    NA

    EXAMPLE:

    LT4L

    MORE DETAILS:

    Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0. 
    The file 3FA0 must be in the current working directory. 


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 3fa0.pdb;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white;show cartoon;color red, ss H;color yellow, ss S;color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); ray 1500,1600; 


    PYTHON CODE:

def LT4L():
    cmd.reinitialize()
    cmd.load('3fa0.pdb')
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

    '''

    cmd.reinitialize()
    cmd.load('3fa0.pdb')
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
    DESCRIPTION:

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.


    USAGE:

    LU8

    Arguments:

    NA

    EXAMPLE:

    LU8

    MORE DETAILS:

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.
    Has one strand in the asymmetric unit. Uses quat.py to generate
    the second strand. Cartoon with filled rings and bases cartoon.
    The file 3nd3.pdb needs to be in the current working directory.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 3nd3.pdb;run $HOME/mg18OU/quat.py;quat 3nd3;hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon;set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 

    PYTHON CODE:

def LU8():
    cmd.reinitialize()
    cmd.load('3nd3.pdb')
    #cmd.do('run $HOME/mg18OU/quat.py')
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

    '''

    cmd.reinitialize()
    cmd.load('3nd3.pdb')
    #cmd.do('run $HOME/mg18OU/quat.py')
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
    DESCRIPTION:

    16-mer dsRNA, Watson-Crick helix RNA, 3nd4.

    USAGE:

    LWC8

    Arguments:

    NA

    EXAMPLE:

    LWC8

    MORE DETAILS:

    16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom 
    resolution: 3nd4.  Has one strand in the asymmetric unit. 
    Needs quat.py to generate the second strand. 
    Cartoon with filled rings and bases cartoon.
    The file 3nd4.pdb must be in the current working directory.

    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;load 3nd4.pdb;hide everything;run $HOME/mg18OU/quat.py; quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; set nb_spheres_size, 0.25; show nb_spheres; set stick_ball, on; set stick_ball_ratio, 1.8;set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);hide everything, name H*;rock 

    PYTHON CODE:

def LWC8():
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

def MA(searchTerm='pymol'):
    ''' 
    DESCRIPTION:

    Send search term to all searchable websites in pymolshortcuts.

    USAGE:

    MA

    Arguments:

    search term(s), N

    EXAMPLE:

    MA pymol plugin

    MORE DETAILS:

    Send search term to all searchable websites in pymolshortcuts.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MA(searchTerm='pymol'):
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
    DESCRIPTION:

    Send search term to multiple sites that contain book content.

    USAGE:

    MB

    Arguments:

    search term(s)

    EXAMPLE:

    MB pymol plugin

    MORE DETAILS:

    Send search term to search multiple sites that contain book content.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MB(searchTerm='pymol'):
    GB(searchTerm)
    SD(searchTerm)
    SP(searchTerm)

    '''

    GB(searchTerm)
    SD(searchTerm)
    SP(searchTerm)

cmd.extend('MB',MB)

def MC(searchTerm='pymol'):
    ''' 
    DESCRIPTION:

    Send search term to search ten core websites in pymolshortcuts:

    USAGE:

    MC

    Arguments:

    searchTerm

    EXAMPLE:

    MC pymol

    MORE DETAILS:

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


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MC(searchTerm='pymol'):
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

def MCL():
    ''' 
    DESCRIPTION:

    Open website of Macromolecular Crystallography Laboratory at the University of Oklahoma.

    USAGE:

    MCL

    Arguments:

    NA

    EXAMPLE:

    MCL

    MORE DETAILS:

    Open website of Macromolecular Crystallography Laboratory at the University of Oklahoma.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MCL():
    webbrowser.open('http://structuralbiology.ou.edu/mcl')
    '''

    webbrowser.open('http://structuralbiology.ou.edu/mcl')
cmd.extend('MCL',MCL)

def MG():
    ''' 
    DESCRIPTION:

    Open website of the OUHSC molecular graphics course.

    USAGE:

    MG

    Arguments:

    NA

    EXAMPLE:

    MG

    MORE DETAILS:

    Open website of the OUHSC molecular graphics course.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MG():
    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/moleculargraphicslinks.html')
    '''

    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/moleculargraphicslinks.html')
cmd.extend('MG',MG)

def MM(searchTerm='pymol'):
    ''' 
    DESCRIPTION:

    Send search term to search for manuscripts in pymolshortcuts.

    USAGE:

    MM

    Arguments:

    searchTerm

    EXAMPLE:

    MM pymol

    MORE DETAILS:

    Send search term to search for manuscripts in pymolshortcuts:

    arXiv
    bioRxiv
    Google Scholar
    PubMed
    Research Gate
    Science Direct
    Springer


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def MM(searchTerm='pymol'):
    AX(searchTerm)
    BX(searchTerm)
    GS(searchTerm)
    PM(searchTerm)
    RG(searchTerm)
    SD(searchTerm)
    SP(searchTerm)

    '''

    AX(searchTerm)
    BX(searchTerm)
    GS(searchTerm)
    PM(searchTerm)
    RG(searchTerm)
    SD(searchTerm)
    SP(searchTerm)

cmd.extend('MM',MM)

def N9():
    ''' 
    DESCRIPTION:

    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr.


    USAGE:

    N9

    Arguments:

    None

    EXAMPLE:

    N9

    MORE DETAILS:

    Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr.
    The biological unit has four copies of the asymmetric unit.
    View is down the four-fold axis. Requires the quat.py script by
    Thomas Holder that is available at the PyMOL Wiki page. 

    Type 'N9' to activate. Type 'help N9' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 4dgr, type=pdb, async=0;run $HOME/mg18OU/quat.py; quat 4dgr;as cartoon; bg_color white;color red, 4dgr_1 and ss H;color yellow,4dgr_1 and ss S;color green, 4dgr_1 and ss L+;color cyan, (not 4dgr_1 and ss H);color magenta, (not 4dgr_1 and ss S);color orange, (not 4dgr_1 and ss L+);set_view (0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0); draw 


    PYTHON CODE:

def N9():
    cmd.reinitialize()
    cmd.fetch('4dgr', type='pdb1')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.do('split_state 4dgr')	
    cmd.color('red', '4dgr_0002 and ss H')
    cmd.color('yellow', '4dgr_0002 and ss S')
    cmd.color('green', '4dgr_0002 and ss L+')
    cmd.color('cyan', '(not 4dgr_0002 and ss H)')
    cmd.color('magenta', '(not 4dgr_0002 and ss S)')
    cmd.color('orange', '(not 4dgr_0002 and ss L+)')
    cmd.set_view('(0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0)')
    cmd.draw()

    '''

    cmd.reinitialize()
    cmd.fetch('4dgr', type='pdb1')
    cmd.show_as('cartoon')
    cmd.bg_color('white')
    cmd.do('split_state 4dgr')	
    cmd.color('red', '4dgr_0002 and ss H')
    cmd.color('yellow', '4dgr_0002 and ss S')
    cmd.color('green', '4dgr_0002 and ss L+')
    cmd.color('cyan', '(not 4dgr_0002 and ss H)')
    cmd.color('magenta', '(not 4dgr_0002 and ss S)')
    cmd.color('orange', '(not 4dgr_0002 and ss L+)')
    cmd.set_view('(0.98,-0.22,0.01,0.22,0.98,0.02,-0.01,-0.02,1.0,-0.0,0.0,-323.44,1.46,5.33,56.19,274.72,372.15,-20.0)')
    cmd.draw()

cmd.extend('N9',N9)

def NA():
    ''' 
    DESCRIPTION:

    Hydrated sodium cation bound in major groove of a 16-mer RNA of Watson-Crick base pairs.


    USAGE:

    NA

    Arguments:

    None

    EXAMPLE:

    NA

    MORE DETAILS:

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


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;viewport 900,600;fetch 3nd4, type=pdb,async=0;run ~/mg18OU/quat.py;quat 3nd4; show sticks;set stick_radius=0.125;hide everything, name H*;bg_color white;create coorCov, (3nd4_1 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)));bond (coorCov//A/NA`19/NA),(coorCov//A/A`3/N7); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`119/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`219/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`319/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`419/O); bond (coorCov//A/NA`19/NA),(coorCov//A/HOH`519/O);distance (3nd4_1 and chain Aand resi 19 and name NA), (3nd4_1 and chain A and resi 519);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 419);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 119);distance (3nd4_1 and chain A and resi 19 and name NA),(3nd4_1 and chain A and resi 319);distance (3nd4_1 and chain A and resi 19 and name NA), (3nd4_1 and chain A and resi 219);show nb_spheres; set nb_spheres_size, .35;distance hbond1,/3nd4_1/1/A/HOH`119/O, /3nd4_1/1/A/A`3/OP2;distance hbond2,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/A`3/OP2;distance hbond3,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/HOH`119/O;distance hbond4,/3nd4_1/1/A/G`4/N7,/3nd4_1/1/A/HOH`91/O;distance hbond5,/3nd4_1/1/A/G`4/O6, /3nd4_1/1/A/HOH`419/O;distance hbond6,/3nd4_1/1/A/HOH`91/O, /3nd4_1/1/A/G`4/OP2;distance hbond7,/3nd4_1/1/A/HOH`319/O, /3nd4_1/1/A/G`2/OP2;distance  hbond9,/3nd4_1/1/A/HOH`419/O,/3nd4_2/2/A/HOH`74/O;distance hbond10,/3nd4_2/2/A/C`15/O2,/3nd4_1/1/A/G`2/N2;distance hbond11, /3nd4_2/2/A/C`15/N3,/3nd4_1/1/A/G`2/N1;distance hbond12,/3nd4_2/2/A/C`15/N4,/3nd4_1/1/A/G`2/O6;distance hbond13, /3nd4_2/2/A/U`14/N3,/3nd4_1/1/A/A`3/N1;distance hbond14,3nd4_2/2/A/U`14/O4,/3nd4_1/1/A/A`3/N6;distance hbond15, /3nd4_2/2/A/C`13/N4,/3nd4_1/1/A/G`4/O6;distance hbond16,/3nd4_2/2/A/C`13/N3, /3nd4_1/1/A/G`4/N1;distance hbond17, /3nd4_1/1/A/G`4/N2,/3nd4_2/2/A/C`13/O2;distance hbond18,/3nd4_1/1/A/G`2/N2,/3nd4_2/2/A/C`15/O2;distance hbond19,/3nd4_1/1/A/HOH`91/O,/3nd4_1/1/A/G`4/OP2;set depth_cue=0;set ray_trace_fog=0;set dash_color, black;set label_font_id, 5;set label_size, 36;set label_position, (0.5, 1.0, 2.0);set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;set label_color, black;set dash_gap, 0.2;set dash_width, 2.0;set dash_length, 0.2;select carbon, element C; color yellow, carbon;disable carbon;set_view (-0.9,0.34,-0.26,0.33,0.18,-0.93,-0.27,-0.92,-0.28,-0.07,-0.23,-27.83,8.63,19.85,13.2,16.0,31.63,-20.0); 


    PYTHON CODE:

def NA():
    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.fetch('3nd4', type='pdb1');
    cmd.do('split_states 3nd4')
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'name H*');
    cmd.bg_color('white');
    cmd.create('coorCov', '(3nd4_0001 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)))');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/A`3/N7)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`119/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`219/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`319/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`419/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`519/O)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 519)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 419)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 119)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 319)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 219)');
    cmd.show('nb_spheres');
    cmd.set('nb_spheres_size', '.35');
    cmd.distance('hbond1', '/3nd4_0001/1/A/HOH`119/O', '/3nd4_0001/1/A/A`3/OP2');
    cmd.distance('hbond2', '/3nd4_0001/1/A/HOH`319/O', '/3nd4_0001/1/A/A`3/OP2');
    cmd.distance('hbond3', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/HOH`119/O');
    cmd.distance('hbond4', '/3nd4_0001/1/A/G`4/N7', '/3nd4_0001/1/A/HOH`91/O');
    cmd.distance('hbond5', '/3nd4_0001/1/A/G`4/O6', '/3nd4_0001/1/A/HOH`419/O');
    cmd.distance('hbond6', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/G`4/OP2');
    cmd.distance('hbond7', '/3nd4_0001/1/A/HOH`319/O', '/3nd4_0001/1/A/G`2/OP2');
    cmd.distance('hbond9', '/3nd4_0001/1/A/HOH`419/O', '/3nd4_0002/2/A/HOH`74/O');
    cmd.distance('hbond10', '/3nd4_0002/2/A/C`15/O2', '/3nd4_0001/1/A/G`2/N2');
    cmd.distance('hbond11', '/3nd4_0002/2/A/C`15/N3', '/3nd4_0001/1/A/G`2/N1');
    cmd.distance('hbond12', '/3nd4_0002/2/A/C`15/N4', '/3nd4_0001/1/A/G`2/O6');
    cmd.distance('hbond13', '/3nd4_0002/2/A/U`14/N3', '/3nd4_0001/1/A/A`3/N1');
    cmd.distance('hbond14', '/3nd4_0002/2/A/U`14/O4', '/3nd4_0001/1/A/A`3/N6');
    cmd.distance('hbond15', '/3nd4_0002/2/A/C`13/N4', '/3nd4_0001/1/A/G`4/O6');
    cmd.distance('hbond16', '/3nd4_0002/2/A/C`13/N3', '/3nd4_0001/1/A/G`4/N1');
    cmd.distance('hbond17', '/3nd4_0001/1/A/G`4/N2', '/3nd4_0002/2/A/C`13/O2');
    cmd.distance('hbond18', '/3nd4_0001/1/A/G`2/N2', '/3nd4_0002/2/A/C`15/O2');
    cmd.distance('hbond19', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/G`4/OP2');
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

    '''

    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.fetch('3nd4', type='pdb1');
    cmd.do('split_states 3nd4')
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'name H*');
    cmd.bg_color('white');
    cmd.create('coorCov', '(3nd4_0001 and (resi 19 or resi 119 or resi 219 or resi 319 or resi 419 or resi 519 or (resi 3 and name N7)))');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/A`3/N7)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`119/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`219/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`319/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`419/O)');
    cmd.bond('(coorCov//A/NA`19/NA)','(coorCov//A/HOH`519/O)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 519)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 419)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 119)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 319)');
    cmd.distance('(3nd4_0001 and chain A and resi 19 and name NA)','(3nd4_0001 and chain A and resi 219)');
    cmd.show('nb_spheres');
    cmd.set('nb_spheres_size', '.35');
    cmd.distance('hbond1', '/3nd4_0001/1/A/HOH`119/O', '/3nd4_0001/1/A/A`3/OP2');
    cmd.distance('hbond2', '/3nd4_0001/1/A/HOH`319/O', '/3nd4_0001/1/A/A`3/OP2');
    cmd.distance('hbond3', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/HOH`119/O');
    cmd.distance('hbond4', '/3nd4_0001/1/A/G`4/N7', '/3nd4_0001/1/A/HOH`91/O');
    cmd.distance('hbond5', '/3nd4_0001/1/A/G`4/O6', '/3nd4_0001/1/A/HOH`419/O');
    cmd.distance('hbond6', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/G`4/OP2');
    cmd.distance('hbond7', '/3nd4_0001/1/A/HOH`319/O', '/3nd4_0001/1/A/G`2/OP2');
    cmd.distance('hbond9', '/3nd4_0001/1/A/HOH`419/O', '/3nd4_0002/2/A/HOH`74/O');
    cmd.distance('hbond10', '/3nd4_0002/2/A/C`15/O2', '/3nd4_0001/1/A/G`2/N2');
    cmd.distance('hbond11', '/3nd4_0002/2/A/C`15/N3', '/3nd4_0001/1/A/G`2/N1');
    cmd.distance('hbond12', '/3nd4_0002/2/A/C`15/N4', '/3nd4_0001/1/A/G`2/O6');
    cmd.distance('hbond13', '/3nd4_0002/2/A/U`14/N3', '/3nd4_0001/1/A/A`3/N1');
    cmd.distance('hbond14', '/3nd4_0002/2/A/U`14/O4', '/3nd4_0001/1/A/A`3/N6');
    cmd.distance('hbond15', '/3nd4_0002/2/A/C`13/N4', '/3nd4_0001/1/A/G`4/O6');
    cmd.distance('hbond16', '/3nd4_0002/2/A/C`13/N3', '/3nd4_0001/1/A/G`4/N1');
    cmd.distance('hbond17', '/3nd4_0001/1/A/G`4/N2', '/3nd4_0002/2/A/C`13/O2');
    cmd.distance('hbond18', '/3nd4_0001/1/A/G`2/N2', '/3nd4_0002/2/A/C`15/O2');
    cmd.distance('hbond19', '/3nd4_0001/1/A/HOH`91/O', '/3nd4_0001/1/A/G`4/OP2');
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

def NDB():
    ''' 
    DESCRIPTION:

    Open website of the Nucleic Acid Database.

    USAGE:

    NDB

    Arguments:

    NA

    EXAMPLE:

    NDB

    MORE DETAILS:

    Open website of the Nucleic Acid Database.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def NDB():
    webbrowser.open('https://en.wikipedia.org/wiki/List_of_molecular_graphics_systems')
    '''

    webbrowser.open('https://en.wikipedia.org/wiki/List_of_molecular_graphics_systems')
cmd.extend('NDB',NDB)

def NSLSII():
    ''' 
    DESCRIPTION:

    Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory.

    USAGE:

    NSLSII

    Arguments:

    NA

    EXAMPLE:

    NSLSII

    MORE DETAILS:

    Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def NSLSII():
    webbrowser.open('https://www.bnl.gov/ps/')
    '''

    webbrowser.open('https://www.bnl.gov/ps/')
cmd.extend('NSLSII',NSLSII)

def PDB(searchTerm="3fa0",numHits="5"):
    ''' 
    DESCRIPTION:

    Submit a search term to the Protein Data Bank.

    USAGE:

    PDB

    Arguments:

    searchTerm

    EXAMPLE:

    PBB 3fa0

    MORE DETAILS:

    Submit a search term to the Protein Data Bank.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PDB(searchTerm="3fa0",numHits="5"):
    webbrowser.open('https://www.rcsb.org/structure/'+searchTerm)
    '''

    webbrowser.open('https://www.rcsb.org/structure/'+searchTerm)
cmd.extend('PDB',PDB)

def PE(selection):
    ''' 
    DESCRIPTION:

    Apply pearl effect about cations.

    USAGE:

    PE

    Arguments:

    selection

    EXAMPLE:

    PE selection

    MORE DETAILS:

    Apply pearl effect about cations. 
    This effect shows a transparent sphere around an inner opaque sphere.
    Must supply the <selection> of a cation or anion. 

    VERTICAL PML SCRIPT:

    select magnesium1, <selection>
    create magnesium2, magnesium1
    show spheres, magnesium1
    show spheres, magnesium2
    set spehrical transparency, 0.4, magnesium2
    set spehrescale, 1.05, magnesium2

    HORIZONTAL PML SCRIPT:

    select magnesium1, <selection>;create magnesium2, magnesium1;show spheres, magnesium1;show spheres, magnesium2;set spehrical transparency, 0.4, magnesium2;set spehrescale, 1.05, magnesium2

    PYTHON CODE:

def PE(selection):
    cmd.select('magnesium1',selection)
    cmd.create('magnesium2', 'magnesium1')
    cmd.show(spheres, 'magnesium1')
    cmd.show(spheres, 'magnesium2')
    cmd.set('spehrical transparency, 0.4, magnesium2')
    cmd.set('spehrescale, 1.05, magnesium2')

    '''

    cmd.select('magnesium1',selection)
    cmd.create('magnesium2', 'magnesium1')
    cmd.show(spheres, 'magnesium1')
    cmd.show(spheres, 'magnesium2')
    cmd.set('spehrical transparency, 0.4, magnesium2')
    cmd.set('spehrescale, 1.05, magnesium2')

cmd.extend('PE',PE)

def PM(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Send search term or phrase to PubMed.

    USAGE:

    PM

    Arguments:

    searchTerm

    EXAMPLE:

    PM

    MORE DETAILS:

    Send search term or phrase to PubMed.
    The default web browser is used.
    The multi word search terms do not need to be enclosed in quotes. 
    Takes one search term but multiple commands can be submitted at once (see below).


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PM(searchTerm="pymol"):
    webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+searchTerm)

    '''

    webbrowser.open('https://www.ncbi.nlm.nih.gov/pubmed/?term='+searchTerm)

cmd.extend('PM',PM)

def PML(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:

    Submit a search term to the PyMOL Users Mail Service.

    USAGE:

    PML

    Arguments:

    searchTerm

    EXAMPLE:

    PML session file

    MORE DETAILS:

    Submit a search term to the PyMOL Users Mail Service.

    Single term search (multi word searches do NOT have to be inside quotes):
    PML session file

    Multiple term search: 
    PML text editor; PML 3d pdf; PML black and white cartoon;


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PML(searchTerm="3d_pdf"):
    webbrowser.open('https://sourceforge.net/p/pymol/mailman/search/?q='+searchTerm)

    '''

    webbrowser.open('https://sourceforge.net/p/pymol/mailman/search/?q='+searchTerm)

cmd.extend('PML',PML)

def PPC():
    ''' 
    DESCRIPTION:

    Open the website of the Protein Production Facility at the University of Oklahoma in Norman.

    USAGE:

    PPC

    Arguments:

    NA

    EXAMPLE:

    PPC

    MORE DETAILS:

    Open the website of the Protein Production Facility at the University of Oklahoma in Norman.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PPC():
    webbrowser.open('http://www.ou.edu/cas/chemistry/research/research-support-services/protein-production-core')

    '''

    webbrowser.open('http://www.ou.edu/cas/chemistry/research/research-support-services/protein-production-core')

cmd.extend('PPC',PPC)

def PS():
    ''' 
    DESCRIPTION:

    Open the home page of the Protein Soceity.

    USAGE:

    PS

    Arguments:

    NA

    EXAMPLE:

    PS

    MORE DETAILS:

    Open the home page of the Protein Soceity.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PS():
    webbrowser.open('https://www.proteinsociety.org/')
    '''

    webbrowser.open('https://www.proteinsociety.org/')
cmd.extend('PS',PS)

def PW(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:

    Submit search of the PyMOL Wiki. 

    USAGE:

    PW

    Arguments:

    NA

    EXAMPLE:

    PW

    MORE DETAILS:

    Submit search of the PyMOL Wiki.    

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def PW(searchTerm="3d_pdf"):
    webbrowser.open('https://pymolwiki.org/index.php/'+searchTerm)
    '''

    webbrowser.open('https://pymolwiki.org/index.php/'+searchTerm)
cmd.extend('PW',PW)

def RG(searchTerm='best molecular graphics program'):
    ''' 
    DESCRIPTION:

    Submit a search query of Research Gate. 

    USAGE:

    RG

    Arguments:

    searchTerm

    EXAMPLE:

    RG best molecular graphics program

    MORE DETAILS:

    Submit a search query of Research Gate. 

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def RG(searchTerm='best molecular graphics program'):
    webbrowser.open('https://www.researchgate.net/search.Search.html?type=researcher&query='+searchTerm)
    '''

    webbrowser.open('https://www.researchgate.net/search.Search.html?type=researcher&query='+searchTerm)
cmd.extend('RG',RG)

def RS():
    ''' 
    DESCRIPTION:

    Open the homepage of the RNA Society.

    USAGE:

    RS

    Arguments:

    NA

    EXAMPLE:

    RS

    MORE DETAILS:

     Open the homepage of the RNA Society.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def RS():
    webbrowser.open('https://www.rnasociety.org/')

    '''

    webbrowser.open('https://www.rnasociety.org/')

cmd.extend('RS',RS)

def SAXS():
    ''' 
    DESCRIPTION:

    Open the webpage of SAXS links at OUHSC. 

    USAGE:

    SAXS

    Arguments:

    NA

    EXAMPLE:

    SAXS

    MORE DETAILS:

    Open the webpage of SAXS links at OUHSC. 

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SAXS():
    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/small-angle-scattering-links-27aug2014.html?sfvrsn=0')
    '''

    webbrowser.open('https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/small-angle-scattering-links-27aug2014.html?sfvrsn=0')
cmd.extend('SAXS',SAXS)

def SB():
    ''' 
    DESCRIPTION:

    Open the webpage of SSRL Biological SAXS at BL 4-2.

    USAGE:

    SB

    Arguments:

    NA

    EXAMPLE:

    SB

    MORE DETAILS:

    Open the webpage of SSRL Biological SAXS at BL 4-2.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SB():
    webbrowser.open('https://www-ssrl.slac.stanford.edu/~saxs/')
    '''

    webbrowser.open('https://www-ssrl.slac.stanford.edu/~saxs/')
cmd.extend('SB',SB)

def SBGRID():
    ''' 
    DESCRIPTION:

    Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel.

    USAGE:

    SBGRID

    Arguments:

    NA

    EXAMPLE:

    SBGRID

    MORE DETAILS:

    Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SBGRID():
    webbrowser.open('https://www.youtube.com/user/SBGridTV/videos')
    '''

    webbrowser.open('https://www.youtube.com/user/SBGridTV/videos')
cmd.extend('SBGRID',SBGRID)

def SC():
    ''' 
    DESCRIPTION:

    Print to screen list of the shortcuts that are available in the script pymolshortcuts.py. 

    USAGE:

    SC

    Arguments:

    NA

    EXAMPLE:

    SC

    MORE DETAILS:

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

    Print shortcuts and their descriptions:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    SC        Print to screen list of the shortcuts that are available in the script pymolshortcuts.py.   

    Show many models (NMR and crystal packing):
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    nmr          Show all of the models in nmr structure.                         
    nmroff      Hide all but first model in a nmr structure. 
    rmd          Remove all measurement objects.                    
    rmsc         Remove supercell and the symmetry mates.                         
    sc111       Make a lattice of 1 x 1 x 1 unit cells.
    sc211       Make a lattice of 2 x 1 x 1 unit cells 
    sc121       Make a lattice of 1 x 2 x 1 unit cells                        
    sc112       Make a lattice of 1 x 1 x 2 unit cells.                                    
    sc221       Make a lattice of 2 x 2 x 1 unit cells.        
    sc122       Make a lattice of 2 x 2 x 1 unit cells.        
    sc212       Make a lattice of 2 x 2 x 1 unit cells.        
    sc222       Make a lattice of 2 x 2 x 2 unit cells.   
    sc113       Make a lattice of 1 x 1 x 3 unit cells.                            
    sc311       Make a lattice of 3 x 1 x 1 unit cells.                            
    sc131       Make a lattice of 3 x 1 x 1 unit cells.                            
    sc133       Make a lattice of 1 x 3 x 3 unit cells.                            
    sc313       Make a lattice of 3 x 1 x 3 unit cells.                            
    sc331       Make a lattice of 3 x 3 x 1 unit cells. 
    sc123       Make a lattice of 1 x 2 x 3 unit cells. 
    sc132       Make a lattice of 1 x 3 x 2 unit cells.               
    sc321       Make a lattice of 3 x 2 x 1 unit cells. 
    sc231       Make a lattice of 2 x 3 x 1 unit cells. 
    sc213       Make a lattice of 2 x 1 x 3 unit cells.
    sc312       Make a lattice of 3 x 1 x 2 unit cells.
    sc321       Make a lattice of 3 x 1 x 2 unit cells.                         
    sc133       Make a lattice of 1 x 3 x 3 unit cells.                                       
    sc333       Make a lattice of 3 x 3 x 3 unit cells.                          

    Save files with date and time in filename:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    saln          Save a aln file (alignment file) with a time stamp included in the filename to avoid overwriting work.
    scif           Save a cif file (Crystallographic Information File) with a time stamp included in the filename to avoid overwriting work. 
    sccp4       Save a ccp4 file (CCP4 electron density map file) with a time stamp included in the filename to avoid overwriting work.
    sdae         Save a dae file (Collada File) with a time stamp included in the filename to avoid overwriting work.
    sdat          Save dat file (output data file) with a time stamp included in the filename to avoid overwriting work. 
    sfasta       Save a fasta file (sequence file) with a time stamp included in the filename to avoid overwriting work.
    sidtf         Save a idtf file (Intermediate Data Text Format) with a time stamp included in the filename to avoid overwriting work.. 
    smae        Save mae file (Maestro file) with a time stamp included in the filename to avoid overwriting work. 
    smmd       Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work. 
    smmod     Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work. 
    spmo        Save pmo file (XYZ, binary format file) with a time stamp included in the filename to avoid overwriting work. 
    smoe        Save moe file (Molecular Operating Environment) with a time stamp included in the filename to avoid overwriting work. 
    smol         Save mol file with a time stamp included in the filename to avoid overwriting work. 
    smol2       Save mol2 (Sybyl file format) file with a time stamp included in the filename to avoid overwriting work. 
    smtl          Save mtl (Wavefront Material file format) file with a time stamp included in the filename to avoid overwriting work. 
    sobj          Save obj file (Wavefront mesh file) with a time stamp included in the filename to avoid overwriting work. 
    sout          Save out file (output data file) with a time stamp included in the filename to avoid overwriting work. 
    spdb         Save pdb file with a time stamp included in the filename to avoid overwriting work.
    spkl          Save a pkl file (Python pickle file) with a time stamp included in the filename to avoid overwriting work.
    spkla        Save a pkla file (Python pickle file) with a time stamp included in the filename to avoid overwriting work.
    spng         Save png file with a time stamp included in the filename to avoid overwriting work. 
    spov         Save pov (POV-ray tracing file format) file with a time stamp included in the filename to avoid overwriting work. 
    spqr         Save pqr file with a time stamp included in the filename to avoid overwriting work. 
    spse         Save a session file with a time stamp included in the filename to avoid overwriting work. 
    ssdf          Save sdf file with a time stamp included in the filename to avoid overwriting work. 
    swrl          Save wrl (VRML 2 file format) file with a time stamp included in the filename to avoid overwriting work. 

    Make molecular representations that are not available in PyMOL:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    AO            Commands to make ambient occlusion image like those in Qutemole.  
    AOD          Make ambient occlusion image of any with dark carbon atoms.      
    BW             Commands to make black-and white-ribbon cartoon on a white background. 
    BU             Commands to make biological unit. Requires a pdb file. There are 
    CB             Loads Jared Sampson's script "colorblindfriendly.py"   
    CR             Commands to make colored filled-ring cartoon of nucleic acids. May 
    CSS            Commands to color ribbon or cartoon representations of proteins by 
    CBSS          Apply colorblind-friendly coloring to ribbon or cartoon representations. 
    DU             Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins.  
    FR              Make filled-ring cartoon of nucleic acids. May need to enter 'hide everything' first.  
    HH             Hide hydrogen atoms of currently visible molecular objects.      
    PE              Apply pearl effect about cations. Must supply selection.         
    PU             Make putty cartoon of main chain of nucleic acids and proteins.  
    SE              Commands to make SAXS envelope from a bead model.    
    cav             Show buried cavities and pockets as molecular surfaces.            
    getchem     Create selections based on the biophysical properties of each residue. 
    timcolor      Use Tim Mather's coloring scheme applied to the selections defined in getchem().  
    
    Launch a full-featured text editor from PyMOL:
    (You may have to edit the path to the applicaiton.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    atom           Open file with the text editor Atom from within PyMOL.           
    bbedit         Open file with the text editor bbedit from within PyMOL.         
    code            Open file with Visual Studio Code from within PyMOL.             
    emacs          Open file with emacs from within PyMOL.                          
    gedit            Open file with gedit from within PyMOL.                          
    jedit             Open file with jedit from within PyMOL.                          
    mate            Open file with Textmate (Mac OS only) from within PyMOL.         
    notepadpp   Open file with notepadpp from within PyMOL.        
    nv                Open file with neovim from within PyMOL.                         
    oni               Open the editor Oni from within PyMOL.                           
    pdbed          Open PDBEditor.jar from within PyMOL.                            
    st3               Open sublime text 3 from within PyMOL.                           
    vim              Open vim from within PyMOL. 
                                   

    Open word processor:
    (You may have to edit the path to the Word applicaiton.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    word         Open word from within PyMOL.                                     
    (You may have to edit the path to the applicaiton in the corresponding function.) 
    
    Open data analysis programs:
    (You may have to edit the path to the applicaiton in the corresponding function.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    cranR         Open the Cran R from within PyMOL.                               
    ddb            Open the DBBrowserSQLite.                                        
    excel          Open the excel from within PyMOL.                                
    JASP           Open the JASP from within PyMOL.                                 
    JMP            Open the JMP from within PyMOL.                                  
    jabref         Open the jabref from within PyMOL.                               
    julia           Open the jabref from within PyMOL.                               
    oc              Open the jabref from within PyMOL.                               
    ppt             Open the powerpoint from within PyMOL.

    Open terminal windows:
    (You may have to edit the path to the applicaiton in the corresponding function.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    iterm     Open iTerm2 window on MacOS.                                     
    term      Open a Terminal window on MacOS. 
  
    Open other molecular graphics programs:
    (You may have to edit the path to the applicaiton in the corresponding function.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    ccp4mg    Open ccp4mg from within PyMOL.                                   
    chimera    Open Chimera from within PyMOL.                                  
    coot          Open coot from within PyMOL.                                     
    jmol          Open Jmol from within PyMOL.                                     
    vmd          Open vmd from within PyMOL.                                      
    yasara       Open the molecular graphics prograom YASASRA from within PyMOL.  


    Image manipulation programs:
    (You may have to edit the path to the applicaiton in the corresponding function.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    gimp         Open the molecular graphics program with gimp from within PyMOL.  
    inkscape    Open the molecular graphics program with gimp from within PyMOL.
   
    Open certain webapps:
    (You may have to edit the url in the corresponding function.)  
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    gcal            Open Google Calendar.                                            
    GM             Open gmail.                                                      
    WM             Open Web Mail in defualt browser. Adjust url for your institution. 
    WS              Open National Weather Service website for locale.                
      
    Samples:
    (You may have to edit the path to the applicaiton in the corresponding function.)
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    GGT          WT human gamma glutamyl transpeptidase at 1.67 Angstrom          
    GU            10-mer dsRNA with 8 contiguous Us. U-helix RNA.                  
    N9            Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr. 
    T4L           WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0.              
    U8            16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. 
    WC8         16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom
                  

    Commands to display complex scenes:
    Shortcuts  Description                                                          
    --------  ---------------------------------------------------------------
    BST           G2G3/U9U8 base step , PDB code 4PCO.                             
    LG             Nine sugar glycan in influenza N9 neuraminidase at 1.55 Ang.              
    NA            Hydrated sodium cation bound in major groove of a                

    Commands to display complex scenes with pdb files on computer.:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    LGGT         WT human gamma glutamyl transpeptidase at 1.67 Angstrom          
    LGU           10-mer dsRNA.                                                    
    LN9           Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 
    LT4L          Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0.  
    LU8           16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. 
    LWC8         16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom              
    LBST           G2G3/U9U8 base step , PDB code 4PCO.                             
    LLG            Nine sugar glycan in influenza N9 neuraminidase at               
    LNA            Hydrated sodium cation bound in major groove of a                
      
    Re-orient molecule:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    oy             Align long axis of molecule along z-axis.                        
    omxy        Align long axis of molecule along minus x-y axis.                
    oxy           Align long axis of molecule along x-y axis.                      
    oz             Align long axis of molecule along y-axis.                        

    Horizontal scripting:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    cntfiles      Count number of files in current directory.                      
    cntpdb      Count number of pdb files in current directory.                  
    rline          Enter "help(rline)" to refresh memory of the readline commands.  
    rv              Get the view settings in a compact format on one line.           

    Print commands for using git for version control:
    Shortcuts Description                                                          
    --------   ---------------------------------------------------------------
    gitAdd        Enter help(gitAdd) to print steps for adding a file to version control. 
    gitCommit  Enter help(gitInit) to print steps for saving updates to a file under version control. 
    gitInit         Enter help(gitInit) to print steps for creating a git repository. 
    gitPull         Enter help(gitPush) to print steps to send to updates to a repository on github.com.  
    gitPush       Enter help(gitPush) to print steps to send to updates to a repository on github.com.  

    Send search term(s) to websites with search boxes.:
    (You may have to edit the path to your default webbrowser.) 
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    AB             Send search term or phrase to Amazon.com Books in default browser. 
    GB             Send search term or phrase to Google Books in default browser.     
    GH             Send search term or phrase to GitHub in default browser.         
    GHN           Send search term or phrase to GitHub in default browser.         
    GO             Send search term or phrase Google in default browser.            
    GON           Send search term or phrase Google in default browser and opens the top N results in N new tabs. 
    GS              Send search term or phrase to Google Scholar in default browser. 
    GSN            Send search term or phrase to Google Scholar in default browser. 
    GV              Send search term or phrase to Google Videos in default browser.  
    GVN           Send search term or phrase to Google Videos in default browser.  
    MA             Send search term to all searchable websites in pymolshortcuts:   
    MB             Send search term to search multiple sites for term in books:     
    MC            Send search term to search ten core websites in pymolshortcuts:  
    MM            Send search term to search for manuscripts in pymolshortcuts:    
    PDB            Submit a search term to the Protein Data Bank.                   
    PDBN         Submit a search term to the Protein Data Bank and open the top N hits in separate tabs. 
    PML           Submit a search term to the PyMOL Users Mail Service.            
    PMLN         Submit a search term to the PyMOL Users Mail Service.            
    PM            Send search term or phrase to PubMed.                            
    PMN          Send search term or phrase to PubMed and open top N hits in separate tabs. 
    IPM           Read list of search terms and submit each term to PubMed in a separate browser tab. 
    IPMN         Read list of search terms and submit each term to PubMed in a separate browser tab. 
    RG            Submit a search query of Research Gate.                          
    RGN          Submit a search query of Research Gate and open the top N hits in sepearte tabs.  
    SD             Submit a search term to Science Direct.                          
    SDN          Submit a search term to Science Direct and open the top N hits in sepearte tabs. 
    SF             Send search term to sourceforge.                                 
    SFN           Send search term to sourceforge and open the top N hits in sepearte tabs. 
    SP             Submit a search term to Springer Books                           
    SPN           Submit a search term to Springer Books and open the top N hits in sepearte tabs.
             

    Open static web sites:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    ACA          Open the American Crystallographic Association Annual Meeting webpage. 
    ALS           Open website of the Advanced Light Source.                       
    APS           Open website of the Advanced Photon Source.                      
    AX             Send search term or phrase to arXiv.                             
    BC             Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source. 
    BD             Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB).  
    BX             Send search term or phrase to bioRxiv                            
    CH            Open the webste of UCSF Chimera.                                 
    CHESS        Open the website of CHESS.                                       
    EMDB         Open the website of the Electron Microscopy Data Bank.           
    EP              EasyPyMOL github site.                                           
    JM             Open the Jmol wiki.                                              
    IUCR          Open website of the IUCr Journals.                               
    LBSF          Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC. 
    MCL          Open website of Macromolecular Crystallography Laboratory at the University of Oklahoma.  
    MG            Open website of the OUHSC molecular graphics course.             
    NDB           Open website of the Nucleic Acid Database.                       
    notPyMOL  Open website with list of other molecular graphics programs.     
    NSLSII        Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory. 
    PPC            Open the website of the Protein Production Facility at the University of Oklahoma in Norman. 
    PS              Open the home page of the Protein Soceity.                       
    PW             Submit search of the PyMOL Wiki.                                 
    RS              Open the homepage of the RNA Society.                            
    SAXS          Open the webpage of SAXS links at OUHSC.                         
    SB              Open the webpage of SSRL Biological SAXS at BL 4-2.              
    SBGRID       Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel. 
    SciPy18      Open the SciPy 2018 YouTube Channel.                             
    SSRL           Open the webpage of SSRL Structural Molecular Biology.           
    SSURF         Open the webpage of the Society for Science at User Research Facilities (SSURF). 
    SO              Submit a search term to Stackoverflow. 
    (You may have to edit the path to your default webbrowser.)                          

    3D-PDFs:
    Shortcuts Description                                                          
    --------  ---------------------------------------------------------------
    ms2pdf    Send molecular surface or ribbon cartoon from PyMOL to 3dpdf.    
    topdf       Send stick models as pse file from PyMOL through Jmol to 3DPDF.  
    (You may have to edit the path to JMol.) 

    Type 'help <ShortCutName>' (e.g., help spng) for a description of the
    function and for two sets of commands. The first set of commands has
    line breaks for easy selection single lines of commands. The second
    set of commands is one one line for easy copying and pasting of the
    entire horizontal script. The commands can be copied from the
    command history window and pasted onto the command line for code
    reuse. Some aliases require additional scripts. 
    
    Type 'SC' to refresh the list of shortcuts.
    Type 'help rline' to see commands for moving cursor on the command line.



    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA
    Scroll up to see the list of shortcuts and their descriptions.

    PYTHON CODE:

def SC():
    print(SC.__doc__)

    '''

    print(SC.__doc__)

cmd.extend('SC',SC)

def SD(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Submit a search term to Science Direct.

    USAGE:

    SD

    Arguments:

    searchTerm

    EXAMPLE:

    SD session file

    MORE DETAILS:

    Submit a search term to Science Direct.

    Single term search (multi word searches do NOT have to be inside quotes):
    SD session file

    Multiple term search: 
    SD text editor; SD 3d pdf; SD black and white cartoon;


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SD(searchTerm="pymol"):
    url1 = 'https://www.sciencedirect.com/search/advanced?qs='
    url2 = '&show=100&sortBy=relevance'
    webbrowser.open(url1+searchTerm+url2)

    '''

    url1 = 'https://www.sciencedirect.com/search/advanced?qs='
    url2 = '&show=100&sortBy=relevance'
    webbrowser.open(url1+searchTerm+url2)

cmd.extend('SD',SD)

def SF(searchTerm='pymol'):
    ''' 
    DESCRIPTION:

    Send search term to sourceforge.

    USAGE:

    SF

    Arguments:

    searchTerm

    EXAMPLE:

    SF pymol

    MORE DETAILS:

    Send search term to sourceforge.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SF(searchTerm='pymol'):
    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)

    '''

    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)

cmd.extend('SF',SF)

def SO(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:

    Submit a search term to Stackoverflow.

    USAGE:

    SO

    Arguments:

    NA

    EXAMPLE:

    SO

    MORE DETAILS:

    Submit a search term to Stackoverflow.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SO(searchTerm="3d_pdf"):
    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)

    '''

    url = "https://stackoverflow.com/search?q="
    webbrowser.open(url+searchTerm)

cmd.extend('SO',SO)

def SP(searchTerm="pymol"):
    ''' 
    DESCRIPTION:

    Submit a search term to Springer Books

    USAGE:

    SP

    Arguments:

    searchTerm

    EXAMPLE:

    SP

    MORE DETAILS:

    Submit a search term to Springer Books.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SP(searchTerm="pymol"):
    url1 = 'https://www.springer.com/gp/search?query='
    url2 = '&submit=Submit+Query'
    webbrowser.open(url1+searchTerm+url2)

    '''

    url1 = 'https://www.springer.com/gp/search?query='
    url2 = '&submit=Submit+Query'
    webbrowser.open(url1+searchTerm+url2)

cmd.extend('SP',SP)

def SSRL():
    ''' 
    DESCRIPTION:

    Open the webpage of SSRL Structural Molecular Biology.

    USAGE:

    SSRL

    Arguments:

    NA

    EXAMPLE:

    SSRL

    MORE DETAILS:

    Open the webpage of SSRL Structural Molecular Biology.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SSRL():
    webbrowser.open('http://ssrl.slac.stanford.edu/smb/index.html')

    '''

    webbrowser.open('http://ssrl.slac.stanford.edu/smb/index.html')

cmd.extend('SSRL',SSRL)

def SSURF():
    ''' 
    DESCRIPTION:

    Open the webpage of the Society for Science at User Research Facilities (SSURF).


    USAGE:

    SSURF

    Arguments:

    NA

    EXAMPLE:

    SSURF

    MORE DETAILS:

    Open the webpage of the Society for Science at User Research Facilities (SSURF).
    SSURF is nonprofit organization in the US that serves as the
    nexus for users and user executive committees at
    national laboratories. SUURF is not a lobbying organization, but
    it helps organize visits to Congress to educate legislators about
    the importance of national laboratories in science. Membership
    is free of students. The annual fee is nominal for PIs. 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SSURF():
    webbrowser.open('http://www.ssurf.org/')
    '''

    webbrowser.open('http://www.ssurf.org/')
cmd.extend('SSURF',SSURF)

def SciPy19():
    ''' 
    DESCRIPTION:

    Open the SciPy 2019 YouTube Channel.

    USAGE:

    SciPy19

    Arguments:

    NA

    EXAMPLE:

    SciPy19

    MORE DETAILS:

    Open the SciPy 2019 YouTube Channel.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def SciPy19():
    webbrowser.open('https://www.scipy2019.scipy.org')
    '''

    webbrowser.open('https://www.scipy2019.scipy.org')
cmd.extend('SciPy19',SciPy19)

def T4L():
    ''' 
    DESCRIPTION:

    WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0. 

    USAGE:

    T4L

    Arguments:

    None

    EXAMPLE:

    T4L

    MORE DETAILS:

    WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0. 
    Type 'T4L' to activate. Type 'help T4L' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 3fa0,type=pdb,async=0;orient;turn z,-90;turn y,-5;turn x,10; hide everything; bg_color white;show cartoon;color red, ss H;color yellow, ss S;color green, ss L+;set_view (-0.18,-0.69,-0.7,0.98,-0.17,-0.09,-0.06,-0.7,0.71,0.0,0.0,-165.67,34.77,11.27,9.52,132.07,199.27,-20.0); ray 1500,1600; 


    PYTHON CODE:

def T4L():
    cmd.reinitialize()
    cmd.fetch('3fa0', type='pdb')
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

    '''

    cmd.reinitialize()
    cmd.fetch('3fa0', type='pdb')
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
    DESCRIPTION:

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.


    USAGE:

    U8

    Arguments:

    None

    EXAMPLE:

    U8

    MORE DETAILS:

    16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3.
    Has one strand in the asymmetric unit. Uses quat.py to generate
    the second strand. Cartoon with filled rings and bases cartoon.
    Type 'U8' to activate. Type 'help U8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.


   >>>  Edit file path to quat.py in Python code below.

    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all;fetch 3nd3,type=pdb,async=0;run $HOME/mg18OU/quat.py;quat 3nd3;hide everything;bg_color white; show sticks;set cartoon_ring_mode, 3;set cartoon_ring_finder, 1;set cartoon_ladder_mode, 1;set cartoon_nucleic_acid_mode, 4;set cartoon_ring_transparency, 0.5;as cartoon;set_view (-1.0,-0.03,0.06,-0.06,0.01,-1.0,0.04,-1.0,-0.01,-0.09,-0.02,-168.02,7.85,15.56,-0.21,137.38,199.33,-20.0);draw; 


    PYTHON CODE:

def U8():
    cmd.reinitialize()
    cmd.fetch('3nd3', type='pdb1')
    cmd.do('split_states 3nd3')
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

    '''

    cmd.reinitialize()
    cmd.fetch('3nd3', type='pdb1')
    cmd.do('split_states 3nd3')
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
    DESCRIPTION:

    16-mer dsRNA, Watson-Crick helix RNA: 3nd4.

    USAGE:

    WC8

    Arguments:

    None

    EXAMPLE:

    WC8

    MORE DETAILS:

    16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom 
    resolution: 3nd4.  Has one strand in the asymmetric unit. 
    Needs quat.py to generate the second strand. Use the 
    BU alias. Cartoon with filled rings and bases cartoon.

    Type 'WC8' to activate. Type 'help WC8' to see this documentation
    printed to the command history window. Select from the command
    history individual lines of code to build a new script. Select the
    hortizontal script at the bottom if retaining most of the commands
    in your new script. Copy and paste onto the command line below.
    Works only with the command line immediately under the command
    history window at the top of the gui.

   >>>  Edit file path to quat.py in Python code below.

    VERTICAL PML SCRIPT:

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


    HORIZONTAL PML SCRIPT:

    delete all; fetch 3nd4,type=pdb,async=0;hide everything; run $HOME/mg18OU/quat.py; quat 3nd4;bg_color white; show sticks; set stick_radius, 0.12; set nb_spheres_size, 0.25; show nb_spheres; set stick_ball, on; set stick_ball_ratio, 1.8;set_view (-0.99,-0.03,0.17,-0.18,0.02,-0.98,0.03,-1.0,-0.03,0.0,0.0,-169.97,8.1,15.62,-1.69,139.24,200.7,-20.0);hide everything, name H*;rock 

    PYTHON CODE:

def WC8():
    cmd.reinitialize()
    cmd.fetch('3nd4', type='pdb1')
    cmd.remove('name H*')
    cmd.hide('everything')
    cmd.do('split_states 3nd4')
    cmd.bg_color('white')
    cmd.do('show stick')
    cmd.do('set stick_radius, 0.12') 
    cmd.do('set nb_spheres_size, 0.25')
    cmd.do('show nb_spheres')
    cmd.do('set stick_ball, on')
    cmd.do('set stick_ball_ratio, 1.8')
    cmd.set_view('(-0.98,-0.03,0.22,-0.23,0.02,-0.97,0.03,-1.0,-0.03,0.0,0.0,-175.3,8.16,15.68,-1.66,144.53,206.07,-20.0)')
    cmd.rock()

    '''

    cmd.reinitialize()
    cmd.fetch('3nd4', type='pdb1')
    cmd.remove('name H*')
    cmd.hide('everything')
    cmd.do('split_states 3nd4')
    cmd.bg_color('white')
    cmd.do('show stick')
    cmd.do('set stick_radius, 0.12') 
    cmd.do('set nb_spheres_size, 0.25')
    cmd.do('show nb_spheres')
    cmd.do('set stick_ball, on')
    cmd.do('set stick_ball_ratio, 1.8')
    cmd.set_view('(-0.98,-0.03,0.22,-0.23,0.02,-0.97,0.03,-1.0,-0.03,0.0,0.0,-175.3,8.16,15.68,-1.66,144.53,206.07,-20.0)')
    cmd.rock()

cmd.extend('WC8',WC8)

def atom(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open the text editor Atom from within PyMOL. 

    USAGE:

    atom

    Arguments:

    None

    EXAMPLE:

    atom script.pml

    MORE DETAILS:

    Open the text editor Atom from within PyMOL. Adjust the path to Atom as needed for your system.
     
    >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

     Open file with the text editor Atom from within PyMOL. Adjust the path as needed for your system.

    HORIZONTAL PML SCRIPT:

    arg = ("/usr/local/bin/atom " + fileName)


    PYTHON CODE:

def atom(fileName="test.pml"):
    arg = ("/usr/local/bin/atom " + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = ("/usr/local/bin/atom " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('atom',atom)

def bbedit(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with the text editor bbedit from within PyMOL. 

    USAGE:

    bbedit

    Arguments:

    None

    EXAMPLE:

    bbedit script.pml

    MORE DETAILS:

    Open file with the text editor bbedit from within PyMOL. 
    Adjust the path as needed for your system.
    Only available for the Mac OS.

    >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (bbeditPath + fileName)
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = (bbeditPath + fileName); subprocess.call(arg,shell=True)


    PYTHON CODE:

def bbedit(fileName="test.pml"):
    arg = (bbeditPath + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (bbeditPath + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('bbedit',bbedit)

def bs(selection='all'):
    ''' 
    DESCRIPTION:

    BallnStick

    USAGE:

    bs selection

    Arguments:

    selection


    EXAMPLE:

    bs 3nd3

    MORE DETAILS:

    # BallnStick creates a ball and stick representation of an object 
# Bondi VDW values added below to override default Pymol settings
# From https://gist.githubusercontent.com/bobbypaton/1cdc4784f3fc8374467bae5eb410edef/raw/9995d51d6a64b8bcf01590c944cc38059b2f8d7f/pymol_style.py


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def bs(selection='all'):
    # Bondi VDW values 
    cmd.alter("elem Ac", "vdw=2.00")
    cmd.alter("elem Al", "vdw=2.00")
    cmd.alter("elem Am", "vdw=2.00")
    cmd.alter("elem Sb", "vdw=2.00")
    cmd.alter("elem Ar", "vdw=1.88")
    cmd.alter("elem As", "vdw=1.85")
    cmd.alter("elem At", "vdw=2.00")
    cmd.alter("elem Ba", "vdw=2.00")
    cmd.alter("elem Bk", "vdw=2.00")
    cmd.alter("elem Be", "vdw=2.00")
    cmd.alter("elem Bi", "vdw=2.00")
    cmd.alter("elem Bh", "vdw=2.00")
    cmd.alter("elem B ", "vdw=2.00")
    cmd.alter("elem Br", "vdw=1.85")
    cmd.alter("elem Cd", "vdw=1.58")
    cmd.alter("elem Cs", "vdw=2.00")
    cmd.alter("elem Ca", "vdw=2.00")
    cmd.alter("elem Cf", "vdw=2.00")
    cmd.alter("elem C ", "vdw=1.70")
    cmd.alter("elem Ce", "vdw=2.00")
    cmd.alter("elem Cl", "vdw=1.75")
    cmd.alter("elem Cr", "vdw=2.00")
    cmd.alter("elem Co", "vdw=2.00")
    cmd.alter("elem Cu", "vdw=1.40")
    cmd.alter("elem Cm", "vdw=2.00")
    cmd.alter("elem Ds", "vdw=2.00")
    cmd.alter("elem Db", "vdw=2.00")
    cmd.alter("elem Dy", "vdw=2.00")
    cmd.alter("elem Es", "vdw=2.00")
    cmd.alter("elem Er", "vdw=2.00")
    cmd.alter("elem Eu", "vdw=2.00")
    cmd.alter("elem Fm", "vdw=2.00")
    cmd.alter("elem F ", "vdw=1.47")
    cmd.alter("elem Fr", "vdw=2.00")
    cmd.alter("elem Gd", "vdw=2.00")
    cmd.alter("elem Ga", "vdw=1.87")
    cmd.alter("elem Ge", "vdw=2.00")
    cmd.alter("elem Au", "vdw=1.66")
    cmd.alter("elem Hf", "vdw=2.00")
    cmd.alter("elem Hs", "vdw=2.00")
    cmd.alter("elem He", "vdw=1.40")
    cmd.alter("elem Ho", "vdw=2.00")
    cmd.alter("elem In", "vdw=1.93")
    cmd.alter("elem I ", "vdw=1.98")
    cmd.alter("elem Ir", "vdw=2.00")
    cmd.alter("elem Fe", "vdw=2.00")
    cmd.alter("elem Kr", "vdw=2.02")
    cmd.alter("elem La", "vdw=2.00")
    cmd.alter("elem Lr", "vdw=2.00")
    cmd.alter("elem Pb", "vdw=2.02")
    cmd.alter("elem Li", "vdw=1.82")
    cmd.alter("elem Lu", "vdw=2.00")
    cmd.alter("elem Mg", "vdw=1.73")
    cmd.alter("elem Mn", "vdw=2.00")
    cmd.alter("elem Mt", "vdw=2.00")
    cmd.alter("elem Md", "vdw=2.00")
    cmd.alter("elem Hg", "vdw=1.55")
    cmd.alter("elem Mo", "vdw=2.00")
    cmd.alter("elem Nd", "vdw=2.00")
    cmd.alter("elem Ne", "vdw=1.54")
    cmd.alter("elem Np", "vdw=2.00")
    cmd.alter("elem Ni", "vdw=1.63")
    cmd.alter("elem Nb", "vdw=2.00")
    cmd.alter("elem N ", "vdw=1.55")
    cmd.alter("elem No", "vdw=2.00")
    cmd.alter("elem Os", "vdw=2.00")
    cmd.alter("elem O ", "vdw=1.52")
    cmd.alter("elem Pd", "vdw=1.63")
    cmd.alter("elem P ", "vdw=1.80")
    cmd.alter("elem Pt", "vdw=1.72")
    cmd.alter("elem Pu", "vdw=2.00")
    cmd.alter("elem Po", "vdw=2.00")
    cmd.alter("elem K ", "vdw=2.75")
    cmd.alter("elem Pr", "vdw=2.00")
    cmd.alter("elem Pm", "vdw=2.00")
    cmd.alter("elem Pa", "vdw=2.00")
    cmd.alter("elem Ra", "vdw=2.00")
    cmd.alter("elem Rn", "vdw=2.00")
    cmd.alter("elem Re", "vdw=2.00")
    cmd.alter("elem Rh", "vdw=2.00")
    cmd.alter("elem Rb", "vdw=2.00")
    cmd.alter("elem Ru", "vdw=2.00")
    cmd.alter("elem Rf", "vdw=2.00")
    cmd.alter("elem Sm", "vdw=2.00")
    cmd.alter("elem Sc", "vdw=2.00")
    cmd.alter("elem Sg", "vdw=2.00")
    cmd.alter("elem Se", "vdw=1.90")
    cmd.alter("elem Si", "vdw=2.10")
    cmd.alter("elem Ag", "vdw=1.72")
    cmd.alter("elem Na", "vdw=2.27")
    cmd.alter("elem Sr", "vdw=2.00")
    cmd.alter("elem S ", "vdw=1.80")
    cmd.alter("elem Ta", "vdw=2.00")
    cmd.alter("elem Tc", "vdw=2.00")
    cmd.alter("elem Te", "vdw=2.06")
    cmd.alter("elem Tb", "vdw=2.00")
    cmd.alter("elem Tl", "vdw=1.96")
    cmd.alter("elem Th", "vdw=2.00")
    cmd.alter("elem Tm", "vdw=2.00")
    cmd.alter("elem Sn", "vdw=2.17")
    cmd.alter("elem Ti", "vdw=2.00")
    cmd.alter("elem W ", "vdw=2.00")
    cmd.alter("elem U ", "vdw=1.86")
    cmd.alter("elem V ", "vdw=2.00")
    cmd.alter("elem Xe", "vdw=2.16")
    cmd.alter("elem Yb", "vdw=2.00")
    cmd.alter("elem Y ", "vdw=2.00")
    cmd.alter("elem Zn", "vdw=1.39")
    cmd.alter("elem Zr", "vdw=2.00")
    cmd.rebuild()

    # workspace settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 0)
    cmd.set("transparency", 0.5)
    cmd.set("dash_gap", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_texture", 2)
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", .1)
    cmd.space("cmyk")

    # defines BallnStick settings
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85","elem C and "+arg1)
    cmd.color("gray98","elem H and "+arg1)
    cmd.color("slate","elem N and "+arg1)
    cmd.set("stick_radius",0.07, arg1)
    cmd.set("sphere_scale",0.18, arg1)
    cmd.set("sphere_scale",0.13, arg1+" and elem H")
    cmd.set("dash_gap",0.01, arg1)
    cmd.set("dash_radius",0.07, arg1)
    cmd.set("stick_color","black", arg1)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.zoom(arg1)
    cmd.hide("labels")
    '''

    # Bondi VDW values 
    cmd.alter("elem Ac", "vdw=2.00")
    cmd.alter("elem Al", "vdw=2.00")
    cmd.alter("elem Am", "vdw=2.00")
    cmd.alter("elem Sb", "vdw=2.00")
    cmd.alter("elem Ar", "vdw=1.88")
    cmd.alter("elem As", "vdw=1.85")
    cmd.alter("elem At", "vdw=2.00")
    cmd.alter("elem Ba", "vdw=2.00")
    cmd.alter("elem Bk", "vdw=2.00")
    cmd.alter("elem Be", "vdw=2.00")
    cmd.alter("elem Bi", "vdw=2.00")
    cmd.alter("elem Bh", "vdw=2.00")
    cmd.alter("elem B ", "vdw=2.00")
    cmd.alter("elem Br", "vdw=1.85")
    cmd.alter("elem Cd", "vdw=1.58")
    cmd.alter("elem Cs", "vdw=2.00")
    cmd.alter("elem Ca", "vdw=2.00")
    cmd.alter("elem Cf", "vdw=2.00")
    cmd.alter("elem C ", "vdw=1.70")
    cmd.alter("elem Ce", "vdw=2.00")
    cmd.alter("elem Cl", "vdw=1.75")
    cmd.alter("elem Cr", "vdw=2.00")
    cmd.alter("elem Co", "vdw=2.00")
    cmd.alter("elem Cu", "vdw=1.40")
    cmd.alter("elem Cm", "vdw=2.00")
    cmd.alter("elem Ds", "vdw=2.00")
    cmd.alter("elem Db", "vdw=2.00")
    cmd.alter("elem Dy", "vdw=2.00")
    cmd.alter("elem Es", "vdw=2.00")
    cmd.alter("elem Er", "vdw=2.00")
    cmd.alter("elem Eu", "vdw=2.00")
    cmd.alter("elem Fm", "vdw=2.00")
    cmd.alter("elem F ", "vdw=1.47")
    cmd.alter("elem Fr", "vdw=2.00")
    cmd.alter("elem Gd", "vdw=2.00")
    cmd.alter("elem Ga", "vdw=1.87")
    cmd.alter("elem Ge", "vdw=2.00")
    cmd.alter("elem Au", "vdw=1.66")
    cmd.alter("elem Hf", "vdw=2.00")
    cmd.alter("elem Hs", "vdw=2.00")
    cmd.alter("elem He", "vdw=1.40")
    cmd.alter("elem Ho", "vdw=2.00")
    cmd.alter("elem In", "vdw=1.93")
    cmd.alter("elem I ", "vdw=1.98")
    cmd.alter("elem Ir", "vdw=2.00")
    cmd.alter("elem Fe", "vdw=2.00")
    cmd.alter("elem Kr", "vdw=2.02")
    cmd.alter("elem La", "vdw=2.00")
    cmd.alter("elem Lr", "vdw=2.00")
    cmd.alter("elem Pb", "vdw=2.02")
    cmd.alter("elem Li", "vdw=1.82")
    cmd.alter("elem Lu", "vdw=2.00")
    cmd.alter("elem Mg", "vdw=1.73")
    cmd.alter("elem Mn", "vdw=2.00")
    cmd.alter("elem Mt", "vdw=2.00")
    cmd.alter("elem Md", "vdw=2.00")
    cmd.alter("elem Hg", "vdw=1.55")
    cmd.alter("elem Mo", "vdw=2.00")
    cmd.alter("elem Nd", "vdw=2.00")
    cmd.alter("elem Ne", "vdw=1.54")
    cmd.alter("elem Np", "vdw=2.00")
    cmd.alter("elem Ni", "vdw=1.63")
    cmd.alter("elem Nb", "vdw=2.00")
    cmd.alter("elem N ", "vdw=1.55")
    cmd.alter("elem No", "vdw=2.00")
    cmd.alter("elem Os", "vdw=2.00")
    cmd.alter("elem O ", "vdw=1.52")
    cmd.alter("elem Pd", "vdw=1.63")
    cmd.alter("elem P ", "vdw=1.80")
    cmd.alter("elem Pt", "vdw=1.72")
    cmd.alter("elem Pu", "vdw=2.00")
    cmd.alter("elem Po", "vdw=2.00")
    cmd.alter("elem K ", "vdw=2.75")
    cmd.alter("elem Pr", "vdw=2.00")
    cmd.alter("elem Pm", "vdw=2.00")
    cmd.alter("elem Pa", "vdw=2.00")
    cmd.alter("elem Ra", "vdw=2.00")
    cmd.alter("elem Rn", "vdw=2.00")
    cmd.alter("elem Re", "vdw=2.00")
    cmd.alter("elem Rh", "vdw=2.00")
    cmd.alter("elem Rb", "vdw=2.00")
    cmd.alter("elem Ru", "vdw=2.00")
    cmd.alter("elem Rf", "vdw=2.00")
    cmd.alter("elem Sm", "vdw=2.00")
    cmd.alter("elem Sc", "vdw=2.00")
    cmd.alter("elem Sg", "vdw=2.00")
    cmd.alter("elem Se", "vdw=1.90")
    cmd.alter("elem Si", "vdw=2.10")
    cmd.alter("elem Ag", "vdw=1.72")
    cmd.alter("elem Na", "vdw=2.27")
    cmd.alter("elem Sr", "vdw=2.00")
    cmd.alter("elem S ", "vdw=1.80")
    cmd.alter("elem Ta", "vdw=2.00")
    cmd.alter("elem Tc", "vdw=2.00")
    cmd.alter("elem Te", "vdw=2.06")
    cmd.alter("elem Tb", "vdw=2.00")
    cmd.alter("elem Tl", "vdw=1.96")
    cmd.alter("elem Th", "vdw=2.00")
    cmd.alter("elem Tm", "vdw=2.00")
    cmd.alter("elem Sn", "vdw=2.17")
    cmd.alter("elem Ti", "vdw=2.00")
    cmd.alter("elem W ", "vdw=2.00")
    cmd.alter("elem U ", "vdw=1.86")
    cmd.alter("elem V ", "vdw=2.00")
    cmd.alter("elem Xe", "vdw=2.16")
    cmd.alter("elem Yb", "vdw=2.00")
    cmd.alter("elem Y ", "vdw=2.00")
    cmd.alter("elem Zn", "vdw=1.39")
    cmd.alter("elem Zr", "vdw=2.00")
    cmd.rebuild()

    # workspace settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 0)
    cmd.set("transparency", 0.5)
    cmd.set("dash_gap", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_texture", 2)
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", .1)
    cmd.space("cmyk")

    # defines BallnStick settings
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85","elem C and "+arg1)
    cmd.color("gray98","elem H and "+arg1)
    cmd.color("slate","elem N and "+arg1)
    cmd.set("stick_radius",0.07, arg1)
    cmd.set("sphere_scale",0.18, arg1)
    cmd.set("sphere_scale",0.13, arg1+" and elem H")
    cmd.set("dash_gap",0.01, arg1)
    cmd.set("dash_radius",0.07, arg1)
    cmd.set("stick_color","black", arg1)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.zoom(arg1)
    cmd.hide("labels")
cmd.extend('bs', bs)

def ccp4mg(fileName="test.pdb"):
    ''' 
    DESCRIPTION:

    Open ccp4mg from within PyMOL. 



    USAGE:

    ccp4mg

    Arguments:

    None

    EXAMPLE:

    ccp4mg

    MORE DETAILS:

    Open ccp4mg from within PyMOL. 
    Adjust url for your location.


    VERTICAL PML SCRIPT:

    arg = ( ccp4mgCommand + fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (ccp4mgCommand + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def ccp4mg(fileName="test.pdb"):
    arg = (ccp4mgCommand + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (ccp4mgCommand + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('ccp4mg',ccp4mg)

def chimera(fileName="test.pdb"):
    ''' 
    DESCRIPTION:

    Open Chimera from within PyMOL. 
 

    USAGE:

    chimera

    Arguments:

    None

    EXAMPLE:

    chimera

    MORE DETAILS:

    Open Chimera from within PyMOL. 
    Adjust url for your location.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (chimeriaPath + fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (chimeriaPath + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def chimera(fileName="test.pdb"):
    arg = (chimeriaPath+ fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (chimeriaPath+ fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('chimera',chimera)

def cntfiles():
    ''' 
    DESCRIPTION:

    Count number of files in current directory.

    USAGE:

    cntfiles

    Arguments:

    None

    EXAMPLE:

    cntfiles

    MORE DETAILS:

    Count number of files in current directory.

    VERTICAL PML SCRIPT:

    arg = 'echo "Count the files in the directory." && echo "Usage: cntfiles." && find . -type f | wc -l'
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:

    arg = 'echo "Count the files in the directory." && echo "Usage: cntfiles." && find . -type f | wc -l';subprocess.call(arg,shell=True);return

    PYTHON CODE:

def cntfiles():
    arg = 'echo "Count the files in the directory." && echo "Usage: cntfiles." && find . -type f | wc -l';
    subprocess.call(arg,shell=True);
    return
    '''

    arg = 'echo "Count the files in the directory." && echo "Usage: cntfiles." && find . -type f | wc -l';
    subprocess.call(arg,shell=True);
    return
cmd.extend('cntfiles',cntfiles)

def cntpdb():
    ''' 
    DESCRIPTION:

    Count number of pdb files in current directory.

    USAGE:

    cntpdb

    Arguments:

    None

    EXAMPLE:

    cntpdb

    MORE DETAILS:

    Count number of pdb files in current directory.

    VERTICAL PML SCRIPT:

    arg = "echo 'Count the number of pdb files in subfolders.' && echo 'Usage: cntpdb' && find ./ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e '{}/*_([1-9999]|****).pdb' ';' -print | wc -l"
    subprocess.call(arg,shell=True)
    return


    HORIZONTAL PML SCRIPT:

    arg = "echo 'Count the number of pdb files in subfolders.' && echo 'Usage: cntpdb' && find ./ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e '{}/*_([1-9999]|****).pdb' ';' -print | wc -l"
subprocess.call(arg,shell=True);return


    PYTHON CODE:

def cntpdb():
    arg = "echo 'Count the number of pdb files in subfolders.' && echo 'Usage: cntpdb' && find ./ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e '{}/*_([1-9999]|****).pdb' ';' -print | wc -l"
    subprocess.call(arg,shell=True)
    return

    '''

    arg = "echo 'Count the number of pdb files in subfolders.' && echo 'Usage: cntpdb' && find ./ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e '{}/*_([1-9999]|****).pdb' ';' -print | wc -l"
    subprocess.call(arg,shell=True)
    return

cmd.extend('cntpdb',cntpdb)

def code(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with Visual Studio Code from within PyMOL.

    USAGE:

    code

    Arguments:

    None

    EXAMPLE:

    code script.pml

    MORE DETAILS:

    Open file with Visual Studio Code from within PyMOL. 
    Install the bioSyntax extension (free) from the Visual Studio Marketplace
    to get color syntax highlighting of pml files along with Fasta and other
    sequence files and PDB coordinate files. 

    https://marketplace.visualstudio.com/items?itemName=reageyao.biosyntax

    >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (codePath + fileName);
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = {codePath + fileName);subprocess.call(arg,shell=True)


    PYTHON CODE:

def code(fileName="test.pml"):
    arg = (codePath + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (codePath + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('code',code)

def coot(fileName="test.pdb"):
    ''' 
    DESCRIPTION:

    Open coot from within PyMOL. 



    USAGE:

    coot

    Arguments:

    None

    EXAMPLE:

    coot

    MORE DETAILS:

    Open coot from within PyMOL. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (cootPath + fileName)
    subprocess.call(arg,shell=True)
    return


    HORIZONTAL PML SCRIPT:

    arg = (cootPath  + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def coot(fileName="test.pdb"):
    arg = (cootPath + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (cootPath + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('coot',coot)

def cranR():
    ''' 
    DESCRIPTION:

    Open Cran R from within PyMOL.

    USAGE:

    cranR

    Arguments:

    None

    EXAMPLE:

    cranR script.R

    MORE DETAILS:

    Open Cran R from within PyMOL.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (RPath);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (RPath);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def cranR():
    arg = (RPath)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (RPath)
    subprocess.call(arg,shell=True)
    return

cmd.extend('cranR',cranR)

def ddb():
    ''' 
    DESCRIPTION:

    Open DBBrowserSQLite. 

    USAGE:

    ddb

    Arguments:

    None

    EXAMPLE:

    ddb database.db

    MORE DETAILS:

    Open DBBrowserSQLite. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = dbbrowserPath;
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = dbbrowserPath;subprocess.call(arg,shell=True);return

    PYTHON CODE:

def ddb():
    arg = dbbrowserPath;
    subprocess.call(arg,shell=True);
    return
    '''

    arg = dbbrowserPath;
    subprocess.call(arg,shell=True);
    return
cmd.extend('ddb',ddb)

def emacs(fileName="testme.pml"):
    ''' 
    DESCRIPTION:

    Open file with emacs from within PyMOL.

    USAGE:

    emacs

    Arguments:

    None

    EXAMPLE:

    emacs script.pml

    MORE DETAILS:

    Open file with emacs from within PyMOL. 
    Adjust path to emacs on your computer as needed.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg2 = ('--file ' + fileName);
    subprocess.call(emacsCommand)


    HORIZONTAL PML SCRIPT:

    arg2 = ('--file ' + fileName);subprocess.call(emacsCommand)

    PYTHON CODE:

def emacs(fileName="testme.pml"):
    arg = (emacsPath+ fileName)
    arg2 = ('--file ' + fileName)
    subprocess.call(emacsCommand)
    return
    '''

    arg = (emacsPath+ fileName)
    arg2 = ('--file ' + fileName)
    subprocess.call(emacsCommand)
    return
cmd.extend('emacs',emacs)

def excel():
    ''' 
    DESCRIPTION:

    Open excel from within PyMOL. 

    USAGE:

    excel

    Arguments:

    None

    EXAMPLE:

    excel

    MORE DETAILS:

    Open excel from within PyMOL. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = excelCommand
    subprocess.call(arg,shell=True)
    return


    HORIZONTAL PML SCRIPT:

    arg = excelCommand;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def excel():
    arg = excelCommand
    subprocess.call(arg,shell=True)
    return

    '''

    arg = excelCommand
    subprocess.call(arg,shell=True)
    return

cmd.extend('excel',excel)

def gcal():
    ''' 
    DESCRIPTION:

    Open Google Calendar.

    USAGE:

    gcal

    Arguments:

    None

    EXAMPLE:

    gcal

    MORE DETAILS:

    Open Google Calendar.

   >>>  Edit url in python code below

    VERTICAL PML SCRIPT:

    webbrowser.open('https://calendar.google.com/calendar/r')

    HORIZONTAL PML SCRIPT:

    webbrowser.open('https://calendar.google.com/calendar/r')

    PYTHON CODE:

def gcal():
    webbrowser.open('https://calendar.google.com/calendar/r')
    '''

    webbrowser.open('https://calendar.google.com/calendar/r')
cmd.extend('gcal',gcal)

def gedit(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with gedit from within PyMOL.

    USAGE:

    gedit

    Arguments:

    None

    EXAMPLE:

    gedit script.pml

    MORE DETAILS:

    Open file with gedit from within PyMOL. 
    Adjust the filepath for location of your executable.
    Can be installed via macports on the mac.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (geditPath +" -w " + fileName);
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = ( geditPath + " -w " + fileName);    subprocess.call(arg,shell=True)


    PYTHON CODE:

def gedit(fileName="test.pml"):
    arg = (geditPath + " -w " + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (geditPath + " -w " + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('gedit',gedit)

def gitAdd():
    ''' 
    DESCRIPTION:

    Enter help(gitAdd) to print steps for adding a file to version control.

    USAGE:

    gitAdd

    Arguments:

    None

    EXAMPLE:

    gitAdd

    MORE DETAILS:

    Enter help(gitAdd) to print steps for adding a file to version control.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def gitAdd():
    #NA
    '''

    #NA
cmd.extend("gitAdd",gitAdd)

def gitCommit():
    ''' 
    DESCRIPTION:

    Enter help(gitCommit) to print steps for saving updates to a file under version control.


    USAGE:

    gitCommit

    Arguments:

    None

    EXAMPLE:

    gitCommit

    MORE DETAILS:

    Enter help(gitInit) to print steps for saving updates to a file under version control.

    git add filename
    git commit -m "Message about changes"


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def gitCommit():
    #NA
    '''

    #NA
cmd.extend("gitCommit",gitCommit)

def gitInit():
    ''' 
    DESCRIPTION:

    Enter help(gitInit) to print steps for creating a git repository.

    USAGE:

    gitInit

    Arguments:

    None

    EXAMPLE:

    gitInit

    MORE DETAILS:

    Enter help(gitInit) to print steps for creating a git repository.
        Plain text files can be put under version control.
        Binary and bitmap files should not be put under version control.

      Step 1: initialize repository in current directory
        git init 

      Step 2: make list of binary file types to be ignored (e.g. 
    *.pse
    *.ccp4
    *.mtz
    *.png
    *.dmg
    *.pdf
    *.tiff
    *.jpeg
    *.jpg
    *.zip
    *.tar
    *.rar
    *.jar
    *.iso
    *.7z
    *.o
    *.so
    *.pyc
    *.doc
    *.docx
    *.idtf
    *.u3d
    *.aux
     )
        touch .gitignore
        git add .gitignore
        git commit -m "message" .gitignore

      Step 3: add all other files to repository
        git add .
      Later add new files one at a time
        git add new.pml

      Step 4: commit new files or changes with message
        git commit -m "Edited new.pml."

        or to commit changes to a group of changed files
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
        git commit -m "Added new text in $1.pml" 
        }


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def gitInit():
    print(gitInit.__doc__)
    '''

    print(gitInit.__doc__)
cmd.extend("gitInit",gitInit)

def gitPull():
    ''' 
    DESCRIPTION:

    Enter help(gitPush) to print steps to supdate a repository on github.com.

    USAGE:

    gitPull

    Arguments:

    None

    EXAMPLE:

    gitPull

    MORE DETAILS:

    Enter help(gitPull) to print steps to send to updates to a repository on github.com.
      Step 1: pull file from an existing repository. 

        git pull 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def gitPull():
    #NA
    '''

    #NA
cmd.extend("gitPull",gitPull)

def gitPush():
    ''' 
    DESCRIPTION:

    Enter help(gitPush) to print steps update a repository on github.com. 


    USAGE:

    gitPush

    Arguments:

    None

    EXAMPLE:

    gitPush

    MORE DETAILS:

    Enter help(gitPush) to print steps to send to updates to a repository on github.com. 

      Step 1: push updated file to an existing repository. 

        git push 


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def gitPush():
    #NA
    '''

    #NA
cmd.extend("gitPush",gitPush)

def GM():
    ''' 
    DESCRIPTION:

    Open gmail. 

    USAGE:

    gmail

    Arguments:

    None

    EXAMPLE:

    gmail

    MORE DETAILS:

    Open gmail. Edit url in python code below

    VERTICAL PML SCRIPT:

    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')

    HORIZONTAL PML SCRIPT:

    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')

    PYTHON CODE:

def GM():
    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')
    '''

    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')
cmd.extend('GM',GM)

def iterm():
    ''' 
    DESCRIPTION:

    Open iTerm2 window on MacOS. 

    USAGE:

    iterm

    Arguments:

    None

    EXAMPLE:

    iterm

    MORE DETAILS:

    Open iTerm2 window on MacOS. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    subprocess.call(itermCommand);
    return


    HORIZONTAL PML SCRIPT:

    subprocess.call(itermPath);return


    PYTHON CODE:

def iterm():
    subprocess.call(itermPath)
    return

    '''

    subprocess.call(itermPath)
    return

cmd.extend('iterm',iterm)

def jabref():
    ''' 
    DESCRIPTION:

    Open the jabref from within PyMOL.

    USAGE:

    jabref

    Arguments:

    None

    EXAMPLE:

    jabref

    MORE DETAILS:

    Open the jabref from within PyMOL.

   >>>  Edit file path in python code below.

    VERTICAL PML SCRIPT:

    arg = (jabrefPath);
    subprocess.call(arg,shell=True);
    return;


    HORIZONTAL PML SCRIPT:

    arg = jabrefCommand; subprocess.call(arg,shell=True);return


    PYTHON CODE:

def jabref():
    arg = jabrefCommand;
    subprocess.call(arg,shell=True);
    return

    '''

    arg = jabrefCommand;
    subprocess.call(arg,shell=True);
    return

cmd.extend('jabref',jabref)

def jedit(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with jedit from within PyMOL. 

    USAGE:

    jedit

    Arguments:

    None

    EXAMPLE:

    jedit script.pml

    MORE DETAILS:

    Open file with jedit from within PyMOL. 
    Adjust file path to your executable.
    Can be installed via macports on the Mac.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (jeditCommand+ fileName);
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = (jeditCommand + fileName);    subprocess.call(arg,shell=True)


    PYTHON CODE:

def jedit(fileName="test.pml"):
    arg = (jeditCommand + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (jeditCommand + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('jedit',jedit)

def jmol(fileName="test.pdb"):
    ''' 
    DESCRIPTION:

    Open Jmol from within PyMOL.

    USAGE:

    jmol

    Arguments:

    None

    EXAMPLE:

    jmol

    MORE DETAILS:

    Open Jmol from within PyMOL.
    Adjust file path to your location of Jmol.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (jmolPath + fileName)
    subprocess.call(arg,shell=True)
    return


    HORIZONTAL PML SCRIPT:

    arg = (jmolPath + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def jmol(fileName="test.pdb"):
    arg = (jmolPath + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (jmolPath + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('jmol',jmol)

def julia():
    ''' 
    DESCRIPTION:

    Open the julia from within PyMOL.

    USAGE:

    julia

    Arguments:

    None

    EXAMPLE:

    julia

    MORE DETAILS:

    Open the julia from within PyMOL.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = juliaPath;
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg =juliaPath;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def julia():
    arg = juliaPath;
    subprocess.call(arg,shell=True);
    return

    '''

    arg = juliaPath;
    subprocess.call(arg,shell=True);
    return

cmd.extend('julia',julia)

def mate(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with mate from within PyMOL. 

    USAGE:

    mate

    Arguments:

    None

    EXAMPLE:

    mate script.pml

    MORE DETAILS:

    Open file with Textmate (Mac OS only) from within PyMOL. 
    Adjust path to Textmate on your computer as needed.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (matePath + fileName);
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = (matePath + fileName);    subprocess.call(arg,shell=True)


    PYTHON CODE:

def mate(fileName="test.pml"):
    arg = (matePath + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (matePath + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('mate',mate)

def nmr():

    ''' 
    DESCRIPTION:

    Show all models in a nmr structure. 

    USAGE:

    nmr

    Arguments:

    None

    EXAMPLE:

    nmr

    MORE DETAILS:

    Show all models in a nmr structure. 

    VERTICAL PML SCRIPT:

    set all_states, on

    HORIZONTAL PML SCRIPT:

    set all_states, on

    PYTHON CODE:

def nmr():

    cmd.do('set all_states, on')
    '''

    cmd.do('set all_states, on')
cmd.extend("nmr", nmr)


def nmroff():

    ''' 
    DESCRIPTION:

    Hide all but first model in a nmr structure. 

    USAGE:

    nmroff

    Arguments:

    None

    EXAMPLE:

    nmroff

    MORE DETAILS:

    Hide all but first model in a nmr structure. 

    VERTICAL PML SCRIPT:

    set all_states, off

    HORIZONTAL PML SCRIPT:

    set all_states, off

    PYTHON CODE:

def nmroff():

    cmd.do('set all_states, off')
    '''

    cmd.do('set all_states, off')
cmd.extend("nmr", nmr)


def notPyMOL():
    ''' 
    DESCRIPTION:

    Open website with list of other molecular graphics programs.

    USAGE:

    notPyMOL

    Arguments:

    NA

    EXAMPLE:

    notPyMOL

    MORE DETAILS:

    Open website with list of other molecular graphics programs.

    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def notPyMOL():
    webbrowser.open('https://www.bnl.gov/ps/')
    '''

    webbrowser.open('https://www.bnl.gov/ps/')
cmd.extend('notPyMOL',notPyMOL)

def notepadpp(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open file with notepadpp from within PyMOL. 

    USAGE:

    notepadpp

    Arguments:

    None

    EXAMPLE:

    notepadpp script.pml

    MORE DETAILS:

    Open file with notepadpp from within PyMOL. 
    Adjust path to notepadpp on your computer as needed.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (notepadppPath + fileName);
    subprocess.call(arg,shell=True)


    HORIZONTAL PML SCRIPT:

    arg = (notepadppPath + fileName);    subprocess.call(arg,shell=True)


    PYTHON CODE:

def notepadpp(fileName="test.pml"):
    arg = (notepadppPath + fileName)
    subprocess.call(arg,shell=True)
    return
    '''

    arg = (notepadppPath + fileName)
    subprocess.call(arg,shell=True)
    return
cmd.extend('notepadpp',notepadpp)

def nv(fileName="testme.pml"):
    ''' 
    DESCRIPTION:

    Open file with neovim from within PyMOL.

    USAGE:

    nv

    Arguments:

    None

    EXAMPLE:

    nv script.pml

    MORE DETAILS:

    Open file with neovim from within PyMOL. 
    Adjust path to neovim on your computer as needed.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    subprocess.call(nvimCommand)


    HORIZONTAL PML SCRIPT:

        subprocess.call(nvimCommand)


    PYTHON CODE:

def nv(fileName="testme.pml"):
    subprocess.call(nvimCommand)

    '''

    subprocess.call(nvimCommand)

cmd.extend('nv',nv)

def oc():
    ''' 
    DESCRIPTION:

    Open the octave from within PyMOL.

    USAGE:

    octave

    Arguments:

    None

    EXAMPLE:

    oc

    MORE DETAILS:

    Open the octave from within PyMOL.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = octaveCommand;
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = octaveCommand;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def oc():
    arg = octaveCommand
    subprocess.call(arg,shell=True)
    return

    '''

    arg = octaveCommand
    subprocess.call(arg,shell=True)
    return

cmd.extend('oc',oc)

def oni(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open the editor Oni from within PyMOL. 

    USAGE:

    oni

    Arguments:

    None

    EXAMPLE:

    oni script.pml

    MORE DETAILS:

    Open the editor Oni from within PyMOL. 
    The is an editor based on neovim.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (oniPath + fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (oniPath + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def oni(fileName="test.pml"):
    arg = (oniPath + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (oniPath + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('oni',oni)

def ppt():
    ''' 
    DESCRIPTION:

    Open the powerpoint from within PyMOL. 

    USAGE:

    ppt

    Arguments:

    None

    EXAMPLE:

    ppt

    MORE DETAILS:

    Open the powerpoint from within PyMOL. 

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = powerpointCommand
    subprocess.call(arg,shell=True)
    return


    HORIZONTAL PML SCRIPT:

    arg = powerpointCommand;subprocess.call(arg,shell=True);return


    PYTHON CODE:

def ppt():
    arg =powerpointCommand
    subprocess.call(arg,shell=True)
    return

    '''

    arg =powerpointCommand
    subprocess.call(arg,shell=True)
    return

cmd.extend('ptt',ppt)

def rline():
    ''' 
    DESCRIPTION:

    Prints cheat sheet for the readline commands. 


    USAGE:

    rline

    Arguments:

    None

    EXAMPLE:

    rline

    MORE DETAILS:

    Prints cheat sheet for the readline commands.

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


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def rline():
    print(readline.__doc__)
    '''

    print(readline.__doc__)
cmd.extend("rline",rline)

def rmd():
    ''' 
    DESCRIPTION:

    Remove all measurement objects in the interal GUI.

    USAGE:

    rmd

    Arguments:

    None

    EXAMPLE:

    rmd

    MORE DETAILS:

    Remove all measurement objects in the interal GUI.
    Note that there is a "delete all measurements" toggle in the internal gui.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    delete measure*;
    delete m*_*

    HORIZONTAL PML SCRIPT:

    delete measure*; delete m*_*

    PYTHON CODE:

def rmd():
    cmd.do('delete measure*')
    cmd.do('delete m*_*')
    '''

    cmd.do('delete measure*')
    cmd.do('delete m*_*')
cmd.extend("rmd", rmd)

def rmsc():

    ''' 
    DESCRIPTION:

    Remove supercell objects.

    USAGE:

    rmsc

    Arguments:

    None

    EXAMPLE:

    rmsc

    MORE DETAILS:

    Use 'rmsc' to remove supercell objects.


    VERTICAL PML SCRIPT:

    delete supercell*;delete m*_*

    HORIZONTAL PML SCRIPT:

    delete supercell*;delete m*_*

    PYTHON CODE:

def rmsc():

    cmd.do('delete supercell*;delete m*_*')
    '''

    cmd.do('delete supercell*;delete m*_*')
cmd.extend("rmsc", rmsc)


def rv(StoredView=0, decimal_places=2, outname="roundedview.txt"):
    ''' 
    DESCRIPTION:

    Get the view settings in a compact format on one line.


    USAGE:

    rv

    Arguments:

    [view, decimal_places, outname]
Note that the values in the [] are optional.


    EXAMPLE:

    rv

    MORE DETAILS:

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


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def rv(StoredView=0, decimal_places=2, outname="roundedview.txt"):
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
    print(x.format(*myRoundedList))

    #print 'set_view ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},\
    #{10},{11},{12},{13},{14},{15},{16},{17})'.format(*myRoundedList)

    #Write to a text file.
    myFile = open("roundedview.txt", "a")
    myFile.write(x.format(*myRoundedList) + "\n")
    myFile.close()
    return

    '''

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
    print(x.format(*myRoundedList))

    #print 'set_view ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},\
    #{10},{11},{12},{13},{14},{15},{16},{17})'.format(*myRoundedList)

    #Write to a text file.
    myFile = open("roundedview.txt", "a")
    myFile.write(x.format(*myRoundedList) + "\n")
    myFile.close()
    return

cmd.extend("rv", rv)

def saln(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a aln file (alignment file) with a time stamp.

    USAGE:

    saln

    Arguments:

    alignmentFileNameStem

    EXAMPLE:

    saln  testFile

    MORE DETAILS:

    Save a aln file (alignment file) with a time stamp included in the filename to avoid
    overwriting work. Read as a commandline argument, a string as the filename 
    stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");
    s = str(DT);
    cmd.save(stemName+s+".aln")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".aln")

    PYTHON CODE:

def saln(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".aln") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".aln") 
cmd.extend('saln',saln)

def sc111():

    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 1 x 1 unit cells. 

    USAGE:

    sc111

    Arguments:

    None

    EXAMPLE:

    sc111

    MORE DETAILS:

    Make a lattice of 1 x 1 x 1 unit cells.
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\nsupercell 1, 1, 1, , orange, supercell111, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;supercell 1, 1, 1, , orange, supercell111, 1

    PYTHON CODE:

def sc111():

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 1, , orange, supercell111, 1')

    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 1, , orange, supercell111, 1')

cmd.extend("sc111", sc111)


def sc112():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 1 x 2 unit cells

    USAGE:

    sc112

    Arguments:

    None

    EXAMPLE:

    sc112

    MORE DETAILS:

    Make a lattice of 1 x 1 x 2 unit cells.
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.

     >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 1, 2, , orange, supercell112, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;supercell 1, 1, 2, , orange, supercell112, 1

    PYTHON CODE:

def sc112():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 2, , orange, supercell112, 1')

    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 2, , orange, supercell112, 1')

cmd.extend("sc112", sc112)


def sc113():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 1 x 3 unit cells.

    USAGE:

    sc113

    Arguments:

    None

    EXAMPLE:

    sc113

    MORE DETAILS:

    Make a lattice of 1 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 1, 3, , orange, supercell113, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 1, 3, , orange, supercell113, 1

    PYTHON CODE:

def sc113():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 3, , orange, supercell113, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 1, 3, , orange, supercell113, 1')
cmd.extend("sc113", sc113)


def sc121():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 2 x 1 unit cells. 

    USAGE:

    sc121

    Arguments:

    None

    EXAMPLE:

    sc121

    MORE DETAILS:

    Make a lattice of 1 x 2 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 2, 1, , orange, supercell121, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 1, , orange, supercell121, 1

    PYTHON CODE:

def sc121():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 1, , orange, supercell121, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 1, , orange, supercell121, 1')
cmd.extend("sc121", sc121)


def sc122():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 2 x 2 unit cells. 

    USAGE:

    sc122

    Arguments:

    None

    EXAMPLE:

    sc122

    MORE DETAILS:

    Make a lattice of 1 x 2 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 2, 2, , orange, supercell122, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 2, , orange, supercell122, 1

    PYTHON CODE:

def sc122():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 2, , orange, supercell122, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 2, , orange, supercell122, 1')
cmd.extend("sc122", sc122)


def sc123():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 2 x 3 unit cells. 

    USAGE:

    sc123

    Arguments:

    None

    EXAMPLE:

    sc123

    MORE DETAILS:

    Make a lattice of 1 x 2 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 3, , orange, supercell123, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 3, , orange, supercell123, 1

    PYTHON CODE:

def sc123():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 3, , orange, supercell123, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 2, 3, , orange, supercell123, 1')
cmd.extend("sc123", sc123)


def sc131():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 3 x 1 unit cells. 

    USAGE:

    sc131

    Arguments:

    None

    EXAMPLE:

    sc131

    MORE DETAILS:

    Make a lattice of 1 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 1, , orange, supercell131, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 1, , orange, supercell131, 1

    PYTHON CODE:

def sc131():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 1, , orange, supercell131, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 1, , orange, supercell131, 1')
cmd.extend("sc131", sc131)


def sc132():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 3 x 2 unit cells. 

    USAGE:

    sc132

    Arguments:

    None

    EXAMPLE:

    sc132

    MORE DETAILS:

    Make a lattice of 1 x 3 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 3, 2, , orange, supercell132, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 2, , orange, supercell132, 1

    PYTHON CODE:

def sc132():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 2, , orange, supercell132, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 2, , orange, supercell132, 1')
cmd.extend("sc132", sc132)


def sc133():
    ''' 
    DESCRIPTION:

    Make a lattice of 1 x 3 x 3 unit cells. 

    USAGE:

    sc133

    Arguments:

    None

    EXAMPLE:

    sc133

    MORE DETAILS:

    Make a lattice of 1 x 3 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/mg18OU/supercell.py;
    supercell 1, 3, 3, , orange, supercell133, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 3, , orange, supercell133, 1

    PYTHON CODE:

def sc133():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 3, , orange, supercell133, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 1, 3, 3, , orange, supercell133, 1')
cmd.extend("sc133", sc133)


def sc211():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 1 x 1 unit cells.

    USAGE:

    sc211

    Arguments:

    None

    EXAMPLE:

    sc211

    MORE DETAILS:

    Make a lattice of 2 x 1 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 1, , orange, supercell211, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 1, , orange, supercell211, 1

    PYTHON CODE:

def sc211():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 1, 1, , orange, supercell211, 1')

    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 1, 1, , orange, supercell211, 1')

cmd.extend("sc211", sc211)


def sc212():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 1 x 2 unit cells.

    USAGE:

    sc212

    Arguments:

    None

    EXAMPLE:

    sc212

    MORE DETAILS:

    Make a lattice of 2 x 1 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 2, , orange, supercell212, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 2, , orange, supercell121, 1

    PYTHON CODE:

def sc212():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 1, 2, , orange, supercell212, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 1, 2, , orange, supercell212, 1')
cmd.extend("sc212", sc212)


def sc213():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 1 x 3 unit cells. 

    USAGE:

    sc213

    Arguments:

    None

    EXAMPLE:

    sc213

    MORE DETAILS:

    Make a lattice of 2 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 3, , orange, supercell213, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 3, , orange, supercell213, 1

    PYTHON CODE:

def sc213():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    Savecmd.do('supercell 2, 1, 3, , orange, supercell213, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    Savecmd.do('supercell 2, 1, 3, , orange, supercell213, 1')
cmd.extend("sc213", sc213)


def sc221():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 2 x 1 unit cells. 

    USAGE:

    sc221

    Arguments:

    None

    EXAMPLE:

    sc221

    MORE DETAILS:

    Make a lattice of 2 x 2 x 1 unit cells. 
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 2, 1, , orange, supercell221, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;supercell 2, 2, 1, , orange, supercell221, 1

    PYTHON CODE:

def sc221():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 2, 1, , orange, supercell221, 1')

    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 2, 1, , orange, supercell221, 1')

cmd.extend("sc221", sc221)


def sc222():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 2 x 2 unit cells

    USAGE:

    sc222

    Arguments:

    None

    EXAMPLE:

    sc222

    MORE DETAILS:

    Make a lattice of 2 x 2 x 2 unit cells. 
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 2, 2, , orange, supercell222, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 2, 2, , orange, supercell222, 1

    PYTHON CODE:

def sc222():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 2, 2, , orange, supercell222, 1')

    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 2, 2, , orange, supercell222, 1')

cmd.extend("sc222", sc222)


def sc231():
    ''' 
    DESCRIPTION:

    Make a lattice of 2 x 3 x 1 unit cells. 

    USAGE:

    sc231

    Arguments:

    None

    EXAMPLE:

    sc231

    MORE DETAILS:

    Make a lattice of 2 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScriptssupercell.py;\n supercell 2, 3, 1, , orange, supercell231, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 3, 1, , orange, supercell231, 1

    PYTHON CODE:

def sc231():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 3, 1, , orange, supercell231, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 2, 3, 1, , orange, supercell231, 1')
cmd.extend("sc231", sc231)

def sc311():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 1 x 1 unit cells. 

    USAGE:

    sc311

    Arguments:

    None

    EXAMPLE:

    sc311

    MORE DETAILS:

    Make a lattice of 3 x 1 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 1, , orange, supercell311, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 1, , orange, supercell311, 1

    PYTHON CODE:

def sc311():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 1, , orange, supercell311, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 1, , orange, supercell311, 1')
cmd.extend("sc311", sc311)


def sc312():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 1 x 2 unit cells. 

    USAGE:

    sc312

    Arguments:

    None

    EXAMPLE:

    sc312

    MORE DETAILS:

    Make a lattice of 3 x 1 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 1, 2, , orange, supercell312, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 2, , orange, supercell312, 1

    PYTHON CODE:

def sc312():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 2, , orange, supercell312, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 2, , orange, supercell312, 1')
cmd.extend("sc312", sc312)

def sc313():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 1 x 3 unit cells.

    USAGE:

    sc313

    Arguments:

    None

    EXAMPLE:

    sc313

    MORE DETAILS:

    Make a lattice of 3 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py
    supercell 3, 1, 3, , orange, supercell313, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 3, , orange, supercell313, 1

    PYTHON CODE:

def sc313():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 3, , orange, supercell313, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 1, 3, , orange, supercell313, 1')
cmd.extend("sc313", sc313)


def sc321():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 2 x 1 unit cells. 

    USAGE:

    sc321

    Arguments:

    None

    EXAMPLE:

    sc321

    MORE DETAILS:

    Make a lattice of 3 x 2 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 2, 1, , orange, supercell321, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 2, 1, , orange, supercell321, 1

    PYTHON CODE:

def sc321():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 2, 1, , orange, supercell321, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 2, 1, , orange, supercell321, 1')
cmd.extend("sc321", sc321)

def sc331():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 3 x 1 unit cells. 

    USAGE:

    sc331

    Arguments:

    None

    EXAMPLE:

    sc331

    MORE DETAILS:

    Make a lattice of 3 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 3, 1, , orange, supercell331, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 3, 1, , orange, supercell331, 1

    PYTHON CODE:

def sc331():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 3, 1, , orange, supercell331, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 3, 1, , orange, supercell331, 1')
cmd.extend("sc331", sc331)


def sc333():
    ''' 
    DESCRIPTION:

    Make a lattice of 3 x 3 x 3 unit cells. 

    USAGE:

    sc333

    Arguments:

    None

    EXAMPLE:

    sc333

    MORE DETAILS:

    Make a lattice of 3 x 3 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 3, 3, , orange, supercell333, 1

    HORIZONTAL PML SCRIPT:

    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 3, 3, , orange, supercell333, 1

    PYTHON CODE:

def sc333():
    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 3, 3, , orange, supercell333, 1')
    '''

    cmd.do('run $HOME/Scripts/PyMOLScripts/supercell.py')
    cmd.do('supercell 3, 3, 3, , orange, supercell333, 1')
cmd.extend("sc333", sc333)


def sccp4(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a ccp4 electron density map with a time stamp.

    USAGE:

    sccp4

    Arguments:

    electronDensityMapFileNameStem

    EXAMPLE:

    sccp4  testFile

    MORE DETAILS:

    Save a ccp4 electron density map with a time stamp included in the 
    filename to avoid overwriting an existing file.Read as a commandline 
    argument, a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".ccp4")

    PYTHON CODE:

def sccp4(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4") 
cmd.extend('sccp4',sccp4)

def scif(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a cif file (Crystallographic Information File) with a time stamp.

    USAGE:

    scif

    Arguments:

    currentModelFileNameStem

    EXAMPLE:

    scif  testFile

    MORE DETAILS:

    Save a cif file (Crystallographic Information File) with a time stamp 
    included in the filename to avoid overwriting an existing file.Read 
    as a commandline argument, a string as the filename stem or use 
    the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".cif")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".cif")

    PYTHON CODE:

def scif(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".cif") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".cif") 
cmd.extend('scif',scif)

def sdae(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a dae file (Collada File) with a time stamp.

    USAGE:

    sdae

    Arguments:

    ColladaFileNameStem

    EXAMPLE:

    sdae  testFile

    MORE DETAILS:

    Save a dae file (Collada File) with a time stamp included in the filename 
    to avoid overwriting an existing file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".ccp4")

    PYTHON CODE:

def sdae(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".ccp4") 
cmd.extend('sdae',sdae)

def sdat(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a dat file (output data file) with a time stamp.

    USAGE:

    sdat

    Arguments:

    outputDataFileStemName

    EXAMPLE:

    sdat  testFile

    MORE DETAILS:

    Save a dat file (output data file) with a time stamp included in the filename 
    to avoid overwriting an existing dat file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".dat")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".dat")

    PYTHON CODE:

def sdat(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".dat") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".dat") 
cmd.extend('sdat',sdat)

def sfasta(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a fasta file (output data file) with a time stamp.

    USAGE:

    sfasta

    Arguments:

    fastaFileStemName

    EXAMPLE:

    sfasta  testFile

    MORE DETAILS:

    Save a fasta (sequence) file with a time stamp included in the filename 
    to avoid overwriting an existing dat file.\n Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".fasta")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".fasta")

    PYTHON CODE:

def sfasta(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".fasta") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".fasta") 
cmd.extend('sfasta',sfasta)

def sidtf(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a idtf file (Intermediate Data Text Format) with a time stamp.

    USAGE:

    sidtf

    Arguments:

    idtfFileStemName

    EXAMPLE:

    sidtf  testFile

    MORE DETAILS:

    Save a idtf (Intermediate Data Text Format) file with a time stamp 
    included in the filename to avoid overwriting an existing dat file. 
    Read as a commandline argument, a string as the filename stem, 
    or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".idtf")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".idtf")

    PYTHON CODE:

def sidtf(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".idtf") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".idtf") 
cmd.extend('sidtf',sidtf)

def smae(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a mae (Maestro) file with a time stamp.

    USAGE:

    smae

    Arguments:

    maeFileStemName

    EXAMPLE:

    smae  testFile

    MORE DETAILS:

    Save a mae (Maestro) file with a time stamp included in the filename 
    to avoid overwriting an existing dat file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".mae")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mae")

    PYTHON CODE:

def smae(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mae") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mae") 
cmd.extend('smae',smae)

def smmd(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a mmd (Macromodel) file with a time stamp.

    USAGE:

    smmd

    Arguments:

    MacromodelFileStemName

    EXAMPLE:

    smmd  testFile

    MORE DETAILS:

    Save a mmd (Macromodel) file with a time stamp included in the filename to avoid
    overwriting an existing Macromodel file. Read as a commandline argument, a string 
    as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".mmd")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mmd")

    PYTHON CODE:

def smmd(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mmd") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mmd") 
cmd.extend('smmd',smmd)

def smmod(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a mmod (Macromodel) file with a time stamp.

    USAGE:

    smmod

    Arguments:

    MacromodelFileStemName

    EXAMPLE:

    smmod  testFile

    MORE DETAILS:

    Save a mmod (Macromodel) file with a time stamp included in the filename 
    to avoid overwriting an existing Macromodel file.  Read as a commandline 
    argument, a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".mmod")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mmod")

    PYTHON CODE:

def smmod(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mmod") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mmod") 
cmd.extend('smmod',smmod)

def smoe(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a moe file (Molecular Operating Environment) with a time stamp.

    USAGE:

    smoe MOEFileStemName

    Arguments:

    MOEFileStemName

    EXAMPLE:

    smoe  testFile

    MORE DETAILS:

    Save moe file (Molecular Operating Environment) with a time stamp included 
    in the filename to avoid overwriting an existing moe file. Read as a commandline
    argument, a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".moe")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".moe")

    PYTHON CODE:

def smoe(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".moe") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".moe") 
cmd.extend('smoe',smoe)

def smol(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a mol file with a time stamp.

    USAGE:

    smol MOEFileStemName

    Arguments:

    MOEFileStemName

    EXAMPLE:

    smol  testFile

    MORE DETAILS:

    Save mol file with a time stamp included in the filename to avoid overwriting 
    an existing pmo file. Read as a commandline argument, a string as the 
    filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
   s = str(DT)
   cmd.save(stemName+s+".mol")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mol")

    PYTHON CODE:

def smol(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mol") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mol") 
cmd.extend('smol',smol)

def smol2(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a mol2 file with a time stamp.

    USAGE:

    smol2 mol2FileStemName

    Arguments:

    mol2FileStemName

    EXAMPLE:

    smol2  testFile

    MORE DETAILS:

    Save mol2 (Sybyl file format) file with a time stamp included in the filename 
    to avoid overwriting an existing mol2 file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".mol2")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mol2")

    PYTHON CODE:

def smol2(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mol2") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mol2") 
cmd.extend('smol2',smol2)

def smtl(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save mtl (Wavefront Material) file format with a time stamp.

    USAGE:

    smtl mtlFileStemName

    Arguments:

    mtlFileStemName

    EXAMPLE:

    smtl  testFile

    MORE DETAILS:

    Save mtl (Wavefront Material) file format with a time stamp included 
    in the filename to avoid overwriting an existing mtl file. Read as a 
    commandline argument, a string as the filename stem, or use the default 
    filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".mtl")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".mtl")

    PYTHON CODE:

def smtl(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mtl") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".mtl") 
cmd.extend('smtl',smtl)

def sobj(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save obj file (Wavefront mesh file format) with a time stamp.

    USAGE:

    sobj sobjFileStemName

    Arguments:

    sobjFileStemName

    EXAMPLE:

    sobj  testFile

    MORE DETAILS:

    Save obj (Wavefront mesh file format) file format with a time stamp 
    included in the filename to avoid overwriting an existing obj file. Read 
    as a commandline argument, a string as the filename stem, or use the 
    default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".obj")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".obj")

    PYTHON CODE:

def sobj(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".obj") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".obj") 
cmd.extend('sobj',sobj)

def sout(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save output data file with a time stamp.

    USAGE:

    sout outFileStemName

    Arguments:

    outFileStemName

    EXAMPLE:

    sout  testFile

    MORE DETAILS:

    Save output data file format with a time stamp included in the filename 
    to avoid overwriting an existing out file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".out")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".out")

    PYTHON CODE:

def sout(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".out") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".out") 
cmd.extend('sout',sout)

def spdb(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save pdb data file with a time stamp.

    USAGE:

    spdb pdbFileStemName

    Arguments:

    pdbFileStemName

    EXAMPLE:

    spdb testFile

    MORE DETAILS:

    Save pdb (Protein Data Bank) file format with a time stamp included in the 
    filename to avoid overwriting an existing out file. Read as a commandline 
    argument, a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pdb")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pdb")

    PYTHON CODE:

def spdb(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pdb") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pdb") 
cmd.extend('spdb',spdb)

def spkl(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a pkl file (Python pickle file) with a time stamp.

    USAGE:

    spkl pklFileStemName

    Arguments:

    pklFileStemName

    EXAMPLE:

    spkl testFile

    MORE DETAILS:

    Save a pkl file (Python pickle file) with a time stamp included in the filename 
    to avoid overwriting an existing out file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)\
    cmd.save(stemName+s+".pkl")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pkl")

    PYTHON CODE:

def spkl(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pkl") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pkl") 
cmd.extend('spkl',spkl)

def spkla(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a pkla file (Python pickle file) with a time stamp.

    USAGE:

    spkla pklaFileStemName

    Arguments:

    pklaFileStemName

    EXAMPLE:

    spkla testFile

    MORE DETAILS:

    Save a pkla file (Python pickle file) with a time stamp included in 
    the filename to avoid overwriting an existing out file. Read as a 
    commandline argument, a string as the filename stem, or use the 
    default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pkla")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pkla")

    PYTHON CODE:

def spkla(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pkla") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pkla") 
cmd.extend('spkla',spkla)

def spmo(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a pmo file (XYZ, binary format file) with a time stamp.

    USAGE:

    spmo MacromodelFileStemName

    Arguments:

    MacromodelFileStemName

    EXAMPLE:

    spmo testFile

    MORE DETAILS:

    Save pmo file (XYZ, binary format file) with a time stamp included in the filename 
    to avoid overwriting an existing pmo file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pmo")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pmo")

    PYTHON CODE:

def spmo(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pmo") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pmo") 
cmd.extend('spmo',spmo)

def spng(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save a png file (Python pickle file) with a time stamp.

    USAGE:

    spng pngFileStemName

    Arguments:

    pngFileStemName

    EXAMPLE:

    spng testFile

    MORE DETAILS:

    Save a png file (Python pickle file) with a time stamp included in the filename to avoid     
    overwriting an existing out file.\n Read as a commandline argument, a string as the 
    ilename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".png")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".png")

    PYTHON CODE:

def spng(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".png") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".png") 
cmd.extend('spng',spng)

def spov(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save pov (POV-ray tracing file format) file with a time stamp.

    USAGE:

    spov povFileStemName

    Arguments:

    povFileStemName

    EXAMPLE:

    spov testFile

    MORE DETAILS:

    Save pov (POV-ray tracing file format) file with a time stamp included in the filename 
    to avoid overwriting an existing out file. Read as a commandline argument, a string 
    as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pov")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pov")

    PYTHON CODE:

def spov(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pov") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pov") 
cmd.extend('spov',spov)

def spqr(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save pqr with file with timestamp.

    USAGE:

    spqr pqrFileStemName

    Arguments:

    pqrFileStemName

    EXAMPLE:

    spqr testFile

    MORE DETAILS:

    Save pqr with file with timestamp.

A pqr (PDB file with the temperature and occupancy columns replaced
    by columns containing the per-atom charge (Q) and radius (R)) file with 
    a time stamp included in the filename to avoid overwriting an existing 
    out file. Read as a commandline argument, a string as the filename stem, 
    or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pqr")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pqr")

    PYTHON CODE:

def spqr(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pqr") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pqr") 
cmd.extend('spqr',spqr)

def spse(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save session file with a time stamp.

    USAGE:

    spse pseFileStemName

    Arguments:

    pseFileStemName

    EXAMPLE:

    spse testFile

    MORE DETAILS:

    Save session file with a time stamp included in the filename to avoid 
    overwriting an existing pse file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".pse")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".pse")

    PYTHON CODE:

def spse(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pse") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pse") 
cmd.extend('spse',spse)

def ssdf(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save session file with a time stamp.

    USAGE:

    ssdf sdfFileStemName

    Arguments:

    sdfFileStemName

    EXAMPLE:

    ssdf testFile

    MORE DETAILS:

    Save session file with a time stamp included in the filename to avoid 
    overwriting an existing sdf file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".sdf")

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".sdf")

    PYTHON CODE:

def ssdf(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".sdf") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".sdf") 
cmd.extend('ssdf',ssdf)

def st3(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open Sublime Text 3 from within PyMOL.

    USAGE:

    st3

    Arguments:

    None

    EXAMPLE:

    st3 script.pml

    MORE DETAILS:

    Open Open Sublime Text 3 from within PyMOL.
    Adjust file path to executable.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (sublimetext3Path  + "-w " + fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (sublimetext3Path  + "-w " + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def st3(fileName="test.pml"):
    arg = (sublimetext3Path  + "-w " + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (sublimetext3Path  + "-w " + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('st3',st3)

def swrl(stemName="saved"):
    ''' 
    DESCRIPTION:

    Save wrl (VRML 2 file format) file with a time stamp.

    USAGE:

    swrl wrlFileStemName

    Arguments:

    wrlFileStemName

    EXAMPLE:

    swrl testFile

    MORE DETAILS:

    Save wrl (VRML 2 file format) file with a time stamp included in the filename 
    to avoid overwriting an existing sdf file. Read as a commandline argument, 
    a string as the filename stem, or use the default filename stem "saved".

    VERTICAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S")
    s = str(DT)
    cmd.save(stemName+s+".wrl")'

    HORIZONTAL PML SCRIPT:

    DT =datetime.datetime.now().strftime("y%Ym%md%dh%Hm%Ms%S");s = str(DT);cmd.save(stemName+s+".wrl")

    PYTHON CODE:

def swrl(stemName="saved"):
    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".wrl") 
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".wrl") 
cmd.extend('swrlf',swrl)

def term():
    ''' 
    DESCRIPTION:

    Open a Terminal window on MacOS.

    USAGE:

    term

    Arguments:

    None

    EXAMPLE:

    term

    MORE DETAILS:

    Open a Terminal window on MacOS.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    subprocess.call(terminalCommand)
    return


    HORIZONTAL PML SCRIPT:

    subprocess.call(terminalCommand);return


    PYTHON CODE:

def term():
    subprocess.call(terminalCommand)
    return

    '''

    subprocess.call(terminalCommand)
    return

cmd.extend('term',term)

def vdw(selection='all'):
    ''' 
    DESCRIPTION:

    Transparent vdw surface

    USAGE:

    vdw selection

    Arguments:

    selection


    EXAMPLE:

    vdw 3nd3

    MORE DETAILS:

    # vdw creates a copy of an object with full-sized, transparent spheres
# Bondi VDW values added below to override default Pymol settings
# From https://gist.githubusercontent.com/bobbypaton/1cdc4784f3fc8374467bae5eb410edef/raw/9995d51d6a64b8bcf01590c944cc38059b2f8d7f/pymol_style.py


    VERTICAL PML SCRIPT:

    NA

    HORIZONTAL PML SCRIPT:

    NA

    PYTHON CODE:

def vdw(selection='all'):
    # Bondi VDW values 
    cmd.alter("elem Ac", "vdw=2.00")
    cmd.alter("elem Al", "vdw=2.00")
    cmd.alter("elem Am", "vdw=2.00")
    cmd.alter("elem Sb", "vdw=2.00")
    cmd.alter("elem Ar", "vdw=1.88")
    cmd.alter("elem As", "vdw=1.85")
    cmd.alter("elem At", "vdw=2.00")
    cmd.alter("elem Ba", "vdw=2.00")
    cmd.alter("elem Bk", "vdw=2.00")
    cmd.alter("elem Be", "vdw=2.00")
    cmd.alter("elem Bi", "vdw=2.00")
    cmd.alter("elem Bh", "vdw=2.00")
    cmd.alter("elem B ", "vdw=2.00")
    cmd.alter("elem Br", "vdw=1.85")
    cmd.alter("elem Cd", "vdw=1.58")
    cmd.alter("elem Cs", "vdw=2.00")
    cmd.alter("elem Ca", "vdw=2.00")
    cmd.alter("elem Cf", "vdw=2.00")
    cmd.alter("elem C ", "vdw=1.70")
    cmd.alter("elem Ce", "vdw=2.00")
    cmd.alter("elem Cl", "vdw=1.75")
    cmd.alter("elem Cr", "vdw=2.00")
    cmd.alter("elem Co", "vdw=2.00")
    cmd.alter("elem Cu", "vdw=1.40")
    cmd.alter("elem Cm", "vdw=2.00")
    cmd.alter("elem Ds", "vdw=2.00")
    cmd.alter("elem Db", "vdw=2.00")
    cmd.alter("elem Dy", "vdw=2.00")
    cmd.alter("elem Es", "vdw=2.00")
    cmd.alter("elem Er", "vdw=2.00")
    cmd.alter("elem Eu", "vdw=2.00")
    cmd.alter("elem Fm", "vdw=2.00")
    cmd.alter("elem F ", "vdw=1.47")
    cmd.alter("elem Fr", "vdw=2.00")
    cmd.alter("elem Gd", "vdw=2.00")
    cmd.alter("elem Ga", "vdw=1.87")
    cmd.alter("elem Ge", "vdw=2.00")
    cmd.alter("elem Au", "vdw=1.66")
    cmd.alter("elem Hf", "vdw=2.00")
    cmd.alter("elem Hs", "vdw=2.00")
    cmd.alter("elem He", "vdw=1.40")
    cmd.alter("elem Ho", "vdw=2.00")
    cmd.alter("elem In", "vdw=1.93")
    cmd.alter("elem I ", "vdw=1.98")
    cmd.alter("elem Ir", "vdw=2.00")
    cmd.alter("elem Fe", "vdw=2.00")
    cmd.alter("elem Kr", "vdw=2.02")
    cmd.alter("elem La", "vdw=2.00")
    cmd.alter("elem Lr", "vdw=2.00")
    cmd.alter("elem Pb", "vdw=2.02")
    cmd.alter("elem Li", "vdw=1.82")
    cmd.alter("elem Lu", "vdw=2.00")
    cmd.alter("elem Mg", "vdw=1.73")
    cmd.alter("elem Mn", "vdw=2.00")
    cmd.alter("elem Mt", "vdw=2.00")
    cmd.alter("elem Md", "vdw=2.00")
    cmd.alter("elem Hg", "vdw=1.55")
    cmd.alter("elem Mo", "vdw=2.00")
    cmd.alter("elem Nd", "vdw=2.00")
    cmd.alter("elem Ne", "vdw=1.54")
    cmd.alter("elem Np", "vdw=2.00")
    cmd.alter("elem Ni", "vdw=1.63")
    cmd.alter("elem Nb", "vdw=2.00")
    cmd.alter("elem N ", "vdw=1.55")
    cmd.alter("elem No", "vdw=2.00")
    cmd.alter("elem Os", "vdw=2.00")
    cmd.alter("elem O ", "vdw=1.52")
    cmd.alter("elem Pd", "vdw=1.63")
    cmd.alter("elem P ", "vdw=1.80")
    cmd.alter("elem Pt", "vdw=1.72")
    cmd.alter("elem Pu", "vdw=2.00")
    cmd.alter("elem Po", "vdw=2.00")
    cmd.alter("elem K ", "vdw=2.75")
    cmd.alter("elem Pr", "vdw=2.00")
    cmd.alter("elem Pm", "vdw=2.00")
    cmd.alter("elem Pa", "vdw=2.00")
    cmd.alter("elem Ra", "vdw=2.00")
    cmd.alter("elem Rn", "vdw=2.00")
    cmd.alter("elem Re", "vdw=2.00")
    cmd.alter("elem Rh", "vdw=2.00")
    cmd.alter("elem Rb", "vdw=2.00")
    cmd.alter("elem Ru", "vdw=2.00")
    cmd.alter("elem Rf", "vdw=2.00")
    cmd.alter("elem Sm", "vdw=2.00")
    cmd.alter("elem Sc", "vdw=2.00")
    cmd.alter("elem Sg", "vdw=2.00")
    cmd.alter("elem Se", "vdw=1.90")
    cmd.alter("elem Si", "vdw=2.10")
    cmd.alter("elem Ag", "vdw=1.72")
    cmd.alter("elem Na", "vdw=2.27")
    cmd.alter("elem Sr", "vdw=2.00")
    cmd.alter("elem S ", "vdw=1.80")
    cmd.alter("elem Ta", "vdw=2.00")
    cmd.alter("elem Tc", "vdw=2.00")
    cmd.alter("elem Te", "vdw=2.06")
    cmd.alter("elem Tb", "vdw=2.00")
    cmd.alter("elem Tl", "vdw=1.96")
    cmd.alter("elem Th", "vdw=2.00")
    cmd.alter("elem Tm", "vdw=2.00")
    cmd.alter("elem Sn", "vdw=2.17")
    cmd.alter("elem Ti", "vdw=2.00")
    cmd.alter("elem W ", "vdw=2.00")
    cmd.alter("elem U ", "vdw=1.86")
    cmd.alter("elem V ", "vdw=2.00")
    cmd.alter("elem Xe", "vdw=2.16")
    cmd.alter("elem Yb", "vdw=2.00")
    cmd.alter("elem Y ", "vdw=2.00")
    cmd.alter("elem Zn", "vdw=1.39")
    cmd.alter("elem Zr", "vdw=2.00")
    cmd.rebuild()

    # workspace settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 0)
    cmd.set("transparency", 0.5)
    cmd.set("dash_gap", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_texture", 2)
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", .1)
    cmd.space("cmyk")
   
    # ball and stick settings
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85","elem C and "+arg1)
    cmd.color("gray98","elem H and "+arg1)
    cmd.color("slate","elem N and "+arg1)
    cmd.set("stick_radius",0.07, arg1)
    cmd.set("sphere_scale",0.18, arg1)
    cmd.set("sphere_scale",0.13, arg1+" and elem H")
    cmd.set("dash_gap",0.01, arg1)
    cmd.set("dash_radius",0.07, arg1)
    cmd.set("stick_color","black", arg1)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.zoom(arg1)
    cmd.hide("labels")

    '''

    # Bondi VDW values 
    cmd.alter("elem Ac", "vdw=2.00")
    cmd.alter("elem Al", "vdw=2.00")
    cmd.alter("elem Am", "vdw=2.00")
    cmd.alter("elem Sb", "vdw=2.00")
    cmd.alter("elem Ar", "vdw=1.88")
    cmd.alter("elem As", "vdw=1.85")
    cmd.alter("elem At", "vdw=2.00")
    cmd.alter("elem Ba", "vdw=2.00")
    cmd.alter("elem Bk", "vdw=2.00")
    cmd.alter("elem Be", "vdw=2.00")
    cmd.alter("elem Bi", "vdw=2.00")
    cmd.alter("elem Bh", "vdw=2.00")
    cmd.alter("elem B ", "vdw=2.00")
    cmd.alter("elem Br", "vdw=1.85")
    cmd.alter("elem Cd", "vdw=1.58")
    cmd.alter("elem Cs", "vdw=2.00")
    cmd.alter("elem Ca", "vdw=2.00")
    cmd.alter("elem Cf", "vdw=2.00")
    cmd.alter("elem C ", "vdw=1.70")
    cmd.alter("elem Ce", "vdw=2.00")
    cmd.alter("elem Cl", "vdw=1.75")
    cmd.alter("elem Cr", "vdw=2.00")
    cmd.alter("elem Co", "vdw=2.00")
    cmd.alter("elem Cu", "vdw=1.40")
    cmd.alter("elem Cm", "vdw=2.00")
    cmd.alter("elem Ds", "vdw=2.00")
    cmd.alter("elem Db", "vdw=2.00")
    cmd.alter("elem Dy", "vdw=2.00")
    cmd.alter("elem Es", "vdw=2.00")
    cmd.alter("elem Er", "vdw=2.00")
    cmd.alter("elem Eu", "vdw=2.00")
    cmd.alter("elem Fm", "vdw=2.00")
    cmd.alter("elem F ", "vdw=1.47")
    cmd.alter("elem Fr", "vdw=2.00")
    cmd.alter("elem Gd", "vdw=2.00")
    cmd.alter("elem Ga", "vdw=1.87")
    cmd.alter("elem Ge", "vdw=2.00")
    cmd.alter("elem Au", "vdw=1.66")
    cmd.alter("elem Hf", "vdw=2.00")
    cmd.alter("elem Hs", "vdw=2.00")
    cmd.alter("elem He", "vdw=1.40")
    cmd.alter("elem Ho", "vdw=2.00")
    cmd.alter("elem In", "vdw=1.93")
    cmd.alter("elem I ", "vdw=1.98")
    cmd.alter("elem Ir", "vdw=2.00")
    cmd.alter("elem Fe", "vdw=2.00")
    cmd.alter("elem Kr", "vdw=2.02")
    cmd.alter("elem La", "vdw=2.00")
    cmd.alter("elem Lr", "vdw=2.00")
    cmd.alter("elem Pb", "vdw=2.02")
    cmd.alter("elem Li", "vdw=1.82")
    cmd.alter("elem Lu", "vdw=2.00")
    cmd.alter("elem Mg", "vdw=1.73")
    cmd.alter("elem Mn", "vdw=2.00")
    cmd.alter("elem Mt", "vdw=2.00")
    cmd.alter("elem Md", "vdw=2.00")
    cmd.alter("elem Hg", "vdw=1.55")
    cmd.alter("elem Mo", "vdw=2.00")
    cmd.alter("elem Nd", "vdw=2.00")
    cmd.alter("elem Ne", "vdw=1.54")
    cmd.alter("elem Np", "vdw=2.00")
    cmd.alter("elem Ni", "vdw=1.63")
    cmd.alter("elem Nb", "vdw=2.00")
    cmd.alter("elem N ", "vdw=1.55")
    cmd.alter("elem No", "vdw=2.00")
    cmd.alter("elem Os", "vdw=2.00")
    cmd.alter("elem O ", "vdw=1.52")
    cmd.alter("elem Pd", "vdw=1.63")
    cmd.alter("elem P ", "vdw=1.80")
    cmd.alter("elem Pt", "vdw=1.72")
    cmd.alter("elem Pu", "vdw=2.00")
    cmd.alter("elem Po", "vdw=2.00")
    cmd.alter("elem K ", "vdw=2.75")
    cmd.alter("elem Pr", "vdw=2.00")
    cmd.alter("elem Pm", "vdw=2.00")
    cmd.alter("elem Pa", "vdw=2.00")
    cmd.alter("elem Ra", "vdw=2.00")
    cmd.alter("elem Rn", "vdw=2.00")
    cmd.alter("elem Re", "vdw=2.00")
    cmd.alter("elem Rh", "vdw=2.00")
    cmd.alter("elem Rb", "vdw=2.00")
    cmd.alter("elem Ru", "vdw=2.00")
    cmd.alter("elem Rf", "vdw=2.00")
    cmd.alter("elem Sm", "vdw=2.00")
    cmd.alter("elem Sc", "vdw=2.00")
    cmd.alter("elem Sg", "vdw=2.00")
    cmd.alter("elem Se", "vdw=1.90")
    cmd.alter("elem Si", "vdw=2.10")
    cmd.alter("elem Ag", "vdw=1.72")
    cmd.alter("elem Na", "vdw=2.27")
    cmd.alter("elem Sr", "vdw=2.00")
    cmd.alter("elem S ", "vdw=1.80")
    cmd.alter("elem Ta", "vdw=2.00")
    cmd.alter("elem Tc", "vdw=2.00")
    cmd.alter("elem Te", "vdw=2.06")
    cmd.alter("elem Tb", "vdw=2.00")
    cmd.alter("elem Tl", "vdw=1.96")
    cmd.alter("elem Th", "vdw=2.00")
    cmd.alter("elem Tm", "vdw=2.00")
    cmd.alter("elem Sn", "vdw=2.17")
    cmd.alter("elem Ti", "vdw=2.00")
    cmd.alter("elem W ", "vdw=2.00")
    cmd.alter("elem U ", "vdw=1.86")
    cmd.alter("elem V ", "vdw=2.00")
    cmd.alter("elem Xe", "vdw=2.16")
    cmd.alter("elem Yb", "vdw=2.00")
    cmd.alter("elem Y ", "vdw=2.00")
    cmd.alter("elem Zn", "vdw=1.39")
    cmd.alter("elem Zr", "vdw=2.00")
    cmd.rebuild()

    # workspace settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 0)
    cmd.set("transparency", 0.5)
    cmd.set("dash_gap", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_texture", 2)
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", .1)
    cmd.space("cmyk")
   
    # ball and stick settings
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85","elem C and "+arg1)
    cmd.color("gray98","elem H and "+arg1)
    cmd.color("slate","elem N and "+arg1)
    cmd.set("stick_radius",0.07, arg1)
    cmd.set("sphere_scale",0.18, arg1)
    cmd.set("sphere_scale",0.13, arg1+" and elem H")
    cmd.set("dash_gap",0.01, arg1)
    cmd.set("dash_radius",0.07, arg1)
    cmd.set("stick_color","black", arg1)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.zoom(arg1)
    cmd.hide("labels")

cmd.extend('vdw', vdw)

def vim(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open vim from within PyMOL.

    USAGE:

    vim

    Arguments:

    None

    EXAMPLE:

    vim script.pml

    MORE DETAILS:

    Open vim from within PyMOL. 
    Adjust file path to vim on your computer.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    subprocess.call(vimCommand)
    return


    HORIZONTAL PML SCRIPT:

    subprocess.call(vimCommand);return


    PYTHON CODE:

def vim(fileName="test.pml"):
    subprocess.call(vimCommand)
    return

    '''

    subprocess.call(vimCommand)
    return

cmd.extend('vim',vim)

def vmd(fileName="test.pdb"):
    ''' 
    DESCRIPTION:

    Open vmd from within PyMOL.

    USAGE:

    vmd

    Arguments:

    None

    EXAMPLE:

    vmd

    MORE DETAILS:

    Open vmd from within PyMOL.
    Adjust file path for your location of vmd.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (vmdCommand+ fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (vmdCommand + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def vmd(fileName="test.pdb"):
    arg = (vmdCommand+ fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (vmdCommand+ fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('vmd',vmd)

def weather():
    ''' 
    DESCRIPTION:

    Open National Weather Service website for locale. 

    USAGE:

    weather

    Arguments:

    None

    EXAMPLE:

    weather

    MORE DETAILS:

    Open National Weather Service website for locale. 
    Adjust url for your location.


    VERTICAL PML SCRIPT:

    webbrowser.open('https://www.weather.gov/oun/')

    HORIZONTAL PML SCRIPT:

    webbrowser.open('https://www.weather.gov/oun/')

    PYTHON CODE:

def weather():
    webbrowser.open('https://www.weather.gov/oun/')
    '''

    webbrowser.open('https://www.weather.gov/oun/')
cmd.extend('weather',weather)

def webmail():
    ''' 
    DESCRIPTION:

    Open Web Mail in defualt browser.

    USAGE:

    webmail

    Arguments:

    None

    EXAMPLE:

    webmail

    MORE DETAILS:

    Open Web Mail in defualt browser. Adjust url for your institution.


    VERTICAL PML SCRIPT:

    webbrowser.open('https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f')

    HORIZONTAL PML SCRIPT:

    webbrowser.open('https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f')

    PYTHON CODE:

def webmail():
    webbrowser.open('https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f')
    '''

    webbrowser.open('https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f')
cmd.extend('webmail',webmail)

def word(fileName="Script.docx"):
    ''' 
    DESCRIPTION:

    Open word from within PyMOL.

    USAGE:

    word

    Arguments:

    None

    EXAMPLE:

    word paper.docx

    MORE DETAILS:

    Open MS Word from within PyMOL. 
    Adjust file path to MS Word on your computer.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    subprocess.call(wordCommand, '--args', fileName);
    return


    HORIZONTAL PML SCRIPT:

    subprocess.call(wordCommand, '--args', fileName);return


    PYTHON CODE:

def word(fileName="Script.docx"):
    subprocess.call(wordCommand, '--args', fileName)
    return

    '''

    subprocess.call(wordCommand, '--args', fileName)
    return

cmd.extend('word',word)

def yasara(fileName="test.pml"):
    ''' 
    DESCRIPTION:

    Open the molecular graphics prograom YASASRA from within PyMOL.

    USAGE:

    yasara

    Arguments:

    None

    EXAMPLE:

    yasara

    MORE DETAILS:

    Open the molecular graphics prograom YASASRA from within PyMOL.

   >>>  Edit file path in python code below

    VERTICAL PML SCRIPT:

    arg = (yasaraPath + fileName);
    subprocess.call(arg,shell=True);
    return


    HORIZONTAL PML SCRIPT:

    arg = (yasaraPath + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:

def yasara(fileName="test.pml"):
    arg = (yasaraPath + fileName)
    subprocess.call(arg,shell=True)
    return

    '''

    arg = (yasaraPath + fileName)
    subprocess.call(arg,shell=True)
    return

cmd.extend('yasara',yasara)


 
""" Print the shortcuts on startup of PyMOL"""
print(SC.__doc__)
