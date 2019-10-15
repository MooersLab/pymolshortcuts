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
import datetime, time, webbrowser, random, glob
import os, os.path

from pymol import cmd, stored, math, cgo, xray
import numpy

__author__ = "Blaine Mooers"
__copyright__ = "2019 Board of Regents for the University of Oklahoma"
__license__ = "MIT Licencse"
__version__ = "1.0.0" 
# Versioning follows follows MAJOR.MINOR[.PATCH] where major releases are not backward compatable. 
__credits__ = ["Miriam Shakir", "Safra Shakir"] 
# Credits are for people who tested the code, reported bug fixes, made suggestions, etc. 
__date__ = "2 October 2019"
__maintainer__ = "Blaine Mooers"
__email__ = "blaine-mooers@ouhsc.edu"
__status__ = "Production" 

cmd.set('ray_opaque_background','on')


##################################  Edit PATHS to local directories of files ###############################################
# $HOME/ is the user's main directory: /Users/username or /home/username. r'$HOME' is equivalent to ~. 
localPDBfilePath = r'$HOME/pdbFiles/'
localEMAPfilePath = r'$HOME/emapFiles/'
localHKLfilePath = r'$HOME/hklFiles/'


AppPaths='''You may have to edit the file paths to your applications around line 477 in pymolshortcut.py.'''
print(AppPaths)

##################################  Edit PATHS to Applications ###############################################
#
# Find your operating system from the three sets of paths or commands below.
# The default OS is Mac.
# If your computer is not a Mac, comment out the Mac set and uncomment the set for your operating system.
# Comment out the paths to the applications that you lack. 
#
# The paths are for my Mac OS 10.14 with some applications stored in /usr/local 
# and some in /opt/local/bin thanks to macports.
# Some paths can be found by entering "which <program>" in a terminal window.
#
# You can replace the file path with a command to open the application.
# For example, you can use the command "open -a bbedit" on the Mac. 
# This approach is the preferred one: it avoids the spinning ball.
# The analogous command is "start" on Windows.
#
##################################  Mac Open Commands and PATHS to Applications ###############################################

# Data analysis
DBBrowserSQLiteOpen = ['open','-a','DBBrowserForSQLite']
DBBrowserSQLitePath = '/Applications/DBBrowserForSQLite.app/Contents/MacOS/DB\ Browser\ for\ SQLite'
excelOpen = ['open','-a','Microsoft Excel']
excelPath = r'/Applications/Microsoft\ Excel.app/Contents/MacOS/Microsoft\ Excel'
JabRefOpen = ['open','-a','JabRef.app']
JabRefPath = r'/Applications/JabRef.app/Contents/MacOS/JavaApplicationStub'
jaspOpen = ['open','-a','JASP.app']
jaspPath = r'Applications/JASP.app/Contents/MacOS/JASP'
jmpOpen = ['open','-a','JMP Pro 14.app']
jmpPath = r'/Applications/JMP\ Pro\ 14.app/Contents/MacOS/JMP'
juliaOpen = ['open','-a','Julia-1.2.app']
juliaPath = r'/Applications/Julia-1.2.app/Contents/MacOS/applet'
juliaproOpen = ['open','-a','JuliaPro-1.2.0-1.app']
juliaproPath = r'/Applications/JuliaPro-1.2.0-1.app/Contents/MacOS/applet'
jupyterPath = r'/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter-notebook'
octaveCommand = '/Applications/MacPorts/Octave.app/Contents/MacOS/Octave'
octaveOpen = ['open', '-n','-a','/Applications/MacPorts/Octave.app/Contents/MacOS/Octave']
octavePath = '/Applications/MacPorts/Octave.app/Contents/MacOS/Octave'
Rpath = r'/usr/local/bin/R'
ROpen = ['open','-a','R'] 
RStudioOpen = ['open','-a','RStudio'] 
RStudioPath = r'/Applications/RStudio.app/Contents/MacOS/Rstudio'


# Image manipulation programs
gimpOpen=['/opt/local/bin/gimp','-n']
gimpPath = r'/opt/local/bin/gimp'

inkscapeCommand = ['/opt/local/bin/inkscape','--file=','fileName']
inkscapeOpen = ['/opt/local/bin/inkscape']
inkscapePath = '/opt/local/bin/inkscape'

pptOpen = ['open','-a','Microsoft PowerPoint']
pptPath =  r'/Applications/Microsoft\ Excel.app/Contents/MacOS/Microsoft\ PowerPoint'


# Molecular graphics programs
ccp4mgCommand = ['/Applications/ccp4-7.0/QtMG.app/Contents/MacOS/QtMG','shell=True']
ccp4mgOpen = ['open','-n','-a','QtMG.app']
ccp4mgPath = r'/Applications/ccp4-7.0/QtMG.app/Contents/MacOS/QtMG'
chimeraOpen = ['open','-a','Chimera.app']
chimeraPath = r'/Applications/Chimera.app/Contents/MacOS/chimera'
cootOpen =['open','-a','coot.app']
cootPath = '/usr/local/bin/coot'
jmolCommand = 'java -jar /Applications/jars/jmol-14.29.54/Jmol.jar'
jmolOpen = ['java','-jar','/Applications/jars/jmol-14.29.54/Jmol.jar']
jmolPath = '/Applications/jars/jmol-14.29.54/Jmol.jar'
yasaraOpen = ['open','-a','YASARA.app']
yasaraPath = r'/Applications/YASARA.app/Contents/MacOS/yasara.app'
vmdOpen = ['open','-a','VMD194.app']
vmdPath = r'/Applications/VMD194.app/Contents/MacOS/startup.command'

# Text editors
atomOpen = ['open','-a','atom']
atomPath = r'/usr/local/bin/atom'
atomStart = ['start','atom']
bbeditOpen = ['open','-a','bbedit']
bbeditPath = r'/usr/local/bin/bbedit'
codeOpen = ['/usr/local/bin/code']
codePath = r'/usr/local/bin/code'
emacsCommand = ['open','-a','Emacs.app', '"$@"']
emacsOpen = ['open','-a','Emacs.app']
emacsPath = r'/opt/local/bin/emacs'
geditOpen = ['open','-a','gedit2.30.2.app']
geditPath = r'/Applications/gedit2.30.2.app/Contents/MacOS/gedit'
jeditOpen = ['open','-a','jEdit.app']
jeditPath = r'/Applications/jEdit.app/Contents/MacOS/jedit'
neovimPath = r'/Users/blaine/software/nvim-osx64/bin/nvim'
neovimOpen = ['open','-n','-a','/Users/blaine/software/nvim-osx64/bin/nvim']
nppOpen = ['open','-a','Notepad++.app']
nppPath = r'/Applications/Notepad++.app/Contents/MacOS/startwine'
oniOpen = ['open','-a','Oni.app']
oniPath = r'/Applications/Oni.app/Contents/MacOS/Oni'
pdbeditorOpen = ['java','-jar','/Applications/jars/PDB_Editor_FIX090203.jar']
pdbeditorPath = r'/Applications/jars/PDB_Editor_FIX090203.jar'
sublimeText3Open = ['/usr/local/bin/subl','test.pml']
sublimeText3Path = r'/usr/local/bin/subl'
textmatePath  = r'/usr/local/bin/mate'
textMateOpen = ['open','-a','TextMate']
vimPath = r'/opt/local/bin/vim'
vimOpen = ['open','-a','MacVim.app']

#Terminals
itermOpen = ['open','-n','-a','iTerm.app']
terminalOpen = ['open','-n','-a','Terminal.app']
xquartzOpen = ['open','-n','-a','XQuartz.app']
x11Open = ['open','-n','-a','X11.app']


# Web sites 
gcalURL = 'https://calendar.google.com/calendar/r'
gmailURL = 'https://mail.google.com/mail/u/0/#inbox'
webmailURL = 'https://webmail.ouhsc.edu/owa/auth/logon.aspx?replaceCurrent=1&url=http%3a%2f%2fwebmail.ouhsc.edu%2fowa%2f'
weatherServiceRadarURL = 'https://radar.weather.gov/radar.php?rid=TLX'

# Web static sites
atsasURL = 'https://www.embl-hamburg.de/biosaxs/download.html'
scatterURL = 'http://www.bioisis.net/tutorial/9'
bioxtasrawURL = 'https://bioxtas-raw.readthedocs.io/en/latest/manual/Introduction_to_RAW_and_this_documentation.html'
acaURL = 'http://www.amercrystalassn.org'
alsURL = 'https://als.lbl.gov/'
apsURL = 'https://www.aps.anl.gov/'
biocatURL = 'https://www.bio.aps.anl.gov/pages/links.html'
sasbdbURL = 'https://www.sasbdb.org/'
chimeriaURL = 'https://www.cgl.ucsf.edu/chimera/'
chessURL = 'https://www.chess.cornell.edu'
emdbURL = 'https://www.ebi.ac.uk/pdbe/emdb/'
epURL = 'https://github.com/MooersLab/EasyPyMOL'
jmURL = 'http://wiki.jmol.org/index.php/Main_Page'
lbsfURL = 'https://research.ouhsc.edu/CoreFacilities/LaboratoryofBiomolecularStructureandFunction.aspx'
ouMCLURL = 'http://structuralbiology.ou.edu/mcl'
molgrURL = 'https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/moleculargraphicslinks.html'
molgrwikiURL = 'https://en.wikipedia.org/wiki/List_of_molecular_graphics_systems'
nslsIIURL = 'https://www.bnl.gov/ps/'
ndbURL = 'http://ndbserver.rutgers.edu'

ppcURL = 'http://www.ou.edu/cas/chemistry/research/research-support-services/protein-production-core'
psURL = 'https://www.proteinsociety.org/'
rsURL = 'https://www.rnasociety.org/'
saxsURL = 'https://www.oumedicine.com/docs/default-source/ad-biochemistry-workfiles/small-angle-scattering-links-27aug2014.html?sfvrsn=0'
saxierURL = 'https://www.saxier.org/forum/'
ssrlbl42URL = 'https://www-ssrl.slac.stanford.edu/~saxs/'
sbgridURL= 'https://www.youtube.com/user/SBGridTV/videos'
ssrlsmbURL = 'http://smb.slac.stanford.edu'
ssurfURL = 'http://www.ssurf.org'
scipyURL = 'https://www.scipy2019.scipy.org'

#### Web search sites
# Send search term to Amazon.com
abURL = 'https://www.amazon.com/s/ref=nb_sb_noss_2?url=search-alias%3Dstripbooks&field-keywords='


anacondaURL = 'https://anaconda.org/search?q='

# Send search term to arxiv
axURL = 'https://arxiv.org/search/?query='

# Send search term to Biorxiv 
bxURL = 'https://www.biorxiv.org/search/'


gbURL = 'https://www.google.com/search?tbm=bks&q='


ghURL = 'https://www.github.com/search?q='

# Send search terms to Google.com. 
goURL = 'http://www.google.com/search?btnG=1&q='

# Had to use Sweden's Google Scholar.
gsURL = 'https://scholar.google.se/scholar?hl=en&q='
gsnURL = 'https://scholar.google.se/scholar?hl=en&q='

# Send search terms to Google Videos.
gvURL = 'https://www.google.com/videohp/search?btnG=1&q='

# Send search terms to PubMed. 
ipmURL = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=' 

iucrURL = 'https://journals.iucr.org/'



stackoverflowURL = 'https://stackoverflow.com/search?q='
spnURL = 'https://www.springer.com/gp/search?query='
spURL = 'https://www.springer.com/gp/search?query='
scienceDirectURL1 = 'https://www.sciencedirect.com/search/advanced?qs='
scienceDirectURL2 = '&show=100&sortBy=relevance'

rgnURL = 'https://www.researchgate.net/search.Search.html?type=researcher&query='
rgURL = 'https://www.researchgate.net/search.Search.html?type=researcher&query='
pwURL = 'https://pymolwiki.org/index.php/'
pmlURL = 'https://sourceforge.net/p/pymol/mailman/search/?q='
pmlnURL = 'https://sourceforge.net/p/pymol/mailman/search/?q='
pmURL = 'https://www.ncbi.nlm.nih.gov/pubmed/?term='
pdbURL = 'https://www.rcsb.org/structure/'
pdbnURL = 'https://www.rcsb.org/structure/'
pymolURL = 'https://pymolwiki.org/index.php/'

researchGateURL = 'https://www.researchgate.net/search.Search.html?type=researcher&query='

scienceDirectURL1 = 'https://www.sciencedirect.com/search/advanced?qs='
scienceDirectURL2 = '&show=100&sortBy=relevance'

sourceForgeURL = "https://sourceforge.net/directory/os:mac/?q="

springerBooksURL1 = 'https://www.springer.com/gp/search?query='
springerBooksURL2 = '&submit=Submit+Query'


#wordProcessor
wordOpen = ['open','-a','Microsoft Word.app']

# Optional local mirror of the PDB
local_mirror_divided = '/mnt/bio/db/pdb.divided'

# ##################################  Linux run commands and PATHS to Applications ###############################################
#
# # Data analysis
# DBBrowserSQLiteStart = ['DBBrowserForSQLite']
# DBBrowserSQLitePath = '/Applications/DBBrowserForSQLite/Contents/MacOS/DB\ Browser\ for\ SQLite'
# excelOpen = ['Microsoft Excel']
# excelPath = '/Applications/Microsoft\ Excel/Contents/MacOS/Microsoft\ Excel'
# jabrefOpen = ['JabRef'] '/Applications/JabRef/Contents/MacOS/JavaApplicationStub'
# jabrefPath = '/Applications/JabRef/Contents/MacOS/JavaApplicationStub'
# jaspOpen = ['JASP']
# jaspPath = 'Applications/JASP/Contents/MacOS/JASP'
# jmpOpen = ['JMP Pro 14']
# jmpPath = '/Applications/JMP\ Pro\ 14/Contents/MacOS/JMP'
# juliaOpen = ['Julia-1.2']
# juliaPath = '/Applications/Julia-1.2/Contents/MacOS/applet'
# jupyterPath = '/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter-notebook'
# RstudioPathOpen = ['RStudio']
# RstudioPath = '/Applications/RStudio/Contents/MacOS/Rstudio'
#
#
# # Local files
# localPDBfilePath = '$HOME/pdbFiles/'
# localEMAPfilePath = '$HOME/emapFiles/'
# localHKLfilePath = '$HOME/hklFiles/'
#
# # Image manipulation programs
# gimp = '/usr/local/bin/gimp'
# inkscapePath = '/opt/local/bin/inkscape'
# pptOpen = ['Microsoft PowerPoint']
# pptPath =  '/Applications/Microsoft\ Excel/Contents/MacOS/Microsoft\ PowerPoint'
#
#
# # Molecular graphics programs
# ccp4mgPath = '/Applications/ccp4-7.0/ccp4i2/Contents/MacOS/ccp4mg'
# chimeraOpen = ['Chimera']
# chimeraPath = '/Applications/Chimera/Contents/MacOS/chimera'
# cootPath = '/usr/local/bin/coot'
# jmolPath = 'java -jar /Applications/jars/jmol-14.29.52/Jmol.jar'
# yasaraOpen = ['YASARA']
# yasaraPath = '/Applications/YASARA/Contents/MacOS/yasara'
# vmdOpen = ['VMD194']
# vmdPath = '/Applications/VMD194/Contents/MacOS/startup.command'
#
#
# # Text editors
# atomOpen = ['atom']
# atomPath = '/usr/local/bin/atom'
# atomStart = ['start','atom']
# codeOpen = ['code']
# codePath = '/usr/local/bin/code'
# emacsOpen = ['emacs']
# emacsPath = '/opt/local/bin/emacs'
# geditOpen = ['gedit2.30.2']
# geditPath = '/Applications/gedit2.30.2/Contents/MacOS/gedit'
# jeditOpen = ['jEdit']
# jeditPath = '/usr/local/bin/jEdit'
# neovimPath = '$HOME/software/nvim-osx64/bin/nvim'
# nppOpen = ['Notepad++']
# nppPath = '/usr/local//Notepad++/Contents/MacOS/startwine'
# oniOpen = ['Oni']
# oniPath = '/Applications/Oni/Contents/MacOS/Oni'
# pdbedPath = 'java -jar /Applications/jars/PDB_Editor_FIX090203.jar'
# sublOpen = ['subl']
# sublimeText3Path = '/usr/local/bin/subl'
# textmatePath  = '/usr/local/bin/mate'
# textmateOpen = ['mate']
# vimPath = '/usr/local/bin/vim'
#
#
# #Terminals
# itermOpen = ['iTerm','-n','"`pwd`"']
# terminalOpen = ['Terminal','-n','"`pwd`"']
# x11termOpen = ['XQuartz','-n','"`pwd`"']
#
#
# # Web sites
# gcalPath = 'firefox https://calendar.google.com/calendar/r'
# gmailPath = 'firefox https://mail.google.com/mail/u/0/#inbox'
# webmailPath = 'firefox https://webmail.ouhsc.edu'
# weatherServicePath = 'firefox https://radar.weather.gov/radar.php?rid=TLX'
#
#
# wordProcessor
# wordOpen = ['open','-a','Microsoft Word']
#
# # Optional local mirror of the PDB
# local_mirror_divided = '/mnt/bio/db/pdb.divided'



# ###################### Windows Start Commands and PATHS to Applications###########################
#
# # Data analysis
# DBBrowserSQLiteStart = ['start','DBBrowserForSQLite']
# DBBrowserSQLitePath = '/Applications/DBBrowserForSQLite.exe/Contents/MacOS/DB\ Browser\ for\ SQLite'
# excelStart = ['start','Microsoft Excel']
# excelPath = '/Applications/Microsoft\ Excel.exe/Contents/MacOS/Microsoft\ Excel'
# jabrefStart = ['start','JabRef.exe'] '/Applications/JabRef.exe/Contents/MacOS/JavaApplicationStub'
# jabrefPath = '/Applications/JabRef.exe/Contents/MacOS/JavaApplicationStub'
# jaspStart = ['start','JASP.exe']
# jaspPath = 'Applications/JASP.exe/Contents/MacOS/JASP'
# jmpStart = ['start','JMP Pro 14.exe']
# jmpPath = '/Applications/JMP\ Pro\ 14.exe/Contents/MacOS/JMP'
# juliaStart = ['start','Julia-1.2.exe']
# juliaPath = '/Applications/Julia-1.2.exe/Contents/MacOS/applet'
# jupyterPath = '/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter-notebook'
# RstudioPathStart = ['start','RStudio']
# RstudioPath = '/Applications/RStudio.exe/Contents/MacOS/Rstudio'
#
#
# # Local files
# localPDBfilePath = '$HOME/pdbFiles/'
# localEMAPfilePath = '$HOME/emapFiles/'
# localHKLfilePath = '$HOME/hklFiles/'
#
# # Image manipulation programs
# gimp = '/usr/local/bin/gimp'
# inkscapePath = '/opt/local/bin/inkscape'
# pptStart = ['start','Microsoft PowerPoint']
# pptPath =  '/Applications/Microsoft\ Excel.exe/Contents/MacOS/Microsoft\ PowerPoint'
#
#
# # Molecular graphics programs
# ccp4mgPath = '/Applications/ccp4-7.0/ccp4i2.exe/Contents/MacOS/ccp4mg'
# chimeraStart = ['start','Chimera.exe']
# chimeraPath = '/Applications/Chimera.exe/Contents/MacOS/chimera'
# cootPath = '/usr/local/bin/coot'
# jmolPath = 'java -jar /Applications/jars/jmol-14.29.52/Jmol.jar'
# yasaraStart = ['start','YASARA.exe']
# yasaraPath = '/Applications/YASARA.exe/Contents/MacOS/yasara.exe'
# vmdStart = ['start','VMD194.exe']
# vmdPath = '/Applications/VMD194.exe/Contents/MacOS/startup.command'
#
# # Text editors
# atomStart = ['start','atom']
# atomPath = '/usr/local/bin/atom'
# atomStart = ['start','atom']
# codeStart = ['start','code']
# codePath = '/usr/local/bin/code'
# codeStart = ['start','code']
# emacsStart = ['start','emacs']
# emacsPath = '/opt/local/bin/emacs'
# geditStart = ['start','gedit2.30.2.exe']
# geditPath = '/Applications/gedit2.30.2.exe/Contents/MacOS/gedit'
# jeditStart = ['start','jEdit.exe']
# jeditPath = '/Applications/jEdit.exe/Contents/MacOS/jedit'
# neovimPath = '/Users/blaine/software/nvim-osx64/bin/nvim'
# nppStart = ['start','Notepad++.exe']
# nppPath = '/Applications/Notepad++.exe/Contents/MacOS/startwine'
# oniStart = ['start','Oni.exe']
# oniPath = '/Applications/Oni.exe/Contents/MacOS/Oni'
# pdbedPath = 'java -jar /Applications/jars/PDB_Editor_FIX090203.jar'
# sublStart = ['start','subl']
# sublimeText3Path = '/usr/local/bin/subl'
# textmatePath  = '/usr/local/bin/mate'
# textmateStart = ['start','mate']
# vimPath = '/opt/local/bin/vim'
#
#
# #Terminals
# itermStart = ['start','iTerm.exe','-n','"`pwd`"']
# terminalStart = ['start','Terminal','-n','"`pwd`"']
# x11termStart = ['start','XQuartz','-n','"`pwd`"']
#
# # Web sites
# gcalPath = 'open -a Safari.exe https://calendar.google.com/calendar/r'
# gmailPath = 'open -a Safari.exe https://mail.google.com/mail/u/0/#inbox'
# webmailPath = 'open -a Safari.exe https://webmail.ouhsc.edu'
# weatherServicePath = 'open -a Safari.exe https://radar.weather.gov/radar.php?rid=TLX'
#
# #wordProcessor
# wordStart = ['start','Microsoft Word.exe']

# # Optional local mirror of the PDB
# local_mirror_divided = '/mnt/bio/db/pdb.divided'


########################################################################################################
def AB(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Send search term or phrase to Amazon.com Books in default browser.

    USAGE:
    AB

    ARGUMENTS:
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
    url = abURL
    try:
        print("Sending",  searchTerm, " to Amazon.com Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to Amazon.com Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AB',AB)
    '''

    url = abURL
    try:
        print("Sending",  searchTerm, " to Amazon.com Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to Amazon.com Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AB',AB)


def AC(searchTerm='pymol'):
    ''' 
    DESCRIPTION:
    Send search term to Anaconda Cloud.



    USAGE:
    AC



    ARGUMENTS:
    NA


    EXAMPLE:
    AC



    MORE DETAILS:
    Send search term to Anaconda Cloud.
    May have to login first on default browser.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def AC(searchTerm='pymol'):
    url=anacondaURL
    try:
        print("Sending ",  searchTerm, " to  the Anaconda Cloud webpage in default browser. May have to register and login in first.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  Anaconda Cloud  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AC',AC)
    '''

    url=anacondaURL
    try:
        print("Sending ",  searchTerm, " to  the Anaconda Cloud webpage in default browser. May have to register and login in first.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  Anaconda Cloud  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AC',AC)


def ACA():
    ''' 
    DESCRIPTION:
    Open the American Crystallographic Association Annual Meeting webpage.

    USAGE:
    ACA

    ARGUMENTS:
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
    url=acaURL
    try:
        print("Trying to open American Crystallographic Association (ACA) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening ACA homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('ACA',ACA)
    '''

    url=acaURL
    try:
        print("Trying to open American Crystallographic Association (ACA) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening ACA homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('ACA',ACA)


def ALS():
    ''' 
    DESCRIPTION:
    Open website of the Advanced Light Source.

    USAGE:
    ALS

    ARGUMENTS:
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
    url=alsURL
    try:
        print("Trying to open Advanced Light Source (ALS) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening ALS homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('ALS',ALS)
    '''

    url=alsURL
    try:
        print("Trying to open Advanced Light Source (ALS) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening ALS homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('ALS',ALS)


def AO():
    ''' 
    DESCRIPTION:
    Commands to make ambient occlusion image like those in Qutemole.

    USAGE:
    AO

    ARGUMENTS:
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
cmd.extend('AO',AO)
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

    ARGUMENTS:
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
    cmd.color("gray20", "symbol c")
    cmd.color("gray90", "symbol h")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
cmd.extend("AOD",AOD)
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
    cmd.color("gray20", "symbol c")
    cmd.color("gray90", "symbol h")
    cmd.set("ray_opaque_background","on")
    cmd.ray()
cmd.extend("AOD",AOD)


def APS():
    ''' 
    DESCRIPTION:
    Open website of the Advanced Photon Source.

    USAGE:
    APS

    ARGUMENTS:
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
    url=apsURL
    try:
        print("Trying to open Advanced Photon Source (APS) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening APS homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('APS',APS)
    '''

    url=apsURL
    try:
        print("Trying to open Advanced Photon Source (APS) homepage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening APS homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('APS',APS)


def AX(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Send search term or phrase to arXiv.

    USAGE:
    AX

    ARGUMENTS:
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
    url = axURL
    searchType= '&searchtype=all&order=-announced_date_first&size=50'
    try:
        print("Sending ",  searchTerm, " to  Amazon.com Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm+searchType)
        print("Sent ",  searchTerm, " to  Amazon.com Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AX',AX)
    '''

    url = axURL
    searchType= '&searchtype=all&order=-announced_date_first&size=50'
    try:
        print("Sending ",  searchTerm, " to  Amazon.com Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm+searchType)
        print("Sent ",  searchTerm, " to  Amazon.com Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('AX',AX)


def BST():
    ''' 
    DESCRIPTION:
    G2G3/U9U8 base step , PDB code 4PCO. 

    USAGE:
    BST

    ARGUMENTS:
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
    delete all;fetch 4PCO, type=pdb, async=0;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));hide cartoon;set valence, off;remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view(-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. dd24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png
    PYTHON CODE:
def BST():
    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb')
    cmd.select('G2G3', '( ((resi 2 or resi 3) and chain A)or ((resi 8 or resi 9) and chain B) )')
    cmd.hide('cartoon')
    cmd.set('valence', 'off')

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
    '''

    cmd.reinitialize()
    cmd.fetch('4PCO', type='pdb')
    cmd.select('G2G3', '( ((resi 2 or resi 3) and chain A)or ((resi 8 or resi 9) and chain B) )')
    cmd.hide('cartoon')
    cmd.set('valence', 'off')

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

    ARGUMENTS:
    None
    EXAMPLE:
    BU

    MORE DETAILS:
    Commands to make biological unit. 
    Requires a pdb file rather than the cif file returned by the fetch command. 
    There are other ways of displaying the biological unit in PyMOL including downloading *.pdb1 from the PDB. 
    Depends on the quat3.py script by Thomas Holder.

    Type 'help BU' to see this documentation
    printed to the command history window. 

    Select from the command history individual 
    lines of code to build a new script. 

    Select the hortizontal script at the bottom if retaining 
    most of the commands in your new script. 

    Copy and paste onto the command line below.
    Works only with the command line immediately under the     command history window at the top of the gui.

    VERTICAL PML SCRIPT:
    run ~/Scripts/PyMOLScripts/quat3.py; 
    quat 

    HORIZONTAL PML SCRIPT:
    run ~/Scripts/PyMOLScripts/quat3.py; quat 
    PYTHON CODE:
def BU():
    cmd.do('run $HOME/Scripts/PyMOLScripts/quat3.py')
    cmd.do('quat')
cmd.extend('BU',BU)
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


    ARGUMENTS:
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
    cmd.set("opaque_background", "on")
    cmd.set('ray_trace_mode', '2'); 
    cmd.bg_color('white'); 
    cmd.set('antialias', '2'); 
    cmd.ray('600','600'); 
    cmd.png('test.png')
cmd.extend('BW', BW)
    '''

    cmd.show_as("cartoon", "all"); 
    cmd.hide('lines'); 
    cmd.hide('nonbonded'); 
    # black and white cartoon; 
    cmd.set("opaque_background", "on")
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

    ARGUMENTS:
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
    url = bxURL
    try:
        print("Sending", searchTerm,"to bioRxiv in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent", searchTerm," to bioRxiv in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('BX',BX)
    '''

    url = bxURL
    try:
        print("Sending", searchTerm,"to bioRxiv in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent", searchTerm," to bioRxiv in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('BX',BX)


def CB():
    ''' 
    DESCRIPTION:
    Runs Jared Sampson's script "colorblindfriendly.py". 

    USAGE:
    CB

    ARGUMENTS:
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

cmd.extend('CB',CB)
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

    ARGUMENTS:
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

cmd.extend('CBSS',CBSS)
    '''

    cmd.do('CB')
    cmd.show_as('cartoon')
    cmd.color('cb_red', 'ss H')
    cmd.color('cb_yellow', 'ss S')
    cmd.color('cb_green', 'ss L+')

cmd.extend('CBSS',CBSS)


def CHESS():
    ''' 
    DESCRIPTION:
    Open the website of CHESS.


    USAGE:
    CHESS


    ARGUMENTS:
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
    url=chessURL
    try:
        print("Opening the website of CHESS.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of CHESS.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('CHESS',CHESS)
    '''

    url=chessURL
    try:
        print("Opening the website of CHESS.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of CHESS.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('CHESS',CHESS)


def CR():
    ''' 
    DESCRIPTION:
    Commands to make colored filled-ring cartoon of nucleic acids.

    USAGE:
    CR

    ARGUMENTS:
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

cmd.extend('CR',CR)
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
    Commands to color ribbon or cartoon representations of proteins by secondary structure. 

    USAGE:
    CSS

    ARGUMENTS:
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

cmd.extend('CSS',CSS)
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

    ARGUMENTS:
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

cmd.extend('DU',DU) 
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

    ARGUMENTS:
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
    url=emdbURL
    try:
        print("Opening the website of the Electron Microscopy Data Bank.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of the Electron Microscopy Data Bank.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('EMDB',EMDB)
    '''

    url=emdbURL
    try:
        print("Opening the website of the Electron Microscopy Data Bank.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of the Electron Microscopy Data Bank.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('EMDB',EMDB)


def EP():
    ''' 
    DESCRIPTION:
    EasyPyMOL github site.

    USAGE:
    EP

    ARGUMENTS:
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
    url=epURL
    try:
        print("Opening the EasyPyMOL github page in new tab of default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the EasyPyMOL github page in new tab of default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('EP',EP)
    '''

    url=epURL
    try:
        print("Opening the EasyPyMOL github page in new tab of default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the EasyPyMOL github page in new tab of default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('EP',EP)


def FR():
    ''' 
    DESCRIPTION:
    Make filled-ring cartoon of nucleic acids.

    USAGE:
    FR

    ARGUMENTS:
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

cmd.extend('FR',FR)
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

    ARGUMENTS:
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
    url = gbURL
    try:
        print("Sending",  searchTerm, " to Google Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to Google Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GB',GB)
    '''

    url = gbURL
    try:
        print("Sending",  searchTerm, " to Google Books in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to Google Books in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GB',GB)


def GGT():
    ''' 
    DESCRIPTION:
    WT human gamma glutamyl transpeptidase at 1.67 Angstrom resolution as cartoon. PDB Code 4gdx.


    USAGE:
    GGT

    ARGUMENTS:
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

cmd.extend('GGT',GGT)
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

    ARGUMENTS:
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
    url = ghURL
    try:
        print("Sending",  searchTerm, " to GitHub in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to GitHubs in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GH',GH)
    '''

    url = ghURL
    try:
        print("Sending",  searchTerm, " to GitHub in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to GitHubs in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GH',GH)


def GO(searchTerm="pymol",numHits="200"):
    ''' 
    DESCRIPTION:
    Send search term or phrase Google in default browser.

    USAGE:
    GO

    ARGUMENTS:
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
    url = goURL
    nhits= '&num='
    try:
        print("Sending ", searchTerm," search term or phrase to Google in default browser along with number of hits to be returned.");
        # URL must be in single quotes
        webbrowser.open_new_tab(goURL+searchTerm + nhits + str(numHits))
        print("Sent ", searchTerm," to Google in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GO',GO)
    '''

    url = goURL
    nhits= '&num='
    try:
        print("Sending ", searchTerm," search term or phrase to Google in default browser along with number of hits to be returned.");
        # URL must be in single quotes
        webbrowser.open_new_tab(goURL+searchTerm + nhits + str(numHits))
        print("Sent ", searchTerm," to Google in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GO',GO)


def GS(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Send search term or phrase to Google Scholar in default browser.

    USAGE:
    GS

    ARGUMENTS:
    searchTerm, searchTerm, searchTerm

    EXAMPLE:
    GS molecular replacement software

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
    url = gsURL
    try:
        print("Sending ", searchTerm," to Google Sholar in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ", searchTerm," to Google Sholar in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GS',GS)
    '''

    url = gsURL
    try:
        print("Sending ", searchTerm," to Google Sholar in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ", searchTerm," to Google Sholar in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GS',GS)


def GU():
    ''' 
    DESCRIPTION:
    10-mer dsRNA with 8 GU wobble base pairs.


    USAGE:
    GU

    ARGUMENTS:
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

cmd.extend('GU',GU)
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
    GV searchTerm

    ARGUMENTS:
    searchTerm  (Do not have to be in quotes. Can be multiple terms.)
    EXAMPLE:
    GV pymol movie

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
    url = gvURL
    try:
        print("Sending",  searchTerm, "to Google Video in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(goURL+searchTerm)
        print("Sent ",  searchTerm, " to Google Video in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GV',GV)
    '''

    url = gvURL
    try:
        print("Sending",  searchTerm, "to Google Video in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(goURL+searchTerm)
        print("Sent ",  searchTerm, " to Google Video in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('GV',GV)


def HH():
    ''' 
    DESCRIPTION:
    Hide hydrogen atoms of currently visible molecular objects.

    USAGE:
    HH

    ARGUMENTS:
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
cmd.extend('HH',HH)
    '''

    cmd.hide('everything', 'name H*')
cmd.extend('HH',HH)


def IPM(searchTerms = [], *args):
    ''' 
    DESCRIPTION:
    Read list of search terms and submit each term to PubMed in a separate browser tab.

    USAGE:
    IPM

    ARGUMENTS:
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
    url = ipmURL
    try:
        print('Sending',  searchTerm, ' to Pubmed and display list of search results in separate tabs of the default brower.')
        for term in searchTerms:
            t0 = time.time()
            sterm = str(term)
            webbrowser.open_new_tab(url+sterm)
            response_delay = time.time() - t0
            time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
            print('Finished searching PubMed for ', sterm, '.')
        print('Finished searching PubMed for  ',  searchTerms, '.') 
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('IPM',IPM)
    '''

    url = ipmURL
    try:
        print('Sending',  searchTerm, ' to Pubmed and display list of search results in separate tabs of the default brower.')
        for term in searchTerms:
            t0 = time.time()
            sterm = str(term)
            webbrowser.open_new_tab(url+sterm)
            response_delay = time.time() - t0
            time.sleep(10*response_delay)  # wait 10x longer than it took them to respond
            print('Finished searching PubMed for ', sterm, '.')
        print('Finished searching PubMed for  ',  searchTerms, '.') 
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('IPM',IPM)


def IUCR(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Open website of the IUCr Journals.

    USAGE:
    IUCR

    ARGUMENTS:
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
    url=iucrURL
    try:
        print("Opening the IUCr journals webpage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the IUCr journals webpage.");
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('IUCR',IUCR)
    '''

    url=iucrURL
    try:
        print("Opening the IUCr journals webpage.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the IUCr journals webpage.");
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('IUCR',IUCR)


def JASP():
    ''' 
    DESCRIPTION:
    Open JASP from within PyMOL.

    USAGE:
    JASP

    ARGUMENTS:
    None
    EXAMPLE:
    JASP

    MORE DETAILS:
    Open JASP from within PyMOL.  The is a data analysis program 
    that can do Bayesian and frequentist statistics in parallel.

    VERTICAL PML SCRIPT:
        subprocess.call(jaspOpen);
    return

    HORIZONTAL PML SCRIPT:
        subprocess.call(jaspOpen);return
    PYTHON CODE:
def JASP():
    try:
        print("Opening the JASP.");
        subprocess.check_output(jaspOpen)
        print("Success opening JASP.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jaspOpen'. \n  Or use 'jaspPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('JASP',JASP)
    '''

    try:
        print("Opening the JASP.");
        subprocess.check_output(jaspOpen)
        print("Success opening JASP.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jaspOpen'. \n  Or use 'jaspPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('JASP',JASP)


def JM():
    ''' 
    DESCRIPTION:
    Open the Jmol wiki.


    USAGE:
    JM

    ARGUMENTS:
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
    url=jmURL
    try:
        print("Opening the Jmol wiki.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the Jmol wiki.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found


cmd.extend('JM',JM)
    '''

    url=jmURL
    try:
        print("Opening the Jmol wiki.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the Jmol wiki.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found


cmd.extend('JM',JM)


def JMP():
    ''' 
    DESCRIPTION:
    Open the JMP from within PyMOL. 

    USAGE:
    JMP

    ARGUMENTS:
    None
    EXAMPLE:
    JMP

    MORE DETAILS:
    Open the JMP from within PyMOL. 


    VERTICAL PML SCRIPT:
    arg = jmpPath;
    subprocess.call(arg,shell=True);
    return

    HORIZONTAL PML SCRIPT:
    arg = jmpPath;subprocess.call(arg,shell=True);return

    PYTHON CODE:
def JMP():
    try:
        print("Opening the JMP.");
        subprocess.check_output(jmpOpen)
        print("Success opening JMP.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jmpOpen'. \n  Or use 'jmpPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('JMP',JMP)
    '''

    try:
        print("Opening the JMP.");
        subprocess.check_output(jmpOpen)
        print("Success opening JMP.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jmpOpen'. \n  Or use 'jmpPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('JMP',JMP)


def LBSF():
    ''' 
    DESCRIPTION:
    Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.

    USAGE:
    LBSF

    ARGUMENTS:
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
    url=lbsfURL
    try:
        print("Opening the website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('LBSF',LBSF)
    '''

    url=lbsfURL
    try:
        print("Opening the website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('LBSF',LBSF)


def LBST():
    ''' 
    DESCRIPTION:
    G2G3/U9U8 base step , PDB code 4PCO. 

    USAGE:
    LBST

    ARGUMENTS:
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
    hide cartoon;
    set valence, off;
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
    delete all;load 4PCO.pdb;hide cartoon;set valence, off;select G2G3, ( ((resi 2 or resi 3) and chain A) or ((resi 8 or resi 9) and chain B));remove not G2G3;bg_color white;show sticks;set stick_radius=0.14;set stick_ball, on;set stick_ball_ratio,1.9;set_view (-0.75,0.09,0.66,-0.2,0.92,-0.35,-0.64,-0.39,-0.67,-0.0,-0.0,-43.7,7. 24,9.55,11.78,29.46,57.91,-20.0);remove name H*;select carbon1, element C and (resi 3 or resi 8);select carbon2, element C and (resi 2 or resi 9);color gray70, carbon1;color gray10, carbon2;show sticks;space cmyk;distance hbond1, /4PCO//B/U`9/N3,/4PCO//A/G`2/O6;distance hbond2, /4PCO//B/U`9/O2,/4PCO//A/G`2/N1;distance hbond3, /4PCO//A/U`3/N3,/4PCO//B/G`8/O6;distance hbond4, /4PCO//A/U`3/O2,/4PCO//B/G`8/N1;color black, hbond1;color black, hbond2;color gray70, hbond3;color gray70, hbond4;show nb_spheres;set nb_spheres_size, 0.35;hide labels;ray 1600,1000;png 4PCO.png
    PYTHON CODE:
def LBST():
    cmd.reinitialize()
    cmd.load(localPDBfilePath + '4pco.pdb')
    cmd.hide('cartoon')
    cmd.set('valence', 'off')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '4pco.pdb')
    cmd.hide('cartoon')
    cmd.set('valence', 'off')
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

    ARGUMENTS:
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

cmd.extend('LG',LG)
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '4gdx.pdb')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '4gdx.pdb')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '4PCO.pdb')
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
    '''

    cmd.reinitialize();
    cmd.load(localPDBfilePath + '4PCO.pdb')
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

    ARGUMENTS:
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
    cmd.load(localEMAPfilePath + '4dgr2FoFc.ccp4')
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
    '''

    cmd.reinitialize()
    cmd.load('4dgr.pdb')
    cmd.load(localEMAPfilePath + '4dgr2FoFc.ccp4')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '4dgr.pdb')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '4dgr.pdb')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '3nd4.pdb');
    cmd.hide('cartoon');
    cmd.set('valence','off');
    cmd.show('sticks');
    cmd.do('run $HOME/Scripts/PyMOLScripts/quat3.py')
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
    '''

    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.load(localPDBfilePath + '3nd4.pdb');
    cmd.hide('cartoon');
    cmd.set('valence','off');
    cmd.show('sticks');
    cmd.do('run $HOME/Scripts/PyMOLScripts/quat3.py')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '3fa0.pdb')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '3fa0.pdb')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '3nd3.pdb')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '3nd3.pdb')
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

    ARGUMENTS:
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
    cmd.load(localPDBfilePath + '3nd4.pdb')
    cmd.remove('name H*')
    cmd.hide('everything')
    # cmd.do('run $HOME/mg18OU/quat.py')
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
    '''

    cmd.reinitialize()
    cmd.load(localPDBfilePath + '3nd4.pdb')
    cmd.remove('name H*')
    cmd.hide('everything')
    # cmd.do('run $HOME/mg18OU/quat.py')
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

    ARGUMENTS:
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

cmd.extend('MA',MA)
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

    ARGUMENTS:
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

cmd.extend('MB',MB)
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

    ARGUMENTS:
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

cmd.extend('MC',MC)
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

    ARGUMENTS:
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
    url=ouMCLURL
    try:
        print("Opening the website of Macromolecular Crystallography Laboratory at the University of Oklahoma.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of Macromolecular Crystallography Laboratory at the University of Oklahoma.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MCL',MCL)
    '''

    url=ouMCLURL
    try:
        print("Opening the website of Macromolecular Crystallography Laboratory at the University of Oklahoma.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of Macromolecular Crystallography Laboratory at the University of Oklahoma.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MCL',MCL)


def MG():
    ''' 
    DESCRIPTION:
    Open website of the OUHSC molecular graphics course.

    USAGE:
    MG

    ARGUMENTS:
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
    url=molgrURL
    try:
        print("Opening the website with molecular graphic links.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened website with molecular graphic links..")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MG',MG)
    '''

    url=molgrURL
    try:
        print("Opening the website with molecular graphic links.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened website with molecular graphic links..")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MG',MG)


def MGW():
    ''' 
    DESCRIPTION:
    Open Wikipedia webpage about molecular graphics.

    USAGE:
    MGW

    ARGUMENTS:
    NA
    EXAMPLE:
    MGW

    MORE DETAILS:
    Open Wikipedia webpage about molecular graphics.

    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def MGW():
    url=molgrwikiURL
    try:
        print("Opening the Wikipedia webpage about molecular graphics.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the Wikipedia webpage about molecular graphics.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MGW',MGW)
    '''

    url=molgrwikiURL
    try:
        print("Opening the Wikipedia webpage about molecular graphics.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the Wikipedia webpage about molecular graphics.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('MGW',MGW)


def MM(searchTerm='pymol'):
    ''' 
    DESCRIPTION:
    Send search term to search for manuscripts in pymolshortcuts.

    USAGE:
    MM

    ARGUMENTS:
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

cmd.extend('MM',MM)
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

    ARGUMENTS:
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

cmd.extend('N9',N9)
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

    ARGUMENTS:
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
    cmd.do('split_states 3nd4');
    cmd.hide('cartoon');
    cmd.hide('spheres');
    cmd.set('valence','off');
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'elem H*');
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
    cmd.distance('hbond1', '3nd4_0001 and resi 119 and name O', '3nd4_0001 and resi 1 and name OP2');
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
    '''

    cmd.reinitialize();
    cmd.viewport('900','600');
    cmd.fetch('3nd4', type='pdb1');
    cmd.do('split_states 3nd4');
    cmd.hide('cartoon');
    cmd.hide('spheres');
    cmd.set('valence','off');
    cmd.show('sticks');
    cmd.set('stick_radius', '0.125');
    cmd.hide('everything', 'elem H*');
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
    cmd.distance('hbond1', '3nd4_0001 and resi 119 and name O', '3nd4_0001 and resi 1 and name OP2');
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



    ARGUMENTS:
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
    url=ndbURL
    try:
        print("Opening the website of the Nucleic Acid Database.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened website of the Nucleic Acid Database.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('NDB',NDB)
    '''

    url=ndbURL
    try:
        print("Opening the website of the Nucleic Acid Database.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened website of the Nucleic Acid Database.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('NDB',NDB)


def NSLSII():
    ''' 
    DESCRIPTION:
    Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory.

    USAGE:
    NSLSII

    ARGUMENTS:
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
    url=nslsIIURL
    try:
        print("Opening the website of the National Synchrotron Light Source II ,NSLSII, at Brookhaven National Laboratory..");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of the National Synchrotron Light Source II ,NSLSII, at Brookhaven National Laboratory.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('NSLSII',NSLSII)
    '''

    url=nslsIIURL
    try:
        print("Opening the website of the National Synchrotron Light Source II ,NSLSII, at Brookhaven National Laboratory..");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of the National Synchrotron Light Source II ,NSLSII, at Brookhaven National Laboratory.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('NSLSII',NSLSII)


def PDB(searchTerm="3fa0",numHits="5"):
    ''' 
    DESCRIPTION:
    Submit a search term to the Protein Data Bank.

    USAGE:
    PDB

    ARGUMENTS:
    searchTerm
    EXAMPLE:
    PBB 3fa0



    MORE DETAILS:
    Submit a pdbcode to the Protein Data Bank and get back the webpage for the orrespondong structure.

    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def PDB(searchTerm="3fa0",numHits="5"):
    url = pdbURL
    try:
        print("Sending",  searchTerm, " to the PBD webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to the PBD webpage in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PDB',PDB)
    '''

    url = pdbURL
    try:
        print("Sending",  searchTerm, " to the PBD webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to the PBD webpage in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PDB',PDB)


def PE(selection):
    ''' 
    DESCRIPTION:
    Apply pearl effect about selection cation or anion.

    USAGE:
    PE selection



    ARGUMENTS:
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
    cmd.create('magnesium2', selection)
    cmd.show('spheres', 'magnesium1')
    cmd.show('spheres', 'magnesium2')
    cmd.set('sphere_transparency', '0.4', 'magnesium2')
    cmd.set('sphere_scale', '1.1', 'magnesium2')

cmd.extend('PE',PE)
    '''

    cmd.select('magnesium1',selection)
    cmd.create('magnesium2', selection)
    cmd.show('spheres', 'magnesium1')
    cmd.show('spheres', 'magnesium2')
    cmd.set('sphere_transparency', '0.4', 'magnesium2')
    cmd.set('sphere_scale', '1.1', 'magnesium2')

cmd.extend('PE',PE)


def PE2(selection):
    ''' 
    DESCRIPTION:
    Apply alternative pearl effect about selected cation or anion.

    USAGE:
    PE2 selection



    ARGUMENTS:
    selection of one catoin or anion
    EXAMPLE:
    PE2 /3nd4//A/NA`19/NA

    MORE DETAILS:
    Apply alternative pearl effect about selected cation or anion.
    This effect shows a transparent sphere (1.0 sphere_scale) around an inner opaque sphere (.35 sphere_scale).
    Must supply the <selection> of a cation or anion. 


    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def PE2(selection):
    cmd.select('magnesium1',selection)
    cmd.create('magnesium2',selection)
    cmd.show('spheres', 'magnesium1')
    cmd.show('spheres', 'magnesium2')
    cmd.set('sphere_transparency', '0.0', 'magnesium1')
    cmd.set('sphere_transparency', '0.5', 'magnesium2')
    cmd.set('sphere_scale', '0.35', 'magnesium1')
    cmd.set('sphere_scale', '1.0', 'magnesium2')

cmd.extend('PE2',PE2)
    '''

    cmd.select('magnesium1',selection)
    cmd.create('magnesium2',selection)
    cmd.show('spheres', 'magnesium1')
    cmd.show('spheres', 'magnesium2')
    cmd.set('sphere_transparency', '0.0', 'magnesium1')
    cmd.set('sphere_transparency', '0.5', 'magnesium2')
    cmd.set('sphere_scale', '0.35', 'magnesium1')
    cmd.set('sphere_scale', '1.0', 'magnesium2')

cmd.extend('PE2',PE2)


def PM(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Send search term or phrase to PubMed.

    USAGE:
    PM

    ARGUMENTS:
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
    url = pmURL
    try:
        print("Sending",  searchTerm, " to the PubMed webpage in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  the PubMed webpagein default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PM',PM)
    '''

    url = pmURL
    try:
        print("Sending",  searchTerm, " to the PubMed webpage in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  the PubMed webpagein default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PM',PM)


def PML(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:
    Submit a search term to the PyMOL Users Mail Service.

    USAGE:
    PML

    ARGUMENTS:
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
    url = pmlURL
    try:
        print("Sending",  searchTerm, " to the PyMOL mailing list in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to the PubMed webpagein default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PML',PML)
    '''

    url = pmlURL
    try:
        print("Sending",  searchTerm, " to the PyMOL mailing list in default browser.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to the PubMed webpagein default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PML',PML)


def PPC():
    ''' 
    DESCRIPTION:
    Open the website of the Protein Production Facility at the University of Oklahoma in Norman.

    USAGE:
    PPC

    ARGUMENTS:
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
    url=ppcURL
    try:
        print("Opening the website of the Protein Production Facility at the University of Oklahoma in Norman.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of the Protein Production Facility at the University of Oklahoma in Norman.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('PPC',PPC)
    '''

    url=ppcURL
    try:
        print("Opening the website of the Protein Production Facility at the University of Oklahoma in Norman.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the website of the Protein Production Facility at the University of Oklahoma in Norman.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('PPC',PPC)


def PS():
    ''' 
    DESCRIPTION:
    Open the home page of the Protein Soceity.

    USAGE:
    PS

    ARGUMENTS:
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
    url=psURL
    try:
        print("Opening the home page of the Protein Soceity.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the home page of the Protein Soceity.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('PS',PS)
    '''

    url=psURL
    try:
        print("Opening the home page of the Protein Soceity.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the home page of the Protein Soceity.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('PS',PS)


def PW(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:
    Submit search of the PyMOL Wiki. 

    USAGE:
    PW

    ARGUMENTS:
    searchTerm
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
    url = pymolURL
    try:
        print("Sending ",  searchTerm, " to  the PyMOL Wiki webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  PyMOL Wiki  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PW',PW)
    '''

    url = pymolURL
    try:
        print("Sending ",  searchTerm, " to  the PyMOL Wiki webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to  PyMOL Wiki  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('PW',PW)


def RG(searchTerm='best molecular graphics program'):
    ''' 
    DESCRIPTION:
    Submit a search query of Research Gate. 

    USAGE:
    RG

    ARGUMENTS:
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
    url = researchGateURL
    try:
        print("Sending",  searchTerm, " to  the Research Gate webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to Research Gate  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('RG',RG)
    '''

    url = researchGateURL
    try:
        print("Sending",  searchTerm, " to  the Research Gate webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to Research Gate  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('RG',RG)


def RS():
    ''' 
    DESCRIPTION:
    Open the homepage of the RNA Society.

    USAGE:
    RS

    ARGUMENTS:
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
    url=rsURL
    try:
        print("Opening the homepage of the RNA Society.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the homepage of the RNA Society.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('RS',RS)
    '''

    url=rsURL
    try:
        print("Opening the homepage of the RNA Society.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the homepage of the RNA Society.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('RS',RS)


def RStudio():
    ''' 
    DESCRIPTION:
    Open Rstudio GUI.

    USAGE:
    RStudio

    ARGUMENTS:
    None
    EXAMPLE:
    RStudio

    MORE DETAILS:
    Open RStudio from within PyMOL.


    VERTICAL PML SCRIPT:
    subprocess.call(RStudioOpen); 
return   

    HORIZONTAL PML SCRIPT:
    subprocess.call(RStudioOpen); return
    PYTHON CODE:
def RStudio():
    try:
        print("Opening the RStudio.");
        subprocess.check_output(RStudioOpen)
        print("Success opening RStudio.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'RStudioOpen'. \n  Or use ''RStudioPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('RStudio',RStudio)
    '''

    try:
        print("Opening the RStudio.");
        subprocess.check_output(RStudioOpen)
        print("Success opening RStudio.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'RStudioOpen'. \n  Or use ''RStudioPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('RStudio',RStudio)


def SAXS():
    ''' 
    DESCRIPTION:
    Open the webpage of SAXS links at OUHSC. 

    USAGE:
    SAXS

    ARGUMENTS:
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
    url=saxsURL
    try:
        print("Opening the webpage of SAXS links at OUHSC.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SAXS links at OUHSC.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('SAXS',SAXS)
    '''

    url=saxsURL
    try:
        print("Opening the webpage of SAXS links at OUHSC.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SAXS links at OUHSC.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 

cmd.extend('SAXS',SAXS)


def SC():
    ''' 
    DESCRIPTION:
    Print to screen list of the shortcuts that are available in the script pymolshortcuts.py. 

    USAGE:
    SC

    ARGUMENTS:
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

    PYTHON CODE:
def SC():
    print(SC.__doc__)

cmd.extend('SC',SC)
    '''

    print(SC.__doc__)

cmd.extend('SC',SC)


def SD(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Submit a search term to Science Direct.

    USAGE:
    SD

    ARGUMENTS:
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
    url1 = scienceDirectURL1
    url2 = scienceDirectURL2
    try:
        print("Sending",  searchTerm, " to  the Science Direct webpage in the default browser.");
        webbrowser.open_new_tab(url1+searchTerm+url2)
        print('Sent', searchTerm, 'to the Science Direct webpage in the default browser.')
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SD',SD)
    '''

    url1 = scienceDirectURL1
    url2 = scienceDirectURL2
    try:
        print("Sending",  searchTerm, " to  the Science Direct webpage in the default browser.");
        webbrowser.open_new_tab(url1+searchTerm+url2)
        print('Sent', searchTerm, 'to the Science Direct webpage in the default browser.')
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SD',SD)


def SF(searchTerm='pymol'):
    ''' 
    DESCRIPTION:
    Send search term to sourceforge.

    USAGE:
    SF

    ARGUMENTS:
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
    url = sourceForgeURL
    try:
        print("Sending ",  searchTerm, " to  the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to sourceforge  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SF',SF)
    '''

    url = sourceForgeURL
    try:
        print("Sending ",  searchTerm, " to  the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, " to sourceforge  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SF',SF)


def SH():
    ''' 
    DESCRIPTION:
    Show hydrogen atoms of currently visible molecular objects.

    USAGE:
    SH

    ARGUMENTS:
    None
    EXAMPLE:
    SH

    MORE DETAILS:
    Show hydrogen atoms of currently visible molecular objects.
    VERTICAL PML SCRIPT:
    show sticks, element H
    HORIZONTAL PML SCRIPT:
    show sticks, element H
    PYTHON CODE:
def SH():
    cmd.show('sticks', 'element H')
cmd.extend('SH',SH)
    '''

    cmd.show('sticks', 'element H')
cmd.extend('SH',SH)


def SO(searchTerm="3d_pdf"):
    ''' 
    DESCRIPTION:
    Submit a search term to Stackoverflow.


    USAGE:
    SO



    ARGUMENTS:
    searchTerm


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
    url = stackoverflowURL
    try:
        print("Sending",  searchTerm, " to the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to  stackoverflow  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SO',SO)
    '''

    url = stackoverflowURL
    try:
        print("Sending",  searchTerm, " to the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url+searchTerm)
        print("Sent ",  searchTerm, "to  stackoverflow  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SO',SO)


def SP(searchTerm="pymol"):
    ''' 
    DESCRIPTION:
    Submit a search term to Springer Books

    USAGE:
    SP

    ARGUMENTS:
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
    url1 = springerBooksURL1
    url2 = springerBooksURL2
    try:
        print("Sending ",  searchTerm, " to  the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url1+searchTerm+url2)
        print("Sent ",  searchTerm, " to sourceforge  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SP',SP)
    '''

    url1 = springerBooksURL1
    url2 = springerBooksURL2
    try:
        print("Sending ",  searchTerm, " to  the sourceforge webpage in default browser.");
        webbrowser.open_new_tab(url1+searchTerm+url2)
        print("Sent ",  searchTerm, " to sourceforge  in default browser.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('SP',SP)


def SSRLSMB():
    ''' 
    DESCRIPTION:
    Open the webpage of SSRL Structural Molecular Biology.



    USAGE:
    SSRLSMB



    ARGUMENTS:
    NA
    EXAMPLE:
    SSRLSMB



    MORE DETAILS:
    Open the webpage of SSRL Structural Molecular Biology.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def SSRLSMB():
    url=ssrlsmbURL
    try:
        print("Opening the webpage of SSRL Structural Molecular Biology.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SSRL Structural Molecular Biology.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 



cmd.extend('SSRLSMB',SSRLSMB)
    '''

    url=ssrlsmbURL
    try:
        print("Opening the webpage of SSRL Structural Molecular Biology.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SSRL Structural Molecular Biology.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 



cmd.extend('SSRLSMB',SSRLSMB)


def SSURF():
    ''' 
    DESCRIPTION:
    Open the webpage of the Society for Science at User Research Facilities (SSURF).


    USAGE:
    SSURF



    ARGUMENTS:
    NA
    EXAMPLE:
    SSURF



    MORE DETAILS:
    Open the webpage of the Society for Science at User Research Facilities (SSURF).
    SSURF is nonprofit organization in the US that serves as the
    nexus for users and user executive committees at
    national laboratories. SSURF is not a lobbying organization, but
    it helps organize visits to Congress to educate legislators about
    the importance of national laboratories in science. Membership
    is free of students. The annual fee is nominal for PIs. 


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def SSURF():
    url=ssurfURL
    try:
        print("Opening the webpage of the Society for Science at User Research Facilities (SSURF).");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Society for Science at User Research Facilities (SSURF).")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('SSURF',SSURF)
    '''

    url=ssurfURL
    try:
        print("Opening the webpage of the Society for Science at User Research Facilities (SSURF).");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Society for Science at User Research Facilities (SSURF).")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('SSURF',SSURF)


def SciPy19():
    ''' 
    DESCRIPTION:
    Open the SciPy 2019 YouTube Channel.

    USAGE:
    SciPy19

    ARGUMENTS:
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
    url=scipyURL
    try:
        print("Opening the webpage of the SciPy19 meeting.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the SciPy19 meeting.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('SciPy19',SciPy19)
    '''

    url=scipyURL
    try:
        print("Opening the webpage of the SciPy19 meeting.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the SciPy19 meeting.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('SciPy19',SciPy19)


def T4L():
    ''' 
    DESCRIPTION:
    WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0. 

    USAGE:
    T4L

    ARGUMENTS:
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

cmd.extend('T4L',T4L)
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

    ARGUMENTS:
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

cmd.extend('U8',U8)
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

    ARGUMENTS:
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

cmd.extend('WC8',WC8)
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

    ARGUMENTS:
    None
    EXAMPLE:
    atom

    MORE DETAILS:
    Open the text editor Atom from within PyMOL.


    VERTICAL PML SCRIPT:
    subprocess.call(atomOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(atomOpen);return

    PYTHON CODE:
def atom(fileName="test.pml"):
    try:
        print("Opeing the text editor atom. Please be patient. It starts slowly.")
        subprocess.check_output(atomOpen)
        print("Success opening atom")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'atomOpen'. \n Or use 'atomPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('atom',atom)
    '''

    try:
        print("Opeing the text editor atom. Please be patient. It starts slowly.")
        subprocess.check_output(atomOpen)
        print("Success opening atom")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'atomOpen'. \n Or use 'atomPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('atom',atom)


def bbedit(fileName="test.pml"):
    ''' 
    DESCRIPTION:
    Open file with the text editor bbedit from within PyMOL. 

    USAGE:
    bbedit

    ARGUMENTS:
    None
    EXAMPLE:
    bbedit

    MORE DETAILS:
    Open file with the text editor bbedit from within PyMOL. 
    Adjust the path as needed for your system.
    Only available for Mac OS.


    VERTICAL PML SCRIPT:
    subprocess.call(bbeditOpen);
    return
    HORIZONTAL PML SCRIPT:
    subprocess.call(bbeditOpen);return

    PYTHON CODE:
def bbedit(fileName="test.pml"):
    try:
        print("Opening the molecular graphics program bbedit.");
        subprocess.check_output(bbeditOpen)
        print("Success opening bbedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the bbeditOpen'. \n  Or use 'bbeditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('bbedit',bbedit)
    '''

    try:
        print("Opening the molecular graphics program bbedit.");
        subprocess.check_output(bbeditOpen)
        print("Success opening bbedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the bbeditOpen'. \n  Or use 'bbeditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('bbedit',bbedit)


def biocat():
    ''' 
    DESCRIPTION:
    Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.



    USAGE:
    biocat


    ARGUMENTS:
    NA

    EXAMPLE:
    biocat


    MORE DETAILS:
    Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def biocat():
    url=biocatURL
    try:
        print("Opening the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening the BIOCAT biological SAXS beamline at the Advanced Photon Source homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('biocat',biocat)
    '''

    url=biocatURL
    try:
        print("Opening the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening the BIOCAT biological SAXS beamline at the Advanced Photon Source homepage.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('biocat',biocat)


def bs(selection='all'):
    ''' 
    DESCRIPTION:
    bs creates a ball and stick representation of an object. 

    USAGE:
    bs selection

    ARGUMENTS:
    selection

    EXAMPLE:
    bs 3nd3

    MORE DETAILS:
    bs creates a ball and stick representation of an object. 
    Bondi VDW values added below to override default Pymol settings.
    From https://gist.githubusercontent.com/bobbypaton/1cdc4784f3fc8374467bae5eb410edef/raw/9995d51d6a64b8bcf01590c944cc38059b2f8d7f/pymol_style.py



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
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("slate","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")

cmd.extend('bs', bs)
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
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("slate","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")

cmd.extend('bs', bs)


def bsbw(selection='all'):
    ''' 
    DESCRIPTION:
    bs creates a gray-scaled ball-and-stick representation of an object. 

    USAGE:
    bsbw

    ARGUMENTS:
    None
    EXAMPLE:
    bsbw

    MORE DETAILS:
    bs creates a gray-scaled ball-and-stick representation of an object. 


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def bsbw(selection='all'):
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
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("gray55","elem O and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem CL and "+selection)
    cmd.color("gray40","elem K and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")
cmd.extend('bsbw', bsbw)
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
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("gray55","elem O and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem CL and "+selection)
    cmd.color("gray40","elem K and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")
cmd.extend('bsbw', bsbw)


def bsbwsc(selection='all'):
    ''' 
    DESCRIPTION:
    bsbwsc creates a gray-scaled ball-and-stick representation of an object or selection.
Only the side chains are shown as ball and stick when used with a cartoon (ribbon diagram).

    USAGE:
    bsbwsc
    bsbwsc selectedResidues

    ARGUMENTS:
    None or a named selection
    EXAMPLE:
    bsbwsc

    MORE DETAILS:
    bsbwsc creates a gray-scaled ball-and-stick representation of an object or selection.
    Only the side chains are shown as ball and stick when used with a cartoon (ribbon diagram).


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def bsbwsc(selection='all'):
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
    cmd.show("cartoon", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("gray55","elem O and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem CL and "+selection)
    cmd.color("gray40","elem K and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.do("set cartoon_side_chain_helper, on")
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")
cmd.extend('bsbwsc', bsbwsc)
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
    cmd.show("cartoon", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)
    cmd.color("gray85","elem C and "+selection)
    cmd.color("gray98","elem H and "+selection)
    cmd.color("gray55","elem O and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem CL and "+selection)
    cmd.color("gray40","elem K and "+selection)
    cmd.color("gray70","elem S and "+selection)
    cmd.color("gray40","elem N and "+selection)
    cmd.set("stick_radius",0.07, selection)
    cmd.do("set cartoon_side_chain_helper, on")
    cmd.set("sphere_scale",0.18, selection)
    cmd.set("sphere_scale",0.13, selection+" and elem H")
    cmd.set("dash_gap",0.01, selection)
    cmd.set("dash_radius",0.07, selection)
    cmd.set("stick_color","black", selection)
    cmd.set("dash_gap",0.01)
    cmd.set("dash_radius",0.035)
    cmd.hide("nonbonded", selection)
    cmd.hide("lines", selection)
    cmd.zoom(selection)
    cmd.hide("labels")
cmd.extend('bsbwsc', bsbwsc)


def bsvdw(arg1="all"):
    ''' 
    DESCRIPTION:
    Transparent vdw surface over ball and stick representation by Bobby Patton at Colorato State University. 

    USAGE:
    bsvdw selection

    ARGUMENTS:
    selection
    EXAMPLE:
    bsvdw 1lw9

    MORE DETAILS:
    Transparent vdw surface over ball and stick representation by Bobby Patton at Colorato State University. 
    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def bsvdw(arg1="all"):
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

    # defines VDW Sphere settings
    cmd.copy(arg1+"_vdw", arg1)
    cmd.set("sphere_scale",1.0, arg1+"_vdw and elem H")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1+"_vdw")
    cmd.hide("nonbonded", arg1+"_vdw")
    cmd.hide("lines", arg1+"_vdw")
    cmd.hide("sticks", arg1+"_vdw")
    cmd.hide("cartoon", arg1+"_vdw")
    cmd.show("spheres", arg1+"_vdw")
    cmd.set("sphere_transparency", 0.7, arg1+"_vdw")
    print("Note that the selection of 'all' does not work when applying this function")
    print("to multiple models as are found in pdb1 files and when multiple chains")
    print("make up the biological unit.")
    print("The shortcut has to be applied separately to each model.")


cmd.extend('bsvdw', bsvdw)
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

    # defines VDW Sphere settings
    cmd.copy(arg1+"_vdw", arg1)
    cmd.set("sphere_scale",1.0, arg1+"_vdw and elem H")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1+"_vdw")
    cmd.hide("nonbonded", arg1+"_vdw")
    cmd.hide("lines", arg1+"_vdw")
    cmd.hide("sticks", arg1+"_vdw")
    cmd.hide("cartoon", arg1+"_vdw")
    cmd.show("spheres", arg1+"_vdw")
    cmd.set("sphere_transparency", 0.7, arg1+"_vdw")
    print("Note that the selection of 'all' does not work when applying this function")
    print("to multiple models as are found in pdb1 files and when multiple chains")
    print("make up the biological unit.")
    print("The shortcut has to be applied separately to each model.")


cmd.extend('bsvdw', bsvdw)


def buriedW(sele='all', cutoff=-1, state=1, quiet=1, _self=cmd):
    ''' 
    DESCRIPTION:
    Return a selection of buried waters. 

    USAGE:
    buriedW

    ARGUMENTS:
    sele = string: atom selection {default: all}
cutoff = float: threshold on what one considers an "exposed"
atom (in A**2) {default: surface_residue_cutoff}
state = int: object state {default: 1}

    EXAMPLE:
    buried 1lw9


    MORE DETAILS:
    Return a selection of buried waters. 

Source: https://pymolwiki.org/index.php/Find_buried_waters

Author:  Unknown

License: GNU Free Documentation License 1.2
License Link: http://www.gnu.org/licenses/fdl-1.2.html


    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def buriedW(sele='all', cutoff=-1, state=1, quiet=1, _self=cmd):
    cutoff, state, quiet = float(cutoff), int(state), int(quiet)

    tmpObj=_self.get_unused_name("__tmp")
    _self.create(tmpObj, sele, state, 1, zoom=0)

    _self.set("dot_solvent", 1, tmpObj);
    _self.get_area(tmpObj, state=1, load_b=1)

    if cutoff < 0:
        cutoff = _self.get("surface_residue_cutoff")
    _self.remove(tmpObj + " and not solvent")
    _self.remove(tmpObj + " and b > %s" % cutoff)

    exposed = set()
    _self.iterate(tmpObj, "exposed.add((chain,resv))", space=locals())

    selName = _self.get_unused_name("buried")
    _self.select(selName, "(%s) in %s" % (sele, tmpObj))
	
    cmd.show("spheres","buried01")
    
    # clean up
    _self.delete(tmpObj)

    if not quiet: print('Found %d buried water atoms' % (len(exposed)))

    return sorted(exposed)


cmd.extend('buriedW', buriedW)
    '''

    cutoff, state, quiet = float(cutoff), int(state), int(quiet)

    tmpObj=_self.get_unused_name("__tmp")
    _self.create(tmpObj, sele, state, 1, zoom=0)

    _self.set("dot_solvent", 1, tmpObj);
    _self.get_area(tmpObj, state=1, load_b=1)

    if cutoff < 0:
        cutoff = _self.get("surface_residue_cutoff")
    _self.remove(tmpObj + " and not solvent")
    _self.remove(tmpObj + " and b > %s" % cutoff)

    exposed = set()
    _self.iterate(tmpObj, "exposed.add((chain,resv))", space=locals())

    selName = _self.get_unused_name("buried")
    _self.select(selName, "(%s) in %s" % (sele, tmpObj))
	
    cmd.show("spheres","buried01")
    
    # clean up
    _self.delete(tmpObj)

    if not quiet: print('Found %d buried water atoms' % (len(exposed)))

    return sorted(exposed)


cmd.extend('buriedW', buriedW)


def cartoonbw(arg1='all'):
    ''' 
    DESCRIPTION:
    Grayscale by secondary structure.

    USAGE:
    cartoonbw

    ARGUMENTS:
    selection

    EXAMPLE:
    cartoonbw

    MORE DETAILS:
    Grayscale by secondary structure.
    Helices were converted from red using formula in gscale().
    Strands were converted from yellow using formula in gscale().
    Loops were converted from green using formula in gscale().
    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def cartoonbw(arg1='all'):
    cmd.show_as('cartoon')
    cmd.color('grey21', 'ss H')
    cmd.color('grey92', 'ss S')
    cmd.color('grey71', 'ss L+')
    cmd.set('cartoon_discrete_colors','on')


cmd.extend('cartoonbw', cartoonbw)
    '''

    cmd.show_as('cartoon')
    cmd.color('grey21', 'ss H')
    cmd.color('grey92', 'ss S')
    cmd.color('grey71', 'ss L+')
    cmd.set('cartoon_discrete_colors','on')


cmd.extend('cartoonbw', cartoonbw)


def ccp4mg():
    ''' 
    DESCRIPTION:
    Open ccp4mg from within PyMOL. 



    USAGE:
    ccp4mg

    ARGUMENTS:
    None
    EXAMPLE:
    ccp4mg

    MORE DETAILS:
    Open ccp4mg from within PyMOL. 
    Adjust url for your location.


    VERTICAL PML SCRIPT:
    NA

    HORIZONTAL PML SCRIPT:
    NA

    PYTHON CODE:
def ccp4mg():
    try:
        print("Opening the molecular graphics program ccp4mg.");
        subprocess.check_output(ccp4mgCommand)
        print("Success opening ccp4mg.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the ccp4mgOpen'. \n  Or use 'ccp4mgPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('ccp4mg',ccp4mg)
    '''

    try:
        print("Opening the molecular graphics program ccp4mg.");
        subprocess.check_output(ccp4mgCommand)
        print("Success opening ccp4mg.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the ccp4mgOpen'. \n  Or use 'ccp4mgPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('ccp4mg',ccp4mg)


def cellbasis(angles, edges):
    ''' 
    DESCRIPTION:
    For the unit cell with given angles and edge lengths, calculate the basis transformation (vectors) as a 4x4 numpy.array. Used by the function supercell.

    USAGE:
    NA

    ARGUMENTS:
    NA
    EXAMPLE:
    NA

    MORE DETAILS:
    For the unit cell with given angles and edge lengths calculate the basis transformation 
     (vectors) as a 4x4 numpy.array. Used by the function supercell.
    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def cellbasis(angles, edges):
    rad = [radians(i) for i in angles]
    basis = numpy.identity(4)
    basis[0][1] = cos(rad[2])
    basis[1][1] = sin(rad[2])
    basis[0][2] = cos(rad[1])
    basis[1][2] = (cos(rad[0]) - basis[0][1]*basis[0][2])/basis[1][1]
    basis[2][2] = sqrt(1 - basis[0][2]**2 - basis[1][2]**2)
    edges.append(1.0)
    return basis * edges # numpy.array multiplication!

cmd.extend('cellbasis',cellbasis)
    '''

    rad = [radians(i) for i in angles]
    basis = numpy.identity(4)
    basis[0][1] = cos(rad[2])
    basis[1][1] = sin(rad[2])
    basis[0][2] = cos(rad[1])
    basis[1][2] = (cos(rad[0]) - basis[0][1]*basis[0][2])/basis[1][1]
    basis[2][2] = sqrt(1 - basis[0][2]**2 - basis[1][2]**2)
    edges.append(1.0)
    return basis * edges # numpy.array multiplication!

cmd.extend('cellbasis',cellbasis)


def checkParams(needle, haystack, selName, het, firstOnly):
    ''' 
    DESCRIPTION:
    Checks user params for the findSeq function.

    USAGE:
    NA

    ARGUMENTS:
    NA
    EXAMPLE:
    NA

    MORE DETAILS:
    Checks user params for the findSeq function.

    Source: https://pymolwiki.org/index.php/Findseq
    Author: Jason Vertrees
    Modificaiton: rewrote furnction name in camelcase
    License:  GNU Free Documentation License 1.2
    License Link: 
    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def checkParams(needle, haystack, selName, het, firstOnly):
    """
    This is just a helper function for checking the user input
    """
    # check Needle
    if len(needle) == 0 or not cmd.is_string(needle):
        print("Error: Please provide a string 'needle' to search for.")
        print("Error: For help type 'help motifFinder'.")
    return False

    # check Haystack
    if len(haystack) == 0 or not cmd.is_string(haystack):
        print("Error: Please provide valid PyMOL object or selection name")
        print("Error: in which to search.")
        print("Error: For help type 'help motifFinder'.")
    return False

    # check het
    try:
        het = bool(int(het))
    except ValueError:
        print("Error: The 'het' parameter was not 0 or 1.")
    return False

    # check first Only
    try:
        firstOnly = bool(int(het))
    except ValueError:
        print("Error: The 'firstOnly' parameter was not 0 or 1.")
    return False

    # check selName
    if not cmd.is_string(selName):
        print("Error: selName was not a string.")
        return False
    return True

    '''

    """
    This is just a helper function for checking the user input
    """
    # check Needle
    if len(needle) == 0 or not cmd.is_string(needle):
        print("Error: Please provide a string 'needle' to search for.")
        print("Error: For help type 'help motifFinder'.")
    return False

    # check Haystack
    if len(haystack) == 0 or not cmd.is_string(haystack):
        print("Error: Please provide valid PyMOL object or selection name")
        print("Error: in which to search.")
        print("Error: For help type 'help motifFinder'.")
    return False

    # check het
    try:
        het = bool(int(het))
    except ValueError:
        print("Error: The 'het' parameter was not 0 or 1.")
    return False

    # check first Only
    try:
        firstOnly = bool(int(het))
    except ValueError:
        print("Error: The 'firstOnly' parameter was not 0 or 1.")
    return False

    # check selName
    if not cmd.is_string(selName):
        print("Error: selName was not a string.")
        return False
    return True



def chimera(fileName="test.pdb"):
    ''' 
    DESCRIPTION:
    Open Chimera from within PyMOL. 
 

    USAGE:
    chimera

    ARGUMENTS:
    None
    EXAMPLE:
    chimera

    MORE DETAILS:
    Open Chimera from within PyMOL. 
    Adjust url for your location.


    VERTICAL PML SCRIPT:
        subprocess.call(chimeraOpen)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(chimeraOpen);return
    PYTHON CODE:
def chimera(fileName="test.pdb"):
    try:
        print("Opening the molecular graphics program CHIMERA.");
        subprocess.check_output(chimeraOpen)
        print("Success opening CHIMERA.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'chimeraOpen'. \n  Or use 'chimeraPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('chimera',chimera)
    '''

    try:
        print("Opening the molecular graphics program CHIMERA.");
        subprocess.check_output(chimeraOpen)
        print("Success opening CHIMERA.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'chimeraOpen'. \n  Or use 'chimeraPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('chimera',chimera)


def chimeraWeb():
    ''' 
    DESCRIPTION:
    Open the webste of UCSF Chimera.


    USAGE:
    chimeraWeb

    ARGUMENTS:
    NA
    EXAMPLE:
    chimeraWeb

    MORE DETAILS:
    Open the webste of UCSF Chimera.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def chimeraWeb():
    url=chimeriaURL
    try:
        print("Opening the website of UCSF Chimera.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of UCSF Chimera.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('chimeraWeb',chimeraWeb)
    '''

    url=chimeriaURL
    try:
        print("Opening the website of UCSF Chimera.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the website of UCSF Chimera.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('chimeraWeb',chimeraWeb)


def cntccp4s():
    ''' 
    DESCRIPTION:
    Count number of *.ccp4 (electron density map) files in current directory.

    USAGE:
    cntccp4s

    ARGUMENTS:
    None
    EXAMPLE:
    cntccp4s

    MORE DETAILS:
    Count number of *.ccp4 (electron density map) files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of ccp4 electron density files in current directory.");
    print("Usage: cntccp4s");
    myPath = os.getcwd();
    ccp4Counter = len(glob.glob1(myPath,"*.pse"));
    print("Number of number of ccp4 electron density files in the current directory: ", ccp4Counter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of ccp4 electron density files in current directory");print("Usage: cntpses");myPath = os.getcwd();pseCounter = len(glob.glob1(myPath,"*.pdb"));print("Number of pdb files in the current directory: ", pseCounter);

    PYTHON CODE:
def cntccp4s():
    print("Count the number of ccp4 electron density files in current directory.");
    print("Usage: cntccp4s");
    myPath = os.getcwd();
    ccp4Counter = len(glob.glob1(myPath,"*.pse"));
    print("Number of number of ccp4 electron density files in the current directory: ", ccp4Counter);


cmd.extend('cntccp4s',cntccp4s)
    '''

    print("Count the number of ccp4 electron density files in current directory.");
    print("Usage: cntccp4s");
    myPath = os.getcwd();
    ccp4Counter = len(glob.glob1(myPath,"*.pse"));
    print("Number of number of ccp4 electron density files in the current directory: ", ccp4Counter);


cmd.extend('cntccp4s',cntccp4s)


def cntfiles():
    ''' 
    DESCRIPTION:
    Count number of files in current directory.

    USAGE:
    cntfiles

    ARGUMENTS:
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
    print("Count the files in the directory.")
    print("Usage: cntfiles.")
    # simple version for working with CWD
    print("Number of files in current working directory: ", len([name for name in os.listdir('.') if os.path.isfile(name)]))
    return
cmd.extend('cntfiles',cntfiles)
    '''

    print("Count the files in the directory.")
    print("Usage: cntfiles.")
    # simple version for working with CWD
    print("Number of files in current working directory: ", len([name for name in os.listdir('.') if os.path.isfile(name)]))
    return
cmd.extend('cntfiles',cntfiles)


def cntlogs():
    ''' 
    DESCRIPTION:
    Count number of *.log files in current directory.

    USAGE:
    cntlogs

    ARGUMENTS:
    None
    EXAMPLE:
    cntlogs

    MORE DETAILS:
    Count number of *.log files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of log files in current directory.");
    print("Usage: cntlogs");
    myPath = os.getcwd();
    pngCounter = len(glob.glob1(myPath,"*.log"));
    print("Number of number of log image files in the current directory: ", logCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of log image files in current directory.");print("Usage: cntlogs");myPath = os.getcwd();logCounter = len(glob.glob1(myPath,"*.log"));print("Number of number of log files in the current directory: ", logCounter);


    PYTHON CODE:
def cntlogs():
    print("Count the number of log image files in current directory.");
    print("Usage: cntlogs");
    myPath = os.getcwd();
    logCounter = len(glob.glob1(myPath,"*.log"));
    print("Number of number of log image files in the current directory: ", logCounter);


cmd.extend('cntlogs',cntlogs)
    '''

    print("Count the number of log image files in current directory.");
    print("Usage: cntlogs");
    myPath = os.getcwd();
    logCounter = len(glob.glob1(myPath,"*.log"));
    print("Number of number of log image files in the current directory: ", logCounter);


cmd.extend('cntlogs',cntlogs)


def cntmtzs():
    ''' 
    DESCRIPTION:
    Count number of *.mtz (structure factor) files in current directory.

    USAGE:
    cntmtzs

    ARGUMENTS:
    None
    EXAMPLE:
    cntmtzs

    MORE DETAILS:
    Count number of *.mtz (structure factor) files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of mtz structure factor files in current directory.");
    print("Usage: cntmtzs");
    myPath = os.getcwd();
    mtzCounter = len(glob.glob1(myPath,"*.mtz"));
    print("Number of number of mtz structure factor  files in the current directory: ", mtzCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of mtz structure factor files in current directory.");
print("Usage: cntmtzs");myPath = os.getcwd();mtzCounter = len(glob.glob1(myPath,"*.mtz"));print("Number of number of mtz structure factor  files in the current directory: ", mtzCounter);


    PYTHON CODE:
def cntmtzs():
    print("Count the number of mtz structure factor files in current directory.");
    print("Usage: cntmtzs");
    myPath = os.getcwd();
    mtzCounter = len(glob.glob1(myPath,"*.mtz"));
    print("Number of number of mtz structure factor  files in the current directory: ", mtzCounter);


cmd.extend('cntmtzs',cntmtzs)
    '''

    print("Count the number of mtz structure factor files in current directory.");
    print("Usage: cntmtzs");
    myPath = os.getcwd();
    mtzCounter = len(glob.glob1(myPath,"*.mtz"));
    print("Number of number of mtz structure factor  files in the current directory: ", mtzCounter);


cmd.extend('cntmtzs',cntmtzs)


def cntpdbs():
    ''' 
    DESCRIPTION:
    Count number of pdb files in current directory.

    USAGE:
    cntpdbs

    ARGUMENTS:
    None
    EXAMPLE:
    cntpdbs

    MORE DETAILS:
    Count number of pdb files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of pdb files in current working directory.");
    print("Usage: cntpdbs");
    myPath = os.getcwd();
    pdbCounter = len(glob.glob1(myPath,"*.pdb"));
    print("Number of pdb files in the current directory: ", pdbCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of pdb files in current directory.");print("Usage: cntpdb");myPath = os.getcwd();pdbCounter = len(glob.glob1(myPath,"*.pdb"));print("Number of pdb files in the current directory: ", pdbCounter);

    PYTHON CODE:
def cntpdbs():
    print("Count the number of pdb files in the current directory.")
    print("Usage: cntpdb")
    myPath = os.getcwd()
    pdbCounter = len(glob.glob1(myPath,"*.pdb"))
    print("Number of pdb files in the current directory: ", pdbCounter)
    return

cmd.extend('cntpdbs',cntpdbs)
    '''

    print("Count the number of pdb files in the current directory.")
    print("Usage: cntpdb")
    myPath = os.getcwd()
    pdbCounter = len(glob.glob1(myPath,"*.pdb"))
    print("Number of pdb files in the current directory: ", pdbCounter)
    return

cmd.extend('cntpdbs',cntpdbs)


def cntpmls():
    ''' 
    DESCRIPTION:
    Count number of pml (Pymol macro language) files in current directory.

    USAGE:
    cntpmls

    ARGUMENTS:
    None
    EXAMPLE:
    cntpmls

    MORE DETAILS:
    Count number of pml (Pymol macro language) files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of pml (PyMOL macro language) files in current directory.");
    print("Usage: cntpmls");
    myPath = os.getcwd();
    pmlCounter = len(glob.glob1(myPath,"*.pml"));
    print("Number of pml (PyMOL macro language)  files in the current directory: ", pmlCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of pml  (Pymol macro language) files in current directory.");print("Usage: cntpmls");myPath = os.getcwd();pmlCounter = len(glob.glob1(myPath,"*.pdb"));print("Number of pml  files in the current directory: ", pmlCounter);

    PYTHON CODE:
def cntpmls():
    print("Count the number of pml (Pymol macro language) files in current directory.");
    print("Usage: cntpmls");
    myPath = os.getcwd();
    pmlCounter = len(glob.glob1(myPath,"*.pml"));
    print("Number of pml files in the current directory: ", pmlCounter);
    return	

cmd.extend('cntpmls',cntpmls)
    '''

    print("Count the number of pml (Pymol macro language) files in current directory.");
    print("Usage: cntpmls");
    myPath = os.getcwd();
    pmlCounter = len(glob.glob1(myPath,"*.pml"));
    print("Number of pml files in the current directory: ", pmlCounter);
    return	

cmd.extend('cntpmls',cntpmls)


def cntpngs():
    ''' 
    DESCRIPTION:
    Count number of *.png image files in current directory.

    USAGE:
    cntpngs

    ARGUMENTS:
    None
    EXAMPLE:
    cntpngs

    MORE DETAILS:
    Count number of *.png image files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of png image files in current directory.");
    print("Usage: cntpngs");
    myPath = os.getcwd();
    pngCounter = len(glob.glob1(myPath,"*.png"));
    print("Number of number of png image files in the current directory: ", pngCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of png image files in current directory.");print("Usage: cntpngs");myPath = os.getcwd();pngCounter = len(glob.glob1(myPath,"*.png"));print("Number of number of png image files in the current directory: ", pngCounter);


    PYTHON CODE:
def cntpngs():
    print("Count the number of png image files in current directory.");
    print("Usage: cntpngs");
    myPath = os.getcwd();
    pngCounter = len(glob.glob1(myPath,"*.png"));
    print("Number of number of png image files in the current directory: ", pngCounter);


cmd.extend('cntpngs',cntpngs)
    '''

    print("Count the number of png image files in current directory.");
    print("Usage: cntpngs");
    myPath = os.getcwd();
    pngCounter = len(glob.glob1(myPath,"*.png"));
    print("Number of number of png image files in the current directory: ", pngCounter);


cmd.extend('cntpngs',cntpngs)


def cntpses():
    ''' 
    DESCRIPTION:
    Count number of *.pse (session) files in current directory.

    USAGE:
    cntpses

    ARGUMENTS:
    None
    EXAMPLE:
    cntpses

    MORE DETAILS:
    Count number of *.pse (session) files in current directory.


    VERTICAL PML SCRIPT:
    print("Count the number of pmls files in current directory.");
    print("Usage: cntpses");
    myPath = os.getcwd();
    pseCounter = len(glob.glob1(myPath,"*.pse"));
    print("Number of pml files in the current directory: ", pseCounter);


    HORIZONTAL PML SCRIPT:
    print("Count the number of *.pse files (session files) in current directory.");print("Usage: cntpses");myPath = os.getcwd();pseCounter = len(glob.glob1(myPath,"*.pdb"));print("Number of pdb files in the current directory: ", pseCounter);

    PYTHON CODE:
def cntpses():
    print("Count the number of *.pse (session) files in current directory.");
    print("Usage: cntpses");
    myPath = os.getcwd();
    pseCounter = len(glob.glob1(myPath,"*.pse"));
    print("Number of *.pse (session) files in the current directory: ", pseCounter);
    return	

cmd.extend('cntpses',cntpses)
    '''

    print("Count the number of *.pse (session) files in current directory.");
    print("Usage: cntpses");
    myPath = os.getcwd();
    pseCounter = len(glob.glob1(myPath,"*.pse"));
    print("Number of *.pse (session) files in the current directory: ", pseCounter);
    return	

cmd.extend('cntpses',cntpses)


def code():
    ''' 
    DESCRIPTION:
    Open file with Visual Studio Code from within PyMOL.

    USAGE:
    code

    ARGUMENTS:
    None
    EXAMPLE:
    code script.pml

    MORE DETAILS:
    Open file with Visual Studio Code from within PyMOL. 
    Install the bioSyntax extension (free) from the Visual Studio Marketplace
    to get color syntax highlighting of pml files along with Fasta and other
    sequence files and PDB coordinate files. 

    https://marketplace.visualstudio.com/items?itemName=reageyao.biosyntax



    VERTICAL PML SCRIPT:
    subprocess.call(codeOpen);
     return

    HORIZONTAL PML SCRIPT:
    subprocess.call(codeOpen);return

    PYTHON CODE:
def code():
    try:
        print("Opening the molecular graphics program Virtual Studio Code.");
        subprocess.check_output(codeOpen)
        print("Success opening Virtual Studio Code.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the codeOpen'. \n  Or use 'codePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('code',code)
    '''

    try:
        print("Opening the molecular graphics program Virtual Studio Code.");
        subprocess.check_output(codeOpen)
        print("Success opening Virtual Studio Code.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the codeOpen'. \n  Or use 'codePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('code',code)


def colorh1(selection='all'):
    ''' 
    DESCRIPTION:
    Color protein molecules according to the Eisenberg hydrophobicity scale. Uses scheme 1.

    USAGE:
    colorh1

    ARGUMENTS:
     <selection>
    EXAMPLE:
    colorh1
colorh1 1lw9

    MORE DETAILS:
    Color protein molecules according to the Eisenberg hydrophobicity scale. Uses scheme 1 according to https://pymolwiki.org/index.php/Color_h.


    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def colorh1(selection='all'):
    s = str(selection)
    print(s)
    cmd.set_color('color_ile',[0.996,0.062,0.062])
    cmd.set_color('color_phe',[0.996,0.109,0.109])
    cmd.set_color('color_val',[0.992,0.156,0.156])
    cmd.set_color('color_leu',[0.992,0.207,0.207])
    cmd.set_color('color_trp',[0.992,0.254,0.254])
    cmd.set_color('color_met',[0.988,0.301,0.301])
    cmd.set_color('color_ala',[0.988,0.348,0.348])
    cmd.set_color('color_gly',[0.984,0.394,0.394])
    cmd.set_color('color_cys',[0.984,0.445,0.445])
    cmd.set_color('color_tyr',[0.984,0.492,0.492])
    cmd.set_color('color_pro',[0.980,0.539,0.539])
    cmd.set_color('color_thr',[0.980,0.586,0.586])
    cmd.set_color('color_ser',[0.980,0.637,0.637])
    cmd.set_color('color_his',[0.977,0.684,0.684])
    cmd.set_color('color_glu',[0.977,0.730,0.730])
    cmd.set_color('color_asn',[0.973,0.777,0.777])
    cmd.set_color('color_gln',[0.973,0.824,0.824])
    cmd.set_color('color_asp',[0.973,0.875,0.875])
    cmd.set_color('color_lys',[0.899,0.922,0.922])
    cmd.set_color('color_arg',[0.899,0.969,0.969])
    cmd.color("color_ile","("+s+" and resn ile)")
    cmd.color("color_phe","("+s+" and resn phe)")
    cmd.color("color_val","("+s+" and resn val)")
    cmd.color("color_leu","("+s+" and resn leu)")
    cmd.color("color_trp","("+s+" and resn trp)")
    cmd.color("color_met","("+s+" and resn met)")
    cmd.color("color_ala","("+s+" and resn ala)")
    cmd.color("color_gly","("+s+" and resn gly)")
    cmd.color("color_cys","("+s+" and resn cys)")
    cmd.color("color_tyr","("+s+" and resn tyr)")
    cmd.color("color_pro","("+s+" and resn pro)")
    cmd.color("color_thr","("+s+" and resn thr)")
    cmd.color("color_ser","("+s+" and resn ser)")
    cmd.color("color_his","("+s+" and resn his)")
    cmd.color("color_glu","("+s+" and resn glu)")
    cmd.color("color_asn","("+s+" and resn asn)")
    cmd.color("color_gln","("+s+" and resn gln)")
    cmd.color("color_asp","("+s+" and resn asp)")
    cmd.color("color_lys","("+s+" and resn lys)")
    cmd.color("color_arg","("+s+" and resn arg)")

cmd.extend('colorh1',colorh1)
    '''

    s = str(selection)
    print(s)
    cmd.set_color('color_ile',[0.996,0.062,0.062])
    cmd.set_color('color_phe',[0.996,0.109,0.109])
    cmd.set_color('color_val',[0.992,0.156,0.156])
    cmd.set_color('color_leu',[0.992,0.207,0.207])
    cmd.set_color('color_trp',[0.992,0.254,0.254])
    cmd.set_color('color_met',[0.988,0.301,0.301])
    cmd.set_color('color_ala',[0.988,0.348,0.348])
    cmd.set_color('color_gly',[0.984,0.394,0.394])
    cmd.set_color('color_cys',[0.984,0.445,0.445])
    cmd.set_color('color_tyr',[0.984,0.492,0.492])
    cmd.set_color('color_pro',[0.980,0.539,0.539])
    cmd.set_color('color_thr',[0.980,0.586,0.586])
    cmd.set_color('color_ser',[0.980,0.637,0.637])
    cmd.set_color('color_his',[0.977,0.684,0.684])
    cmd.set_color('color_glu',[0.977,0.730,0.730])
    cmd.set_color('color_asn',[0.973,0.777,0.777])
    cmd.set_color('color_gln',[0.973,0.824,0.824])
    cmd.set_color('color_asp',[0.973,0.875,0.875])
    cmd.set_color('color_lys',[0.899,0.922,0.922])
    cmd.set_color('color_arg',[0.899,0.969,0.969])
    cmd.color("color_ile","("+s+" and resn ile)")
    cmd.color("color_phe","("+s+" and resn phe)")
    cmd.color("color_val","("+s+" and resn val)")
    cmd.color("color_leu","("+s+" and resn leu)")
    cmd.color("color_trp","("+s+" and resn trp)")
    cmd.color("color_met","("+s+" and resn met)")
    cmd.color("color_ala","("+s+" and resn ala)")
    cmd.color("color_gly","("+s+" and resn gly)")
    cmd.color("color_cys","("+s+" and resn cys)")
    cmd.color("color_tyr","("+s+" and resn tyr)")
    cmd.color("color_pro","("+s+" and resn pro)")
    cmd.color("color_thr","("+s+" and resn thr)")
    cmd.color("color_ser","("+s+" and resn ser)")
    cmd.color("color_his","("+s+" and resn his)")
    cmd.color("color_glu","("+s+" and resn glu)")
    cmd.color("color_asn","("+s+" and resn asn)")
    cmd.color("color_gln","("+s+" and resn gln)")
    cmd.color("color_asp","("+s+" and resn asp)")
    cmd.color("color_lys","("+s+" and resn lys)")
    cmd.color("color_arg","("+s+" and resn arg)")

cmd.extend('colorh1',colorh1)


def colorh2(selection='all'):
    ''' 
    DESCRIPTION:
    Color protein molecules according to the Eisenberg hydrophobicity scale. Uses scheme 2.

    USAGE:
    colorh2

    ARGUMENTS:
     <selection>
    EXAMPLE:
    colorh2
colorh2 1lw9

    MORE DETAILS:
    Color protein molecules according to the Eisenberg hydrophobicity scale. Uses scheme 2 according to  https://pymolwiki.org/index.php/Color_h.


    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def colorh2(selection='all'):
    s = str(selection)
    print(s)
    cmd.set_color("color_ile2",[0.938,1,0.938])
    cmd.set_color("color_phe2",[0.891,1,0.891])
    cmd.set_color("color_val2",[0.844,1,0.844])
    cmd.set_color("color_leu2",[0.793,1,0.793])
    cmd.set_color("color_trp2",[0.746,1,0.746])
    cmd.set_color("color_met2",[0.699,1,0.699])
    cmd.set_color("color_ala2",[0.652,1,0.652])
    cmd.set_color("color_gly2",[0.606,1,0.606])
    cmd.set_color("color_cys2",[0.555,1,0.555])
    cmd.set_color("color_tyr2",[0.508,1,0.508])
    cmd.set_color("color_pro2",[0.461,1,0.461])
    cmd.set_color("color_thr2",[0.414,1,0.414])
    cmd.set_color("color_ser2",[0.363,1,0.363])
    cmd.set_color("color_his2",[0.316,1,0.316])
    cmd.set_color("color_glu2",[0.27,1,0.27])
    cmd.set_color("color_asn2",[0.223,1,0.223])
    cmd.set_color("color_gln2",[0.176,1,0.176])
    cmd.set_color("color_asp2",[0.125,1,0.125])
    cmd.set_color("color_lys2",[0.078,1,0.078])
    cmd.set_color("color_arg2",[0.031,1,0.031])   
    cmd.color("color_ile2","("+s+" and resn ile)")
    cmd.color("color_phe2","("+s+" and resn phe)")
    cmd.color("color_val2","("+s+" and resn val)")
    cmd.color("color_leu2","("+s+" and resn leu)")
    cmd.color("color_trp2","("+s+" and resn trp)")
    cmd.color("color_met2","("+s+" and resn met)")
    cmd.color("color_ala2","("+s+" and resn ala)")
    cmd.color("color_gly2","("+s+" and resn gly)")
    cmd.color("color_cys2","("+s+" and resn cys)")
    cmd.color("color_tyr2","("+s+" and resn tyr)")
    cmd.color("color_pro2","("+s+" and resn pro)")
    cmd.color("color_thr2","("+s+" and resn thr)")
    cmd.color("color_ser2","("+s+" and resn ser)")
    cmd.color("color_his2","("+s+" and resn his)")
    cmd.color("color_glu2","("+s+" and resn glu)")
    cmd.color("color_asn2","("+s+" and resn asn)")
    cmd.color("color_gln2","("+s+" and resn gln)")
    cmd.color("color_asp2","("+s+" and resn asp)")
    cmd.color("color_lys2","("+s+" and resn lys)")
    cmd.color("color_arg2","("+s+" and resn arg)")

cmd.extend('colorh2',colorh2)
    '''

    s = str(selection)
    print(s)
    cmd.set_color("color_ile2",[0.938,1,0.938])
    cmd.set_color("color_phe2",[0.891,1,0.891])
    cmd.set_color("color_val2",[0.844,1,0.844])
    cmd.set_color("color_leu2",[0.793,1,0.793])
    cmd.set_color("color_trp2",[0.746,1,0.746])
    cmd.set_color("color_met2",[0.699,1,0.699])
    cmd.set_color("color_ala2",[0.652,1,0.652])
    cmd.set_color("color_gly2",[0.606,1,0.606])
    cmd.set_color("color_cys2",[0.555,1,0.555])
    cmd.set_color("color_tyr2",[0.508,1,0.508])
    cmd.set_color("color_pro2",[0.461,1,0.461])
    cmd.set_color("color_thr2",[0.414,1,0.414])
    cmd.set_color("color_ser2",[0.363,1,0.363])
    cmd.set_color("color_his2",[0.316,1,0.316])
    cmd.set_color("color_glu2",[0.27,1,0.27])
    cmd.set_color("color_asn2",[0.223,1,0.223])
    cmd.set_color("color_gln2",[0.176,1,0.176])
    cmd.set_color("color_asp2",[0.125,1,0.125])
    cmd.set_color("color_lys2",[0.078,1,0.078])
    cmd.set_color("color_arg2",[0.031,1,0.031])   
    cmd.color("color_ile2","("+s+" and resn ile)")
    cmd.color("color_phe2","("+s+" and resn phe)")
    cmd.color("color_val2","("+s+" and resn val)")
    cmd.color("color_leu2","("+s+" and resn leu)")
    cmd.color("color_trp2","("+s+" and resn trp)")
    cmd.color("color_met2","("+s+" and resn met)")
    cmd.color("color_ala2","("+s+" and resn ala)")
    cmd.color("color_gly2","("+s+" and resn gly)")
    cmd.color("color_cys2","("+s+" and resn cys)")
    cmd.color("color_tyr2","("+s+" and resn tyr)")
    cmd.color("color_pro2","("+s+" and resn pro)")
    cmd.color("color_thr2","("+s+" and resn thr)")
    cmd.color("color_ser2","("+s+" and resn ser)")
    cmd.color("color_his2","("+s+" and resn his)")
    cmd.color("color_glu2","("+s+" and resn glu)")
    cmd.color("color_asn2","("+s+" and resn asn)")
    cmd.color("color_gln2","("+s+" and resn gln)")
    cmd.color("color_asp2","("+s+" and resn asp)")
    cmd.color("color_lys2","("+s+" and resn lys)")
    cmd.color("color_arg2","("+s+" and resn arg)")

cmd.extend('colorh2',colorh2)


def coot(fileName="test.pdb"):
    ''' 
    DESCRIPTION:
    Open coot from within PyMOL. 



    USAGE:
    coot

    ARGUMENTS:
    None
    EXAMPLE:
    coot

    MORE DETAILS:
    Open coot from within PyMOL. 


    VERTICAL PML SCRIPT:
    arg = (cootPath + fileName)
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    arg = (cootPath  + fileName);subprocess.call(arg,shell=True);return

    PYTHON CODE:
def coot(fileName="test.pdb"):
    try:
        print("Opening the molecular graphics program COOT.");
        subprocess.check_output(cootOpen)
        print("Success opening COOT.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'cootOpen'. \n  Or use 'cootPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('coot',coot)
    '''

    try:
        print("Opening the molecular graphics program COOT.");
        subprocess.check_output(cootOpen)
        print("Success opening COOT.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'cootOpen'. \n  Or use 'cootPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('coot',coot)


def cranR():
    ''' 
    DESCRIPTION:
    Open Cran R from within PyMOL.

    USAGE:
    cranR

    ARGUMENTS:
    None
    EXAMPLE:
    ROpen   

    MORE DETAILS:
    Open Cran R from within PyMOL.



    VERTICAL PML SCRIPT:
    subprocess.call(ROpen); 
return      
    HORIZONTAL PML SCRIPT:
    subprocess.call(ROpen); return      

    PYTHON CODE:
def cranR():
    try:
        print("Opening the Cran R.");
        subprocess.check_output(ROpen)
        print("Success opening R.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'RlOpen'. \n  Or use 'RPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('cranR',cranR)
    '''

    try:
        print("Opening the Cran R.");
        subprocess.check_output(ROpen)
        print("Success opening R.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'RlOpen'. \n  Or use 'RPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('cranR',cranR)


def ddb():
    ''' 
    DESCRIPTION:
    Open DBBrowserSQLite. 

    USAGE:
    ddb

    ARGUMENTS:
    None
    EXAMPLE:
    ddb database.db

    MORE DETAILS:
    Open DBBrowserSQLite. 


    VERTICAL PML SCRIPT:
    arg = dbbrowserPath;
    subprocess.call(arg,shell=True);
    return

    HORIZONTAL PML SCRIPT:
    arg = dbbrowserPath;subprocess.call(arg,shell=True);return
    PYTHON CODE:
def ddb():
    try:
        print("Opening the DBBrowserSQLite.");
        subprocess.check_output(DBBrowserSQLiteOpen)
        print("Success opening DBBrowserSQLite.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'excelOpen'. \n  Or use 'excelPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('ddb',ddb)
    '''

    try:
        print("Opening the DBBrowserSQLite.");
        subprocess.check_output(DBBrowserSQLiteOpen)
        print("Success opening DBBrowserSQLite.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'excelOpen'. \n  Or use 'excelPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('ddb',ddb)


def emacs(fileName="test.pml"):
    ''' 
    DESCRIPTION:
    Open file with emacs from within PyMOL.

    USAGE:
    emacs

    ARGUMENTS:
    None
    EXAMPLE:
    emacs

    MORE DETAILS:
    Open file with emacs from within PyMOL. 
    Adjust path to emacs on your computer as needed.


    VERTICAL PML SCRIPT:
    subprocess.call(emacsOpen);
    return
    HORIZONTAL PML SCRIPT:
    subprocess.call(emacsOpen);return
    PYTHON CODE:
def emacs(fileName="test.pml"):
    try:
        print("Opening the text editor emacs.");
        subprocess.check_output(emacsOpen)
        print("Success opening emacs.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the emacsOpen'. \n  Or use 'emacsPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('emacs',emacs)
    '''

    try:
        print("Opening the text editor emacs.");
        subprocess.check_output(emacsOpen)
        print("Success opening emacs.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the emacsOpen'. \n  Or use 'emacsPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('emacs',emacs)


def excel():
    ''' 
    DESCRIPTION:
    Open excel from within PyMOL. 

    USAGE:
    excel

    ARGUMENTS:
    None
    EXAMPLE:
    excel

    MORE DETAILS:
    Open excel from within PyMOL. 


    VERTICAL PML SCRIPT:
    arg = excelCommand
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    arg = excelCommand;subprocess.call(arg,shell=True);return

    PYTHON CODE:
def excel():
    try:
        print("Opening the Microsoft Excel.");
        subprocess.check_output(excelOpen)
        print("Success opening Excel.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'excelOpen'. \n  Or use 'excelPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('excel',excel)
    '''

    try:
        print("Opening the Microsoft Excel.");
        subprocess.check_output(excelOpen)
        print("Success opening Excel.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'excelOpen'. \n  Or use 'excelPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('excel',excel)


def gcal():
    ''' 
    DESCRIPTION:
    Open Google Calendar.

    USAGE:
    gcal

    ARGUMENTS:
    None

    EXAMPLE:
    gcal

    MORE DETAILS:
    Open Google Calendar.


    VERTICAL PML SCRIPT:
    webbrowser.open('https://calendar.google.com/calendar/r')
    HORIZONTAL PML SCRIPT:
    webbrowser.open('https://calendar.google.com/calendar/r')
    PYTHON CODE:
def gcal():
    url=gcalURL
    try:
        print("Trying to open gcal.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening gcal.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('gcal',gcal)
    '''

    url=gcalURL
    try:
        print("Trying to open gcal.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Success opening gcal.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('gcal',gcal)


def gedit():
    ''' 
    DESCRIPTION:
    Open file with gedit from within PyMOL.

    USAGE:
    gedit

    ARGUMENTS:
    None

    EXAMPLE:
    gedit

    MORE DETAILS:
    Open file with gedit from within PyMOL. 
    Adjust the filepath for location of your executable.
    Can be installed via macports on Mac OS.


    VERTICAL PML SCRIPT:
    subprocess.call(emacsOpen)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(geditOpen);return

    PYTHON CODE:
def gedit():
    try:
        print("Opening the molecular graphics program gedit.");
        subprocess.check_output(geditOpen)
        print("Success opening gedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the geditOpen'. \n  Or use 'geditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('gedit',gedit)
    '''

    try:
        print("Opening the molecular graphics program gedit.");
        subprocess.check_output(geditOpen)
        print("Success opening gedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the geditOpen'. \n  Or use 'geditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('gedit',gedit)


def gimp():
    ''' 
    DESCRIPTION:
    Open the image editing program gimp from within PyMOL.


    USAGE:
    gimp

    ARGUMENTS:
    None
    EXAMPLE:
    gimp

    MORE DETAILS:
    Open the image editing program gimp from within PyMOL. 


    VERTICAL PML SCRIPT:
    arg = (gimpPath + fileName)
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    arg = (gimpPath + fileName);subprocess.call(arg,shell=True);return

    PYTHON CODE:
def gimp():
    try:
        print("Opening the gimp.");
        subprocess.check_output(gimpOpen)
        print("Success opening gimp.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'gimpOpen'. \n  Or use 'gimpPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('gimp',gimp)
    '''

    try:
        print("Opening the gimp.");
        subprocess.check_output(gimpOpen)
        print("Success opening gimp.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'gimpOpen'. \n  Or use 'gimpPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('gimp',gimp)


def gitAdd():
    ''' 
    DESCRIPTION:
    Enter help(gitAdd) to print steps for adding a file to version control.

    USAGE:
    gitAdd

    ARGUMENTS:
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
cmd.extend("gitAdd",gitAdd)
    '''

    #NA
cmd.extend("gitAdd",gitAdd)


def gitCommit():
    ''' 
    DESCRIPTION:
    Enter help(gitCommit) to print steps for saving updates to a file under version control.


    USAGE:
    gitCommit

    ARGUMENTS:
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
cmd.extend("gitCommit",gitCommit)
    '''

    #NA
cmd.extend("gitCommit",gitCommit)


def gitInit():
    ''' 
    DESCRIPTION:
    Enter help(gitInit) to print steps for creating a git repository.

    USAGE:
    gitInit

    ARGUMENTS:
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

cmd.extend("gitInit",gitInit)
    '''

    print(gitInit.__doc__)

cmd.extend("gitInit",gitInit)


def gitPull():
    ''' 
    DESCRIPTION:
    Enter help(gitPush) to print steps to supdate a repository on github.com.


    USAGE:
    gitPull


    ARGUMENTS:
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
cmd.extend("gitPull",gitPull)

    '''

    #NA
cmd.extend("gitPull",gitPull)



def gitPush():
    ''' 
    DESCRIPTION:
    Enter help(gitPush) to print steps update a repository on github.com. 


    USAGE:
    gitPush

    ARGUMENTS:
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
cmd.extend("gitPush",gitPush)
    '''

    #NA
cmd.extend("gitPush",gitPush)


def gmail():
    ''' 
    DESCRIPTION:
    Open gmail. 

    USAGE:
    gmail

    ARGUMENTS:
    None
    EXAMPLE:
    gmail

    MORE DETAILS:
    Open gmail. Edit url in python code below.


    VERTICAL PML SCRIPT:
    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')
    HORIZONTAL PML SCRIPT:
    webbrowser.open('https://mail.google.com/mail/u/0/#inbox')
    PYTHON CODE:
def gmail():
    url=gmailURL
    try:
        print("Trying to open gmail.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening gmail.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('gmail',gmail)
    '''

    url=gmailURL
    try:
        print("Trying to open gmail.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening gmail.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('gmail',gmail)


def gscale(selection="all"):
    ''' 
    DESCRIPTION:
    Apply grayscale to all atoms by element. 

    USAGE:
    gs <selection> # the default selection is 'all'


    ARGUMENTS:
    <selection>
    EXAMPLE:
    gs 3nd4


    MORE DETAILS:
    Apply grayscale to all atoms by element. 
    Made a dictionary of elements and their RGB values on the 0 to 1 scale.
    Converted the RGB values to grayscale with this formula:

    Y = 0.2126 * R + 0.7152 * G + 0.0722 * B

    where Y is the grayscale. 
    Wrote out the colors and elements as pml and python commands. 
    VERTICAL PML SCRIPT:
    color grey64, (elem Ac)
color grey67, (elem Al)
color grey39, (elem Am)
color grey46, (elem Sb)
color grey75, (elem Ar)
color grey58, (elem As)
color grey33, (elem At)
color grey56, (elem Ba)
color grey40, (elem Bk)
color grey87, (elem Be)
color grey40, (elem Bi)
color grey20, (elem Bh)
color grey77, (elem B)
color grey26, (elem Br)
color grey86, (elem Cd)
color grey76, (elem Ca)
color grey34, (elem Cf)
color grey77, (elem C)
color grey98, (elem Ce)
color grey17, (elem Cs)
color grey70, (elem Cl)
color grey60, (elem Cr)
color grey64, (elem Co)
color grey54, (elem Cu)
color grey42, (elem Cm)
color grey89, (elem D)
color grey19, (elem Db)
color grey79, (elem Dy)
color grey29, (elem Es)
color grey67, (elem Er)
color grey85, (elem Eu)
color grey28, (elem Fm)
color grey93, (elem F)
color grey8, (elem Fr)
color grey82, (elem Gd)
color grey60, (elem Ga)
color grey52, (elem Ge)
color grey80, (elem Au)
color grey68, (elem Hf)
color grey20, (elem Hs)
color grey96, (elem He)
color grey75, (elem Ho)
color grey89, (elem H)
color grey49, (elem In)
color grey16, (elem I)
color grey29, (elem Ir)
color grey48, (elem Fe)
color grey65, (elem Kr)
color grey76, (elem La)
color grey19, (elem Lr)
color grey34, (elem Pb)
color grey60, (elem Li)
color grey48, (elem Lu)
color grey83, (elem Mg)
color grey52, (elem Mn)
color grey20, (elem Mt)
color grey23, (elem Md)
color grey72, (elem Hg)
color grey62, (elem Mo)
color grey93, (elem Nd)
color grey85, (elem Ne)
color grey43, (elem Np)
color grey67, (elem Ni)
color grey69, (elem Nb)
color grey25, (elem N)
color grey23, (elem No)
color grey36, (elem Os)
color grey44, (elem O)
color grey33, (elem Pd)
color grey57, (elem P)
color grey82, (elem Pt)
color grey37, (elem Pu)
color grey40, (elem Po)
color grey35, (elem K)
color grey95, (elem Pr)
color grey90, (elem Pm)
color grey52, (elem Pa)
color grey35, (elem Ra)
color grey46, (elem Rn)
color grey43, (elem Re)
color grey39, (elem Rh)
color grey27, (elem Rb)
color grey47, (elem Ru)
color grey19, (elem Rf)
color grey89, (elem Sm)
color grey90, (elem Sc)
color grey20, (elem Sg)
color grey66, (elem Se)
color grey80, (elem Si)
color grey75, (elem Ag)
color grey46, (elem Na)
color grey71, (elem Sr)
color grey76, (elem S)
color grey60, (elem Ta)
color grey53, (elem Tc)
color grey51, (elem Te)
color grey81, (elem Tb)
color grey39, (elem Tl)
color grey59, (elem Th)
color grey61, (elem Tm)
color grey48, (elem Sn)
color grey75, (elem Ti)
color grey50, (elem W)
color grey47, (elem U)
color grey65, (elem V)
color grey54, (elem Xe)
color grey55, (elem Yb)
color grey91, (elem Y)
color grey51, (elem Zn)
color grey81, (elem Zr)


    HORIZONTAL PML SCRIPT:
    color grey64, (elem Ac);color grey67, (elem Al);color grey39, (elem Am);color grey46, (elem Sb);color grey58, (elem As);color grey33, (elem At);color grey56, (elem Ba);color grey40, (elem Bk);color grey87, (elem Be);color grey40, (elem Bi);color grey20, (elem Bh);color grey77, (elem B);color grey26, (elem Br);color grey86, (elem Cd);color grey76, (elem Ca);color grey34, (elem Cf);color grey77, (elem C);color grey98, (elem Ce);color grey70, (elem Cl);color grey60, (elem Cr);color grey64, (elem Co);color grey54, (elem Cu);color grey42, (elem Cm);color grey89, (elem D);color grey19, (elem Db);color grey79, (elem Dy);color grey29, (elem Es);color grey67, (elem Er);color grey85, (elem Eu);color grey28, (elem Fm);color grey93, (elem F);color grey8, (elem Fr);color grey82, (elem Gd);color grey60, (elem Ga);color grey52, (elem Ge);color grey80, (elem Au);color grey68, (elem Hf);color grey20, (elem Hs);color grey96, (elem He);color grey75, (elem Ho);color grey89, (elem H);color grey49, (elem In);color grey16, (elem I);color grey29, (elem Ir);color grey48, (elem Fe);color grey65, (elem Kr);color grey76, (elem La);color grey19, (elem Lr);color grey34, (elem Pb);color grey60, (elem Li);color grey48, (elem Lu);color grey83, (elem Mg);color grey52, (elem Mn);color grey20, (elem Mt);color grey23, (elem Md);color grey72, (elem Hg);color grey62, (elem Mo);color grey93, (elem Nd);color grey85, (elem Ne);color grey43, (elem Np);color grey67, (elem Ni);color grey69, (elem Nb);color grey25, (elem N);color grey23, (elem No);color grey36, (elem Os);color grey44, (elem O);color grey33, (elem Pd);color grey57, (elem P);color grey82, (elem Pt);color grey37, (elem Pu);color grey40, (elem Po);color grey35, (elem K);color grey95, (elem Pr);color grey90, (elem Pm);color grey52, (elem Pa);color grey35, (elem Ra);color grey46, (elem Rn);color grey43, (elem Re);color grey39, (elem Rh);color grey27, (elem Rb);color grey47, (elem Ru);color grey19, (elem Rf);color grey89, (elem Sm);color grey90, (elem Sc);color grey20, (elem Sg);color grey66, (elem Se);color grey80, (elem Si);color grey75, (elem Ag);color grey46, (elem Na);color grey71, (elem Sr);color grey76, (elem S);color grey60, (elem Ta);color grey53, (elem Tc);color grey51, (elem Te);color grey81, (elem Tb);color grey39, (elem Tl);color grey59, (elem Th);color grey61, (elem Tm);color grey48, (elem Sn);color grey75, (elem Ti);color grey50, (elem W);color grey47, (elem U);color grey65, (elem V);color grey54, (elem Xe);color grey55, (elem Yb);color grey91, (elem Y);color grey51, (elem Zn);color grey81, (elem Zr);

    PYTHON CODE:
def gscale(selection="all"):
    cmd.color('grey64', 'elem Ac')
    cmd.color('grey67', 'elem Al')
    cmd.color('grey39', 'elem Am')
    cmd.color('grey46', 'elem Sb')
    cmd.color('grey75', 'elem Ar')
    cmd.color('grey58', 'elem As')
    cmd.color('grey33', 'elem At')
    cmd.color('grey56', 'elem Ba')
    cmd.color('grey40', 'elem Bk')
    cmd.color('grey87', 'elem Be')
    cmd.color('grey40', 'elem Bi')
    cmd.color('grey20', 'elem Bh')
    cmd.color('grey77', 'elem B')
    cmd.color('grey26', 'elem Br')
    cmd.color('grey86', 'elem Cd')
    cmd.color('grey76', 'elem Ca')
    cmd.color('grey34', 'elem Cf')
    cmd.color('grey77', 'elem C')
    cmd.color('grey98', 'elem Ce')
    cmd.color('grey17', 'elem Cs')
    cmd.color('grey70', 'elem Cl')
    cmd.color('grey60', 'elem Cr')
    cmd.color('grey64', 'elem Co')
    cmd.color('grey54', 'elem Cu')
    cmd.color('grey42', 'elem Cm')
    cmd.color('grey89', 'elem D')
    cmd.color('grey19', 'elem Db')
    cmd.color('grey79', 'elem Dy')
    cmd.color('grey29', 'elem Es')
    cmd.color('grey67', 'elem Er')
    cmd.color('grey85', 'elem Eu')
    cmd.color('grey28', 'elem Fm')
    cmd.color('grey93', 'elem F')
    cmd.color('grey8', 'elem Fr')
    cmd.color('grey82', 'elem Gd')
    cmd.color('grey60', 'elem Ga')
    cmd.color('grey52', 'elem Ge')
    cmd.color('grey80', 'elem Au')
    cmd.color('grey68', 'elem Hf')
    cmd.color('grey20', 'elem Hs')
    cmd.color('grey96', 'elem He')
    cmd.color('grey75', 'elem Ho')
    cmd.color('grey89', 'elem H')
    cmd.color('grey49', 'elem In')
    cmd.color('grey16', 'elem I')
    cmd.color('grey29', 'elem Ir')
    cmd.color('grey48', 'elem Fe')
    cmd.color('grey65', 'elem Kr')
    cmd.color('grey76', 'elem La')
    cmd.color('grey19', 'elem Lr')
    cmd.color('grey34', 'elem Pb')
    cmd.color('grey60', 'elem Li')
    cmd.color('grey48', 'elem Lu')
    cmd.color('grey83', 'elem Mg')
    cmd.color('grey52', 'elem Mn')
    cmd.color('grey20', 'elem Mt')
    cmd.color('grey23', 'elem Md')
    cmd.color('grey72', 'elem Hg')
    cmd.color('grey62', 'elem Mo')
    cmd.color('grey93', 'elem Nd')
    cmd.color('grey85', 'elem Ne')
    cmd.color('grey43', 'elem Np')
    cmd.color('grey67', 'elem Ni')
    cmd.color('grey69', 'elem Nb')
    cmd.color('grey25', 'elem N')
    cmd.color('grey23', 'elem No')
    cmd.color('grey36', 'elem Os')
    cmd.color('grey44', 'elem O')
    cmd.color('grey33', 'elem Pd')
    cmd.color('grey57', 'elem P')
    cmd.color('grey82', 'elem Pt')
    cmd.color('grey37', 'elem Pu')
    cmd.color('grey40', 'elem Po')
    cmd.color('grey35', 'elem K')
    cmd.color('grey95', 'elem Pr')
    cmd.color('grey90', 'elem Pm')
    cmd.color('grey52', 'elem Pa')
    cmd.color('grey35', 'elem Ra')
    cmd.color('grey46', 'elem Rn')
    cmd.color('grey43', 'elem Re')
    cmd.color('grey39', 'elem Rh')
    cmd.color('grey27', 'elem Rb')
    cmd.color('grey47', 'elem Ru')
    cmd.color('grey19', 'elem Rf')
    cmd.color('grey89', 'elem Sm')
    cmd.color('grey90', 'elem Sc')
    cmd.color('grey20', 'elem Sg')
    cmd.color('grey66', 'elem Se')
    cmd.color('grey80', 'elem Si')
    cmd.color('grey75', 'elem Ag')
    cmd.color('grey46', 'elem Na')
    cmd.color('grey71', 'elem Sr')
    cmd.color('grey76', 'elem S')
    cmd.color('grey60', 'elem Ta')
    cmd.color('grey53', 'elem Tc')
    cmd.color('grey51', 'elem Te')
    cmd.color('grey81', 'elem Tb')
    cmd.color('grey39', 'elem Tl')
    cmd.color('grey59', 'elem Th')
    cmd.color('grey61', 'elem Tm')
    cmd.color('grey48', 'elem Sn')
    cmd.color('grey75', 'elem Ti')
    cmd.color('grey50', 'elem W')
    cmd.color('grey47', 'elem U')
    cmd.color('grey65', 'elem V')
    cmd.color('grey54', 'elem Xe')
    cmd.color('grey55', 'elem Yb')
    cmd.color('grey91', 'elem Y')
    cmd.color('grey51', 'elem Zn')
    cmd.color('grey81', 'elem Zr')
cmd.extend('gscale',gscale)
    '''

    cmd.color('grey64', 'elem Ac')
    cmd.color('grey67', 'elem Al')
    cmd.color('grey39', 'elem Am')
    cmd.color('grey46', 'elem Sb')
    cmd.color('grey75', 'elem Ar')
    cmd.color('grey58', 'elem As')
    cmd.color('grey33', 'elem At')
    cmd.color('grey56', 'elem Ba')
    cmd.color('grey40', 'elem Bk')
    cmd.color('grey87', 'elem Be')
    cmd.color('grey40', 'elem Bi')
    cmd.color('grey20', 'elem Bh')
    cmd.color('grey77', 'elem B')
    cmd.color('grey26', 'elem Br')
    cmd.color('grey86', 'elem Cd')
    cmd.color('grey76', 'elem Ca')
    cmd.color('grey34', 'elem Cf')
    cmd.color('grey77', 'elem C')
    cmd.color('grey98', 'elem Ce')
    cmd.color('grey17', 'elem Cs')
    cmd.color('grey70', 'elem Cl')
    cmd.color('grey60', 'elem Cr')
    cmd.color('grey64', 'elem Co')
    cmd.color('grey54', 'elem Cu')
    cmd.color('grey42', 'elem Cm')
    cmd.color('grey89', 'elem D')
    cmd.color('grey19', 'elem Db')
    cmd.color('grey79', 'elem Dy')
    cmd.color('grey29', 'elem Es')
    cmd.color('grey67', 'elem Er')
    cmd.color('grey85', 'elem Eu')
    cmd.color('grey28', 'elem Fm')
    cmd.color('grey93', 'elem F')
    cmd.color('grey8', 'elem Fr')
    cmd.color('grey82', 'elem Gd')
    cmd.color('grey60', 'elem Ga')
    cmd.color('grey52', 'elem Ge')
    cmd.color('grey80', 'elem Au')
    cmd.color('grey68', 'elem Hf')
    cmd.color('grey20', 'elem Hs')
    cmd.color('grey96', 'elem He')
    cmd.color('grey75', 'elem Ho')
    cmd.color('grey89', 'elem H')
    cmd.color('grey49', 'elem In')
    cmd.color('grey16', 'elem I')
    cmd.color('grey29', 'elem Ir')
    cmd.color('grey48', 'elem Fe')
    cmd.color('grey65', 'elem Kr')
    cmd.color('grey76', 'elem La')
    cmd.color('grey19', 'elem Lr')
    cmd.color('grey34', 'elem Pb')
    cmd.color('grey60', 'elem Li')
    cmd.color('grey48', 'elem Lu')
    cmd.color('grey83', 'elem Mg')
    cmd.color('grey52', 'elem Mn')
    cmd.color('grey20', 'elem Mt')
    cmd.color('grey23', 'elem Md')
    cmd.color('grey72', 'elem Hg')
    cmd.color('grey62', 'elem Mo')
    cmd.color('grey93', 'elem Nd')
    cmd.color('grey85', 'elem Ne')
    cmd.color('grey43', 'elem Np')
    cmd.color('grey67', 'elem Ni')
    cmd.color('grey69', 'elem Nb')
    cmd.color('grey25', 'elem N')
    cmd.color('grey23', 'elem No')
    cmd.color('grey36', 'elem Os')
    cmd.color('grey44', 'elem O')
    cmd.color('grey33', 'elem Pd')
    cmd.color('grey57', 'elem P')
    cmd.color('grey82', 'elem Pt')
    cmd.color('grey37', 'elem Pu')
    cmd.color('grey40', 'elem Po')
    cmd.color('grey35', 'elem K')
    cmd.color('grey95', 'elem Pr')
    cmd.color('grey90', 'elem Pm')
    cmd.color('grey52', 'elem Pa')
    cmd.color('grey35', 'elem Ra')
    cmd.color('grey46', 'elem Rn')
    cmd.color('grey43', 'elem Re')
    cmd.color('grey39', 'elem Rh')
    cmd.color('grey27', 'elem Rb')
    cmd.color('grey47', 'elem Ru')
    cmd.color('grey19', 'elem Rf')
    cmd.color('grey89', 'elem Sm')
    cmd.color('grey90', 'elem Sc')
    cmd.color('grey20', 'elem Sg')
    cmd.color('grey66', 'elem Se')
    cmd.color('grey80', 'elem Si')
    cmd.color('grey75', 'elem Ag')
    cmd.color('grey46', 'elem Na')
    cmd.color('grey71', 'elem Sr')
    cmd.color('grey76', 'elem S')
    cmd.color('grey60', 'elem Ta')
    cmd.color('grey53', 'elem Tc')
    cmd.color('grey51', 'elem Te')
    cmd.color('grey81', 'elem Tb')
    cmd.color('grey39', 'elem Tl')
    cmd.color('grey59', 'elem Th')
    cmd.color('grey61', 'elem Tm')
    cmd.color('grey48', 'elem Sn')
    cmd.color('grey75', 'elem Ti')
    cmd.color('grey50', 'elem W')
    cmd.color('grey47', 'elem U')
    cmd.color('grey65', 'elem V')
    cmd.color('grey54', 'elem Xe')
    cmd.color('grey55', 'elem Yb')
    cmd.color('grey91', 'elem Y')
    cmd.color('grey51', 'elem Zn')
    cmd.color('grey81', 'elem Zr')
cmd.extend('gscale',gscale)


def hb(selection='all'):
    ''' 
    DESCRIPTION:
    Creates an object of all H-bonds found by PyMOL.


    USAGE:
    hb <selection>

    ARGUMENTS:
    The selection is optional. It is "all" by default.


    EXAMPLE:
    hbond 1lw9

    MORE DETAILS:
    Creates an object of all polar contacts found by PyMOL. 
    They all may not be H-bonds (e.g., a water cannot partake in 
    five H-bonds at once). User be aware!


    VERTICAL PML SCRIPT:
    hbonds, all, all, 3.2, mode=2 
    HORIZONTAL PML SCRIPT:
    hbonds, all, all, 3.2, mode=2 
    PYTHON CODE:
def hb(selection='all'):
    cmd.distance("hbonds", "all", "all", "3.2", mode="2")
    cmd.set("dash_gap","0.4")
    cmd.set("dash_color","grey30")
    cmd.set("dash_width","1.5")
    cmd.set("dash_length",".25")
    print("Enter 'rmhb' to remove the hbonds.")
cmd.extend('hb',hb)
    '''

    cmd.distance("hbonds", "all", "all", "3.2", mode="2")
    cmd.set("dash_gap","0.4")
    cmd.set("dash_color","grey30")
    cmd.set("dash_width","1.5")
    cmd.set("dash_length",".25")
    print("Enter 'rmhb' to remove the hbonds.")
cmd.extend('hb',hb)


def inkscape():
    ''' 
    DESCRIPTION:
    Open the image editing program inkscape from within PyMOL. 

    USAGE:
    inkscape


    ARGUMENTS:
    None

    EXAMPLE:
    inkscape


    MORE DETAILS:
    Open the image editing program inkscape from within PyMOL. 


    VERTICAL PML SCRIPT:
    arg = (inkscapePath + fileName)
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    arg = (inkscapePath + fileName);subprocess.call(arg,shell=True);return

    PYTHON CODE:
def inkscape():
    try:
        print("Opening the inkscape.");
        subprocess.check_output(inkscapeOpen)
        print("Success opening inkscape.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'inkscapeOpen'. \n  Or use 'inkscapePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('inkscape',inkscape)
    '''

    try:
        print("Opening the inkscape.");
        subprocess.check_output(inkscapeOpen)
        print("Success opening inkscape.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'inkscapeOpen'. \n  Or use 'inkscapePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('inkscape',inkscape)


def interface(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    ''' 
    DESCRIPTION:
    Returns a selection of interface residues named according to what you passed into selName.

    USAGE:
    objname, selection1, selection2, cutoff=1.0, selectionNameForInterface


    ARGUMENTS:
    PARAMS
    cmpx
        The complex containing cA and cB

    cA
        The first chain in which we search for residues at an interface
        with cB

    cB
        The second chain in which we search for residues at an interface
        with cA

    cutoff
        The difference in area OVER which residues are considered
        interface residues.  Residues whose dASA from the complex to
        a single chain is greater than this cutoff are kept.  Zero
        keeps all residues.

    selName
        The name of the selection to return.

    EXAMPLE:
    interface 1BRS, c. C, c. F

    MORE DETAILS:
    interfaceResidues -- finds 'interface' residues between two chains in a complex.

    PARAMS
	cmpx
		The complex containing cA and cB

	cA
		The first chain in which we search for residues at an interface
		with cB

	cB
		The second chain in which we search for residues at an interface
		with cA

	cutoff
		The difference in area OVER which residues are considered
		interface residues.  Residues whose dASA from the complex to
		a single chain is greater than this cutoff are kept.  Zero
		keeps all residues.

	selName
		The name of the selection to return.

    RETURNS
	* A selection of interface residues is created and named
		depending on what you passed into selName
	* An array of values is returned where each value is:
		( modelName, residueNumber, dASA )

    NOTES
	If you have two chains that are not from the same PDB that you want
	to complex together, use the create command like:
		create myComplex, pdb1WithChainA or pdb2withChainX
	then pass myComplex to this script like:
		interfaceResidues myComlpex, c. A, c. X

	This script calculates the area of the complex as a whole.  Then,
	it separates the two chains that you pass in through the arguments
	cA and cB, alone.  Once it has this, it calculates the difference
	and any residues ABOVE the cutoff are called interface residues.

    AUTHOR:
	Jason Vertrees, 2009.
     
    SOURCE:
                Source: https://pymolwiki.org/index.php/InterfaceResidues

    LICENSE:
                GNU Free Documentation License 1.2
                https://www.gnu.org/licenses/old-licenses/fdl-1.2.en.html

			

    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def interface(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)

    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"

    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)

    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))

    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')

    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)

    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )

    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

    cmd.enable(cmpx)
    cmd.select(selName1, None)
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)

    # reset users settings
    cmd.set("dot_solvent", oldDS)

    return rVal

cmd.extend('interface', interface)
    '''

    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)

    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"

    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)

    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))

    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')

    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)

    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )

    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

    cmd.enable(cmpx)
    cmd.select(selName1, None)
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)

    # reset users settings
    cmd.set("dot_solvent", oldDS)

    return rVal

cmd.extend('interface', interface)


def iterm():
    ''' 
    DESCRIPTION:
    Open iTerm2 window on MacOS. 

    USAGE:
    iterm

    ARGUMENTS:
    None
    EXAMPLE:
    iterm

    MORE DETAILS:
    Open iTerm2 window on MacOS. 



    VERTICAL PML SCRIPT:
    subprocess.call(itermOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(itermOpen);return

    PYTHON CODE:
def iterm():
    try:
        print("Opening an iTerm window.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.Popen(itermOpen)
        print("Opened an iTerm window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('iterm',iterm)
    '''

    try:
        print("Opening an iTerm window.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.Popen(itermOpen)
        print("Opened an iTerm window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('iterm',iterm)


def jabref():
    ''' 
    DESCRIPTION:
    Open the jabref from within PyMOL.

    USAGE:
    jabref

    ARGUMENTS:
    None
    EXAMPLE:
    jabref

    MORE DETAILS:
    Open the jabref from within PyMOL.


    VERTICAL PML SCRIPT:
        subprocess.call(jabrefOpen);
    return

    HORIZONTAL PML SCRIPT:
        subprocess.call(jabrefOpen);return
    PYTHON CODE:
def jabref():
    try:
        print("Opening the bibliography manager JabRef.");
        subprocess.check_output(JabRefOpen)
        print("Success opening JabRef.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jabrefOpen'. \n  Or use 'jabrefPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jabref',jabref)
    '''

    try:
        print("Opening the bibliography manager JabRef.");
        subprocess.check_output(JabRefOpen)
        print("Success opening JabRef.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jabrefOpen'. \n  Or use 'jabrefPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jabref',jabref)


def jedit(fileName="test.pml"):
    ''' 
    DESCRIPTION:
    Open file with jedit from within PyMOL. 

    USAGE:
    jedit

    ARGUMENTS:
    None
    EXAMPLE:
    jedit script.pml

    MORE DETAILS:
    Open file with jedit from within PyMOL. 
    Adjust file path to your executable.
    Can be installed via macports on the Mac.


    VERTICAL PML SCRIPT:
    subprocess.call(jeditOpen)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(jeditOpen);return

    PYTHON CODE:
def jedit(fileName="test.pml"):
    try:
        print("Opening the molecular graphics program jedit.");
        subprocess.check_output(jeditOpen)
        print("Success opening jedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the jeditOpen'. \n  Or use 'jeditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jedit',jedit)
    '''

    try:
        print("Opening the molecular graphics program jedit.");
        subprocess.check_output(jeditOpen)
        print("Success opening jedit.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the jeditOpen'. \n  Or use 'jeditPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jedit',jedit)


def jmol():
    ''' 
    DESCRIPTION:
    Open Jmol from within PyMOL.


    USAGE:
    jmol

    ARGUMENTS:
    None

    EXAMPLE:
    jmol


    MORE DETAILS:
    Open Jmol from within PyMOL.
    Adjust file path to your location of Jmol.


    VERTICAL PML SCRIPT:
    arg = (jmolPath + fileName)
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    arg = (jmolPath + fileName);subprocess.call(arg,shell=True);return


    PYTHON CODE:
def jmol():
    try:
        print("Opening the molecular graphics program JMOL.");
        subprocess.check_output(jmolOpen)
        print("Success opening JMOL.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jmolOpen'. \n  Or use 'jmolPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jmol',jmol)
    '''

    try:
        print("Opening the molecular graphics program JMOL.");
        subprocess.check_output(jmolOpen)
        print("Success opening JMOL.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'jmolOpen'. \n  Or use 'jmolPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('jmol',jmol)


def julia():
    ''' 
    DESCRIPTION:
    Open the julia from within PyMOL.

    USAGE:
    julia

    ARGUMENTS:
    None

    EXAMPLE:
    julia

    MORE DETAILS:
    Open the julia from within PyMOL.


    VERTICAL PML SCRIPT:
    arg = juliaPath;
    subprocess.call(arg,shell=True);
    return

    HORIZONTAL PML SCRIPT:
    arg =juliaPath;subprocess.call(arg,shell=True);return

    PYTHON CODE:
def julia():
    try:
        print("Opening the REPL of the programming language julia REPL.");
        subprocess.check_output(juliaOpen)
        print("Success opening julia.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'juliaOpen'. \n  Or use 'juliaPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('julia',julia)
    '''

    try:
        print("Opening the REPL of the programming language julia REPL.");
        subprocess.check_output(juliaOpen)
        print("Success opening julia.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'juliaOpen'. \n  Or use 'juliaPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('julia',julia)


def juliapro():
    ''' 
    DESCRIPTION:
    Open the juliapro from within PyMOL.


    USAGE:
    juliapro

    ARGUMENTS:
    None

    EXAMPLE:
    juliapro

    MORE DETAILS:
    Open the juliapro from within PyMOL.


    VERTICAL PML SCRIPT:
    arg = juliaproPath;
    subprocess.call(arg,shell=True);
    return

    HORIZONTAL PML SCRIPT:
    arg =juliaproPath;subprocess.call(arg,shell=True);return

    PYTHON CODE:
def juliapro():
    try:
        print("Please be patient. Juliapro depends on atom which starts slowly.");
        subprocess.check_output(juliaproOpen)
        print("Success opening juliapro")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'juliaproOpen'. \n  Or use 'juliaproPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('juliapro',juliapro)

    '''

    try:
        print("Please be patient. Juliapro depends on atom which starts slowly.");
        subprocess.check_output(juliaproOpen)
        print("Success opening juliapro")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'juliaproOpen'. \n  Or use 'juliaproPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('juliapro',juliapro)



def mate():
    ''' 
    DESCRIPTION:
    Open textmate from within PyMOL. 



    USAGE:
    mate

    ARGUMENTS:
    None
    EXAMPLE:
    mate

    MORE DETAILS:
    Open file with Textmate (Mac OS only) from within PyMOL.  
    Adjust path to Textmate on your computer as needed.


    VERTICAL PML SCRIPT:
    subprocess.call(mateOpen)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(mateOpen);return

    PYTHON CODE:
def mate():
    try:
        print("Opening the molecular graphics program mate.");
        subprocess.check_output(textMateOpen)
        print("Success opening mate.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the mateOpen'. \n  Or use 'matePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('mate',mate)


    '''

    try:
        print("Opening the molecular graphics program mate.");
        subprocess.check_output(textMateOpen)
        print("Success opening mate.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the mateOpen'. \n  Or use 'matePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('mate',mate)




def nmr():

    ''' 
    DESCRIPTION:
    Show all models in a nmr structure. 

    USAGE:
    nmr

    ARGUMENTS:
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
cmd.extend("nmr", nmr)

    '''

    cmd.do('set all_states, on')
cmd.extend("nmr", nmr)



def nmroff():

    ''' 
    DESCRIPTION:
    Hide all but first model in a nmr structure.


    USAGE:
    nmroff


    ARGUMENTS:
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

cmd.extend("nmr", nmr)

    '''

    cmd.do('set all_states, off')

cmd.extend("nmr", nmr)



def npp():
    ''' 
    DESCRIPTION:
    Open notepadpp from within PyMOL. 

    USAGE:
    npp

    ARGUMENTS:
    None
    EXAMPLE:
    npp

    MORE DETAILS:
    Open notepadpp from within PyMOL. 
    Notepadpp can be installed on Mac OS and Linux via wine.


    VERTICAL PML SCRIPT:
    subprocess.call(nppOpen)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(nppOpen);return

    PYTHON CODE:
def npp():
    try:
        print("Opening the molecular graphics program notepad++.");
        subprocess.check_output(nppOpen)
        print("Success opening notepad++.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the nppOpen'. \n  Or use 'nppPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('npp',npp)
    '''

    try:
        print("Opening the molecular graphics program notepad++.");
        subprocess.check_output(nppOpen)
        print("Success opening notepad++.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the nppOpen'. \n  Or use 'nppPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('npp',npp)


def nv():
    ''' 
    DESCRIPTION:
    Open neovim from within PyMOL.

    USAGE:
    nv

    ARGUMENTS:
    None
    EXAMPLE:
    nv

    MORE DETAILS:
    Open neovim from within PyMOL. 
    Adjust path to neovim on your computer as needed.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA


    PYTHON CODE:
def nv():
    try:
        print("Opening the text editor neovim.");
        subprocess.check_output(neovimOpen)
        print("Success opening neovim.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the neovimOpen'. \n  Or use 'neovimPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('nv',nv)
    '''

    try:
        print("Opening the text editor neovim.");
        subprocess.check_output(neovimOpen)
        print("Success opening neovim.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the neovimOpen'. \n  Or use 'neovimPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('nv',nv)


def oc():
    ''' 
    DESCRIPTION:
    Open the data analysis program octave (open source analog of matlab) from within PyMOL.


    USAGE:
    oc

    ARGUMENTS:
    None
    EXAMPLE:
    oc

    MORE DETAILS:
    Open the data analysis program octave (open source analog of matlab) from within PyMOL.


    VERTICAL PML SCRIPT:
    subprocess.call(octaveOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(octaveOpen);return

    PYTHON CODE:
def oc():
    try:
        print("Opening octave.");
        subprocess.check_output(octaveOpen)
        print("Success opening octave.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'octaveOpen'. \n  Or use 'octavePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('oc',oc)
    '''

    try:
        print("Opening octave.");
        subprocess.check_output(octaveOpen)
        print("Success opening octave.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'octaveOpen'. \n  Or use 'octavePath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('oc',oc)


def omx(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of molecule along the x-axis of the viewport. 

    USAGE:
    omx <selection>

    ARGUMENTS:
    optional selection
    EXAMPLE:
    omx 1lw9

    MORE DETAILS:
    Align long axis of molecule along the x-axis of the viewport in the negative direction. 


    VERTICAL PML SCRIPT:
    orient selection; 
    rotate selection,z,180;

    HORIZONTAL PML SCRIPT:
    orient selection;rotate selection,z,180;

    PYTHON CODE:
def omx(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',180,selection); 

cmd.extend('omx',omx)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',180,selection); 

cmd.extend('omx',omx)


def omxy(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of molecule along minus x-y axis.

    USAGE:
    omxy

    ARGUMENTS:
    optional selection
    EXAMPLE:
    omxy

    MORE DETAILS:
    Align long axis of molecule along the minus x*y axis of the viewport.


    VERTICAL PML SCRIPT:
    cmd.orient(); 
    cmd.turn('z',315) 

    HORIZONTAL PML SCRIPT:
    cmd.orient();cmd.turn('z',315) 

    PYTHON CODE:
def omxy(selection='all'):
    cmd.orient(); 
    cmd.turn('z',315) 

cmd.extend('omxy',omxy)
    '''

    cmd.orient(); 
    cmd.turn('z',315) 

cmd.extend('omxy',omxy)


def omxyz(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of the selection along the mxyz axis. 

    USAGE:
    omxyz

    or 

    omxyz selection

    ARGUMENTS:
    optional selection
    EXAMPLE:
    omxyz 1lw9

    MORE DETAILS:
    Align long axis of the selection along the mxyz axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient selection; 
    rotate selection,z,315;
    rotate selection,y,315;
    HORIZONTAL PML SCRIPT:
    orient selection; rotate selection,z,135;rotate selection,y,135);
    PYTHON CODE:
def omxyz(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',315,selection); 
    cmd.rotate('y',315,selection);
cmd.extend('omxyz',omxyz)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',315,selection); 
    cmd.rotate('y',315,selection);
cmd.extend('omxyz',omxyz)


def omy(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of the selection along the y-axis of the viewport in the negative direction.

    USAGE:
    omy <selection>

    ARGUMENTS:
    selection
    EXAMPLE:
    omy;
    omy 1lw9

    MORE DETAILS:
    Align long axis of the selection along the y-axis of the viewport in the negative direction.


    VERTICAL PML SCRIPT:
    orient; 
    rotate <selection>, z,270

    HORIZONTAL PML SCRIPT:
    orient; rotate <selection>,z,270

    PYTHON CODE:
def omy(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',270,selection) 

cmd.extend('omy',omy)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',270,selection) 

cmd.extend('omy',omy)


def omz(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of selection along the z-axis of the viewport in the negative z direction. 

    USAGE:
    omz <selection>

    ARGUMENTS:
    selection
    EXAMPLE:
    omz;
    omz 1lw9

    MORE DETAILS:
    Align long axis of selection along the negative z-axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient; 
    rotate y,270

    HORIZONTAL PML SCRIPT:
    orient; rotate z,270

    PYTHON CODE:
def omz(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('y',270,selection) 

cmd.extend('omz',omz)
    '''

    cmd.orient(selection); 
    cmd.rotate('y',270,selection) 

cmd.extend('omz',omz)


def oni():
    ''' 
    DESCRIPTION:
    Open the editor Oni from within PyMOL. 

    USAGE:
    oni

    ARGUMENTS:
    None
    EXAMPLE:
    oni

    MORE DETAILS:
    Open the editor Oni from within PyMOL. 
    The is an editor based on neovim.


    VERTICAL PML SCRIPT:
    subprocess.call(oniOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(oniOpen);return

    PYTHON CODE:
def oni():
    try:
        print("Opening the text editor oni.");
        subprocess.check_output(oniOpen)
        print("Success opening oni.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'oniOpen'. \n  Or use 'oniPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('oni',oni)
    '''

    try:
        print("Opening the text editor oni.");
        subprocess.check_output(oniOpen)
        print("Success opening oni.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'oniOpen'. \n  Or use 'oniPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('oni',oni)


def ox(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of molecule along x-axis. 

    USAGE:
    ox <selection>

    ARGUMENTS:
    selection
    EXAMPLE:
    ox;
    ox 1lw9


    MORE DETAILS:
    Align long axis of molecule along the x-axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient selection
    HORIZONTAL PML SCRIPT:
    orient selection
    PYTHON CODE:
def ox(selection='all'):
    cmd.orient(selection) 

cmd.extend('ox',ox)
    '''

    cmd.orient(selection) 

cmd.extend('ox',ox)


def oxy(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of the selection along the x-y axis. 

    USAGE:
    oxy

     or

    oxy selection

    ARGUMENTS:
    optional selection
    EXAMPLE:
    oxy

    MORE DETAILS:
    Align long axis of the selection along the x-y axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient selection
    rotate selection,z,270 

    HORIZONTAL PML SCRIPT:
    orient selection;rotate selection,z,270 

    PYTHON CODE:
def oxy(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',135,selection) 

cmd.extend('oxy',oxy)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',135,selection) 

cmd.extend('oxy',oxy)


def oxyz(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of the selection along the xyz axis. 

    USAGE:
    oxyz

    or 

    oxyz selection

    ARGUMENTS:
    optional selection
    EXAMPLE:
    oxyz 1lw9

    MORE DETAILS:
    Align long axis of the selection along the xyz axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient selection; 
    rotate selection,z,135;
    rotate selection,y,135;
    HORIZONTAL PML SCRIPT:
    orient selection; rotate selection,z,315;rotate selection,y,315);
    PYTHON CODE:
def oxyz(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',135,selection); 
    cmd.rotate('y',135,selection);
cmd.extend('oxyz',oxyz)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',135,selection); 
    cmd.rotate('y',135,selection);
cmd.extend('oxyz',oxyz)


def oy(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of the selection along the y-axis of the viewport. 

    USAGE:
    oy <selection>

    ARGUMENTS:
    selection
    EXAMPLE:
    oy;
    oy 1lw9

    MORE DETAILS:
    Align long axis of the selection along the y-axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient; 
    rotate <selection>, z,90

    HORIZONTAL PML SCRIPT:
    orient; rotate <selection>,z,90

    PYTHON CODE:
def oy(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('z',90,selection) 

cmd.extend('oy',oy)
    '''

    cmd.orient(selection); 
    cmd.rotate('z',90,selection) 

cmd.extend('oy',oy)


def oz(selection='all'):
    ''' 
    DESCRIPTION:
    Align long axis of selection along the z-axis of the viewport. 

    USAGE:
    oz <selection>

    ARGUMENTS:
    selection
    EXAMPLE:
    oz;
    oz 1lw9

    MORE DETAILS:
    Align long axis of selection along the z-axis of the viewport. 


    VERTICAL PML SCRIPT:
    orient; 
    rotate y,90

    HORIZONTAL PML SCRIPT:
    orient; rotate z,90

    PYTHON CODE:
def oz(selection='all'):
    cmd.orient(selection); 
    cmd.rotate('y',90,selection) 

cmd.extend('oz',oz)
    '''

    cmd.orient(selection); 
    cmd.rotate('y',90,selection) 

cmd.extend('oz',oz)


def pairD(sel1, sel2, max_dist, output="N", sidechain="N", show="N"):
    ''' 
    DESCRIPTION:
    Find the pairwise distances between two selections.

    USAGE:
    pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]

    ARGUMENTS:
    sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]
    EXAMPLE:
    pairD chA., ch.B, 3.5, output=S, sidechain=Y, show=Y
pairD i. 80, i. 89, 20, output=P, sidechain=Y, show=Y

    MORE DETAILS:
    usage: pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]

    sel1 and sel2 can be any to pre-existing or newly defined selections

    max_dist: maximum distance in Angstrom between atoms in the two selections

    --optional settings:
    output: accepts Screen/Print/None (default N)
    sidechain: limits (Y) results to sidechain atoms (default N)
    show: shows (Y) individual distances in pymol menu (default=N)


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA

    PYTHON CODE:
def pairD(sel1, sel2, max_dist, output="N", sidechain="N", show="N"):
    print("")
    cmd.delete ("dist*")
    extra=""
    if sidechain=="Y": extra=" and not name c+o+n"

    #builds models
    m1=cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
    m1o=cmd.get_object_list(sel1)
    m2=cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
    m2o=cmd.get_object_list(sel2)

    #defines selections
    cmd.select("__tsel1a", sel1+" around "+str(max_dist)+" and "+sel2+extra)
    cmd.select("__tsel1", "__tsel1a and "+sel2+extra)
    cmd.select("__tsel2a", sel2+" around "+str(max_dist)+" and "+sel1+extra)
    cmd.select("__tsel2", "__tsel2a and "+sel1+extra)
    cmd.select("IntAtoms_"+max_dist, "__tsel1 or __tsel2")
    cmd.select("IntRes_"+max_dist, "byres IntAtoms_"+max_dist)

    #controlers-1
    if len(m1o)==0:
        print("warning, '"+sel1+extra+"' does not contain any atoms.")
        return
    if len(m2o)==0:
        print("warning, '"+sel2+extra+"' does not contain any atoms.")
        return

    #measures distances
    s=""
    counter=0
    for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
            distance=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord))))
            if distance<float(max_dist):
                s+="%s/%s/%s/%s/%s to %s/%s/%s/%s/%s: %.3f\n" % (m1o[0],m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2o[0],m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)
                counter+=1
                if show=="Y": cmd.distance (m1o[0]+" and "+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name, m2o[0]+" and "+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)

    #controler-2
    if counter==0:
        print("warning, no distances were measured! Check your selections/max_dist value")
        return

    #outputs
    if output=="S": print(s)
    if output=="P":
        f=open('IntAtoms_'+max_dist+'.txt','w')
        f.write("Number of distances calculated: %s\n" % (counter))
        f.write(s)
        f.close()
        print("Results saved in IntAtoms_%s.txt" % max_dist)
    print("Number of distances calculated: %s" % (counter))
    cmd.hide("lines", "IntRes_*")
    if show=="Y": cmd.show("lines","IntRes_"+max_dist)
    cmd.deselect()
cmd.extend('pairD', pairD)
    '''

    print("")
    cmd.delete ("dist*")
    extra=""
    if sidechain=="Y": extra=" and not name c+o+n"

    #builds models
    m1=cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
    m1o=cmd.get_object_list(sel1)
    m2=cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
    m2o=cmd.get_object_list(sel2)

    #defines selections
    cmd.select("__tsel1a", sel1+" around "+str(max_dist)+" and "+sel2+extra)
    cmd.select("__tsel1", "__tsel1a and "+sel2+extra)
    cmd.select("__tsel2a", sel2+" around "+str(max_dist)+" and "+sel1+extra)
    cmd.select("__tsel2", "__tsel2a and "+sel1+extra)
    cmd.select("IntAtoms_"+max_dist, "__tsel1 or __tsel2")
    cmd.select("IntRes_"+max_dist, "byres IntAtoms_"+max_dist)

    #controlers-1
    if len(m1o)==0:
        print("warning, '"+sel1+extra+"' does not contain any atoms.")
        return
    if len(m2o)==0:
        print("warning, '"+sel2+extra+"' does not contain any atoms.")
        return

    #measures distances
    s=""
    counter=0
    for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
            distance=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord))))
            if distance<float(max_dist):
                s+="%s/%s/%s/%s/%s to %s/%s/%s/%s/%s: %.3f\n" % (m1o[0],m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2o[0],m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)
                counter+=1
                if show=="Y": cmd.distance (m1o[0]+" and "+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name, m2o[0]+" and "+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)

    #controler-2
    if counter==0:
        print("warning, no distances were measured! Check your selections/max_dist value")
        return

    #outputs
    if output=="S": print(s)
    if output=="P":
        f=open('IntAtoms_'+max_dist+'.txt','w')
        f.write("Number of distances calculated: %s\n" % (counter))
        f.write(s)
        f.close()
        print("Results saved in IntAtoms_%s.txt" % max_dist)
    print("Number of distances calculated: %s" % (counter))
    cmd.hide("lines", "IntRes_*")
    if show=="Y": cmd.show("lines","IntRes_"+max_dist)
    cmd.deselect()
cmd.extend('pairD', pairD)


def pdbed():
    ''' 
    DESCRIPTION:
    Open PDBEditor.jar from within PyMOL. 

    USAGE:
    edpdb


    ARGUMENTS:
    None
    EXAMPLE:
    pdbed pdb filename

    MORE DETAILS:
    Open PDBEditor.jar from within PyMOL. 
    Adjust file path to location of the jar file.

    https://sourceforge.net/projects/pdbeditorjl/


    VERTICAL PML SCRIPT:
    print("Please wait. The PDBEditor is slow to start.")
    arg = ("java -jar" + pdbeditorPath + fileName)
    subprocess.call(arg,shell=True)
    return

    HORIZONTAL PML SCRIPT:
    print("Please wait. The PDBEditoris slow to start.");arg = ("java -jar " + pdbeditorPath + fileName);subprocess.call(arg,shell=True);return

    PYTHON CODE:
def pdbed():
    try:
        print("Opening PDB_Editor.");
        subprocess.check_output(pdbeditorOpen)
        print("Success opening PDBeditor.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'pdbeditorOpen'. \n  Or use 'pdbeditorPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('pdbed',pdbed)
    '''

    try:
        print("Opening PDB_Editor.");
        subprocess.check_output(pdbeditorOpen)
        print("Success opening PDBeditor.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'pdbeditorOpen'. \n  Or use 'pdbeditorPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('pdbed',pdbed)


def pdbremarks(filename):
    ''' 
    DESCRIPTION:
    Read REMARK lines from PDB file. 
Return dictionary with remarkNum as key and list of lines as value.
Called by the function quat().

    USAGE:
    pdbremarks(filename)

    ARGUMENTS:
    filename
    EXAMPLE:
    pdbremarks(filename)


    MORE DETAILS:
    Read REMARK lines from PDB file. 
    Return dictionary with remarkNum as key and list of lines as value.
    Called by the function quat().


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def pdbremarks(filename):
    remarks = dict()
    if not isinstance(filename, basestring):
        f = filename
    elif filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    return remarks

cmd.extend('pdbremarks', pdbremarks)
    '''

    remarks = dict()
    if not isinstance(filename, basestring):
        f = filename
    elif filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    return remarks

cmd.extend('pdbremarks', pdbremarks)


def ppt():
    ''' 
    DESCRIPTION:
    Open the powerpoint from within PyMOL. 

    USAGE:
    ppt

    ARGUMENTS:
    None
    EXAMPLE:
    ppt

    MORE DETAILS:
    Open the powerpoint from within PyMOL. 


    VERTICAL PML SCRIPT:
    subprocess.call(pptOpen); 
    return


    HORIZONTAL PML SCRIPT:
    subprocess.call(pptOpen); return


    PYTHON CODE:
def ppt():
    try:
        print("Opening the MS powerpoint.");
        subprocess.check_output(pptOpen)
        print("Success opening ppt.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'pptOpen'. \n  Or use 'pptPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('ppt',ppt)
    '''

    try:
        print("Opening the MS powerpoint.");
        subprocess.check_output(pptOpen)
        print("Success opening ppt.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'pptOpen'. \n  Or use 'pptPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('ppt',ppt)


def quat(name=None, filename=None, prefix=None, quiet=0):
    ''' 
    DESCRIPTION:
    Runs Thomas Holder's quat.py script to generate a biological unit using crystallographic symmetry.
    Requires a pdb file. The defult file type with the fetch command is *.cif.
    When using fetch, add the optional parameter type=pdb.
    Of course, you can also use type=pdb1 to retrieve the biological unit form the PDB. 
    Reads REMARK 350 from the pdb file  `filename` and creates the biological unit (quaternary structure).

    USAGE:
    quat [name [, filename [, prefix]]]

    ARGUMENTS:
    The code for a molecular object with symmetry information. 
    name = string: name of object and basename of PDB file, if filename is not given {default: first loaded object}
    filename = string: file path {default: <name>.pdb}
    prefix = string: prefix for new objects {default: <name>}

    EXAMPLE:
    fetch 4dgr
    quat 4dgr



    MORE DETAILS:
    Copyright 2010 to 2011 Thomas Holder, MPI for Developmental Biology
 
    Module for reading REMARK records from PDB files and in particular
    generate quaterny structure from REMARK 350.

    These notes below were added by Blaine Mooers on 12 September 2019.

    This content was made available under the GNU Free Documentation License 1.2

    See https://pymolwiki.org/index.php/BiologicalUnit/Quat for more information.

    quat is equavalent to biomolecule in the psico package.
 

    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def quat(name=None, filename=None, prefix=None, quiet=0):
    quiet = int(quiet)
    if name is None:
        name = cmd.get_object_list()[0]
    if prefix is None:
        prefix = name
    if filename is None:
        candidates = [
        '%s.pdb' % (name),
        '%s/%s.pdb' % (cmd.get('fetch_path'), name),
        '%s/%s/pdb%s.ent.gz' % (local_mirror_divided, name[1:3], name),
        ]
        for filename in candidates:
            if os.path.exists(filename):
                break
        else:
            print('please provide filename')
            return
        if not quiet:
            print('loading from %s' % (filename))
    remarks = pdbremarks(filename)
    if 350 not in remarks:
            print('There is no REMARK 350 in', filename)
            return
    quat = quat350(remarks[350])
    for chains in quat:
            matrices = quat[chains]
            for num in matrices:
                mat = matrices[num][0:12]
                mat.extend([0,0,0,1])
                copy = '%s_%d' % (prefix, num)
                if not quiet:
                    print('creating %s' % (copy))
                cmd.create(copy, '/%s//%s' % (name, '+'.join(chains)))
                cmd.alter(copy, 'segi="%d"' % (num))
                cmd.transform_object(copy, mat)
    cmd.disable(name)
    cmd.group('%s_quat' % (prefix), '%s_*' % (prefix))
cmd.extend('quat', quat)
    '''

    quiet = int(quiet)
    if name is None:
        name = cmd.get_object_list()[0]
    if prefix is None:
        prefix = name
    if filename is None:
        candidates = [
        '%s.pdb' % (name),
        '%s/%s.pdb' % (cmd.get('fetch_path'), name),
        '%s/%s/pdb%s.ent.gz' % (local_mirror_divided, name[1:3], name),
        ]
        for filename in candidates:
            if os.path.exists(filename):
                break
        else:
            print('please provide filename')
            return
        if not quiet:
            print('loading from %s' % (filename))
    remarks = pdbremarks(filename)
    if 350 not in remarks:
            print('There is no REMARK 350 in', filename)
            return
    quat = quat350(remarks[350])
    for chains in quat:
            matrices = quat[chains]
            for num in matrices:
                mat = matrices[num][0:12]
                mat.extend([0,0,0,1])
                copy = '%s_%d' % (prefix, num)
                if not quiet:
                    print('creating %s' % (copy))
                cmd.create(copy, '/%s//%s' % (name, '+'.join(chains)))
                cmd.alter(copy, 'segi="%d"' % (num))
                cmd.transform_object(copy, mat)
    cmd.disable(name)
    cmd.group('%s_quat' % (prefix), '%s_*' % (prefix))
cmd.extend('quat', quat)


def quat350(rem350):
    ''' 
    DESCRIPTION:
    Get transformation matrices for biomolecule 1 from REMARK 350.



    USAGE:
    quat350(rem350)



    ARGUMENTS:
    rem350


    EXAMPLE:
    quat350(rem350)



    MORE DETAILS:
    Get transformation matrices for biomolecule 1 from REMARK 350.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def quat350(rem350):
    biomt = dict()
    chains = tuple()
    seenbiomolecule = False
    for line in rem350:
        if line.startswith('BIOMOLECULE:'):
            if seenbiomolecule:
                break
            seenbiomolecule = True
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains = tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('                   AND CHAINS:'):
            chains += tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('  BIOMT'):
            row = int(line[7])
            num = int(line[8:12])
            vec = line[12:].split()
            vec = map(float, vec)
            biomt.setdefault(chains, dict()).setdefault(num, []).extend(vec)
    return biomt
cmd.extend('quat350', quat350)
    '''

    biomt = dict()
    chains = tuple()
    seenbiomolecule = False
    for line in rem350:
        if line.startswith('BIOMOLECULE:'):
            if seenbiomolecule:
                break
            seenbiomolecule = True
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains = tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('                   AND CHAINS:'):
            chains += tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('  BIOMT'):
            row = int(line[7])
            num = int(line[8:12])
            vec = line[12:].split()
            vec = map(float, vec)
            biomt.setdefault(chains, dict()).setdefault(num, []).extend(vec)
    return biomt
cmd.extend('quat350', quat350)


def rline():
    ''' 
    DESCRIPTION:
    Prints cheat sheet for the readline commands. 


    USAGE:
    rline



    ARGUMENTS:
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
    print(rline.__doc__)


cmd.extend("rline",rline)
    '''

    print(rline.__doc__)


cmd.extend("rline",rline)


def rmd():
    ''' 
    DESCRIPTION:
    Remove all measurement objects in the interal GUI.

    USAGE:
    rmd

    ARGUMENTS:
    None
    EXAMPLE:
    rmd

    MORE DETAILS:
    Remove all measurement objects in the interal GUI.
    Note that there is a "delete all measurements" toggle in the internal gui.


    VERTICAL PML SCRIPT:
    delete measure*;
    delete m*_*
    delete dist*
    HORIZONTAL PML SCRIPT:
    delete measure*; delete m*_*; delete dist*
    PYTHON CODE:
def rmd():
    cmd.do('delete measure*')
    cmd.do('delete m*_*')
    cmd.do('delete dist*')
cmd.extend("rmd", rmd)
    '''

    cmd.do('delete measure*')
    cmd.do('delete m*_*')
    cmd.do('delete dist*')
cmd.extend("rmd", rmd)


def rmhb(selection='all'):
    ''' 
    DESCRIPTION:
    Delete all H-bonds in the selection, which is all by default.




    USAGE:
    rmhb <selection>

    ARGUMENTS:
    The selection is optional. It is "all" by default.


    EXAMPLE:
    rmhb

 or

rmhb 1lw9



    MORE DETAILS:
    Delete all H-bonds in the selection, which is all by default.


    VERTICAL PML SCRIPT:
    delete hbonds


    HORIZONTAL PML SCRIPT:
    cmd.delete('hbonds')


    PYTHON CODE:
def rmhb(selection='all'):
    cmd.delete('hbonds')


cmd.extend('rmhb',rmhb)
    '''

    cmd.delete('hbonds')


cmd.extend('rmhb',rmhb)


def rmsc():

    ''' 
    DESCRIPTION:
    Remove supercell objects.

    USAGE:
    rmsc

    ARGUMENTS:
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
cmd.extend("rmsc", rmsc)

    '''

    cmd.do('delete supercell*;delete m*_*')
cmd.extend("rmsc", rmsc)



def rv(StoredView=0, decimal_places=2, outname="roundedview.txt"):
    ''' 
    DESCRIPTION:
    Get the view settings in a compact format on one line.


    USAGE:
    rv

or 

rv 1, 2, bestview

    ARGUMENTS:
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

cmd.extend("rv", rv)
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

    ARGUMENTS:
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
cmd.extend('saln',saln)
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".aln") 
cmd.extend('saln',saln)


def sasbdb():
    ''' 
    DESCRIPTION:
    Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB). 

    USAGE:
    sasbdb


    ARGUMENTS:
    NA
    EXAMPLE:
    sasbdb


    MORE DETAILS:
    Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB). 


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def sasbdb():
    url=sasbdbURL
    try:
        print("Opening the webpage of the Small Angle Scattering Biological Data Bank (SASBDB): a curated repository for small angle scattering data and models.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Small Angle Scattering Biological Data Bank (SASBDB): a curated repository for small angle scattering data and models.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('sasbdb',sasbdb)
    '''

    url=sasbdbURL
    try:
        print("Opening the webpage of the Small Angle Scattering Biological Data Bank (SASBDB): a curated repository for small angle scattering data and models.");
        # URL must be in single quotes
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Small Angle Scattering Biological Data Bank (SASBDB): a curated repository for small angle scattering data and models.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('sasbdb',sasbdb)


def sb(selection='(all)', name='bump_check', quiet=1):
    ''' 
    DESCRIPTION:
    Show van der Waals clashes as red discs. 

    USAGE:
    show_bumps [ selection [, name ]]

    ARGUMENTS:
    show_bumps [ selection [, name ]]

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create {default: bump_check}

    EXAMPLE:
    sb
	
    sb 1lw9

    sb i. 60:80

    MORE DETAILS:
    Show van der Waals clashes as red discs. 

    Source: http://pymolwiki.org/index.php/show_bumps
    (c) 2011 Thomas Holder, MPI for Developmental Biology

    License: BSD-2-Clause:
    https://opensource.org/licenses/BSD-2-Clause


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def sb(selection='(all)', name='bump_check', quiet=1):
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        if not int(quiet):
            print('VDW Strain in state %d: %f' % (state, strain))
    cmd.show_as('cgo', name)

cmd.extend('sb',sb)
    '''

    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        if not int(quiet):
            print('VDW Strain in state %d: %f' % (state, strain))
    cmd.show_as('cgo', name)

cmd.extend('sb',sb)


def sbgrid():
    ''' 
    DESCRIPTION:
    Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel.

    USAGE:
    sbgrid

    ARGUMENTS:
    NA
    EXAMPLE:
    sbgrid

    MORE DETAILS:
    Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def sbgrid():
    url=sbgridURL
    try:
        print("Opening the webpage of the Structural Biology Grid,SBGRID, YouTube Channel.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Structural Biology Grid,SBGRID, YouTube Channel.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('sbgrid',sbgrid)
    '''

    url=sbgridURL
    try:
        print("Opening the webpage of the Structural Biology Grid,SBGRID, YouTube Channel.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of the Structural Biology Grid,SBGRID, YouTube Channel.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('sbgrid',sbgrid)


def sc111():

    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 1 x 1 unit cells. 

    USAGE:
    sc111

    ARGUMENTS:
    None
    EXAMPLE:
    sc111

    MORE DETAILS:
    Make a lattice of 1 x 1 x 1 unit cells.
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\nsupercell 1, 1, 1, , orange, supercell111, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;supercell 1, 1, 1, , orange, supercell111, 1
    PYTHON CODE:
def sc111():

    supercell(1, 1, 1, object=None, color='orange', name='supercell111', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")

cmd.extend("sc111", sc111)

    '''

    supercell(1, 1, 1, object=None, color='orange', name='supercell111', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")

cmd.extend("sc111", sc111)



def sc112():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 1 x 2 unit cells

    USAGE:
    sc112

    ARGUMENTS:
    None
    EXAMPLE:
    sc112

    MORE DETAILS:
    Make a lattice of 1 x 1 x 2 unit cells.
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 1, 2, , orange, supercell112, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;supercell 1, 1, 2, , orange, supercell112, 1
    PYTHON CODE:
def sc112():
    supercell(1, 1, 2, object=None, color='orange', name='supercell112', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc112", sc112)

    '''

    supercell(1, 1, 2, object=None, color='orange', name='supercell112', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc112", sc112)



def sc113():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 1 x 3 unit cells.

    USAGE:
    sc113

    ARGUMENTS:
    None
    EXAMPLE:
    sc113

    MORE DETAILS:
    Make a lattice of 1 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 1, 3, , orange, supercell113, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 1, 3, , orange, supercell113, 1
    PYTHON CODE:
def sc113():
    supercell(1, 1, 3, object=None, color='orange', name='supercell113', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc113", sc113)

    '''

    supercell(1, 1, 3, object=None, color='orange', name='supercell113', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc113", sc113)



def sc121():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 2 x 1 unit cells. 

    USAGE:
    sc121

    ARGUMENTS:
    None
    EXAMPLE:
    sc121

    MORE DETAILS:
    Make a lattice of 1 x 2 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

 
    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 2, 1, , orange, supercell121, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 1, , orange, supercell121, 1
    PYTHON CODE:
def sc121():
    supercell(1, 2, 1, object=None, color='orange', name='supercell121', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc121", sc121)

    '''

    supercell(1, 2, 1, object=None, color='orange', name='supercell121', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc121", sc121)



def sc122():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 2 x 2 unit cells. 

    USAGE:
    sc122

    ARGUMENTS:
    None
    EXAMPLE:
    sc122

    MORE DETAILS:
    Make a lattice of 1 x 2 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 2, 2, , orange, supercell122, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 2, , orange, supercell122, 1
    PYTHON CODE:
def sc122():
    supercell(1, 2, 2, object=None, color='orange', name='supercell122', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc122", sc122)

    '''

    supercell(1, 2, 2, object=None, color='orange', name='supercell122', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc122", sc122)



def sc123():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 2 x 3 unit cells. 

    USAGE:
    sc123

    ARGUMENTS:
    None
    EXAMPLE:
    sc123

    MORE DETAILS:
    Make a lattice of 1 x 2 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 3, , orange, supercell123, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 2, 3, , orange, supercell123, 1
    PYTHON CODE:
def sc123():
    supercell(1, 2, 3, object=None, color='orange', name='supercell123', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc123", sc123)

    '''

    supercell(1, 2, 3, object=None, color='orange', name='supercell123', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc123", sc123)



def sc131():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 3 x 1 unit cells. 

    USAGE:
    sc131

    ARGUMENTS:
    None
    EXAMPLE:
    sc131

    MORE DETAILS:
    Make a lattice of 1 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 1, , orange, supercell131, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 1, , orange, supercell131, 1
    PYTHON CODE:
def sc131():
    supercell(1, 3, 1, object=None, color='orange', name='supercell131', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc131", sc131)

    '''

    supercell(1, 3, 1, object=None, color='orange', name='supercell131', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc131", sc131)



def sc132():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 3 x 2 unit cells. 

    USAGE:
    sc132

    ARGUMENTS:
    None
    EXAMPLE:
    sc132

    MORE DETAILS:
    Make a lattice of 1 x 3 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 1, 3, 2, , orange, supercell132, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 2, , orange, supercell132, 1
    PYTHON CODE:
def sc132():
    supercell(1, 3, 2, object=None, color='orange', name='supercell132', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc132", sc132)

    '''

    supercell(1, 3, 2, object=None, color='orange', name='supercell132', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc132", sc132)



def sc133():
    ''' 
    DESCRIPTION:
    Make a lattice of 1 x 3 x 3 unit cells. 

    USAGE:
    sc133

    ARGUMENTS:
    None
    EXAMPLE:
    sc133

    MORE DETAILS:
    Make a lattice of 1 x 3 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/mg18OU/supercell.py;
    supercell 1, 3, 3, , orange, supercell133, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 1, 3, 3, , orange, supercell133, 1
    PYTHON CODE:
def sc133():
    supercell(1, 3, 3, object=None, color='orange', name='supercell133', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc133", sc133)

    '''

    supercell(1, 3, 3, object=None, color='orange', name='supercell133', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc133", sc133)



def sc211():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 1 x 1 unit cells.

    USAGE:
    sc211

    ARGUMENTS:
    None
    EXAMPLE:
    sc211

    MORE DETAILS:
    Make a lattice of 2 x 1 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 1, , orange, supercell211, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 1, , orange, supercell211, 1
    PYTHON CODE:
def sc211():
    supercell(2, 1, 1, object=None, color='orange', name='supercell211', withmates=1)
    print('When finished, remove the symmetry mates with the command rmsc')

cmd.extend("sc211", sc211)

    '''

    supercell(2, 1, 1, object=None, color='orange', name='supercell211', withmates=1)
    print('When finished, remove the symmetry mates with the command rmsc')

cmd.extend("sc211", sc211)



def sc212():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 1 x 2 unit cells.

    USAGE:
    sc212

    ARGUMENTS:
    None
    EXAMPLE:
    sc212

    MORE DETAILS:
    Make a lattice of 2 x 1 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

 
    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 2, , orange, supercell212, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 2, , orange, supercell121, 1
    PYTHON CODE:
def sc212():
    supercell(2, 1, 2, object=None, color='orange', name='supercell212', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc212", sc212)

    '''

    supercell(2, 1, 2, object=None, color='orange', name='supercell212', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc212", sc212)



def sc213():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 1 x 3 unit cells. 

    USAGE:
    sc213

    ARGUMENTS:
    None
    EXAMPLE:
    sc213

    MORE DETAILS:
    Make a lattice of 2 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 1, 3, , orange, supercell213, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 1, 3, , orange, supercell213, 1
    PYTHON CODE:
def sc213():
    supercell(2, 1, 3, object=None, color='orange', name='supercell213', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc213", sc213)

    '''

    supercell(2, 1, 3, object=None, color='orange', name='supercell213', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc213", sc213)



def sc221():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 2 x 1 unit cells. 

    USAGE:
    sc221

    ARGUMENTS:
    None
    EXAMPLE:
    sc221

    MORE DETAILS:
    Make a lattice of 2 x 2 x 1 unit cells. 
    Use 'rmsc' to remove supercell objects.
    Requires Thomas Holder's supercell.py script.


    VERTICAL PML SCRIPT:
    supercell(2, 2, 1, object=None, color='orange', name='supercell221', withmates=1)

    HORIZONTAL PML SCRIPT:
    supercell(2, 2, 1, object=None, color='orange', name='supercell221', withmates=1)

    PYTHON CODE:
def sc221():
    supercell(2, 2, 1, object=None, color='orange', name='supercell221', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc221", sc221)

    '''

    supercell(2, 2, 1, object=None, color='orange', name='supercell221', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc221", sc221)



def sc222():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 2 x 2 unit cells

    USAGE:
    sc222

    ARGUMENTS:
    None
    EXAMPLE:
    sc222

    MORE DETAILS:
    Make a lattice of 2 x 2 x 2 unit cells. 
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 2, 2, 2, , orange, supercell222, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 2, 2, , orange, supercell222, 1
    PYTHON CODE:
def sc222():
    supercell(2, 2, 2, object=None, color='orange', name='supercell222', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc222", sc222)

    '''

    supercell(2, 2, 2, object=None, color='orange', name='supercell222', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc222", sc222)



def sc231():
    ''' 
    DESCRIPTION:
    Make a lattice of 2 x 3 x 1 unit cells. 

    USAGE:
    sc231

    ARGUMENTS:
    None
    EXAMPLE:
    sc231

    MORE DETAILS:
    Make a lattice of 2 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScriptssupercell.py;\n supercell 2, 3, 1, , orange, supercell231, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 2, 3, 1, , orange, supercell231, 1
    PYTHON CODE:
def sc231():
    supercell(2, 3, 1, object=None, color='orange', name='supercell231', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc231", sc231)
    '''

    supercell(2, 3, 1, object=None, color='orange', name='supercell231', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc231", sc231)


def sc311():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 1 x 1 unit cells. 

    USAGE:
    sc311

    ARGUMENTS:
    None
    EXAMPLE:
    sc311

    MORE DETAILS:
    Make a lattice of 3 x 1 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 1, , orange, supercell311, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 1, , orange, supercell311, 1
    PYTHON CODE:
def sc311():
    supercell(3, 1, 1, object=None, color='orange', name='supercell311', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc311", sc311)

    '''

    supercell(3, 1, 1, object=None, color='orange', name='supercell311', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc311", sc311)



def sc312():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 1 x 2 unit cells. 

    USAGE:
    sc312

    ARGUMENTS:
    None
    EXAMPLE:
    sc312

    MORE DETAILS:
    Make a lattice of 3 x 1 x 2 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 1, 2, , orange, supercell312, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 2, , orange, supercell312, 1
    PYTHON CODE:
def sc312():
    supercell(3, 2, 1, object=None, color='orange', name='supercell321', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc312", sc312)
    '''

    supercell(3, 2, 1, object=None, color='orange', name='supercell321', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc312", sc312)


def sc313():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 1 x 3 unit cells.

    USAGE:
    sc313

    ARGUMENTS:
    None
    EXAMPLE:
    sc313

    MORE DETAILS:
    Make a lattice of 3 x 1 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py
    supercell 3, 1, 3, , orange, supercell313, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 1, 3, , orange, supercell313, 1
    PYTHON CODE:
def sc313():
    supercell(3, 1, 3, object=None, color='orange', name='supercell313', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc313", sc313)

    '''

    supercell(3, 1, 3, object=None, color='orange', name='supercell313', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc313", sc313)



def sc321():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 2 x 1 unit cells. 

    USAGE:
    sc321

    ARGUMENTS:
    None
    EXAMPLE:
    sc321

    MORE DETAILS:
    Make a lattice of 3 x 2 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 2, 1, , orange, supercell321, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 2, 1, , orange, supercell321, 1
    PYTHON CODE:
def sc321():
    supercell(3, 2, 1, object=None, color='orange', name='supercell321', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc321", sc321)
    '''

    supercell(3, 2, 1, object=None, color='orange', name='supercell321', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc321", sc321)


def sc331():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 3 x 1 unit cells. 

    USAGE:
    sc331

    ARGUMENTS:
    None
    EXAMPLE:
    sc331

    MORE DETAILS:
    Make a lattice of 3 x 3 x 1 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.


    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 3, 1, , orange, supercell331, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 3, 1, , orange, supercell331, 1
    PYTHON CODE:
def sc331():
    supercell(3, 3, 1, object=None, color='orange', name='supercell331', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc331", sc331)

    '''

    supercell(3, 3, 1, object=None, color='orange', name='supercell331', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc331", sc331)



def sc333():
    ''' 
    DESCRIPTION:
    Make a lattice of 3 x 3 x 3 unit cells. 

    USAGE:
    sc333

    ARGUMENTS:
    None
    EXAMPLE:
    sc333

    MORE DETAILS:
    Make a lattice of 3 x 3 x 3 unit cells.
    Use "rmsc" to remove supercell objects.
    Requires supercell.py script by Thomas Holder.

 
    VERTICAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;
    supercell 3, 3, 3, , orange, supercell333, 1
    HORIZONTAL PML SCRIPT:
    run $HOME/Scripts/PyMOLScripts/supercell.py;\n supercell 3, 3, 3, , orange, supercell333, 1
    PYTHON CODE:
def sc333():
    supercell(3, 3, 3, object=None, color='orange', name='supercell333', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc333", sc333)

    '''

    supercell(3, 3, 3, object=None, color='orange', name='supercell333', withmates=1)
    print("Use 'rmsc' to remove supercell objects.")
cmd.extend("sc333", sc333)



def sccp4(stemName="saved"):
    ''' 
    DESCRIPTION:
    Save a ccp4 electron density map with a time stamp.

    USAGE:
    sccp4

    ARGUMENTS:
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
cmd.extend('sccp4',sccp4)
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

    ARGUMENTS:
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
cmd.extend('scif',scif)
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

    ARGUMENTS:
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
cmd.extend('sdae',sdae)
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

    ARGUMENTS:
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
cmd.extend('sdat',sdat)
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

    ARGUMENTS:
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
cmd.extend('sfasta',sfasta)
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

    ARGUMENTS:
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
cmd.extend('sidtf',sidtf)
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

    ARGUMENTS:
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
cmd.extend('smae',smae)
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

    ARGUMENTS:
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
cmd.extend('smmd',smmd)
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

    ARGUMENTS:
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
cmd.extend('smmod',smmod)
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

    ARGUMENTS:
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
cmd.extend('smoe',smoe)
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

    ARGUMENTS:
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
cmd.extend('smol',smol)
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

    ARGUMENTS:
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
cmd.extend('smol2',smol2)
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

    ARGUMENTS:
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
cmd.extend('smtl',smtl)
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

    ARGUMENTS:
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
cmd.extend('sobj',sobj)
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

    ARGUMENTS:
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
cmd.extend('sout',sout)
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

    ARGUMENTS:
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
cmd.extend('spdb',spdb)
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

    ARGUMENTS:
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
cmd.extend('spkl',spkl)
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

    ARGUMENTS:
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
cmd.extend('spkla',spkla)
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

    ARGUMENTS:
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
cmd.extend('spmo',spmo)
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

    ARGUMENTS:
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
cmd.extend('spng',spng)
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

    ARGUMENTS:
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
cmd.extend('spov',spov)
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

    ARGUMENTS:
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
cmd.extend('spqr',spqr)
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

    ARGUMENTS:
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
cmd.extend('spse',spse)
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".pse") 
cmd.extend('spse',spse)


def srv(StoredView=0, decimal_places=2, fileStemName="roundedview"):
    ''' 
    DESCRIPTION:
    Get the view settings in a compact format on one line. Save to file with timestamp appended to the stem of the filename.


    USAGE:
    srv

    ARGUMENTS:
    [view, decimal_places, outname]
Note that the values in the [] are optional.

    EXAMPLE:
    srv

or 

srv 1, 2, bestview

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
def srv(StoredView=0, decimal_places=2, fileStemName="roundedview"):
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
    DT=datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    myFile = open(fileStemName+s+".txt", "a")
    myFile.write(x.format(*myRoundedList) + "\n")
    myFile.close()
    return

cmd.extend("srv", srv)
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
    DT=datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    myFile = open(fileStemName+s+".txt", "a")
    myFile.write(x.format(*myRoundedList) + "\n")
    myFile.close()
    return

cmd.extend("srv", srv)


def ssdf(stemName="saved"):
    ''' 
    DESCRIPTION:
    Save session file with a time stamp.

    USAGE:
    ssdf sdfFileStemName

    ARGUMENTS:
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
cmd.extend('ssdf',ssdf)
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".sdf") 
cmd.extend('ssdf',ssdf)


def ssrlbl42():
    ''' 
    DESCRIPTION:
    Open the webpage of SSRL Biological SAXS at BL 4-2.



    USAGE:
    ssrlbl42



    ARGUMENTS:
    NA
    EXAMPLE:
    ssrlbl42



    MORE DETAILS:
    Open the webpage of SSRL Biological SAXS at BL 4-2.


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def ssrlbl42():
    url=ssrlbl42URL
    try:
        print("Opening the webpage of SSRL Biological SAXS at BL 4-2.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SSRL Biological SAXS at BL 4-2.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('ssrlbl42',ssrlbl42)
    '''

    url=ssrlbl42URL
    try:
        print("Opening the webpage of SSRL Biological SAXS at BL 4-2.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        webbrowser.open_new_tab(url)
        print("Opened the webpage of SSRL Biological SAXS at BL 4-2.")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found 


cmd.extend('ssrlbl42',ssrlbl42)


def st3(fileName="test.pml"):
    ''' 
    DESCRIPTION:
    Open Sublime Text 3 from within PyMOL.

    USAGE:
    st3

    ARGUMENTS:
    None
    EXAMPLE:
    st3

    MORE DETAILS:
    Open Open Sublime Text 3 from within PyMOL.
    Adjust file path to executable.


    VERTICAL PML SCRIPT:
    subprocess.call(sublOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(sublOpen);return

    PYTHON CODE:
def st3(fileName="test.pml"):
    try:
        print("Opening the text editor Sublime Text 3.");
        subprocess.Popen(sublimeText3Open)
        print("Success opening Sublime Text 3.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the sublpen'. \n  Or use 'sublimeText3Path' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('st3',st3)
    '''

    try:
        print("Opening the text editor Sublime Text 3.");
        subprocess.Popen(sublimeText3Open)
        print("Success opening Sublime Text 3.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the sublpen'. \n  Or use 'sublimeText3Path' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('st3',st3)


def supercell(a=1, b=1, c=1, object=None, color='blue', name='supercell', withmates=1):

    ''' 
    DESCRIPTION:
    Generate symmetry mates to fill unit cell. 
    (c) 2010 Thomas Holder




    USAGE:
    supercell a, b, c [, object [, color [, name [, withmates]]]]

    ARGUMENTS:
    a, b, c = integer: repeat cell in x,y,z direction a,b,c times {default: 1,1,1}
 
object = string: name of object to take cell definition from
 
color = string: color of cell {default: blue}
 
name = string: name of the cgo object to create {default: supercell}
 
withmates = bool: also create symmetry mates in displayed cells {default: 1}

    EXAMPLE:
    NA

    MORE DETAILS:
    Draw a supercell, as requested by Nicolas Bock on the pymol-users
    mailing list (Subject: [PyMOL] feature request: supercell construction
    Date: 04/12/2010 10:12:17 PM (Mon, 12 Apr 2010 14:12:17 -0600))

    source: https://pymolwiki.org/index.php/Supercell

    Calls two other functions: cellbasis and symexpcell.
    Both of these functions are also in this collection.

    (c) 2010 Thomas Holder

    GNU Free Documentation License 1.2


    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def supercell(a=1, b=1, c=1, object=None, color='blue', name='supercell', withmates=1):

    if object is None:
        object = cmd.get_object_list()[0]
    withmates = int(withmates)
    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
 
    basis = cellbasis(cell_angles, cell_edges)
    assert isinstance(basis, numpy.ndarray)
 
    ts = list()
    for i in range(int(a)):
        for j in range(int(b)):
            for k in range(int(c)):
                ts.append([i,j,k])
    obj = [
            cgo.BEGIN,
            cgo.LINES,
            cgo.COLOR,
    ]
    obj.extend(cmd.get_color_tuple(color))
 
    for t in ts:
        shift = basis[0:3,0:3] * t
        shift = shift[:,0] + shift[:,1] + shift[:,2]
 
        for i in range(3):
            vi = basis[0:3,i]
            vj = [
                numpy.array([0.,0.,0.]),
                basis[0:3,(i+1)%3],
                basis[0:3,(i+2)%3],
                basis[0:3,(i+1)%3] + basis[0:3,(i+2)%3]
                ]
            for j in range(4):
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j]).tolist())
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j] + vi).tolist())
 
        if withmates:
            symexpcell('m%d%d%d_' % tuple(t), object, *t)
 
    obj.append(cgo.END)
    cmd.delete(name)
    cmd.load_cgo(obj, name)

cmd.extend('supercell', supercell)
    '''

    if object is None:
        object = cmd.get_object_list()[0]
    withmates = int(withmates)
    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
 
    basis = cellbasis(cell_angles, cell_edges)
    assert isinstance(basis, numpy.ndarray)
 
    ts = list()
    for i in range(int(a)):
        for j in range(int(b)):
            for k in range(int(c)):
                ts.append([i,j,k])
    obj = [
            cgo.BEGIN,
            cgo.LINES,
            cgo.COLOR,
    ]
    obj.extend(cmd.get_color_tuple(color))
 
    for t in ts:
        shift = basis[0:3,0:3] * t
        shift = shift[:,0] + shift[:,1] + shift[:,2]
 
        for i in range(3):
            vi = basis[0:3,i]
            vj = [
                numpy.array([0.,0.,0.]),
                basis[0:3,(i+1)%3],
                basis[0:3,(i+2)%3],
                basis[0:3,(i+1)%3] + basis[0:3,(i+2)%3]
                ]
            for j in range(4):
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j]).tolist())
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j] + vi).tolist())
 
        if withmates:
            symexpcell('m%d%d%d_' % tuple(t), object, *t)
 
    obj.append(cgo.END)
    cmd.delete(name)
    cmd.load_cgo(obj, name)

cmd.extend('supercell', supercell)


def swrl(stemName="saved"):
    ''' 
    DESCRIPTION:
    Save wrl (VRML 2 file format) file with a time stamp.



    USAGE:
    swrl wrlFileStemName

    ARGUMENTS:
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
cmd.extend('swrlf',swrl)
    '''

    DT =datetime.datetime.now().strftime("yr%Ymo%mday%dhr%Hmin%Msec%S")
    s = str(DT)
    cmd.save(stemName+s+".wrl") 
cmd.extend('swrlf',swrl)


def symexpcell(prefix='mate', object=None, a=0, b=0, c=0):
    ''' 
    DESCRIPTION:
    Creates all symmetry-related objects for the specified object that occur with their bounding box center within the unit cell.
    (c) 2010 Thomas Holder



    USAGE:
    symexpcell prefix, object, [a, b, c]

    ARGUMENTS:
    prefix = string: prefix of new objects

object = string: object for which to create symmetry mates

a, b, c = integer: create neighboring cell {default: 0,0,0}

    EXAMPLE:
    NA

    MORE DETAILS:
    Creates all symmetry-related objects for the specified object that occur with their 
    bounding box center within the unit cell.

    Used by the function supercell.

    Also see symexp, http://www.pymolwiki.org/index.php/SuperSym.

    (c) 2010 Thomas Holder

    VERTICAL PML SCRIPT:
    NotYet
    HORIZONTAL PML SCRIPT:
    NotYet
    PYTHON CODE:
def symexpcell(prefix='mate', object=None, a=0, b=0, c=0):
    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
    spacegroup = sym[6]

    basis = cellbasis(cell_angles, cell_edges)
    basis = numpy.matrix(basis)

    extent = cmd.get_extent(object)
    center = sum(numpy.array(extent)) * 0.5
    center = numpy.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center

    extra_shift = [[float(i)] for i in (a,b,c)]

    i = 0
    matrices = xray.sg_sym_to_mat_list(spacegroup)
    for mat in matrices:
        i += 1

        mat = numpy.matrix(mat)
        shift = numpy.floor(mat * center_cell)
        mat[0:3,3] -= shift[0:3,0]
        mat[0:3,3] += extra_shift

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = '%s%d' % (prefix, i)
        cmd.create(name, object)
        cmd.transform_object(name, mat_list, 0)
        cmd.color(i+1, name)
cmd.extend('symexpcell', symexpcell)

    '''

    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
    spacegroup = sym[6]

    basis = cellbasis(cell_angles, cell_edges)
    basis = numpy.matrix(basis)

    extent = cmd.get_extent(object)
    center = sum(numpy.array(extent)) * 0.5
    center = numpy.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center

    extra_shift = [[float(i)] for i in (a,b,c)]

    i = 0
    matrices = xray.sg_sym_to_mat_list(spacegroup)
    for mat in matrices:
        i += 1

        mat = numpy.matrix(mat)
        shift = numpy.floor(mat * center_cell)
        mat[0:3,3] -= shift[0:3,0]
        mat[0:3,3] += extra_shift

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = '%s%d' % (prefix, i)
        cmd.create(name, object)
        cmd.transform_object(name, mat_list, 0)
        cmd.color(i+1, name)
cmd.extend('symexpcell', symexpcell)



def term():
    ''' 
    DESCRIPTION:
    Open a Terminal window on MacOS.



    USAGE:
    term

    ARGUMENTS:
    None
    EXAMPLE:
    term

    MORE DETAILS:
    Open a Terminal window on MacOS.


    VERTICAL PML SCRIPT:
    subprocess.call(terminalCommand)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(terminalCommand);return

    PYTHON CODE:
def term():
    try:
        print("Opening a terminal.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.check_output(terminalOpen)
        print("Opened a terminal.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('term',term)
    '''

    try:
        print("Opening a terminal.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.check_output(terminalOpen)
        print("Opened a terminal.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('term',term)


def timcolor(selection='all'):
    ''' 
    DESCRIPTION:
    Tim Mather's biophysical coloring scheme for proteins.




    USAGE:
    timcolor <selection>

    ARGUMENTS:
     <selection>
    EXAMPLE:
    timcolor 1lw9

    MORE DETAILS:
    Tim Mather's biophyscial coloring scheme for proteins.
    The scheme is defined in the dictionary below.
    This scheme applies only to proteins. 
    See the yrb shortcut for an alternate coloring scheme. 

    code={'acid'    :  'red'    ,
          'basic'   :  'blue'   ,
          'nonpolar':  'orange' ,
          'polar'   :  'green'  ,
          'cys'     :  'yellow'}

    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def timcolor(selection='all'):
    print(selection)
    code={'acid'    :  'red'    ,
          'basic'   :  'blue'   ,
          'nonpolar':  'orange' ,
          'polar'   :  'green'  ,
          'cys'     :  'yellow'}

    cmd.select('calcium', 'resn ca or resn cal')
    cmd.select('acid', 'resn asp or resn glu or resn cgu')
    cmd.select('basic', 'resn arg or resn lys or resn his')
    cmd.select('nonpolar', 'resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
    cmd.select('polar', 'resn ser or resn thr or resn asn or resn gln or resn tyr')
    cmd.select('cys', 'resn cys or resn cyx')
    cmd.select('backbone', 'name ca or name n or name c or name o')

    cmd.select ('none')
    for elem in code:
        line='color '+code[elem]+','+elem+'&'+selection
        print(line)
        cmd.do (line)
    word='color white,backbone &'+selection
    print(word)
    cmd.do (word)  #Used to be in code, but looks like
                              #dictionnaries are accessed at random
    cmd.hide ('everything','resn HOH')
    cmd.show('surface',selection)

cmd.extend('timcolor',timcolor)
    '''

    print(selection)
    code={'acid'    :  'red'    ,
          'basic'   :  'blue'   ,
          'nonpolar':  'orange' ,
          'polar'   :  'green'  ,
          'cys'     :  'yellow'}

    cmd.select('calcium', 'resn ca or resn cal')
    cmd.select('acid', 'resn asp or resn glu or resn cgu')
    cmd.select('basic', 'resn arg or resn lys or resn his')
    cmd.select('nonpolar', 'resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
    cmd.select('polar', 'resn ser or resn thr or resn asn or resn gln or resn tyr')
    cmd.select('cys', 'resn cys or resn cyx')
    cmd.select('backbone', 'name ca or name n or name c or name o')

    cmd.select ('none')
    for elem in code:
        line='color '+code[elem]+','+elem+'&'+selection
        print(line)
        cmd.do (line)
    word='color white,backbone &'+selection
    print(word)
    cmd.do (word)  #Used to be in code, but looks like
                              #dictionnaries are accessed at random
    cmd.hide ('everything','resn HOH')
    cmd.show('surface',selection)

cmd.extend('timcolor',timcolor)


def tvdw(arg1='all'):
    ''' 
    DESCRIPTION:
    Transparent vdw surface by Bobby Patton at Colorato State University. 

    USAGE:
    vdwTrans selection

    ARGUMENTS:
    selection

    EXAMPLE:
    vdwTrans 3nd3

    MORE DETAILS:
    vdwTrans creates a copy of an object with full-sized, transparent spheres.
    This looks cool when combined with the bs shortcut. 
    Bondi VDW values added below to override default Pymol settings
    The code is from  Bobby Paton at Colorado State University
     https://gist.githubusercontent.com/bobbypaton/\
    1cdc4784f3fc8374467bae5eb410edef/raw/\
    9995d51d6a64b8bcf01590c944cc38059b2f8d7f/pymol_style.py


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def tvdw(arg1='all'):
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

   #van der Waals settings
    cmd.copy(arg1+"_vdw", arg1)
    cmd.set("sphere_scale",1.0, arg1+"_vdw and elem H")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1+"_vdw")
    cmd.hide("nonbonded", arg1+"_vdw")
    cmd.hide("lines", arg1+"_vdw")
    cmd.hide("sticks", arg1+"_vdw")
    cmd.hide("cartoon", arg1+"_vdw")
    cmd.show("spheres", arg1+"_vdw")
    cmd.set("sphere_transparency", 0.7, arg1+"_vdw")

    print("Note that the selection of 'all' does not work when applying this function")
    print("to multiple models as are found in pdb1 files and when multiple chains")
    print("make up the biological unit.")
    print("The shortcut has to be applied separately to each model.")


cmd.extend('tvdw', tvdw)
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

   #van der Waals settings
    cmd.copy(arg1+"_vdw", arg1)
    cmd.set("sphere_scale",1.0, arg1+"_vdw and elem H")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1+"_vdw")
    cmd.hide("nonbonded", arg1+"_vdw")
    cmd.hide("lines", arg1+"_vdw")
    cmd.hide("sticks", arg1+"_vdw")
    cmd.hide("cartoon", arg1+"_vdw")
    cmd.show("spheres", arg1+"_vdw")
    cmd.set("sphere_transparency", 0.7, arg1+"_vdw")

    print("Note that the selection of 'all' does not work when applying this function")
    print("to multiple models as are found in pdb1 files and when multiple chains")
    print("make up the biological unit.")
    print("The shortcut has to be applied separately to each model.")


cmd.extend('tvdw', tvdw)


def tvdwbw(arg1='all'):
    ''' 
    DESCRIPTION:
    Transparent vdw surface by Bobby Patton at Colorato State University and grayscale by gscale function.

    USAGE:
    tvdwbw

    ARGUMENTS:
    selection

    EXAMPLE:
    tvdwbw

    MORE DETAILS:
    Transparent vdw surface by Bobby Patton at Colorato State University and grayscale by gscale function.
    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def tvdwbw(arg1='all'):
    tvdw(arg1)
    gscale(ag1)

cmd.extend('tvdwbw', tvdwbw)
    '''

    tvdw(arg1)
    gscale(ag1)

cmd.extend('tvdwbw', tvdwbw)


def vim():
    ''' 
    DESCRIPTION:
    Open vim from within PyMOL.

    USAGE:
    vim

    ARGUMENTS:
    None
    EXAMPLE:
    vim script.pml

    MORE DETAILS:
    Open vim from within PyMOL. 
    Adjust file path to vim on your computer.

 
    VERTICAL PML SCRIPT:
    subprocess.call(vimCommand)
    return

    HORIZONTAL PML SCRIPT:
    subprocess.call(vimCommand);return

    PYTHON CODE:
def vim():
    try:
        print("Opening the molecular graphics program vim.");
        subprocess.check_output(vimOpen)
        print("Success opening vim.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the vimOpen'. \n  Or use 'vimPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('vim',vim)
    '''

    try:
        print("Opening the molecular graphics program vim.");
        subprocess.check_output(vimOpen)
        print("Success opening vim.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the vimOpen'. \n  Or use 'vimPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('vim',vim)


def vmd():
    ''' 
    DESCRIPTION:
    Open vmd from within PyMOL.

    USAGE:
    vmd

    ARGUMENTS:
    None
    EXAMPLE:
    vmd

    MORE DETAILS:
    Open vmd from within PyMOL.
    Adjust file path for your location of vmd.


    VERTICAL PML SCRIPT:
    subprocess.call(vmdOpen);
return

    HORIZONTAL PML SCRIPT:
    subprocess.call(vmdOpen);return
    PYTHON CODE:
def vmd():
    try:
        print("Opening the molecular graphics program VMD.");
        subprocess.check_output(vmdOpen)
        print("Success opening VMD.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'vmdOpen'. \n  Or use 'vmdPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('vmd',vmd)
    '''

    try:
        print("Opening the molecular graphics program VMD.");
        subprocess.check_output(vmdOpen)
        print("Success opening VMD.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'vmdOpen'. \n  Or use 'vmdPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('vmd',vmd)


def weather():
    ''' 
    DESCRIPTION:
    Open National Weather Service website for locale. 

    USAGE:
    weather

    ARGUMENTS:
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
    url=weatherServiceRadarURL
    try:
        print("Trying to open National Weather Service website.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening National Weather Service website")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('weather',weather)
    '''

    url=weatherServiceRadarURL
    try:
        print("Trying to open National Weather Service website.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening National Weather Service website")
    except Exception as e:
        print("Webbrowser error: " % e) # prints error if browser is not found

cmd.extend('weather',weather)


def webmail():
    ''' 
    DESCRIPTION:
    Open Web Mail in defualt browser.

    USAGE:
    webmail

    ARGUMENTS:
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
    url=webmailURL
    try:
        print("Trying to open webmail.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening webmail.")
    except Exception as e:
        print("Webbroser error: " % e) # prints error if browser is not found

cmd.extend('webmail',webmail)
    '''

    url=webmailURL
    try:
        print("Trying to open webmail.");
       # Try to open the default webrower
        client = webbrowser.get()
        # URL must be in single quotes
        client.open_new_tab(url)
        print("Success opening webmail.")
    except Exception as e:
        print("Webbroser error: " % e) # prints error if browser is not found

cmd.extend('webmail',webmail)


def word():
    ''' 
    DESCRIPTION:
    Open word from within PyMOL.

    USAGE:
    word

    ARGUMENTS:
    None
    EXAMPLE:
    word paper.docx

    MORE DETAILS:
    Open MS Word from within PyMOL. 
    Adjust file path to MS Word on your computer.
	

    VERTICAL PML SCRIPT:
    subprocess.call(wordOpen);
return

    HORIZONTAL PML SCRIPT:
    subprocess.call(wordOpen);return

    PYTHON CODE:
def word():
    try:
        print("Opening the program Microsoft Word.");
        subprocess.check_output(wordOpen)
        print("Success opening Microsoft Word.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the wordopen'. \n  Or use 'wordPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('word',word)
    '''

    try:
        print("Opening the program Microsoft Word.");
        subprocess.check_output(wordOpen)
        print("Success opening Microsoft Word.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the wordopen'. \n  Or use 'wordPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found

cmd.extend('word',word)


def x11():
    ''' 
    DESCRIPTION:
    Open x11 terminal.

    USAGE:
    x11

    ARGUMENTS:
    None
    EXAMPLE:
    x11

    MORE DETAILS:
    Open x11 terminal.


    VERTICAL PML SCRIPT:
    subprocess.call(x11Open)
return
    HORIZONTAL PML SCRIPT:
    subprocess.call(x11Open);return
    PYTHON CODE:
def x11():
    try:
        print("Opening an X11 window.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.check_output(x11Open)
        print("Opened an X11 window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('x11',x11)
    '''

    try:
        print("Opening an X11 window.");                               
        # URL must be in single quotes                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        subprocess.check_output(x11Open)
        print("Opened an X11 window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('x11',x11)


def xquartz():
    ''' 
    DESCRIPTION:
    Open new XQuartz terminal on MacOS. 


    USAGE:
    xquartz

    ARGUMENTS:
    None
    EXAMPLE:
    xquartz

    MORE DETAILS:
    Open new XQuartz terminal on MacOS. 
    Does not work on Mohave. 
    Recommend installing and using iterm or term. 

    VERTICAL PML SCRIPT:
    subprocess.call(xquartzOpen);
    return

    HORIZONTAL PML SCRIPT:
    subprocess.callxquartzOpen);return

    PYTHON CODE:
def xquartz():
    try:
        print("Opening a new XQuartz window.");                               
        subprocess.Popen(xquartzOpen)
        print("Opened a new XQuartz window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('xquartz',xquartz)
    '''

    try:
        print("Opening a new XQuartz window.");                               
        subprocess.Popen(xquartzOpen)
        print("Opened a new XQuartz window.")
    except Exception as e:
        print("Subprocess error: " % e) # prints error if browser is not found 


cmd.extend('xquartz',xquartz)


def yasara(fileName="test.pml"):
    ''' 
    DESCRIPTION:
    Open the molecular graphics prograom YASASRA from within PyMOL.

    USAGE:
    yasara

    ARGUMENTS:
    None
    EXAMPLE:
    yasara

    MORE DETAILS:
    Open the molecular graphics prograom YASASRA from within PyMOL.


    VERTICAL PML SCRIPT:
    subprocess.call(yasaraOpen);
return


    HORIZONTAL PML SCRIPT:
        subprocess.call(yasaraOpen);return

    PYTHON CODE:
def yasara(fileName="test.pml"):
    try:
        print("Opening the molecular graphics program Yasara.");
        subprocess.check_output(yasaraOpen)
        print("Success opening Yasara.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'yasaraOpen'. \n  Or use 'yasaraPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('yasara',yasara)
    '''

    try:
        print("Opening the molecular graphics program Yasara.");
        subprocess.check_output(yasaraOpen)
        print("Success opening Yasara.")
    except subprocess.CalledProcessError:
        print("Executable not found! \n  Check syntax of the 'yasaraOpen'. \n  Or use 'yasaraPath' as the argument of check_output().")
        pass # handle errors in the called executable
    except OSError:
        pass # executable not found
cmd.extend('yasara',yasara)


def yrb(selection='all'):
    ''' 
    DESCRIPTION:
    A script to highlight hydrophobicity and charge on protein surfaces (Hagemans et al. 2015).

    USAGE:
    "yrb" to colour all structures
"yrb 'selection'" to colour only that specific structure


    ARGUMENTS:
    None to apply to all or a selection
    EXAMPLE:
    yrb
yrb 1lw9

    MORE DETAILS:
    A script to highlight hydrophobicity and charge on protein surfaces.

    @Article{Hagemans2015AScriptToHighlightHydrophobicityAndChargeOnProteinSurfaces,
      author   = {Hagemans, Dominique and van Belzen, Ianthe A. E. M. and Morán Luengo, Tania and Rüdiger, Stefan G. D.},
      title    = {A script to highlight hydrophobicity and charge on protein surfaces},
      journal  = {Frontiers in Molecular Biosciences},
      year     = {2015},
      volume   = {2},
      pages    = {56},
      issn     = {2296-889X},
      abstract = {The composition of protein surfaces determines both affinity and specificity of protein-protein interactions. Matching of hydrophobic contacts and charged groups on both sites of the interface are crucial to ensure specificity. Here, we propose a highlighting scheme, YRB, which highlights both hydrophobicity and charge in protein structures. YRB highlighting visualizes hydrophobicity by highlighting all carbon atoms that are not bound to nitrogen and oxygen atoms. The charged oxygens of glutamate and aspartate are highlighted red and the charged nitrogens of arginine and lysine are highlighted blue. For a set of representative examples, we demonstrate that YRB highlighting intuitively visualizes segments on protein surfaces that contribute to specificity in protein-protein interfaces, including Hsp90/co-chaperone complexes, the SNARE complex and a transmembrane domain. We provide YRB highlighting in form of a script that runs using the software PyMOL.},
      doi      = {10.3389/fmolb.2015.00056},
      url      = {https://www.frontiersin.org/article/10.3389/fmolb.2015.00056},
    }

    Created by Dominique Hagemans and Ianthe A.E.M. van Belzen, July 2015
    Rudiger group CPC, Utrecht University

    yellow: C, CH, CH2, CH3 groups that are not bound to N or O groups.
    red: negatively charged atoms
    blue: positively charged atoms
    grey: backbone, polar groups and remaining atoms


    VERTICAL PML SCRIPT:
    NA
    HORIZONTAL PML SCRIPT:
    NA
    PYTHON CODE:
def yrb(selection='all'):
    cmd.remove("hydro")
    cmd.set_color('yellow',[0.950,0.78,0.0])
    cmd.set_color('grey',[0.95,0.95,0.95])
    cmd.set_color('red',[1.0,0.4,0.4])  
    cmd.set_color('blue',[0.2,0.5,0.8]) 

    mapping = {}
    mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
    mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
    mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
    mapping['cys'] = [ ('SG', 'grey') ] 
    mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
    mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
    mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ] 
    mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
    mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
    mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
    mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
    mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
    mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
    mapping['ser'] = [ ('CB,OG', 'grey') ]
    mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
    mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
    mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
    mapping['val'] = [ ('CG1,CG2', 'yellow') ]

    obj_list = cmd.get_names('objects')
    for obj in obj_list:
        if (obj == selection or selection == 'all'):
            cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
            cmd.color('yellow','(n. CB and ' + obj + ')')

            for key in mapping:
                for (atom, color) in mapping[key]:
                    cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

cmd.extend('yrb',yrb)
    '''

    cmd.remove("hydro")
    cmd.set_color('yellow',[0.950,0.78,0.0])
    cmd.set_color('grey',[0.95,0.95,0.95])
    cmd.set_color('red',[1.0,0.4,0.4])  
    cmd.set_color('blue',[0.2,0.5,0.8]) 

    mapping = {}
    mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
    mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
    mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
    mapping['cys'] = [ ('SG', 'grey') ] 
    mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
    mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
    mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ] 
    mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
    mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
    mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
    mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
    mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
    mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
    mapping['ser'] = [ ('CB,OG', 'grey') ]
    mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
    mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
    mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
    mapping['val'] = [ ('CG1,CG2', 'yellow') ]

    obj_list = cmd.get_names('objects')
    for obj in obj_list:
        if (obj == selection or selection == 'all'):
            cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
            cmd.color('yellow','(n. CB and ' + obj + ')')

            for key in mapping:
                for (atom, color) in mapping[key]:
                    cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

cmd.extend('yrb',yrb)



 
""" Print the shortcuts on startup of PyMOL"""
print(SC.__doc__)
