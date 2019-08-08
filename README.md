# pymolshortcuts

This repository for ***pymolschortucts.py*** which contains 192 functions in 20 categories mapped to short-names that work like aliases (the table below is incomplete).
These shortcuts include many convenience functions that make work in PyMOL more productive and fun!
Some shortcuts save the many hours of work required to assemble a new script file while other shortcuts save only a few minutes but lower motivational barriers.
You do not need to understand Python to be able run the script.
To get a quick overview of this project, see the slides in the pdf file *mooers23jul2019ACAPyMOLshortcuts.pdf*.

To download the *pymolshortcuts.py* file, git clone the repository. 
Otherwise, left-click on  **pymolshortcuts.py** above and left-click "Raw".
The raw file will be displayed in your browser.
Select save under the file pull-down of your browser
Then save the file **pymolshortcuts.py** to your home direcotry.

To have the shortcuts always available in PyMOL, add the command 'run ~/pymolshortcuts.py' to your *.pymolrc*, *pymolrc* or *pymolrc.pml* file in your home directory to load the functions in *pymolshortcuts.py* on startup of PyMOL.
The pymolrc file is an optional file.
You may have to creat it with a text editor if you have done so already.
If you do not have text editor, you can use PyMOL's built-in text editor.
Go to File --> Edit pymolrc.
In spite of the name *run*, the functions will be loaded into memory but will not be executed.
You may want to store the script **pymolshortcut.py** in a safer place than your home directory.
I store mine in */Users/blaine/Scripts/PyMOLScripts/*

Enter **SC** at the upper **PyMOL>** prompt the to get a list of shortcuts printed to the command history window.
Use the **help** function to see the documentation for each shortcut. 
For examples, enter **help PW** to print the documentation for the shortcut **PW** to the command history window. The documentation has four sections: a description of what the shortcut does, an example of running the shortcut, 
the corresponding pml code with one command per line for easy reuse in a script file, and all of the commands on a single line for re-use on the command line as hortizontal script. 

While most shortcuts are ready to use, some shortcuts require the file paths to the executables on your computer.
The sites where edits are required are marked with ** >>> ** in the script.
You can search for ** >>> ** in a good text editor like Microsoft's free Visual Studio Code (avialable for Mac, Windows, and Linux).

Some shortcuts require additional scripts that available on the PyMOL Wiki.
Error messages will indicate the missing scripts. 
Otherwise, search **pymolshortcut.py** for ** .py ** to find the names of the missing scripts.

A few shortcuts require access to Python modules that were not packaged with the incentive version of PyMOL.
Paste the following command at the upper PyMOL> prompt.

```bash
conda install requests beautifulsoup4 datetime
```

Be patient. The prompt can appear to hang for five to ten minutes while the installation occurs. 


The shortcut **PW** takes one or more search terms and the sends them to the PyMOL Wiki.
A browser tab opens for each search term, so multiple searches are run in parallel while you continue your work in PyMOL.
You can inspect the results of the searches when there is a natural break in your work in PyMOL.
Other search functions can submit parallel searches of PubMed, Google, bioRxiv, Research Gate, GitHub, and more.
See the tables below.

Another class of shortcuts launch your favorite full-featured text editor from within PyMOL.
You have to install the text editor and edit the file path to the executable.

Another class of shortcuts opens the manuscript or grant application that you are working on in Overleaf.
You have to edit the script by pasting in the appropriate link to the specific document.

Another class of shortcuts saves files with timestamps embedded in the filename to avoid overwriting png, pdb, pse, and other types of files written out from PyMOL.
These save function names begin with **s**, (e.g., **spse filename** saves the current session with date and time to nearest second embedded in the file name).
You can delete the unwanted version(s) at a latter time. 
These functions are useful if you do not have these files under version control.


## Print shortcuts and their descriptions

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|SC        |Print to screen list of the shortcuts that are available in the script pymolshortcuts.py.   |
|github    |Print the url of the README file for the pymolshortcuts repository. |

## Show many models (NMR and crystal packing)

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|nmr       |Show all of the models in nmr structure.                       |
|nmroff    |Hide all but first model in a nmr structure.                   |
|rmsc      |Remove supercell and the symmetry mates.                       |
|sc111     |Make a lattice of 1 x 1 x 1 unit cells.                        |
|sc221     |Make a lattice of 2 x 2 x 1 unit cells.                        |
|sc112     |Make a lattice of 1 x 1 x 2 unit cells.                        |
|sc222     |Make a lattice of 2 x 2 x 2 unit cells.                        |
|sc331     |Make a lattice of 3 x 3 x 1 unit cells.                        |
|sc313     |Make a lattice of 3 x 1 x 3 unit cells.                        |
|sc133     |Make a lattice of 1 x 3 x 3 unit cells.                        |
|sc333     |Make a lattice of 3 x 3 x 3 unit cells.                        |

## Save files with date and time in filename

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|saln      |Save a aln file (alignment file) with a time stamp included in the filename to avoid overwriting work.. |
|scif      |Save a cif file (Crystallographic Information File) with a time stamp included in the filename to avoid overwriting work.. |
|sccp4     |Save a ccp4 file (CCP4 electron density map file) with a time stamp included in the filename to avoid overwriting work.. |
|sdae      |Save a dae file (Collada File) with a time stamp included in the filename to avoid overwriting work.. |
|sdat      |Save dat file (output data file) with a time stamp included in the filename to avoid overwriting work. |
|sfasta    |Save a fasta file (sequence file) with a time stamp included in the filename to avoid overwriting work.. |
|sidtf     |Save a idtf file (Intermediate Data Text Format) with a time stamp included in the filename to avoid overwriting work.. |
|smae      |Save mae file (Maestro file) with a time stamp included in the filename to avoid overwriting work. |
|smmd      |Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work. |
|smmod     |Save mmd file (Macromodel file) with a time stamp included in the filename to avoid overwriting work. |
|spmo      |Save pmo file (XYZ, binary format file) with a time stamp included in the filename to avoid overwriting work. |
|smoe      |Save moe file (Molecular Operating Environment) with a time stamp included in the filename to avoid overwriting work. |
|smol      |Save mol file with a time stamp included in the filename to avoid overwriting work. |
|smol2     |Save mol2 (Sybyl file format) file with a time stamp included in the filename to avoid overwriting work. |
|smtl      |Save mtl (Wavefront Material file format) file with a time stamp included in the filename to avoid overwriting work. |
|sobj      |Save obj file (Wavefront mesh file) with a time stamp included in the filename to avoid overwriting work. |
|sout      |Save out file (output data file) with a time stamp included in the filename to avoid overwriting work. |
|spdb      |Save pdb file with a time stamp included in the filename to avoid overwriting work.. |
|spkl      |Save a pkl file (Python pickle file) with a time stamp included in the filename to avoid overwriting work.. |
|spkla     |Save a pkla file (Python pickle file) with a time stamp included in the filename to avoid overwriting work.. |
|spng      |Save png file with a time stamp included in the filename to avoid overwriting work. |
|spov      |Save pov (POV-ray tracing file format) file with a time stamp included in the filename to avoid overwriting work. |
|spqr      |Save pqr file with a time stamp included in the filename to avoid overwriting work. |
|spse      |Read as a commandline argument, a string as the filename stem or  |
|ssdf      |Save sdf file with a time stamp included in the filename to avoid overwriting work. |
|swrl      |Save wrl (VRML 2 file format) file with a time stamp included in the filename to avoid overwriting work. |

## Make molecular representations that are not available in PyMOL

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|AO        |Commands to make ambient occlusion image like those in Qutemole.  |
|AOD       |Make ambient occlusion image of any with dark carbon atoms.    |
|BW        |Commands to make black-and white-ribbon cartoon on a white background. |
|BU        |Commands to make biological unit. Requires a pdb file. There are |
|CB        |Loads Jared Sampson's script "colorblindfriendly.py" from the  |
|CR        |Commands to make colored filled-ring cartoon of nucleic acids. May |
|CSS       |Commands to color ribbon or cartoon representations of proteins by |
|CBSS      |Apply colorblind-friendly coloring to ribbon or cartoon representations. |
|DU        |Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins.  |
|FR        |Make filled-ring cartoon of nucleic acids. May need to enter 'hide everything' first.  |
|HH        |Hide hydrogen atoms of currently visible molecular objects.    |
|PE        |Apply pearl effect about cations. Must supply selection.       |
|PU        |Make putty cartoon of main chain of nucleic acids and proteins. |
|SE        |Commands to make SAXS envelope from a bead model.              |
|getchem   |Create selections based on the biophysical properties of each residue. |
|timcolor  |Use Tim Mather's coloring scheme applied to the selections defined in getchem().  |

## Launch a full-featured text editor from PyMOL

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|atom      |Open file with the text editor Atom from within PyMOL.         |
|bbedit    |Open file with the text editor bbedit from within PyMOL.       |
|code      |Open file with Visual Studio Code from within PyMOL.           |
|emacs     |Open file with emacs from within PyMOL.                        |
|gedit     |Open file with gedit from within PyMOL.                        |
|jedit     |Open file with jedit from within PyMOL.                        |
|mate      |Open file with Textmate (Mac OS only) from within PyMOL.       |
|notepadpp |Open file with notepadpp (Mac OS only) from within PyMOL.      |
|nv        |Open file with neovim from within PyMOL.                       |
|oni       |Open the editor Oni from within PyMOL.                         |
|pdbed     |Open PDBEditor.jar from within PyMOL.                          |
|st3       |Open sublime text 3 from within PyMOL.                         |
|vim       |Open vim from within PyMOL.                                    |

## Open word processord

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|word      |Open word from within PyMOL.                                   |

## Open data analysis programs

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|cranR     |Open the Cran R from within PyMOL.                             |
|ddb       |Open the DBBrowserSQLite.                                      |
|excel     |Open the excel from within PyMOL.                              |
|JASP      |Open the JASP from within PyMOL.                               |
|JMP       |Open the JMP from within PyMOL.                                |
|jabref    |Open the jabref from within PyMOL.                             |
|julia     |Open the jabref from within PyMOL.                             |
|oc        |Open the jabref from within PyMOL.                             |
|ppt       |Open the powerpoint from within PyMOL.                         |

## Open terminal windows

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|iterm     |Open iTerm2 window on MacOS.                                   |
|term      |Open a Terminal window on MacOS.                               |

## Open other molecular graphics programs

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|ccp4mg    |Open ccp4mg from within PyMOL.                                 |
|chimera   |Open Chimera from within PyMOL.                                |
|coot      |Open coot from within PyMOL.                                   |
|jmol      |Open Jmol from within PyMOL.                                   |
|vmd       |Open vmd from within PyMOL.                                    |
|yasara    |Open the molecular graphics prograom YASASRA from within PyMOL.  |

## Image manipulation programs

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|gimp      |Open the molecular graphics program with gimp from within PyMOL.  |
|inkscape  |Open the molecular graphics program with gimp from within PyMOL.  |

## Open certain webapps

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|gcal      |Open Google Calendar.                                          |
|GM        |Open gmail.                                                    |
|WM        |Open Web Mail in defualt browser. Adjust url for your institution. |
|WS        |Open National Weather Service website for locale.              |

## Samples

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|GGT       |WT human gamma glutamyl transpeptidase at 1.67 Angstrom        |
|GU        |10-mer dsRNA with 8 contiguous Us. U-helix RNA.                |
|N9        |Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code 4dgr. |
|T4L       |WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0.            |
|U8        |16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. |
|WC8       |16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom            |

## Commands to display complex scenes.

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|BST       |G2G3/U9U8 base step , PDB code 4PCO.                           |
|LG        |Nine sugar glycan in influenza N9 neuraminidase at             |
|NA        |Hydrated sodium cation bound in major groove of a              |

## Commands to display complex scenes with pdb files on computer.

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|LGGT      |WT human gamma glutamyl transpeptidase at 1.67 Angstrom        |
|LGU       |10-mer dsRNA.                                                  |
|LN9       |Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code |
|LT4L      |Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0.  |
|LU8       |16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. |
|LWC8      |16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom            |
|LBST      |G2G3/U9U8 base step , PDB code 4PCO.                           |
|LLG       |Nine sugar glycan in influenza N9 neuraminidase at             |
|LNA       |Hydrated sodium cation bound in major groove of a              |

## Re-orient molecule

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|oy        |Align long axis of molecule along z-axis.                      |
|omxy      |Align long axis of molecule along minus x-y axis.              |
|oxy       |Align long axis of molecule along x-y axis.                    |
|oz        |Align long axis of molecule along y-axis.                      |

## Horizontal scripting

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|cntfiles  |Count number of files in current directory.                    |
|cntpdb    |Count number of pdb files in current directory.                |
|rline     |Enter "help(rline)" to refresh memory of the readline commands. |
|rv        |Get the view settings in a compact format on one line.         |

## Print commands for using git for version control

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|gitAdd    |Enter help(gitAdd) to print steps for adding a file to version control. |
|gitCommit |Enter help(gitInit) to print steps for saving updates to a file under version control. |
|gitInit   |Enter help(gitInit) to print steps for creating a git repository. |
|gitPull   |Enter help(gitPush) to print steps to send to updates to a repository on github.com.  |
|gitPush   |Enter help(gitPush) to print steps to send to updates to a repository on github.com.  |

## Send search term(s) to websites with search boxes.

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|AB        |Send search term or phrase to Amazon.com Books in default browser. |
|def       |#     Send search term or phrase to anaconda.com books in default browser. |
|GB        |Send search term or phrase to Amazon.com in default browser.   |
|GH        |Send search term or phrase to GitHub in default browser.       |
|GHN       |Send search term or phrase to GitHub in default browser.       |
|GO        |Send search term or phrase Google in default browser.          |
|GON       |Send search term or phrase Google in default browser and opens the top N results in N new tabs. |
|GS        |Send search term or phrase to Google Scholar in default browser. |
|GSN       |Send search term or phrase to Google Scholar in default browser. |
|GV        |Send search term or phrase to Google Videos in default browser. |
|GVN       |Send search term or phrase to Google Videos in default browser. |
|MA        |Send search term to all searchable websites in pymolshortcuts: |
|MB        |Send search term to search multiple sites for term in books:   |
|MC        |Send search term to search ten core websites in pymolshortcuts: |
|MM        |Send search term to search for manuscripts in pymolshortcuts:  |
|PDB       |Submit a search term to the Protein Data Bank.                 |
|PDBN      |Submit a search term to the Protein Data Bank and open the top N hits in separate tabs. |
|PML       |Submit a search term to the PyMOL Users Mail Service.          |
|PMLN      |Submit a search term to the PyMOL Users Mail Service.          |
|PM        |Send search term or phrase to PubMed.                          |
|PMN       |Send search term or phrase to PubMed and open top N hits in separate tabs. |
|IPM       |Read list of search terms and submit each term to PubMed in a separate browser tab. |
|IPMN      |Read list of search terms and submit each term to PubMed in a separate browser tab. |
|RG        |Submit a search query of Research Gate.                        |
|RGN       |Submit a search query of Research Gate and open the top N hits in sepearte tabs.  |
|SD        |Submit a search term to Science Direct.                        |
|SDN       |Submit a search term to Science Direct and open the top N hits in sepearte tabs. |
|SF        |Send search term to sourceforge.                               |
|SFN       |Send search term to sourceforge and open the top N hits in sepearte tabs. |
|SP        |Submit a search term to Springer Books                         |
|SPN       |Submit a search term to Springer Books and open the top N hits in sepearte tabs. |

## Open static web sites

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|ACA       |Open the American Crystallographic Association Annual Meeting webpage. |
|ALS       |Open website of the Advanced Light Source.                     |
|APS       |Open website of the Advanced Photon Source.                    |
|AX        |Send search term or phrase to arXiv.                           |
|BC        |Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source. |
|BD        |Open the webpage of the Small Angle Scattering Biological Data Bank (SASBDB).  |
|BX        |Send search term or phrase to bioRxiv                          |
|CH        |Open the webste of UCSF Chimera.                               |
|CHESS     |Open the website of CHESS.                                     |
|EMDB      |Open the website of the Electron Microscopy Data Bank.         |
|EP        |EasyPyMOL github site.                                         |
|JM        |Open the Jmol wiki.                                            |
|IUCR      |Open website of the IUCr Journals.                             |
|LBSF      |Open website of Laboratory of Biomolecular Structure and Function, the X-ray diffraction core facility at OUHSC. |
|MCL       |Open website of Macromolecular Crystallography Laboratory at the University of Oklahoma.  |
|MG        |Open website of the OUHSC molecular graphics course.           |
|NDB       |Open website of the Nucleic Acid Database.                     |
|notPyMOL  |Open website with list of other molecular graphics programs.   |
|NSLSII    |Open the website of the National Synchrotron Light Source II (NSLSII) at Brookhaven National Laboratory. |
|PPC       |Open the website of the Protein Production Facility at the University of Oklahoma in Norman. |
|PS        |Open the home page of the Protein Soceity.                     |
|PW        |Submit search of the PyMOL Wiki.                               |
|RS        |Open the homepage of the RNA Society.                          |
|SAXS      |Open the webpage of SAXS links at OUHSC.                       |
|SB        |Open the webpage of SSRL Biological SAXS at BL 4-2.            |
|SBGRID    |Open the webpage of the Structural Biology Grid (SBGRID) YouTube Channel. |
|SciPy18   |Open the SciPy 2018 YouTube Channel.                           |
|SSRL      |Open the webpage of SSRL Structural Molecular Biology.         |
|SSURF     |Open the webpage of the Society for Science at User Research Facilities (SSURF). |
|SO        |Submit a search term to Stackoverflow.                         |

## 3D-PDFs

|Shortcuts |Description                                                      |
|:--------|:---------------------------------------------------------------|
|ms2pdf    |Send molecular surface or ribbon cartoon from PyMOL to 3dpdf.  |
|topdf     |Send stick models as pse file from PyMOL through Jmol to 3DPDF. |