# pymolshortcuts

This repository for ***pymolschortucts.py*** which contains 131 functions mapped to short-names that work like aliases (the table below is incomplete).
These shortcuts include many convenience functions that make work in PyMOL more productive and fun!
Some shortcuts save the many hours of work required to assemble a new script file while other shortcuts save only a few minutes but lower motivational barriers.
You do not need to understand Python to be able run the script.
To have the shortcuts always available, add the command 'run ~/pymolshortcuts.py' to your pymolrc or pymolrc.pml file in your home directory to load the functions in pymolshortcuts.py on startup of PyMOL.
In spite of the command *run*, the functions will be loaded into memory but not executed.
Enter **SC** to get a list of shortcuts printed to the command history window.
Enter **help PW** to print the documentation for the shortcut **PW** to the command history window.
While most shortcuts are ready to use, some shortcuts require pasting of path names to the executables on your computer.

For example, the **BW** command makes cartoon or sphere representations of the currently loaded molecule into black and white figures suitable for classroom handouts and coloring books for aspiring scientists.

The shortcut **PW** takes one or more search terms and the sends them to the PyMOL Wiki.
A browser tab opens for each search term, so multiple searches are run in parallel while you continue your work in PyMOL.
You can inspect the results of the searches when there is a natural break in your work in PyMOL.
Other search functions can submit parallel searches of PubMed, Google, bioRxiv, Research Gate, GitHub, and more.
See the table below.

Another class of shortcuts launch your favorite full-featured text editor from within PyMOL.
You have to install the text editor and edit the file path to the executable.

Another class of shortcuts opens the manuscript or grant application that you are working on in Overleaf.
You have to edit the script by pasting in the appropriate link to the specific document.

Another class of shortcuts saves files with timestamps embedded in the filename to avoid overwriting png, pdb, pse, and other types of files written out from PyMOL.
These save function names begin with **s**, (e.g., **spse filename** saves the current session with date and time to nearest second embedded in the file name).
You can delete the unwanted version(s) at a latter time. 
These functions are useful if you do not have these files under version control.

The functions are stored in a single file: pymolshortcuts.py

Edit the file paths to outside executables like text editors.
The shortcuts work best on the command line below the command history window because you can copy and paste text onto this command line.
 
Videos that demonstrate representatives from each class of shortcut are planned.

## List available shortcuts

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|       SC|Print to screen list of the shortcuts that are available in the script pymolshortcuts.py.   |

## Show many models (NMR and crystal packing)

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|      nmr|                         Show all of the models in nmr structure.  |
|   nmroff|                     Hide all but first model in a nmr structure.  |
|     rmsc|                          Remove supercell and the symmetry mates. |
|    sc111|                          Make a lattice of 1 x 1 x 1 unit cells.  |
|    sc221|                          Make a lattice of 2 x 2 x 1 unit cells.  |
|    sc112|                          Make a lattice of 1 x 1 x 2 unit cells.  |
|    sc222|                          Make a lattice of 2 x 2 x 2 unit cells.  |
|    sc333|                          Make a lattice of 3 x 3 x 3 unit cells.  |


## Save files with date and time in filename

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|     spdb| Save pdb file with a time stamp in the filename to avoid overwriting previous work. |
|     spng| Save png file with a time stamp in the filename to avoid overwriting previous work. |
|     spse| Save pse file with a time stamp in the filename to avoid overwriting previous work.  |

## Make molecular representations that are not available in PyMOL

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|       AO| Commands to make ambient occlusion image like those in Qutemole.  |
|      AOD|      Make ambient occlusion image of any with dark carbon atoms.  |
|       BW|Commands to make black-and white-ribbon cartoon on a white background. |
|       BU|  Commands to make biological unit. Requires a pdb file. There are |
|       CB|     Loads Jared Sampson's script "colorblindfriendly.py" from the |
|     CBSS|Apply colorblind-friendly coloring to ribbon or cartoon representations. |
|       CR|Commands to make colored filled-ring cartoon of nucleic acids. May |
|      CSS|Commands to color ribbon or cartoon representations of proteins by |
|       DU|Make dumbbell (ribbons with rolled edges) cartoon of the main chains of nucleic acids and proteins.  |
|       FR|Make filled-ring cartoon of nucleic acids. May need to enter 'hide everything' first.  |
|       HH|      Hide hydrogen atoms of currently visible molecular objects.  |
|       PU|   Make putty cartoon of main chain of nucleic acids and proteins. |
|       SE|                 Commands to make SAXS envelope from a bead model. |
|  getchem|Create selections based on the biophysical properties of each residue. |
| timcolor|Use Tim Mather's coloring scheme applied to the selections defined in getchem().  |

## Launch a full-featured text editor from PyMOL

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|     atom|           Open file with the text editor Atom from within PyMOL.  |
|     code|             Open file with Visual Studio Code from within PyMOL.  |
|    emacs|                          Open file with emacs from within PyMOL.  |
|    gedit|                          Open file with gedit from within PyMOL.  |
|    jedit|                          Open file with jedit from within PyMOL.  |
|     mate|         Open file with Textmate (Mac OS only) from within PyMOL.  |
|       nv|                         Open file with neovim from within PyMOL.  |
|      oni|                           Open the editor Oni from within PyMOL.  |
|    pdbed|                            Open PDBEditor.jar from within PyMOL.  |
|      st3|                           Open sublime text 3 from within PyMOL.  |
|      vim|                         Open file with vim from within PyMOL.  |


## Submit search terms using default web browser

Send multiple searches with multiple commands on one line separated by semicolons.

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|       AX|                              Send search term or phrase to arXiv. |
|       BX|                            Send search term or phrase to bioRxiv  |
|       GO|             Send search term or phrase Google in default browser. |
|       GS|  Send search term or phrase to Google Scholar in default browser. |
|      PDB|                Submit a term for searching the Protein Data Bank. |
|       PM|                             Send search term or phrase to PubMed. |
|      PML|             Submit a search term to the PyMOL Users Mail Service. |
|       PW|                            Submit search term to the PyMOL Wiki.  |
|       RG|                         Submit a search query to Research Gate.  |


## Open static websites with default web browser

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|      ACA|Open the American Crystallographic Association Annual Meeting webpage. |
|      ALS|                        Open website of the Advanced Light Source. |
|      APS|                       Open website of the Advanced Photon Source. |
|       BC|Open the webpage of the BIOCAT biological SAXS beamline at the Advanced Photon Source. |
|       BD|Open the webpage of the  Small Angle Scattering Biological Data Bank (SASBDB).|
|       CH|                                  Open the website of UCSF Chimera. |
|    CHESS|                                       Open the website of CHESS.  |
|     EMDB|            Open the website of the Electron Microscopy Data Bank. |
|       EP|                                  Open the  EasyPyMOL github site. |
|       GM|                                                      Open gmail.  |
|       JM|                                               Open the Jmol wiki. |
|     IUCR|                                Open website of the IUCr Journals. |
|     LBSF|Open website of the Laboratory of Biomolecular Structure and Function (LBSF), the X-ray diffraction facility at OUHSC. |
|      MCL|    Open website of Macromolecular Crystallography Laboratory at OU |
|       MG|               Open website of the OUHSC molecular graphics course |
|      NDB|                        Open website of the Nucleic Acid Database. |
| notPyMOL|      Open website with list of other molecular graphics programs. |
|   NSLSII|Open the website of the National Synchrotron Light Source (NSLSII) at Brookhaven National Laboratory. |
|      PPC|        Open the website of the Protein Production Facility at OU. |
|       PS|                       Open the home page of the Protein Soceity.  |
|       RS|                             Open the homepage of the RNA Society. |
|     SAXS|                         Open the webpage of SAXS links at OUHSC.  |
|       SB|               Open the webpage of SSRL Biological SAXS at BL 4-2. |
|   SBGRID|                   Open the webpage of the SBGRID YouTube Channel. |
|  SciPy18|               Open the webpage of the SciPy 2018 YouTube Channel. |
|     SSRL|            Open the webpage of SSRL Structural Molecular Biology. |
|    SSURF|Open the webpage of the Society for Science at User Research Facilities. |
|      VSC|Open Visual Studio Code Market to obtain extensions. Use Safari on the Mac. |
|       WM|Open Web Mail in default browser. Adjust url for your institution. |
|       WS|Open National Weather Service website. Adjust url for your locale. |


## Launch other molecular graphics programs from inside PyMOL

Edit file paths to your executable files.

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|  chimera|                         Open Chimera from within PyMOL.  |
|     jmol|                            Open Jmol from within PyMOL.  |
|      vmd|                             Open vmd from within PyMOL.  |

## Molecules in standard orientation

Require internet connection to retrieve the coordinate and map files.

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|      GGT|           WT human gamma glutamyl transpeptidase at 1.67 Angstrom |
|       GU|                  10-mer dsRNA with 8 contiguous Us. U-helix RNA.  |
|       N9|  Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code |
|      T4L|              WT T4 lysozyme as ribbon diagram (1.08 Ang):  3FA0.  |
|       U8| 16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. |
|      WC8|              16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom  |

## Molecules in standard orientation 

Require local copies of pdb and map files.

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|     LGGT|           WT human gamma glutamyl transpeptidase at 1.67 Angstrom |
|      LGU|                                                    10-mer dsRNA.  |
|      LN9|  Influenza N9 neuraminidase at 1.55 Angstrom resolution, PDB code |
|     LT4L|Display WT T4 lysozyme as ribbon diagram (resoluton 1.08 Ang):  3FA0.  |
|      LU8| 16-mer dsRNA with 8 contiguous Us. U-helix RNA (1.37 Ang):  3nd3. |
|     LWC8|              16-mer dsRNA, Watson-Crick helix RNA. 1.55 Angstrom  |

## Complex figures, require internet connection

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|      BST|                             G2G3/U9U8 base step , PDB code 4PCO.  |
|       LG|               Nine sugar glycan in influenza N9 neuraminidase at  |
|       NA|                Hydrated sodium cation bound in major groove of a  |


## Complex figures, internet free (local copies of pdb files required)

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|     LBST|                             G2G3/U9U8 base step , PDB code 4PCO.  |
|      LLG|               Nine sugar glycan in influenza N9 neuraminidase at  |
|      LNA|                Hydrated sodium cation bound in major groove of a  |


## PyMOL to 3D PDFs

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|   ms2pdf|    Send molecular surface or ribbon cartoon from PyMOL to 3D-PDF.  |
|    topdf|   Send stick models as pse file from PyMOL through Jmol to 3D-PDF. |


## Horizontal scripting related

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|    rline|   Enter "help(rline)" to refresh memory of the readline commands. |
|       rv|   Get the view settings in a compact format on one line suitable for use on the command line. |


## File backups

| shortcut | Description |
|:--------|:---------------------------------------------------------------|
|  gitInit| Enter **help(gitInit)** to print the steps for creating a git repository. |
