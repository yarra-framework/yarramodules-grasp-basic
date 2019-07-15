# YarraModules-GRASP-Basic

Yarra integrated Matlab implementation of the basic GRASP reconstruction technique. Developed by Li Feng, Ricardo Otazo, and Kai Tobias Block.

This module uses several open-source Matlab packages developed by other contributors, including the mapVBVD library by Philip Ehses, the NUFFT library developed by Jeff Fessler, and the inifile reader developed by Primoz Cermelj. Please give credit to these developers when using the module.


##Reference

This method implements the "GRASP" reconstruction for free-breathing DCE-MRI acquisitions, which has been described in the following publication:

Feng L, Grimm R, Block KT, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. 
Golden-angle radial sparse parallel MRI: combination of compressed sensing, parallel imaging, and golden-angle radial sampling for fast and flexible dynamic volumetric MRI. 
Magn Reson Med. 2014 Sep;72(3):707-17.

https://www.ncbi.nlm.nih.gov/pubmed/24142845

Please cite this publication when using the module for your work.


##License Information
The Yarra framework is provided free of charge for use in research applications. It must not be used for diagnostic applications. The author takes no responsibility of any kind for the accuracy or integrity of the created data sets. No data created using the framework should be used to make diagnostic decisions.

The Yarra framework, including all of its software components, comes without any warranties. Operation of the software is solely on the user's own risk. The author takes no responsibility for damage or misfunction of any kind caused by the software, including permanent damages to MR systems, temporary unavailability of MR systems, and loss or corruption of research or patient data. This applies to all software components included in the Yarra software package.

The source is released under the GNU General Public License, GPL (http://www.gnu.org/copyleft/gpl.html). The source code is provided without warranties of any kind.