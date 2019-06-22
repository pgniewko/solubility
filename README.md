## Solubility Challange

### Background
Intrinsic solubility (water solubility): solubility of non-charged molecules, i.e. free acid and base free form. It is required that the solubility is determined in the presence of solid substance.     

Put here an example (table) where I show how solubility changes with pH.     

The project is motivated by the [challange](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00345). 

### Potential sources of bulk data
1. [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
2. [OCHEM](https://ochem.eu/home/show.do)

### This work
To be written.

### Results

For the comaprison I provide data from [PotentialNet paper](https://pubs.acs.org/doi/full/10.1021/acscentsci.8b00507)

| Network       | Valid RMSE | Test RMSE |
| ------------- |:----------:| ---------:|
| PotentialNet  | 0.517      |  0.490    |
| Weave         | 0.549      |  0.553    |
| GraphConv     | 0.721      |  0.648    |
| XGBoost       | 1.182      |  0.912    |


[Datasets](#datasets)    
[Papers](#papers)    

### Datasets

| Dataset                   | Do I trust it? | Comments                                 |
|---------------------------|:--------------:|:-----------------------------------------|
| AB.2001.EJPS              | (+/-)          | Given in weird unists log(1/S0)[mol/L]?  |
| ABB.2000.PR               | (+/-)          | Given in weird unists log(1/S0)[mol/L]?  |
| BOM.2017.JC               | (+)            |                                          |
| D.2008.JCIC               | (+)            |                                          |
| H.2000.test1              | (+)            | Downloaded from the website              |
| H.2000.test2              | (+)            | Downloaded from the website              |
| H.2000.train              | (+)            | Downloaded from the website              |
| HXZ.2004.JCIC.data_set    | (+)            | Downloaded from the website              |
| HXZ.2004.JCIC.test_set1   | (+)            | Downloaded from the website              |
| LGG.2008.JCIM.32          | (+)            |                                          |
| LGG.2008.JCIM.100         | (+)            |                                          |
| LPB.2013.JCIC [all]       | (-)            | Can't understand the format of the data! |
| POG.2007.JCIM.test        | (+)            | Data obtained from authors               |
| POG.2007.JCIM.train       | (+)            | Data obtained from authors               |
| WKH.2007.JCIM.solubility  | (+)            | ADME website data                        |    
| OCHEM.WaterSolubility     | (+/-)          | Lots of repeats, some sign error         |
| PubChem                   | (+/-)          | No logS0 data, Measurements at pH=7.4    |

### Papers          


1. **Can You Predict Solubilities of Thirty-Two Molecules Using a Database of One Hundred Reliable Measurements?**     
Antonio Llinàs, Robert C. Glen and Jonathan M. Goodman    
*J. Chem. Inf. Modeling 2008, 48, 1289-1303*      
[[paper 1]](https://pubs.acs.org/doi/10.1021/ci800058v)      
[[paper 2]](https://pubs.acs.org/doi/10.1021/ci800436c)    
[[website]](http://www-jmg.ch.cam.ac.uk/data/solubility/)      
**Note 0: This is the reference for the original Solubility Challange**      
**Note 1: In the test set, SMILES strings for probenecid and pseudoephedrine were swapped. Use only soldataswap.xls file.**      
**Note 2: Solubility for 32 compounds taken from HEL.2009.JCIM.pdf**       
**Note 3: Data was downloaded from the original website, but the numbers are dubious (IMO) - use CAREFULLY!**      


2. **ESOL: Estimating Aqueous Solubility Directly from Molecular Structure**    
John S. Delaney           
*J. Chem. Inf. Comput. Sci. 2004, 44, 1000-1005*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci034243x)    
**Note: There are two files D.2008.JCIC.solubility.v[1-2].txt. These files are the same but come from two 
different sources: (i) [Pat Walters Blog](https://github.com/PatWalters/solubility) (ii) [ChemDB](ftp://ftp.ics.uci.edu/pub/baldig/learning/Delaney)**      

3. **Can You Predict Solubilities of Thirty-Two Molecules Using a Database of One Hundred Reliable Measurements?**     
Jarmo Huuskonen
*J. Chem. Inf. Comput. Sci. 2000, 40, 773-777*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci9901338)     
[[website]](http://cheminformatics.org/datasets/huuskonen/index.html)     
**Note: Quite a few repeats from Delaney Set. Different measurements though.**     


4. **ADME evaluation in drug discovery. 4. Prediction of aqueous solubility based on atom contribution approach**
Tingjun Hou, Ke Xia, Wei Zhang, Xiaojie Xu     
*Journal of Chemical Information and Computer Sciences, 2004, 44, 266-275*       
[[paper]](https://pubs.acs.org/doi/abs/10.1021/ci034184n)     
[[website]](http://modem.ucsd.edu/adme/databases/databases_logS.htm)      


5. **Development of reliable aqueous solubility models and their application in drug-like analysis**     
Junmei Wang, George Krudy, Tingjun Hou, George Holland, Xiaojie Xu      
*Journal of Chemical Information and Modeling, 2007, 47, 1395-1404*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci700096r)     
[[website]](http://modem.ucsd.edu/adme/databases/databases_logS.htm)     
**Note: In logS database, the aqueous solubility was expressed as logS, where S is the solubility at a temperature of 20-25°C in mol/L. These are two databases for our modeling. In reference [4], the data afforded by Tetko was used. This database includes 1290 organic compounds. The data set was converted from the SMILES flat file representation to the MACCS/sdf structured data file. In reference [5], some new molecules collected from literature were added. This database includes 1708 molecules.**


6. **Can human experts predict solubility better than computers?**      
Samuel Boobier, Anne Osbourn and John B. O. Mitchell    
*Journal of Cheminformatics, 2017, 9:63*       
[[paper]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0250-y)      
[[website]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0250-y)     
**Note: Source codes accompany the paper.**      


7. **pH-metric solubility. 3. Dissolution titration template method for solubility determination**      
Alex Avdeef, Cynthia M. Berger     
*European Journal of Pharmaceutical Sciences 14 (2001) 281–29*       
[[paper]](https://www.ncbi.nlm.nih.gov/pubmed/11684402)      


8. **pH-Metric Solubility. 2: Correlation Between the Acid-Base Titration and the Saturation Shake-Flask Solubility-pH Methods**    
Alex Avdeef, Cynthia M. Berger, and Charles Brownell         
*Pharmaceutical Research, Vol. 17, No. 1, 2000*      
[[paper]](https://www.ncbi.nlm.nih.gov/pubmed/10714613)      


9. **Online chemical modeling environment (OCHEM): web platform for data storage, model development and publishing of chemical information**
Iurii Sushko et al.,     
*J Comput Aided Mol Des (2011) 25:533–554*      
[[paper]](https://www.ncbi.nlm.nih.gov/pubmed/21660515)      
[[server]](https://ochem.eu/home/show.do)      
**Note: Reference for OCHEM. Not yet sure how to download the data :-( **    


10. **Solubility Challenge revisited after 10 years, with multi-lab shake- flask data, using tight (SD 0.17 log) and loose (SD 0.62 log) test sets**   
Antonio Llinas, and Alex Avdeef       
*J. Chem. Inf. Model., 2019*
[[paper]](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00345)     
**Note: The reference for the new challange.**      


11. **Random Forest Models To Predict Aqueous Solubility**      
David S. Palmer, Noel M. O’Boyle, Robert C. Glen, and John B. O. Mitchell       
*J. Chem. Inf. Model. 2007,471, 150-158*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci060164k)    
**Note: Data extracted from pdfs**     


12. **Deep Architectures and Deep Learning in Chemoinformatics**      
Alessandro Lusci, Gianluca Pollastri, and Pierre Baldi     
*J. Chem. Inf. Model. 2013,537, 1563-1575*       
[[paper]](https://pubs.acs.org/doi/abs/10.1021/ci400187y)      
**Note: Some of the files/data are duplicates**       


13. **Is Experimental Data Quality the Limiting Factor in Predicting the
Aqueous Solubility of Druglike Molecules?**    
David S. Palmer and John B. O. Mitchell      
*Mol. Pharmaceutics 2014, 11, 2962−2972*      
[[paper]](https://pubs.acs.org/doi/10.1021/mp500103r)      
**Note:  Good overview of the sources of the errors in solubility prediction.**      


14. **Convolutional Networks on Graphs for Learning Molecular Fingerprints**    
David Duvenaud, Dougal Maclaurin, Jorge Aguilera-Iparraguirre, Rafael Gómez-Bombarelli, Timothy Hirzel, Alán Aspuru-Guzik, and Ryan P. Adams.      
*arXiv, 2015*      
[[paper]](https://arxiv.org/pdf/1509.09292.pdf)      
[[code]](https://github.com/HIPS/neural-fingerprint)       
**Note1: Original code in Python 2. In order to make it work use futurize to convert to Python 3**      
**Note2: install with `python setup.py install`**      

15. **Random Forest: A Classification and Regression Tool for Compound Classification and QSAR Modeling**      
Vladimir Svetnik, Andy Liaw, Christopher Tong, J. Christopher Culberson, Robert P. Sheridan, and Bradley P. Feuston     
*J. Chem. Inf. Comput. Sci. 2003, 43, 1947-1958*       
[[paper]](https://pubs.acs.org/doi/10.1021/ci034160g)      


16. **Binary Classification of Aqueous Solubility Using Support Vector Machines with Reduction and Recombination Feature Selection**            
Cheng, T., Li, Q., Wang, Y., and Bryant, S.H.    
*Journal of Chemical Information and Modeling, 2011, 51, 229-236*
[[paper]](https://pubs.acs.org/doi/10.1021/ci100364a)    
**Note: The measurements come from BioAssay AID:1996, and are done at pH=7.4. Not very useful for logS0 estimations.**              







