# Fragariyo: Hunting for Protein Fragments
###### Update: May 13th, 2024

## What is Top-Down Mass Spectrometry?
In Top-Down Mass Spectrometry protein or protein complexes are fragmented intact inside a mass spectrometer. Therefore, if there are specific combination of Post-Translationla Modifications or amino acid substitutions fragmentation will be obtined for each one (unfrotunatly, in bottom the resultatn digested mix of peptides precludes the determination of an specific protein from - a.k.a proteoform). For more details, see the publications below.

  - 2016 [Top-Down Review](https://www.annualreviews.org/content/journals/10.1146/annurev-anchem-071015-041550)
  - [Top-Down Proteomics and the Human Proteoform Project](https://www.science.org/doi/full/10.1126/sciadv.abk0734)
  
### Package Installation

#### 1. Download .exe file from [the release page](https://github.com/RuotoloLab/Fragariyo/releases) for the latests versions 
#### 2. Follow the installation prompts
#### 3. Press on icon and a Graphical User Interface (GUI) will appeared. 

### Terminal Fragments

###### Note: If the data does not have good mass accuracy, please internally calibrate the dataset before searching....otherwise internal fragment false discovery rate will reach the moon!

#### 1. On the GUI select **Change Output Directory**
#### 2. Load Modifications Library **(template name: ModificationsRepository.txt)**
#### 2.5 (optional Waters IM-MS files) **Run IMTBX output extraction and renaming**
#### 3. Select **Terminal Run**. Formats accepted: .isotopes files, mMass files, CSV files (m/z, int), as well as unmatched files (produced by Fragariyo after running a terminal fragment analysis). **(template file: NISTmAb_TerminalFragments_template.csv)**
#### 3.5 There is a dropdown menu with sequence coverage graphing capabilities **Data Analysis Options**.

### 4. Internal Fragments 
SEARCH TERMINAL FRAGMENTS FIRST!!! Searching for internal fragments  first it is found to turnup a lot of false positives as the search space can be huge. Unless you have removed terminal fragmetns first by Fragariyo or other methods or somehow your fragmentation expeirment produced 100% internal fragmetns *gasp*), search terminal fragment first. 

#### 4.1 Running **Input Generator**. It requires unmatched files (left over ions from the terminal search) or ions in CSV format, and also requires xy files (m/x, int) files to obtain isotope information for each experimental ions. 
#### 4.2.1 Select Mass Ion Resolution, required for proper isotope theoretical production
#### 4.2.2 Select Mass Error (recommend to be 1 ppm). 
#### 4.2.3 Select **Internal fragment**. A batch file will be required (template: Internalfrag_batchmode_NISTmAb.csv), where files and param template files (template: NISTmAb_InternalFragments_template.csv) can be run in a high throughput manner.










Authors:
Carolina Rojas Ramirez 
(author of the Disulfinator, reconstruction of the IMAnnotator into the Fragmentor and new data analysis tools)

Dan Polasky 
(author of the initial version of some data analysis tools and of the first draft of IMAnnotator)

##### Fragariyo is a collection of Python 3 scripts to enable peak matching 
and some data analysis tools for native top-down experiments. See the 
SOP for more information.


Upon isntallation the user can find template files inside the main Fragariyo folder under the C drive.
