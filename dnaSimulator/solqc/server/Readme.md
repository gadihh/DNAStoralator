# SOLQC Server
The easiest way to operate the SOLQC is by using the SOLQC server.<br>
The server uses a GUI to expose and operate the different functionalities of the tool.

# Contents
* [Setup](#setup)
* [Preparation (What you need)](#preparation)
* [Usage](#usage)
* [Output](#output)
* [Configuration File](#configuration)
* [Example using the toy data](#example)

## Setup
In it's current state we assume the user as some familiarity with python. <br>
You'll need to run the tool with python *3.6.5*.<br>
Start by cloning the repository to a local directory.<br>
Next, open a command line tool, go to the root folder and run:<br>
`pip install -r requirements.txt`* <br>
This will install all the necessary modules to run the tool.<br>
*If you are familiar with python virtual enviorments, it is advised
to run this command in a new virtual env.

## Preparation
In order to use the tool you'll need the following:
* Design, could one of 2 options:
  * A design file, in a csv format containing 2 columns : [barcode, variant]
    * barcode - a sequence identifier for the variant. [Needed for matching between a read and a variant].
    * variant - the complete variant sequence. [Needed for the alignment to analyse missmatches and indel's.
  * IUPAC string
* A reads text file containing all the fasta/q files names of the sequenced read (one row for each file).
* A config.json file containing different possible configuration, see - [configuration](#configuration)

Here is an example for each of those files:
* [design.csv](https://github.com/yoavo1984/SOLQC/blob/master/data/toy_data/design.csv)

  <img src="img/desing.png" alt="drawing" width="224" height="224"/>

* [reads.txt](https://github.com/yoavo1984/SOLQC/blob/master/data/toy_data/reads.txt)
```txt
data/my_data/reads_1.fastq
data/my_data/reads_2.fastq
```
* [config.json](https://github.com/yoavo1984/SOLQC/blob/master/data/toy_data/config.json)
```json
{
    "prefix" : "ACAACGCTTTCTGTGTCGTG",
    "suffix" : "",
    "length" : 0,
    "barcode_start" : 20,
    "barcode_end" : 32,
}
```

## Usage
Open a command line and to go the root folder and run one of the following:<br>
`$ python qc_server.py` [You will navigate to the folder contatining the data using a GUI]<br>
or <br>
`$ python qc_server.py <path_to_folder>` <br>
Where path to folder is a the path to a folder containing the design.csv, reads.txt and config.json.<br>
If all went well this should open a web browser looking like this : 

<img src="https://raw.githubusercontent.com/yoavo1984/SOLQC/master/img/web_interface.png" alt="web" width="812" height="500"/>

### Opeartion
Now, choose your design, reads and config file. <br>
Next, choose your method of matching and what type of analysis you want to preform. 

### Matching
#### Barcdoe Aligner 
A read will be matched to a variant only if they have the exact same sequence in the barcode area (defined in the *config.json* file.
#### N Edit barcode aligner
A read will be matched to a variant only if the edit distance between the variant barcode and the read barcode area is smaller than a given tolerance(in the config.json file)
