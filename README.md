# DNASimulator
DNA Simulator for Technion

# [Please check the wiki for the documentation](https://github.com/NiliStein/DNASimulator/wiki)

# DNASimulator

DNASimulator is a tool developed by Gadi Chaykin and Nili Furman and supervised by Omer Sabary and Eithan Yaakovi from the Technion - Israel Institute of Technology, in order to provide a convenient simulator environment for testing methods of data encoding and decoding in DNA storage.

The GUI of the tool was developed with pyQt.

The simulator modules were tested by checking the outputs' results and statistics, and the testing scripts or parts of them can be found commented out in the end of the relevant .py files.

# Running The App

First, install **python 3.8.6**. This is the version best compatible with some of the other libraries we use.

Install the following:
* PyQt5
* scipy
* Pillow (PyPI)
* matplotlib
* edlib

Then, to run the app, run:

```
python app.py
```
or
```
python3 app.py
```
(depends on your environment)

# Project Hierarchy

The package consists of the following directories and files, in the following hierarchy:
* 📂 input
    * 📄 strands_in.txt
* 📂 reconstruction_algs
    * ▶️ BMALookahead.exe
    * ▶️ BMALookaheadLinux
    * ▶️ BMALookaheadMac
    * ▶️ DivBMA.exe
    * ▶️ DivBMALinux
    * ▶️ DivBMAMac
    * ▶️ Hybrid.exe
    * ▶️ HybridLinux
    * ▶️ HybridMac
    * ▶️ Iterative.exe
    * ▶️ IterativeLinux
    * ▶️ IterativeMac
    * 📚 libgcc_s_seh-1.dll
    * 📚 libstdc++-6.dll
    * 📚 libwinpthread-1.dll
* 📂 shuffle_prog
    * ▶️ shuf_windows.exe
    * ▶️ shuf_mac
    * 📚 msys-2.0.dll
    * 📚 msys-iconv-2.dll
    * 📚 msys-intl-8.dll
* 🐍 app.py
* 🐍 simulator.py
* 🐍 strand_error_sim.py
* 🐍 custom_random_variable.py
* 🐍 dnaSimulator_ui2.py
* 🐍 SpinBoxCustom.py
* 📄 dnaSimulator2.ui

In the **input** folder you can find an input file for example for the error simulation.

The **reconstruction_algs** folder holds all runnables of the reconstruction algorithms for each of the 3 platforms the app can run on: Windows, Linux and OSX. In addition, it also contains a few dlls that are required for the Windows executables to run on certain Windows environments.

The **shuffle_prog** folder contains the shuffling programs (and the dlls required to run in on Windows) that are used in the code, for the platforms that need them: Windows and OSX. (Linux has a shuffle command and therefore doesn't need a program)

The rest of the files are python code files of the app, and a .ui file compatible with Qt Creator.

In addition, the following will be generated **after running the errors simulator or the reconstruction**:
* 📂  output
    * 📄 evyat.txt
    * 📄 errors_shuffled.txt
    * 📄 histogram.txt
    * 🖼️ histogram.png
    * 📄 output.txt
    * 📄 output-results-fail.txt
    * 📄 output-results-success.txt

The output folder holds all outputs from the app: evyat.txt & errors_shuffled.txt from the error simulator, and the rest from the reconstruction. The clustering modifies the evyat.txt from this directory.
