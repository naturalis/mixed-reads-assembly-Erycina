Mixed-reads-assembly
====================
Config files and scripts for assembly of mixed short (Illumina) and long (PacBio) reads

Description
-----------
What is this project trying to accomplish?
This project was started as a collaboration between Hogeschool Leiden students and Naturalis Biodiversity Center in order to try and make a first draft assembly of the *Erycina pusilla* genome using Illumina reads and PacBio contigs. In order to do this we have created a work environment. This work environment can be found on the OpenStack cloud server of Naturalis Biodiversity Center. 

Dependencies
------------
Which dependencies are used and how are they managed?
The work environment has been created on an Ubuntu operating system. This project utilizes several software packages, namely:

1. Burrows-Wheeler Aligner
2. SAMtools
3. Seqtk
4. ncbi-blast+

Requirements
------------
How much RAM/disk space/processor cores are needed for this work environment?
For this project we had been allocated 64 Gb RAM, 8 VCPU's and 160 Gb disk storage. For this project we made use of 6 of 8 cores. We had chosen to use 6 cores to speed up the assembly. It is possible to use 7 cores but we wanted to play it safe and allow the operating system to use 2 cores for its other functions. 

The RAM usage never came above 10 Gb so the allocated 64 Gb is more than enough.

The 160 Gb disk storage wasn't enough to store all of the initial data so we also used an additional volume (disk storage) of 1000 Gb which was more than enough to store our data and results and also makes the work environment somewhat future proof.

Architecture
------------
How do the work environment and the data volume(s) fit together?
The work environment we have created is located on the instance and our data and results are located on the volume. This is better explained in the figure below.

![alt text](https://github.com/naturalis/mixed-reads-assembly-Erycina/blob/master/doc/work_environment_layout.png "Work environment layout")

In this figure you can see the general layout of the work environment and it's connection with the volume. We used the doc folder to store all our small documents and papers. In the src folder we stored the code that we have written and in the bin folder we stored the script we used as initialization script for the OpenStack instance.

The broken lines display the connection of the instance and the volume. On the volume we have the data folder which houses all of the initial data and the results folder contains the results we have gathered during this project.

For more information please refer to the project report Draft_assembly_Erycina_pusilla_report.docx in the doc folder.
