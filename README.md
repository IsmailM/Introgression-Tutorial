# Introgression in cichlid fishes

######A tutorial on genomic tests for introgression

#### Summary

This tutorial demonstrates how phylogenies based on whole-genome sequence data can be used to detect introgression between closely related species. The data set used in this tutorial comes from a recent study on five species of Princess cichlid fishes (Lamprologini) from Lake Tanganyika, and the methodology used here reflects the one used in this study. In brief, the analysis will include the use of a hidden Markov model to identify alignment regions that support the same topology, Bayesian phylogenetic analyses of alignment blocks within these regions, and the application of a maximum likelihood framework to detect introgression based on a set of local phylogenies.

#### Background

Hybridization between closely related species, followed by back-crossing with the parental species, can lead to transfer of genetic material between established species, so-called introgression (e.g. [Mallet et al. 2016](http://onlinelibrary.wiley.com/doi/10.1002/bies.201500149/abstract)). This genetic transfer between species is an important evolutionary process that has facilitated adaptation in several species. For example, genes responsible for mimicry-related wing patterns were transferred between species of *Heliconius* butterflies ([The Heliconius Genome Consortium 2012](http://www.nature.com/nature/journal/v487/n7405/full/nature11041.html)), and Tibetans owe their altitude adaptations to the archaic human lineage of Denisovans ([Huerta-Sañchez et al. 2014](Altitude adaptation in Tibetans caused by introgression of Denisovan-like DNA)).

If the transferred genetic material becomes fixed in the recipient species, the true phylogeny of the genomic regions affected by introgression may differ from the species tree. Thus, by identifying regions in which phylogenies differ from each other, it is in principle possible to identify past hybridization that led to introgression. Unfortunately, the inference is complicated by several factors: First, the true species tree is usually unknown and will need to be inferred jointly with the introgression events. Second, even genomic regions not affected by introgression often differ in their phylogenetic histories due to incomplete lineage sorting (ILS; see e.g. [Suh et al. 2016](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002224)). Third, these phylogenetic histories of individual genomic regions are not known and can only be estimated based on sequence data. And fourth, the boundaries of these genomic regions characterized by different phylogenetic histories are unknown, and bias is introduced to phylogenetic inference if the alignment block used for the inference contains such boundaries ([Edwards et al. 2016](http://www.sciencedirect.com/science/article/pii/S1055790315003309)). In addition, the latter two issues are confounded, as the reliability of phylogenetic estimates increases with the amount of data used for inference, however, longer alignments are more likely to include regions with multiple phylogenetic histories and may therefore be more biased.

To account for all these complications, the methodology used here first aims to identify boundaries between genomic regions characterized by different phylogenetic histories. This is done using a hidden Markov model, implemented in the software Saguaro ([Zamani et al. 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)). Once these boundaries are identified, a large number of long alignment blocks are selected for reliable phylogenetic inference so that none of these blocks crosses any boundaries. Time-calibrated phylogenies are produced for each of these alignment blocks, using the Bayesian phylogenetic software BEAST 2 ([Bouckaert et al. 2014](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003537)). Finally, the inferred set of time-calibrated phylogenies is used as input for the software PhyloNet ([Than et al. 2008](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-322)), and statistical support for hypotheses of introgression is assessed with the maximum likelihood framework implemented in PhyloNet's CalGTProb function ([Yu et al. 2014](http://www.pnas.org/content/111/46/16448.abstract)).

In this tutorial, the above methods will be used to identify introgression between five closely related species of Princess cichlid fishes (Lamprologini) of Lake Tanganyika. Species from this group have previously been suggested to hybridize, based on mitochondrial sequences ([Salzburger et al. 2002](http://onlinelibrary.wiley.com/doi/10.1046/j.0962-1083.2001.01438.x/abstract)) and AFLP data ([Sturmbauer et al. 2010](http://www.sciencedirect.com/science/article/pii/S1055790310002897)). The sequence alignment used here is based on Illumina sequence data (with around 20x coverage) for the four lamprologine species *Neolamprologus gracilis*, *N. marunguensis*, *N. olivaceous*, and *N. pulcher*, mapped against the [BROAD institute's version 1.1 of the genome of *Oreochromis niloticus* (tilapia)](https://www.broadinstitute.org/ftp/pub/assemblies/fish/tilapia/), which was used as an outgroup. In addition, a fifth species of Lake Tanganyika Lamprologini, *N. brichardi*, as well as a second outgroup species from Lake Malawi, *Metriaclima zebra*, were included by also mapping available sequence data for these species against the genome of *Oreochromis niloticus* ([Brawand et al. 2014](http://www.nature.com/nature/journal/v513/n7518/full/nature13726.html)). Per species, consensus sequences from all reads were produced with BCFtools ([Li 2011](http://bioinformatics.oxfordjournals.org/content/27/21/2987.abstract)), VCFtools ([Danecek et al. 2011](http://bioinformatics.oxfordjournals.org/content/27/15/2156)), and Seqtk ([https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)). The resulting chromosome-length alignments thus included sequence data for seven species, including two outgroups (*O. niloticus* and *M. zebra*) and five ingroup species (*N. gracilis*, *N. marunguensis*, *N. olivaceous*, *N. pulcher*, and *N. brichardi*).

#### Requirements for this tutorial

The following software needs to be installed on your machine in order to run all analyses included in this tutorial:

* **Saguaro:** If Saguaro is not yet installed on your system, please follow the installation instructions at [http://saguarogw.sourceforge.net](http://saguarogw.sourceforge.net). Note that Saguaro runs on Linux, and while in principle the installation should also be possible on a Mac, I have been unable to get it to run on Mac OS X Yosemite. If you have access to a Linux server with a Saguaro installation, but you would like to run the rest of the tutorial on your own machine, you can do so by transferring input and ouput files of Saguaro via scp between your machine and the Linux server. To check whether the Saguaro installation succeeded, just type `Saguaro` in a console window.

* **Ruby:** This tutorial comes with a number of scripts for file conversion and visualization, and these are written in the programming language Ruby. Please make sure that a recent version of Ruby (2.0 and newer) is installed on your machine, which you can check by typing `ruby -v` in a console window. On most Linux and Mac systems, this should already be the case, and Ruby can be installed on all systems if it isn't. Dowloads and installation instructions are available at [https://www.ruby-lang.org/en/](https://www.ruby-lang.org/en/).

* **Firefox or another web browser:** For the visualization of Saguaro output, one of the Ruby scripts produces vector graphic figures in SVG format. These figures can be opened with a range of tools, including [Adobe Illustrator](http://www.adobe.com/products/illustrator.html) and [Inkscape](https://inkscape.org/en/), but any modern web browser ([Firefox](https://www.mozilla.org/en-US/firefox/new/), [Safari](http://www.apple.com/safari/), [Chrome](https://www.google.com/chrome/)) will also do the job.

* **AliView:** To visualize sequence alignments, the software AliView is recommended. The installation of AliView is described at [http://www.ormbunkar.se/aliview/](http://www.ormbunkar.se/aliview/) and should be possible on all operating systems. While AliView has many options for alignment editing, it is here used only for visualization, and can be replaced by other alignment viewers (or even a text editor), if for some reason the installation fails. A more extensive tutorial on AliView is available from the [Workshop on Molecular Evolution](http://evomics.org/learning/bioinformatics/multiple-sequence-alignment-ami-version/).

* **Java SDK 8:** In order to run BEAST, version 8 of the Java SE Development Kit is required. To check whether your machine already has this version, type `java -version` in a console window. If you then see something like "java version 1.8.0_31", your version is recent enough (as for some reason version 8 is synonymous with version 1.8). If this is not the case, the latest version for your operating system can be downloaded from [http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html), after accepting the license agreement.

* **BEAST 2:** The Bayesian phylogenetic software BEAST 2 comes as a package combining BEAST itself with BEAUti, TreeAnnotator, and other tools. Of these, BEAUti will be used to generate input for BEAST in XML format, BEAST will be used to run the Bayesian phylogenetic analysis, and the results will be summarized with TreeAnnotator. Each of these tools runs on Linux, Mac, and Windows machines. The latest version of the BEAST 2 package can be downloaded from [http://beast2.org](http://beast2.org).

* **PhyloNet:** The latest binary jar file of the software PhyloNet can be downloaded from [http://bioinfo.cs.rice.edu/phylonet](http://bioinfo.cs.rice.edu/phylonet), which should run on all operating systems. After download, the file should be placed in a directory that can easily be accessed from the command line.

#### Identification of genomic regions for phylogenetic inference

In this part of the tutorial, the sofware Saguaro will be used to detect boundaries between genomic regions that are characterized by different phylogenetic histories. However, for computational reasons, Saguaro does not infer these phylogenetic histories directly. Instead, the analysis performed by Saguaro is based on what the authors call "cacti", sets of distance matrices that describe how different each pair of genomes is relative to all others. For the purpose of this analysis, these cacti can be considered as proxies for phylogenetic histories, as the difference between pairs of genomes is obviously linked to their phylogenetic relatedness. However, to reconstruct local phylogenetic histories more accurately, Bayesian phylogenetic analysis will be performed subsequently for the genomic regions identified with Saguaro.

At the beginning of the analysis, Saguaro will calculate a single cactus for the entire alignment, and a score is calcuated for each variable alignment position, describing the fit between this site and the first cactus. Based on these scores, genomic regions with a poor fit to the current cactus are identied with the hidden Markov model implemented in Saguaro, and a new cactus is defined for these. This process is repeated multiple times, thus further partitioning the alignment into segments, and at the same time assigning one out of an increasing set of cacti to each segment. Details of this procedure are described in Zamani et al. ([2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)).

As a full Saguaro analysis can take several hours to days depending on the size of the sequence alignment, only a preliminary analysis of a single chromosome (linkage group 5 of the tilapia genome assembly) will be done as part of this tutorial.

* As a first step, please **download the sequence alignment** for linkage group 5 from [https://github.com/mmatschiner/Introgression-Tutorial/blob/master/LG05.fasta.zip?raw=true](https://github.com/mmatschiner/Introgression-Tutorial/blob/master/LG05.fasta.zip?raw=true) (the download may take a few moments to start). Save it in a directory, in which you'ld like to run this tutorial.

*  **Unzip the alignment file** by double-clicking, or by typing

		unzip LG05.fasta.zip
in a console window (after navigating to the directory in which you saved it).

Saguaro expects its input in a particular binary format, which the authors call a "feature". To allow the conversion of commonly used file formats into this feature format, the Saguaro installation includes several conversion tools like Fasta2HMMFeature, Maf2HMMFeature, and VCF2HMMFeature. These should be located in the same directory as the Saguaro executable.

* Type the following in a console window to **convert the alignment into a Saguaro feature**:

		Fasta2HMMFeature -i LG05.fasta -o LG05.feature -n LG05 -m 4 > LG05.out
Here, the input file is specified with the `-i` option, the name of the output file in "feature" format is given with the `-o` option, and a name for the linkage group is specified with `-n`. By using `-m 4`, only sites with a minimum coverage of 4 (i.e. missing data for maximally 3 out of the 7 sequences in the alignment) are included in the "feature" file. To avoid too much output on the screen, the output is sent to a file that we name `LG05.out`.

* **Run Saguaro** on the "feature" file for linkage group 5 by typing

		Saguaro -f LG05.feature -o saguaro_results -iter 16 -cycle 1 -neurons 100
This specifies that file `LG05.feature` is used as input for Saguaro and that all results are written to a directory called `saguaro_results`, which is created by Saguaro. If more than one alignment was used as the input, these could be specified something like `-l all_features.txt`, where the file `all_features.txt` would have to contain only a list of the file names of all input "feature" files. We here use the options `-iter 16`, `-cycle 1`, and `-neurons 100` only in order to reduce the run time for the purpose of this tutorial. With these options, Saguaro runs 16 iterations in which cacti are added, and each iteration contains a single cycle in which cacti are optimized and assigned to sites. The 100 "neurons" are used to generate a new cactus in each iteration based on the sites that have the lowest fit to any of the cacti already included in the model. If we would not specify these options, Saguaro would by default perform 40 iterations with 2 cycles per iteration, and it would use 800 neurons. Each of these changes would increase the run duration, but probably also the accuracy of the results, and it is recommend to keep the default options whenever the resulting run time is acceptable. With the settings used here, the Saguaro analysis should take 10-15 minutes, depending on the speed of your computer.

* Once the Saguaro analysis has finished, results should have been written to directory `saguaro_results`. **Have a look at the content of this directory**. You'll see the following files:

		HMMTrain.out.0
		HMMTrain.out.1
		HMMTrain.out.2
		...
		LocalTrees.out
		saguaro.cactus
		saguaro.config
		saguaro.garbage
		saguaro.garbage.vec
Of these, only `LocalTrees.out` and `saguaro.cactus` are interesting for us.

* **Open file** `saguaro.cactus` in a text editor, or in a console window (e.g. by typing `less saguaro_results/saguaro.cactus`). You should see more or less the following content:

		cactus0
		Br	Gr	Ma	Mz	Ol	Pu	On
		Br	0.021704	0.236808	0.262940	1.824387	0.240841	0.220309	1.722597
		Gr	0.236808	0.020616	0.208407	1.680392	0.236869	0.245299	1.825267
		Ma	0.262940	0.208407	0.009898	1.679830	0.227633	0.246533	1.799380
		Mz	1.824387	1.680392	1.679830	0.001310	1.879314	1.994391	0.116698
		Ol	0.240841	0.236869	0.227633	1.879314	0.020897	0.145890	2.067251
		Pu	0.220309	0.245299	0.246533	1.994391	0.145890	0.003647	2.138905
		On	1.722597	1.825267	1.799380	0.116698	2.067251	2.138905	-0.000000
		cactus1
		...
This is the distance matrix for the first cactus (cactus0), followed by distance matrices for all other cacti.

The two-character identifiers "Br", "Gr", etc. are the codes used here for the seven species included in the analysis. They will be used throughout this tutorial.

| Code | Species                        | Group    |
|------|--------------------------------|----------|
| Br   | *Neolamprologus brichardi*     | Ingroup  |
| Gr   | *Neolamprologus gracilis*      | Ingroup  |
| Ma   | *Neolamprologus marunguensis*  | Ingroup  |
| Mz   | *Metriaclima zebra*            | Outgroup |
| Ol   | *Neolamprologus olivaceous*    | Ingroup  |
| Pu   | *Neolamprologus pulcher*       | Ingroup  |
| On   | *Oreochromis niloticus*        | Outgroup |

* Next, **open file** `LocalTrees.out` in a text editor or a console window. You'll see something like this:

		Reading features...
		done!
		Reading models.
		Dynprog'ing...
		Setting up HMM...
		Adding word cactus0 as # 0
		...
Scroll down a bit to this part:

		...
		Processed: 380000 (98.559 %)
		REPORTING Traceback and Update
		cactus14	LG05: 1278 - 487927	length: 486649	(frames 1-3313 l=3312) 	score=20.7077
		Br	Gr	Ma	Mz	Ol	Pu	On
		Br	0.02	0.59	0.47	1.36	0.46	0.45	1.14
		Gr	0.59	0.03	0.50	1.04	0.55	0.54	1.23
		Ma	0.47	0.50	0.03	1.34	0.35	0.31	1.52
		Mz	1.36	1.04	1.34	0.01	1.43	1.42	0.29
		Ol	0.46	0.55	0.35	1.43	0.03	0.27	1.63
		Pu	0.45	0.54	0.31	1.42	0.27	0.03	1.64
		On	1.14	1.23	1.52	0.29	1.63	1.64	-0.00
		cactus0	LG05: 487940 - 499860	length: 11920	(frames 3314-3504 l=190) 	score=378.688
		...
These lines contain information on the first segment identified by Saguaro. It is assigned to cactus14, and is located between position 1278 and position 487927 on linkage group 5 (LG05 - not too surprisingly, as this is the only linkage group used here). In contrast to distance matrices given in file `saguaro.cactus`, the distance matrices given in this in `LocalTrees.out` are calculated per segment, not per cactus. Nevertheless, the distance matrix of the first segment between position 1278 and position 487927 is likely very similar to that of cactus14, otherwise, this segment would not have been assigned to this cactus.

As we're not particularly interested in the distance matrices of segments, but more in the placement of segment boundaries so that we can select alignment regions for phylogenetic analyses that are not broken up by boundaries, the most imporant information for us is in the header lines for each segment.

* To **see only header information** for each cactus, type

		cat saguaro_results/LocalTrees.out | grep length
You should see output like this:

		cactus14	LG05: 1278 - 487927	length: 486649	(frames 1-3313 l=3312) 	score=20.7077
		cactus0	LG05: 487940 - 499860	length: 11920	(frames 3314-3504 l=190) 	score=378.688
		cactus9	LG05: 499901 - 706017	length: 206116	(frames 3505-4791 l=1286) 	score=77.6632
		cactus0	LG05: 706087 - 733860	length: 27773	(frames 4792-5075 l=283) 	score=365.399
		cactus14	LG05: 736928 - 1158725	length: 421797	(frames 5076-8167 l=3091) 	score=54.6317
		cactus0	LG05: 1158770 - 1187064	length: 28294	(frames 8168-8338 l=170) 	score=1003.93
		...
		
* To **find out how many segments** Saguaro has identified, you could type

 		cat saguaro_results/LocalTrees.out | grep length | wc -l
 		
* In order to visualize segment boundaries and the cacti assigned to segments, **download the Ruby script paint_chromosomes.rb** from [https://github.com/mmatschiner/Introgression-Tutorial/blob/master/paint_chromosomes.rb?raw=true](https://github.com/mmatschiner/Introgression-Tutorial/blob/master/paint_chromosomes.rb?raw=true), and save it in the directory that you're using for this tutorial.

* **See whether the script can run** on your machine by typing

		ruby paint_chromosomes.rb
If you see this short help text, everything seems fine:

		paint_chromosomes.rb

		This script uses output from the software Saguaro to paint
		chromosomes according to Saguaro cacti. The output will be
		written in svg format.

		This script should be run e.g. with
		ruby paint_chromosomes.rb LocalTrees.out LocalTrees.svg
		where 'LocalTrees.out' should be replaced with the actual path
		to the Saguaro output file.

* Now run the script with your Saguaro output file

		ruby paint_chromosomes.rb saguaro_results/LocalTrees.out
This will generate a vector graphic file in SVG format that will be written to the same directory in which `LocalTrees.out` is placed (i.e. `saguaro_results`), and it will be named `LocalTrees.svg`.

* Open the vector graphic file `LocalTrees.svg` in Firefox or another web browser. This is what you should see:<br>
![LocalTrees.svg](https://rawgit.com/mmatschiner/Introgression-Tutorial/master/LocalTrees.svg "LocalTrees.svg")<br>
In this image, segments assigned to the most common cactus are drawn in dark gray, and segments assigned to other cacti are shown in red, orange, cyan, and light green, purple (in decreasing frequency). With more than six different cacti, all remaining cacti are shown in light gray. As you can see, only four cacti are common (dark gray, red, orange, and cyan). Also, you'll notice that the most frequent cactus (in dark gray) is mostly found in the center of the linage group, while other cacti dominate towards the ends of the linkage group. These results may not be particularly accurate, as we've only performed a preliminary Saguaro analysis. For comparison, this is how the results from a much longer analysis with default settings would look like:<br>
![LocalTrees_full.svg](https://rawgit.com/mmatschiner/Introgression-Tutorial/master/LocalTrees_full.svg "LocalTrees_full.svg")<br>
You'll notice that some details are different in these two images, however, the overall pattern is the same.