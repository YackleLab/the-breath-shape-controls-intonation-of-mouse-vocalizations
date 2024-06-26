# the-breath-shape-controls-intonation-of-mouse-vocalizations
MATLAB code for analyses described in MacDonald et al 2024 (https://elifesciences.org/reviewed-preprints/93079)
This repository contains the code for analysis of mouse vocalization, breathing +/- optogenetic data collected and reported in our paper titled 'the breath shape controls intonation of mouse vocalizations'
The data are available on the Dryad repository (https://datadryad.org/stash/dataset/doi:10.5061/dryad.n8pk0p34d)

**Overview**
The basic workflow of these scripts is to have one script that batch processes the recordings and saves relevant information in a .mat file and then have a second script that opens these .mat file and groups the variables in a relevant way (e.g. segmenting all the breaths and USVs in a recording, saving them in a .mat and then a second script that opens them and groups by USV class)

The first two figures of the paper contain analysis of recordings of mice vocalizing and breathing in a plethysmography chamber.
To analyze the features of breaths with and without USVs (e.g. Figure 1): first run 'batchProcessStats.m' which runs through each recording (specified in a metadata table) and outputs summary statistics as a .mat file
To group summary statistics and make plots run 'compileSummaryData.m'

To analyze the features of breaths containing different classes of USV (e.g. Figure 2, Figure 5): first create pitch vectors from the spectrogram with 'peakFreqWriter.m' and analyze the recording with vocalMat (https://github.com/ahof1704/VocalMat) for classification of USV types
Then use 'callVarsPreProcessing.m' to gather each vocal breath, pitch vector, onset/offset and vocal class and save in a .mat.
Finally run the various compliler scripts to generate summary data and/or raster plots

To analyze the abaility of optogenetic stimulation to produce USVs run 'writeOptoVarsBatch'.

Note, the basic breathing analysis uses functions from Bachmutsky et al 2020. Availble in this respository: https://github.com/YackleLab/Opioids-depress-breathing-through-two-small-brainstem-sites
