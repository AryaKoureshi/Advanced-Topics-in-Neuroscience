![Logo](logo.png)
###### [Arya Koureshi](https://aryakoureshi.github.io)
###### 401204008
---

#### Introduction
Irregular spontaneous activity, i.e., activity that is not related in any obvious way to external stimulation, and trial-to-trial variations in neuronal responses are often considered as noise. The origin of this irregularity of neuronal dynamics in vivo is poorly understood. For instance, trial-to-trial fluctuations in response strength are shared between neurons, and spikes often occur synchronously. Understanding the properties and mechanisms that generate these forms of correlation is critical for determining their role in cortical processing. Matthew Smith and Adam Kohn have investigated the spatial extent and functional specificity of correlated spontaneous and evoked activity. In this exercise, we are going to replicate some part of their research and manipulate neuronal data.

---

#### Data description
Information about the data including cells, stimuli and format of the data are in crcns_pvc11_data_description file. For further details on how the data were obtained, see the reference 2.

---

#### Tuning curve and noise correlation
Tuning curve represents the average response of a neuron to a set of stimuli, with the average taken across many presentations of each stimulus, and the noise refers to the trial-to-trial variability in the responses. In panel a, tuning curves are shown for two neurons that have slightly different preferred stimuli. In panel b, we show two hypothetical scatter plots of the single trial responses for this pair of neurons, in response to the repeated presentation of a single stimulus s1 (arrow in panel a). Ellipses represent 95% confidence intervals. The example on the left illustrates positive noise correlation and the example on the right illustrates negative noise correlation. Responses also show a second sort of correlation known as signal correlation. These are correlations in the average response. Neurons with similar tuning curves (panel a) typically have positive signal correlations, because when s increases, the mean responses of both neurons tend to increase, or decrease, together. Conversely, neurons with dissimilar tuning curves typically have negative signal correlations.
![Tuning Curves](TuningCurves.png)
![Examples](Examples.png)
![Uncorrelated - Correlated](Uncorrelated_Correlated.png)

---

#### Questions
First you have to read the reference 2 ( just read the abstracts and go through the figures and read the captions to get an idea) and download their data

1.  As mentioned in data_description file, evoked data set (‘spikes_gratings’) comprises spiking activity from 44 to 104 neurons collected from 3 arrays (data_monkey1_gratings, up to 3). So you need to find data_monkey1_gratings, 1, 2 and 3 files; and then find the most active neuron in each array. Then plot their Tuning curves.

##### PSTHs (V1 neurons)
<img width="410" height="315" src="Q1_PSTH_0°.png">
<img width="410" height="315" src="Q1_PSTH_30°.png">
<img width="410" height="315" src="Q1_PSTH_60°.png">
<img width="410" height="315" src="Q1_PSTH_90°.png">
<img width="410" height="315" src="Q1_PSTH_120°.png">
<img width="410" height="315" src="Q1_PSTH_150°.png">
<img width="410" height="315" src="Q1_PSTH_180°.png">
<img width="410" height="315" src="Q1_PSTH_210°.png">
<img width="410" height="315" src="Q1_PSTH_240°.png">
<img width="410" height="315" src="Q1_PSTH_270°.png">
<img width="410" height="315" src="Q1_PSTH_300°.png">
<img width="410" height="315" src="Q1_PSTH_330°.png">

PSTHs are plotted for the most active neuron in each monkey.

##### Tuning Curves
<img width="410" height="315" src="Q1_TC_Monkey1.png">
<img width="410" height="315" src="Q1_TC_Monkey2.png">
<img width="410" height="315" src="Q1_TC_Monkey3.png">

The most active neuron exhibits the highest firing rate in response to each stimulus. Therefore, I identified these neurons by determining the maximum responses to various orientations.

---

2. Use MAP file and CHANNELS files; to Plot 3(10×10) color mesh for these 3 arrays and use same color for neurons with similar preferred orientation. (as shown below). Is your results similar to pinwheel organization of orientation in the cortex? Why or why not?

![Q2Exm](Q2Exm.png)

The preferred orientation of each neuron is determined by identifying the neuron's maximum response to various orientations.

![Spatial-distribution-of-orientation-preferences-in-monkey-V1-a-Orientation-preference](Spatial-distribution-of-orientation-preferences-in-monkey-V1-a-Orientation-preference.png)

Spatial distribution of orientation preferences in monkey V1. (a) Orientation preference map obtained by optical imaging. Color scale at the bottom left shows preferred orientations in degrees. White and black dots indicate clockwise and counterclockwise pinwheel centers (PWCs), respectively. The white circle represents an example of a region of interest (ROI), and the plus sign indicates the center of the ROI. (b) Magnified view of ROI, consisting of 20 concentric annular regions whose radii vary at 0.5 mm intervals (black circles). Color scale at bottom right shows the range of preferred orientation relative to the preferred orientation at the ROI center (plus sign). The top right panel shows the orientation population within the ROI as a function of annular region radius. (c) Variation of orientation populations in ROIs across the PWC. The imaged area (4.5 × 3.9 mm2) was scanned by shifting a circular ROI, and the data were pooled into five groups according to the distance between the ROI center and the nearest PWC: 0–0.04 to 0.16–0.20 mm from the PWC. Top row, orientation populations obtained from the median of pooled data in each group plotted as a function of annular region radius. Bottom row, the orientation distribution indices (ODIs) plotted as a function of annular region radius. Solid lines and gray areas indicate the median ODI and the range from 25th to 75th percentiles of ODIs, respectively. The ODIs in many regions differ significantly from zero (red asterisks, P < 0.0025). Data for all monkeys are in Supplementary Fig. S1.


##### Preferred Orientations (0° - 330°)
<img width="410" height="315" src="Q2_PO1(0°-330°)_Monkey1.png">
<img width="410" height="315" src="Q2_PO1(0°-330°)_Monkey2.png">
<img width="410" height="315" src="Q2_PO1(0°-330°)_Monkey3.png">
<img width="410" height="315" src="Q2_PO2(0°-330°)_Monkey1.png">
<img width="410" height="315" src="Q2_PO2(0°-330°)_Monkey2.png">
<img width="410" height="315" src="Q2_PO2(0°-330°)_Monkey3.png">

In the figures, it is evident that only Monkey 3 exhibits this pattern to a certain degree. However, it is important to note that spatial distribution of orientation preferences in monkey V1 represents orientations between 0-180 degrees, suggesting that the paper's authors either wrapped the stimulus orientations within this range or solely used orientations between 0-180 degrees. I have followed the same approach in below. While the pinwheel pattern appears slightly more apparent in this figure compared to (0° - 330°) figures, it still does not display a distinct pinwheel pattern.

##### Preferred Orientations (0° - 150°)
<img width="410" height="315" src="Q2_PO3(0°-150°)_Monkey1.png">
<img width="410" height="315" src="Q2_PO3(0°-150°)_Monkey2.png">
<img width="410" height="315" src="Q2_PO3(0°-150°)_Monkey3.png">
<img width="410" height="315" src="Q2_PO4(0°-150°)_Monkey1.png">
<img width="410" height="315" src="Q2_PO4(0°-150°)_Monkey2.png">
<img width="410" height="315" src="Q2_PO4(0°-150°)_Monkey3.png">

---

3. You need to find the dependence of $r_{sc}$ (noise correlation) on distance for pairs grouped based on their orientation tuning similarity. In other word, you should find most populated group of neurons with similar orientation preferences, and investigate the relation of correlation and distance of pair neurons. Is your answer similar with paper's conclusion?

##### No Shuffling vs Shuffled
The results depicted in the figures are similar to the results in the paper. However, to validate these findings, I shuffled the trials and generated plots for each monkey.

* Monkey 1

<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey1_Shuffle=0.png">
<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey1_Shuffle=1.png">

<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey1_Shuffle=0.png">
<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey1_Shuffle=1.png">

<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey1_Shuffle=0.png">
<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey1_Shuffle=1.png">

* Monkey 2

<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey2_Shuffle=0.png">
<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey2_Shuffle=1.png">

<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey2_Shuffle=0.png">
<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey2_Shuffle=1.png">

<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey2_Shuffle=0.png">
<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey2_Shuffle=1.png">

* Monkey 3

<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey3_Shuffle=0.png">
<img width="600" height="475" src="Q3_1_(rsc_distance)_Monkey3_Shuffle=1.png">

<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey3_Shuffle=0.png">
<img width="600" height="475" src="Q3_2_(rsc_rs)_Monkey3_Shuffle=1.png">

<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey3_Shuffle=0.png">
<img width="600" height="475" src="Q3_3_(rs_rsc_distance)_Monkey3_Shuffle=1.png">

The difference between shuffled and non-shuffled trials is more evident.

The figures show that the downward trend in the $r_{sc}$ vs $Distance$ plot and the upward trend in the $r_s$ vs $r_{sc}$ plot have been removed.

---

4. As mentioned earlier in the paper, the spatial structure of correlation was similar between spontaneous and evoked activity, there was a striking difference in its strength: the average $r_{sc}$ value for spontaneous activity was 0.299, nearly twofold higher than the average correlation of evoked activity (0.176) in dataset. How can you explain his difference in correlation between evoked and spontaneous periods?

In the paper, it is stated that there is a difference in correlation between spontaneous and evoked activity, indicating that sensory input can significantly decrease the correlation of ongoing activity. This suggests that input can reduce noise correlation, which is logical, as neurons should focus on the new input rather than spontaneous activity. In summary, the paper states that the spatial structure of correlated activity is similar for both spontaneous and evoked activities, suggesting that they might originate from the same mechanisms and circuits. However, when spontaneous activity occurs for extended periods without interruption, the correlation is roughly twice as strong as during evoked activity. The decrease in correlation due to stimulus is quick, prominently observed at the beginning and end of the stimulus. Following the stimulus conclusion, the correlation gradually increases over the 1.5 s ISI. These findings suggest that the mechanisms causing correlated variability are influenced by network inputs, with effects lasting several seconds after the stimulus has ceased.

---

5. [Additional Score question not mandatory] So we know that the neurons with similar orientation preferences are highly correlated with each other. Can we use this fact and compute color mesh (in question2) from spontaneous activity of neurons? Explain your approach

No, because there is only one condition (not 12), which makes it impossible to compute $r_s$.

---

#### References
1. Averbeck, Bruno B., Peter E. Latham, and Alexandre Pouget. "Neural correlations, population coding and computation." Nature reviews neuroscience 7.5 (2006): 358.
2. Smith, Matthew A., and Adam Kohn. "Spatial and temporal scales of neuronal correlation in primary visual cortex." Journal of Neuroscience 28.48 (2008): 12591-12603.
3. Gerstner, Wulfram, et al. Neuronal dynamics: From single neurons to networks and models of cognition. Cambridge University Press, 2014.
4. Mainen, Zachary F., and Terrence J. Sejnowski. "Reliability of spike timing in neocortical neurons." Science 268.5216 (1995): 1503-1506.
5. Kohn, A., Smith, M.A. (2016) Utah array extracellular recordings of spontaneous and visually evoked activity from anesthetized macaque primary visual cortex (V1). CRCNS.org   
http://dx.doi.org/10.6080/K0NC5Z4X

