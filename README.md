# Advanced-Topics-in-Neuroscience
This repository contains homework assignments and projects completed for the course "Advanced Topics in Neuroscience" instructed by Dr. Ali Ghazizadeh at Sharif University of Technology.

---

## [Simulation #1 (Neural Encoding)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/1)
This assignment is divided into two main sections, each dealing with a different type of neuron model.
### Integrate and Fire Neuron
Here we are going to simulate spike trains based on [Softky & Koch, 1993] and compare its statistics with real neural data. The simplest neuron model in this paper was a kind of integrate and fire neuron where it’s inter-spike intervals (ISI) distribution explained by equation [9].
### Leaky Integrate and Fire Neuron
In this section we are going to implement a more realistic model of neuron, Leaky Integrate and fire (LIF) neuron, which consider the leakage of postsynaptic inputs. 
The leaky integrate-and-fire (LIF) neuron is probably one of the simplest spiking neuron models, but it is still very popular due to the ease with which it can be analyzed and simulated.
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/1/Papers)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/1/AdvNeuro_HW1_AryaKoureshi.ipynb)

---

## [Simulation #2 (Neural Encoding)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/2)
### Aim
Studying the population response structure
### Task Description
Uncertainty is a ubiquitous component of our environment, such that humans and animals are regularly confronted with conditions that vary in degrees of uncertainty. It is therefore a fundamental requirement of our brain systems to accurately process uncertain information so that we may function appropriately, both in our mundane day-to-day activities and in more profound moments. Many brain regions including the frontal cortex, basal ganglia, amygdala, parietal cortex, cingulate cortex, and insular cortex have been identified as key areas in processing information about uncertain rewards. In this assignment we are going to analyze the activity of a population of single units recorded with multi-electrode array in Parietal cortex. The task is designed to study the encoding of reward expected value in area 7a. This area encodes the spatial location of cue.
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/2/Papers)
### [Data](https://drive.google.com/file/d/1FdO7RF8IJpGe76sMRqa2oDxob_IUnzYs/view?usp=share_link)
### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/2/Codes)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Neural-Encoding/2/AdvNeuro_HW2_AryaKoureshi.html)

---

## [Simulation #3 (Noise Correlation)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Noise-Correlation)
### Introduction
Irregular spontaneous activity, i.e., activity that is not related in any obvious way to external stimulation, and trial-to-trial variations in neuronal responses are often considered as noise. The origin of this irregularity of neuronal dynamics in vivo is poorly understood. For instance, trial-to-trial fluctuations in response strength are shared between neurons, and spikes often occur synchronously. Understanding the properties and mechanisms that generate these forms of correlation is critical for determining their role in cortical processing. Matthew Smith and Adam Kohn have investigated the spatial extent and functional specificity of correlated spontaneous and evoked activity. In this exercise, we are going to replicate some part of their research and manipulate neuronal data.
### Data Description
Information about the data including cells, stimuli and format of the data are in crcns_pvc-11_data_description file. For further details on how the data were obtained, see the reference 2.
### [Data](https://drive.google.com/file/d/10ZDjIfTTnw_BUO3UzsCGDLYa-_2wyBYn/view?usp=share_link)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Noise-Correlation/Papers)
### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Noise-Correlation/Codes)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Noise-Correlation/AdvNeuro_HW3_AryaKoureshi.html)

---

## [Simulation #4 (Traveling Waves)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Traveling-Waves)
### Introduction
In this assignment we are going to analyze the activity of local field potential (LFP) recorded with multielectrode array in motor cortex. The task is designed to study the encoding of working memory in Premotor Area F5. Although searching for significant correlations between neural signals (including spikes and LFP) and task variables might be interesting itself, here in this assignment we are going to investigate existence of mesoscopic traveling waves and their properties. Accordingly, you only need to know task timings of rather than task details.
### Task Description
In this task, a monkey learns to follow instruction to receive tasty juice as reward. At first, the monkey should keep its eyes on a fixation point in center of screen (red/green circle). To make sure monkey is looking at the desired point, we follow its eyes using an eye tracker. After a random duration between 300 to 500 ms fixating at fixation point, a cue will be presented in the periphery. Monkey must remember the spatial position of cue to either saccade or reach to that point after waiting for center light offset. At the end, based on the color of target monkey reach/saccade to the expected location. If the monkey saccade/reach in the correct position, he will receive tasty juice!
### Data Description
After the monkey gets fully trained in the task, we start recording from the Premotor Area F5, using an electrode array. This array includes 49 single electrodes distributed in an area of $12𝑚𝑚^2$. Distance between neighboring electrodes are 400 𝜇𝑚. In assignment folder you can find a session of electrophysiological recordings including timing of task.
### [Data](https://drive.google.com/file/d/1PalKnMMMLB7FbtcpTIyX80JtU_bxoHky/view?usp=sharing)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Traveling-Waves/Papers)
### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Traveling-Waves/Code)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Traveling-Waves/AdvNeuro_HW4_AryaKoureshi.html)
### [Video](https://drive.google.com/file/d/1qfLWo0i4vuWrkFz_yl4cc-hHGifxVBX9/view?usp=sharing)
### [Results](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Traveling-Waves/Results)

---

## [Simulation #5 (Motivation and Classical Conditioning)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Motivation-and-Classical-Conditioning)
### Aim
Modelling learning in classical conditioning paradigms
### Description
We often need to rapidly learn about the value of new stimuli that we encounter or be ready for changes to familiar stimulus values that we were familiar without prior notice. Classical conditioning includes a large group of paradigms where different combination of reward histories with stimuli can form our judgement about their current values.

A well know model in classical conditioning and reinforcement learning is the Rescola-Wagner. This model predicts that violations of our expected reward for each stimuli or combination of stimuli causes incremental changes in our belief about their values.

### [Notebook (Code)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Motivation-and-Classical-Conditioning/AdvNeuro_HW5_AryaKoureshi.ipynb)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Motivation-and-Classical-Conditioning/Papers)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Motivation-and-Classical-Conditioning/AdvNeuro_HW5_AryaKoureshi.html)

---

## [Simulation #6 (Reinforcement Learning)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Reinforcement-Learning)
### Learning the Water Maze
As an example of generalized reinforcement learning, we consider the water maze task. This is a navigation problem in which rats are placed in a large pool of milky water and have to swim around until they find a small platform that is submerged slightly below the surface of the water. The opaqueness of the water prevents them from seeing the platform directly, and their natural aversion to water (although they are competent swimmers) motivates them to find the platform. After several trials, the rats learn the location of the platform and swim directly to it when placed in the water.
We are going to simulate a simple model of navigation problem.

### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Reinforcement-Learning/Code)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Reinforcement-Learning/Papers)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Reinforcement-Learning/AdvNeuro_HW6_AryaKoureshi.html)
### [Results](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Reinforcement-Learning/Results)

---

## [Simulation #7 (Evidence Accumulation)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Evidence-Accumulation)
### Aim
Simulation of evidence accumulation model and studying the relationship between accuracy and reaction time (RT) in decision making.

### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Evidence-Accumulation/Code)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Evidence-Accumulation/Papers)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Evidence-Accumulation/AdvNeuro_HW7_AryaKoureshi.html)
### [Results](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Evidence-Accumulation/Results)

---

## [Simulation #8 (Visual Attention)](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Visual-Attention)
### Learning to predict where humans look
For many applications in graphics, design, and human computer interaction, it is essential to understand where humans look in a scene. Where eye tracking devices are not a viable option, models of saliency can be used to predict fixation locations. Most saliency approaches are based on bottom-up computation that does not consider top-down image semantics and often does not match actual eye movements. To address this problem, we used eye tracking data of 15 viewers on 1003 images and use this database as training and testing examples to learn a model of saliency based on low, middle and high-level image features. This large database of eye tracking data is publicly available with the paper (Tilke Judd).
### [Database](http://people.csail.mit.edu/tjudd/WherePeopleLook/index.html)

### [Code](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Visual-Attention/Code)
### [Papers](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Visual-Attention/Papers)
### [Report](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Visual-Attention/AdvNeuro_HW8_AryaKoureshi.html)
### [Results](https://github.com/AryaKoureshi/Advanced-Topics-in-Neuroscience/tree/main/Visual-Attention/Results)

---
