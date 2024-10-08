# This is archived. Please use the upstream instead.

# Microwave Radar-based Imaging Toolbox: Efficient Reconstruction Software

MERIT provides free software algorithms for Delay-and-Sum reconstruction of
microwave imaging signals in medical and industrial settings.
MERIT allows users to easily and configurably test different algorithms,
easily switching between time and frequency domain signal representations.
All inbuilt algorithms can be configured to run in parallel or on GPU without
changing the code.
Features include:
 - visualize signals: view and compare signals;
 - manage signals: reorder antenna numbering, exclude channels etc., limit to
    monostatic/bistatic signals etc.;
 - estimate propagation paths based on transmit and receive locations;
 - delay signals based on propagation paths through multiple media
    with dispersive dielectric properties;
 - image using a highly configurable and extensible set of beamforming
     algorithms such as delay-and-sum;
 - analyze resulting images using a variety of metrics such as
    signal-to-clutter and signal-to-mean ratios;
 - and visualize image results in two and three dimensions.

![example_3D_scan](/.showcase/example_3D_scan.png)

# Examples

MERIT is designed to make the imaging code short, clear and efficient. For
example:

```matlab
%% Load sample data (antenna locations, frequencies and signals)
frequencies = dlmread('example_data/frequencies.csv');
antenna_locations = dlmread('example_data/antenna_locations.csv');
channel_names = dlmread('example_data/channel_names.csv');

scan1 = dlmread('example_data/B0_P3_p000.csv');
scan2 = dlmread('example_data/B0_P3_p036.csv');

%% Perform rotation subtraction
signals = scan1-scan2;

%% Generate imaging domain
[points, axes_] = merit.domain.hemisphere('radius', 7e-2, 'resolution', 2.5e-3);

%% Calculate delays for synthetic focusing
delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  'relative_permittivity', 8);

%% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

%% Plot image using MATLAB functions
im_slice = merit.visualize.get_slice(img, points, axes_, 'z', 35e-3);
imagesc(axes_{1:2}, im_slice);
```

In a few lines of code, radar-based images can be efficiently created.
MERIT allows the user to change the beamformer, imaging domain and other features easily and simply.
Functions are designed to accept options allowing the user to easily change the imaging procedure.

# Getting started

To try MERIT:

Follow the [Getting Started Guide](https://github.com/EMFMed/MERIT/wiki/Getting-Started) which shows, using a step-by-step guide, how to load, process, image and visualise microwave breast images. Bug reports, feature requests, or code or documentation contributions are welcome.

If you have found MERIT useful and publish your work, we would be grateful if you could cite us using:

D. O’Loughlin, M. A. Elahi, E. Porter, et al., "Open-source Software for Microwave Radar-based Image Reconstruction", Proceedings of the 12th European Conference on Antennas and Propagation (EuCAP), London, the UK, 9-13 April.

# Publications using MERIT

If you have used and cited MERIT, please consider adding your publication here (using the Github Issue Tracker: [how-to guide](https://help.github.com/articles/creating-an-issue/)).

Research publications which have used MERIT:

## Journal Publications

* [Sensitivity and Specificity Estimation Using Patient-Specific Microwave Imaging in Diverse Experimental Breast Phantoms](https://ieeexplore.ieee.org/abstract/document/8428660)

* [Microwave Breast Imaging: experimental tumour phantoms for the evaluation of new breast cancer diagnosis systems](http://iopscience.iop.org/article/10.1088/2057-1976/aaaaff/meta)

* [Parameter Search Algorithms for Microwave Radar-Based Breast Imaging: Focal Quality Metrics as Fitness Functions](http://www.mdpi.com/1424-8220/17/12/2823)

* [A 16‐modified antipodal Vivaldi antenna array for microwave‐based breast tumor imaging applications](https://onlinelibrary.wiley.com/doi/full/10.1002/mop.31873)

## Conference Publications

* [Effects of Interpatient Variance on Microwave Breast Images: Experimental Evaluation](https://ieeexplore.ieee.org/abstract/document/8513673)

* [Evaluation of Experimental Microwave Radar-Based Images: Evaluation Criteria](https://ieeexplore.ieee.org/abstract/document/8608682)

# Changelog

The most recent version is 0.1.0.
Notable changes to MERIT are recorded in CHANGELOG.
The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

# License

MERIT is available under the Apache 2.0 license. See LICENSE for more information.
