#' Create an integrative multiassay container
#' 
#' An integrative package for coordinated representation of multiple experiments
#' on partially overlapping samples, with associated metadata at the level of
#' an entire study and at the level of the "biological unit". The biological
#' unit may be a patient, plant, yeast strain, etc. This package is designed
#' around the following hierarchy of information: 
#' study (highest level): The study can encompass several different types of
#' experiments performed on one set of biological units, for example cancer
#' patients. A MultiAssayExperiment represents a whole study, containing:
#' 1. metadata about the study as a whole
#' 2. metadata about each biological unit: for example, age, grade, stage
#' for cancer patients
#' 3. results from a set of experiments performed on the biological units
#' 4. a map for matching data from the experiments back to the corresponding
#' biological units.
#' experiment: A set of assays of a single type performed on some or all of 
#' the biological units. It is permissible that an experiment may be performed 
#' only on a subset of the biological units, and may be performed in duplicate 
#' on some of the biological units. For example, an experiment could be somatic 
#' mutation calls for some or all of the biological units. 
#' Experiments may be ID-based, where measurements are indexed identifiers of 
#' genes, microRNA, proteins, microbes, etc. Alternatively, experiments may be 
#' range-based, where measurements correspond to genomic ranges that can be 
#' represented as GRanges objects, such as gene expression or copy number. Note
#' that for ID-based experiments, there is no requirement that the same IDs be 
#' present for different experiments. For range-based experiments, there is 
#' also no requirement that the same ranges be present for different 
#' experiments; furthermore, it is possible for different samples within an 
#' experiment to be represented by different ranges.
#' samples (lowest level): An individual set of measurements performed on a 
#' single biological unit. These measurements are indexed either by identifiers 
#' (IDs) or ranges, and these IDs or ranges are referred to here as features.
"_PACKAGE"