#!/bin/bash

# Project simulated parental origins onto markers to produce F1/F2 truth sets.
python python/04_build_simulation_truth.py --truth-dir simulation_truth --panel simulation_panel/markers_sim.tsv --f1-list pseudo_f1_pairs.tsv

# Lift the MM-coordinate marker, haplotype-block and breakpoint truth sets to SL6.
python python/05_lift_truth_coordinates.py --truth-dir simulation_truth --panel-dir simulation_panel --liftover-dir liftover
