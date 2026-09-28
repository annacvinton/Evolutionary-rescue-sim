R pipeline for the v3 dataset (all base R, no packages)

Inputs (from the Drive folder top level and Supplement/):
  run_summary_v3.csv    one row per run
  controls_v3.csv       unperturbed control arm
  traj_subset_v3.csv    trajectory slice for Fig 1 (slope 0.8, disp 3, mut on, sev 11)
  kernel_check2.csv     extended kernel check (Supplement/; optional, for Fig 6 panel 2)

Step 1 — analysis and every figure input:
  Rscript analyse_v3.R
  writes: viability_final.csv, viability_labeled.csv, phase_profiles_v3.csv,
          cells_persist_base_v3.csv, fig3_panels.csv, fig5_data.csv, fig6_draws.csv,
          fig6_kernel.csv, interactions_final.csv, maineffects_eta2.csv, models_v3.txt

Step 2 — figures:
  Rscript plot_trajectories.R traj_subset_v3.csv   -> fig1_trajectories_v3.pdf
  Rscript plot_phase_profiles.R                    -> fig2_phase_profiles_v3.pdf
  Rscript plot_fig3_v3.R                           -> fig3_abundance_mediation_v3.pdf
  Rscript plot_landscape_v3.R                      -> fig4_viability_map_v3.pdf (also an alt fig3)
  Rscript plot_fig56.R                             -> fig5_homogeneous_v3.pdf, fig6_robustness_v3.pdf

Every number in methods_and_results_summary_v3 and main_messages_v3 comes from
models_v3.txt or the CSVs above. The simulator itself (simulation/) is C++; the
raw-output reduction (combine_results_v3.sh, reduce_sweep_v3.py) is bash/Python.
Everything downstream of run_summary_v3.csv is R.
