# Validation launch status

Observed **2026-09-18 04:28 UTC**; this is a launch record, not an accuracy result.

- DōgoHLA 0.1.0 is documented and pushed on `hla-experiments`.
- Reference preparation 20663742 completed successfully.
- All four fold-specific DRB1 graphs (20663743) completed successfully in 9–10 minutes.
- Six paired donor jobs are running; the remaining 26 are queued under 20663744.
- Native long-indel feasibility 20663769 passed on pilot HG02155.
- All 32 supplemental native long-indel controls are queued under 20664057.
- Four-way scoring covers native default, native long-indel, DōgoHLA without
  graph reconstruction, and full DōgoHLA, on identical input reads.
- Local monitoring service `dogohla-validation-20260918.service` is active.
- A deadline publisher is scheduled for **09:02:30 UTC / 12:02:30 Riyadh**,
  immediately after the five-hour reporting cutoff. It will push the compact
  completed-or-incomplete report, with all planned donors retained.

Current local report: `hla-spechla-pg/results/validation-20260918/analysis/REPORT.md`.
Expected GitHub report after publication: the same path on `hla-experiments`.
The monitor stops when all four methods finish or at 09:02:07 UTC. Completion
within that window is a target, not a guarantee from the scheduler. No method
superiority is inferred from job submission or the earlier selected pilot.

Validation: 35 tests passed, shell syntax checked, remote reference preparation
and graph builds completed with exit 0, and persistent monitoring plus GitHub
SSH access were verified outside the interactive session.
