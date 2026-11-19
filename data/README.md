# Data (external)

This folder is a placeholder. Point `TOPMED_OUTLIERS_ROOT` at an analysis
tree that already contains `data/`, or symlink:

```bash
# from a clone of this code repo sitting beside a data checkout:
export TOPMED_OUTLIERS_ROOT=/path/to/pre_Watershed_scripts
# scripts in this repo can be used via PATH / by copying into that tree
```

Never commit expression matrices, underlier RDS files, or VEP caches.
