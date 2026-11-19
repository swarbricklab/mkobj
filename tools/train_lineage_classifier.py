#!/usr/bin/env python3
"""
tools/train_lineage_classifier.py
CLI for training a lineage classifier offline (#44). The code lives in
workflow/scripts/lineage_train.py so the workflow's train_lineage_classifier rule
runs exactly the same thing; see that file or --help for the arguments.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from lineage_train import main  # noqa: E402

if __name__ == "__main__":
    main()
