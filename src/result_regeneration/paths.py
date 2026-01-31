# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve()
REPO_ROOT = PACKAGE_ROOT.parents[2]

EXTRACTED_DATA = REPO_ROOT / "extracted_data"
FIGURES = REPO_ROOT / "figures" / "regenerated"
