# Code style and maintenance

This is an archive of MATLAB research experiments. Keep numerical changes separate
from formatting so a reviewer can tell when the historical computation changes.

## Source conventions

- Use UTF-8, LF line endings, four spaces per indentation level and a final newline.
- Remove trailing whitespace. Use at most one blank line between logical sections;
  avoid blank lines immediately after control-flow headers.
- Keep MATLAB function bodies at column one; indent loops and conditional bodies.
- Use spaces around binary operators and after commas. Preserve MATLAB's distinction
  between transpose (`'`) and nonconjugate transpose (`.'`).
- Prefer one executable statement per line. Keep output-suppression semicolons and
  command-form arguments unchanged when making cosmetic edits.
- Aim for 120 columns in MATLAB and 100 in C++. Wrap MATLAB expressions with `...`.
- Start function help with its name and purpose, then document array orientation,
  file inputs/outputs, side effects and relevant dependencies.
- Describe what the active code does. Label historical alternatives, incomplete
  branches and newly reconstructed adapters explicitly. Preserve original credits.
- Preserve public function names, MAT variable names, filenames and argument order.
  The legacy code uses string-based calls and `eval`; renaming requires a separate
  dependency review.
- Keep `thirdparty/` and prebuilt MEX files unchanged. Vendor upgrades should identify
  the upstream version and license and be reviewed separately.

`.editorconfig`, `miss_hit.cfg` and `.clang-format` encode these conventions.
MISS_HIT naming and copyright-generation rules are disabled to retain the historical
API and existing attribution. Its statement-ending rewrite is disabled to preserve
command-window output. Syntax and indentation checks remain enabled.

## Check and format

Install the pinned tools in a virtual environment:

```sh
python3 -m venv /tmp/lrfb-dev
/tmp/lrfb-dev/bin/python -m pip install -r requirements-dev.txt
/tmp/lrfb-dev/bin/mh_style --input-encoding utf-8 .
/tmp/lrfb-dev/bin/mh_lint --input-encoding utf-8 --json /tmp/lrfb-lint.json .
/tmp/lrfb-dev/bin/clang-format --dry-run --Werror mbs/src/*.cpp mbs/src/*.h
```

Use `mh_style --fix` or `clang-format -i` on the same targets to format changes.
Inspect the diff afterwards. The lint command currently reports the three historical
chained comparisons in the archived sources; they have not been suppressed or
silently changed.

To validate a formatting-only change against a tar.gz backup whose top-level folder
is named `LRFb`:

```sh
/tmp/lrfb-dev/bin/python tools/check_polish.py /path/to/LRFb-before-polish.tar.gz
```

This compares MATLAB syntax trees and significant tokens, C++ raw tokens, vendor
bytes and binary bytes. It requires `clang++`; it does not require MATLAB headers.
These checks establish preservation of source structure, not numerical correctness
or MEX runtime compatibility. Execute the actual MATLAB experiments before claiming
reproduction of the paper.

## Provenance

`tools/provenance.tsv` records the initial consolidation and original source MD5s.
It is a historical manifest, not a checksum of the polished working tree.
`tools/polish_manifest.tsv` records the before/after SHA-256 hashes of this polish.
Keep the original source provenance when making later changes.
